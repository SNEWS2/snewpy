"""
The module ``snewpy.snowglobes_interface`` contains a low-level Python interface for SNOwGLoBES v1.2.

The class: class: `SNOwGLoBES` manages the input and output files of the SNOwGLoBES application,
and allows running the simulation for one or several setups::

    from snewpy.snowglobes_interface import SNOwGLoBES
    sng = SNOwGLoBES()

The method :meth:`SNOwGLoBES.run` performs the simulation for one or more flux files,
and returns the resulting tables as a list containing a `pandas.DataFrame`_ for each input file::

    flux_files = ['fluxes/fluence_timeBin1.dat', 'fluxes/fluence_timeBin2.dat']
    result = sng.run(flux_files, detector='icecube')

    # get results, summed over all energies and all channels:
    Ntotal_0 = results[0].smeared.weighted.sum().sum()
    Ntotal_1 = results[1].smeared.weighted.sum().sum()

    # get only sum of ibd channel
    Nibd_1 = results[1].smeared.weighted.ibd.sum()

Reading the detector and configurations, used by SNOwGLoBES::

    sng.detectors  # a table of detectors known to SNOwGLoBES
    sng.channels  # a dictionary: list of channels for each detector
    sng.efficiencies  # channel detection efficiencies for each detector

.. _pandas.DataFrame: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html

Reference:

.. autoclass:: SNOwGLoBES
   :members:

"""
from pathlib import Path
import jinja2
import pandas as pd
import numpy as np
import os

import logging
logger = logging.getLogger(__name__)

from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
import threading
import subprocess
from tqdm.auto import tqdm

def guess_material(detector):
    if detector.startswith('wc') or detector.startswith('ice'):
        mat = 'water'
    elif detector.startswith('d2O'):
        mat = 'heavywater'
    elif detector.startswith('ar'):
        mat = 'argon'
    elif detector.startswith('nova'):
        mat = 'nova_soup'
    elif detector.startswith('halo'):
        mat = 'lead'
    elif detector.startswith('scint'):
        mat = 'scint'
    else: 
        raise ValueError(f'Please provide material for {detector}')
    return mat

class SNOwGLoBES:
    def __init__(self, base_dir:Path=''):
        """ SNOwGLoBES interface 

        Parameters
        ----------
        base_dir: Path or None
            Path to the SNOwGLoBES installation.
            If empty, try to get it from $SNOWGLOBES environment var


        On construction SNOwGLoBES will read: 

        * detectors from `<base_dir>/detector_configurations.dat`,
        * channels  from `<base_dir>/channels/channel_*.dat`
        * efficiencies from `<base_dir>/effic/effic_*.dat`

        After that use :meth:`SNOwGLoBES.run` method to run the simulation for specific detector and flux file.
        """
        if not base_dir:
            base_dir = os.environ['SNOWGLOBES']
        self.base_dir = Path(base_dir)
        self._load_detectors(self.base_dir/'detector_configurations.dat')
        self._load_channels(self.base_dir/'channels')
        self._load_efficiencies(self.base_dir/'effic')
        self._load_smearing(self.base_dir/'smear')

        env = jinja2.Environment(loader=jinja2.PackageLoader('snewpy'), enable_async=True)
        self.template = env.get_template('supernova.glb')

    def _load_detectors(self, path:Path):
        df = pd.read_table(path,names=['name','mass','factor'], delim_whitespace=True, comment='#')
        df['tgt_mass']=df.mass*df.factor
        self.detectors=df.set_index('name').T
        logger.info(f'read masses for detectors {list(self.detectors)}')
        logger.debug(f'detectors: {self.detectors}')
       
    def _load_channels(self, chan_dir):

        def _read_binning(fname):
            with open(fname) as f:
                l = f.readline().strip()
                if not l.startswith('%'):
                    l = '% 200 0.0005 0.100 200 0.0005 0.100'
                tokens = l.split(' ')[1:]
                return dict(zip(['nsamples','smin','smax','nbins','emin','emax'],tokens))

        self.channels = {}
        self.binning = {}
        all_channels = []
        for f in chan_dir.glob('channels_*.dat'):
            material = f.stem[len('channels_'):]
            df = pd.read_table(f,delim_whitespace=True, index_col=1, comment='%', names=['name','n','parity','flavor','weight'])
            self.channels[material] = df
            self.binning[material] = _read_binning(f)
        self.materials = list(self.channels.keys())
        self.chan_dir = chan_dir
        logger.info(f'read channels for  materials: {self.materials}')
        logger.debug(f'channels: {self.channels}')
        
    def _load_efficiencies(self, path):
        result = {}
        for detector in self.detectors:
            res_det = {}
            for file in path.glob(f'effic_*_{detector}.dat'):
                channel =  file.stem[len('effic_'):-len(detector)-1]
                logger.debug(f'Reading file ({detector},{channel}): {file}')
                with open(file) as f:
                    res_det[channel]= f.read()
            result[detector]=res_det
        self.efficiencies = result 
        logger.info(f'read efficiencies for materials: {list(self.efficiencies.keys())}')
        logger.debug(f'efficiencies: {self.efficiencies}')

    def _load_smearing(self, path):
        result = {}
        for detector in self.detectors:
            res_det = {}
            for file in path.glob(f'smear*_{detector}.dat'):
                channel =  file.stem[len('smear_'):-len(detector)-1]
                res_det[channel]= file
            result[detector]=res_det
        self.smearings = result 
        logger.info(f'read efficiencies for materials: {list(self.efficiencies.keys())}')
        logger.debug(f'efficiencies: {self.efficiencies}')

       
    def run(self, flux_files, detector:str, material:str=None):
        """ Run the SNOwGLoBES simulation for given configuration,
        collect the resulting data and return it in `pandas.DataFrame`

        Parameters
        -----------
        flux_files: list(str) or str
            An iterable of flux table filenames to process, or a single filename
        detector: str
            Detector name, known to SNOwGLoBES
        material: str or None
            Material name, known to SNOwGLoBES. If None, we'll try to guess it

        Returns
        --------
        list(pd.DataFrame or Exception)
            List with the data table for each flux_file, keeping the order.
            Each table containing Energy (GeV) as index values, 
            and number of events for each energy bin, for all interaction channels.
            Columns are hierarchical: (is_weighted, is_smeared, channel),
            so one can easily access the desired final table. 
            If run failed with exception, this exception will be returned (not raised).

        Raises
        ------
        ValueError
            if material or detector value is invalid
        RuntimeError
            if SNOwGLoBES run has failed

        """
        if not  detector in self.detectors:
            raise ValueError(f'Detector "{detector}" is not in {list(self.detectors)}')
        if material is None:
            material = guess_material(detector)
        if not material in self.materials:
            raise ValueError(f'Material "{material}" is not in {self.materials}')
        if not self.efficiencies[detector]:
            logger.warning(f'Missing efficiencies for detector={detector}! Results will assume 100% efficiency')
        if not self.smearings[detector]:
            logger.warning(f'Missing smearing for detector={detector}! Results will not be smeared')
        if isinstance(flux_files,str):
            flux_files = [flux_files]
        
        with tqdm(total=len(flux_files), leave=False, desc='Flux files') as progressbar:
            def do_run(file):
                result =  Runner(self,Path(file),detector,material).run()
                progressbar.update()
                return result

            self.lock = threading.Lock() #global lock, ensuring that snowglobes files aren't mixed!
            with ThreadPoolExecutor() as executor:
                result = executor.map(do_run, flux_files)
                return list(result)

@dataclass
class Runner:
    sng: SNOwGLoBES
    flux_file: Path
    detector: str
    material: str
    
    def __post_init__(self):
        self.channels=self.sng.channels[self.material]
        self.binning=self.sng.binning[self.material]
        self.efficiency=self.sng.efficiencies[self.detector]
        self.smearing=self.sng.smearings[self.detector]
        self.det_config=self.sng.detectors[self.detector]
        self.base_dir=self.sng.base_dir
        self.out_dir=self.base_dir/'out'
            
    def _generate_globes_config(self):
        cfg =  self.sng.template.render(flux_file=self.flux_file.resolve(),
                                    detector=self.detector,
                                    target_mass=self.det_config.tgt_mass,
                                    smearing=self.smearing,
                                    xsec_dir =self.base_dir/'xscns',
                                    channels =list(self.channels.itertuples()),
                                    efficiency =self.efficiency,
                                    **self.binning
                                    )
        return cfg

    def _parse_output(self, output):
        data = {}
        for l in output:
            #read the generated files from output
            if l.endswith('weighted.dat'):
                channum, fname = l.split()
                fname = self.base_dir/fname
                #load the data from file
                try:
                    E,N = np.loadtxt(fname, comments=['--','Total'], unpack=True)
                    channel=self.channels.loc[int(channum)]
                    smeared= 'smeared'  if '_smeared' in fname.stem else 'unsmeared'
                    data[(channel['name'],smeared,'unweighted')] = N
                    data[(channel['name'],smeared,'weighted')] = N*channel['weight']
                except ValueError:
                    logger.error(f'Failed reading data from file {fname}')
                finally:
                    fname.unlink() #cleanup file
        #collect everything to pandas DataFrame
        df = pd.DataFrame(data, index = E)
        df.index.rename('E', inplace=True)
        df.columns.rename(['channel','is_smeared','is_weighted'], inplace=True)
        return df.reorder_levels([2,1,0], axis='columns')

    def run(self):
        """write configuration file and run snowglobes"""
        cfg = self._generate_globes_config()
        chan_file = self.sng.chan_dir/f'channels_{self.material}.dat'
        #this section is exclusive to one process at a time, 
        # because snowglobes must  read the "$SNOGLOBES/supernova.glb" file
        with self.sng.lock:
            #write configuration file:
            with open(self.base_dir/'supernova.glb','w') as f:
                f.write(cfg)
            #run the snowglobes process:
            p = subprocess.run(['bin/supernova', self.flux_file.stem, str(chan_file), self.detector],
                capture_output=True,
                cwd=self.base_dir)
        #process the output
        stdout = p.stdout.decode('utf_8')
        stderr = p.stderr.decode('utf_8')
        if(stderr):
            logger.error('Run errors: \n'+stderr)
        if(p.returncode==0):
            tables = self._parse_output(stdout.split('\n'))
            return tables
        else:
            raise RuntimeError('SNOwGLoBES run failed:\n'+stderr)

class SimpleRate(SNOwGLoBES):
    def __init__(self, base_dir:Path=''):
        """ Simple rate calculation interface 
        Computes expected rate for a perfect detector (100% efficiencies, no smearing)
        without using GLOBES. The formula for the rate is
                Rate = [cross-section in 10^-38 cm^2] x 10^-38 x [fluence in cm^-2] x [target mass in kton] 
                    x [Dalton per kton] x [energy bin size in GeV]
        with [target mass in kton] x [Dalton per kton] = number of reference targets in experiment.

        Parameters
        ----------
        base_dir: Path or None
            Path to the directory where the cross-section, detector, and channel files are located
            If empty, try to get it from $SNOWGLOBES environment var

        On construction the code will read: 

        * detectors from `<base_dir>/detector_configurations.dat`,
        * channels  from `<base_dir>/channels/channel_*.dat`

        After that use :meth:`SimpleRate.run` method to run the simulation for specific detector and flux file.
        """
        if not base_dir:
            base_dir = os.environ['SNOWGLOBES']
        self.base_dir = Path(base_dir)
        self._load_detectors(self.base_dir/'detector_configurations.dat')
        self._load_channels(self.base_dir/'channels')

    def _compute_rates(self, detector, material, flux_file:Path):
        flux_file = flux_file.resolve()
        fluxes = np.loadtxt(flux_file)
        TargetMass = self.detectors[detector].tgt_mass
        data = {}
        energies = np.linspace(7.49e-4, 9.975e-2, 200) # Use the same energy grid as SNOwGLoBES
        for chan_num,channel in enumerate(self.channels[material].itertuples()):
            xsec_path = f"xscns/xs_{channel.name}.dat"
            xsec = np.loadtxt(self.base_dir/xsec_path)
            flavor_index = 0 if 'e' in channel.flavor else (1 if 'm' in channel.flavor else 2)
            flavor = flavor_index + (3 if channel.parity == '-' else 0)
            flux = fluxes[:, (0,1+flavor)]
            binsize = energies[1] - energies[0]
            # Cross-section in 10^-38 cm^2
            xsecs = np.interp(np.log(energies)/np.log(10), xsec[:, 0], xsec[:, 1+flavor], left=0, right=0) * energies
            # Fluence (flux integrated over time bin) in cm^-2 
            # (must be divided by 0.2 MeV to compensate the multiplication in generate_time_series)
            fluxs = np.interp(energies, flux[:, 0], flux[:, 1], left=0, right=0)/2e-4
            # Rate computation
            rates = xsecs * 1e-38 * fluxs * float(TargetMass) * 1./1.661e-33 * binsize
            # Weighting
            weighted_rates = rates * channel.weight
            # Write to dictionary
            data[(channel.name,'unsmeared','unweighted')] = rates
            data[(channel.name,'unsmeared','weighted')] = weighted_rates
        #collect everything to pandas DataFrame
        df = pd.DataFrame(data, index = energies)
        df.index.rename('E', inplace=True)
        df.columns.rename(['channel','is_smeared','is_weighted'], inplace=True)
        return df.reorder_levels([2,1,0], axis='columns')
       
    def run(self, flux_files, detector:str, material:str=None):
        """ Compute expected rates for given configuration,
        collect the resulting data and return it in `pandas.DataFrame`

        Parameters
        -----------
        flux_files: list(str) or str
            An iterable of flux table filenames to process, or a single filename
        detector: str
            Detector name, SNOwGLoBES style
        material: str or None
            Material name, SNOwGLoBES style. If None, we'll try to guess it

        Returns
        --------
        list(pd.DataFrame or Exception)
            List with the data table for each flux_file, keeping the order.
            Each table containing Energy (GeV) as index values, 
            and number of events for each energy bin, for all interaction channels.
            Columns are hierarchical: (is_weighted, channel),
            so one can easily access the desired final table. 
            If run failed with exception, this exception will be returned (not raised).

        Raises
        ------
        ValueError
            if material or detector value is invalid

        """
        if not  detector in self.detectors:
            raise ValueError(f'Detector "{detector}" is not in {list(self.detectors)}')
        if material is None:
            material = guess_material(detector)
        if not material in self.materials:
            raise ValueError(f'Material "{material}" is not in {self.materials}')
        if isinstance(flux_files,str):
            flux_files = [flux_files]
        
        with tqdm(total=len(flux_files), leave=False, desc='Flux files') as progressbar:
            results = []
            for flux_file in flux_files:
                result = self._compute_rates(detector,material,Path(flux_file))
                progressbar.update()
                results.append(result)
            return results
