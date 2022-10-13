"""
The module ``snewpy.snowglobes_interface`` contains a low-level Python interface for SNOwGLoBES v1.2.

.. note::
    Users should only use the high-level interface described above.
    This low-level interface is not guaranteed to be stable and may change at
    any time without warning, e.g. to support new SNOwGLoBES versions.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import os

import logging
logger = logging.getLogger(__name__)

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

    if '_he' in detector:
        mat = mat+'_he'
    return mat


class SimpleRate():
    def __init__(self, detectors:str="all", detector_effects=True, base_dir:Path=''):
        """Simple rate calculation interface.
        Computes expected rate for a detector without using GLOBES. The formula for the rate is

                Rate = Sum_i ([cross-section in 10^-38 cm^2] x 10^-38 x [fluence in cm^-2])_i x [smearing matrix]_{ij}
                       x [target mass in kton] x [Dalton per kton] x [energy bin size in GeV] x [efficiency]

        with [target mass in kton] x [Dalton per kton] = number of reference targets in experiment.

        On construction the code will read: 

        * detectors from `<base_dir>/detector_configurations.dat`,
        * channels from `<base_dir>/channels/channel_*.dat`

        if the detector_effects option is on, also read efficiencies and smearing:

        * efficiencies from `<base_dir>/effic/effic_*.dat`,
        * smearing matrices from `<base_dir>/smear/smear_*.dat`

        After that use :meth:`SimpleRate.run` method to run the simulation for specific detector and flux file.

        Parameters
        ----------
        detectors: str
            Name of detector. If ``"all"``, will use all detectors supported by SNOwGLoBES.

        detector_effects: bool
            If true, account for efficiency and smearing. If false, consider a perfect detector.

        base_dir:         Path or None
            Path to the directory where the cross-section, detector, and channel files are located
            If empty, try to get it from ``$SNOWGLOBES`` environment var
        """
        if not base_dir:
            base_dir = os.environ['SNOWGLOBES']
        self.base_dir = Path(base_dir)
        self._load_detectors(self.base_dir/'detector_configurations.dat', detectors)
        self._load_channels(self.base_dir/'channels')
        self.efficiencies = None
        self.smearings = None
        if detector_effects:
            self._load_efficiency_vectors(self.base_dir/'effic')
            self._load_smearing_matrices(self.base_dir/'smear')

    def _load_detectors(self, path:Path, detectors:str):
        df = pd.read_table(path,names=['name','mass','factor'], delim_whitespace=True, comment='#')
        df['tgt_mass']=df.mass*df.factor

        if detectors == 'all':
            detectors = list(df['name'])
            detectors.remove('d2O')
        elif isinstance(detectors, str):
            detectors = [detectors]

        self.detectors = df.set_index('name').T.filter(items=detectors)
        logger.info(f'read masses for detectors {list(self.detectors)}')
        logger.debug(f'detectors: {self.detectors}')
       
    def _load_channels(self, chan_dir):

        def _read_binning(fname):
            with open(fname) as f:
                l = f.readline().strip()
                if not l.startswith('%'):
                    l = '% 200 0.0005 0.100 200 0.0005 0.100'
                tokens = l.split(' ')[1:]
                nsamples, smin,smax, nbins,emin,emax = [float(t) for t in tokens]
                return {'e_true' :np.linspace(smin,smax,int(nsamples)+1),
                        'e_smear':np.linspace(emin,emax,int(nbins)+1)
                        }

        self.channels = {}
        self.binning = {}
        for f in chan_dir.glob('channels_*.dat'):
            material = f.stem[len('channels_'):]
            df = pd.read_table(f,delim_whitespace=True, index_col=1, comment='%', names=['name','n','parity','flavor','weight'])
            self.channels[material] = df
            self.binning[material] = _read_binning(f)
        self.materials = list(self.channels.keys())
        self.chan_dir = chan_dir
        logger.info(f'read channels for materials: {self.materials}')
        logger.debug(f'channels: {self.channels}')

    def _load_efficiency_vectors(self, path):
        result = {}
        for detector in self.detectors:
            res_det = {}
            for file in path.glob(f'effic_*_{detector}.dat'):
                channel =  file.stem[len('effic_'):-len(detector)-1]
                logger.debug(f'Reading file ({detector},{channel}): {file}')
                with open(file) as f:
                    effs = np.fromiter(f.readlines()[0].split("{")[-1].split("}")[0].split(","), float)
                    res_det[channel]= effs
            result[detector]=res_det
        self.efficiencies = result 
        logger.info(f'read efficiencies for detectors: {list(self.efficiencies.keys())}')
        logger.debug(f'efficiencies: {self.efficiencies}')

    def _load_smearing_matrices(self, path):
        result = {}
        for detector in self.detectors:
            res_det = {}
            for file in path.glob(f'smear*_{detector}.dat'):
                channel =  file.stem[len('smear_'):-len(detector)-1]
                with open(file) as f:
                    lines = f.readlines()[1:-1]
                    while not "{" in lines[0]: lines = lines[1:]
                    while not "{" in lines[-1]: lines = lines[:-1]
                    matrix = np.zeros((len(lines),len(lines)))
                    for i,l in enumerate(lines):
                        elements = np.fromiter(l.split("{")[-1].split("}")[0].split(","), float)
                        matrix[i, int(elements[0]+0.1):int(elements[1]+0.1)+1] = elements[2:]
                    res_det[channel]= matrix
            result[detector]=res_det
        self.smearings = result 
        logger.info(f'read smearing matrices for detectors: {list(self.smearings.keys())}')
        logger.debug(f'smearing matrices: {self.smearings}')

    def compute_rates(self, detector:str, material:str, fluxes:np.ndarray, flux_energies:np.ndarray):
        """ Calculate the rates for the given neutrino fluxes interacting in the given detector.

        Parameters
        ----------
        detector: str
            Detector name in SNOwGLoBES. Check `SimpleRate.detectors` for the list of options
        material: str
            Material name in SNOwGLoBES. Check `SimpleRate.materials` for the list of options
        fluxes:
            2d array of the neutrino fluxes.
            First array dimension corresponds to flavors [nu_e, nu_mu, anti_nu_e, anti_nu_mu]
            Second array dimension corresponds to the `energies` bins,

        flux_energies:
            1d array of the neutrino energy bins in GeV, corresponding to the `fluxes`

        Returns
        -------
        pd.DataFrame
            Table with Energy (GeV) as index values,
            and number of events for each energy bin, for all interaction channels.
            Columns are hierarchical: (is_weighted, channel),
            so one can easily access the desired final table. 
        """
        TargetMass = self.detectors[detector].tgt_mass
        data = {}

        #load the binning for the smearing
        binning= self.binning[material]
        #calculate bin centers 
        energies_t = 0.5*(binning['e_true'][1:]+ binning['e_true'][:-1] )
        energies_s= 0.5*(binning['e_smear'][1:]+binning['e_smear'][:-1])
        binsize = np.diff(binning['e_true'])

        for channel in self.channels[material].itertuples():
            xsec_path = f"xscns/xs_{channel.name}.dat"
            xsec = np.loadtxt(self.base_dir/xsec_path)
            flavor_index = 0 if 'e' in channel.flavor else (1 if 'm' in channel.flavor else 2)
            flavor = flavor_index + (3 if channel.parity == '-' else 0)
            flux = fluxes[flavor]
            # Cross-section in 10^-38 cm^2
            xsecs = np.interp(np.log(energies_t)/np.log(10), xsec[:, 0], xsec[:, 1+flavor], left=0, right=0) * energies_t
            # Fluence (flux integrated over time bin) in cm^-2 
            # (must be divided by 0.2 MeV to compensate the multiplication in generate_time_series)
            fluxs = np.interp(energies_t, flux_energies, flux, left=0, right=0)/2e-4
            # Rate computation
            rates = xsecs * 1e-38 * fluxs * float(TargetMass) * 1./1.661e-33 * binsize
            # Weighting
            weighted_rates = rates * channel.weight
            # Write to dictionary
            data[(channel.name,'unsmeared','unweighted')] = rates
            data[(channel.name,'unsmeared','weighted')] = weighted_rates
            # Add detector effects
            if self.smearings and self.efficiencies:
                smear = self.smearings[detector].get(channel.name, np.eye(len(rates)))
                effic = self.efficiencies[detector].get(channel.name, np.ones(len(rates)))
                rates = np.dot(smear,rates) * effic
                weighted_rates = rates * channel.weight
                # Write to dictionary
                data[(channel.name,'smeared','unweighted')] = rates
                data[(channel.name,'smeared','weighted')] = weighted_rates
        #collect everything to pandas DataFrame
        df = pd.DataFrame(data, index = energies_s)
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
        if not detector in self.detectors:
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
                e_flux = np.loadtxt(Path(flux_file).resolve()).T
                energies = e_flux[0]
                fluxes = e_flux[1:]
                result = self.compute_rates(detector, material, fluxes, energies)
                progressbar.update()
                results.append(result)
            return results
