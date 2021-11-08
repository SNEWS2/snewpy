"""
Python interface for `SNOwGLoBES` v1.2

Usage:
>>> from snewpy.snowglobes_interface import SNOwGLoBES
>>> sng = SNOwGLoBES() 
>>> result = sng.run('./Bollig_2016_s11.2c_AdiabaticMSW_NMO.dat', detector='icecube', material='water')
>>> results.smeared.weighted.sum().sum() #get results, summed over all energies and all channels:
320622.97449880163
"""
from pathlib import Path
import jinja2
import pandas as pd
import numpy as np
import os

import logging
logger = logging.getLogger(__name__)

from dataclasses import dataclass
import subprocess
import asyncio

class SNOwGLoBES:
    def __init__(self, base_dir:Path=''):
        """ SNOwGLoBES interface

        Args:
            base_dir
                Path to the SNOwGLoBES installation
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

        env = jinja2.Environment(loader=jinja2.PackageLoader('snewpy'))
        self.template = env.get_template('supernova.glb')

    def _load_detectors(self, path:Path):
        df = pd.read_table(path,names=['name','mass','factor'], delim_whitespace=True, comment='#')
        df['tgt_mass']=df.mass*df.factor
        self.detectors=df.set_index('name').T
        logger.info(f'read masses for detectors {list(self.detectors)}')
        logger.debug(f'detectors: {self.detectors}')
       
    def _load_channels(self, chan_dir):
        self.channels = {}
        all_channels = []
        for f in chan_dir.glob('channels_*.dat'):
            material = f.stem[len('channels_'):]
            df = pd.read_table(f,delim_whitespace=True, index_col=1, comment='%', names=['name','n','parity','flavor','weight'])
            self.channels[material] = df
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
       
    def run(self, flux_files, detector:str, material:str) -> pd.DataFrame:
        """ Run the SNOwGLoBES simulation for given configuration,
        collect the resulting data and return it in `pandas.DataFrame`

        Args:
            flux_files
                An iterable of flux table filenames to process, or a single filename
            detector
                Detector name, known to SNOwGLoBES
            material
                Material name, known to SNOwGLoBES

        Returns:
            pd.DataFrame
                The table, containing Energy (GeV) as index values, 
                and number of events for each energy bin, for all interaction channels.
                Columns are hierarchical: (is_weighted, is_smeared, channel),
                so one can easily access :code:`data.weighted.unsmeared.ibd` 
        """
        if not material in self.materials:
            raise ValueError(f'material "{material}" is not in {self.materials}')
        if not  detector in self.detectors:
            raise ValueError(f'detector "{detector}" is not in {list(self.detectors)}')
        if isinstance(flux_files,str):
            flux_files = [flux_files]

        #make a function to collect results from asyncio
        async def do_run():
            self.lock = asyncio.Lock() #global lock, ensuring that snowglobes files aren't mixed!
            return await self.run_async(flux_files,detector,material)

        return asyncio.run(do_run())
 
    async def run_async(self, flux_files, detector:str, material:str):
        tasks = [Runner(self,Path(f),detector,material).run() for f in flux_files]
        results = await asyncio.gather(*tasks, return_exceptions=False)
        return results
        
@dataclass
class Runner:
    sng: SNOwGLoBES
    flux_file: Path
    detector: str
    material: str
    
    def __post_init__(self):
        self.channels=self.sng.channels[self.material]
        self.efficiency=self.sng.efficiencies[self.detector]
        self.det_config=self.sng.detectors[self.detector]
        self.base_dir=self.sng.base_dir
        self.out_dir=self.base_dir/'out'
        if not self.efficiency:
            logger.warning(f'Missing efficiencies for detector={self.detector}!')
            
    def _generate_globes_config(self):
        cfg =  self.sng.template.render(flux_file=self.flux_file.resolve(),
                                    detector=self.detector,
                                    target_mass=self.det_config.tgt_mass,
                                    smear_dir=self.base_dir/'smear',
                                    xsec_dir =self.base_dir/'xscns',
                                    channels =list(self.channels.itertuples()),
                                    efficiency =self.efficiency)
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

    async def run(self):
        """write configuration file and run snowglobes"""
        cfg = self._generate_globes_config()
        chan_file = self.sng.chan_dir/f'channels_{self.material}.dat'
        #this section is exclusive to one process at a time, 
        # because snowglobes must  read the "$SNOGLOBES/supernova.glb" file
        async with self.sng.lock:
            #write configuration file:
            with open(self.base_dir/'supernova.glb','w') as f:
                f.write(cfg)
            #run the snowglobes process:
            p = await asyncio.create_subprocess_exec('bin/supernova', self.flux_file.stem, chan_file, self.detector,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=self.base_dir)
            #collect the output
            stdout,stderr = await p.communicate()
        #process the output
        stdout = stdout.decode('utf_8')
        stderr = stderr.decode('utf_8')
        if(stderr):
            logger.error('Run failed: \n'+stderr)
        if(p.returncode==0):
            tables = self._parse_output(stdout.split('\n'))
            return tables
        
 
