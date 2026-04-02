# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.typeIa_loaders`` contains classes to load type Ia supernova
models from files stored on disk.
"""

# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.pisn_loaders`` contains classes to load pair-instability supernova
models from files stored on disk.
"""
import logging
import os
import numpy as np
import pandas as pd
import tarfile

from astropy import units as u
from astropy.table import Table

from scipy import interpolate 

from snewpy.models.base import SupernovaModel
from snewpy.flavor import ThreeFlavor
from snewpy.flux import Spectrum
from snewpy import _model_downloader

class TypeIa(SupernovaModel):
    """The DDT model is described in 'Neutrinos from type Ia supernovae: The deflagration-to-detonation transition scenario', by Warren P. Wright et al.,
    [Phys. Rev. D94 (2016) 025026](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.025026), [arXiv:1605.01408](https://arxiv.org/abs/1605.01408).  
    There are 30 snapshots in time and the format of each data file is the SNOwGLoBES format.  

    The GCD model is described in 'Neutrinos from type Ia supernovae: The gravitationally confined detonation scenario' by Warren P. Wright et al.,
    [Phys. Rev. D95 (2017) 043006](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.043006), [arXiv:1609.07403](https://arxiv.org/abs/1609.07403).  
    There are 64 snapshots in time and the format of each data file is the SNOwGLoBES format.
    """

    def  __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        tarfilename: str
            Absolute or relative path to tar archive
        """
        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)

        tf = tarfile.open(datafile)
        
        # Find the "NoOsc" files.
        datafiles = sorted([f.name for f in tf if '.dat' in f.name])
        nooscfiles = [df for df in datafiles if 'NoOsc' in df]
        nooscfiles.sort(key=len)

        # Loop through the NoOsc files and pull out the number fluxes.
        self.time = []
        self.energy = None   
        self.initial_spectra = {}
        self.interpolation = {}         
        
        self._flavorkeys = {ThreeFlavor.NU_E: 'NuE',
                            ThreeFlavor.NU_E_BAR: 'aNuE',
                            ThreeFlavor.NU_MU: 'NuMu',
                            ThreeFlavor.NU_MU_BAR: 'aNuMu',
                            ThreeFlavor.NU_TAU: 'NuTau',
                            ThreeFlavor.NU_TAU_BAR: 'aNuTau'}        

        for nooscfile in nooscfiles:
            with tf.extractfile(nooscfile) as f:
                logging.debug('Reading {}'.format(nooscfile))
                meta = f.readline()
                metatext = meta.decode('utf-8')
                t = float(metatext.split('TBinMid=')[-1].split('sec')[0]) 
                dt = float(metatext.split('tBinWidth=')[-1].split('s')[0]) 
                dE = float(metatext.split('eBinWidth=')[-1].split('MeV')[0])

                data = Table.read(f, format='ascii.commented_header', header_start=-1)
                data.meta['t'] = t
                data.meta['dt'] = dt
                data.meta['dE'] = dE

                self.time.append(t)
                
                if self.energy is None:
                    self.energy = (data['E(GeV)'].data*1000).tolist()

            for flavor in ThreeFlavor:
                key = self._flavorkeys[flavor]
                # convert from flux back to initial spectra: number per /s/erg                
                spectrum = data[key].data.tolist() * (4*np.pi*(10*u.kpc)**2) /dt/dE 
                if flavor in self.initial_spectra:
                    self.initial_spectra[flavor].append(spectrum)
                else:
                    self.initial_spectra[flavor] = [spectrum]                

        for flavor in ThreeFlavor:
            self.interpolation[flavor] = interpolate.RegularGridInterpolator((self.time, self.energy), self.initial_spectra[flavor], method='cubic')
            
        self.time *= u.s
        self.energy *= u.MeV
                            
        self.filename = os.path.basename(filename)            

    def _get_initial_spectra_dict(self, t, E, flavors=ThreeFlavor):
        """Get neutrino spectra/luminosity curves after oscillation.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity or ndarray of astropy.Quantity
            Energies to evaluate the initial spectra.
        flavors: iterable of snewpy.neutrino.Flavor
            Return spectra for these flavors only (default: all)

        Returns
        -------
        spectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """   
        #convert input arguments to 1D arrays
        t = u.Quantity(t, ndmin=1)
        E = u.Quantity(E, ndmin=1)

        initial_spectra = {}
        for flavor in ThreeFlavor:
            initial_spectra[flavor] = self.interpolation[flavor]((t, E)) / (u.MeV * u.s)

        return initial_spectra





