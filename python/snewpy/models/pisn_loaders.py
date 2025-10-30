# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.pisn_loaders`` contains classes to load pair-instability supernova
models from files stored on disk.
"""

import numpy as np
import pandas as pd
import tarfile

from scipy.interpolate import interp1d

from astropy import units as u
from astropy.table import Table

from snewpy.models.base import SupernovaModel
from pathlib import Path

class Wright_2017(SupernovaModel):
    """PISN model described in the paper `Neutrino signal from pair-instability supernovae' by Wright et al.,  
    Phys. Rev. D96 (2017) 103008  <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.103008>`_
    """

    def  __init__(self, tarfilename, metadata={}):
        """
        Parameters
        ----------
        tarfilename: str
            Absolute or relative path to tar archive
        """
        self.tfname = tarfilename
        tf = tarfile.open(self.tfname)
        
        # Pull out the "NoOsc" files.
        datafiles = sorted([f.name for f in tf if '.dat' in f.name])
        nooscfiles = [df for df in datafiles if 'NoOsc' in df]
        nooscfiles.sort(key=len)

        # Loop through the noosc files and pull out the number fluxes.
        self.time = []
        self.energy = None
        flux = {}
        self.initialspectra = {}
        self.fmin = 1e99
        self.fmax = -1e99

        # Extract the data: The ordering is rows=times and columns=energies.
        for nooscfile in nooscfiles:
            with tf.extractfile(nooscfile) as f:
                meta = f.readline()
                metatext = meta.decode('utf-8')
                t = float(metatext.split('TBinMid=')[-1].split('sec')[0])
                dt = float(metatext.split('tBinWidth=')[-1].split('s')[0])
                dE = float(metatext.split('eBinWidth=')[-1].split('MeV')[0])

                # the data in the table is number of neutrinos at Earth per cm^2 in a given time and energy bin 
                data = Table.read(f, format='ascii.commented_header', header_start=-1)
                data.meta['t'] = t
                data.meta['dt'] = dt
                data.meta['dE'] = dE

                self.time.append(t)
                
                if self.energy is None:
                    self.energy = (data['E(GeV)'].data*1000).tolist()

            for flavor in ['NuE', 'NuMu', 'NuTau', 'aNuE', 'aNuMu', 'aNuTau']:
                if flavor in flux:
                    flux[flavor].append(data[flavor].data.tolist())
                else:
                    flux[flavor] = [data[flavor].data.tolist()]
                    
            # convert from flux back to initial spectra: number per /s/erg
            self.initialspectra[flavor] = flux[flavor] * (4*np.pi*(10*u.kpc)**2) /dt /(dE*u.MeV) 

        # Transpose so that rows=energy and cols=time.
        for k, v in self.initialspectra.items():
            self.initialspectra[k] = np.transpose(self.initialspectra[k])


    def get_initial_spectra(self, t, E):
        E = E.to_value('MeV')
        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to_value('s')    
        
        for flavor in flavors:
        # Use np.interp rather than scipy.interpolate.interp1d because it
        # can handle dimensional units (astropy.Quantity).
            result = get_value(np.interp(t, self.time, self.initialspectrum[flavor]))  
            initialspectra[flavor] = result
            
        return initialspectra





