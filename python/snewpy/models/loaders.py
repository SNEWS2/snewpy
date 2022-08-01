# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.loaders`` contains classes to load core-collapse
supernova models from files stored on disk.
"""

import os
import re

from astropy import units as u
from astropy.table import Table
from astropy.io import ascii
import h5py
import numpy as np

from snewpy.models.base import PinchedModel, _GarchingArchiveModel
from snewpy.neutrino import Flavor

class Nakazato_2013(PinchedModel):
    def __init__(self, filename, metadata={}):
        """Model initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        """
        # Read FITS table using the astropy reader.
        simtab = Table.read(filename)
        self.filename = os.path.basename(filename)
        super().__init__(simtab, metadata)


class Sukhbold_2015(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """
        # Read FITS table using the astropy unified Table reader.
        simtab = Table.read(filename)
        self.filename = os.path.basename(filename)
        super().__init__(simtab, metadata)


class Tamborra_2014(_GarchingArchiveModel):
    pass


class Bollig_2016(_GarchingArchiveModel):
    pass


class Walk_2018(_GarchingArchiveModel):
    pass


class Walk_2019(_GarchingArchiveModel):
    pass


class OConnor_2015(PinchedModel):
    """Model based on the black hole formation simulation in `O'Connor (2015) <https://arxiv.org/abs/1411.7058>`_.
    """

    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        simtab = Table.read(filename,
                            names=['TIME', 'L_NU_E', 'L_NU_E_BAR', 'L_NU_X',
                                   'E_NU_E', 'E_NU_E_BAR', 'E_NU_X',
                                   'RMS_NU_E', 'RMS_NU_E_BAR', 'RMS_NU_X'],
                            format='ascii')

        header = ascii.read(simtab.meta['comments'], delimiter='=', format='no_header', names=['key', 'val'])
        tbounce = float(header['val'][0])
        simtab['TIME'] -= tbounce

        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2) / \
            (simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2) / \
            (simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2) / \
            (simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        # SYB: double-check on this factor of 4. Should be factor of 2?
        simtab['L_NU_X'] /= 4.0

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Zha_2021(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        """
        simtab = Table.read(filename,
                            names=['TIME', 'L_NU_E', 'L_NU_E_BAR', 'L_NU_X',
                                   'E_NU_E', 'E_NU_E_BAR', 'E_NU_X',
                                   'RMS_NU_E', 'RMS_NU_E_BAR', 'RMS_NU_X'],
                            format='ascii')

        header = ascii.read(simtab.meta['comments'], delimiter='=', format='no_header', names=['key', 'val'])
        tbounce = float(header['val'][0])
        simtab['TIME'] -= tbounce

        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2) / \
            (simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2) / \
            (simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2) / \
            (simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        # SYB: double-check on this factor of 4. Should be factor of 2?
        simtab['L_NU_X'] /= 4.0

        #prevent neagative lums
        simtab['L_NU_E'][simtab['L_NU_E'] < 0] = 1
        simtab['L_NU_E_BAR'][simtab['L_NU_E_BAR'] < 0] = 1
        simtab['L_NU_X'][simtab['L_NU_X'] < 0] = 1

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Warren_2020(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        """
        # Read data from HDF5 files, then store.
        f = h5py.File(filename, 'r')
        simtab = Table()

        for i in range(len(f['nue_data']['lum'])):
            if f['sim_data']['shock_radius'][i][1] > 0.00001:
                bounce = f['sim_data']['shock_radius'][i][0]
                break

        simtab['TIME'] = f['nue_data']['lum'][:, 0] - bounce
        simtab['L_NU_E'] = f['nue_data']['lum'][:, 1] * 1e51
        simtab['L_NU_E_BAR'] = f['nuae_data']['lum'][:, 1] * 1e51
        simtab['L_NU_X'] = f['nux_data']['lum'][:, 1] * 1e51
        simtab['E_NU_E'] = f['nue_data']['avg_energy'][:, 1]
        simtab['E_NU_E_BAR'] = f['nuae_data']['avg_energy'][:, 1]
        simtab['E_NU_X'] = f['nux_data']['avg_energy'][:, 1]
        simtab['RMS_NU_E'] = f['nue_data']['rms_energy'][:, 1]
        simtab['RMS_NU_E_BAR'] = f['nuae_data']['rms_energy'][:, 1]
        simtab['RMS_NU_X'] = f['nux_data']['rms_energy'][:, 1]

        simtab['ALPHA_NU_E'] = (2.0 * simtab['E_NU_E'] ** 2 - simtab['RMS_NU_E'] ** 2) / \
            (simtab['RMS_NU_E'] ** 2 - simtab['E_NU_E'] ** 2)
        simtab['ALPHA_NU_E_BAR'] = (2.0 * simtab['E_NU_E_BAR'] ** 2 - simtab['RMS_NU_E_BAR']
                                    ** 2) / (simtab['RMS_NU_E_BAR'] ** 2 - simtab['E_NU_E_BAR'] ** 2)
        simtab['ALPHA_NU_X'] = (2.0 * simtab['E_NU_X'] ** 2 - simtab['RMS_NU_X'] ** 2) / \
            (simtab['RMS_NU_X'] ** 2 - simtab['E_NU_X'] ** 2)

        # Set model metadata.
        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Kuroda_2020(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        """
        # Read ASCII data.
        simtab = Table.read(filename, format='ascii')

        # Get grid of model times.
        simtab['TIME'] = simtab['Tpb[ms]'] << u.ms
        for f in [Flavor.NU_E, Flavor.NU_E_BAR, Flavor.NU_X]:
            fkey = re.sub('(E|X)_BAR', r'A\g<1>', f.name).lower()
            simtab[f'L_{f.name}'] = simtab[f'<L{fkey}>'] * 1e51 << u.erg / u.s
            simtab[f'E_{f.name}'] = simtab[f'<E{fkey}>'] << u.MeV
            # There is no pinch parameter so use alpha=2.0.
            simtab[f'ALPHA_{f.name}'] = np.full_like(simtab[f'E_{f.name}'].value, 2.)

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)
