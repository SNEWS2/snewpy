# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.loaders`` contains classes to load core-collapse
supernova models from files stored on disk.
"""

import os

from astropy import units as u
from astropy.table import Table

from snewpy.models.base import PinchedModel

class Nakazato_2013(PinchedModel):
    """Model based on simulations from Nakazato et al., ApJ S 205:2
    (2013), ApJ 804:75 (2015), PASJ 73:639 (2021). See also http://asphwww.ph.noda.tus.ac.jp/snn/.
    """
    def __init__(self, filename):
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
        # Store model metadata.
        if 't_rev' in filename:
            self.progenitor_mass = float(filename.split('-')[-1].strip('s%.fits')) * u.Msun
            self.revival_time = float(filename.split('-')[-2].strip('t_rev%ms')) * u.ms
            self.metallicity = float(filename.split('-')[-3].strip('z%'))
            self.EOS = filename.split('-')[-4].upper()
        # No revival time because the explosion "failed" (BH formation).
        else:
            self.progenitor_mass = float(filename.split('-')[-1].strip('s%.fits')) * u.Msun
            self.metallicity = float(filename.split('-')[-2].strip('z%'))
            self.revival_time = 0 * u.ms
            self.EOS = filename.split('-')[-4].upper()

        metadata = {
            'Progenitor mass': self.progenitor_mass,
            'EOS': self.EOS,
            'Metallicity': self.metallicity,
            'Revival time': self.revival_time
        }

        # Read FITS table using the astropy reader.
        simtab = Table.read(filename)
        self.filename = os.path.basename(filename)
        super().__init__(simtab, metadata)
