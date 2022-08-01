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
