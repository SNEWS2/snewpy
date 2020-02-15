# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
======
snewpy
======
Front-end for supernova models which provide neutrino luminosity and spectra.
"""

from __future__ import absolute_import
from ._version import __version__
import os

src_path = os.path.realpath(__path__[0])
base_path = os.sep.join(src_path.split(os.sep)[:-2])
