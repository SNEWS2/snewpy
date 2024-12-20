#!/usr/bin/env python
#
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from setuptools import setup, find_packages

# Git-based version info. Remove?
from python.snewpy._git import get_version, SetVersion

#
# Begin setup
#
setup_keywords = dict()
#
setup_keywords['version'] = get_version()
#
# Set other keywords for the setup function.
#
# Use entry_points to let `pip` create executable scripts for each target platform.
# See https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
# setup_keywords['entry_points'] = {'console_scripts': ['to_snowglobes = snewpy.to_snowglobes:generate_time_series', ], },
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages('python')
setup_keywords['package_dir'] = {'': 'python'}
setup_keywords['package_data'] = {'':['templates/*.glb', 'models/*.yml']}
setup_keywords['cmdclass'] = {'version': SetVersion}
setup_keywords['test_suite']='snewpy.test.snewpy_test_suite.snewpy_test_suite'

#
# Internal data directories.
#
#setup_keywords['data_files'] = [('snewpy/data/config', glob('data/config/*')),
#                                ('snewpy/data/spectra', glob('data/spectra/*'))]
#
# Run setup command.
#
setup(**setup_keywords)
