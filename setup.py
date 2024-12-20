#!/usr/bin/env python
#
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from setuptools import setup, find_packages


#
# Begin setup
#
setup_keywords = dict()
#
# Set other keywords for the setup function.
#
# Use entry_points to let `pip` create executable scripts for each target platform.
# See https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
# setup_keywords['entry_points'] = {'console_scripts': ['to_snowglobes = snewpy.to_snowglobes:generate_time_series', ], },
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['packages'] = find_packages('python')
setup_keywords['package_dir'] = {'': 'python'}
setup_keywords['package_data'] = {'':['models/*.yml']}

#
# Run setup command.
#
setup(**setup_keywords)
