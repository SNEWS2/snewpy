#!/usr/bin/env python
#
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages

# Git-based version info. Remove?
from python.snewpy._git import get_version, SetVersion

#
# Begin setup
#
setup_keywords = dict()
#
setup_keywords['name'] = 'snewpy'
setup_keywords['description'] = 'A Python package for working with supernova neutrinos'
setup_keywords['author'] = 'SNEWS Collaboration'
setup_keywords['author_email'] = 'snews2.0@lists.bnl.gov'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/SNEWS2/snewpy'
setup_keywords['version'] = get_version()
#
# Use README.md as a long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
    setup_keywords['long_description_content_type'] = 'text/markdown'
#
# Set other keywords for the setup function.
#
# Use entry_points to let `pip` create executable scripts for each target platform.
# See https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
# setup_keywords['entry_points'] = {'console_scripts': ['to_snowglobes = snewpy.to_snowglobes:generate_time_series', ], },
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.7'
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages('python')
setup_keywords['package_dir'] = {'': 'python'}
setup_keywords['package_data'] = {'':['templates/*.glb', 'models/*.yml']}
setup_keywords['cmdclass'] = {'version': SetVersion, 'sdist': DistutilsSdist}
setup_keywords['test_suite']='snewpy.test.snewpy_test_suite.snewpy_test_suite'

requires = []
with open('requirements.txt', 'r') as f:
    for line in f:
        if line.strip():
            requires.append(line.strip())
setup_keywords['install_requires'] = requires
setup_keywords['extras_require'] = {  # Optional
    'dev': ['pytest', 'pytest-benchmark'],
    'docs':['numpydoc']
}
#
# Internal data directories.
#
#setup_keywords['data_files'] = [('snewpy/data/config', glob('data/config/*')),
#                                ('snewpy/data/spectra', glob('data/spectra/*'))]
#
# Run setup command.
#
setup(**setup_keywords)
