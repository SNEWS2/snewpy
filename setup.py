#!/usr/bin/env python
#
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
#
# Standard imports
#
from glob import glob
import os
import re
import sys
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
#
# Git-based version info. Remove?
#
from python.snewpy._git import get_version, SetVersion
#
# Begin setup
#
setup_keywords = dict()
#
setup_keywords['name'] = 'snewpy'
setup_keywords['description'] = 'SNEWS2 supernova simulation package'
setup_keywords['author'] = 'SNEWS2 Collaboration'
setup_keywords['author_email'] = 'snews2.0@lists.bnl.gov'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/SNEWS2/supernova_models'
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
if os.path.isdir('bin'):
    # Treat everything in bin as a script to be installed.
    setup_keywords['scripts'] = \
    [fname for fname in glob(os.path.join('bin', '*'))]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>3.3.0)']
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages('python')
setup_keywords['package_dir'] = {'': 'python'}
setup_keywords['cmdclass'] = {'version': SetVersion, 'sdist': DistutilsSdist}
setup_keywords['test_suite']='tests'
setup_keywords['tests_require']=['pytest']

requires = []
with open('requirements.txt', 'r') as f:
    for line in f:
        if line.strip():
            requires.append(line.strip())
setup_keywords['install_requires'] = requires
#
# Internal data directories.
#
#setup_keywords['data_files'] = [('snewpy/data/config', glob('data/config/*')),
#                                ('snewpy/data/spectra', glob('data/spectra/*'))]
#
# Run setup command.
#
setup(**setup_keywords)
