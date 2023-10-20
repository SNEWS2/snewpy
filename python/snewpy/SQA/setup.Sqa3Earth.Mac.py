# to use type "sudo python3 setup.Sqa3Earth.py installon the command line

#!/usr/bin/env python
#
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
from setuptools.extension import Extension

#
# Begin setup
#
setup_keywords = dict()
#
setup_keywords['name'] = 'Sqa3Earth'
setup_keywords['description'] = 'A Python module for computing the Earth-matter effect upon neutrinos. Works alone or with snewpy.'
setup_keywords['author'] = 'Jim Kneller and Anne Graf'
setup_keywords['author_email'] = 'jpknelle@ncsu.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/SNEWS2/snewpy'
setup_keywords['version'] = '1.0'
#
# Use README.md as a long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.Sqa3Earth.md'):
    with open('README.Sqa3Earth.md') as readme:
        setup_keywords['long_description'] = readme.read()
    setup_keywords['long_description_content_type'] = 'text/markdown'
#
# Set other keywords for the setup function.
#
# Use entry_points to let `pip` create executable scripts for each target platform.
# See https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
# setup_keywords['entry_points'] = {'console_scripts': ['to_snowglobes = snewpy.to_snowglobes:generate_time_series', ], },
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.9'
setup_keywords['zip_safe'] = False
#setup_keywords['test_suite']='snewpy.test.snewpy_test_suite.snewpy_test_suite'

requires = []
with open('requirements.Sqa3Earth.txt', 'r') as f:
    for line in f:
        if line.strip():
            requires.append(line.strip())
setup_keywords['install_requires'] = requires
setup_keywords['extras_require'] = {  # Optional
    'dev': ['pytest'],
    'docs':['numpydoc']
}

Sqa3Earth = Extension('Sqa3Earth',
                    define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
                    include_dirs = ['/usr/local/lib/python3.9/site-packages/pybind11/include', './src', './src/mstl', './src/mstl/math2', './src/mstl/math2/algebra', './src/mstl/math2/analysis', './src/mstl/math2/spline', './src/mstl/physics'],
                    #libraries = ['stdc++', 'm', 'gomp', 'python3'],
                    library_dirs = ['/usr/lib64'],
                    extra_compile_args = ['-std=c++17', '-fopenmp', '-fPIC', '-nostartfiles'],
                    #extra_link_args = ['-shared'],
                    sources = ['./src/Sqa3Earth.cpp', './src/adiabatic_basis.cpp', './src/eigenvalues.cpp', './src/flavour_basis.cpp', './src/input_class.Sqa3Earth.cpp', './src/jacobians.cpp', './src/mixing_angles.cpp', './src/output.Sqa3Earth.cpp', './src/output_matrix.Sqa3Earth.cpp', './src/parameters.cpp', './src/potentials.cpp', './src/RK.Sqa3Earth.cpp', './src/update.Sqa3Earth.cpp', './src/mstl/errors2.cpp', './src/mstl/messages.cpp', './src/mstl/miscellaneous functions.cpp', './src/mstl/stdarg2.cpp', './src/mstl/math2/algebra/column and row vectors.cpp', './src/mstl/math2/algebra/linear algebra.cpp', './src/mstl/math2/algebra/mmatrix.cpp', './src/mstl/math2/analysis/algorithm3.cpp', './src/mstl/math2/analysis/complex2.cpp', './src/mstl/math2/analysis/derivative.cpp', './src/mstl/math2/analysis/polynomial.cpp', './src/mstl/math2/analysis/roots.cpp', './src/mstl/math2/analysis/runge kutta.cpp', './src/mstl/math2/analysis/special functions.cpp', './src/mstl/math2/spline/discontinuous.cpp', './src/mstl/math2/spline/interpolation data.cpp', './src/mstl/physics/units and constants.cpp'])

setup_keywords['ext_modules'] = [Sqa3Earth]

#
# Run setup command.
#
setup(**setup_keywords)
