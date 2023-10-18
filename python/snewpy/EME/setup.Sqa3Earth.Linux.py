# to use type "sudo python3 setup.Sqa3Earth.py install on the command line

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
                    include_dirs = ['.', './python/snewpy/SQA/src', '/usr/local/lib/python3.9/site-packages/pybind11/include', 
                                    './python/snewpy/SQA/src/mstl', './python/snewpy/SQA/src/mstl/math2', './python/snewpy/SQA/src/mstl/math2/algebra', 
                                    './python/snewpy/SQA/src/mstl/math2/analysis', './python/snewpy/SQA/src/mstl/math2/spline', './python/snewpy/SQA/src/mstl/physics'],
                    libraries = ['stdc++', 'm', 'gomp', 'python3'],
                    library_dirs = ['/usr/lib64'],
                    extra_compile_args = ['-std=c++17', '-fopenmp', '-fPIC', '-nostartfiles'],
                    extra_link_args = ['-shared'],
                    sources = ['python/snewpy/SQA/src/Sqa3Earth.cpp', 'python/snewpy/SQA/src/adiabatic_basis.cpp', 'python/snewpy/SQA/src/eigenvalues.cpp', 
                               'python/snewpy/SQA/src/flavour_basis.cpp', 'python/snewpy/SQA/src/input_class.Sqa3Earth.cpp', 'python/snewpy/SQA/src/jacobians.cpp', 
                               'python/snewpy/SQA/src/mixing_angles.cpp', 'python/snewpy/SQA/src/output.Sqa3Earth.cpp', 'python/snewpy/SQA/src/output_matrix.Sqa3Earth.cpp', 
                               'python/snewpy/SQA/src/parameters.cpp', 'python/snewpy/SQA/src/potentials.cpp', 'python/snewpy/SQA/src/RK.Sqa3Earth.cpp', 
                               'python/snewpy/SQA/src/update.Sqa3Earth.cpp', 'python/snewpy/SQA/src/mstl/errors2.cpp', 'python/snewpy/SQA/src/mstl/messages.cpp', 
                               'python/snewpy/SQA/src/mstl/miscellaneous functions.cpp', 'python/snewpy/SQA/src/mstl/stdarg2.cpp', 
                               'python/snewpy/SQA/src/mstl/math2/algebra/column and row vectors.cpp', 'python/snewpy/SQA/src/mstl/math2/algebra/linear algebra.cpp', 
                               'python/snewpy/SQA/src/mstl/math2/algebra/mmatrix.cpp', 'python/snewpy/SQA/src/mstl/math2/analysis/algorithm3.cpp', 
                               'python/snewpy/SQA/src/mstl/math2/analysis/complex2.cpp', 'python/snewpy/SQA/src/mstl/math2/analysis/derivative.cpp', 
                               'python/snewpy/SQA/src/mstl/math2/analysis/polynomial.cpp', 'python/snewpy/SQA/src/mstl/math2/analysis/roots.cpp', 
                               'python/snewpy/SQA/src/mstl/math2/analysis/runge kutta.cpp', 'python/snewpy/SQA/src/mstl/math2/analysis/special functions.cpp', 
                               'python/snewpy/SQA/src/mstl/math2/spline/discontinuous.cpp', 'python/snewpy/SQA/src/mstl/math2/spline/interpolation data.cpp', 
                               'python/snewpy/SQA/src/mstl/physics/units and constants.cpp'])

setup_keywords['ext_modules'] = [Sqa3Earth]

#
# Run setup command.
#
setup(**setup_keywords)
