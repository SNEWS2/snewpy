#!/usr/bin/env python3

# to use type "python3 setup.Sqa3Earth.py install --install-lib=destination/directory/" as sudo on the command line

from setuptools import setup
from setuptools.extension import Extension

Sqa3Earth = Extension('Sqa3Earth',
                    define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
                    include_dirs = ['.', './src', '/usr/local/lib/python3.9/site-packages/pybind11/include', './src/mstl_lite', './src/mstl_lite/math2', './src/mstl_lite/math2/algebra', './src/mstl_lite/math2/analysis', './src/mstl_lite/math2/spline', './src/mstl_lite/physics'],
                    libraries = ['stdc++', 'm', 'gomp', 'python3'],
                    library_dirs = ['/usr/lib64'],
                    extra_compile_args = ['-std=c++17', '-fopenmp', '-fPIC', '-nostartfiles'],
                    extra_link_args = ['-shared'],
                    sources = ['src/Sqa3Earth.cpp', 'src/adiabatic_basis.cpp', 'src/eigenvalues.cpp', 'src/flavour_basis.cpp', 'src/input_class.Sqa3Earth.cpp', 'src/jacobians.cpp', 'src/mixing_angles.cpp', 'src/output.Sqa3Earth.cpp', 'src/output_matrix.Sqa3Earth.cpp', 'src/parameters.cpp', 'src/potentials.cpp', 'src/RK.Sqa3Earth.cpp', 'src/update.Sqa3Earth.cpp', 'src/mstl_lite/errors2.cpp', 'src/mstl_lite/messages.cpp', 'src/mstl_lite/miscellaneous functions.cpp', 'src/mstl_lite/stdarg2.cpp', 'src/mstl_lite/math2/algebra/column and row vectors.cpp', 'src/mstl_lite/math2/algebra/linear algebra.cpp', 'src/mstl_lite/math2/algebra/mmatrix.cpp', 'src/mstl_lite/math2/analysis/algorithm3.cpp', 'src/mstl_lite/math2/analysis/complex2.cpp', 'src/mstl_lite/math2/analysis/derivative.cpp', 'src/mstl_lite/math2/analysis/polynomial.cpp', 'src/mstl_lite/math2/analysis/roots.cpp', 'src/mstl_lite/math2/analysis/runge kutta.cpp', 'src/mstl_lite/math2/analysis/special functions.cpp', 'src/mstl_lite/math2/spline/discontinuous.cpp', 'src/mstl_lite/math2/spline/interpolation base.cpp', 'src/mstl_lite/math2/spline/interpolation data.cpp', 'src/mstl_lite/physics/units and constants.cpp'])

setup (name = 'Sqa3Earth',
       version = '1.0',
       description = 'Sqa3 for Earth matter effects as a python module',
       author = 'Jim Kneller and Anne Graf',
       author_email = 'jpknelle@ncsu.edu',
       ext_modules = [Sqa3Earth],
      )


