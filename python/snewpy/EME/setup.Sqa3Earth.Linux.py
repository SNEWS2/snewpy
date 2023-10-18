#!/usr/bin/env python3

# to use type "python3 setup.Sqa3Earth.py install --install-lib=destination/directory/" as sudo on the command line

from setuptools import setup
from setuptools.extension import Extension

Sqa3Earth = Extension('Sqa3Earth',
                    define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
                    include_dirs = ['.', '/usr/local/lib/python3.9/site-packages/pybind11/include', '../SQA/src', '../SQA/src/mstl', '../SQA/src/mstl/math2', '../SQA/src/mstl/math2/algebra', '../SQA/src/mstl/math2/analysis', '../SQA/src/mstl/math2/spline', '../SQA/src/mstl/physics'],
                    libraries = ['stdc++', 'm', 'gomp', 'python3'],
                    library_dirs = ['/usr/lib64'],
                    extra_compile_args = ['-std=c++17', '-fopenmp', '-fPIC', '-nostartfiles'],
                    extra_link_args = ['-shared'],
                    sources = ['../SQA/src/Sqa3Earth.cpp', '../SQA/src/adiabatic_basis.cpp', '../SQA/src/eigenvalues.cpp', '../SQA/src/flavour_basis.cpp', '../SQA/src/input_class.Sqa3Earth.cpp', '../SQA/src/jacobians.cpp', '../SQA/src/mixing_angles.cpp', '../SQA/src/output.Sqa3Earth.cpp', '../SQA/src/output_matrix.Sqa3Earth.cpp', '../SQA/src/parameters.cpp', '../SQA/src/potentials.cpp', '../SQA/src/RK.Sqa3Earth.cpp', '../SQA/src/update.Sqa3Earth.cpp', '../SQA/src/mstl/errors2.cpp', '../SQA/src/mstl/messages.cpp', '../SQA/src/mstl/miscellaneous functions.cpp', '../SQA/src/mstl/stdarg2.cpp', '../SQA/src/mstl/math2/algebra/column and row vectors.cpp', '../SQA/src/mstl/math2/algebra/linear algebra.cpp', '../SQA/src/mstl/math2/algebra/mmatrix.cpp', '../SQA/src/mstl/math2/analysis/algorithm3.cpp', '../SQA/src/mstl/math2/analysis/complex2.cpp', '../SQA/src/mstl/math2/analysis/derivative.cpp', '../SQA/src/mstl/math2/analysis/polynomial.cpp', '../SQA/src/mstl/math2/analysis/roots.cpp', '../SQA/src/mstl/math2/analysis/runge kutta.cpp', '../SQA/src/mstl/math2/analysis/special functions.cpp', '../SQA/src/mstl/math2/spline/discontinuous.cpp', '../SQA/src/mstl/math2/spline/interpolation data.cpp', '../SQA/src/mstl/physics/units and constants.cpp'])

setup (name = 'Sqa3Earth',
       version = '1.0',
       description = 'Sqa3 for Earth matter effects as a python module',
       author = 'Jim Kneller and Anne Graf',
       author_email = 'jpknelle@ncsu.edu',
       ext_modules = [Sqa3Earth],
      )


