
#ifndef EMEWS_H
#define EMEWS_H

// *****************************************************************

#include <cmath>

#include <complex>
//using std::complex;
//using std::polar;
//using std::abs;
//using std::arg;
//using std::real;
//using std::imag;
//using std::norm;

#include <cstdarg>
//using std::va_list;

#include <cstdlib>

#include<iostream>
//using::std::cout;

#include <ostream>
//using std::ostream;
//using std::endl;
//using std::flush;

#include <fstream>
//using std::ifstream;
//using std::ofstream;

#include <sstream>
//using std::stringstream;

#include <algorithm>
//using std::min;
//using std::max;
//using std::sort;
//using std::swap;
//using std::lower_bound;
//using std::upper_bound;

#include <string>
//using std::string;

#include <utility>
//using std::pair;

#include <functional>

#include <limits>
//using std::numeric_limits;

#include <vector>
//using std::vector;

#include <array>
//using std::array;

#include <omp.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// ************************

#include "mstl.h"
//using namespace prefixes;
//using interpolation::DISCONTINUOUS;

// ************************

#include "parameters.h"
#include "potentials.h"

#include "eigenvalues.h"
#include "mixing_angles.h"

#include "adiabatic_basis.h"
#include "flavour_basis.h"
#include "jacobians.h"
#include "RK.h"
#include "update.h"

//#include "input.h"
#include "output.h"

#include "input_class.h"
#include "output_matrix.h"

// ********************************************************************** 

std::vector<std::vector<std::vector<std::vector<double> > > > Run(InputDataEMEWS ID);

// ************************ Neutrino Potentials **************************

// DISCONTINUOUS is a cubic spline interpolator based on Akima's algorithm but it can handle discontinuities 
extern interpolation::DISCONTINUOUS rho, lnrho, Ye, v, M;

// *****************************************************************

#endif
