
// the derivative formulae are taken from James, p397

#include <cstdarg>
#include <vector>

#include "mstl.h"

#if !defined(_DERIVATIVE)
#define _DERIVATIVE

// ******************************************************************************
// ******************************************************************************
// ******************** Finite Difference Formulae ******************************
// ******************************************************************************
// ******************************************************************************

// the std::vector this returns comprises the estimate of the function at x0 and all
// the computable derivatives given the value of y at the points x
std::vector<double> FiniteDifference1D(double x0,std::vector<double> x,const std::vector<double> &y);

std::vector<std::vector<double> > FiniteDifference2D(double x0,double y0,std::vector<double> x,std::vector<double> y,const std::vector<std::vector<double> > &z);

// *************************************************************************
// *************************************************************************
// *************************************************************************

#endif

