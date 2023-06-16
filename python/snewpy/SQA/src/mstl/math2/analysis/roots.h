
#include<cstdarg>
#include<complex>
#include<vector>
#include<cmath>

#include "mstl.h"

#if !defined(_ROOTS)
#define _ROOTS

// ********************************************************************************************************************************
// ********************************************************************************************************************************
// ********************************************************************************************************************************

// Brent's algorithm for finding roots of the equation f(x)=C
template <class Functor>
double Root(double xlower,double xupper,const Functor &F,double C,double tolerance);

// ********************************************************************************************************************************
// ********************************************************************************************************************************
// ********************************************************************************************************************************

// Newton's method
// two seperate functions
//std::vector<double> NewtonRoots(double x01,double x02,double(*f1)(std::vector<double>),double(*f2)(std::vector<double>),double accuracy);

std::vector<double> NewtonRoots(std::vector<double> x0,std::vector<double(*)(std::vector<double>)> functions,std::vector<std::vector<double(*)(std::vector<double>)> > derivatives,double accuracy);

std::vector<double> NewtonRoots(std::vector<double> x0,std::vector<double>(*)(std::vector<double>),std::vector<std::vector<double> >(*)(std::vector<double>),double accuracy);

// ********************************************************************************************************************************
// ********************************************************************************************************************************
// ********************************************************************************************************************************

// Broyden's method
// two seperate functions
std::vector<double> BroydenRoots(double x01,double x02,double(*function1)(std::vector<double>),double(*function2)(std::vector<double>),double accuracy);

std::vector<double> BroydenRoots(std::vector<double> x0,std::vector<double(*)(std::vector<double>)> functions,double accuracy);

std::vector<double> BroydenRoots(std::vector<double> x0,std::vector<double(*)(std::vector<double>)> functions,std::vector<std::vector<double(*)(std::vector<double>)> > Jacobian,double accuracy);

std::vector<double> BroydenRoots(std::vector<double> x0,std::vector<double>(*)(std::vector<double>),std::vector<std::vector<double> >(*)(std::vector<double>),double accuracy);

// *************************************************************************
// *************************************************************************
// *************************************************************************

#endif
