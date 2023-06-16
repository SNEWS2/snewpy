
#include<cstdarg>
#include<vector>
#include<iostream>
#include<cmath>

#include "mstl.h"

#if !defined(_ALGORITHM3)
#define _ALGORITHM3

// Line Search Algorithms
double LineSearch(std::vector<double> x0,std::vector<double> deltax,std::vector<double(*)(std::vector<double>)> functions);

double LineSearch(std::vector<double> x0,std::vector<double> deltax,std::vector<double>(*functions)(std::vector<double>));

#endif

