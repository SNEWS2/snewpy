
#include <cmath>
#include <iostream>
#include <complex>

#include "mstl.h"

#if !defined(COMPLEX2)
#define COMPLEX2

std::complex<int> Zero(std::complex<int>);
std::complex<double> Zero(std::complex<double>);

std::complex<int> One(std::complex<int>);
std::complex<double> One(std::complex<double>);

std::complex<int> Two(std::complex<int>);
std::complex<double> Two(std::complex<double>);

extern const std::complex<double> I;

template<typename type> std::complex<type> cbrt(std::complex<type> z);

#endif 
