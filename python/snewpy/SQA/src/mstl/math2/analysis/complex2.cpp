#include "complex2.h"

using std::complex;
using std::ostream;

const complex<double> I(0.,1.);

complex<int> Zero(complex<int>){ return complex<int>(0,0);}
complex<double> Zero(complex<double>){ return complex<double>(0.,0.);}

complex<int> One(complex<int>){ return complex<int>(1,0);}
complex<double> One(complex<double>){ return complex<double>(1.,0.);}

complex<int> Two(complex<int>){ return complex<int>(2,0);}
complex<double> Two(complex<double>){ return complex<double>(2.,0.);}



