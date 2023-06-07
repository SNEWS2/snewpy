
#include <cstdarg>
#include <cstdlib>
#include <cstddef>
#include <climits>
#include <ctime>
#include <cmath>
#include <complex>
#include <values.h>

#include "mstl.h"

#if !defined(_MATH2)
#define _MATH2

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

// additional math functions
template <class Type> Type Theta(const Type &T); // Heaviside Step Function;

int Kronecker(int i,int j); // return 1 if i==j, 0 otherwise

double Heaviside(double x,double x0); 
double tophat(double x,double xmax); 

template <class Type> bool Equality(Type A,Type B,int N=4);
template <class Type> bool Equality(std::complex<Type> A,std::complex<Type> B,int N=4);

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

#endif




