
#include <cstdarg>
#include <cstdlib>
#include <cstddef>
#include <climits>
#include <ctime>
#include <cmath>
#include <complex>
#include <type_traits>
#include <algorithm>
#include <limits>       
//#include <values.h>

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


template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(Type A,Type B,unsigned int N=4);

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(Type A,std::complex<Type> B,unsigned int N=4);

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(std::complex<Type> A,Type B,unsigned int N=4);

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(std::complex<Type> A,std::complex<Type> B,unsigned int N=4);

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

#endif




