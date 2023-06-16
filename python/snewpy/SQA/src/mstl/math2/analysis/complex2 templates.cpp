#if !defined(COMPLEX2_TEMPLATES)
#define COMPLEX2_TEMPLATES

#include "mstl.h"

template<typename type> std::complex<type> cbrt(std::complex<type> z)
       { return pow(z,1./3.);}

#endif

