
#include "adiabatic_basis.h"

// *********************************************************************

using std::complex;
using std::array;

// *********************************************************************

MATRIX<complex<double>,NF,NF> W(array<double,NY> Y)
       { MATRIX<complex<double>,NF,NF> w; 
         w[0][0]=exp(-I*M_2PI*Y[9]); w[1][1]=exp(-I*M_2PI*Y[10]); w[2][2]=exp(-I*M_2PI*Y[11]);
         return w;
        }






