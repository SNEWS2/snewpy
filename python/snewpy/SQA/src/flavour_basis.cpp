
#include "flavour_basis.h"

// *********************************************************************

using std::complex;

// *********************************************************************

void Evaluate_HfV(void)
      { MATRIX<complex<double>,NF,NF> KV;

	for(int i=0;i<=NE-1;i++)
           { for(int j=0;j<=NF-1;j++){ KV[j][j]=kV[i][j];}
	     HfV[nu][i]=UV[nu]*KV*Adjoint(UV[nu]);
	     HfV[antinu][i]=UV[antinu]*KV*Adjoint(UV[antinu]);
	    }

       }


