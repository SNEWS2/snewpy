
#include "output_matrix.h"

// *******************************************************************************

using std::complex;
using std::array;
using std::vector;
using std::cout;
using std::norm;

// ********************************************************************************

void Pmf(double r,vector<vector<array<double,NY> > > &Y,vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > &C0,vector<vector<array<array<double,NF>,NF> > > A0,vector<vector<MATRIX<complex<double>,NF,NF> > > &Scumulative,vector<vector<vector<vector<double> > > > &PPmf)
      { double rrho, YYe;

        // ******

        rrho=exp(lnrho(log(r)));
        YYe=Ye(r); 

        array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;

        VfMSW[nu][e][e]=Ve(rrho,YYe);
        VfMSW[nu][mu][mu]=Vmu(rrho,YYe);
        VfMSW[nu][tau][tau]=Vtau(rrho,YYe);
	VfMSW[antinu]=-VfMSW[nu];

	vector<vector<MATRIX<complex<double>,NF,NF> > > Hf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
	vector<vector<MATRIX<complex<double>,NF,NF> > > UU(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));	
	vector<vector<MATRIX<complex<double>,NF,NF> > > Sa(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Smf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)); 

        vector<vector<array<double,NF> > > kk(NM,vector<array<double,NF> >(NE));
        vector<vector<array<double,NF> > > dkk(NM,vector<array<double,NF> >(NE));

        int i;
        #pragma omp parallel for schedule(static)
	for(i=0;i<=NE-1;i++)
           { Hf[nu][i]=HfV[nu][i] + VfMSW[nu];
             kk[nu][i]=k(Hf[nu][i]);
	     dkk[nu][i]=deltak(kk[nu][i]);
             UU[nu][i] = U(dkk[nu][i],C0[nu][i],A0[nu][i]);

	     Sa[nu][i] = W(Y[nu][i]) * B(Y[nu][i]);

	     Sm[nu][i] = Sa[nu][i] * Scumulative[nu][i];
             Smf[nu][i]= Sm[nu][i] * Adjoint(U0[nu][i]);

	     // *********
	     Hf[antinu][i]=HfV[antinu][i] + VfMSW[antinu];
	     kk[antinu][i]=kbar(Hf[antinu][i]);
	     dkk[antinu][i]=deltakbar(kk[antinu][i]);
	     UU[antinu][i]=MixingMatrix(dkk[antinu][i],C0[antinu][i],A0[antinu][i]);
       
	     Sa[antinu][i] = W(Y[antinu][i]) * B(Y[antinu][i]);

	     Sm[antinu][i] = Sa[antinu][i] * Scumulative[antinu][i];
             Smf[antinu][i]= Sm[antinu][i] * Adjoint(U0[antinu][i]);
	    }

        // *******

	for(i=0;i<=NE-1;i++)
           { PPmf[nu][i][e][0]=norm(Smf[nu][i][e][0]); PPmf[nu][i][e][1]=norm(Smf[nu][i][e][1]); PPmf[nu][i][e][2]=norm(Smf[nu][i][e][2]);
	     PPmf[nu][i][mu][0]=norm(Smf[nu][i][mu][0]); PPmf[nu][i][mu][1]=norm(Smf[nu][i][mu][1]); PPmf[nu][i][mu][2]=norm(Smf[nu][i][mu][2]);
	     PPmf[nu][i][tau][0]=norm(Smf[nu][i][tau][0]); PPmf[nu][i][tau][1]=norm(Smf[nu][i][tau][1]); PPmf[nu][i][tau][2]=norm(Smf[nu][i][tau][2]);

	     PPmf[antinu][i][e][0]=norm(Smf[antinu][i][e][0]); PPmf[antinu][i][e][1]=norm(Smf[antinu][i][e][1]); PPmf[antinu][i][e][2]=norm(Smf[antinu][i][e][2]);
	     PPmf[antinu][i][mu][0]=norm(Smf[antinu][i][mu][0]); PPmf[antinu][i][mu][1]=norm(Smf[antinu][i][mu][1]); PPmf[antinu][i][mu][2]=norm(Smf[antinu][i][mu][2]);	     
             PPmf[antinu][i][tau][0]=norm(Smf[antinu][i][tau][0]); PPmf[antinu][i][tau][1]=norm(Smf[antinu][i][tau][1]); PPmf[antinu][i][tau][2]=norm(Smf[antinu][i][tau][2]);
	    }

         //return PPmf;
        }
