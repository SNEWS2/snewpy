
#include "output_matrix.h"

// *******************************************************************************

using std::complex;
using std::array;
using std::vector;
using std::cout;

// ********************************************************************************

void Pfm(double lambda,vector<vector<array<double,NY> > > &Y,vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > &C0,vector<vector<array<array<double,NF>,NF> > > A0,vector<vector<MATRIX<complex<double>,NF,NF> > > &Scumulative,vector<vector<vector<vector<double> > > > &PPfm)
      { double r, rrho, YYe;

        // ******

        r=sqrt(RE*RE+lambda*lambda+2.*RE*lambda*sin(altitude));
        rrho=rho(r);
        YYe=Ye(r); 

        array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;

        VfMSW[nu][e][e]=Ve(rrho,YYe);
        VfMSW[nu][mu][mu]=Vmu(rrho,YYe);
        VfMSW[nu][tau][tau]=Vtau(rrho,YYe);
	VfMSW[antinu]=-VfMSW[nu];

	vector<vector<MATRIX<complex<double>,NF,NF> > > Hf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
	vector<vector<MATRIX<complex<double>,NF,NF> > > UU(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));	
	vector<vector<MATRIX<complex<double>,NF,NF> > > Sa(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), Sfm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)); 

        vector<vector<array<double,NF> > > kk(NM,vector<array<double,NF> >(NE));
        vector<vector<array<double,NF> > > dkk(NM,vector<array<double,NF> >(NE));

        int i;
        #pragma omp parallel for schedule(static)
	for(i=0;i<=NE-1;i++)
           { Hf[nu][i]=HfV[nu][i] + VfMSW[nu];
             kk[nu][i]=k(Hf[nu][i]);
	     dkk[nu][i]=deltak(kk[nu][i]);
             UU[nu][i] = MixingMatrix(dkk[nu][i],C0[nu][i],A0[nu][i]);

	     Sa[nu][i] = W(Y[nu][i]) * B(Y[nu][i]);

             // take into account the density jump from Earth matter back to vacuum
	     Sm[nu][i] = Adjoint(UV[nu])*UU[nu][i] * Sa[nu][i] * Scumulative[nu][i];
             Sfm[nu][i]= UV[nu] * Sm[nu][i];

	     // *********
	     Hf[antinu][i]=HfV[antinu][i] + VfMSW[antinu];
	     kk[antinu][i]=kbar(Hf[antinu][i]);
	     dkk[antinu][i]=deltakbar(kk[antinu][i]);
	     UU[antinu][i]=MixingMatrix(dkk[antinu][i],C0[antinu][i],A0[antinu][i]);
       
	     Sa[antinu][i] = W(Y[antinu][i]) * B(Y[antinu][i]);

             // take into account the density jump from Earth matter back to vacuum
	     Sm[antinu][i] = Adjoint(UV[antinu])*UU[antinu][i] * Sa[antinu][i] * Scumulative[antinu][i];
             Sfm[antinu][i]= UV[antinu] * Sm[antinu][i];
	    }

        // *******

	for(i=0;i<=NE-1;i++)
           { PPfm[nu][i][e][0]=norm(Sfm[nu][i][e][0]); PPfm[nu][i][e][1]=norm(Sfm[nu][i][e][1]); PPfm[nu][i][e][2]=norm(Sfm[nu][i][e][2]);
	     PPfm[nu][i][mu][0]=norm(Sfm[nu][i][mu][0]); PPfm[nu][i][mu][1]=norm(Sfm[nu][i][mu][1]); PPfm[nu][i][mu][2]=norm(Sfm[nu][i][mu][2]);
	     PPfm[nu][i][tau][0]=norm(Sfm[nu][i][tau][0]); PPfm[nu][i][tau][1]=norm(Sfm[nu][i][tau][1]); PPfm[nu][i][tau][2]=norm(Sfm[nu][i][tau][2]);

	     PPfm[antinu][i][e][0]=norm(Sfm[antinu][i][e][0]); PPfm[antinu][i][e][1]=norm(Sfm[antinu][i][e][1]); PPfm[antinu][i][e][2]=norm(Sfm[antinu][i][e][2]);
	     PPfm[antinu][i][mu][0]=norm(Sfm[antinu][i][mu][0]); PPfm[antinu][i][mu][1]=norm(Sfm[antinu][i][mu][1]); PPfm[antinu][i][mu][2]=norm(Sfm[antinu][i][mu][2]);	     
             PPfm[antinu][i][tau][0]=norm(Sfm[antinu][i][tau][0]); PPfm[antinu][i][tau][1]=norm(Sfm[antinu][i][tau][1]); PPfm[antinu][i][tau][2]=norm(Sfm[antinu][i][tau][2]);
	    }

         //return PPfm;
        }
