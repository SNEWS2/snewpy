
#include "update.h"

// *****************************************************************************

using std::complex;
using std::array;
using std::vector;

// ***************************** UpdateSm ***************************************

vector<vector<MATRIX<complex<double>,NF,NF> > > UpdateSm(double lambdaminus,double lambdaplus,vector<vector<array<double,NY> > > &Y,vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > &C0,vector<vector<array<array<double,NF>,NF> > > A0,vector<vector<MATRIX<complex<double>,NF,NF> > > &Smprior)
          { array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;
            MATRIX<complex<double>,NF,NF> Hf,Hfbar;
            MATRIX<complex<double>,NF,NF> UU,UUbar; 

            array<double,NF> kk,kkbar,dkk,dkkbar;

            array<MATRIX<complex<double>,NF,NF>,NF> CC;
            array<array<double,NF>,NF> AA;

            vector<vector<MATRIX<complex<double>,NF,NF> > > Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));

            double r, rrho, YYe;

            // multiply SMSW by the mixing matrix at rminus
            r=sqrt(RE*RE+lambdaminus*lambdaminus+2.*RE*lambdaminus*sin(altitude));
            rrho = rho(r);
            YYe = Ye(r); 
    
            VfMSW[nu][e][e] = Ve(rrho,YYe); 
            VfMSW[nu][mu][mu] = Vmu(rrho,YYe); 
            VfMSW[nu][tau][tau] = Vtau(rrho,YYe); 

            VfMSW[antinu] = -VfMSW[nu];

            int i;
	    #pragma omp parallel for schedule(static) private(Hf,Hfbar,kk,kkbar,dkk,dkkbar,UU,UUbar)
	    for(i=0;i<=NE-1;i++)
               { Hf=HfV[nu][i]+VfMSW[nu];
                 kk=k(Hf);
                 dkk=deltak(kk);
                 UU=MixingMatrix(dkk,C0[nu][i],A0[nu][i]);
            	 Sm[nu][i]=UU * W(Y[nu][i]) * B(Y[nu][i]) * Smprior[nu][i];

            	 Hfbar=HfV[antinu][i]+VfMSW[antinu];
            	 kkbar=kbar(Hfbar);
            	 dkkbar=deltakbar(kkbar);
                 UUbar=MixingMatrix(dkkbar,C0[antinu][i],A0[antinu][i]);
            	 Sm[antinu][i]=UUbar * W(Y[antinu][i]) * B(Y[antinu][i]) * Smprior[antinu][i];
		}

            // multiply SMSW by the adjoint of the mixing matrix at rplus
            r=sqrt(RE*RE+lambdaplus*lambdaplus+2.*RE*lambdaplus*sin(altitude));
            rrho = rho(r);
            YYe = Ye(r); 
    
            VfMSW[nu][e][e] = Ve(rrho,YYe); 
            VfMSW[nu][mu][mu] = Vmu(rrho,YYe); 
            VfMSW[nu][tau][tau] = Vtau(rrho,YYe); 

            VfMSW[antinu] = -VfMSW[nu];

	    #pragma omp parallel for schedule(static) private(Hf,Hfbar,kk,kkbar,dkk,dkkbar,UU,UUbar) firstprivate(AA,CC)
	    for(i=0;i<=NE-1;i++)
               { Hf=HfV[nu][i]+VfMSW[nu];
                 kk=k(Hf);
             	 dkk=deltak(kk);
                 CofactorMatrices(Hf,kk,CC);
                 AA=MixingMatrixFactors(CC,C0[nu][i],A0[nu][i]);
                 UU=MixingMatrix(dkk,CC,AA);  
            	 Sm[nu][i]=Adjoint(UU)*MATRIX<complex<double>,NF,NF>(Sm[nu][i]); 

            	 Hfbar=HfV[antinu][i]+VfMSW[antinu];
                 kkbar=kbar(Hfbar);
            	 dkkbar=deltakbar(kkbar);
                 CofactorMatrices(Hfbar,kkbar,CC);
                 AA=MixingMatrixFactors(CC,C0[antinu][i],A0[antinu][i]);
                 UUbar=MixingMatrix(dkkbar,CC,AA);
            	 Sm[antinu][i]=Adjoint(UUbar)*MATRIX<complex<double>,NF,NF>(Sm[antinu][i]);
	        }

            return Sm;
        }

// ********************************************************************************

vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > UpdateC(double lambda)
          { array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;
            MATRIX<complex<double>,NF,NF> Hf,Hfbar;
            array<double,NF> kk,kkbar;

            vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > CC(NM,vector<array<MATRIX<complex<double>,NF,NF>,NF> >(NE));

            double r=sqrt(RE*RE+lambda*lambda+2.*RE*lambda*sin(altitude));
            double rrho = rho(r);
            double YYe=Ye(r); 
    
            VfMSW[nu][e][e] = Ve(rrho,YYe); 
            VfMSW[nu][mu][mu] = Vmu(rrho,YYe); 
            VfMSW[nu][tau][tau] = Vtau(rrho,YYe); 

            VfMSW[antinu] = -VfMSW[nu];

            int i;
            #pragma omp parallel for schedule(static) private(Hf,Hfbar,kk,kkbar)
	    for(i=0;i<=NE-1;i++)
               { Hf=HfV[nu][i]+VfMSW[nu]; 
                 kk=k(Hf);
                 CofactorMatrices(Hf,kk,CC[nu][i]);

            	 Hfbar=HfV[antinu][i]+VfMSW[antinu];
            	 kkbar=kbar(Hfbar);
                 CofactorMatrices(Hfbar,kkbar,CC[antinu][i]);
	        }

            return CC;
        }

// ********************************************************************************

vector<vector<array<array<double,NF>,NF> > > UpdateA(vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > &C,vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > &C0,vector<vector<array<array<double,NF>,NF> > > A0)
      { vector<vector<array<array<double,NF>,NF> > > A(NM,vector<array<array<double,NF>,NF> >(NE));

        int i;
        #pragma omp parallel for schedule(static)
        for(i=0;i<=NE-1;i++)
           { A[nu][i]=MixingMatrixFactors(C[nu][i],C0[nu][i],A0[nu][i]);
             A[antinu][i]=MixingMatrixFactors(C[antinu][i],C0[antinu][i],A0[antinu][i]);
            }

        return A;
       }


