
#include "eigenvalues.h"

// ********************************************************************

using std::complex;

using std::max;
using std::sort;

using std::array;
using std::vector;

// *********************************************************************

double Q(MATRIX<complex<double>,NF,NF> Hf)
       { return -( norm(Hf[e][e]-Hf[mu][mu]) +norm(Hf[e][e]-Hf[tau][tau]) +norm(Hf[mu][mu]-Hf[tau][tau]) )/18.
                -( norm(Hf[e][mu]) +norm(Hf[e][tau]) +norm(Hf[mu][tau]) )/3.;
        }

double R(MATRIX<complex<double>,NF,NF> Hf)
       { return real( (2.*Hf[e][e]-Hf[mu][mu]-Hf[tau][tau])*(2.*Hf[mu][mu]-Hf[e][e]-Hf[tau][tau])*(2.*Hf[tau][tau]-Hf[e][e]-Hf[mu][mu]) )/54.
               +real( Hf[e][mu]*Hf[mu][tau]*Hf[tau][e] )
               -real( norm(Hf[mu][tau])*(2.*Hf[e][e]-Hf[mu][mu]-Hf[tau][tau]) +norm(Hf[e][tau])*(2.*Hf[mu][mu]-Hf[e][e]-Hf[tau][tau]) +norm(Hf[e][mu])*(2.*Hf[tau][tau]-Hf[e][e]-Hf[mu][mu]) )/6.;
        }

double omega(double Q,double sqrtminusQ,double R)
       { return acos(-R/Q/sqrtminusQ);}

double piminusomegabar(double Qbar,double sqrtminusQbar,double Rbar)
       { return acos(Rbar/Qbar/sqrtminusQbar);}

// *********************************************************************

double k1(double T,double sqrtminusQ,double comega_3,double somega_3)
       { return T/3. +2.*sqrtminusQ*(comega_3*comega1-somega_3*somega1);}

double k1bar(double Tbar,double sqrtminusQbar,double cpiminusomegabar_3,double spiminusomegabar_3)
       { return Tbar/3. +2.*sqrtminusQbar*(cpiminusomegabar_3*comega1p60+spiminusomegabar_3*somega1p60);}

// **************************

double k2(double T,double sqrtminusQ,double comega_3,double somega_3)
       { return T/3. +2.*sqrtminusQ*(comega_3*comega2-somega_3*somega2);}

double k2bar(double Tbar,double sqrtminusQbar,double cpiminusomegabar_3,double spiminusomegabar_3)
       { return Tbar/3. +2.*sqrtminusQbar*(cpiminusomegabar_3*comega2p60+spiminusomegabar_3*somega2p60);}

// **********************************

double k3(double T,double sqrtminusQ,double comega_3,double somega_3)
       { return T/3. +2.*sqrtminusQ*(comega_3*comega3-somega_3*somega3);}

double k3bar(double Tbar,double sqrtminusQbar,double cpiminusomegabar_3,double spiminusomegabar_3)
       { return Tbar/3. +2.*sqrtminusQbar*(cpiminusomegabar_3*comega3p60+spiminusomegabar_3*somega3p60);}

// ****************************

double deltak12(double sqrtminusQ,double comega_3,double somega_3)
       { return 4.*sqrtminusQ*(somega_3*comega12_2+comega_3*somega12_2)*somega12;} 

double deltak13(double sqrtminusQ,double comega_3,double somega_3)
       { return 4.*sqrtminusQ*(somega_3*comega13_2+comega_3*somega13_2)*somega13;}

double deltak23(double sqrtminusQ,double comega_3,double somega_3)
       { return 4.*sqrtminusQ*(somega_3*comega23_2+comega_3*somega23_2)*somega23;}

double deltak12bar(double sqrtminusQbar,double cpiminusomega_3,double spiminusomega_3)
       { return 4.*sqrtminusQbar*(-spiminusomega_3*comega12p120_2+cpiminusomega_3*somega12p120_2)*somega12;} 

double deltak13bar(double sqrtminusQbar,double cpiminusomega_3,double spiminusomega_3)
       { return 4.*sqrtminusQbar*(-spiminusomega_3*comega13p120_2+cpiminusomega_3*somega13p120_2)*somega13;}

double deltak23bar(double sqrtminusQbar,double cpiminusomega_3,double spiminusomega_3)
       { return 4.*sqrtminusQbar*(-spiminusomega_3*comega23p120_2+cpiminusomega_3*somega23p120_2)*somega23;}

// *********************************************************************
// *********************************************************************
// *********************************************************************

array<double,NF> k(MATRIX<complex<double>,NF,NF> Hf)
      { array<double,NF> k;

        double t=real(Trace(Hf)), q=Q(Hf), sqrtminusq=sqrt(-q), r=R(Hf);
        double Omega, cOmega=-r/q/sqrtminusq; 

        Omega=acos(cOmega);
        double cOmega_3=cos(Omega/3.), sOmega_3=sin(Omega/3.);

        k[0]=k1(t,sqrtminusq,cOmega_3,sOmega_3);
        k[1]=k2(t,sqrtminusq,cOmega_3,sOmega_3);
        k[2]=k3(t,sqrtminusq,cOmega_3,sOmega_3);

        return k;
       }

array<double,NF> kbar(MATRIX<complex<double>,NF,NF> Hfbar)
     { array<double,NF> k;

       double t=real(Trace(Hfbar)), q=Q(Hfbar), sqrtminusq=sqrt(-q), r=R(Hfbar);
       double PiminusOmegabar, cPiminusOmegabar=r/q/sqrtminusq;

       PiminusOmegabar=acos(cPiminusOmegabar);
       double cPiminusOmegabar_3=cos(PiminusOmegabar/3.), sPiminusOmegabar_3=sin(PiminusOmegabar/3.);

       k[0]=k1bar(t,sqrtminusq,cPiminusOmegabar_3,sPiminusOmegabar_3);
       k[1]=k2bar(t,sqrtminusq,cPiminusOmegabar_3,sPiminusOmegabar_3);
       k[2]=k3bar(t,sqrtminusq,cPiminusOmegabar_3,sPiminusOmegabar_3);

       return k;
      }

// *********************************************************************

array<double,NF> deltak(MATRIX<complex<double>,NF,NF> Hf)
      { array<double,NF> dk;

        double q=Q(Hf), sqrtminusq=sqrt(-q), Omega, cOmega=-R(Hf)/q/sqrtminusq;

        Omega=acos(cOmega);
        double cOmega_3=cos(Omega/3.), sOmega_3=sin(Omega/3.);

        dk[0]=deltak12(sqrtminusq,cOmega_3,sOmega_3);
        dk[1]=deltak13(sqrtminusq,cOmega_3,sOmega_3);
        dk[2]=deltak23(sqrtminusq,cOmega_3,sOmega_3);

        return dk;
       }

array<double,NF> deltakbar(MATRIX<complex<double>,NF,NF> Hfbar)
      { array<double,NF> dk;

        double qbar=Q(Hfbar), sqrtminusqbar=sqrt(-qbar), PiminusOmegabar, cPiminusOmegabar=R(Hfbar)/qbar/sqrtminusqbar;

        PiminusOmegabar=acos(cPiminusOmegabar);
        double cPiminusOmegabar_3=cos(PiminusOmegabar/3.), sPiminusOmegabar_3=sin(PiminusOmegabar/3.);

        dk[0]=deltak12bar(sqrtminusqbar,cPiminusOmegabar_3,sPiminusOmegabar_3);
        dk[1]=deltak13bar(sqrtminusqbar,cPiminusOmegabar_3,sPiminusOmegabar_3);
        dk[2]=deltak23bar(sqrtminusqbar,cPiminusOmegabar_3,sPiminusOmegabar_3);

        return dk;
       }

array<double,NF> deltak(array<double,NF> k)
      { array<double,NF> dk;
        dk[0]=k[0]-k[1];
        dk[1]=k[0]-k[2];
        dk[2]=k[1]-k[2];
        return dk;
       }

array<double,NF> deltakbar(array<double,NF> kbar)
      { return deltak(kbar);}

// *********************************************************************
// *********************************************************************
// *********************************************************************

double Evaluate_C(MATRIX<complex<double>,NF,NF> Hf,double comegai,double somegai)
       { return ( (1.+2.*comegai)*real(Hf[e][e]) + (1.-comegai)*real(Hf[mu][mu]+Hf[tau][tau]) )/3. - sqrt(4.*norm(Hf[mu][tau])+norm(Hf[mu][mu]-Hf[tau][tau]))/M_SQRT3*somegai;}

double Evaluate_Cbar(MATRIX<complex<double>,NF,NF> Hf,double comegaip60,double somegaip60)
       { return ( (1.-2.*comegaip60)*real(Hf[e][e]) + (1.+comegaip60)*real(Hf[mu][mu]+Hf[tau][tau]) )/3. + sqrt(4.*norm(Hf[mu][tau])+norm(Hf[mu][mu]-Hf[tau][tau]))/M_SQRT3*somegaip60;}

double Evaluate_a(MATRIX<complex<double>,NF,NF> Hf,double comegai,double somegai)
       { return ( 2.*real(Hf[e][mu]*Hf[mu][tau]*Hf[tau][e]+Hf[e][tau]*Hf[tau][mu]*Hf[mu][e]) + (norm(Hf[e][mu])-norm(Hf[e][tau]))*real(Hf[mu][mu]-Hf[tau][tau]) )/M_SQRT3/sqrt(4.*norm(Hf[mu][tau])+norm(Hf[mu][mu]-Hf[tau][tau])) * somegai + (norm(Hf[e][mu])+norm(Hf[e][tau])) * comegai;}

// ***********

void Evaluate_omega1(void)
     { if(ordering[0]==0){ omega1=2.*M_PI/3.;}
       if(ordering[1]==0){ omega1=4.*M_PI/3.;}
       if(ordering[2]==0){ omega1=0.;}  
      }

void Evaluate_omega2(void)
     { if(ordering[0]==1){ omega2=2.*M_PI/3.;}
       if(ordering[1]==1){ omega2=4.*M_PI/3.;}
       if(ordering[2]==1){ omega2=0.;}  
      }

void Evaluate_omega3(void)
     { if(ordering[0]==2){ omega3=2.*M_PI/3.;}
       if(ordering[1]==2){ omega3=4.*M_PI/3.;}
       if(ordering[2]==2){ omega3=0.;}  
      }

void Evaluate_comega1(void)
     { if(ordering[0]==0){ comega1=-0.5;}
       if(ordering[1]==0){ comega1=-0.5;}
       if(ordering[2]==0){ comega1=1.;}
      }

void Evaluate_comega2(void)
     { if(ordering[0]==1){ comega2=-0.5;}
       if(ordering[1]==1){ comega2=-0.5;}
       if(ordering[2]==1){ comega2=1.;}
      }

void Evaluate_comega3(void)
     { if(ordering[0]==2){ comega3=-0.5;}
       if(ordering[1]==2){ comega3=-0.5;}
       if(ordering[2]==2){ comega3=1.;}
      }

void Evaluate_somega1(void)
     { if(ordering[0]==0){ somega1=M_SQRT3/2.;}
       if(ordering[1]==0){ somega1=-M_SQRT3/2.;}
       if(ordering[2]==0){ somega1=0.;}
      }

void Evaluate_somega2(void)
     { if(ordering[0]==1){ somega2=M_SQRT3/2.;}
       if(ordering[1]==1){ somega2=-M_SQRT3/2.;}
       if(ordering[2]==1){ somega2=0.;}
      }

void Evaluate_somega3(void)
     { if(ordering[0]==2){ somega3=M_SQRT3/2.;}
       if(ordering[1]==2){ somega3=-M_SQRT3/2.;}
       if(ordering[2]==2){ somega3=0.;}
      }

// ***********

void Evaluate_comega1p60(void)
     { if(ordering[0]==0){ comega1p60=-1.;}
       if(ordering[1]==0){ comega1p60=0.5;}
       if(ordering[2]==0){ comega1p60=0.5;}
      }

void Evaluate_comega2p60(void)
     { if(ordering[0]==1){ comega2p60=-1.;}
       if(ordering[1]==1){ comega2p60=0.5;}
       if(ordering[2]==1){ comega2p60=0.5;}
      }

void Evaluate_comega3p60(void)
     { if(ordering[0]==2){ comega3p60=-1.;}
       if(ordering[1]==2){ comega3p60=0.5;}
       if(ordering[2]==2){ comega3p60=0.5;}
      }

void Evaluate_somega1p60(void)
     { if(ordering[0]==0){ somega1p60=0.;}
       if(ordering[1]==0){ somega1p60=-M_SQRT3/2.;}
       if(ordering[2]==0){ somega1p60=M_SQRT3/2.;}
      }

void Evaluate_somega2p60(void)
     { if(ordering[0]==1){ somega2p60=0.;}
       if(ordering[1]==1){ somega2p60=-M_SQRT3/2.;}
       if(ordering[2]==1){ somega2p60=M_SQRT3/2.;}
      }

void Evaluate_somega3p60(void)
     { if(ordering[0]==2){ somega3p60=0.;}
       if(ordering[1]==2){ somega3p60=-M_SQRT3/2.;}
       if(ordering[2]==2){ somega3p60=M_SQRT3/2.;}
      }

void Evaluate_somega12(void){ somega12=-sin((omega1-omega2)/2.);}

void Evaluate_somega13(void){ somega13=-sin((omega1-omega3)/2.);}

void Evaluate_somega23(void){ somega23=-sin((omega2-omega3)/2.);}

void Evaluate_comega12_2(void){ comega12_2=cos((omega1+omega2)/2.);}

void Evaluate_comega13_2(void){ comega13_2=cos((omega1+omega3)/2.);}

void Evaluate_comega23_2(void){ comega23_2=cos((omega2+omega3)/2.);}

void Evaluate_comega12p120_2(void){ comega12p120_2=cos((omega1+omega2+2.*M_PI/3.)/2.);}

void Evaluate_comega13p120_2(void){ comega13p120_2=cos((omega1+omega3+2.*M_PI/3.)/2.);}

void Evaluate_comega23p120_2(void){ comega23p120_2=cos((omega2+omega3+2.*M_PI/3.)/2.);}

void Evaluate_somega12_2(void){ somega12_2=sin((omega1+omega2)/2.);}

void Evaluate_somega13_2(void){ somega13_2=sin((omega1+omega3)/2.);}

void Evaluate_somega23_2(void){ somega23_2=sin((omega2+omega3)/2.);}

void Evaluate_somega12p120_2(void){ somega12p120_2=sin((omega1+omega2+2.*M_PI/3.)/2.);}

void Evaluate_somega13p120_2(void){ somega13p120_2=sin((omega1+omega3+2.*M_PI/3.)/2.);}

void Evaluate_somega23p120_2(void){ somega23p120_2=sin((omega2+omega3+2.*M_PI/3.)/2.);}

