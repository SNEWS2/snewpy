#include "special functions.h"

using std::numeric_limits;
using std::complex;
using std::real;
using std::imag;
using std::vector;

// *******************************************************
// *******************************************************
// *******************************************************

// the truncated exponential function, see A&S page 262
double exp(int n,double X)
       { if(n<0){ throw NEGATIVE_NUMBER("exp(int,double");}
         if( n==0){ return 1.;}
         double e=0.;
         if(X>1.){ for(int i=0;i<=n;i++){ e+=pow(X,i)/Factorial(i);} }
         else{ for(int i=n;i>=0;i--){ e+=pow(X,i)/Factorial(i);} }
         return e;
        }

// *******************************************************

// see Arfken p558 on Numerical Computation of the Gamma function and G&R p941
double Gamma(double N){ return tgamma(N);}
     /*{ if(N<0. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("Gamma(double)");}
         if(Equality(N,0.5)==true){ return M_SQRTPI;}
         if(N<1.){ return Gamma(N+1.)/N;}

         if(N>2.){ return (N-1.)*Gamma(N-1.);}
         if(N>=1.){ return 1. +(N-1.)*( -0.577191652 +(N-1.)*( 0.988205891 
                              +(N-1.)*( -0.897056937 +(N-1.)*( 0.918206857 
                              +(N-1.)*( -0.756704078 +(N-1.)*( 0.482199394 
                              +(N-1.)*( -0.193527818 +(N-1.)*0.035868343 )))))));         
                   }
        }*/

// see Numerical Recipes for the algorithm
double logGamma(double N){ return gamma(N);}
       /*{ if(N<0.){ throw NEGATIVE_NUMBER("logGamma(double)");}
         double ser; double n=N;
         double c[6]={ 76.18009172947146,-86.50532032941677, 24.01409824083091,
                      -1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5 };
         ser=1.000000000190015;
         for(int i=0;i<=5;i++){ ser+=c[i]/(++n);}
         return -(N+5.5) + (N+0.5)*log(N+5.5) + log(M_SQRT2*M_SQRTPI*ser/N);
        }*/

complex<double> Gamma(complex<double> Z)
       { if(real(Z)<0. && Equality(fmod(real(Z),1.),0.)==true){ throw NEGATIVE_INTEGER("Gamma(complex<double>)");}
         return exp(logGamma(Z));
        }

// see efg2.com, this is just Lanczos approximation
complex<double> logGamma(complex<double> Z)
       { if(real(Z)<0. && Equality(fmod(real(Z),1.),0.)==true){ throw NEGATIVE_INTEGER("logGamma(complex<double>)");}
         if(imag(Z)<0.){ return conj(logGamma(conj(Z)));}
         if(real(Z)<9.){ return logGamma(Z+1.)-log(Z);}
 
         double C[]={1./12.,-1./360.,1./1260.,-1./1680.,1./1188.,-691./360360.,1./156.,-3617./122400.};
         complex<double> Z2(Z*Z);  
         return (Z-0.5)*log(Z) -Z +log(M_2PI)/2.
               +( C[0] +(C[1] +(C[2] +(C[3] +(C[4] +(C[5] +(C[6] +C[7]/Z2)/Z2)/Z2)/Z2)/Z2)/Z2)/Z2)/Z;
         
        }

// *******************************************************

// Gamma(N,X) see Numerical Recipes for the algorithm
double Gamma(double N,double X)
       { if(N<1. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("Gamma(double,double)");}
         if(X<0.){ NEGATIVE_NUMBER("Gamma(double,double)");}
         if(Equality(X,0.)==true){ return Gamma(N);}
         if(N>2.){ return (N-1.)*Gamma(N-1.,X)+exp(-X+(N-1.)*log(X));}
         if(N<1.){ return ( Gamma(N+1.,X)-exp(-X+N*log(X)) )/N;}

         if(X<N+1.){ return Gamma(N)-IncompleteGamma(N,X);}

         double a,b,j,C,D,F,T=1e-30;

         j=1.; a=1.; b=X-N+1.; C=a/numeric_limits<double>::epsilon(); D=1./b; F=D;

         do{ a=j*(N-j);
             b+=2.;
             C=b+a/C; if(Equality(fabs(C),0.)==true){ C=numeric_limits<double>::epsilon();}
             D=b+a*D; if(Equality(fabs(D),0.)==true){ D=numeric_limits<double>::epsilon();}
             D=1./D;
             F*=C*D;
             j++;
            }
         while(1.-fabs(C*D)>T);

         return exp(-X+N*log(X))*F;
        }

double logGamma(double N,double X)
       { if(Equality(X,0.)==true){ return logGamma(N);}
         return log(Gamma(N,X));
        }

// Gamma(N,X1)-Gamma(N,X2)
double Gamma(double N,double X1,double X2)
       { if(N<1. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("Gamma(double,double,double)");}
         if(X1<0. || X2<0.){ NEGATIVE_NUMBER("Gamma(double,double,double)");}
         if(Equality(X1,X2)==true){ return 0.;}

         return -IncompleteGamma(N,X1,X2);
        }

complex<double> Gamma(double N,complex<double> Z)
       { if(N<1. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("Gamma(double,complex<double>)");}
         if(real(Z)<0.){ NEGATIVE_NUMBER("Gamma(double,complex<double>)");}
         if(Equality(real(Z),0.)==true && Equality(imag(Z),0.)==true){ return complex<double>(Gamma(N));}
         if(N>2.){ return (N-1.)*Gamma(N-1.,Z)+exp(-Z+(N-1.)*log(Z));}
         if(N<1.){ return ( Gamma(N+1.,Z)-exp(-Z+N*log(Z)) )/N;}

         if(real(Z)<N+1.){ return Gamma(N)-IncompleteGamma(N,Z);}

         double j,T=1e-30;
         complex<double> a,b,C,D,F;

         j=1.; a=1.; b=Z-N+1.; C=a/numeric_limits<double>::epsilon(); D=1./b; F=D;

         do{ a=j*(N-j);
             b+=2.;
             C=b+a/C; if(Equality(abs(C),0.)==true){ C=numeric_limits<double>::epsilon();}
             D=b+a*D; if(Equality(abs(D),0.)==true){ D=numeric_limits<double>::epsilon();}
             D=1./D;
             F*=C*D;
             j++;
            }
         while(1.-abs(C*D)>T);// && arg(C*D)>T);

         return exp(-Z+N*log(Z))*F;
        }

// *******************************************************

double IncompleteGamma(double N,double X)
       { if(N<1. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("IncompleteGamma(double,double)");}
         if(X<0.){ NEGATIVE_NUMBER("IncompleteGamma(double,double)");}
         if(Equality(X,0.)==true){ return 0.;}
         if(N>2.){ return (N-1.)*IncompleteGamma(N-1.,X)-exp(-X+(N-1.)*log(X));}
         if(N<1.){ return ( IncompleteGamma(N+1.,X)+exp(-X+N*log(X)) )/N;}

         if(X>N+1.){ return Gamma(N)-Gamma(N,X);}

         double dG=1./N, G=0., A=1.;
         do{ G+=dG; dG*=X/(N+A); A++;}
         while(fabs(dG)>=numeric_limits<double>::epsilon()*fabs(G));
         return exp(-X+N*log(X))*G;
        }

double logIncompleteGamma(double N,double X)
       { if(Equality(X,0.)==true){ throw ZERO_NUMBER("logIncompleteGamma(double,double)");}
         return log(IncompleteGamma(N,X));
        }

double IncompleteGamma(double N,double X1,double X2) //=gamma(N,X2) - gamma(N,X1)
       { if(N<1. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("IncompleteGamma(double,double,double)");}
         if(X1<0. || X2<0.){ NEGATIVE_NUMBER("IncompleteGamma(double,double,double)");}
         if(Equality(X1,X2)==true){ return 0.;}
         if(Equality(X2,0.)==true){ return IncompleteGamma(N,X1);}

         if(fabs(X1-X2)>fabs(X1)){ return IncompleteGamma(N,X1,(X1+X2)/2.)+IncompleteGamma(N,(X1+X2)/2.,X2);}

         double G=0., dG, Y=X1-X2; int A=0;
         do{ dG=exp(-X2+(N-1.-A)*log(X2))*pow(-1.,A)*(1.-exp(-Y)*exp(A,Y))*exp( logGamma(1.-N+A)-logGamma(1.-N) );
             G+=dG;
             A++;
            }
         while(fabs(dG)>=numeric_limits<double>::epsilon()*fabs(G));
         return G;
        }

complex<double> IncompleteGamma(double N,complex<double> Z)
       { if(N<1. && Equality(fmod(N,1.),0.)==true){ throw NEGATIVE_INTEGER("IncompleteGamma(double,complex<double>)");}
         if(real(Z)<0.){ NEGATIVE_NUMBER("IncompleteGamma(double,complex<double>)");}
         if(Equality(abs(Z),0.)==true){ return Zero(Z);}
         if(N>2.){ return (N-1.)*IncompleteGamma(N-1.,Z)-exp(-Z+(N-1.)*log(Z));}
         if(N<1.){ return ( IncompleteGamma(N+1.,Z)+exp(-Z+N*log(Z)) )/N;}

         if(real(Z)>N+1.){ return Gamma(N)-Gamma(N,Z);}

         complex<double> dG=1./N, G=Zero(Z);
         double A=1.;
         do{ G+=dG; dG*=Z/(N+A); A++;}
         while(abs(dG)>=numeric_limits<double>::epsilon()*abs(G));
         return exp(-Z+N*log(Z))*G;
        }

// *******************************************************

double Factorial(int N)
       { const int Nmax=30;
         if(N<0){ throw NEGATIVE_INTEGER("Factorial(int)");}
         if(N>30-1){ return N*Factorial(N-1);}         
         static int Ntop=5;
         static double F[Nmax]={ 1., 1., 2., 6., 24.};
         while(Ntop<=N){ F[Ntop]=Ntop*F[Ntop-1]; Ntop++;};
         return F[N]; 
        }

double Factorial(double N){ return Gamma(N+1.);}

double DoubleFactorial(int N)
       { if(N%2==0){ return pow(2.,N)*Factorial(N);}
         else{ return pow(2.,N)/sqrt(M_PI)*Gamma(N+0.5);}
        }

double RisingFactorial(double N,double X){ return exp( logGamma(X+N) -logGamma(X) );}

double Pochhammer(double N,double X){ return RisingFactorial(N,X);}

double FallingFactorial(double N,double X){ return exp( logGamma(X+1.) -logGamma(X-N+1.) );}

// ********************************************************************

double Beta(double a,double b)
//       { return std::beta(a,b);}
       { return exp( logGamma(a)+logGamma(b)-logGamma(a+b) );}

double IncompleteBeta(double a,double b,double X)
       { if(X<0. || X>1.){ throw OUT_OF_RANGE<double>(X,0.,1.,"IncompleteBeta");}
         if(Equality(X,0.)==true){ return 0.;}
         if(Equality(X,1.)==true){ return Beta(a,b);}
         if( X>(a+1.)/(a+b+2.) ){ return Beta(b,a)-IncompleteBeta(b,a,1.-X);}

         double C=1., D=1.-(a+b)*X/(a+1.); 
         if(fabs(D)<numeric_limits<double>::epsilon()){ D=numeric_limits<double>::epsilon();}
         D=1./D;
         double B=D;

         int m=1;
         do{ D=1.+ D*m*(b-m)*X/(a+2.*m-1.)/(a+2.*m);
             if(fabs(D)<numeric_limits<double>::epsilon()){ D=numeric_limits<double>::epsilon();}
             D=1./D;
             C=1.+m*(b-m)*X/(a+2.*m-1.)/(a+2.*m)/C;
             if(fabs(C)<numeric_limits<double>::epsilon()){ C=numeric_limits<double>::epsilon();}
             B*=C*D;

             D=1.-D*(a+m)*(a+b+m)*X/(a+2.*m)/(a+2.*m+1.);
             if(fabs(D)<numeric_limits<double>::epsilon()){ D=numeric_limits<double>::epsilon();}
             D=1./D;
             C=1.-(a+m)*(a+b+m)*X/(a+2.*m)/(a+2.*m+1.)/C;
             if(fabs(C)<numeric_limits<double>::epsilon()){ C=numeric_limits<double>::epsilon();}
             B*=C*D;

             ++m;
            }while( fabs(C*D-1.)>10.*numeric_limits<double>::epsilon() );

         return B*exp( a*log(X)+b*log(1.-X) )/a;
        }

double IncompleteBeta(double a,double b,double X1,double X2){ return IncompleteBeta(a,b,X1)-IncompleteBeta(a,b,X2);}

double RegularizedBeta(double a,double b,double X){ return IncompleteBeta(a,b,X)/Beta(a,b);}

// *******************************************************

//double BinomialCoefficient(int n,int k){ return Factorial(n)/Factorial(k)/Factorial(n-k);}

double BinomialCoefficient(double n,double k){ return Gamma(n+1.)/Gamma(k+1.)/Gamma(n-k+1.);}

double BinomialSeries(double X,double n)
       { if(fabs(X)>1.){ throw OUT_OF_RANGE<double>(X,-1.,1.,"BinomialSeries");}
         if(Equality(fabs(X),0.)==true){ return X;}
         if(Equality(n,1.)==true){ return X;}

         double ipart, Y;
         if(Equality(modf(fabs(n),&ipart),0.)==true)
           { double k=1.,dY=n*X; 
             Y=dY;             
             do{ k++; dY*=(n+1.-k)/k*X; Y+=dY;}
             while(fabs(dY)>=numeric_limits<double>::epsilon()*fabs(Y));
            }
         else{ double k=1.,dY=BinomialCoefficient(n,k)*X;
               Y=dY;
               do{ k++; dY=BinomialCoefficient(n,k)*pow(X,k); Y+=dY;}
               while(fabs(dY)>=numeric_limits<double>::epsilon()*fabs(Y));
              }
         return Y;
        }








                  
        



