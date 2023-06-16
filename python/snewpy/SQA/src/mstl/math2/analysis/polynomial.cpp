#include "polynomial.h"

using std::complex;
using std::vector;
using std::numeric_limits;

void POLYNOMIAL::Copy(vector<double> v){ coefficients=v;}

void POLYNOMIAL::Initialize(void){ for(int i=0;i<=(int)N()-1;i++){ coefficients[i]=0.;} }

void POLYNOMIAL::Fill(double C0,va_list &val){ coefficients[0]=C0; for(int i=1;i<=(int)N()-1;i++){ coefficients[i]=va_arg(val,double);} }

POLYNOMIAL::POLYNOMIAL(vector<double> Xpoints,vector<double> Ypoints)
            { Create(Xpoints.size());
              try{ DCVECTOR C(VandermondeMatrixInverse(Xpoints)*DCVECTOR(Ypoints));
                   for(int i=0;i<=(int)N()-1;i++){ coefficients[i]=C[i];}
                  }
              catch(DIFFERENT_LENGTHS &DL){ DL.Change("POLYNOMIAL(vector<double>,vector<double>)","POLYNOMIAL"); throw DL;}
             }

double POLYNOMIAL::Evaluate(double X)
       { double Y=coefficients.back(); for(int i=(int)N()-2;i>=0;i--){ Y*=X; Y+=coefficients[i];} return Y;}

double POLYNOMIAL::Derivative(double X)
       { double dYdX=coefficients.back()*N(); for(int i=(int)N()-2;i>=1;i--){ dYdX*=X; dYdX+=i*coefficients[i];} return dYdX;}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

double PolynomialSum(int N,double *C,double X)
       { return PolynomialSum(vector<double>(C,C+N),X);}

double PolynomialSum(vector<double> C,double X)
       { double Y=C.back();;
         for(int i=(int)C.size()-2;i>=0;i--){ Y*=X; Y+=C[i];}
         return Y;
        }

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

// the Cs are the coefficients of the polynomial in increasing order, (X+A) is the root and R is the remainder
POLYNOMIAL SyntheticDivision(vector<double> C,double A,double &R)
       { return SyntheticDivision(POLYNOMIAL(C),A,R);}

POLYNOMIAL SyntheticDivision(vector<double> Xpoints,vector<double> Ypoints,double A,double &R)
       { if(Xpoints.size()!=Ypoints.size()){ throw DIFFERENT_LENGTHS("SyntheticDivision");}
         return SyntheticDivision(PolynomialCoefficients(Xpoints,Ypoints),A,R);
        }

POLYNOMIAL SyntheticDivision(POLYNOMIAL P,double A,double &R)
       { POLYNOMIAL Q(P.N()-1);
         R=P[P.N()-1];
         for(int i=P.N()-2;i>=0;i--){ Q[i]=R; R=P[i]-Q[i]*A;}
         return Q;
        }

//************************************************

POLYNOMIAL SyntheticDivision(vector<double> C,vector<double> A,vector<double> &R)
       { POLYNOMIAL RR(R);
         POLYNOMIAL Q=SyntheticDivision(POLYNOMIAL(C),POLYNOMIAL(A),RR);
         R=RR.Coefficients();
         return Q;
        }

POLYNOMIAL SyntheticDivision(vector<double> Xpoints,vector<double> Ypoints,vector<double> A,vector<double> &R)
       { if(Xpoints.size()!=Ypoints.size()){ throw DIFFERENT_LENGTHS("SyntheticDivision");}
         return SyntheticDivision(PolynomialCoefficients(Xpoints,Ypoints),A,R);
        }

POLYNOMIAL SyntheticDivision(POLYNOMIAL P,POLYNOMIAL A,POLYNOMIAL &R)
       { int NP=P.N(), NA=A.N(); 
         if(NA>NP){ R=P; return POLYNOMIAL();}
         int NR=NA-1, NQ=NP-NR;
         POLYNOMIAL Q(NQ);
         R=POLYNOMIAL(NR);        
         for(int j=NR-1;j>=0;j--){ R[j]=P[NP-NR+j];}
         for(int i=NQ-1;i>=0;i--)
            { Q[i]=R[NR-1]/A[NA-1]; 
              for(int j=NR-1;j>=1;j--){ R[j]=R[j-1]-Q[i]*A[j];} 
              R[0]=P[i]-Q[i]*A[0];
             } 
         return Q;
        }

//*******************************************************************************************************************

vector<double> PolynomialCoefficients(int N,double *Xpoints,double *Ypoints)
       { return PolynomialCoefficients(vector<double>(Xpoints,Xpoints+N),vector<double>(Ypoints,Ypoints+N));}

vector<double> PolynomialCoefficients(vector<double> Xpoints,vector<double> Ypoints)
       { if(Xpoints.size()!=Ypoints.size()){ throw DIFFERENT_LENGTHS("PolynomialCoefficients");}
         if(Xpoints.size()==1){ return Ypoints;}

         DCVECTOR C(VandermondeMatrixInverse(Xpoints)*DCVECTOR(Ypoints));
         return vector<double>(&C[0],&C[C.N()-1]);
        }

//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************

double Hermite(int N,double X)
//       { return std::hermite(N,X);}
       { if(N==0){ return 1.;}
         if(N==1){ return 2.*X;}

         double Hnplus1, Hn=2.*X, Hnminus1=1.;
         for(int n=1;n<=N-1;n++){ Hnplus1=2.*(X*Hn-n*Hnminus1); Hnminus1=Hn; Hn=Hnplus1;}

         return Hn;
        }

double UnitNormalizedHermite(int N,double X)
       { if(N==0){ return 1./pow(M_PI,0.25);}
         if(N==1){ return M_SQRT2*X/pow(M_PI,0.25);}

         double Hnplus1, Hn=M_SQRT2*X/pow(M_PI,0.25), Hnminus1=1./pow(M_PI,0.25);
         for(int n=1;n<=N-1;n++){ Hnplus1=(M_SQRT2*X*Hn-sqrt(n)*Hnminus1)/sqrt(n+1.); Hnminus1=Hn; Hn=Hnplus1;}

         return Hn;
        }

//*********************************

double Laguerre(int N,double X)
//       { return std::laguerre(N,X);}
       { if(fabs(X)<0.){ throw NEGATIVE_NUMBER("Laguerre");}

         if(N==0){ return 1.;}
         if(N==1){ return 1.-X;}

         double logplus1, log=1.-X, logminus1=1.;
         for(int n=1;n<=N-1;n++){ logplus1=2.*log-logminus1-((1.+X)*log-logminus1)/(n+1.); logminus1=log; log=logplus1;}

         return log;
        }

//*********************************

double Legendre(int L,double X)
//       { return std::legendre(N,X);}
       { if(fabs(X)>1.){ throw OUT_OF_RANGE<double>(X,-1.,1.,"Legendre");}

         if(L==0){ return 1.;}
         if(L==1){ return X;}

         double Plplus1, Pl=X, Plminus1=1.;
         for(int l=1;l<=L-1;l++){ Plplus1=2.*X*Pl-Plminus1-(X*Pl-Plminus1)/(l+1.); Plminus1=Pl; Pl=Plplus1;}

         return Pl;
        }

double UnitNormalizedLegendre(int L,double X)
       { try{ return Legendre(L,X)*sqrt(L+0.5);}
         catch(OUT_OF_RANGE<double> &OOR){ OOR.ChangeFunction("UnitNormalizedLegendre"); throw OOR;}
        }

double AssociatedLegendre(int L,int M,double X)
//       { return std::assoc_legendre(L,M,X);}
       { if(abs(M)>(int)L){ throw OUT_OF_RANGE<int>(M,-(int)L,L,"AssociatedLegendre");}
         if(fabs(X)>1.){ throw OUT_OF_RANGE<double>(X,-1.,1.,"AssociatedLegendre");}

         if(M<0){ return pow(-1.,M)*Factorial(L-M)/Factorial(L+M)*AssociatedLegendre(L,-M,X);}

         if(L==0){ return 1.;}
         double Plm,Pmm,Pmplus1m;

         Pmm=pow(-1.,M)*pow(1.-X*X,M/2.); for(int i=1;i<=2*M-1;i+=2){ Pmm*=i;} if((int)L==M){ return Pmm;}
         Pmplus1m=X*(2.*M+1.)*Pmm; if((int)L==M+1){ return Pmplus1m;}

         for(int l=M+2;l<=L;l++){ Plm=( (2*l+1)*X*Pmplus1m-(l+M)*Pmm)/(l+1-M); Pmm=Pmplus1m; Pmplus1m=Plm;}

         return Plm;
        }

//*********************************

double ChebyshevI(int N,double X)
       { if(fabs(X)>1.){ throw OUT_OF_RANGE<double>(X,-1.,1.,"ChebyshevI");}
         if(N==0){ return 1.;}
         if(N==1){ return X;}

         double Tnplus1, Tn=X, Tnminus1=1.;
         for(int n=1;n<=N-1;n++){ Tnplus1=2.*X*Tn-Tnminus1; Tnminus1=Tn; Tn=Tnplus1;}

         return Tn;
        }

double UnitNormalizedChebyshevI(int N,double X)
       { try{ if(N==0){ return ChebyshevI(N,X)*M_1_SQRTPI;} else{ return ChebyshevI(N,X)*M_SQRT2*M_1_SQRTPI;} }
         catch(OUT_OF_RANGE<double> &OOR){ OOR.ChangeFunction("UnitNormalizedChebyshevI"); throw OOR;}
        }

double ChebyshevII(int N,double X)
       { if(fabs(X)>1.){ throw OUT_OF_RANGE<double>(X,-1.,1.,"ChebyshevII");}
         if(N==0){ return 1.;}
         if(N==1){ return 2.*X;}

         double Unplus1, Un=2.*X, Unminus1=1.;
         for(int n=1;n<=N-1;n++){ Unplus1=2.*X*Un-Unminus1; Unminus1=Un; Un=Unplus1;}

         return Un;
        }

double DiscreteChebyshev(int N,int n,double X)
       { if(n>N-1){ throw OUT_OF_RANGE<int>(n,0,N-1,"DiscreteChebyshev");}
         if(X<0. || X>N-1.){ throw OUT_OF_RANGE<double>(X,0.,N-1.,"DiscreteChebyshev");}
         if(n==0){ return 1.;}
         if(n==1){ return 2.*X-N+1.;}

         double Tlplus1, Tl=2.*X-N+1., Tlminus1=1.;
         for(int l=1;l<=n-1;l++)
            { Tlplus1=(2.*X-N+1.)*Tl - (N-l)*(N+l)*Tlminus1 + ( l*(2.*X-N+1.)*Tl + (N-l)*(N+l)*Tlminus1 )/(l+1.);
              Tlminus1=Tl;
              Tl=Tlplus1;
             }

         return Tl;
        }

double UnitNormalizedDiscreteChebyshev(int N,int n,double X)
       { return sqrt(2.*n+1.)*exp( 0.5*logGamma(N-n)-0.5*logGamma(N+n+1.) )*DiscreteChebyshev(N,n,X);}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************

// solutions of the equation Y=A3 X^3 + A2 X^2 + A1 X + A0
vector<complex<double> > QuarticRoots(double A4,double A3,double A2,double A1,double A0,double Y)
      { if(Equality(A4,0.)==true){ return CubicRoots(A3,A2,A1,A0,Y);}

        A3/=A4; A2/=A4; A1/=A4; A0=(A0-Y)/A4;

        double p=-3.*A3*A3/8.+A2; 
        double q=A3*A3*A3/8.-A3*A2/2.+A1;
        double r=-3.*A3*A3*A3*A3/256.+A2*A3*A3/16. -A3*A1/4.+A0;

        complex<double> y=CubicRoots(1.,2.*p,p*p-4.*r,-q*q)[0];
        vector<complex<double> > Z(4,-complex<double>(A3/4.)),dZ;

        dZ=QuadraticRoots(1.,+sqrt(y),0.5*(p+y-q/sqrt(y)));
        Z[0]+=dZ[0]; Z[1]+=dZ[1];
        dZ=QuadraticRoots(1.,-sqrt(y),0.5*(p+y+q/sqrt(y)));
        Z[2]+=dZ[0]; Z[3]+=dZ[1];

        return Z;
       }

vector<double> RealQuarticRoots(double A4,double A3,double A2,double A1,double A0,double Y)
      { if(Equality(A4,0.)==true){ return RealCubicRoots(A3,A2,A1,A0,Y);}
        vector<complex<double> > Z(QuarticRoots(A4,A3,A2,A1,A0,Y));
        vector<double> X;

        for(int i=0;i<=static_cast<int>(Z.size())-1;i++){ if(Equality(imag(Z[i]),0.)==true){ X.push_back(real(Z[i]));} }

        return X;
       }

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************

// solutions of the equation Y=A3 X^3 + A2 X^2 + A1 X + A0
vector<complex<double> > CubicRoots(double A3,double A2,double A1,double A0,double Y)
      { if(Equality(A3,0.)==true){ return QuadraticRoots(A2,A1,A0,Y);}

        A2/=A3; A1/=A3; A0=(A0-Y)/A3;
        double Q=(3.*A1-A2*A2)/9., R=(9.*A2*A1-27.*A0-2.*pow(A2,3.))/54.;

        // A=S+T, B=S-T
        vector<complex<double> > Z(3);
        complex<double> A,B;
        double p,q,r,omega;

        if(Equality(R,0.)==true)
          { if(Q>=0.){ B=sqrt(Q);} else{ B=I*sqrt(-Q);} 
            Z[0]=-A2/3.;
            Z[1]=-A2/3. -I*M_SQRT3*B;
            Z[2]=-A2/3. +I*M_SQRT3*B;
           }
        else{ q=pow(Q,3.)/R/R;
              if(q<=-1.){ omega=acos(Sign(R)/sqrt(-q)); 
                          Z[0]=-A2/3. +2.*sqrt(-Q)*cos(omega/3.);
                          Z[1]=-A2/3. +2.*sqrt(-Q)*cos((omega+2.*M_2PI)/3.);
                          Z[2]=-A2/3. +2.*sqrt(-Q)*cos((omega+M_2PI)/3.);
                         }
             else{ r=cbrt(R);
                   if(fabs(q)<1e-7){ p=BinomialSeries(q,0.5);
                                     A=r*( M_CBRT2*(1.+BinomialSeries(p/2.,1./3.)) + cbrt(-p) );
                                     B=r*( M_CBRT2*(1.+BinomialSeries(p/2.,1./3.)) + cbrt(p) );
                                    }
                   if(q>1e7){ p=(1.+BinomialSeries(1./q,-0.5))/sqrt(q);
                              A=r*pow(1.+q,1./6.)*( BinomialSeries(p,1./3.) -BinomialSeries(-p,1./3.) );
                              B=r*pow(1.+q,1./6.)*( 2.+BinomialSeries(p,1./3.)+BinomialSeries(-p,1./3.) );
                             }
                   if( (q<=-1e-7 && q>=-1.) || (q>=1e-7 && q<=1e7) )
                     { A=r*( cbrt(1.+sqrt(1.+q)) + cbrt(1.-sqrt(1.+q)) );
                       B=r*( cbrt(1.+sqrt(1.+q)) - cbrt(1.-sqrt(1.+q)) );
                      }

                   Z[0]=-A2/3. +A;
                   Z[1]=-A2/3. -(A +I*M_SQRT3*B)/2.;
                   Z[2]=-A2/3. -(A -I*M_SQRT3*B)/2.;
                  }
             }

        return Z;
       }

vector<double> RealCubicRoots(double A3,double A2,double A1,double A0,double Y)
      { if(Equality(A3,0.)==true){ return RealQuadraticRoots(A2,A1,A0,Y);}
        vector<complex<double> > Z(CubicRoots(A3,A2,A1,A0,Y));
        vector<double> X;

        for(int i=0;i<=static_cast<int>(Z.size())-1;i++){ if(Equality(imag(Z[i]),0.)==true){ X.push_back(real(Z[i]));} }

        return X;
       }

double SingleRealCubicRoot(double A3,double A2,double A1,double A0,double Y)
      { A2/=A3; A1/=A3; A0=(A0-Y)/A3;
        double Q=(3.*A1-A2*A2)/9., R=(9.*A2*A1-27.*A0-2.*pow(A2,3.))/54.;

        // A=S+T
        double Z;
        double A;
        double p,q,r;

        if(Equality(R,0.)==true)
          { Z=-A2/3.;}
        else{ q=pow(Q,3.)/R/R;
              r=cbrt(R);
              if(fabs(q)<1e-7){ p=BinomialSeries(q,0.5);
                                A=r*( M_CBRT2*(1.+BinomialSeries(p/2.,1./3.)) + cbrt(-p) );
                               }
              if(q>1e7){ p=(1.+BinomialSeries(1./q,-0.5))/sqrt(q);
                         A=r*pow(1.+q,1./6.)*( BinomialSeries(p,1./3.) -BinomialSeries(-p,1./3.) );
                        }
              if( (q<=-1e-7 && q>=-1.) || (q>=1e-7 && q<=1e7) )
                { A=r*( cbrt(1.+sqrt(1.+q)) + cbrt(1.-sqrt(1.+q)) );}

              Z=-A2/3. +A;
             }

        return Z;
       }

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************

// solutions of the equation Y=A2 Z^2 + A1 Z + A0
vector<complex<double> > QuadraticRoots(complex<double> A2,complex<double> A1,complex<double> A0,complex<double> Y)
      { if(Equality(norm(A2),0.)==true){ return LinearRoots(A1,A0,Y);}

        A0-=Y;
        complex<double> q=4.*A2*A0/(A1*A1);

        vector<complex<double> > Z(2);
        Z[0]=-A1*(1.+sqrt(1.-q))/2./A2;
        Z[1]=-A1*(1.-sqrt(1.-q))/2./A2;

        return Z;
       }

// solutions of the equation Y=A2 X^2 + A1 X + A0
vector<complex<double> > QuadraticRoots(double A2,double A1,double A0,double Y)
      { if(Equality(A2,0.)==true){ return LinearRoots(A1,A0,Y);}

        A0-=Y;
        double q=4.*A2*A0/(A1*A1);

        vector<complex<double> > Z(2);

        if(fabs(q)<1e-7){ Z[0]=A1*BinomialSeries(-q,0.5)/2./A2;
                          Z[1]=-A1*( 2.+BinomialSeries(-q,0.5) )/2./A2;
                          return Z;
                         }
        if(q<0. && fabs(q)>1e7){ Z[0]=-A1*( 1.+sqrt(-q)+sqrt(-q)*BinomialSeries(-1./q,0.5) )/2./A2;
                                 Z[1]=-A1*( 1.-sqrt(-q)-sqrt(-q)*BinomialSeries(-1./q,0.5) )/2./A2;
                                 return Z;
                                }
        if(q>0. && fabs(q)>1e7){ Z[0]=-A1*( 1. + I*sqrt(q) +I*sqrt(q)*BinomialSeries(-1./q,0.5) )/2./A2;
                                 Z[1]=-A1*( 1. - I*sqrt(q) -I*sqrt(q)*BinomialSeries(-1./q,0.5) )/2./A2;
                                 return Z;
                                }
        if(q<1.){ Z[0]=-A1*(1.+sqrt(1.-q))/2./A2;
                  Z[1]=-A1*(1.-sqrt(1.-q))/2./A2;
                  return Z;
                 }
        else{ Z[0]=-A1*(1. + I*sqrt(q-1.) )/2./A2;
              Z[1]=-A1*(1. - I*sqrt(q-1.) )/2./A2;
              return Z;
             }
       }

vector<complex<double> > QuadraticRoots(POLYNOMIAL P,double Y)
      { return QuadraticRoots(P[2],P[1],P[0],Y);}

//******************************

vector<double> RealQuadraticRoots(double A2,double A1,double A0,double Y)
      { if(Equality(A2,0.)==true){ return RealLinearRoots(A1,A0,Y);}

        A0-=Y;
        double q=4.*A2*A0/(A1*A1);

        if(q>1.){ return vector<double>();}
        
        vector<double> X(2);
        // |q| << 1
        if(fabs(q)<1e-7){ X[0]=A1*BinomialSeries(-q,0.5)/2./A2;
                          X[1]=-A1*( 2.+BinomialSeries(-q,0.5) )/2./A2;
                          return X;
                         }
        // q < 0 and |q|>>1 
        if(q<0. && fabs(q)>1e7){ X[0]=-A1*( 1.+sqrt(-q)+sqrt(-q)*BinomialSeries(-1./q,0.5) )/2./A2;
                                 X[1]=-A1*( 1.-sqrt(-q)-sqrt(-q)*BinomialSeries(-1./q,0.5) )/2./A2;
                                 return X;
                                }
        // usual case 
        X[0]=-A1*(1.+sqrt(1.-q))/2./A2;
        X[1]=-A1*(1.-sqrt(1.-q))/2./A2;
        return X;
       }

vector<double> RealQuadraticRoots(POLYNOMIAL P,double Y)
      { return RealQuadraticRoots(P[2],P[1],P[0],Y);}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************


// solutions of the equation Y=A1 X + A0
vector<complex<double> > LinearRoots(complex<double> A1,complex<double> A0,complex<double> Y)
      { return vector<complex<double> >(1,(Y-A0)/A1);}

//******************************

vector<double> RealLinearRoots(double A1,double A0,double Y)
      { return vector<double>(1,(Y-A0)/A1);}

vector<double> RealLinearRoots(POLYNOMIAL P,double Y)
      { return vector<double>(1,(Y-P[0])/P[1]);}

//************************************************************************************************
//************************************************************************************************
//************************************************************************************************

vector<double> RealPolynomialRoots(vector<double> A,double Y)
      { return RealPolynomialRoots(POLYNOMIAL(A),Y);}

vector<double> RealPolynomialRoots(POLYNOMIAL P0,double Y)
      { vector<double> roots,temproots;

        if(P0.N()==3){ return RealQuadraticRoots(P0,Y);}
        if(P0.N()==2){ return RealLinearRoots(P0,Y);}
        if(P0.N()==1){ return vector<double>();}

        POLYNOMIAL P(P0); 
        P[0]-=Y;

        POLYNOMIAL U(3); 
        U[2]=1.;
        U[1]=P[P.N()-2]/P[P.N()-1];
        U[0]=P[P.N()-3]/P[P.N()-1];

        vector<double> dU(2); 
        
        POLYNOMIAL R1,R2;
        double lambda;
        while(P.N()>=4)
             { do{ SyntheticDivision(SyntheticDivision(P,U,R1),U,R2); 

                   lambda=U[0]*R2[1]*R2[1] +R2[0]*(R2[0]-U[1]*R2[1]);

                   dU[0]=( R2[1]*R1[1]*U[0] +(R2[0]-U[1]*R2[1])*R1[0] )/lambda;
                   dU[1]=( R2[0]*R1[1] -R2[1]*R1[0] )/lambda; 
                   U[0]+=dU[0];
                   U[1]+=dU[1];
                  }while( fabs(dU[0])>4.*fabs(U[0])*numeric_limits<double>::epsilon() && fabs(dU[1])>4.*fabs(U[1])*numeric_limits<double>::epsilon() );  

               temproots=RealQuadraticRoots(U);
               if(temproots.empty()==false)
                 { if(roots.empty()==true){ roots=temproots;} 
                   else{ roots.push_back(temproots[0]); roots.push_back(temproots[1]);}
                  }
               P=SyntheticDivision(P,U,R1);
              }
        if(P.N()==3){ temproots=RealQuadraticRoots(P); if(temproots.empty()==false){ roots.push_back(temproots[0]); roots.push_back(temproots[1]);} }
        if(P.N()==2){ temproots=RealLinearRoots(P); roots.push_back(temproots[0]);}

        return roots;
       }