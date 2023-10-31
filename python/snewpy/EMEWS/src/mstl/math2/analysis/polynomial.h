
#include <cstdarg>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "mstl.h"

#if !defined(_POLYNOMIAL)
#define _POLYNOMIAL

// *******************************************************
// *******************************************************
// *******************************************************

class POLYNOMIAL;

// *******************************************************
// *******************************************************
// *******************************************************

// add the terms in a polynomial with coefficients C
double PolynomialSum(int N,double *C,double X);
double PolynomialSum(std::vector<double> C,double X);

// the Cs are the coefficients of the polynomial in increasing order, (X+A) is the root and R is the remainder
POLYNOMIAL SyntheticDivision(std::vector<double> C,double A,double &R);
POLYNOMIAL SyntheticDivision(std::vector<double> Xpoints,std::vector<double> Ypoints,double A,double &R);
POLYNOMIAL SyntheticDivision(POLYNOMIAL P,double A,double &R);

// the Cs are the coefficients of the polynomial in increasing order, A the coefficients of the divisor, and R is the remainder
POLYNOMIAL SyntheticDivision(std::vector<double> C,std::vector<double> A,std::vector<double> &R);
POLYNOMIAL SyntheticDivision(std::vector<double> Xpoints,std::vector<double> Ypoints,std::vector<double> A,std::vector<double> &R);
POLYNOMIAL SyntheticDivision(POLYNOMIAL P,POLYNOMIAL A,POLYNOMIAL &R);

// Neville's alogrithm for polynomial interpolation and extrapolation
//template <typename Type> Type PolynomialInterpolation(Type X,std::vector<Type> Xpoints,std::vector<Type> Ypoints,double ORDER=1.);
//template <typename Type> Type PolynomialInterpolation(Type X,int N,Type *Xpoints,Type *Ypoints,double ORDER=1.);

std::vector<double> PolynomialCoefficients(std::vector<double> Xpoints,std::vector<double> Ypoints);
std::vector<double> PolynomialCoefficients(int N,double *Xpoints,double *Ypoints);

// ***********************************************************
/*
double Hermite(int N,double X);
double UnitNormalizedHermite(int N,double X);

double Laguerre(int N,double X);

double Legendre(int L,double X);
double UnitNormalizedLegendre(int L,double X);
double AssociatedLegendre(int L,int M,double X);

double ChebyshevI(int N,double X);
double UnitNormalizedChebyshevI(int N,double X);
double ChebyshevII(int N,double X);
double DiscreteChebyshev(int N,int n,double X);
double UnitNormalizedDiscreteChebyshev(int N,int n,double X);
*/
// ***********************************************************
// ***********************************************************
// ***********************************************************

std::vector<std::complex<double> > QuarticRoots(double A4,double A3,double A2,double A1,double A0,double Y=0.);
std::vector<double> RealQuarticRoots(double A4,double A3,double A2,double A1,double A0,double Y=0.);

std::vector<std::complex<double> > CubicRoots(double A3,double A2,double A1,double A0,double Y=0.);
std::vector<double> RealCubicRoots(double A3,double A2,double A1,double A0,double Y=0.);
double SingleRealCubicRoot(double A3,double A2,double A1,double A0,double Y=0.);

// ***************************************************************

std::vector<std::complex<double> > QuadraticRoots(std::complex<double> A2,std::complex<double> A1,std::complex<double> A0,std::complex<double> Y=std::complex<double>(0.));
std::vector<std::complex<double> > QuadraticRoots(double A2,double A1,double A0,double Y=0.);
std::vector<std::complex<double> > QuadraticRoots(POLYNOMIAL P,double Y=0.);

std::vector<double> RealQuadraticRoots(double A2,double A1,double A0,double Y=0.);
std::vector<double> RealQuadraticRoots(POLYNOMIAL P,double Y=0.);

// ***************************************************************

std::vector<std::complex<double> > LinearRoots(std::complex<double> A1,std::complex<double> A0,std::complex<double> Y=std::complex<double>(0.)); // coefficients may be complex

std::vector<double> RealLinearRoots(double A1,double A0,double Y=0.); 
std::vector<double> RealLinearRoots(POLYNOMIAL P,double Y=0.);

// ***************************************************************

// roots of polynomials using Bairstow's method
std::vector<double> RealPolynomialRoots(std::vector<double> A,double Y=0.);
std::vector<double> RealPolynomialRoots(POLYNOMIAL P,double Y=0.);

// *******************************************************
// *******************************************************
// *******************************************************

class POLYNOMIAL // the coefficients are held as P = C0 + C1.X + C2.X^2 ...
      { protected : std::vector<double> coefficients;

                    void Create(int NC=0){ coefficients=std::vector<double>(NC);}
                    void Copy(std::vector<double> C);
                    void Initialize(void);
                    void Fill(double C0,va_list &val);

          public : POLYNOMIAL(void){ Create();}
                   POLYNOMIAL(int N,double C0,...){ Create(N); va_list val; va_start(val,C0); Fill(C0,val); va_end(val);}
                   POLYNOMIAL(int N){ Create(N); Initialize();}
                   POLYNOMIAL(std::vector<double> C){ Create(C.size()); Copy(C);}
                   POLYNOMIAL(std::vector<double> Xpoints,std::vector<double> Ypoints);

                   int N(void) const { return coefficients.size();}
                   int Order(void) const { return N()-1;}

                   double operator[](int i) const { return coefficients[i];}
                   double& operator[](int i) { return coefficients[i];}
                   std::vector<double> Coefficients(void) const { return coefficients;}

                   double& SetCoefficient(int i) { return coefficients[i];}
                   void SetCoefficient(int i,double Ci){ coefficients[i]=Ci;}

                   double operator()(double X){ return Evaluate(X);}
                   double Evaluate(double X);

                   double Derivative(double X);
         };

#endif
