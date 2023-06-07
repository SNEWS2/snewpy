
#include<cmath>
#include<cstdlib>

#include<limits>
#include<vector>
#include<complex>

#include "mstl.h"

#if !defined(_SPECIALFUNCTIONS99)
#define _SPECIALFUNCTIONS99

double exp(int n,double X); // the truncated exponential function

double Gamma(double N); // use the standard C function
double logGamma(double N); //natural log // use the standard C function
std::complex<double> Gamma(std::complex<double> N);
std::complex<double> logGamma(std::complex<double> N); //natural log

double Gamma(double N,double X); // This is the upper incomplete Gamma function
double logGamma(double N,double X);
double Gamma(double N,double X1,double X2); // Gamma(N,X1)-Gamma(N,X2)
std::complex<double> Gamma(double N,std::complex<double> Z); // This is the upper incomplete Gamma function

double IncompleteGamma(double N,double X); // This is the lower incomplete gamma function
double logIncompleteGamma(double N,double X);
double IncompleteGamma(double N,double X1,double X2); // gamma(N,X1)-gamma(N,X2)
std::complex<double> IncompleteGamma(double N,std::complex<double> Z); // This is the lower incomplete gamma function

double Factorial(int N);
double Factorial(double N);

double DoubleFactorial(int N);
double RisingFactorial(double N,double X);
double Pochhammer(double N,double X);
double FallingFactorial(double N,double X);

double Beta(double a,double b);
double IncompleteBeta(double a,double b,double X); 
double IncompleteBeta(double a,double b,double X1,double X2); // Beta(A,B,X1)-Beta(A,B,X2)
double RegularizedBeta(double a,double b,double X);

double BinomialCoefficient(double n,double k);
double BinomialSeries(double x,double n); // This returns the expansion of (1+x)^n - 1


#endif
