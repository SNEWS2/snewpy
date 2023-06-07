
#include<cstdarg>
#include<cmath>
#include<vector>

#include "mstl.h"

#if !defined(_RUNGE_KUTTA)
#define _RUNGE_KUTTA

void RungeKuttaCashKarpParameters(int &NRK,int &NOrder,const double* &A,const double** &B,const double* &C,const double* &D);


template <typename XType,typename YType,typename Functor> 
YType RungeKuttaCashKarp(const Functor &f,YType &error,XType x0,YType y0,XType dx);

// *******************************************************************

// These 4 functions are the workhorses of the RK routines. x is the integration variable, y the dependent variables
// f the list of derivative functions or function that returns a list of derivaties, x0 and y0 the starting values
// dx the interval, A is the x increments for each step, B the K update, C the final combination of Ks and
// D the final combination of Ks for the embedded routine (if there is one). the error estimates are inserted
// into error and the functions return the increments for y.

template <typename XType,typename YType,typename Functor>
YType RungeKutta(const Functor &f,YType &error,XType x0,YType y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D);

// ******************

template <typename XType,typename YType,typename Functor>
std::vector<YType> RungeKutta(std::vector<Functor> &f,std::vector<YType> &error,
                              XType x0,std::vector<YType> y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D);

// ******************

// for the case where the function returns a list rather than a list of functions
template <typename XType,typename YType,typename Functor>
std::vector<YType> RungeKutta(const Functor &f,std::vector<YType> &error,XType x0,std::vector<YType> y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D);

// ******************

// for the case where a function returns a multicolumn list rather than a list of functions
template <typename XType,typename YType,typename Functor>
std::vector<std::vector<YType> > RungeKutta(const Functor &f,std::vector<std::vector<YType> > &error,XType x0,std::vector<std::vector<YType> > y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D);

// ******************

// for the case where a function returns a multicolumn list rather than a list of functions
template <typename XType,typename YType,typename Functor>
std::vector<std::vector<std::vector<YType> > > RungeKutta(const Functor &f,std::vector<std::vector<std::vector<YType> > > &error,XType x0,std::vector<std::vector<std::vector<YType> > > y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D);

// *****************

#endif



