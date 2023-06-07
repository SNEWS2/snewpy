
#include<cmath>
#include<vector>

#include "mstl.h"

#if !defined(_INTERPOLATION)
#define _INTERPOLATION

// Neville's alogrithm for polynomial interpolation and extrapolation
// The XPOWER parameter controls the powers of x in the polynomial
template <typename XType,typename YType> YType PolynomialInterpolation(XType X,std::vector<XType> Xpoints,std::vector<YType> Ypoints,double XPOWER=1.);
template <typename XType,typename YType> YType PolynomialInterpolation(XType X,int N,XType *Xpoints,YType *Ypoints,double XPOWER=1.);

// The XPOWER parameter controls the powers of x in the polynomial
// The DENEXP paramater determines the highest term in the denominator's polynomial
// This routine is based upon the formulae in Stoer and Burlisch p69,70
template <typename XType,typename YType> YType RationalInterpolation(XType X,std::vector<XType> Xpoints,std::vector<YType> Ypoints,int DENEXP,double XPOWER=1.);
template <typename XType,typename YType> YType RationalInterpolation(XType X,int N,XType *Xpoints,YType *Ypoints,int DENEXP,double XPOWER=1.);

// The XPOWER parameter controls the powers of x in the polynomial
template <typename XType,typename YType> YType DiagonalRationalInterpolation(XType X,std::vector<XType> Xpoints,std::vector<YType> Ypoints,double XPOWER=1.);
template <typename XType,typename YType> YType DiagonalRationalInterpolation(XType X,int N,XType *Xpoints,YType *Ypoints,double XPOWER=1.);

// *************************************************************************
// *************************************************************************
// *************************************************************************

#endif
