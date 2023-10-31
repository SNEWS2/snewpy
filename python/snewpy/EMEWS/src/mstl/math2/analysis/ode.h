
#include<cstdarg>
#include<cmath>
#include<vector>

#include "mstl.h"

#if !defined(_ODE_INTEGRATORS)
#define _ODE_INTEGRATORS

template<typename XType,typename YType> YType SecondOrderFirstReduction(XType,YType,YType dYdx);

template<typename XType,typename YType> YType ThirdOrderFirstReduction(XType,YType,YType dYdx,YType);
template<typename XType,typename YType> YType ThirdOrderSecondReduction(XType,YType,YType,YType d2Ydx2);

template<typename XType,typename YType> YType NthReduction(XType,std::vector<YType> dYdx,int N);

// *************************************************************************
// *************************************************************************
// *************************************************************************

#endif



