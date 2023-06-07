#if !defined(_ODE_INTEGRATORS_TEMPLATES)
#define _ODE_INTEGRATORS_TEMPLATES

#include "mstl.h"

template<typename XType,typename YType> YType SecondOrderFirstReduction(XType,YType,YType dYdx){ return dYdx;}

template<typename XType,typename YType> YType ThirdOrderFirstReduction(XType,YType,YType dYdx,YType){ return dYdx;}
template<typename XType,typename YType> YType ThirdOrderSecondReduction(XType,YType,YType,YType d2Ydx2){ return d2Ydx2;}

template<typename XType,typename YType> YType NthReduction(XType,std::vector<YType> dYdx,int N){ return dYdx[N];}

#endif


