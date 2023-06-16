#if !defined(_INTERPOLATION_BASE)
#include "interpolation base.h"
#endif

#if !defined(_INTERPOLATION_BASE_FUNCTIONS)
#define _INTERPOLATION_BASE_FUNCTIONS

namespace interpolation{

template <typename IType> EXPRESSION_NEGATEEXPRESSION<EXPRESSION<IType> > operator-(const EXPRESSION<IType> &C)
      { return EXPRESSION_NEGATEEXPRESSION<EXPRESSION<IType> >(C.X(),C);}

// ************************

template <typename IType> EXPRESSION_EXPRESSIONSCALAR_ADDITION<EXPRESSION<IType> > operator+(const EXPRESSION<IType> &C,const double &D)
      { return EXPRESSION_EXPRESSIONSCALAR_ADDITION<EXPRESSION<IType> >(C.X(),C,D);}

template <typename IType> EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION<EXPRESSION<IType> > operator-(const EXPRESSION<IType> &C,const double &D)
      { return EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION<EXPRESSION<IType> >(C.X(),C,D);}

template <typename IType> EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION<EXPRESSION<IType> > operator*(const EXPRESSION<IType> &C,const double &D)
      { return EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION<EXPRESSION<IType> >(C.X(),C,D);}

template <typename IType> EXPRESSION_EXPRESSIONSCALAR_DIVISION<EXPRESSION<IType> > operator/(const EXPRESSION<IType> &C,const double &D)
      { return EXPRESSION_EXPRESSIONSCALAR_DIVISION<EXPRESSION<IType> >(C.X(),C,D);}

// ************************

template <typename IType> EXPRESSION_SCALAREXPRESSION_ADDITION<EXPRESSION<IType> > operator+(const double &D,const EXPRESSION<IType> &C)
      { return EXPRESSION_SCALAREXPRESSION_ADDITION<EXPRESSION<IType> >(C.X(),D,C);}

template <typename IType> EXPRESSION_SCALAREXPRESSION_SUBTRACTION<EXPRESSION<IType> > operator-(const double &D,const EXPRESSION<IType> &C)
      { return EXPRESSION_SCALAREXPRESSION_SUBTRACTION<EXPRESSION<IType> >(C.X(),D,C);}

template <typename IType> EXPRESSION_SCALAREXPRESSION_MULTIPLICATION<EXPRESSION<IType> > operator*(const double &D,const EXPRESSION<IType> &C)
      { return EXPRESSION_SCALAREXPRESSION_MULTIPLICATION<EXPRESSION<IType> >(C.X(),D,C);}

// ************************

template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION<EXPRESSION<IType> > operator+(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double> &NF)
      { return EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION<EXPRESSION<IType> >(C.X(),C,NF);}

template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION<EXPRESSION<IType> > operator-(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double> &NF)
      { return EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION<EXPRESSION<IType> >(C.X(),C,NF);}

template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION<EXPRESSION<IType> > operator*(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double> &NF)
      { return EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION<EXPRESSION<IType> >(C.X(),C,NF);}

template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION<EXPRESSION<IType> > operator/(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double> &NF)
      { return EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION<EXPRESSION<IType> >(C.X(),C,NF);}

// ************************

template <typename IType> EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION<EXPRESSION<IType> > operator+(const NULLARYFUNCTOR<double> &NF,const EXPRESSION<IType> &C)
      { return EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION<EXPRESSION<IType> >(C.X(),NF,C);}

template <typename IType> EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION<EXPRESSION<IType> > operator-(const NULLARYFUNCTOR<double> &NF,const EXPRESSION<IType> &C)
      { return EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION<EXPRESSION<IType> >(C.X(),NF,C);}

template <typename IType> EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION<EXPRESSION<IType> > operator*(const NULLARYFUNCTOR<double> &NF,const EXPRESSION<IType> &C)
      { return EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION<EXPRESSION<IType> >(C.X(),NF,C);}

// ************************

template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION<EXPRESSION<IType> > operator+(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double> &UF)
      { return EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION<EXPRESSION<IType> >(C.X(),C,UF);}

template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION<EXPRESSION<IType> > operator-(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double> &UF)
      { return EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION<EXPRESSION<IType> >(C.X(),C,UF);}

template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION<EXPRESSION<IType> > operator*(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double> &UF)
      { return EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION<EXPRESSION<IType> >(C.X(),C,UF);}

template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION<EXPRESSION<IType> > operator/(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double> &UF)
      { return EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION<EXPRESSION<IType> >(C.X(),C,UF);}

// ************************

template <typename IType> EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION<EXPRESSION<IType> > operator+(const UNARYFUNCTOR<double> &UF,const EXPRESSION<IType> &C)
      { return EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION<EXPRESSION<IType> >(C.X(),UF,C);}

template <typename IType> EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION<EXPRESSION<IType> > operator-(const UNARYFUNCTOR<double> &UF,const EXPRESSION<IType> &C)
      { return EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION<EXPRESSION<IType> >(C.X(),UF,C);}

template <typename IType> EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION<EXPRESSION<IType> > operator*(const UNARYFUNCTOR<double> &UF,const EXPRESSION<IType> &C)
      { return EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION<EXPRESSION<IType> >(C.X(),UF,C);}

// *******************************************************
// *******************************************************
// *******************************************************

template <typename IType1,typename IType2> 
std::vector<double> OverlapX(const EXPRESSION<IType1> &C1,const EXPRESSION<IType2> &C2)
     { double xmin=max(C1.XMin(),C2.XMin()), xmax=min(C1.XMax(),C2.XMax()); 
       if(xmin>xmax){ throw NO_OVERLAP("OverlapX");}
       int a1min=C1.XInterval(xmin),a1max=C1.XInterval(xmax), a2min=C2.XInterval(xmin),a2max=C2.XInterval(xmax);

       std::vector<double> points(1,xmin); 

       int a1=a1min+1, a2=a2min+1;

       while(C1.X(a1)<xmax && C2.X(a2)<xmax)
            { if(C1.X(a1)==C2.X(a2)){ points.push_back(C1.X(a1)); a1++; a2++;}
              else{ if(C1.X(a1)<C2.X(a2)){ points.push_back(C1.X(a1)); a1++;}
                    else{ if(C1.X(a1)>C2.X(a2)){ points.push_back(C2.X(a2)); a2++;} }
                   }
             };
       points.push_back(xmax);

       return points;
      }

// *******************************************************
// *******************************************************
// *******************************************************

template <typename IType1,typename IType2> 
EXPRESSION_EXPRESSIONEXPRESSION_ADDITION<EXPRESSION<IType1>,EXPRESSION<IType2> > operator+(const EXPRESSION<IType1> &C1,const EXPRESSION<IType2> &C2)
      { return EXPRESSION_EXPRESSIONEXPRESSION_ADDITION<EXPRESSION<IType1>,EXPRESSION<IType2> >(OverlapX(C1,C2),C1,C2);}

template <typename IType1,typename IType2> 
EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION<EXPRESSION<IType1>,EXPRESSION<IType2> > operator-(const EXPRESSION<IType1> &C1,const EXPRESSION<IType2> &C2)
      { return EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION<EXPRESSION<IType1>,EXPRESSION<IType2> >(OverlapX(C1,C2),C1,C2);}

template <typename IType1,typename IType2> 
EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION<EXPRESSION<IType1>,EXPRESSION<IType2> > operator*(const EXPRESSION<IType1> &C1,const EXPRESSION<IType2> &C2)
      { return EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION<EXPRESSION<IType1>,EXPRESSION<IType2> >(OverlapX(C1,C2),C1,C2);}

template <typename IType1,typename IType2> 
EXPRESSION_EXPRESSIONEXPRESSION_DIVISION<EXPRESSION<IType1>,EXPRESSION<IType2> > operator/(const EXPRESSION<IType1> &C1,const EXPRESSION<IType2> &C2)
      { return EXPRESSION_EXPRESSIONEXPRESSION_DIVISION<EXPRESSION<IType1>,EXPRESSION<IType2> >(OverlapX(C1,C2),C1,C2);}

} // end of namespace interpolation

#endif

