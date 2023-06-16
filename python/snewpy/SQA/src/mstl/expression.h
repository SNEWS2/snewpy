
#include<cstdarg>

#include "mstl.h"

#if !defined(_EXPRESSION)
#define _EXPRESSION

// ******************************************************************
// ******************************************************************
// ******************************************************************

template<typename Type> class VARIABLE;
template<typename Type> class CONSTANT;

template <typename RType,typename XType,typename AType,typename OType> class UNARYEXPRESSION;

template <typename RType,typename XType1,typename XType2,typename AType,typename OType> class BINARYEXPRESSION;

template <typename RType,typename AType,typename EType> class EXPRESSION;

// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************

template<typename Type>
class VARIABLE{ Type operator()(const Type &X){ return X;} };

template<typename Type>
class CONSTANT
      { private : const Type &c; 
        public: CONSTANT(const Type &C) : c(C) {;}
                CONSTANT(Type &C) : c(C) {;}

                Type operator()(void){ return c;} 
       }; 

// ******************************************************************
// ******************************************************************
// ******************************************************************

template <typename RType,typename XType,typename AType,typename OType>
class UNARYEXPRESSION
      { private : const XType &x;
                  template <typename YType> RType operator()(UNARYEXPRESSION<RType,YType,AType,OType>*,const AType &A) const; 
                  RType operator()(UNARYEXPRESSION<RType,CONSTANT<RType>,AType,OType>*,const AType &A) const; 

        public : UNARYEXPRESSION(const XType &X) : x(X) {;}

                 RType operator()(const AType &A) const { return (*this)(this,A);}
       };

// ******************************************************************
// ******************************************************************
// ******************************************************************

template <typename RType,typename XType1,typename XType2,typename AType,typename OType>
class BINARYEXPRESSION 
      { private : const XType1 &x1; const XType2 &x2;
                  template <typename YType1,typename YType2> RType operator()(BINARYEXPRESSION<RType,YType1,YType2,AType,OType>*,const AType &A) const; 
                  template <typename YType2> RType operator()(BINARYEXPRESSION<RType,CONSTANT<RType>,YType2,AType,OType>*,const AType &A) const; 
                  template <typename YType1> RType operator()(BINARYEXPRESSION<RType,YType1,CONSTANT<RType>,AType,OType>*,const AType &A) const; 
                  RType operator()(BINARYEXPRESSION<RType,CONSTANT<RType>,CONSTANT<RType>,AType,OType>*,const AType &A) const; 

        public : BINARYEXPRESSION(const XType1 &X1,const XType2 &X2) : x1(X1), x2(X2) {;}

                 RType operator()(const AType &A) const { return (*this)(this,A);}
       };

// ******************************************************************
// ******************************************************************
// ******************************************************************

template <typename RType,typename AType,typename EType>
class EXPRESSION
         { private : EType e;
           public : EXPRESSION(EType E) : e(E) {;}

                    RType operator()(const AType &A) const { return e(A);}
                    operator RType() const { return (*this)();}
          };

#endif
