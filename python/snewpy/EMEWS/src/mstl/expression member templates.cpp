#if !defined(_EXPRESSIONMEMBERS)
#define _EXPRESSIONMEMBERS

template <typename RType,typename XType,typename AType,typename OType>
template <typename YType> 
RType UNARYEXPRESSION<RType,XType,AType,OType>::operator()(UNARYEXPRESSION<RType,YType,AType,OType>*,const AType &A) const
      { return OType::operator()(x(A));}

template <typename RType,typename XType,typename AType,typename OType>
RType UNARYEXPRESSION<RType,XType,AType,OType>::operator()(UNARYEXPRESSION<RType,CONSTANT<RType>,AType,OType>*,const AType &A) const
      { return OType::operator()(x());}

// ******************************************************************
// ******************************************************************
// ******************************************************************

template <typename RType,typename XType1,typename XType2,typename AType,typename OType>
template <typename YType1,typename YType2> 
RType BINARYEXPRESSION<RType,XType1,XType2,AType,OType>::operator()(BINARYEXPRESSION<RType,YType1,YType2,AType,OType>*,const AType &A) const
      { return OType::operator()(x1(A),x2(A));}

template <typename RType,typename XType1,typename XType2,typename AType,typename OType>
template <typename YType2> 
RType BINARYEXPRESSION<RType,XType1,XType2,AType,OType>::operator()(BINARYEXPRESSION<RType,CONSTANT<RType>,YType2,AType,OType>*,const AType &A) const
      { return OType::operator()(x1(),x2(A));}

template <typename RType,typename XType1,typename XType2,typename AType,typename OType>
template <typename YType1> 
RType BINARYEXPRESSION<RType,XType1,XType2,AType,OType>::operator()(BINARYEXPRESSION<RType,YType1,CONSTANT<RType>,AType,OType>*,const AType &A) const
      { return OType::operator()(x1(A),x2());}

template <typename RType,typename XType1,typename XType2,typename AType,typename OType>
RType BINARYEXPRESSION<RType,XType1,XType2,AType,OType>::operator()(BINARYEXPRESSION<RType,CONSTANT<RType>,CONSTANT<RType>,AType,OType>*,const AType &A) const
      { return OType::operator()(x1(),x2());}

#endif

