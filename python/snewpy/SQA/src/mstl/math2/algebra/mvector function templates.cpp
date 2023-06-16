#if !defined(_MVECTOR_FUNCTIONS)
#define _MVECTOR_FUNCTIONS

#include "mstl.h"

template <typename Type,typename VType> 
MVECTOREXPRESSION_NEGATE<Type,VType> operator-(MVECTOREXPRESSION<Type,VType> const &V)
         { return MVECTOREXPRESSION_NEGATE<Type,VType>(V);}

// ******************************************

template <typename Type,typename VType1,typename VType2>
MVECTOREXPRESSION_VECTORADDITION<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> > operator+(MVECTOREXPRESSION<Type,VType1> const &V1,MVECTOREXPRESSION<Type,VType2> const &V2)
         { return MVECTOREXPRESSION_VECTORADDITION<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> >(V1,V2);}

template <typename Type,typename VType1,typename VType2>
MVECTOREXPRESSION_VECTORSUBTRACTION<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> > operator-(MVECTOREXPRESSION<Type,VType1> const &V1,MVECTOREXPRESSION<Type,VType2> const &V2)
         { return MVECTOREXPRESSION_VECTORSUBTRACTION<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> >(V1,V2);}

template <typename Type,typename VType,typename SType>
MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION<Type,MVECTOREXPRESSION<Type,VType>,SType> operator*(MVECTOREXPRESSION<Type,VType> const &V,SType const &S)
         { return MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION<Type,MVECTOREXPRESSION<Type,VType>,SType>(V,S);}

template <typename Type,typename SType,typename VType>
MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION<Type,SType,MVECTOREXPRESSION<Type,VType> > operator*(SType const &S,MVECTOREXPRESSION<Type,VType> const &V)
         { return MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION<Type,SType,MVECTOREXPRESSION<Type,VType> >(S,V);}

template <typename Type,typename VType,typename SType>
MVECTOREXPRESSION_VECTORSCALARDIVISION<Type,MVECTOREXPRESSION<Type,VType>,SType> operator/(MVECTOREXPRESSION<Type,VType> const &V,SType const &S)
         { return MVECTOREXPRESSION_VECTORSCALARDIVISION<Type,MVECTOREXPRESSION<Type,VType>,SType>(V,S);}

template <typename Type,typename VType1,typename VType2>
MVECTOREXPRESSION_CROSSPRODUCT<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> > operator^(MVECTOREXPRESSION<Type,VType1> const &V1,MVECTOREXPRESSION<Type,VType2> const &V2)
         { 
           #ifdef _LAERRORS
           if(V1.Size()!=3 || V2.Size()!=3){ throw NOT_3D("operator^(MVECTOREXPRESSION<Type,VType1>&,MVECTOREXPRESSION<Type,VType2>&)");}
           #endif
           return MVECTOREXPRESSION_CROSSPRODUCT<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> >(V1,V2);
          }

// ******************************************
// dot product
template <typename Type,typename VType1,typename VType2> 
Type operator*(MVECTOREXPRESSION<Type,VType1> const &V1,MVECTOREXPRESSION<Type,VType2> const &V2)
 	      { 
                #ifdef _LAERRORS
                if(V1.Size()!=V2.Size()){ throw DIFFERENT_LENGTHS("operator*(MVECTOREXPRESSION<Type,VType1>&,MVECTOREXPRESSION<Type,VType2>&)");}
                #endif
                Type T=Zero<Type>();
                int a,amax=static_cast<int>(V1.Size())-1; 
                for(a=0;a<=amax;a++){ T+=V1[a]*V2[a];}
	        return T;
	       }

// ******************************************

template <typename Type,typename VType>
std::ostream& operator<<(std::ostream &os,MVECTOREXPRESSION<Type,VType> const &V)
 	      { os.precision(15);
	        int a,amax=static_cast<int>(V.Size())-1; 
                for(a=0;a<=amax;a++){ os<<V[a]<<",\t"<<std::flush;}
	        os<<"\n";
	        return os;
	       }

#endif

