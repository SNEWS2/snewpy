#if !defined(_LINALG_FUNCTIONS)
#define _LINALG_FUNCTIONS

#include "mstl.h"

// *******************************************************
// ************************ FUNCTION DEFINITIONS *********
// *******************************************************

template <typename Type,typename MType> inline MATRIXEXPRESSION_NEGATE<Type,MATRIXEXPRESSION<Type,MType> > operator-(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_NEGATE<Type,MATRIXEXPRESSION<Type,MType> >(M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXTRANSPOSE<Type,MATRIXEXPRESSION<Type,MType> > Transpose(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_MATRIXTRANSPOSE<Type,MATRIXEXPRESSION<Type,MType> >(M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXANTITRANSPOSE<Type,MATRIXEXPRESSION<Type,MType> > AntiTranspose(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_MATRIXANTITRANSPOSE<Type,MATRIXEXPRESSION<Type,MType> >(M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXADJOINT<Type,MATRIXEXPRESSION<Type,MType> > Adjoint(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_MATRIXADJOINT<Type,MATRIXEXPRESSION<Type,MType> >(M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_CONJUGATE<Type,MATRIXEXPRESSION<Type,MType> > Conjugate(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_CONJUGATE<Type,MATRIXEXPRESSION<Type,MType> >(M);}

// *******************************************************************

template <typename Type,typename MType1,typename MType2> inline MATRIXEXPRESSION_MATRIXADDITION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator+(MATRIXEXPRESSION<Type,MType1> const &M1,MATRIXEXPRESSION<Type,MType2> const &M2)
         { return MATRIXEXPRESSION_MATRIXADDITION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> >(M1,M2);}

template <typename Type,typename MType1,typename MType2> inline MATRIXEXPRESSION_MATRIXSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator-(MATRIXEXPRESSION<Type,MType1> const &M1,MATRIXEXPRESSION<Type,MType2> const &M2)
         { return MATRIXEXPRESSION_MATRIXSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> >(M1,M2);}

template <typename Type,typename MType1,typename MType2> inline MATRIXEXPRESSION_MATRIXMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> > operator*(MATRIXEXPRESSION<Type,MType1> const &M1,MATRIXEXPRESSION<Type,MType2> const &M2)
         { return MATRIXEXPRESSION_MATRIXMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<Type,MType2> >(M1,M2);}

template <typename Type,typename MType1,typename MType2> inline MATRIXEXPRESSION_MATRIXMULTIPLICATION<std::complex<Type>,MATRIXEXPRESSION<std::complex<Type>,MType1>,MATRIXEXPRESSION<Type,MType2> > operator*(MATRIXEXPRESSION<std::complex<Type>,MType1> const &M1,MATRIXEXPRESSION<Type,MType2> const &M2)
         { return MATRIXEXPRESSION_MATRIXMULTIPLICATION<std::complex<Type>,MATRIXEXPRESSION<std::complex<Type>,MType1>,MATRIXEXPRESSION<Type,MType2> >(M1,M2);}

template <typename Type,typename MType1,typename MType2> inline MATRIXEXPRESSION_MATRIXMULTIPLICATION<std::complex<Type>,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<std::complex<Type>,MType2> > operator*(MATRIXEXPRESSION<Type,MType1> const &M1,MATRIXEXPRESSION<std::complex<Type>,MType2> const &M2)
         { return MATRIXEXPRESSION_MATRIXMULTIPLICATION<std::complex<Type>,MATRIXEXPRESSION<Type,MType1>,MATRIXEXPRESSION<std::complex<Type>,MType2> >(M1,M2);}

// *******************************************************************

template <typename Type,typename MType> inline MATRIXEXPRESSION_SCALARMATRIXADDITION<Type,MATRIXEXPRESSION<Type,MType> > operator+(Type const &S,MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_SCALARMATRIXADDITION<Type,MATRIXEXPRESSION<Type,MType> >(S,M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXSCALARADDITION<Type,MATRIXEXPRESSION<Type,MType> > operator+(MATRIXEXPRESSION<Type,MType> const &M,Type const &S)
         { return MATRIXEXPRESSION_MATRIXSCALARADDITION<Type,MATRIXEXPRESSION<Type,MType> >(M,S);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType> > operator-(Type const &S,MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType> >(S,M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType> > operator-(MATRIXEXPRESSION<Type,MType> const &M,Type const &S)
         { return MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION<Type,MATRIXEXPRESSION<Type,MType> >(M,S);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> > operator*(Type const &S,MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> >(S,M);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> > operator*(MATRIXEXPRESSION<Type,MType> const &M,Type const &S)
         { return MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType> >(M,S);}

template <typename Type,typename MType> inline MATRIXEXPRESSION_MATRIXSCALARDIVISION<Type,MATRIXEXPRESSION<Type,MType> > operator/(MATRIXEXPRESSION<Type,MType> const &M,Type const &S)
         { return MATRIXEXPRESSION_MATRIXSCALARDIVISION<Type,MATRIXEXPRESSION<Type,MType> >(M,S);}

// *******************************************************************

template <typename Type,typename CVType,typename RVType> inline MATRIXEXPRESSION_VECTORPRODUCT<Type,CVECTOREXPRESSION<Type,CVType>,RVECTOREXPRESSION<Type,RVType> > operator*(CVECTOREXPRESSION<Type,CVType> const &CV,RVECTOREXPRESSION<Type,RVType> const &RV)
         { return MATRIXEXPRESSION_VECTORPRODUCT<Type,CVECTOREXPRESSION<Type,CVType>,RVECTOREXPRESSION<Type,RVType> >(CV,RV);}

// *******************************************************
// *******************************************************
// *******************************************************

template <typename Type,typename CVType> inline CVECTOREXPRESSION_NEGATE<Type,CVECTOREXPRESSION<Type,CVType> > operator-(CVECTOREXPRESSION<Type,CVType> const &CV)
         { return CVECTOREXPRESSION_NEGATE<Type,CVECTOREXPRESSION<Type,CVType> >(CV);}

template <typename Type,typename CVType> inline RVECTOREXPRESSION_CVECTORTRANSPOSE<Type,CVECTOREXPRESSION<Type,CVType> > Transpose(CVECTOREXPRESSION<Type,CVType> const &CV)
         { return RVECTOREXPRESSION_CVECTORTRANSPOSE<Type,CVECTOREXPRESSION<Type,CVType> >(CV);}

template <typename Type,typename CVType> inline RVECTOREXPRESSION_CVECTORADJOINT<Type,CVECTOREXPRESSION<Type,CVType> > Adjoint(CVECTOREXPRESSION<Type,CVType> const &CV)
         { return RVECTOREXPRESSION_CVECTORADJOINT<Type,CVECTOREXPRESSION<Type,CVType> >(CV);}

template <typename Type,typename CVType> inline CVECTOREXPRESSION_CONJUGATE<Type,CVECTOREXPRESSION<Type,CVType> > Conjugate(CVECTOREXPRESSION<Type,CVType> const &CV)
         { return CVECTOREXPRESSION_CONJUGATE<Type,CVECTOREXPRESSION<Type,CVType> >(CV);}

// *******************************************************************

template <typename Type,typename CVType1,typename CVType2> inline CVECTOREXPRESSION_CVECTORADDITION<Type,CVECTOREXPRESSION<Type,CVType1>,CVECTOREXPRESSION<Type,CVType2> > operator+(CVECTOREXPRESSION<Type,CVType1> const &CV1,CVECTOREXPRESSION<Type,CVType2> const &CV2)
         { return CVECTOREXPRESSION_CVECTORADDITION<Type,CVECTOREXPRESSION<Type,CVType1>,CVECTOREXPRESSION<Type,CVType2> >(CV1,CV2);}

template <typename Type,typename CVType1,typename CVType2> inline CVECTOREXPRESSION_CVECTORSUBTRACTION<Type,CVECTOREXPRESSION<Type,CVType1>,CVECTOREXPRESSION<Type,CVType2> > operator-(CVECTOREXPRESSION<Type,CVType1> const &CV1,CVECTOREXPRESSION<Type,CVType2> const &CV2)
         { return CVECTOREXPRESSION_CVECTORSUBTRACTION<Type,CVECTOREXPRESSION<Type,CVType1>,CVECTOREXPRESSION<Type,CVType2> >(CV1,CV2);}

// *******************************************************************

template <typename Type,typename CVType> inline CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION<Type,CVECTOREXPRESSION<Type,CVType> > operator*(CVECTOREXPRESSION<Type,CVType> const &CV,Type const &S)
         { return CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION<Type,CVECTOREXPRESSION<Type,CVType> >(CV,S);}

template <typename Type,typename CVType> inline CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION<Type,CVECTOREXPRESSION<Type,CVType> > operator*(Type const &S,CVECTOREXPRESSION<Type,CVType> const &CV)
         { return CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION<Type,CVECTOREXPRESSION<Type,CVType> >(S,CV);}

template <typename Type,typename CVType> inline CVECTOREXPRESSION_CVECTORSCALARDIVISION<Type,CVECTOREXPRESSION<Type,CVType> > operator/(CVECTOREXPRESSION<Type,CVType> const &CV,Type const &S)
         { return CVECTOREXPRESSION_CVECTORSCALARDIVISION<Type,CVECTOREXPRESSION<Type,CVType> >(CV,S);}

// *******************************************************************

template <typename Type,typename MType,typename CVType> inline CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType>,CVECTOREXPRESSION<Type,CVType> > operator*(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &CV)
         { return CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION<Type,MATRIXEXPRESSION<Type,MType>,CVECTOREXPRESSION<Type,CVType> >(M,CV);}

// *******************************************************
// *******************************************************
// *******************************************************

template <typename Type,typename RVType> inline RVECTOREXPRESSION_NEGATE<Type,RVECTOREXPRESSION<Type,RVType> > operator-(RVECTOREXPRESSION<Type,RVType> const &RV)
         { return RVECTOREXPRESSION_NEGATE<Type,RVECTOREXPRESSION<Type,RVType> >(RV);}

template <typename Type,typename RVType> inline CVECTOREXPRESSION_RVECTORTRANSPOSE<Type,RVECTOREXPRESSION<Type,RVType> > Transpose(RVECTOREXPRESSION<Type,RVType> const &RV)
         { return CVECTOREXPRESSION_RVECTORTRANSPOSE<Type,RVECTOREXPRESSION<Type,RVType> >(RV);}

template <typename Type,typename RVType> inline CVECTOREXPRESSION_RVECTORADJOINT<Type,RVType> Adjoint(RVECTOREXPRESSION<Type,RVECTOREXPRESSION<Type,RVType> > const &RV)
         { return CVECTOREXPRESSION_RVECTORADJOINT<Type,RVECTOREXPRESSION<Type,RVType> >(RV);}

template <typename Type,typename RVType> inline RVECTOREXPRESSION_CONJUGATE<Type,RVECTOREXPRESSION<Type,RVType> > Conjugate(RVECTOREXPRESSION<Type,RVType> const &RV)
         { return CVECTOREXPRESSION_CONJUGATE<Type,RVECTOREXPRESSION<Type,RVType> >(RV);}

// *******************************************************************

template <typename Type,typename RVType1,typename RVType2> inline RVECTOREXPRESSION_RVECTORADDITION<Type,RVECTOREXPRESSION<Type,RVType1>,RVECTOREXPRESSION<Type,RVType2> > operator+(RVECTOREXPRESSION<Type,RVType1> const &RV1,RVECTOREXPRESSION<Type,RVType2> const &RV2)
         { return RVECTOREXPRESSION_RVECTORADDITION<Type,RVECTOREXPRESSION<Type,RVType1>,RVECTOREXPRESSION<Type,RVType2> >(RV1,RV2);}

template <typename Type,typename RVType1,typename RVType2> inline RVECTOREXPRESSION_RVECTORSUBTRACTION<Type,RVECTOREXPRESSION<Type,RVType1>,RVECTOREXPRESSION<Type,RVType2> > operator-(RVECTOREXPRESSION<Type,RVType1> const &RV1,RVECTOREXPRESSION<Type,RVType2> const &RV2)
         { return RVECTOREXPRESSION_RVECTORSUBTRACTION<Type,RVECTOREXPRESSION<Type,RVType1>,RVECTOREXPRESSION<Type,RVType2> >(RV1,RV2);}

// *******************************************************************

template <typename Type,typename RVType> inline RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> > operator*(RVECTOREXPRESSION<Type,RVType> const &RV,Type const &S)
         { return RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> >(RV,S);}

template <typename Type,typename RVType> inline RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> > operator*(Type const &S,RVECTOREXPRESSION<Type,RVType> const &RV)
         { return RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType> >(S,RV);}

template <typename Type,typename RVType> inline RVECTOREXPRESSION_RVECTORSCALARDIVISION<Type,RVECTOREXPRESSION<Type,RVType> > operator/(RVECTOREXPRESSION<Type,RVType> const &RV,Type const &S)
         { return RVECTOREXPRESSION_RVECTORSCALARDIVISION<Type,RVECTOREXPRESSION<Type,RVType> >(RV,S);}

// *******************************************************************

template <typename Type,typename RVType,typename MType> RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType>,MATRIXEXPRESSION<Type,MType> > operator*(RVECTOREXPRESSION<Type,RVType> const &RV,MATRIXEXPRESSION<Type,MType> const &M)
         { return RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION<Type,RVECTOREXPRESSION<Type,RVType>,MATRIXEXPRESSION<Type,MType> >(RV,M);}

template <typename Type,typename RVType,typename CVType> 
Type operator*(RVECTOREXPRESSION<Type,RVType> const &RV,CVECTOREXPRESSION<Type,CVType> const &CV)
         { int i,imax=static_cast<int>(RV.Size())-1;
           Type T=Zero<Type>();
           for(i=0;i<=imax;i++){ T+=RV[i]*CV[i];}
           return T;
          }

// *************************************************************************
// ********************** << operator **************************************
// *************************************************************************

template <typename Type,typename MType> std::ostream& operator<<(std::ostream &os,MATRIXEXPRESSION<Type,MType> const &M) 
      	 { int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ os<<M(a,b)<<"\t";} os<<"\n"<<std::flush; }
           return os;
          }

template <typename Type,typename CVType> std::ostream& operator<<(std::ostream &os,CVECTOREXPRESSION<Type,CVType> const &CV)
         { os<<"\n(";
           int a,amax=static_cast<int>(CV.N())-1;
	   for(a=0;a<=amax;a++){ os<<"\t"<<CV[a]<<"\n";}
           os<<"\t)"<<std::flush;;
	   return os;
	  }

template <typename Type,typename RVType> std::ostream& operator<<(std::ostream &os,RVECTOREXPRESSION<Type,RVType> const &RV)
         { os<<"\n( ";
           int a,amax=static_cast<int>(RV.N())-1;
           for(a=0;a<=amax;a++){ os<<RV[a]<<"\t";}
           os<<")\n"<<std::flush;;
  	   return os;
   	  }

// *************************************************************************
// *************** operations involving submatrices ************************
// *************************************************************************

template <typename Type,typename MType> 
MATRIX<Type,0,0> SubMatrix(MATRIXEXPRESSION<Type,MType> const &M1,int i0,int i1,int j0,int j1)
         { MATRIX<Type,0,0> M2(i1-i0+1,j1-j0+1);
           int i,j;
           for(i=0;i<=i1-i0;i++){ for(j=0;j<=j1-j0;j++){ M2(i,j)=M1(i0+i,j0+j);} }
           return M2;
          }

template <typename Type,typename MType> 
CVECTOR<Type> SubColumnVector(MATRIXEXPRESSION<Type,MType> const &M,int i0,int i1,int j)
         { CVECTOR<Type> CV(i1-i0+1); 
           int i;
           for(i=0;i<=i1-i0;i++){ CV[i]=M(i0+i,j);} 
           return CV;
          }

template <typename Type,typename MType> 
RVECTOR<Type> SubRowVector(MATRIXEXPRESSION<Type,MType> const &M,int i,int j0,int j1)
         { RVECTOR<Type> RV(j1-j0+1); 
           int j;
           for(j=0;j<=j1-j0;j++){ RV[j]=M(i,j0+j);} 
           return RV;
          }

template <typename Type,typename CVType> 
CVECTOR<Type> SubColumnVector(CVECTOREXPRESSION<Type,CVType> const &CV1,int i0,int i1)
         { CVECTOR<Type> CV2(i1-i0+1); 
           int i;
           for(i=0;i<=i1-i0;i++){ CV2[i]=CV1[i0+i];} 
           return CV2;
          }

template <typename Type,typename RVType> 
RVECTOR<Type> SubRowVector(RVECTOREXPRESSION<Type,RVType> const &RV1,int j0,int j1)
         { RVECTOR<Type> RV2(j1-j0+1); 
           int j;
           for(j=0;j<=j1-j0;j++){ RV2[j]=RV1[j0+j];} 
           return RV2;
          }

// ***************************

template <typename Type,typename MType1,typename MType2> 
MATRIXEXPRESSION_SUBMATRIXADDITION<Type,MType1,MType2> SubMatrixAddition(MATRIXEXPRESSION<Type,MType1> const &M1,MATRIXEXPRESSION<Type,MType2> const &M2,int i0,int j0)
         { return MATRIXEXPRESSION_SUBMATRIXADDITION<Type,MType1,MType2>(M1,M2,i0,j0);}

template <typename Type,typename MType1,typename MType2> 
MATRIXEXPRESSION_SUBMATRIXSUBTRACTION<Type,MType1,MType2> SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType1> const &M1,MATRIXEXPRESSION<Type,MType2> const &M2,int i0,int j0)
         { return MATRIXEXPRESSION_SUBMATRIXSUBTRACTION<Type,MType1,MType2>(M1,M2,i0,j0);}


// ***************************

template <typename Type,typename MType,typename CVType> 
MATRIX<Type,0,0> Deflate(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &CV)
         { MATRIX<Type,0,0> H(HouseholderMatrix(CV));
           return SubMatrix(H*M*Adjoint(H),1,M.N1()-1,1,M.N2()-1); 
          }

template <typename Type,typename MType,typename RVType> 
MATRIX<Type,0,0> Deflate(MATRIXEXPRESSION<Type,MType> const &M,RVECTOREXPRESSION<Type,RVType> const &RV)
         { MATRIX<Type,0,0> H(HouseholderMatrix(RV));
           return SubMatrix(H*M*Adjoint(H),1,M.N1()-1,1,M.N2()-1); 
          }

// *************************************************************************
// ********************** SPECIAL MATRICES AND VECTORS *********************
// *************************************************************************

template <typename Type,std::size_t N1,std::size_t N2> MATRIX<Type,N1,N2> ZeroMatrix(void)
         { 
           #ifdef _LAERRORS
           if(N1==0 || N2==0){ throw ZERO_NUMBER("ZeroMatrix");}
           #endif
           return MATRIX<Type,N1,N2>();
          }

template <typename Type> MATRIX<Type,0,0> ZeroMatrix(std::size_t N1,std::size_t N2)
         { 
           #ifdef _LAERRORS
           if(N1==0 || N2==0){ throw ZERO_NUMBER("ZeroMatrix");}
           #endif
           return MATRIX<Type,0,0>(N1,N2);
          }

template <typename Type> MATRIX<Type,0,0> ZeroMatrix(std::size_t N)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("ZeroMatrix");}
           #endif
           return MATRIX<Type,0,0>(N);
          }

// *******

template <typename Type,std::size_t N> MATRIX<Type,N,N> UnitMatrix(void)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("UnitMatrix");}
           #endif
           MATRIX<Type,N,N> I;
           int i,imax=static_cast<int>(N)-1;
           for(i=0;i<=imax;i++){ I(i,i)=One<Type>();}
           return I;
          }

template <typename Type> MATRIX<Type> UnitMatrix(std::size_t N)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("UnitMatrix");}
           #endif
           MATRIX<Type,0,0> I(N,N);
           int i,imax=static_cast<int>(N)-1;
           for(i=0;i<=imax;i++){ I(i,i)=One<Type>();}
           return I;
          }

// *******

template <typename Type,std::size_t N> MATRIX<Type,N,N> PermutationMatrix(int i,int j)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("PermutationMatrix");}
           #endif
           MATRIX<Type,N,N> P(UnitMatrix<Type,N>());
           std::swap(P(i,i),P(i,j));
           std::swap(P(j,j),P(j,i));
           return P;
          }

template <typename Type> MATRIX<Type,0,0> PermutationMatrix(std::size_t N,int i,int j)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("PermutationMatrix");}
           #endif
           MATRIX<Type,0,0> P(UnitMatrix<Type>(N));
           std::swap(P(i,i),P(i,j));
           std::swap(P(j,j),P(j,i));
           return P;
          }

// *******

template <typename Type,std::size_t N> MATRIX<Type> ProjectionMatrix(int i)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("ProjectionMatrix");}
           #endif
           MATRIX<Type,N,N> P;
           P(i,i)=One<Type>();
           return P;
          }

template <typename Type> MATRIX<Type,0,0> ProjectionMatrix(std::size_t N,int i)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("ProjectionMatrix");}
           #endif
           MATRIX<Type,0,0> P(N,N);
           P(i,i)=One<Type>();
           return P;
          }

// *******

template <typename Type,std::size_t N> MATRIX<Type> RightShiftMatrix(void)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("RightShiftMatrix");}
           #endif
           MATRIX<Type,N,N> RS;
           int i,imax=static_cast<int>(N)-2;
           for(i=0;i<=imax;i++){ RS(i,i+1)=One<Type>();}
           return RS;
          }

template <typename Type> MATRIX<Type,0,0> RightShiftMatrix(std::size_t N)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("RightShiftMatrix");}
           #endif
           MATRIX<Type,0,0> RS(N,N);
           int i,imax=static_cast<int>(N)-2;
           for(i=0;i<=imax;i++){ RS(i,i+1)=One<Type>();}
           return RS;
          }

// ****************************************

template <typename Type,typename CVType> MATRIX<Type,0,0> HouseholderMatrix(CVECTOREXPRESSION<Type,CVType> const &CV)
         { MATRIX<Type,0,0> H(UnitMatrix<Type>(CV.N()));
           Type magnitude=sqrt(Transpose(CV)*CV);
           if(Equality(magnitude,Zero<Type>())==false){ H-=Two<Type>()*CV*Transpose(CV)/magnitude;}
           return H;
          }

template <typename Type,typename CVType> MATRIX<std::complex<Type>,0,0> HouseholderMatrix(CVECTOREXPRESSION<std::complex<Type>,CVType> const &CV)
         { MATRIX<std::complex<Type>,0,0> H(UnitMatrix<std::complex<Type> >(CV.N()));
           Type magnitude=sqrt(std::real(Adjoint(CV)*CV));
           if(Equality(magnitude,Zero<Type>())==false){ H-=Two<std::complex<Type> >()*CV*Adjoint(CV)/magnitude;}
           return H;
          }

template <typename Type,typename RVType> MATRIX<Type,0,0> HouseholderMatrix(RVECTOREXPRESSION<Type,RVType> const &RV)
         { MATRIX<Type,0,0> H(UnitMatrix<Type>(RV.N()));
           Type magnitude=sqrt(RV*Transpose(RV));
           if(Equality(magnitude,Zero<Type>())==false){ H-=Two<Type>()*Transpose(RV)*RV/magnitude;}
           return H;
          }

template <typename Type,typename RVType> MATRIX<std::complex<Type>,0,0> HouseholderMatrix(RVECTOREXPRESSION<std::complex<Type>,RVType> const &RV)
         { MATRIX<std::complex<Type>,0,0> H(UnitMatrix<std::complex<Type> >(RV.N()));
           Type magnitude=sqrt(RV*std::real(Adjoint(RV)));
           if(Equality(magnitude,Zero<Type>())==false){ H-=Two<std::complex<Type> >()*Adjoint(RV)*RV/magnitude;}
           return H;
          }

// ****************************************

template <std::size_t N> MATRIX<double,N,N> Givens(int N1,int N2,double ANGLE)
         { if(N1>N-1){ throw OUT_OF_RANGE<int>(N1,0,N-1,"Givens(ANGLE)");}
           if(N2>N-1){ throw OUT_OF_RANGE<int>(N2,0,N-1,"Givens(ANGLE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"Givens(ANGLE)");}

           MATRIX<double,N,N> R(UnitMatrix<double,N>());
           R(N1,N2)=-R(N1,N1)*sin(ANGLE);   R(N1,N1)*=cos(ANGLE);
           R(N2,N1)=R(N2,N2)*sin(ANGLE);    R(N2,N2)*=cos(ANGLE);
           return R;
          }

// *******

template <std::size_t N> MATRIX<double> Givens(int N1,int N2,double COS,double SIN)
         { if(N1>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,N-1,"Givens(COS,SIN)");}
           if(N2>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,N-1,"Givens(COS,SIN)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"Givens(COS,SIN)");}

           MATRIX<double,N,N> R(UnitMatrix<double,N>());
           R(N1,N2)=-R(N1,N1)*SIN;   R(N1,N1)*=COS;
           R(N2,N1)=R(N2,N2)*SIN;    R(N2,N2)*=COS;
           return R;
          }

// *******

template <std::size_t N> MATRIX<std::complex<double> > ComplexGivens(int N1,int N2,double ANGLE,double PHASE)
         { if(N1>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,N-1,"ComplexGivens(ANGLE,PHASE)");}
           if(N2>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,N-1,"ComplexGivens(ANGLE,PHASE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"ComplexGivens(ANGLE,PHASE)");}

           MATRIX<std::complex<double>,N,N> R(UnitMatrix<std::complex<double>,N>());
           R(N1,N2)=-R(N1,N1)*sin(ANGLE)*std::complex<double>( cos(PHASE),sin(PHASE) );   R(N1,N1)*=cos(ANGLE);
           R(N2,N1)=R(N2,N2)*sin(ANGLE)*std::complex<double>( cos(PHASE),-sin(PHASE) );   R(N2,N2)*=cos(ANGLE);
           return R;
          }

// *******

template <std::size_t N> MATRIX<std::complex<double>,N,N> ComplexGivens(int N1,int N2,double COS,double SIN,double PHASE)
         { if(N1>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,N-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(N2>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,N-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"ComplexGivens(COS,SIN,PHASE)");}

           MATRIX<std::complex<double>,N,N> R(UnitMatrix<std::complex<double>,N>());
           R(N1,N2)=-R(N1,N1)*SIN*std::complex<double>( cos(PHASE),sin(PHASE) );   R(N1,N1)*=COS;
           R(N2,N1)=R(N2,N2)*SIN*std::complex<double>( cos(PHASE),-sin(PHASE) );   R(N2,N2)*=COS;
           return R;
          }

// *******

template <std::size_t N> MATRIX<std::complex<double>,N,N> ComplexGivens(int N1,int N2,std::complex<double> ALPHA,std::complex<double> BETA)
         { if(N1>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,N-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(N2>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,N-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"ComplexGivens(COS,SIN,PHASE)");}

           MATRIX<std::complex<double>,N,N> R(UnitMatrix<std::complex<double>,N>());
           R(N1,N2)=-R(N1,N1)*BETA;             R(N1,N1)*=ALPHA;
           R(N2,N1)=R(N2,N2)*std::conj(BETA);   R(N2,N2)*=ALPHA;
           return R;
          }

// ****************************************
// ****************************************
// ****************************************

template <std::size_t N> MATRIX<double,N,N> IntegerMatrix(void) // a diagonal matrix with entries equal to the integers
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("IntegerMatrix");}
          #endif
          MATRIX<double,N,N> I;
          int i,imax=static_cast<int>(N)-1;
          for(i=1;i<=imax;i++){ I[i][i]=i*1.;}
          return I;
         }

template <std::size_t N> MATRIX<double,N,N> FactorialMatrix(void) // a diagonal matrix with entries equal to the factorial numbers
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("FactorialMatrix");}
          #endif
          MATRIX<double,N,N> M(UnitMatrix<double,N>());
          int i,imax=static_cast<int>(N)-1;
          for(i=2;i<=imax;i++){ M[i][i]*=Factorial(i*1.);}
          return M;
         }

// *******
// *******
// *******

template <typename Type> MATRIX<Type,0,0> DiagonalMatrix(std::vector<Type> const &X)
         { MATRIX<Type,0,0> M(X.size(),X.size());
           int i,imax=static_cast<int>(M.N1())-1;
           for(i=0;i<=imax;i++){ M(i,i)=X[i];}
           return M;
          }

template <typename Type> MATRIX<Type,0,0> VandermondeMatrix(std::vector<Type> const &X)
         { MATRIX<Type,0,0> V(X.size(),X.size()); 
           int i,imax=static_cast<int>(V.N1())-1,j,jmax=static_cast<int>(V.N2())-1;
           for(i=0;i<=imax;i++){ V(i,0)=One<Type>(); for(j=1;j<=jmax;j++){ V(i,j)=V(i,j-1)*X[i];} } 
           return V;
          }

// *******

template <typename Type,std::size_t N> MATRIX<Type,N,N> TranslationMatrix(Type const &X)
         { MATRIX<Type,N,N> S;
           int i,imax=static_cast<int>(N)-1,j;
           for(i=0;i<=imax;i++)
              { S(i,i)=One<Type>();
                 for(j=i-1;j>=0;j--){ S(j,i)=-S(j+1,i)*X*(j+1.)/(i*1.-j);}
               }
           return S;
          }

template <typename Type> MATRIX<Type,0,0> TranslationMatrix(std::size_t N,Type const &X)
         { MATRIX<Type,0,0> S(N,N);
           int i,imax=static_cast<int>(N)-1,j;
           for(i=0;i<=imax;i++)
              { S(i,i)=One<Type>();
                 for(j=i-1;j>=0;j--){ S(j,i)=-S(j+1,i)*X*(j+1.)/(i*1.-j);}
               }
           return S;
          }

// *******
// *******
// *******

template <typename Type,std::size_t M> CVECTOR<Type,M> BasisColumnVector(int N)
         { CVECTOR<Type,M> E; 
           E[N]=One<Type>(); 
           return E;
          }

template <typename Type> CVECTOR<Type,0> BasisColumnVector(std::size_t M,int N)
         { CVECTOR<Type,0> E(M); 
           E[N]=One<Type>(); 
           return E;
          }

template <typename Type,std::size_t M> RVECTOR<Type,M> BasisRowVector(int N)
         { RVECTOR<Type,M> E;
           E[N]=One<Type>();
           return E;
          }

template <typename Type> RVECTOR<Type,0> BasisRowVector(std::size_t M,int N)
         { RVECTOR<Type,0> E(M);
           E[N]=One<Type>();
           return E;
          }

// *************************************************************************
// ****************** MATRIX PROPERTY TESTS ********************************
// *************************************************************************

template <typename Type,typename MType> bool ZeroTest(MATRIXEXPRESSION<Type,MType> const &M)
         {
           #ifdef _LAERRORS 
           if(M.N1()*M.N2()==0){ throw EMPTY("ZeroTest");}
           #endif
           bool zero=true; 
           int a=0,amax=static_cast<int>(M.N1())-1, b=0,bmax=static_cast<int>(M.N2())-1;
           while(zero==true && a<=amax)
                { while(zero==true && b<=bmax)
                       { if(Equality(M(a,b),Zero<Type>())==false){ zero=false;}
                         ++b;
                        }
                  ++a; b=0;
                 }
           return zero;
          }

template <typename Type,typename MType> bool DiagonalTest(MATRIXEXPRESSION<Type,MType> const &M)
         { try{ return BanddiagonalTest(M,0,0);}
           catch(EMPTY &E){ E.ChangeFunction("DiagonalTest"); throw E;}
          }

template <typename Type,typename MType> bool TridiagonalTest(MATRIXEXPRESSION<Type,MType> const &M) // return true is the matrix is tridiagonal
         { try{ return BanddiagonalTest(M,1,1);}
           catch(EMPTY &E){ E.ChangeFunction("TridiagonalTest"); throw E;}
          }

template <typename Type,typename MType> bool BanddiagonalTest(MATRIXEXPRESSION<Type,MType> const &M,int p,int q) // return true is the matrix is banddiagonal
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("BanddiagonalTest");}
           #endif
           bool banddiagonal=true; 
           int a=0,amax=static_cast<int>(M.N1())-1, b=0,bmax=static_cast<int>(M.N2())-1,c=a+q+1;

           while(banddiagonal==true && a<=amax)
                { while(banddiagonal==true && b<=std::min(a-p-1,bmax))
                       { if(Equality(M(a,c),Zero<Type>())==false){ banddiagonal=false;}
                         ++b;
                        }
                  while(banddiagonal==true && c<=bmax)
                       { if(Equality(M(a,c),Zero<Type>())==false){ banddiagonal=false;}
                         ++c;
                        }
                  ++a; b=0; c=a+q+1;
                 }
           return banddiagonal;
          }

template <typename Type,typename MType> bool SpurTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("SpurTest");}
           #endif
           bool spur=true; 
           int a=0,amax=static_cast<int>(M.N1())-1;
           while(spur==true && a<=amax)
                { if(Equality(M(a,a),Zero<Type>())==true){ spur=false;}
                  ++a;
                 }
           return spur;
          }

template <typename Type,typename MType> bool SquareTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("SquareTest");}
           #endif
           if(M.N1()==M.N2()){ return true;} else{ return false;}
          }

template <typename Type,typename MType> bool LowerTriangleTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("LowerTriangleTest");}
           #endif
           bool lower=true; 
           int a=0,amax=static_cast<int>(M.N1())-2, b=a+1,bmax=static_cast<int>(M.N2())-1;
           while(lower==true && a<=amax)
                { while(lower==true && b<=bmax)
                       { if(Equality(M(a,b),Zero<Type>())==false){ lower=false;}
                         ++b;
                        }
                  ++a; b=a+1;
                 }
           return lower;
          }

template <typename Type,typename MType> bool UpperTriangleTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("UpperTriangleTest");}
           #endif
           bool upper=true; 
           int a=1,amax=static_cast<int>(M.N1())-1, b=a-1;
           while(upper==true && a<=amax)
                { while(upper==true && b<=a-1)
                       { if(Equality(M(a,b),Zero<Type>())==false){ upper=false;}
                         ++b;
                        }
                  ++a; b=a-1;
                 }
           return upper;
          }

template <typename Type,typename MType> bool SymmetricTest(MATRIXEXPRESSION<Type,MType> const &M) // return true is the matrix is symmetric
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("SymmetricTest");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("SymmetricTest");}
           #endif
           bool symmetric=true; 
           int a=0,amax=static_cast<int>(M.N1())-1, b=a+1,bmax=static_cast<int>(M.N2())-1;
           while(symmetric==true && a<=amax)
                { while(symmetric==true && b<=bmax)
                       { if( M(a,b)!=M(b,a)){ symmetric=false;}
                         ++b;
                        }
                  ++a; b=a+1;
                 }
           return symmetric;
          }

template <typename Type,typename MType> bool AntiSymmetricTest(MATRIXEXPRESSION<Type,MType> const &M) // return true is the matrix is antisymmetric
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("AntiSymmetricTest");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("AntiSymmetricTest");}
           #endif
           bool antisymmetric=true; 
           int a=0,amax=static_cast<int>(M.N1())-1, b=a+1,bmax=static_cast<int>(M.N2())-1;
           while(antisymmetric==true && a<=amax)
                { while(antisymmetric==true && b<=bmax)
                       { if( M(a,b)!=-M(b,a)){ antisymmetric=false;}
                         ++b;
                        }
                  ++a; b=a+1;
                 }
           return antisymmetric;
          }

template <typename Type,typename MType> bool LowerHessenbergTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("LowerHessenbergTest");}
           #endif
           bool lower=true; 
           int a=0,amax=static_cast<int>(M.N1())-2, b=a+2,bmax=static_cast<int>(M.N2())-1;
           while(lower==true && a<=amax)
                { while(lower==true && b<=bmax)
                       { if(Equality(M(a,b),Zero<Type>())==false){ lower=false;}
                         ++b;
                        }
                  ++a; b=a+2;
                 }
           return lower;
          }

template <typename Type,typename MType> bool LowerHessenbergTest(MATRIXEXPRESSION<std::complex<Type>,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("LowerHessenbergTest");}
           #endif
           bool lower=true; 
           int a=0,amax=static_cast<int>(M.N1())-2, b=a+2,bmax=static_cast<int>(M.N2())-1;
           while(lower==true && a<=amax)
                { while(lower==true && b<=bmax)
                       { if(Equality(M(a,b),Zero<std::complex<Type> >())==false){ lower=false;}
                         ++b;
                        }
                  ++a; b=a+2;
                 }
           return lower;
          }

template <typename Type,typename MType> bool UpperHessenbergTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("UpperHessenbergTest");}
           #endif
           bool upper=true; 
           int a=2,amax=static_cast<int>(M.N1())-1, b=a-2;
           while(upper==true && a<=amax)
                { while(upper==true && b<=a-1)
                       { if(Equality(M(a,b),Zero<Type>())==false){ upper=false;}
                         ++b;
                        }
                  ++a; b=a-2;
                 }
           return upper;
          }

template <typename Type,typename MType> bool UpperHessenbergTest(MATRIXEXPRESSION<std::complex<Type>,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("UpperHessenbergTest");}
           #endif
           bool upper=true; 
           int a=2,amax=static_cast<int>(M.N1())-1, b=a-2;
           while(upper==true && a<=amax)
                { while(upper==true && b<=a-1)
                       { if(Equality(M(a,b),Zero<std::complex<Type> >())==false){ upper=false;}
                         ++b;
                        }
                  ++a; b=a-2;
                 }
           return upper;
          }

template <typename Type,typename MType> bool RowDiagonalDominanceTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("RowDiagonalDominanceTest");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("RowDiagonalDominanceTest");}
           #endif
           bool rowdiagonaldominance=true; 
           int a=0,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           while(rowdiagonaldominance==true && a<=amax)
                { Type sum=Zero<Type>(); 
                  for(b=0;b<=a-1;b++){ sum+=abs(M(a,b));} for(b=a+1;b<=bmax;b++){ sum+=abs(M(a,b));}
                  if(sum>abs(M(a,a))){ rowdiagonaldominance=false;} 
                  ++a; 
                 }
           return rowdiagonaldominance;
          }

template <typename Type,typename MType> bool RowDiagonalDominanceTest(MATRIXEXPRESSION<std::complex<Type>,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("RowDiagonalDominanceTest");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("RowDiagonalDominanceTest");}
           #endif
           bool rowdiagonaldominance=true; 
           int a=0,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           while(rowdiagonaldominance==true && a<=amax)
                { Type sum=Zero<Type>(); 
                  for(b=0;b<=a-1;b++){ sum+=abs(M(a,b));} for(b=a+1;b<=bmax;b++){ sum+=abs(M(a,b));}
                  if(sum>abs(M(a,a))){ rowdiagonaldominance=false;} 
                  ++a; 
                 }
           return rowdiagonaldominance;
          }


template <typename Type,typename MType> bool ColumnDiagonalDominanceTest(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("ColumnDiagonalDominanceTest");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("ColumnDiagonalDominanceTest");}
           #endif
           bool columndiagonaldominance=true; 
           int a,amax=static_cast<int>(M.N1())-1, b=0,bmax=static_cast<int>(M.N2())-1;
           while(columndiagonaldominance==true && b<=bmax)
                { Type sum=Zero<Type>(); 
                  for(a=0;a<=b-1;a++){ sum+=abs(M(a,b));} for(a=b+1;a<=amax;a++){ sum+=abs(M(a,b));}
                  if(sum>abs(M(b,b))){ columndiagonaldominance=false;} 
                  ++b; 
                 }
           return columndiagonaldominance;
          }

template <typename Type,typename MType> bool ColumnDiagonalDominanceTest(MATRIXEXPRESSION<std::complex<Type>,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("ColumnDiagonalDominanceTest");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("ColumnDiagonalDominanceTest");}
           #endif
           bool columndiagonaldominance=true; 
           int a,amax=static_cast<int>(M.N1())-1, b=0,bmax=static_cast<int>(M.N2())-1;
           while(columndiagonaldominance==true && b<=bmax)
                { Type sum=Zero<Type>(); 
                  for(a=0;a<=b-1;a++){ sum+=abs(M(a,b));} for(a=b+1;a<=amax;a++){ sum+=abs(M(a,b));}
                  if(sum>abs(M(b,b))){ columndiagonaldominance=false;} 
                  ++b; 
                 }
           return columndiagonaldominance;
          }

// *************************************************************************
// ****************** MATRIX PROPERTIES ************************************
// *************************************************************************

template <typename Type,typename MType> Type FrobeniusNorm(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("FrobeniusNorm","MATRIX");}
           #endif
           MATRIX<Type,0,0> MM(M.N1(),M.N2());
           int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           Type FN=Zero<Type>();
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ FN+=M(a,b)*M(a,b);} }

           return sqrt(FN);
          }

template <typename Type,typename MType> Type FrobeniusNorm(MATRIXEXPRESSION<std::complex<Type>,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("FrobeniusNorm","MATRIX");}
           #endif
           MATRIX<Type,0,0> MM(M.N1(),M.N2());
           int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           Type FN=Zero<Type>();
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ FN+=M(a,b)*conj(M(a,b));} }

           return sqrt(FN);
          }

template <typename Type,typename MType> MATRIX<Type,0,0> Diagonal(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("Diagonal","MATRIX");}
           #endif
           MATRIX<Type,0,0> MM(M.N1(),M.N2());
           int a,amax=std::min(static_cast<int>(M.N1())-1,static_cast<int>(M.N2())-1);
           for(a=0;a<=amax;a++){ MM(a,a)=M(a,a);}

           return MM;
          }

template <typename Type,typename MType> MATRIX<Type,0,0> LowerTriangle(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("LowerTriangle","MATRIX");}
           #endif
           MATRIX<Type,0,0> MM(M.N1(),M.N2());
           int a,amax=static_cast<int>(M.N1())-1,b;
           for(a=0;a<=amax;a++){ for(b=0;b<=a;b++){ MM(a,b)=M(a,b);} }
           return MM;
          }

template <typename Type,typename MType> MATRIX<Type,0,0> UpperTriangle(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("UpperTriangle","MATRIX");}
           #endif
           MATRIX<Type,0,0> MM(M.N1(),M.N2());
           int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           for(a=0;a<=amax;a++){ for(b=a;b<=bmax;b++){ MM(a,b)=M(a,b);} }
           return MM;
          }

template <typename Type,typename MType> CVECTOR<Type> Column(MATRIXEXPRESSION<Type,MType> const &M,int j)
         { CVECTOR<Type> C(M.N1());
           int i,imax=static_cast<int>(M.N1())-1;
           for(i=0;i<=imax;i++){ C[i]=M(i,j);}
           return C;
          }

template <typename Type,typename MType> RVECTOR<Type> Row(MATRIXEXPRESSION<Type,MType> const &M,int i)
         { RVECTOR<Type> R(M.N2());
           int j, jmax=static_cast<int>(M.N2())-1;
           for(j=0;j<=jmax;j++){ R[j]=M(i,j);}
           return R;
          }

template <typename Type,typename MType> Type Trace(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("Trace","MATRIX");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("Trace");}
           #endif
           Type Tr=Zero<Type>();
           int a,amax=static_cast<int>(M.N1())-1;
           for(a=0;a<=amax;a++){ Tr+=M(a,a);}
           return Tr;
          }

template <typename Type,typename MType> Type Cofactor(MATRIXEXPRESSION<Type,MType> const &M,int i,int j)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("Cofactor","MATRIX");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("Cofactor");}
           #endif

           int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
           MATRIX<Type,0,0> MM(M.N1()-1,M.N2()-1);

                   for(a=0;a<=i-1;a++)
                      { for(b=0;b<=j-1;b++){ MM(a,b)=M(a,b);}
                        for(b=j+1;b<=bmax;b++){ MM(a,b-1)=M(a,b);}
                       }

                   for(a=i+1;a<=amax;a++)
                      { for(b=0;b<=j-1;b++){ MM(a-1,b)=M(a,b);}
                        for(b=j+1;b<=bmax;b++){ MM(a-1,b-1)=M(a,b);}
                       }
           return Determinant(MM)*pow(-One<Type>(),(double)(i+j));
          }

template <typename Type,typename MType> Type Determinant(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("Determinant","MATRIX");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("Determinant");}
           #endif

           if(M.N1()==1){ return M(0,0);}

           int b,bmax=static_cast<int>(M.N2())-1;
           Type T=Zero<Type>();
           for(b=0;b<=bmax;b++){ T+=M(0,b)*Cofactor(M,0,b);}

           return T;
          }

template <typename Type,typename MType> Type Spur(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("Spur","MATRIX");}
           #endif
           Type T=M(0,0);
           int a,amax=std::min(static_cast<int>(M.N1())-1,static_cast<int>(M.N2())-1);
           for(a=1;a<=amax;a++){ T*=M(a,a);}
           return T;
          }

// *************************************************************************
// ************************** MATRIX INVERSIONS ****************************
// *************************************************************************

template <typename Type,typename MType> 
MATRIX<Type,0,0> Invert(MATRIXEXPRESSION<Type,MType> const &M)
         { try{ return Inverse(M);} catch(SINGULAR &S){ S.ChangeFunction("Invert"); throw S;} }


template <typename Type,typename MType> MATRIX<Type,0,0> Inverse(MATRIXEXPRESSION<Type,MType> const &M)
         { if(M.N1()!=M.N2()){ try{ return MPInverse(M);} catch(...){;} }
           else{ if(M.N1()==2){ try{ return LInverse(M);} catch(...){;} }
                 if(SpurTest(M)==true){ if(DiagonalTest(M)==true){ try{ return DInverse(M);} catch(...){;} }
                                        if(LowerTriangleTest(M)==true){ try{ return LTInverse(M);} catch(...){;} }
                                        if(UpperTriangleTest(M)==true){ try{ return UTInverse(M);} catch(...){;} }
                                        try{ return LUInverse(M);} catch(...){;} }
                                       }
           try{ return GJInverse(M);} catch(...){;}
           try{ return LInverse(M);} catch(SINGULAR &S){ S.ChangeFunction("Inverse"); throw S;}
          }

template <typename Type,typename MType> MATRIX<Type,0,0> DInverse(MATRIXEXPRESSION<Type,MType> const &M)
         { try{ return MATRIX<Type,0,0>(M).DInvert();}
           catch(SINGULAR &S){ S.ChangeFunction("DInverse"); throw S;}
          }

template <typename Type,typename MType> MATRIX<Type,0,0> GJInverse(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIX<Type,0,0>(M).GJInvert();}

template <typename Type,typename MType> MATRIX<Type,0,0> LInverse(MATRIXEXPRESSION<Type,MType> const &M,bool message) // Laplace inversion by Cofactors/Determinant
         { try{ return MATRIX<Type,0,0>(M).LInvert(message);}
           catch(SINGULAR &S){ S.ChangeFunction("LInverse"); throw S;}
          }

template <typename Type,typename MType> MATRIX<Type,0,0> LTInverse(MATRIXEXPRESSION<Type,MType> const &M) // Inversion of a Lower Triangular matrix
         { try{ return MATRIX<Type,0,0>(M).LTInvert();}
           catch(SINGULAR &S){ S.ChangeFunction("LTInverse"); throw S;}
          }

template <typename Type,typename MType> MATRIX<Type,0,0> UTInverse(MATRIXEXPRESSION<Type,MType> const &MM) // Inversion of an Upper Triangular matrix
         { return MATRIX<Type>(MM).UTInvert();}

template <typename Type,typename MType> MATRIX<Type,0,0> LUInverse(MATRIXEXPRESSION<Type,MType> const &MM)
         { return MATRIX<Type,0,0>(MM).LUInvert();}

template <typename Type,typename MType> MATRIX<Type,0,0> MPInverse(MATRIXEXPRESSION<Type,MType> const &MM)
         { return MATRIX<Type,0,0>(MM).MPInvert();}

// *************************************************************************
// ****************** SPECIAL MATRIX INVERSIONS ****************************
// *************************************************************************

template <typename Type,typename MType> MATRIX<Type,0,0> TridiagonalInverse(MATRIXEXPRESSION<Type,MType> const &M)
         { return MATRIX<Type,0,0>(M).TridiagonalInvert();}

template <typename Type,typename MType> MATRIX<Type,0,0> VandermondeMatrixInverse(MATRIXEXPRESSION<Type,MType> const &V)
         { if(V.N1()!=V.N2()){ throw NOT_SQUARE("VandermondeMatrixInverse");}
           return VandermondeMatrixInverse(Column(V,1));
          }

template <typename Type> MATRIX<Type,0,0> VandermondeMatrixInverse(std::vector<Type> const &Xoriginal)
         { return VandermondeMatrix(Xoriginal).LUInvert();}

// *************************************************************************
// ****************** MATRIX DECOMPOSITIONS ********************************
// *************************************************************************

template <typename MType> std::vector<MATRIX<double,0,0> > LUDecomposition(MATRIXEXPRESSION<double,MType> const &MM)
         { MATRIX<double,0,0> M(MM); 
           std::vector<MATRIX<double,0,0> > LU(2,MATRIX<double,0,0>(M.N1(),M.N2()));
           double sum;
           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1,k;
           for(i=0;i<=imax;i++){ LU[0](i,i)=1.;} 

           for(j=0;j<=jmax;j++)
              { for(i=0;i<=j;i++)
                   { sum=0.;
                     for(k=0;k<=i-1;k++){ sum+=LU[0](i,k)*LU[1](k,j);}
                     LU[1](i,j)=M(i,j)-sum;
                    }
                for(i=j+1;i<=imax;i++)
                   { sum=0.;
                     for(k=0;k<=j-1;k++){ sum+=LU[0](i,k)*LU[1](k,j);}
                     if(Equality(LU[1](j,j),0.)==true){ throw DIVISION_BY_ZERO("LUDecomposition","MATRIX");}
                     LU[0](i,j)=( M(i,j)-sum )/LU[1](j,j);
                    } 
               }

           return LU;
          }

template <typename MType> std::vector<MATRIX<std::complex<double>,0,0> > LUDecomposition(MATRIXEXPRESSION<std::complex<double>,MType> const &MM)
         { MATRIX<std::complex<double>,0,0> M(MM); 
           std::vector<MATRIX<std::complex<double>,0,0> > LU(2,MATRIX<std::complex<double>,0,0>(M.N1(),M.N2()));
           std::complex<double> sum;
           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1,k;
           for(i=0;i<=imax;i++){ LU[0](i,i)=1.;} 

           for(j=0;j<=jmax;j++)
              { for(i=0;i<=j;i++)
                   { sum=0.;
                     for(k=0;k<=i-1;k++){ sum+=LU[0](i,k)*LU[1](k,j);}
                     LU[1](i,j)=M(i,j)-sum;
                    }
                for(i=j+1;i<=imax;i++)
                   { sum=0.;
                     for(k=0;k<=j-1;k++){ sum+=LU[0](i,k)*LU[1](k,j);}
                     if(Equality(std::real(LU[1](j,j)),0.)==true && Equality(std::imag(LU[1](j,j)),0.)==true){ throw DIVISION_BY_ZERO("LUDecomposition","MATRIX");}
                     LU[0](i,j)=( M(i,j)-sum )/LU[1](j,j);
                    } 
               }

           return LU;
          }

// ******************************************

template <typename MType> 
std::vector<MATRIX<double,0,0> > QRDecomposition(MATRIXEXPRESSION<double,MType> const &MM)
         { MATRIX<double,0,0> M(MM);
           MATRIX<double,0,0> Q=UnitMatrix<double>(M.N1()), Qi, R=M;
           CVECTOR<double,0> u(M.N1());

           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=std::min(imax-1,jmax);i++)
              { for(j=i;j<=imax;j++){ u[j]=R[j][i];}
                u-=CVECTOR<double,0>( sqrt(Transpose(u)*u)*BasisColumnVector<double>(M.N1(),i) ); 
                Qi=UnitMatrix<double>(M.N1()) - 2.*u*Transpose(u)/(Transpose(u)*u);            
                Q*=Transpose(Qi);
                R=MATRIX<double,0,0>(Qi*R);
                u[i]=0.; 
               } 

           std::vector<MATRIX<double,0,0> > QR(2); QR[0]=Q; QR[1]=R;
           return QR;
          }

template <typename MType> 
std::vector<MATRIX<std::complex<double>,0,0> > QRDecomposition(MATRIXEXPRESSION<std::complex<double>,MType> const &MM)
         { MATRIX<std::complex<double>,0,0> M(MM);
           MATRIX<std::complex<double>,0,0> Q=UnitMatrix<std::complex<double> >(M.N1()),Qi, R=M;
           CVECTOR<std::complex<double>,0> x(M.N1()),u(M.N1());

           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=std::min(imax-1,jmax);i++)
              { for(j=i;j<=imax;j++){ x[j]=R[j][i];}
                u=x - CVECTOR<std::complex<double>,0>( std::complex<double>(sqrt(std::real(Adjoint(x)*x)))*BasisColumnVector<std::complex<double> >(M.N1(),i) );
                Qi=UnitMatrix<std::complex<double> >(M.N1()) - (1.+(Adjoint(x)*u)/(Adjoint(u)*x))*u*Adjoint(u)/(Adjoint(u)*u);
                Q*=Adjoint(Qi);
                R=MATRIX<std::complex<double>,0,0>(Qi*R);
                x[i]=std::complex<double>(0.);
               }

           std::vector<MATRIX<std::complex<double>,0,0> > QR(2); QR[0]=Q; QR[1]=R;
           return QR;
          }

// ******************************************

// this is SVD decompostion using Jacobi rotations. 
template <typename MType> 
std::vector<MATRIX<double,0,0> > SVDecomposition(MATRIXEXPRESSION<double,MType> const &M)
         { MATRIX<double,0,0> B(M), U(M.N1(),M.N2()), W(M.N2(),M.N2()), V(UnitMatrix<double>(M.N2())), R;

           CVECTOR<double,0> bi,bj;
           double cosphi,sinphi,p,q,v; 
           bool finish;

           do{ finish=true; 
               int i,imax=static_cast<int>(M.N2())-2, j,jmax=static_cast<int>(M.N2())-1;
               for(i=0;i<=imax;i++)
                  { for(j=i+1;j<=jmax;j++)
                       { bi=Column(B,i);
                         bj=Column(B,j);
     
                         p=Transpose(bi)*bj;
                         q=Transpose(bj)*bj-Transpose(bi)*bi;
                         v=hypot(2.*p,q);

                         if(fabs(p)>=std::numeric_limits<double>::epsilon())
                           { finish=false;
                             if(q>0.){ cosphi=sqrt(1.+q/v)/M_SQRT2; sinphi=-p/v/cosphi;} 
                             else{ sinphi=-sqrt(1.-q/v)/M_SQRT2; cosphi=p/v/sinphi;}
                             R=Givens(M.N2(),i,j,cosphi,sinphi);

                             // this matrix multiplications could be speeded up
                             B*=R; V*=R; 
                            } 
                        }
                   }
              }while(finish==false);

           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1;
           for(j=0;j<=jmax;j++)
              { bj=Column(B,j);
                W[j][j]=sqrt(Transpose(bj)*bj);
                if(Equality(W[j][j],0.)==false){ for(i=0;i<=imax;i++){ U[i][j]=bj[i]/W[j][j];} }
               }

           std::vector<MATRIX<double,0,0> > UWVT(3); UWVT[0]=U; UWVT[1]=W; UWVT[2]=Transpose(V);

           return UWVT;
          }

template <typename MType> 
std::vector<MATRIX<std::complex<double>,0,0> > SVDecomposition(MATRIXEXPRESSION<std::complex<double>,MType> const &M)
         { MATRIX<std::complex<double>,0,0> B=M, U(M.N1(),M.N2()), W(M.N2(),M.N2()), V=UnitMatrix<std::complex<double> >(M.N2()), R;

           CVECTOR<std::complex<double>,0> bi,bj;
           double p,q,v, cosphi, sinphi;

           bool finish;

           do{ finish=true; 
               int i,imax=static_cast<int>(M.N2())-2, j,jmax=static_cast<int>(M.N2())-1;
               for(i=0;i<=imax;i++)
                  { for(j=i+1;j<=jmax;j++)
                       { bi=Column(B,i);
                         bj=Column(B,j);
     
                         p=std::real(Adjoint(bi)*bj);
                         q=std::real(Adjoint(bj)*bj-Adjoint(bi)*bi);
                         v=hypot(2.*p,q);

                         if(fabs(p)>=std::numeric_limits<double>::epsilon())
                           { finish=false;
                             if(q>0.){ cosphi=::sqrt(1.+q/v)/M_SQRT2; sinphi=-p/v/cosphi;} else{ sinphi=-sqrt(1.-q/v)/M_SQRT2; cosphi=p/v/sinphi;}
                             R=Givens(M.N2(),i,j,std::complex<double>(cosphi),std::complex<double>(sinphi));

                             // these matrix multiplications could be sped up given structure of R
                             B*=R; V*=R; 
                            } 
                        }
                   }
              }while(finish==false);

           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1;
           for(j=0;j<=jmax;j++)
              { bj=Column(B,j);
                W[j][j]=sqrt(Adjoint(bj)*bj);
                if(Equality(std::real(W[j][j]),0.)==false && Equality(std::imag(W[j][j]),0.)==false){ for(i=0;i<=imax;i++){ U[i][j]=bj[i]/W[j][j];} }
               }

           std::vector<MATRIX<std::complex<double>,0,0> > UWVT(3); UWVT[0]=U; UWVT[1]=W; UWVT[2]=Transpose(V);

           return UWVT;
          }

// ******************************************

template <typename Type,typename MType> 
MATRIX<Type,0,0> Diagonalize(MATRIXEXPRESSION<Type,MType> const &M) 
         { 
           #ifdef _LAERRORS
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("Diagonalize");}
           #endif
           return SVDecomposition(M)[1];
          }

// *************************************************************************

template <typename MType> 
std::vector<std::complex<double> > QREigenValues(MATRIXEXPRESSION<double,MType> const &ME,bool deflate) 
         { 
           #ifdef _LAERRORS
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("QREigenValues");}
           #endif

           MATRIX<double,0,0> M(UpperHessenberg(ME));
           MATRIX<double,0,0> P1(3,3);
           CVECTOR<double,0> x(3),u(3);

           int b,i,j,k,l,m,n=static_cast<int>(M.N1()), iterations;
           double a1,a0,C, mi0,mi1,mi2, m0j,m1j,m2j;
           std::vector<std::complex<double> > eigenvalues(M.N1()), eigenvaluepair;
           bool rootfound,smallsubdiagonal;

           iterations=0;
           m=l=0;

           if(deflate==true)
             { smallsubdiagonal=true;          
               while(m<=n-2 && smallsubdiagonal==true)
                    { if( abs(M[m+1][m]) <= 4.*std::numeric_limits<double>::epsilon()*abs(M[m][m]+M[m+1][m+1]) && abs(M[m+1][m]*M[m][m+1]) <= 4.*std::numeric_limits<double>::epsilon()*abs(M[m][m]*M[m+1][m+1]) )
                        { eigenvalues[m]=std::complex<double>(M[m][m]);
                          M[m+1][m]=0.;
                          ++m;
                        }
                      else{ smallsubdiagonal=false;}
                     };
              }

           while(n-m>=3)
                { rootfound=false;    
                  if( abs(M[n-1][n-2]) <= 4.*std::numeric_limits<double>::epsilon()*(abs(M[n-1][n-1])+abs(M[n-2][n-2])) )      
                    { eigenvalues[n-1]=std::complex<double>(M[n-1][n-1]);
                      M[n-1][n-2]=0.;
                      n-=1; 
                      rootfound=true;
                     }
                  else{ if( abs(M[n-2][n-3]) <= 4.*std::numeric_limits<double>::epsilon()*(abs(M[n-2][n-2])+abs(M[n-3][n-3])) )
                          { a1=M[n-1][n-1]+M[n-2][n-2];
                            a0=M[n-1][n-1]*M[n-2][n-2]-M[n-1][n-2]*M[n-2][n-1];
                            eigenvaluepair=QuadraticRoots(1.,-a1,a0);
                            eigenvalues[n-1]=eigenvaluepair[0];
                            eigenvalues[n-2]=eigenvaluepair[1];
                            M[n-2][n-3]=0.;
                            n-=2; 
                            rootfound=true;
                           }
                       }

                  if(rootfound==false)
                    { l=m;
                      if(deflate==true)
                        { for(b=n-3;b>=m;b--)
                             { if(b>0 && abs(M[b][b-1]*M[b-1][b]) <= 16.*std::numeric_limits<double>::epsilon()*(abs(M[b][b]*M[b-1][b-1])) && abs(M[b][b-1]) <= 4.*std::numeric_limits<double>::epsilon()*abs(M[b][b]+M[b-1][b-1]) )
                                 { l=std::max(l,b);}
                              }
                         }//end of if(deflate==true) loop

                      if(iterations==20){ a1=1.5*(abs(M[n-1][n-2])+abs(M[n-2][n-3])); a0=pow(abs(M[n-1][n-2])+abs(M[n-2][n-3]),2.); iterations=0;} 
                      else{ a1=M[n-1][n-1]+M[n-2][n-2]; a0=M[n-1][n-1]*M[n-2][n-2]-M[n-1][n-2]*M[n-2][n-1];}

                      x[0]=u[0]=M[l][l]*(M[l][l]-a1) + M[l][l+1]*M[l+1][l] + a0;
                      x[1]=u[1]=M[l+1][l]*(M[l][l]+M[l+1][l+1]-a1);
                      x[2]=u[2]=M[l+1][l]*M[l+2][l+1];
                      u[0]+= Sign(x[0])*sqrt(pow(x[0],2.)+pow(x[1],2.)+pow(x[2],2.));

                      P1=UnitMatrix<double>(3) - ( 1.+(Transpose(x)*u)/(Transpose(u)*x) )*u*Transpose(u)/(Transpose(u)*u);

                      // multiply from the left by (P1,I)                      
                      for(j=l;j<=n-1;j++)
                              { m0j=P1[0][0]*M[l][j] + P1[0][1]*M[l+1][j] + P1[0][2]*M[l+2][j];
                                m1j=P1[1][0]*M[l][j] + P1[1][1]*M[l+1][j] + P1[1][2]*M[l+2][j];
                                m2j=P1[2][0]*M[l][j] + P1[2][1]*M[l+1][j] + P1[2][2]*M[l+2][j];
                                M[l][j]=m0j; 
                                M[l+1][j]=m1j; 
                                M[l+2][j]=m2j;
                               } 

                      // multiply from the right by (P1,I)                      
                      for(i=l;i<=n-1;i++)
                              { mi0=M[i][l]*P1[0][0] + M[i][l+1]*P1[1][0] + M[i][l+2]*P1[2][0];
                                mi1=M[i][l]*P1[0][1] + M[i][l+1]*P1[1][1] + M[i][l+2]*P1[2][1];
                                mi2=M[i][l]*P1[0][2] + M[i][l+1]*P1[1][2] + M[i][l+2]*P1[2][2];
                                M[i][l]=mi0; 
                                M[i][l+1]=mi1; 
                                M[i][l+2]=mi2;
                               } 

                      // chase the bulge - same as in UpperHessenberg function
                      if(l<=n-4)
                        { for(j=l;j<=n-4;j++)//sweep through columns up to four from last
                             { i=j+2;
                               if(abs(M[j+1][j])<abs(M[i][j])){ M.SwapRows(j+1,i); M.SwapColumns(j+1,i);}
                               if(abs(M[j+1][j])<abs(M[i+1][j])){ M.SwapRows(j+1,i+1); M.SwapColumns(j+1,i+1);}
                               if(M[j+1][j]!=0.)
                                 { C=M[i][j]/M[j+1][j];
                                   M[i][j]=0.;
                                   for(k=n-1;k>=j+1;k--){ M[i][k]-=M[j+1][k]*C;}  // subtract rows
                                   for(k=n-1;k>=0;k--){ M[k][j+1]+=M[k][i]*C;}    // add columns

                                   i=j+3;
                                   C=M[i][j]/M[j+1][j];
                                   M[i][j]=0.;
                                   for(k=n-1;k>=j+1;k--){ M[i][k]-=M[j+1][k]*C;}  // subtract rows
                                   for(k=n-1;k>=0;k--){ M[k][j+1]+=M[k][i]*C;}    // add columns
                                  }
                              }
                         }

                      // deal with the column three from the end
                      j=n-3;
                      i=j+2;
                      if(abs(M[j+1][j])<abs(M[i][j])){ M.SwapRows(j+1,i); M.SwapColumns(j+1,i);}
                      if(M[j+1][j]!=0.)
                        { C=M[i][j]/M[j+1][j];
                          M[i][j]=0.;
                          for(k=n-1;k>=j+1;k--){ M[i][k]-=M[j+1][k]*C;}  // subtract rows
                          for(k=n-1;k>=0;k--){ M[k][j+1]+=M[k][i]*C;}    // add columns
                         }

                      iterations++; 
                      //if(iterations==40){ throw NO_SOLUTION("QREigenValues");} 

                     }//end of if(rootfound==false) loop
                  else{ iterations=0;}
 
                 };//end of while(n-m>=3) loop

           if(n-m==2){ a1=M[m+1][m+1]+M[m][m];
                       a0=M[m+1][m+1]*M[m][m]-M[m+1][m]*M[m][m+1];
                       eigenvaluepair=QuadraticRoots(1.,-a1,a0);
                       eigenvalues[m]=eigenvaluepair[0];
                       eigenvalues[m+1]=eigenvaluepair[1];
                      }
           if(n-m==1){ eigenvalues[m]=std::complex<double>(M[m][m]);}
           
           return eigenvalues;
          }

template <typename MType> 
std::vector<std::complex<double> > QREigenValues(MATRIXEXPRESSION<std::complex<double>,MType> const &ME,bool deflate) 
         { 
           #ifdef _LAERRORS
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("QREigenValues");}
           #endif

           MATRIX<std::complex<double>,0,0> M(UpperHessenberg(ME));
           MATRIX<std::complex<double>,0,0> P1(3,3);
           CVECTOR<std::complex<double>,0> x(3),u(3);

           int b,i,j,k,l,m,n=static_cast<int>(M.N1()), iterations;
           std::complex<double> a1,a0,C, mi0,mi1,mi2, m0j,m1j,m2j;
           std::vector<std::complex<double> > eigenvalues(M.N1()), eigenvaluepair;
           bool rootfound,smallsubdiagonal;

           iterations=0;
           m=l=0;

           if(deflate==true)
             { smallsubdiagonal=true;          
               while(m<=n-2 && smallsubdiagonal==true)
                    { if( abs(M[m+1][m]) <= 4.*std::numeric_limits<double>::epsilon()*abs(M[m][m]+M[m+1][m+1]) && abs(M[m+1][m]*M[m][m+1]) <= 4.*std::numeric_limits<double>::epsilon()*abs(M[m][m]*M[m+1][m+1]) )
                        { eigenvalues[m]=M[m][m];
                          M[m+1][m]=0.;
                          ++m;
                        }
                      else{ smallsubdiagonal=false;}
                     };
              }

           while(n-m>=3)
                { rootfound=false;    
                  if( abs(M[n-1][n-2]) <= 4.*std::numeric_limits<double>::epsilon()*(abs(M[n-1][n-1])+abs(M[n-2][n-2])) )      
                    { eigenvalues[n-1]=M[n-1][n-1];
                      M[n-1][n-2]=0.;
                      n-=1; 
                      rootfound=true;
                     }
                  else{ if( abs(M[n-2][n-3]) <= 4.*std::numeric_limits<double>::epsilon()*(abs(M[n-2][n-2])+abs(M[n-3][n-3])) )
                          { a1=M[n-1][n-1]+M[n-2][n-2];
                            a0=M[n-1][n-1]*M[n-2][n-2]-M[n-1][n-2]*M[n-2][n-1];
                            eigenvaluepair=QuadraticRoots(1.,-a1,a0);
                            eigenvalues[n-1]=eigenvaluepair[0];
                            eigenvalues[n-2]=eigenvaluepair[1];
                            M[n-2][n-3]=0.;
                            n-=2; 
                            rootfound=true;
                           }
                       }

                  if(rootfound==false)
                    { l=m;
                      if(deflate==true)
                        { for(b=n-3;b>=m;b--)
                             { if(b>0 && abs(M[b][b-1]*M[b-1][b]) <= 4.*std::numeric_limits<double>::epsilon()*(abs(M[b][b]*M[b-1][b-1])) && abs(M[b][b-1]) <= 4.*std::numeric_limits<double>::epsilon()*abs(M[b][b]+M[b-1][b-1]) )
                                 { l=std::max(l,b);}
                              }
                         }//end of if(deflate==true) loop

                      if(iterations==20){ a1=1.5*(abs(M[n-1][n-2])+abs(M[n-2][n-3])); a0=pow(abs(M[n-1][n-2])+abs(M[n-2][n-3]),2.); iterations=0;} 
                      else{ a1=M[n-1][n-1]+M[n-2][n-2]; a0=M[n-1][n-1]*M[n-2][n-2]-M[n-1][n-2]*M[n-2][n-1];}

                      x[0]=u[0]=M[l][l]*(M[l][l]-a1) + M[l][l+1]*M[l+1][l] + a0;
                      x[1]=u[1]=M[l+1][l]*(M[l][l]+M[l+1][l+1]-a1);
                      x[2]=u[2]=M[l+1][l]*M[l+2][l+1];
                      u[0]+=std::complex<double>( cos(arg(x[0])),sin(arg(x[0])) ) * sqrt(norm(x[0])+norm(x[1])+norm(x[2]));

                      P1=UnitMatrix<std::complex<double> >(3) - ( 1.+(Adjoint(x)*u)/(Adjoint(u)*x) ) * u * Adjoint(u) / (Adjoint(u)*u);

                      // multiply from the left by (P1,I)                      
                      for(j=l;j<=n-1;j++)
                              { m0j=P1[0][0]*M[l][j] + P1[0][1]*M[l+1][j] + P1[0][2]*M[l+2][j];
                                m1j=P1[1][0]*M[l][j] + P1[1][1]*M[l+1][j] + P1[1][2]*M[l+2][j];
                                m2j=P1[2][0]*M[l][j] + P1[2][1]*M[l+1][j] + P1[2][2]*M[l+2][j];
                                M[l][j]=m0j; 
                                M[l+1][j]=m1j; 
                                M[l+2][j]=m2j;
                               } 

                      // multiply from the right by (P1,I)                      
                      for(i=l;i<=n-1;i++)
                              { mi0=M[i][l]*P1[0][0] + M[i][l+1]*P1[1][0] + M[i][l+2]*P1[2][0];
                                mi1=M[i][l]*P1[0][1] + M[i][l+1]*P1[1][1] + M[i][l+2]*P1[2][1];
                                mi2=M[i][l]*P1[0][2] + M[i][l+1]*P1[1][2] + M[i][l+2]*P1[2][2];
                                M[i][l]=mi0; 
                                M[i][l+1]=mi1; 
                                M[i][l+2]=mi2;
                               } 

                      // chase the bulge - same as in UpperHessenberg function
                      if(l<=n-4)
                        { for(j=l;j<=n-4;j++)//sweep through columns up to four from last
                             { i=j+2;
                               if(norm(M[j+1][j])<norm(M[i][j])){ M.SwapRows(j+1,i); M.SwapColumns(j+1,i);}
                               if(norm(M[j+1][j])<norm(M[i+1][j])){ M.SwapRows(j+1,i+1); M.SwapColumns(j+1,i+1);}
                               if(M[j+1][j]!=0.)
                                 { C=M[i][j]/M[j+1][j];
                                   M[i][j]=0.;
                                   for(k=n-1;k>=j+1;k--){ M[i][k]-=M[j+1][k]*C;}  // subtract rows
                                   for(k=n-1;k>=0;k--){ M[k][j+1]+=M[k][i]*C;}    // add columns

                                   i=j+3;
                                   C=M[i][j]/M[j+1][j];
                                   M[i][j]=0.;
                                   for(k=n-1;k>=j+1;k--){ M[i][k]-=M[j+1][k]*C;}  // subtract rows
                                   for(k=n-1;k>=0;k--){ M[k][j+1]+=M[k][i]*C;}    // add columns
                                  }
                              }
                         }

                      // deal with the column three from the end
                      j=n-3;
                      i=j+2;
                      if(norm(M[j+1][j])<norm(M[i][j])){ M.SwapRows(j+1,i); M.SwapColumns(j+1,i);}
                      if(M[j+1][j]!=0.)
                        { C=M[i][j]/M[j+1][j];
                          M[i][j]=0.;
                          for(k=n-1;k>=j+1;k--){ M[i][k]-=M[j+1][k]*C;}  // subtract rows
                          for(k=n-1;k>=0;k--){ M[k][j+1]+=M[k][i]*C;}    // add columns
                         }

                      iterations++; 
                      //if(iterations==40){ throw NO_SOLUTION("QREigenValues");} 

                     }//end of if(rootfound==false) loop
                  else{ iterations=0;}
 
                 };//end of while(n-m>=3) loop

           if(n-m==2){ a1=M[m+1][m+1]+M[m][m];
                       a0=M[m+1][m+1]*M[m][m]-M[m+1][m]*M[m][m+1];
                       eigenvaluepair=QuadraticRoots(1.,-a1,a0);
                       eigenvalues[m]=eigenvaluepair[0];
                       eigenvalues[m+1]=eigenvaluepair[1];
                      }
           if(n-m==1){ eigenvalues[m]=M[m][m];}
           
           return eigenvalues;
          }

// *****************************************************************************************

// for real symmetric matrices
template <typename MType> 
std::vector<double> JacobiEigenValues(MATRIXEXPRESSION<double,MType> const &BB) 
         { 
           #ifdef _LAERRORS
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("JacobiEigenValues");}
           #endif

           MATRIX<double,0,0> B=BB;
           double Bil,Bjl, Bki,Bkj, Bii,Bij,Bji,Bjj;

           double gamma,r,r2,x;
           int N=static_cast<int>(BB.N2());
           bool finish;

           std::vector<double> W(N);

           do{ finish=true;
               for(int i=0;i<=N-2;i++)
                  { for(int j=i+1;j<=N-1;j++)
                       { x=4.*pow(B[i][j]/(B[i][i]-B[j][j]),2.);
                         if(x<1.){ gamma=0.5*(B[i][i]-B[j][j])*BinomialSeries(x,0.5);} else{ gamma=0.5*(B[i][i]-B[j][j])*( sqrt(1.+x)-1. );}

                         r2=1./( gamma*gamma+B[i][j]*B[i][j] );
                         if(r2!=0.){ r=sqrt(r2);} else{ r=0.;}

                         if(fabs(r*gamma) > std::numeric_limits<double>::epsilon() )
                           { finish=false;

                             Bii=B[i][i]; Bij=B[i][j]; Bji=B[j][i]; Bjj=B[j][j];
                             B[i][i]=r2*( Bii*Bij*Bji + gamma*gamma*Bjj + 2.*gamma*Bij*Bji );                           
                             B[j][j]=r2*( Bjj*Bij*Bji + gamma*gamma*Bii - 2.*gamma*Bij*Bji );
                             B[i][j]=B[j][i]=0.;

                             for(int k=0;k<=i-1;k++){ Bki=B[k][i]; Bkj=B[k][j]; B[k][i]=r*( Bki*Bij + gamma*Bkj); B[k][j]=r*( Bkj*Bji - gamma*Bki);}
                             for(int k=i+1;k<=j-1;k++){ Bki=B[k][i]; Bkj=B[k][j]; B[k][i]=r*( Bki*Bij + gamma*Bkj); B[k][j]=r*( Bkj*Bji - gamma*Bki);}
                             for(int k=j+1;k<=N-1;k++){ Bki=B[k][i]; Bkj=B[k][j]; B[k][i]=r*( Bki*Bij + gamma*Bkj); B[k][j]=r*( Bkj*Bji - gamma*Bki);}
 
                             for(int l=0;l<=i-1;l++){ Bil=B[i][l]; Bjl=B[j][l]; B[i][l]=r*( Bji*Bil + gamma*Bjl); B[j][l]=r*( Bij*Bjl - gamma*Bil);}
                             for(int l=i+1;l<=j-1;l++){ Bil=B[i][l]; Bjl=B[j][l]; B[i][l]=r*( Bji*Bil + gamma*Bjl); B[j][l]=r*( Bij*Bjl - gamma*Bil);}
                             for(int l=j+1;l<=N-1;l++){ Bil=B[i][l]; Bjl=B[j][l]; B[i][l]=r*( Bji*Bil + gamma*Bjl); B[j][l]=r*( Bij*Bjl - gamma*Bil);}
                            }
                        }
                   }
              }while(finish==false);

           for(int j=0;j<=N-1;j++){ W[j]=B[j][j];}

           return W;
          }

// for Hermitian matrices
template <typename MType> 
std::vector<std::complex<double> > JacobiEigenValues(MATRIXEXPRESSION<std::complex<double>,MType> const &BB) 

         { 
           #ifdef _LAERRORS
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("JacobiEigenValues");}
           #endif

           MATRIX<std::complex<double>,0,0> B=BB;
           std::complex<double> Bil,Bjl, Bki,Bkj, Bii,Bij,Bji,Bjj;

           double gamma,r,r2,x;
           int N=static_cast<int>(BB.N2());
           bool finish;

           std::vector<std::complex<double> > W(N);

           do{ finish=true;
               for(int i=0;i<=N-2;i++)
                  { for(int j=i+1;j<=N-1;j++)
                       { x=4.*norm(B[i][j])/norm(B[i][i]-B[j][j]);
                         if(x<1.){ gamma=0.5*real(B[i][i]-B[j][j])*BinomialSeries(x,0.5);} else{ gamma=0.5*real(B[i][i]-B[j][j])*( sqrt(1.+x)-1. );}

                         r2=1./( gamma*gamma+norm(B[i][j]) );
                         if(r2!=0.){ r=sqrt(r2);} else{ r=0.;}

                         if(fabs(r*gamma) > std::numeric_limits<double>::epsilon() )
                           { finish=false;

                             Bii=B[i][i]; Bij=B[i][j]; Bji=B[j][i]; Bjj=B[j][j];
                             B[i][i]=r2*( Bii*norm(Bij) + gamma*gamma*Bjj + 2.*gamma*norm(Bij) );                           
                             B[j][j]=r2*( Bjj*norm(Bij) + gamma*gamma*Bii - 2.*gamma*norm(Bij) );
                             B[i][j]=B[j][i]=0.;

                             for(int k=0;k<=i-1;k++){ Bki=B[k][i]; Bkj=B[k][j]; B[k][i]=r*( Bki*Bij + gamma*Bkj); B[k][j]=r*( Bkj*Bji - gamma*Bki);}
                             for(int k=i+1;k<=j-1;k++){ Bki=B[k][i]; Bkj=B[k][j]; B[k][i]=r*( Bki*Bij + gamma*Bkj); B[k][j]=r*( Bkj*Bji - gamma*Bki);}
                             for(int k=j+1;k<=N-1;k++){ Bki=B[k][i]; Bkj=B[k][j]; B[k][i]=r*( Bki*Bij + gamma*Bkj); B[k][j]=r*( Bkj*Bji - gamma*Bki);}
 
                             for(int l=0;l<=i-1;l++){ Bil=B[i][l]; Bjl=B[j][l]; B[i][l]=r*( Bji*Bil + gamma*Bjl); B[j][l]=r*( Bij*Bjl - gamma*Bil);}
                             for(int l=i+1;l<=j-1;l++){ Bil=B[i][l]; Bjl=B[j][l]; B[i][l]=r*( Bji*Bil + gamma*Bjl); B[j][l]=r*( Bij*Bjl - gamma*Bil);}
                             for(int l=j+1;l<=N-1;l++){ Bil=B[i][l]; Bjl=B[j][l]; B[i][l]=r*( Bji*Bil + gamma*Bjl); B[j][l]=r*( Bij*Bjl - gamma*Bil);}
                            }
                        }
                   }
              }while(finish==false);

           for(int j=0;j<=N-1;j++){ W[j]=B[j][j];}

           return W;
          }

// *************************************************************************

template <typename Type,typename MType> 
Type ConditionNumber(MATRIXEXPRESSION<Type,MType> const &M)
         { Type T=InverseConditionNumber(M);
           #ifdef _LAERRORS
           if(Equality(T,Zero<Type>())==true){ throw DIVISION_BY_ZERO("ConditionNumber","MATRIX");}
           #endif
           return One<Type>()/T;
          }

template <typename Type,typename MType> 
Type InverseConditionNumber(MATRIXEXPRESSION<Type,MType> const &M)
         { std::vector<Type> W(MEigenValues(M));
           return Smallest(W)/Largest(W);
          }

// *************************************************************************

template <typename MType> 
CVECTOR<double> EigenVector(MATRIXEXPRESSION<double,MType> const &M0,double lambda)
            { if(M0.N1()!=M0.N2()){ throw NOT_SQUARE("EigenVector");}

              MATRIX<double,0,0> M(M0-lambda);

              int n,i,imax=static_cast<int>(M.N1())-1,j,jmax=static_cast<int>(M.N2())-2,k;
              double C;

              for(j=0;j<=jmax;j++) // go through the columns
                 { C=fabs(M[j][j]); k=j;
                   for(i=j+1;i<=std::min(imax,j);i++){ if(fabs(M[i][j])>C){ k=i; C=fabs(M[i][j]);} }  // find largest element in column
                   if(k>j){ M.SwapRows(j,k);}                                                         // pivot the rows

                   for(i=j+1;i<=std::min(imax,j);i++)
                      { C=M[i][j]/M[j][j];
                        if(Equality(C,0.)==false){ for(k=std::min(jmax+1,j);k>=std::max(0,i);k--){ M[i][k]-=M[j][k]*C;} }// subtract rows 
                       }
                  }

              CVECTOR<double> X(static_cast<int>(M.N2()));
              double sum;

              n=static_cast<int>(X.N())-1;
              while( Equality(M[n][n],0.)==false && n>=1){ X[n]=0.; --n;};
              X[n]=1.;

              for(i=n-1;i>=0;i--)
                 { sum=0;
                   for(j=i+1;j<=n;j++){ sum+=M[i][j]*X[j];}
                   X[i]=-sum/M[i][i];
                  }

              return X/sqrt(Transpose(X)*X);
             }

template <typename MType> 
CVECTOR<std::complex<double> >EigenVector(MATRIXEXPRESSION<std::complex<double>,MType> const &M0,std::complex<double> lambda)
            { if(M0.N1()!=M0.N2()){ throw NOT_SQUARE("EigenVector");}

              MATRIX<std::complex<double>,0,0> M(M0-lambda);

              int n,i,imax=static_cast<int>(M.N1())-1,j,jmax=static_cast<int>(M.N2())-2,k;
              double B;
              std::complex<double> C,zero(0.,0.);

              for(j=0;j<=jmax;j++) // go through the columns
                 { B=abs(M[j][j]); k=j;
                   for(i=j+1;i<=std::min(imax,j);i++){ if(abs(M[i][j])>B){ k=i; B=abs(M[i][j]);} }  // find largest element in column
                   if(k>j){ M.SwapRows(j,k);}                                                       // pivot the rows

                   for(i=j+1;i<=std::min(imax,j);i++)
                      { C=M[i][j]/M[j][j];
                        if(Equality(C,zero)==false){ for(k=std::min(jmax+1,j);k>=std::max(0,i);k--){ M[i][k]-=M[j][k]*C;} }// subtract rows 
                       }
                  }

              CVECTOR<std::complex<double> > X(M.N2());
              std::complex<double> sum;

              n=static_cast<int>(X.N())-1;
              while( Equality(M[n][n],zero)==false && n>=1){ X[n]=0.; --n;};
              X[n]=1.;

              for(i=n-1;i>=0;i--)
                 { sum=zero;
                   for(j=i+1;j<=n;j++){ sum+=M[i][j]*X[j];}
                   X[i]=-sum/M[i][i];
                  }

              return X/sqrt(Adjoint(X)*X);
             }

// *************************************************************************
// ****************** OTHER MATRIX AND VECTOR OPERATIONS *******************
// *************************************************************************

template <typename Type,typename MType> 
MATRIX<Type,0,0> SwapRows(MATRIXEXPRESSION<Type,MType> const &M,int i,int j)
         { return MATRIX<Type,0,0>(M).SwapRows(i,j);}

template <typename Type,typename MType> 
MATRIX<Type,0,0> SwapColumns(MATRIXEXPRESSION<Type,MType> const &M,int i,int j)
         { return MATRIX<Type,0,0>(M).SwapColumns(i,j);}

template <typename Type,typename MType> 
MATRIX<Type,0,0> pow(MATRIXEXPRESSION<Type,MType> const &M,int N)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("pow(MATRIX<Type>,double)");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("pow(MATRIX<Type>,double)");}
           #endif
           if(N==0){ return UnitMatrix<Type>(M.N1());}
           if(N==1){ return M;}
           if(N==-1){ return Inverse(M);}

           if(N>1){ return M*pow(M,N-1);}
           if(N<-1){ return pow(M,N+1)/M;}
          }

// return the lower triangle part of the Cholesky decomposition for positive definite matrices
template <typename Type,typename MType> 
MATRIX<Type,0,0> sqrt(MATRIXEXPRESSION<Type,MType> const &M) 
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("sqrt");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("sqrt");}
           #endif
           MATRIX<Type,0,0> L(M.N1(),M.N2());
           Type Lii,Lji;
           int i,imax=static_cast<int>(M.N1())-1,j,k;

           for(i=0;i<=imax;i++)
              { Lii=M(i,i);
                for(k=0;k<=i-1;k++){ Lii-=pow(L(i,k),2.);}
                if(Lii<Zero<Type>()){ throw NEGATIVE_NUMBER("sqrt","MATRIX");}
                L(i,i)=sqrt(Lii);
                for(j=i+1;j<=imax;j++)
                   { if(Equality(L(i,i),Zero<Type>())==true){ throw DIVISION_BY_ZERO("sqrt","MATRIX");}
                     Lji=M(i,j);
                     for(k=0;k<=i-1;k++){ Lji-=L(i,k)*L(j,k);}
                     L(j,i)=Lji/L(i,i);
                    }
               }
           return L;
          }

// ***********************************

template <typename MType> 
MATRIX<double,0,0> Balance(MATRIXEXPRESSION<double,MType> const &A) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("Balance");}
           if(A.N1()!=A.N2()){ throw NOT_SQUARE("Balance");} 
           #endif

           MATRIX<double,0,0> B(A);
           double rnorm,cnorm,n,f;
           bool finish;

           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           // iterate through rows and columns multiplying/dividing by a factor f such that the norm of the ith row is close to the norm of the ith column
           do{ finish=true;
               for(i=0;i<=imax;i++)
                  { rnorm=cnorm=0;
                    for(j=0;j<=i-1;j++){ rnorm+=B[i][j]*B[i][j]; cnorm+=B[j][i]*B[j][i];} 
                    for(j=i+1;j<=jmax;j++){ rnorm+=B[i][j]*B[i][j]; cnorm+=B[j][i]*B[j][i];} 

                    if(rnorm!=0. && cnorm!=0.)
                      { n=round( log(cnorm/rnorm)/(4.*M_LN2) ); 
                        if(n!=0){ finish=false;
                                  f=pow(2.,n);                                  
                                  for(j=0;j<=i-1;j++){ B[i][j]*=f; B[j][i]/=f;} 
                                  for(j=i+1;j<=jmax;j++){ B[i][j]*=f; B[j][i]/=f;} 
                                 }
                       }
                   }
              }
           while(finish==false);

           return B;
          }

template <typename MType> 
MATRIX<std::complex<double>,0,0> Balance(MATRIXEXPRESSION<std::complex<double>,MType> const &A) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("Balance");}
           if(A.N1()!=A.N2()){ throw NOT_SQUARE("Balance");} 
           #endif

           MATRIX<std::complex<double>,0,0> B(A);
           double rnorm,cnorm,n,f;
           bool finish;

           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           // iterate through rows and columns multiplying/dividing by a factor f such that the norm of the ith row is close to the norm of the ith column
           do{ finish=true;
               for(i=0;i<=imax;i++)
                  { rnorm=cnorm=0;
                    for(j=0;j<=i-1;j++){ rnorm+=norm(B[i][j]); cnorm+=norm(B[j][i]);} 
                    for(j=i+1;j<=jmax;j++){ rnorm+=norm(B[i][j]); cnorm+=norm(B[j][i]);} 

                    if(rnorm!=0. && cnorm!=0.)
                      { n=round( log(cnorm/rnorm)/(4.*M_LN2) );
                        if(n!=0){ finish=false;
                                  f=pow(2.,n);                                  
                                  for(j=0;j<=i-1;j++){ B[i][j]*=f; B[j][i]/=f;} 
                                  for(j=i+1;j<=jmax;j++){ B[i][j]*=f; B[j][i]/=f;} 
                                 }
                       }
                   }
              }
           while(finish==false);

           return B;
          }

template <typename Type,typename MType> 
MATRIX<Type,0,0> Pivot(MATRIXEXPRESSION<Type,MType> const &A,int a,int b) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("Pivot");}
           #endif

           MATRIX<double,0,0> B(A);

           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           for(i=0;i<=imax;i++){ std::swap(B[i][a],B[i][b]);}
           for(j=0;j<=jmax;j++){ std::swap(B[a][j],B[b][j]);}

           return B;
          }

// *****************************************************************************************

template <typename MType> 
MATRIX<double,0,0> SortDiagonalAscending(MATRIXEXPRESSION<double,MType> const &A) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("SortDiagonalAscending");}
           if(A.N1()!=A.N2()){ throw NOT_SQUARE("SortDiagonalAscending");} 
           #endif

           MATRIX<double,0,0> B(A.N1(),A.N2());

           std::vector<int> ordering(B.N1());
           std::vector<double> diagonal(B.N1());
           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           for(i=0;i<=imax;i++){ ordering[i]=i; diagonal[i]=fabs(A[i][i]);}
           Sort(diagonal,ordering,ascending);

           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ B[i][j]=A[ordering[i]][ordering[j]];} }

           return B;
          }

template <typename MType> 
MATRIX<std::complex<double>,0,0> SortDiagonalAscending(MATRIXEXPRESSION<std::complex<double>,MType> const &A) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("SortDiagonalAscending");}
           if(A.N1()!=A.N2()){ throw NOT_SQUARE("SortDiagonalAscending");} 
           #endif

           MATRIX<std::complex<double>,0,0> B(A.N1(),A.N2());

           std::vector<int> ordering(B.N1());
           std::vector<double> diagonal(B.N1());
           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           for(i=0;i<=imax;i++){ ordering[i]=i; diagonal[i]=abs(A[i][i]);}
           Sort(diagonal,ordering,ascending);

           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ B[i][j]=A[ordering[i]][ordering[j]];} }

           return B;
          }

// *****************************************************************************************

template <typename MType> 
MATRIX<double,0,0> SortDiagonalDescending(MATRIXEXPRESSION<double,MType> const &A) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("SortDiagonalDescending");}
           if(A.N1()!=A.N2()){ throw NOT_SQUARE("SortDiagonalDescending");} 
           #endif

           MATRIX<double,0,0> B(A), C(A.N1(),A.N2());

           std::vector<int> ordering(B.N1());
           std::vector<double> diagonal(B.N1());
           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           for(i=0;i<=imax;i++){ ordering[i]=i; diagonal[i]=fabs(B[i][i]);}
           Sort(diagonal,ordering,descending);

           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ C[i][j]=B[ordering[i]][ordering[j]];} }

           return C;
          }

template <typename MType> 
MATRIX<std::complex<double>,0,0> SortDiagonalDescending(MATRIXEXPRESSION<std::complex<double>,MType> const &A) 
         { 
           #ifdef _LAERRORS
           if(A.N1()*A.N2()==0){ throw EMPTY("SortDiagonalDescending");}
           if(A.N1()!=A.N2()){ throw NOT_SQUARE("SortDiagonalDescending");} 
           #endif

           MATRIX<std::complex<double>,0,0> B(A), C(A.N1(),A.N2());

           std::vector<int> ordering(B.N1());
           std::vector<double> diagonal(B.N1());
           int i,imax=static_cast<int>(B.N1())-1, j,jmax=static_cast<int>(B.N2())-1;

           for(i=0;i<=imax;i++){ ordering[i]=i; diagonal[i]=abs(B[i][i]);}
           Sort(diagonal,ordering,descending);

           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ C[i][j]=B[ordering[i]][ordering[j]];} }

           return C;
          }

// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

template <typename Type,typename MType,typename CVType> 
CVECTOR<Type,0> GaussEliminate(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y)
            { try{ return BanddiagonalSolve(M,Y,M.N1()-1,M.N1()-1,false);}
              catch(EMPTY &E){ E.ChangeFunction("GaussEliminate"); throw E;}
              catch(NOT_SQUARE &NS){ NS.ChangeFunction("GaussEliminate"); throw NS;}
              catch(SINGULAR &S){ S.ChangeFunction("GaussEliminate"); throw S;}
              catch(DIFFERENT_SIZES &DS){ DS.ChangeFunction("GaussEliminate"); throw DS;}
              catch(INCORRECT_FORM &IF){ IF.ChangeFunction("GaussEliminate"); throw IF;}
             }

// ***********************************

template <typename Type,typename MType,typename CVType> 
CVECTOR<Type,0> DiagonalSolve(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,bool withtest) 
            { if(withtest==true){ try{ if(DiagonalTest(M)!=true){ throw INCORRECT_FORM("DiagonalSolve");}
                                       if(SpurTest(M)!=true){ throw SINGULAR("DiagonalSolve");}
                                       if(M.N1()!=M.N2()){ throw NOT_SQUARE("DiagonalSolve");}
                                      }
                                  catch(EMPTY &E){ E.ChangeFunction("DiagonalSolve"); throw E;}
                                 }

              try{ return BackSubstitution(M,Y,0,false);}
              catch(DIFFERENT_SIZES &DS){ DS.ChangeFunction("DiagonalSolve"); throw DS;}
             }

// ***********************************

template <typename Type,typename MType,typename CVType> 
CVECTOR<Type,0> TridiagonalSolve(MATRIXEXPRESSION<Type,MType> const &M,CVECTOREXPRESSION<Type,CVType> const &Y,bool withtest) 
            { if(withtest==true){ try{ if(TridiagonalTest(M)!=true){ throw INCORRECT_FORM("TridiagonalSolve");}
                                       if(SpurTest(M)!=true){ throw SINGULAR("TridiagonalSolve");}
                                       if(M.N1()!=M.N2()){ throw NOT_SQUARE("TridiagonalSolve");}
                                      }
                                  catch(EMPTY &E){ E.ChangeFunction("TridiagonalSolve"); throw E;}
                                 }
              try{ return BanddiagonalSolve(M,Y,1,1,false);}
              catch(SINGULAR &S){ S.ChangeFunction("TridiagonalSolve"); throw S;}
              catch(DIFFERENT_SIZES &DS){ DS.ChangeFunction("TridiagonalSolve"); throw DS;}
             }

// ***********************************

template <typename MType,typename CVType> 
CVECTOR<double> BanddiagonalSolve(MATRIXEXPRESSION<double,MType> const &M0,CVECTOREXPRESSION<double,CVType> const &Y0,int p,int q,bool withtest)
            { if(withtest==true){ if(BanddiagonalTest(M0,p,q)!=true){ throw INCORRECT_FORM("BanddiagonalSolve");} }
              if(M0.N2()!=Y0.N()){ throw DIFFERENT_SIZES("BanddiagonalSolve");}
              if(M0.N1()!=M0.N2()){ throw NOT_SQUARE("BanddiagonalSolve");}

              MATRIX<double,0,0> M(M0);
              CVECTOR<double> Y(Y0);

              int i,imax=static_cast<int>(M.N1())-1,j,jmax=static_cast<int>(M.N2())-2,k;
              double C,D;
              for(j=0;j<=jmax;j++) // go through the columns
                 { C=fabs(M[j][j]); k=j;
                   for(i=j+1;i<=std::min(imax,j+p);i++){ if(fabs(M[i][j])>C){ k=i; C=fabs(M[i][j]);} }  // find largest element in column
                   if(k>j){ M.SwapRows(j,k); std::swap(Y[j],Y[k]); q+=k;}                                        // pivot the rows

                   for(i=j+1;i<=std::min(imax,j+p);i++)
                      { D=M[i][j]/M[j][j];
                        if(Equality(D,0.)==false)
                          { for(k=std::min(jmax+1,j+p+q);k>=std::max(0,i-p);k--){ M[i][k]-=M[j][k]*D;}         // subtract rows
                            Y[i]-=D*Y[j];
                           }
                       }
                  }

              try{ return BackSubstitution(M,Y,p+q,false);}
              catch(EMPTY &E){ E.ChangeFunction("BanddiagonalSolve"); throw E;}
              catch(NOT_SQUARE &NS){ NS.ChangeFunction("BanddiagonalSolve"); throw NS;}
              catch(SINGULAR &S){ S.ChangeFunction("BanddiagonalSolve"); throw S;}
              catch(INCORRECT_FORM &IF){ IF.ChangeFunction("BanddiagonalSolve"); throw IF;}
             }

template <typename MType,typename CVType> 
CVECTOR<std::complex<double> > BanddiagonalSolve(MATRIXEXPRESSION<std::complex<double> ,MType> const &M0,CVECTOREXPRESSION<std::complex<double> ,CVType> const &Y0,int p,int q,bool withtest)
            { if(withtest==true){ if(BanddiagonalTest(M0,p,q)!=true){ throw INCORRECT_FORM("BanddiagonalSolve");} }
              if(M0.N2()!=Y0.N()){ throw DIFFERENT_SIZES("BanddiagonalSolve");}
              if(M0.N1()!=M0.N2()){ throw NOT_SQUARE("BanddiagonalSolve");}

              MATRIX<std::complex<double>,0,0> M(M0);
              CVECTOR<std::complex<double> > Y(Y0);

              int i,imax=static_cast<int>(M.N1())-1,j,jmax=static_cast<int>(M.N2())-2,k;
              double C;
              std::complex<double> D;
              for(j=0;j<=jmax;j++) // go through the columns
                 { C=norm(M[j][j]); k=j;
                   for(i=j+1;i<=std::min(imax,j+p);i++){ if(norm(M[i][j])>C){ k=i; C=norm(M[i][j]);} }  // find largest element in column
                   if(k>j){ M.SwapRows(j,k); std::swap(Y[j],Y[k]); q+=k;}                               // pivot the rows

                   for(i=j+1;i<=std::min(imax,j+p);i++)
                      { D=M[i][j]/M[j][j];
                        if(Equality(std::real(D),0.)==false && Equality(std::imag(D),0.)==false)
                          { for(k=std::min(jmax+1,j+p+q);k>=std::max(0,i-p);k--){ M[i][k]-=M[j][k]*D;}        // subtract rows
                            Y[i]-=D*Y[j];
                           }
                       }
                  }

              try{ return BackSubstitution(M,Y,p+q,false);}
              catch(EMPTY &E){ E.ChangeFunction("BanddiagonalSolve"); throw E;}
              catch(NOT_SQUARE &NS){ NS.ChangeFunction("BanddiagonalSolve"); throw NS;}
              catch(SINGULAR &S){ S.ChangeFunction("BanddiagonalSolve"); throw S;}
              catch(INCORRECT_FORM &IF){ IF.ChangeFunction("BanddiagonalSolve"); throw IF;}
             }

//*************************************************************************

template <std::size_t N> CVECTOR<double,N> TridiagonalSolve(MVECTOR<double,N> &L,MVECTOR<double,N> &D,MVECTOR<double,N> &U,CVECTOR<double,N> &Y)
     { CVECTOR<double,N> X, Y0(Y);

       try{ // effectively zero all the subdiagonal terms
            int i,imax=static_cast<int>(D.N())-2;
            for(int i=0;i<=imax;i++){ D[i+1]-=L[i]/D[i]*U[i]; Y[i+1]-=L[i]/D[i]*Y[i]; }

            // back substitution
            X[X.N()-1]=Y[Y.N()-1]/D[D.N()-1];
            for(int i=imax;i>=0;i--){ X[i]=( Y[i] - U[i]*X[i+1] )/D[i];}
            return X;
           }
       catch(...){ MATRIX<double,0,0> M(D.N(),D.N());
                   for(int i=0;i<=static_cast<int>(L.N())-2;i++){ M[i+1][i]=L[i];}
                   for(int i=0;i<=static_cast<int>(D.N())-1;i++){ M[i][i]=D[i];}
                   for(int i=0;i<=static_cast<int>(U.N())-2;i++){ M[i][i+1]=U[i];}
                   return TridiagonalSolve(M,Y0,false);
                  }
      }

//*************************************************************************

// Solve for M*X = Y for X given M is an upper triangular matrix
template <typename Type,typename MType,typename CVType> 
CVECTOR<Type> BackSubstitution(MATRIXEXPRESSION<Type,MType> const &M0,CVECTOREXPRESSION<Type,CVType> const &Y0,int width,bool withtest) 
            { if(withtest==true){ try{ if(UpperTriangleTest(M0)==false){ throw INCORRECT_FORM("BackSubstitution");}
                                       if(SpurTest(M0)==false){ throw SINGULAR("BackSubstitution");}
                                       if(M0.N1()!=M0.N2()){ throw NOT_SQUARE("BackSubstitution");}
                                      }
                                  catch(EMPTY &E){ E.ChangeFunction("BackSubstitution"); throw E;}
                                 }
              if(M0.N2()!=Y0.N()){ throw DIFFERENT_SIZES("BackSubstitution");}

              MATRIX<Type,0,0> M(M0);
              CVECTOR<Type> Y(Y0);
              CVECTOR<Type> X(Y.N());
              Type sum;

              int i,j;
              for(i=(int)X.N()-1;i>=(int)X.N()-width;i--)
                 { sum=Zero<Type>();
                   for(j=i+1;j<=(int)X.N()-1;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              for(i=(int)X.N()-width-1;i>=0;i--)
                 { sum=Zero<Type>();
                   for(j=i+1;j<=i+width;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              return X;
             }

// Solve for M*X = Y for X given M is an lower triangular matrix
template <typename Type,typename MType,typename CVType> 
CVECTOR<Type> ForwardSubstitution(MATRIXEXPRESSION<Type,MType> const &M0,CVECTOREXPRESSION<Type,CVType> const &Y0,int width,bool withtest) 
            { if(withtest==true){ try{ if(LowerTriangleTest(M0)==false){ throw INCORRECT_FORM("ForwardSubstitution");}
                                       if(SpurTest(M0)==false){ throw SINGULAR("ForwardSubstitution");}
                                       if(M0.N1()!=M0.N2()){ throw NOT_SQUARE("ForwardSubstitution");}
                                      }
                                  catch(EMPTY &E){ E.ChangeFunction("ForwardSubstitution"); throw E;}
                                 }
              if(M0.N2()!=Y0.N()){ throw DIFFERENT_SIZES("ForwardSubstitution");}

              MATRIX<Type,0,0> M(M0);
              CVECTOR<Type> Y(Y0);
              CVECTOR<Type> X(Y.N());

              Type sum;
              int i,j;
              for(i=0;i<=width;i++)
                 { sum=Zero<Type>();
                   for(j=0;j<=i-1;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              for(i=width+1;i<=(int)X.N()-1;i++)
                 { sum=Zero<Type>();
                   for(j=i-width;j<=i-1;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              return X;
             }

// *****************************************************************************************

template <typename MType> 
MATRIX<double,0,0> UpperHessenberg(MATRIXEXPRESSION<double,MType> const &M1) 
         { 
           #ifdef _LAERRORS
           if(M1.N1()*M1.N2()==0){ throw EMPTY("UpperHessenberg");}
           #endif

           MATRIX<double,0,0> M2(M1);

           int i,imax=static_cast<int>(M1.N1())-1,j,jmax=static_cast<int>(M1.N2())-3,k;
           double C;
           for(j=0;j<=jmax;j++) // go through the columns
              { C=fabs(M2[j+1][j]); k=j+1;
                for(i=j+2;i<=imax;i++){ if(fabs(M2[i][j])>C){ k=i; C=fabs(M2[i][j]);} }  // find largest element in subdiagonal
                if(k>j+1){ M2.SwapRows(j+1,k); M2.SwapColumns(j+1,k);}                             // swap the rows and columns

                if(fabs(M2[j+1][j])!=0.)
                  { for(i=j+2;i<=imax;i++) // sweep through the rows starting with row 2 below the diagonal
                       { C=-M2[i][j]/M2[j+1][j];
                         M2[i][j]=0.;
                         for(k=jmax+2;k>=j+1;k--){ M2[i][k]-=M2[j+1][k]*C;}      // subtract rows
                         for(k=imax;k>=0;k--){ M2[k][j+1]+=M2[k][i]*C;}        // add columns
                        }
                   }
               }

           return M2;
          }

template <typename MType> 
MATRIX<std::complex<double>,0,0> UpperHessenberg(MATRIXEXPRESSION<std::complex<double>,MType> const &M1) 
         { 
           #ifdef _LAERRORS
           if(M1.N1()*M1.N2()==0){ throw EMPTY("UpperHessenberg");}
           #endif

           MATRIX<std::complex<double>,0,0> M2(M1);

           int i,imax=static_cast<int>(M1.N1())-1,j,jmax=static_cast<int>(M1.N2())-3,k;
           double C;
           std::complex<double> Z; 
           for(j=0;j<=jmax;j++) // go through the columns
              { C=abs(M2[j+1][j]); k=j+1;
                for(i=j+2;i<=imax;i++){ if(abs(M2[i][j])>C){ k=i; C=abs(M2[i][j]);} }  // find largest element in column
                if(k>j+1){ M2.SwapRows(j+1,k); M2.SwapColumns(j+1,k);}                           // swap the rows and columns

                if(abs(M2[j+1][j])!=0.)
                  { for(i=j+2;i<=imax;i++) // sweep through the rows starting with row 2 below the diagonal
                       { Z=M2[i][j]/M2[j+1][j];
                         M2[i][j]=0.;
                         for(k=jmax+2;k>=j+1;k--){ M2[i][k]-=M2[j+1][k]*Z;}      // subtract rows
                         for(k=imax;k>=0;k--){ M2[k][j+1]+=M2[k][i]*Z;}        // add columns
                        }
                    }
               }

           return M2;
          }

// *****************************************************************************************

template <typename MType> 
MATRIX<double,0,0> LowerHessenberg(MATRIXEXPRESSION<double,MType> const &M1) 
         { 
           #ifdef _LAERRORS
           if(M1.N1()*M1.N2()==0){ throw EMPTY("LowerHessenberg");}
           #endif

           MATRIX<double,0,0> M2(M1);

           int i,imax=static_cast<int>(M1.N1())-3,j,jmax=static_cast<int>(M1.N2())-1,k;
           double C;
           for(i=0;i<=imax;i++) // go through the rows
              { C=fabs(M2[i][i+1]); k=i+1;
                for(j=i+2;j<=jmax;j++){ if(abs(M2[i][j])>C){ k=j; C=fabs(M2[i][j]);} }  // find largest element in column
                if(k>i+1){ M2.SwapColumns(i+1,k); M2.SwapRows(i+1,k);}                            // swap the columns and rows

                for(j=i+2;j<=jmax;j++) // sweep through the columns starting with column 2 to the right of the diagonal
                   { C=M2[i][j]/M2[i][i+1];
                     M2[i][j]=0.;                     
                     for(k=imax+2;k>=i+1;k--){ M2[k][j]-=M2[k][i+1]*C;}  // subtract columns
                     for(k=jmax;k>=0;k--){ M2[i+1][k]+=M2[j][k]*C;}    // add rows
                    }
               }

           return M2;
          }

template <typename MType> 
MATRIX<std::complex<double>,0,0> LowerHessenberg(MATRIXEXPRESSION<std::complex<double>,MType> const &M1) 
         { 
           #ifdef _LAERRORS
           if(M1.N1()*M1.N2()==0){ throw EMPTY("LowerHessenberg");}
           #endif

           MATRIX<std::complex<double>,0,0> M2(M1);

           int i,imax=static_cast<int>(M1.N1())-3,j,jmax=static_cast<int>(M1.N2())-1,k;
           double C;
           std::complex<double> Z; 
           for(i=0;i<=imax;i++) // go through the rows
              { C=abs(M2[i][i+1]); k=i+1;
                for(j=i+2;j<=jmax;j++){ if(abs(M2[i][j])>C){ k=j; C=abs(M2[i][j]);} }  // find largest element in column
                if(k>i+1){ M2.SwapColumns(i+1,k); M2.SwapRows(i+1,k);}                 // swap the columns and rows

                for(j=i+2;j<=jmax;j++) // sweep through the columns starting with column 2 to the right of the diagonal
                   { Z=M2[i][j]/M2[i][i+1];
                     M2[i][j]=0.;                     
                     for(k=imax+2;k>=i+1;k--){ M2[k][j]-=M2[k][i+1]*Z;} // subtract columns
                     for(k=jmax;k>=0;k--){ M2[i+1][k]+=M2[j][k]*Z;}     // add rows
                    }
               }

           return M2;
          }

// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

template <typename Type,typename MType> CVECTOR<Type,0> Vectorize(MATRIXEXPRESSION<double,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("Vectorize");}
           #endif

           CVECTOR<Type,0> V(M.N1()*M.N2());

           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1, offset=0;
           for(j=0;j<=jmax;j++){ for(i=0;i<=imax;i++){ V[offset+i]=M[i][j];} offset+=imax+1; }

           return V;
          }

template <typename Type,typename MType> CVECTOR<Type,0> HalfVectorize(MATRIXEXPRESSION<double,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(M.N1()*M.N2()==0){ throw EMPTY("HalfVectorize");}
           if(M.N1()!=M.N2()){ throw NOT_SQUARE("HalfVectorize");}
           #endif

           CVECTOR<Type,0> V(M.N1()*(M.N2()+1)/2);

           int i,imax=static_cast<int>(M.N1())-1, j,jmax=static_cast<int>(M.N2())-1, offset=0;
           for(j=0;j<=jmax;j++){ for(i=j;i<=imax;i++){ V[offset+i]=M[i][j];} offset+=imax-j;}

           return V;
          }

// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

template <std::size_t N> CVECTOR<double,N> LUSolve(MATRIX<double,N,N> &M,CVECTOR<double,N> &Y)
            { std::vector<MATRIX<double,N,N> > LU(LUDecomposition(M));
              return CVECTOR<double,N>(UTInverse(LU[1])*LTInverse(LU[0])*Y);
             }


#endif



