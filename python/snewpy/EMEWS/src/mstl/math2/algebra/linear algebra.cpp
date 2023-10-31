#include "linear algebra.h"

using std::vector;

using std::min;
using std::max;

using std::size_t;

using std::complex;

// **********************************************************************************************

 MATRIX<double,0,0> Givens(std::size_t N,int N1,int N2,double ANGLE)
         { if(static_cast<int>(N1)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,static_cast<int>(N)-1,"Givens(ANGLE)");}
           if(static_cast<int>(N2)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,static_cast<int>(N)-1,"Givens(ANGLE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"Givens(ANGLE)");}

           MATRIX<double,0,0> R(UnitMatrix<double>(N));
           R(N1,N2)=-R(N1,N1)*sin(ANGLE);   R(N1,N1)*=cos(ANGLE);
           R(N2,N1)=R(N2,N2)*sin(ANGLE);    R(N2,N2)*=cos(ANGLE);
           return R;
          }

// *******

MATRIX<double,0,0> Givens(std::size_t N,int N1,int N2,double COS,double SIN)
         { if(static_cast<int>(N1)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,static_cast<int>(N)-1,"Givens(COS,SIN)");}
           if(static_cast<int>(N2)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,static_cast<int>(N)-1,"Givens(COS,SIN)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"Givens(COS,SIN)");}

           MATRIX<double,0,0> R(UnitMatrix<double>(N));
           R(N1,N2)=-R(N1,N1)*SIN;   R(N1,N1)*=COS;
           R(N2,N1)=R(N2,N2)*SIN;    R(N2,N2)*=COS;
           return R;
          }

// *******

MATRIX<complex<double>,0,0> ComplexGivens(std::size_t N,int N1,int N2,double ANGLE,double PHASE)
         { if(static_cast<int>(N1)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,static_cast<int>(N)-1,"ComplexGivens(ANGLE,PHASE)");}
           if(static_cast<int>(N2)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,static_cast<int>(N)-1,"ComplexGivens(ANGLE,PHASE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"ComplexGivens(ANGLE,PHASE)");}

           MATRIX<complex<double>,0,0> R(UnitMatrix<complex<double> >(N));
           R(N1,N2)=-R(N1,N1)*sin(ANGLE)*complex<double>( cos(PHASE),sin(PHASE) );   R(N1,N1)*=cos(ANGLE);
           R(N2,N1)=R(N2,N2)*sin(ANGLE)*complex<double>( cos(PHASE),-sin(PHASE) );   R(N2,N2)*=cos(ANGLE);
           return R;
          }

// *******

MATRIX<complex<double>,0,0> ComplexGivens(std::size_t N,int N1,int N2,double COS,double SIN,double PHASE)
         { if(static_cast<int>(N1)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,static_cast<int>(N)-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(static_cast<int>(N2)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,static_cast<int>(N)-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"ComplexGivens(COS,SIN,PHASE)");}

           MATRIX<complex<double>,0,0> R(UnitMatrix<complex<double> >(N));
           R(N1,N2)=-R(N1,N1)*SIN*complex<double>( cos(PHASE),sin(PHASE) );   R(N1,N1)*=COS;
           R(N2,N1)=R(N2,N2)*SIN*complex<double>( cos(PHASE),-sin(PHASE) );   R(N2,N2)*=COS;
           return R;
          }

// *******

MATRIX<complex<double>,0,0> ComplexGivens(std::size_t N,int N1,int N2,complex<double> ALPHA,complex<double> BETA)
         { if(static_cast<int>(N1)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N1,0,static_cast<int>(N)-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(static_cast<int>(N2)>static_cast<int>(N)-1){ throw OUT_OF_RANGE<int>(N2,0,static_cast<int>(N)-1,"ComplexGivens(COS,SIN,PHASE)");}
           if(N1==N2){ throw EQUAL_VALUES<int>(N1,"ComplexGivens(COS,SIN,PHASE)");}

           MATRIX<complex<double>,0,0> R(UnitMatrix<complex<double> >(N));
           R(N1,N2)=-R(N1,N1)*BETA;             R(N1,N1)*=ALPHA;
           R(N2,N1)=R(N2,N2)*std::conj(BETA);   R(N2,N2)*=ALPHA;
           return R;
          }

// **********************************************************************************************

MATRIX<double,0,0> IntegerMatrix(size_t N) // a diagonal matrix with entries equal to the integers
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("IntegerMatrix");}
          #endif
          MATRIX<double,0,0> I(N,N);
          int i,imax=static_cast<int>(N)-1;
          for(i=1;i<=imax;i++){ I[i][i]=i*1.;}
          return I;
         }

MATRIX<double,0,0> FactorialMatrix(size_t N) // a diagonal matrix with entries equal to the factorial numbers
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("FactorialMatrix");}
          #endif
          MATRIX<double,0,0> M(UnitMatrix<double>(N));
          int i,imax=static_cast<int>(N)-1;
          for(i=2;i<=imax;i++){ M[i][i]*=Factorial(i*1.);}
          return M;
         }

//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************

CVECTOR<double,0> LUSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y)
            { if(M.N2()!=Y.N()){ throw DIFFERENT_SIZES("LUSolve");}
              vector<MATRIX<double,0,0> > LU(LUDecomposition(M));
              return CVECTOR<double,0>(UTInverse(LU[1])*LTInverse(LU[0])*Y);
             }

//*************************************************************************

CVECTOR<double,0> TridiagonalSolve(MVECTOR<double,0> &L,MVECTOR<double,0> &D,MVECTOR<double,0> &U,CVECTOR<double,0> &Y)
     { if(L.Empty()==true || D.Empty()==true || U.Empty()==true || Y.Empty()==true){ throw EMPTY("TridiagonalSolve");}
       if(D.N()!=Y.N()){ throw DIFFERENT_SIZES("TridiagonalSolve");}

       CVECTOR<double,0> X(Y.N()), Y0(Y);

       try{ // effectively zero all the subdiagonal terms
            int imax=static_cast<int>(D.N())-2;
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



