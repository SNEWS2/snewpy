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

MATRIX<double,0,0> PascalLMatrix(size_t N) // the Pascal L matrix
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("PascalLMatrix");}
          #endif
          MATRIX<double,0,0> L(UnitMatrix<double>(N));
          int i,imax=static_cast<int>(N)-1,j;
          for(i=1;i<=imax;i++){ for(j=0;j<=i-1;j++){ L[i][j]=BinomialCoefficient(i*1.,j*1.);} }
          return L;
         }

MATRIX<double,0,0> PascalUMatrix(size_t N) // the Pascal U matrix
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("PascalUMatrix");}
          #endif
          MATRIX<double,0,0> U(UnitMatrix<double>(N));
          int i,imax=static_cast<int>(N)-1,j,jmax=static_cast<int>(N)-1;
          for(i=1;i<=imax;i++){ for(j=i+1;j<=jmax;j++){ U[i][j]=BinomialCoefficient(j*1.,i*1.);} }
          return U;
         }

MATRIX<double,0,0> PascalSMatrix(size_t N)  // the Pascal S matrix
        { 
          #ifdef _LAERRORS
          if(N==0){ throw ZERO_NUMBER("PascalSMatrix");}
          #endif
          MATRIX<double,0,0> S(N,N);
          int i,imax=static_cast<int>(N)-1,j;
          for(i=0;i<=imax;i++){ for(j=0;j<=i;j++){ S[i][j]=S[j][i]=BinomialCoefficient((i+j)*1.,i*1.);} }
          return S;
         }

MATRIX<double,0,0> DiscreteChebyshevMatrix(size_t N)
         { 
           #ifdef _LAERRORS
           if(N==0){ throw ZERO_NUMBER("DiscreteChebyshevMatrix");}
           #endif
           MATRIX<double,0,0> C(N,N);
           int i,imax=static_cast<int>(N)-1,j,jmax=static_cast<int>(N)-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ C[i][j]=UnitNormalizedDiscreteChebyshev(N,j,i*1.);} }
           return C;
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

//CVECTOR<double,0> QRSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y)
//            { if(M.N2()!=Y.N()){ throw DIFFERENT_SIZES("QRSolve");}
//              vector<MATRIX<double,0,0> > QR(QRDecomposition(M));
//              return CVECTOR<double,0>(UTInverse(QR[1])*Transpose(QR[0])*Y);
//             }

//*************************************************************************

// Solve the equation M * X = Y for X
/*CVECTOR<double,0> GaussElimination(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y)
            { try{ return BanddiagonalSolve(M,Y,M.N1()-1,M.N1()-1,false);}
              catch(EMPTY E){ E.ChangeFunction("GaussElimination"); throw E;}
              catch(NOT_SQUARE NS){ NS.ChangeFunction("GaussElimination"); throw NS;}
              catch(SINGULAR S){ S.ChangeFunction("GaussElimination"); throw S;}
              catch(DIFFERENT_SIZES DS){ DS.ChangeFunction("GaussElimination"); throw DS;}
              catch(INCORRECT_FORM IF){ IF.ChangeFunction("GaussElimination"); throw IF;}
             }*/

/*CVECTOR<double,0> DiagonalSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y,bool withtest)
            { if(withtest==true){ try{ if(DiagonalTest(M)!=true){ throw INCORRECT_FORM("DiagonalSolve");}
                                       if(SpurTest(M)!=true){ throw SINGULAR("DiagonalSolve");}
                                       if(M.N1()!=M.N2()){ throw NOT_SQUARE("DiagonalSolve");}
                                      }
                                  catch(EMPTY E){ E.ChangeFunction("DiagonalSolve"); throw E;}
                                 }

              try{ return BackSubstitution(M,Y,0,false);}
              catch(DIFFERENT_SIZES DS){ DS.ChangeFunction("DiagonalSolve"); throw DS;}
             }*/

/*CVECTOR<double,0> TridiagonalSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y,bool withtest)
            { if(withtest==true){ try{ if(TridiagonalTest(M)!=true){ throw INCORRECT_FORM("TridiagonalSolve");}
                                       if(SpurTest(M)!=true){ throw SINGULAR("TridiagonalSolve");}
                                       if(M.N1()!=M.N2()){ throw NOT_SQUARE("TridiagonalSolve");}
                                      }
                                  catch(EMPTY E){ E.ChangeFunction("TridiagonalSolve"); throw E;}
                                 }
              try{ return BanddiagonalSolve(M,Y,1,1,false);}
              catch(SINGULAR S){ S.ChangeFunction("TridiagonalSolve"); throw S;}
              catch(DIFFERENT_SIZES DS){ DS.ChangeFunction("TridiagonalSolve"); throw DS;}
             }*/

/*CVECTOR<double,0> BanddiagonalSolve(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y,int p,int q,bool withtest)
            { if(withtest==true){ if(BanddiagonalTest(M,p,q)!=true){ throw INCORRECT_FORM("BanddiagonalSolve");} }
              if(M.N2()!=Y.N()){ throw DIFFERENT_SIZES("BanddiagonalSolve");}
              if(M.N1()!=M.N2()){ throw NOT_SQUARE("BanddiagonalSolve");}

              int i,j,k;
              bool pivoted=false;
              double C;
              for(j=0;j<=static_cast<int>(M.N2())-2;j++) // go through the columns
                 { C=fabs(M[j][j]); k=j;
                   for(i=j+1;i<=min(static_cast<int>(M.N1())-1,j+p);i++){ if(fabs(M[i][j])>C){ k=i; C=fabs(M[i][j]);} }  // find largest element in column
                   if(k>j){ M.SwapRows(j,k); std::swap(Y[j],Y[k]); pivoted=true; q+=k;}                     // pivot the rows

                   for(i=j+1;i<=min(static_cast<int>(M.N1())-1,j+p);i++)
                      { C=M[i][j]/M[j][j];
                        if(Equality(C,0.)==false)
                          { if(pivoted==false)
                               { for(k=min(static_cast<int>(M.N2())-1,j+q);k>=max(0,i-p);k--)
                                            { M[i][k]-=M[j][k]*C;}            // subtract rows
                                }
                            else{ for(k=min(static_cast<int>(M.N2())-1,j+p+q);k>=max(0,i-p);k--)
                                            { M[i][k]-=M[j][k]*C;}                          // subtract rows
                                 }
                            Y[i]-=C*Y[j];
                           }
                       }
                  }

              try{ if(pivoted==false){ return BackSubstitution(M,Y,q,false);}
                       else{ return BackSubstitution(M,Y,p+q,false);}
                      }
              catch(EMPTY E){ E.ChangeFunction("BanddiagonalSolve"); throw E;}
              catch(NOT_SQUARE NS){ NS.ChangeFunction("BanddiagonalSolve"); throw NS;}
              catch(SINGULAR S){ S.ChangeFunction("BanddiagonalSolve"); throw S;}
              catch(INCORRECT_FORM IF){ IF.ChangeFunction("BanddiagonalSolve"); throw IF;}
             }*/

//*************************************************************************

// Solve for M*X = Y for X given M is an upper triangular matrix
/*CVECTOR<double,0> BackSubstitution(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y,int width,bool withtest)
            { if(withtest==true){ try{ if(UpperTriangleTest(M)==false){ throw INCORRECT_FORM("BackSubstitution");}
                                       if(SpurTest(M)==false){ throw SINGULAR("BackSubstitution");}
                                       if(M.N1()!=M.N2()){ throw NOT_SQUARE("BackSubstitution");}
                                      }
                                  catch(EMPTY E){ E.ChangeFunction("BackSubstitution"); throw E;}
                                 }
              if(M.N2()!=Y.N()){ throw DIFFERENT_SIZES("BackSubstitution");}

              CVECTOR<double,0> X(Y.N());
              double sum;

              int i,j;
              for(i=static_cast<int>(X.N())-1;i>=static_cast<int>(X.N())-width;i--)
                 { sum=Zero<double>();
                   for(j=i+1;j<=static_cast<int>(X.N())-1;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              for(i=static_cast<int>(X.N())-width-1;i>=0;i--)
                 { sum=Zero<double>();
                   for(j=i+1;j<=i+width;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              return X;
             }*/

// Solve for M*X = Y for X given M is an lower triangular matrix
/*CVECTOR<double,0> ForwardSubstitution(MATRIX<double,0,0> &M,CVECTOR<double,0> &Y,int width,bool withtest)
            { if(withtest==true){ try{ if(LowerTriangleTest(M)==false){ throw INCORRECT_FORM("ForwardSubstitution");}
                                       if(SpurTest(M)==false){ throw SINGULAR("ForwardSubstitution");}
                                       if(M.N1()!=M.N2()){ throw NOT_SQUARE("ForwardSubstitution");}
                                      }
                                  catch(EMPTY E){ E.ChangeFunction("ForwardSubstitution"); throw E;}
                                 }
              if(M.N2()!=Y.N()){ throw DIFFERENT_SIZES("ForwardSubstitution");}

              CVECTOR<double,0> X(Y.N());

              double sum;
              int i,j;
              for(i=0;i<=width;i++)
                 { sum=Zero<double>();
                   for(j=0;j<=i-1;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              for(i=width+1;i<=static_cast<int>(X.N())-1;i++)
                 { sum=Zero<double>();
                   for(j=i-width;j<=i-1;j++){ sum+=M[i][j]*X[j];}
                   X[i]=(Y[i]-sum)/M[i][i];
                  }

              return X;
             }*/

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



