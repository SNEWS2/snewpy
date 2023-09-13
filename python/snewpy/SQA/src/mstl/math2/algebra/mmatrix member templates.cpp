#if !defined(_MMATRIX_MEMBERS)
#define _MMATRIX_MEMBERS

#include "mstl.h"

template <typename Type,std::size_t n1,std::size_t n2> std::string MATRIX<Type,n1,n2>::CLASS_NAME("MATRIX - defined size");
template <typename Type,std::size_t n> std::string MATRIX<Type,n,n>::CLASS_NAME("MATRIX - defined square");
template <typename Type> std::string MATRIX<Type,0,0>::CLASS_NAME("MATRIX - undefined size");

//*************************************************************************
//******************** PROTECTED MATRIX MEMBERS ***************************
//*************************************************************************

template <typename Type,std::size_t n1,std::size_t n2> 
void MATRIX<Type,n1,n2>::Copy(std::vector<std::vector<Type> > const &V)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=V[a][b];} } 
          }

template <typename Type,std::size_t n> 
void MATRIX<Type,n,n>::Copy(std::vector<std::vector<Type> > const &V)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=V[a][b];} } 
          }

template <typename Type> 
void MATRIX<Type,0,0>::Copy(std::vector<std::vector<Type> > const &V)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;           
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=V[a][b];} } 
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType>
inline void MATRIX<Type,n1,n2>::Copy(MATRIXEXPRESSION<Type,MType> const &M)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=M(a,b);} } 
          }

template <typename Type,std::size_t n> template <typename MType>
inline void MATRIX<Type,n,n>::Copy(MATRIXEXPRESSION<Type,MType> const &M)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=M(a,b);} } 
          }

template <typename Type> template <typename MType>
inline void MATRIX<Type,0,0>::Copy(MATRIXEXPRESSION<Type,MType> const &M)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=M(a,b);} } 
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> 
void MATRIX<Type,n1,n2>::Initialize(void)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=Zero<Type>();} } 
          }

template <typename Type,std::size_t n> 
void MATRIX<Type,n,n>::Initialize(void)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=Zero<Type>();} } 
          }

template <typename Type> 
void MATRIX<Type,0,0>::Initialize(void)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=Zero<Type>();} } 
          }

//*********

/*template <typename Type> 
void MATRIX<Type,0,0>::Create(const std::size_t N1,const std::size_t N2)
         { n1=N1; n2=N2; 
           if(N1==0 || N2==0){ empty=true; m=NULL;}
           else{ empty=false; 
                 m=new Type*[N1]; 
                 int a,amax=static_cast<int>(N1)-1;
                 for(a=0;a<=amax;a++){ m[a]=new Type[N2];}
                }
          }*/

template <typename Type> 
void MATRIX<Type,0,0>::Create(const std::size_t N1,const std::size_t N2)
         { n1=N1; n2=N2; 
           if(n1==0 || n2==0){ empty=true;} else{ empty=false;}
           m=std::vector<std::vector<Type> >(N1,std::vector<Type>(N2));                  
          }

/*template <typename Type> 
void MATRIX<Type,0,0>::Destroy(void)
         { if(Empty()==false)
               { int a,amax=static_cast<int>(N1())-1;
                 for(a=0;a<=amax;a++){ delete [](m[a]);}
                 delete []m; 
                 m=NULL;
                 empty=true; 
                }
          }*/

template <typename Type> 
void MATRIX<Type,0,0>::Destroy(void)
         { Create(0,0);}

//**********************************************************************
//**************************** PUBLIC MATRIX MEMBERS ******************
//**********************************************************************

template <typename Type,std::size_t n1,std::size_t n2> MATRIX<Type,n1,n2>::MATRIX(Type MM00,...)
      	{ va_list val; va_start(val,MM00);
 	  (*this)(0,0)=MM00; 
          int a,amax=static_cast<int>(n1)-1,b,bmax=static_cast<int>(n2)-1;
          for(b=1;b<=bmax;b++){ (*this)(0,b)=va_arg(val,Type);}
      	  for(a=1;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=va_arg(val,Type);} }
          va_end(val);
	 }

template <typename Type,std::size_t n> MATRIX<Type,n,n>::MATRIX(Type MM00,...)
      	{ va_list val; va_start(val,MM00);

 	  (*this)(0,0)=MM00; 
          int a,amax=static_cast<int>(n)-1,b,bmax=static_cast<int>(n)-1;
          for(b=1;b<=bmax;b++){ (*this)(0,b)=va_arg(val,Type);}
      	  for(a=1;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=va_arg(val,Type);} }
          va_end(val);
	 }

template <typename Type> MATRIX<Type,0,0>::MATRIX(const std::size_t N1,const std::size_t N2,Type MM00,...)
      	{ va_list val; va_start(val,MM00);

          Create(N1,N2);
 	  if(Empty()==false)
            { (*this)(0,0)=MM00; 
              int a,amax=static_cast<int>(N1)-1,b,bmax=static_cast<int>(N2)-1;
              for(b=1;b<=bmax-1;b++){ (*this)(0,b)=va_arg(val,Type);}
      	      for(a=1;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=va_arg(val,Type);} }
             }
          va_end(val);
	 }

//********

template <typename Type,std::size_t n1,std::size_t n2> MATRIX<Type,n1,n2>::MATRIX(std::vector<std::vector<Type> > const &V)
         { Copy(V);}

template <typename Type,std::size_t n> MATRIX<Type,n,n>::MATRIX(std::vector<std::vector<Type> > const &V)
         { Copy(V);}

template <typename Type> MATRIX<Type,0,0>::MATRIX(std::vector<std::vector<Type> > const &V)
         { Create(V.size(),V[0].size()); Copy(V);}

//********

template <typename Type,std::size_t n1,std::size_t n2> MATRIX<Type,n1,n2>::MATRIX(MATRIX<Type,n1,n2> const &M)
         { Copy(M);}

template <typename Type,std::size_t n> MATRIX<Type,n,n>::MATRIX(MATRIX<Type,n,n> const &M)
         { Copy(M);}

template <typename Type> template<std::size_t m1,std::size_t m2> MATRIX<Type,0,0>::MATRIX(MATRIX<Type,m1,m2> const &M)
         { Create(M.N1(),M.N2()); Copy(M);}

//********

template <typename Type,std::size_t n1,std::size_t n2>
MATRIX<Type,n1,n2>::MATRIX(MATRIX<Type,0,0> const &M)
         { if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("MATRIX<Type,0,0> const&");}
           Copy(M);
          }

template <typename Type,std::size_t n>
MATRIX<Type,n,n>::MATRIX(MATRIX<Type,0,0> const &M)
         { if(N()!=M.N1() || N()!=M.N2()){ throw DIFFERENT_SIZES("MATRIX<Type,0,0> const&");}
           Copy(M);
          }

template <typename Type>
MATRIX<Type,0,0>::MATRIX(MATRIX<Type,0,0> const &M)
         { Create(M.N1(),M.N2()); Copy(M);}

//********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
MATRIX<Type,n1,n2>::MATRIX(MATRIXEXPRESSION<Type,MType> const &M)
         { if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("MATRIX(MATRIXEXPRESSION<Type,MType> const &");}
           Copy(M);
          }

template <typename Type,std::size_t n> template <typename MType> 
MATRIX<Type,n,n>::MATRIX(MATRIXEXPRESSION<Type,MType> const &M)
         { if(N()!=M.N1() || N()!=M.N2()){ throw DIFFERENT_SIZES("MATRIX(MATRIXEXPRESSION<Type,MType> const &");}
           Copy(M);
          }

template <typename Type> template <typename MType> 
MATRIX<Type,0,0>::MATRIX(MATRIXEXPRESSION<Type,MType> const &M)
         { Create(M.N1(),M.N2()); Copy(M);}

//**********************************************************************

template <typename Type,std::size_t n1,std::size_t n2> MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::Conjugate(void)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=conj((*this)(a,b));} }
           return (*this);
          }

template <typename Type,std::size_t n> MATRIX<Type,n,n>& MATRIX<Type,n,n>::Conjugate(void)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=conj((*this)(a,b));} }
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::Conjugate(void)
         { if(Empty()==true){ throw EMPTY("Conjugate",CLASS_NAME);}
           int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)=conj((*this)(a,b));} }
           return (*this);
          }

//********

template <typename Type,std::size_t n> MATRIX<Type,n,n>& MATRIX<Type,n,n>::Transpose(void)
         { int a,b, amax=static_cast<int>(N1())-1; 
           for(a=1;a<=amax;a++){ for(b=0;b<=a-1;b++){ std::swap((*this)(a,b),(*this)(b,a));} }
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::Transpose(void)
         { if(Empty()==true){ throw EMPTY("Transpose",CLASS_NAME);}
           if(N1()!=N2())
             { MATRIX<Type,0,0> MT(N2(),N1()); 
               int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
               for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ MT(b,a)=(*this)(a,b);} }
               return (*this)=MT;
              }
           else{ int a,b, amax=static_cast<int>(N1())-1; 
                 for(a=1;a<=amax;a++){ for(b=0;b<=a-1;b++){ std::swap((*this)(a,b),(*this)(b,a));} }
                 return (*this);
                }
          }

//********

template <typename Type,std::size_t n> MATRIX<Type,n,n>& MATRIX<Type,n,n>::AntiTranspose(void)
         { int a,b, amax=static_cast<int>(N1())-2,bmax=static_cast<int>(N2())-2;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax-a;b++){ std::swap((*this)(a,b),(*this)(N1()-1-b,N2()-1-a));} }
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::AntiTranspose(void)
         { if(Empty()==true){ throw EMPTY("AntiTranspose",CLASS_NAME);}
           MATRIX MAT(N2(),N1());
           int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ MAT(b,a)=(*this)(N1()-1-a,N2()-1-b);} }
           return (*this)=MAT;
          }

//********

template <typename Type,std::size_t n> MATRIX<Type,n,n>& MATRIX<Type,n,n>::Adjoint(void)
         { Transpose();
           Conjugate();
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::Adjoint(void)
         { if(Empty()==true){ throw EMPTY("Adjoint",CLASS_NAME);}
           Transpose();
           Conjugate();
           return (*this);
          }

//*************************************************************************

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::Invert(void)
         { if(N1()<=2){ try{ return LInvert();} catch(...){ throw SINGULAR("Invert");} }
           if(SpurTest(*this)==true){ if(DiagonalTest(*this)==true){ try{ return DInvert();} catch(...){ throw SINGULAR("Invert");} }
                                      if(LowerTriangleTest(*this)==true){ try{ return LTInvert();} catch(...){ throw SINGULAR("Invert");} }
                                      if(UpperTriangleTest(*this)==true){ try{ return UTInvert();} catch(...){ throw SINGULAR("Invert");} }
                                      try{ return LUInvert();} catch(...){ throw SINGULAR("Invert");} 
                                     }
           //try{ return QRInvert();} catch(...){;}
           try{ return GJInvert();} catch(...){ throw SINGULAR("Invert");}
           try{ return LInvert();}  catch(...){ throw SINGULAR("Invert");}
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::Invert(void)
         { if(N1()!=N2()){ try{ return MPInvert();} catch(...){ throw SINGULAR("Invert");} }
           else{ if(N1()<=2){ try{ return LInvert();} catch(...){ throw SINGULAR("Invert");} }
                 if(SpurTest(*this)==true){ if(DiagonalTest(*this)==true){ try{ return DInvert();} catch(...){ throw SINGULAR("Invert");} }
                                            if(LowerTriangleTest(*this)==true){ try{ return LTInvert();} catch(...){ throw SINGULAR("Invert");} }
                                            if(UpperTriangleTest(*this)==true){ try{ return UTInvert();} catch(...){ throw SINGULAR("Invert");} }
                                            try{ return LUInvert();} catch(...){ throw SINGULAR("Invert");} 
                                           }
                }
           //try{ return QRInvert();} catch(...){;}
           try{ return GJInvert();} catch(...){ throw SINGULAR("Invert");}
           try{ return LInvert();}  catch(...){ throw SINGULAR("Invert");}
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::DInvert(void) // Inversion of diagonal matrices
         { 
           #ifdef _LAERRORS
           if(DiagonalTest(*this)==false){ throw INCORRECT_FORM("DInvert");}
           if(SpurTest(*this)==false){ throw SINGULAR("DInvert");}
           #endif
           int a,amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ (*this)(a,a)=One<Type>()/(*this)(a,a);}

           return (*this);
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::DInvert(void) // Inversion of diagonal matrices
         { 
           #ifdef _LAERRORS
           if(DiagonalTest(*this)==false){ throw INCORRECT_FORM("DInvert");}
           if(SpurTest(*this)==false){ throw SINGULAR("DInvert");}
           #endif
           int a,amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ (*this)(a,a)=One<Type>()/(*this)(a,a);}

           return (*this);
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::GJInvert(void)
         { int a,amax=static_cast<int>(N1())-1,b,bmax=static_cast<int>(N2())-1,c,cmax=amax;
           double scale;

           MATRIX<Type,n,n> GJ(*this), IGJ(N1(),N2()); for(a=0;a<=amax;a++){ IGJ(a,a)=One<Type>();}

           for(c=0;c<=cmax;c++)
              { if(Equality(GJ(c,c),Zero<Type>())==true){ return LInvert();}
                scale=GJ(c,c); for(b=0;b<=bmax;b++){ GJ(c,b)/=scale; IGJ(c,b)/=scale;}

                for(a=0;a<=amax;a++)                  // go through the rows excepth the c'th
	                { if(a!=c && Equality(GJ(a,c),Zero<Type>())==false)
                            { scale=GJ(a,c);
                              for(b=0;b<=bmax;b++)      // go through the columns in that row
 	                         { GJ(a,b)-=scale*GJ(c,b);    // subtract from element IDM[a][b] so that same manipulation
	                           IGJ(a,b)-=scale*IGJ(c,b);  // on DM would produces a zero
	                          }
                             }
	                 }
               }
   	   return (*this)=IGJ;
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::GJInvert(void)
         { int a,amax=static_cast<int>(N1())-1,b,bmax=static_cast<int>(N2())-1,c,cmax=amax;
           double scale;

           MATRIX<Type,0,0> GJ(*this), IGJ(N1(),N2()); for(a=0;a<=amax;a++){ IGJ(a,a)=One<Type>();}

           for(c=0;c<=cmax;c++)
              { if(Equality(GJ(c,c),Zero<Type>())==true){ return LInvert();}
                scale=GJ(c,c); for(b=0;b<=bmax;b++){ GJ(c,b)/=scale; IGJ(c,b)/=scale;}

                for(a=0;a<=amax;a++)                  // go through the rows excepth the c'th
	                { if(a!=c && Equality(GJ(a,c),Zero<Type>())==false)
                            { scale=GJ(a,c);
                              for(b=0;b<=bmax;b++)      // go through the columns in that row
 	                         { GJ(a,b)-=scale*GJ(c,b);    // subtract from element IDM[a][b] so that same manipulation
	                           IGJ(a,b)-=scale*IGJ(c,b);  // on DM would produces a zero
	                          }
                             }
	                 }
               }
   	   return (*this)=IGJ;
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::LInvert(bool message) // Laplace inversion by Cofactors/Determinant
         { if(N1()>=3){ BASIC_MESSAGE BM("Beginning Matrix Inversion by Laplace Minors","LInvert",CLASS_NAME,message);}

           Type D=Determinant(*this);
           if(Equality(D,Zero<Type>())==true){ throw SINGULAR("LInvert");}
           MATRIX<Type,n,n> IL(N1(),N2());

           if(N1()==1){ if(Equality((*this)(0,0),Zero<Type>())==true){ throw SINGULAR("LInvert");}
                        IL(0,0)=1./(*this)(0,0);
                       }
           else{ int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ IL(b,a)=Cofactor(*this,a,b)/D;} }
                }
           return (*this)=IL;
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::LInvert(bool message) // Laplace inversion by Cofactors/Determinant
         { if(N1()>=3){ BASIC_MESSAGE BM("Beginning Matrix Inversion by Laplace Minors","LInvert",CLASS_NAME,message);}

           Type D=Determinant(*this);
           if(Equality(D,Zero<Type>())==true){ throw SINGULAR("LInvert");}
           MATRIX<Type,0,0> IL(N1(),N2());

           if(N1()==1){ if(Equality((*this)(0,0),Zero<Type>())==true){ throw SINGULAR("LInvert");}
                        IL(0,0)=1./(*this)(0,0);
                       }
           else{ int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ IL(b,a)=Cofactor(*this,a,b)/D;} }
                }
           return (*this)=IL;
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::LTInvert(void) // Inversion of a Lower Triangular matrix
         { 
           #ifdef _LAERRORS
           if(LowerTriangleTest(*this)==false){ throw INCORRECT_FORM("LTInvert");}
           if(SpurTest(*this)==false){ throw SINGULAR("LTInvert");}
           #endif
           MATRIX<Type,n,n> ILT((Diagonal(*this)).DInvert()); 
           int a,b,c, amax=static_cast<int>(N1())-1;
           for(a=1;a<=amax;a++)
              { for(b=a-1;b>=0;b--)
                   { Type C=Zero<Type>(); 
                     for(c=0;c<=a-b-1;c++){ C+=(*this)(a-c,b)*ILT(a,a-c);}
                     ILT(a,b)=-C/(*this)(b,b);
                    }
               }
           return (*this)=ILT;
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::LTInvert(void) // Inversion of a Lower Triangular matrix
         { 
           #ifdef _LAERRORS
           if(LowerTriangleTest(*this)==false){ throw INCORRECT_FORM("LTInvert");}
           if(SpurTest(*this)==false){ throw SINGULAR("LTInvert");}
           #endif
           MATRIX<Type,0,0> ILT((Diagonal(*this)).DInvert()); 
           int a,b,c, amax=static_cast<int>(N1())-1;
           for(a=1;a<=amax;a++)
              { for(b=a-1;b>=0;b--)
                   { Type C=Zero<Type>(); 
                     for(c=0;c<=a-b-1;c++){ C+=(*this)(a-c,b)*ILT(a,a-c);}
                     ILT(a,b)=-C/(*this)(b,b);
                    }
               }
           return (*this)=ILT;
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::UTInvert(void) // Inversion of an Upper Triangular matrix
         { try{ Transpose(); LTInvert(); Transpose(); return (*this);}
           catch(SINGULAR &S){ S.ChangeFunction("UTInvert"); throw S;}
           catch(INCORRECT_FORM &IF){ IF.ChangeFunction("UTInvert"); throw IF;}
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::UTInvert(void) // Inversion of an Upper Triangular matrix
         { try{ Transpose(); LTInvert(); Transpose(); return (*this);}
           catch(SINGULAR &S){ S.ChangeFunction("UTInvert"); throw S;}
           catch(INCORRECT_FORM &IF){ IF.ChangeFunction("UTInvert"); throw IF;}
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::LUInvert(void)
         { std::vector<MATRIX<Type,0,0> > LU(LUDecomposition(*this)); //L is the first, U the second
           return (*this)=LU[1].UTInvert()*LU[0].LTInvert();
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::LUInvert(void)
         { std::vector<MATRIX<Type,0,0> > LU(LUDecomposition(*this)); //L is the first, U the second
           return (*this)=LU[1].UTInvert()*LU[0].LTInvert();
          }

//***********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::QRInvert(void)
         { std::vector<MATRIX<Type,0,0> > QR=QRDecomposition(*this); //Q is the first, R the second
           return (*this)=QR[1].UTInvert()*QR[0].Transpose();
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::QRInvert(void)
         { std::vector<MATRIX<Type,0,0> > QR=QRDecomposition(*this); //Q is the first, R the second
           return (*this)=QR[1].UTInvert()*QR[0].Transpose();
          }

//***********

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::MPInvert(void)
         { MATRIX<Type,0,0> TMM(MATRIX<Type,0,0>(*this).Transpose());
           return (*this)=( MATRIX<Type,0,0>(TMM*(*this)) ).Invert()*TMM;
          }

//*************************************************************************

template <typename Type,std::size_t n1,std::size_t n2> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::SwapRows(int i,int j)
         { int kmax=static_cast<int>(N2())-1;
           for(int k=0;k<=kmax;k++){ std::swap((*this)(i,k),(*this)(j,k));}
           return (*this);
          }

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::SwapRows(int i,int j)
         { int kmax=static_cast<int>(N2())-1;
           for(int k=0;k<=kmax;k++){ std::swap((*this)(i,k),(*this)(j,k));}
           return (*this);
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::SwapRows(int i,int j)
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw EMPTY("SwapRows",CLASS_NAME);}
           #endif
           int kmax=static_cast<int>(N2())-1;
           for(int k=0;k<=kmax;k++){ std::swap((*this)(i,k),(*this)(j,k));}
           return (*this);
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::SwapColumns(int i,int j)
         { int kmax=static_cast<int>(N1())-1;
           for(int k=0;k<=kmax;k++){ std::swap((*this)(k,i),(*this)(k,j));}
           return (*this);
          }

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::SwapColumns(int i,int j)
         { int kmax=static_cast<int>(N1())-1;
           for(int k=0;k<=kmax;k++){ std::swap((*this)(k,i),(*this)(k,j));}
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& 
MATRIX<Type,0,0>::SwapColumns(int i,int j)
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw EMPTY("SwapColumns",CLASS_NAME);}
           #endif
           int kmax=static_cast<int>(N1())-1;
           for(int k=0;k<=kmax;k++){ std::swap((*this)(k,i),(*this)(k,j));}
           return (*this);
          }
//*************************************************************************

template <typename Type,std::size_t n1,std::size_t n2> inline std::size_t MATRIX<Type,n1,n2>::N(void)
         { 
           #ifdef _LAERRORS
           if(N1()!=N2()){ throw NOT_SQUARE("N");}
           #endif 
           return N1();
          }

//*************************************************************************

template <typename Type,std::size_t n1,std::size_t n2> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator=(MATRIX<Type,n1,n2> const &M)
         { Copy(M); 
           return (*this);
          }

template <typename Type,std::size_t n1,std::size_t n2> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator=(MATRIX<Type,0,0> const &M)
         { if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator=(MATRIX<Type,0,0> const &)");} 
           Copy(M); 
           return (*this);
          }

template <typename Type,std::size_t n>
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator=(MATRIX<Type,n,n> const &M)
         { Copy(M);
           return (*this);
          }

template <typename Type,std::size_t n>
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator=(MATRIX<Type,0,0> const &M)
         { if(N()!=M.N1() || N()!=M.N2()){ throw DIFFERENT_SIZES("operator=(MATRIX<Type,0,0> const &)");} 
           Copy(M);
           return (*this);
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator=(MATRIX<Type,0,0> const &M)
         { if(Empty()==true){ Create(M.N1(),M.N2());}
           else{ if(N1()!=M.N1() || N2()!=M.N2()){ Destroy(); Create(M.N1(),M.N2());} }
           Copy(M); 
           return (*this);
          }

template <typename Type> template <std::size_t m1,std::size_t m2> MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator=(MATRIX<Type,m1,m2> const &M)
         { if(Empty()==true){ Create(M.N1(),M.N2());}
           else{ if(N1()!=M.N1() || N2()!=M.N2()){ Destroy(); Create(M.N1(),M.N2());} }
           Copy(M); 
           return (*this);
          }

//**********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator=(MATRIXEXPRESSION<Type,MType> const &M)
         { Copy(M); 
           return (*this);
          }

template <typename Type,std::size_t n> template <typename MType> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator=(MATRIXEXPRESSION<Type,MType> const &M)
         { Copy(M); 
           return (*this);
          }

template <typename Type> template <typename MType> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator=(MATRIXEXPRESSION<Type,MType> const &M)
         { if(Empty()==true){ Create(M.N1(),M.N2());}
           else{ if(N1()!=M.N1() || N2()!=M.N2()){ Destroy(); Create(M.N1(),M.N2());} }
           Copy(M); 
           return (*this);
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator+=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator+=(MATRIXEXPRESSION)");}
           #endif
           int imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1;
           for(int i=0;i<=imax;i++){ for(int j=0;j<=jmax;j++){ (*this)(i,j)+=M(i,j);} }
           return (*this);
          }

template <typename Type,std::size_t n> template <typename MType> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator+=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator+=(MATRIXEXPRESSION)");}
           #endif
           int imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1;
           for(int i=0;i<=imax;i++){ for(int j=0;j<=jmax;j++){ (*this)(i,j)+=M(i,j);} }
           return (*this);
          }

template <typename Type> template <typename MType> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator+=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator+=(MATRIXEXPRESSION)");}
           #endif
           int i,j, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i,j)+=M(i,j);} }
           return (*this);
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator-=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator-=(MATRIX)");}
           #endif
           int i,j, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i,j)-=M(i,j);} }
           return (*this);
          }

template <typename Type,std::size_t n> template <typename MType> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator-=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator-=(MATRIX)");}
           #endif
           int i,j, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i,j)-=M(i,j);} }
           return (*this);
          }

template <typename Type> template <typename MType> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator-=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N1()!=M.N1() || N2()!=M.N2()){ throw DIFFERENT_SIZES("operator-=(MATRIX)");}
           #endif
           int i,j, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i,j)-=M(i,j);} }
           return (*this);
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator*=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N2()!=M.N1()){ throw DIFFERENT_SIZES("operator*=(MATRIXEXPRESSION)");}
           #endif
           MATRIX<Type,0,0> L(N1(),M.N2());
           int i,j,k, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1,kmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ for(k=0;k<=kmax;k++){ L(i,k)+=(*this)(i,j)*M(j,k);} } }
           return (*this)=L;
          }

template <typename Type,std::size_t n> template <typename MType> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator*=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N2()!=M.N1()){ throw DIFFERENT_SIZES("operator*=(MATRIXEXPRESSION)");}
           #endif
           MATRIX<Type,0,0> L(N1(),M.N2());
           int i,j,k, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1,kmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ for(k=0;k<=kmax;k++){ L(i,k)+=(*this)(i,j)*M(j,k);} } }
           return (*this)=L;
          }

template <typename Type> template <typename MType> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator*=(MATRIXEXPRESSION<Type,MType> const &M)
         { 
           #ifdef _LAERRORS
           if(N2()!=M.N1()){ throw DIFFERENT_SIZES("operator*=(MATRIXEXPRESSION)");}
           #endif
           MATRIX<Type,0,0> L(N1(),M.N2());
           int i,j,k, imax=static_cast<int>(N1())-1,jmax=static_cast<int>(N2())-1,kmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ for(k=0;k<=kmax;k++){ L(i,k)+=(*this)(i,j)*M(j,k);} } }
           return (*this)=L;
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator/=(MATRIXEXPRESSION<Type,MType> const &M)
         { try{ (*this)*=MATRIX<Type,0,0>(M).Inverse(); return (*this);}
           catch(SINGULAR &S){ S.ChangeWhich("operator/=(MATRIX)"); throw S;}
           catch(DIFFERENT_SIZES &DS){ DS.ChangeWhich("operator/=(MATRIX)"); throw DS;}
          }

template <typename Type,std::size_t n> template <typename MType> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator/=(MATRIXEXPRESSION<Type,MType> const &M)
         { try{ (*this)*=MATRIX<Type,0,0>(M).Inverse(); return (*this);}
           catch(SINGULAR &S){ S.ChangeWhich("operator/=(MATRIX)"); throw S;}
           catch(DIFFERENT_SIZES &DS){ DS.ChangeWhich("operator/=(MATRIX)"); throw DS;}
          }

template <typename Type> template <typename MType> MATRIX<Type,0,0>& 
MATRIX<Type,0,0>::operator/=(MATRIXEXPRESSION<Type,MType> const &M)
         { try{ (*this)*=MATRIX<Type,0,0>(M).Inverse(); return (*this);}
           catch(SINGULAR &S){ S.ChangeWhich("operator/=(MATRIX)"); throw S;}
           catch(DIFFERENT_SIZES &DS){ DS.ChangeWhich("operator/=(MATRIX)"); throw DS;}
          }

//*************************************************************************

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator=(Type const &T)
         { int i,imax=static_cast<int>(N1())-1;
           for(i=0;i<=imax;i++){ (*this)(i,i)=T;} 
           return (*this);
          }

//*********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator+=(Type const &T)
         { int a, amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ (*this)(a,a)+=T;}
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& 
MATRIX<Type,0,0>::operator+=(Type const &T)
         { 
           #ifdef _LAERRORS
           if(N1()!=N2()){ throw NOT_SQUARE("operator+=(Type)");}
           #endif
           int a, amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ (*this)(a,a)+=T;}
           return (*this);
          }

//*********

template <typename Type,std::size_t n> 
MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator-=(Type const &T)
         { int a, amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ (*this)(a,a)-=T;}
           return (*this);
          }

template <typename Type> 
MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator-=(Type const &T)
         { 
           #ifdef _LAERRORS
           if(N1()!=N2()){ throw NOT_SQUARE("operator-=(Type)");}
           #endif
           int a, amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ (*this)(a,a)-=T;}
           return (*this);
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator*=(Type const &T)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)*=T;} }
           return (*this);
          }

template <typename Type,std::size_t n> MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator*=(Type const &T)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)*=T;} }
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator*=(Type const &T)
         { int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)*=T;} }
           return (*this);
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::operator/=(Type const &T)
         { 
           #ifdef _LAERRORS
           if(Equality(T,Zero(T))==true){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)/=T;} }
           return (*this);
          }

template <typename Type,std::size_t n> MATRIX<Type,n,n>& MATRIX<Type,n,n>::operator/=(Type const &T)
         { 
           #ifdef _LAERRORS
           if(Equality(T,Zero(T))==true){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)/=T;} }
           return (*this);
          }

template <typename Type> MATRIX<Type,0,0>& MATRIX<Type,0,0>::operator/=(Type const &T)
         { 
           #ifdef _LAERRORS
           if(Equality(T,Zero(T))==true){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,b, amax=static_cast<int>(N1())-1,bmax=static_cast<int>(N2())-1; 
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a,b)/=T;} }
           return (*this);
          }

//*************************************************************************

template <typename Type,std::size_t n1,std::size_t n2> 
CVECTOR<Type,n1> MATRIX<Type,n1,n2>::Column(int i) const
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N2())-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N2())-1,"Column",CLASS_NAME);}
           #endif
           CVECTOR<Type,n1> CV;
           int a,amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ CV[a]=(*this)(a,i);}

           return CV;
          }

template <typename Type,std::size_t n> 
CVECTOR<Type,n> MATRIX<Type,n,n>::Column(int i) const
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N2())-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N2())-1,"Column",CLASS_NAME);}
           #endif
           CVECTOR<Type,n> CV;
           int a,amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ CV[a]=(*this)(a,i);}

           return CV;
          }

template <typename Type> 
CVECTOR<Type,0> MATRIX<Type,0,0>::Column(int i) const
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw EMPTY("Column",CLASS_NAME);}
           if(i>static_cast<int>(N2())-1){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N2())-1,"Column",CLASS_NAME);}
           #endif
           CVECTOR<Type,0> CV(N1());
           int a,amax=static_cast<int>(N1())-1;
           for(a=0;a<=amax;a++){ CV[a]=(*this)(a,i);}

           return CV;
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> 
RVECTOR<Type,n2> MATRIX<Type,n1,n2>::Row(int i) const
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N1())-1){ throw OUT_OF_RANGE<int>(i,0,N1()-1,"Row",CLASS_NAME);}
           #endif
           RVECTOR<Type,n2> RV;
           int a, amax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ RV[a]=(*this)(i,a);}

           return RV;
          }

template <typename Type,std::size_t n> 
RVECTOR<Type,n> MATRIX<Type,n,n>::Row(int i) const
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N1())-1){ throw OUT_OF_RANGE<int>(i,0,N1()-1,"Row",CLASS_NAME);}
           #endif
           RVECTOR<Type,n> RV;
           int a, amax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ RV[a]=(*this)(i,a);}

           return RV;
          }

template <typename Type> 
RVECTOR<Type,0> MATRIX<Type,0,0>::Row(int i) const
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw EMPTY("Row",CLASS_NAME);}
           if(i>static_cast<int>(N1())-1){ throw OUT_OF_RANGE<int>(i,0,N1()-1,"Row",CLASS_NAME);}
           #endif
           RVECTOR<Type,0> RV(N2());
           int a, amax=static_cast<int>(N2())-1;
           for(a=0;a<=amax;a++){ RV[a]=(*this)(i,a);}

           return RV;
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> 
void MATRIX<Type,n1,n2>::SetSubMatrix(MATRIXEXPRESSION<Type,MType> const &M,int a0,int b0)
         { int a,b, amax=static_cast<int>(M.N1())-1,bmax=static_cast<int>(M.N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a0+a,b0+b)=M(a,b);} } 
          }

template <typename Type,std::size_t n> template <typename MType> 
void MATRIX<Type,n,n>::SetSubMatrix(MATRIXEXPRESSION<Type,MType> const &M,int a0,int b0)
         { int a,b, amax=static_cast<int>(M.N1())-1,bmax=static_cast<int>(M.N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a0+a,b0+b)=M(a,b);} } 
          }

template <typename Type> template <typename MType> 
void MATRIX<Type,0,0>::SetSubMatrix(MATRIXEXPRESSION<Type,MType> const &M,int a0,int b0)
         { int a,b, amax=static_cast<int>(M.N1())-1,bmax=static_cast<int>(M.N2())-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ (*this)(a0+a,b0+b)=M(a,b);} } 
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <class RVType> 
void MATRIX<Type,n1,n2>::SetSubRow(RVECTOREXPRESSION<Type,RVType> const &RV,int a0,int b0)
         { int b,bmax=static_cast<int>(RV.N())-1;
           for(b=0;b<=bmax;b++){ (*this)(a0,b0+b)=RV[b];} 
          }

template <typename Type,std::size_t n> template <class RVType> 
void MATRIX<Type,n,n>::SetSubRow(RVECTOREXPRESSION<Type,RVType> const &RV,int a0,int b0)
         { int b,bmax=static_cast<int>(RV.N())-1;
           for(b=0;b<=bmax;b++){ (*this)(a0,b0+b)=RV[b];} 
          }

template <typename Type> template <class RVType> 
void MATRIX<Type,0,0>::SetSubRow(RVECTOREXPRESSION<Type,RVType> const &RV,int a0,int b0)
         { int b,bmax=static_cast<int>(RV.N())-1;
           for(b=0;b<=bmax;b++){ (*this)(a0,b0+b)=RV[b];} 
          }

//********

template <typename Type,std::size_t n1,std::size_t n2> template <class CVType> 
void MATRIX<Type,n1,n2>::SetSubColumn(CVECTOREXPRESSION<Type,CVType> const &CV,int a0,int b0)
         { int a, amax=static_cast<int>(CV.N())-1;
           for(a=0;a<=amax;a++){ (*this)(a0+a,b0)=CV[a];} 
          }

template <typename Type,std::size_t n> template <class CVType> 
void MATRIX<Type,n,n>::SetSubColumn(CVECTOREXPRESSION<Type,CVType> const &CV,int a0,int b0)
         { int a, amax=static_cast<int>(CV.N())-1;
           for(a=0;a<=amax;a++){ (*this)(a0+a,b0)=CV[a];} 
          }

template <typename Type> template <class CVType> 
void MATRIX<Type,0,0>::SetSubColumn(CVECTOREXPRESSION<Type,CVType> const &CV,int a0,int b0)
         { int a, amax=static_cast<int>(CV.N())-1;
           for(a=0;a<=amax;a++){ (*this)(a0+a,b0)=CV[a];} 
          }

//*****************************

template <typename Type,std::size_t n1,std::size_t n2> 
inline std::array<Type,n2>& MATRIX<Type,n1,n2>::operator[](int i)
         { return m[i];}

template <typename Type,std::size_t n> 
inline std::array<Type,n>& MATRIX<Type,n,n>::operator[](int i)
         { return m[i];}

//template <typename Type> inline Type*& MATRIX<Type,0,0>::operator[](int i) { return m[i];}
template <typename Type> inline std::vector<Type>& MATRIX<Type,0,0>::operator[](int i) { return m[i];}

//*******

template <typename Type,std::size_t n1,std::size_t n2> 
inline std::array<Type,n2> MATRIX<Type,n1,n2>::operator[](int i) const
         { return m[i];}

template <typename Type,std::size_t n> 
inline std::array<Type,n> MATRIX<Type,n,n>::operator[](int i) const
         { return m[i];}

//template <typename Type> inline const Type* MATRIX<Type,0,0>::operator[](int i) const { return m[i];}
template <typename Type> inline std::vector<Type> MATRIX<Type,0,0>::operator[](int i) const { return m[i];}

//*****************************

template <typename Type,std::size_t n1,std::size_t n2> inline Type& MATRIX<Type,n1,n2>::operator()(int i,int j)
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N1())-1 || 0>i){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N1())-1,"Type& (int,int)",CLASS_NAME);}
           if(j>static_cast<int>(N2())-1 || 0>j){ throw OUT_OF_RANGE<int>(j,0,static_cast<int>(N2())-1,"Type& (int,int)]",CLASS_NAME);}
           #endif
           return m[i][j];
          }

template <typename Type,std::size_t n> inline Type& MATRIX<Type,n,n>::operator()(int i,int j)
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N1())-1 || 0>i){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N1())-1,"Type& (int,int)",CLASS_NAME);}
           if(j>static_cast<int>(N2())-1 || 0>j){ throw OUT_OF_RANGE<int>(j,0,static_cast<int>(N2())-1,"Type& (int,int)]",CLASS_NAME);}
           #endif
           return m[i][j];
          }

template <typename Type> inline Type& MATRIX<Type,0,0>::operator()(int i,int j)
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw EMPTY("Type& (int,int)",CLASS_NAME);}
           if(i>static_cast<int>(N1())-1 || 0>i){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N1())-1,"Type& (int,int)",CLASS_NAME);}
           if(j>static_cast<int>(N2())-1 || 0>j){ throw OUT_OF_RANGE<int>(j,0,static_cast<int>(N2())-1,"Type& (int,int)]",CLASS_NAME);}
           #endif
           return m[i][j];
          }

//*******

template <typename Type,std::size_t n1,std::size_t n2> inline Type MATRIX<Type,n1,n2>::operator()(int i,int j) const
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N1())-1 || 0>i){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N1())-1,"Type (int,int)",CLASS_NAME);}
           if(j>static_cast<int>(N2())-1 || 0>j){ throw OUT_OF_RANGE<int>(j,0,static_cast<int>(N2())-1,"Type (int,int)",CLASS_NAME);}
           #endif
           return m[i][j];
          }

template <typename Type,std::size_t n> inline Type MATRIX<Type,n,n>::operator()(int i,int j) const
         { 
           #ifdef _LAERRORS
           if(i>static_cast<int>(N1())-1 || 0>i){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N1())-1,"Type (int,int)",CLASS_NAME);}
           if(j>static_cast<int>(N2())-1 || 0>j){ throw OUT_OF_RANGE<int>(j,0,static_cast<int>(N2())-1,"Type (int,int)",CLASS_NAME);}
           #endif
           return m[i][j];
          }

template <typename Type> inline Type MATRIX<Type,0,0>::operator()(int i,int j) const
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw EMPTY("Type (int,int)",CLASS_NAME);}
           if(i>static_cast<int>(N1())-1 || 0>i){ throw OUT_OF_RANGE<int>(i,0,static_cast<int>(N1())-1,"Type (int,int)",CLASS_NAME);}
           if(j>static_cast<int>(N2())-1 || 0>j){ throw OUT_OF_RANGE<int>(j,0,static_cast<int>(N2())-1,"Type (int,int)",CLASS_NAME);}
           #endif
           return m[i][j];
          }

//*****************************
         
template <typename Type,std::size_t n1,std::size_t n2> template <typename MType>
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::SubMatrixAddition(MATRIXEXPRESSION<Type,MType> const &M,int i0,int j0)
         { 
           #ifdef _LAERRORS
           if(i0+M.N1()>N1()){ throw INCORRECT_FORM("SubMatrixAddition");}
           if(j0+M.N2()>N2()){ throw INCORRECT_FORM("SubMatrixAddition");}
           #endif
           int i,j, imax=static_cast<int>(M.N1())-1,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i0+i,j0+j)+=M(i,j);} }
           return (*this);
          }

template <typename Type,std::size_t n> template <typename MType>
MATRIX<Type,n,n>& MATRIX<Type,n,n>::SubMatrixAddition(MATRIXEXPRESSION<Type,MType> const &M,int i0,int j0)
         { 
           #ifdef _LAERRORS
           if(i0+M.N1()>N1()){ throw INCORRECT_FORM("SubMatrixAddition");}
           if(j0+M.N2()>N2()){ throw INCORRECT_FORM("SubMatrixAddition");}
           #endif
           int i,j, imax=static_cast<int>(M.N1())-1,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i0+i,j0+j)+=M(i,j);} }
           return (*this);
          }

template <typename Type> template <typename MType>
MATRIX<Type,0,0>& MATRIX<Type,0,0>::SubMatrixAddition(MATRIXEXPRESSION<Type,MType> const &M,int i0,int j0)
         { 
           #ifdef _LAERRORS
           if(i0+M.N1()>N1()){ throw INCORRECT_FORM("SubMatrixAddition");}
           if(j0+M.N2()>N2()){ throw INCORRECT_FORM("SubMatrixAddition");}
           #endif
           int i,j, imax=static_cast<int>(M.N1())-1,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i0+i,j0+j)+=M(i,j);} }
           return (*this);
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType>
MATRIX<Type,n1,n2>& MATRIX<Type,n1,n2>::SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType> const &M,int i0,int j0)
         { 
           #ifdef _LAERRORS
           if(i0+M.N1()>N1()){ throw INCORRECT_FORM("SubMatrixSubtraction");}
           if(j0+M.N2()>N2()){ throw INCORRECT_FORM("SubMatrixSubtraction");}
           #endif
           int i,j, imax=static_cast<int>(M.N1())-1,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i0+i,j0+j)-=M(i,j);} }
           return (*this);
          }

template <typename Type,std::size_t n> template <typename MType>
MATRIX<Type,n,n>& MATRIX<Type,n,n>::SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType> const &M,int i0,int j0)
         { 
           #ifdef _LAERRORS
           if(i0+M.N1()>N1()){ throw INCORRECT_FORM("SubMatrixSubtraction");}
           if(j0+M.N2()>N2()){ throw INCORRECT_FORM("SubMatrixSubtraction");}
           #endif
           int i,j, imax=static_cast<int>(M.N1())-1,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i0+i,j0+j)-=M(i,j);} }
           return (*this);
          }

template <typename Type> template <typename MType>
MATRIX<Type,0,0>& MATRIX<Type,0,0>::SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType> const &M,int i0,int j0)
         { 
           #ifdef _LAERRORS
           if(i0+M.N1()>N1()){ throw INCORRECT_FORM("SubMatrixSubtraction");}
           if(j0+M.N2()>N2()){ throw INCORRECT_FORM("SubMatrixSubtraction");}
           #endif
           int i,j, imax=static_cast<int>(M.N1())-1,jmax=static_cast<int>(M.N2())-1;
           for(i=0;i<=imax;i++){ for(j=0;j<=jmax;j++){ (*this)(i0+i,j0+j)-=M(i,j);} }
           return (*this);
          }

//*************************************************************************

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> bool MATRIX<Type,n1,n2>::operator==(MATRIXEXPRESSION<Type,MType> const &M) const
         { bool equal=true;
           if(N1()!=M.N1() || N2()!=M.N2()){ equal=false;}
           else{ int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ if(Equality((*this)(a,b),M(a,b))==false){ equal=false;}}}
                }
           return equal;
          }

template <typename Type,std::size_t n> template <typename MType> bool MATRIX<Type,n,n>::operator==(MATRIXEXPRESSION<Type,MType> const &M) const
         { bool equal=true;
           if(N1()!=M.N1() || N2()!=M.N2()){ equal=false;}
           else{ int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ if(Equality((*this)(a,b),M(a,b))==false){ equal=false;}}}
                }
           return equal;
          }

template <typename Type> template <typename MType> bool MATRIX<Type,0,0>::operator==(MATRIXEXPRESSION<Type,MType> const &M) const
         { bool equal=true;
           if(N1()!=M.N1() || N2()!=M.N2()){ equal=false;}
           else{ int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ if(Equality((*this)(a,b),M(a,b))==false){ equal=false;}}}
                }
           return equal;
          }

//*********

template <typename Type,std::size_t n1,std::size_t n2> template <typename MType> bool MATRIX<Type,n1,n2>::operator!=(MATRIXEXPRESSION<Type,MType> const &M) const
         { bool notequal=false;
           if(N1()!=M.N1() || N2()!=M.N2()){ notequal=true;}
           else{ int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ if(Equality((*this)(a,b),M(a,b))==false){ notequal=true;}}}
                }
           return notequal;
          }

template <typename Type,std::size_t n> template <typename MType> bool MATRIX<Type,n,n>::operator!=(MATRIXEXPRESSION<Type,MType> const &M) const
         { bool notequal=false;
           if(N1()!=M.N1() || N2()!=M.N2()){ notequal=true;}
           else{ int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ if(Equality((*this)(a,b),M(a,b))==false){ notequal=true;}}}
                }
           return notequal;
          }

template <typename Type> template <typename MType> bool MATRIX<Type,0,0>::operator!=(MATRIXEXPRESSION<Type,MType> const &M) const
         { bool notequal=false;
           if(N1()!=M.N1() || N2()!=M.N2()){ notequal=true;}
           else{ int a,amax=static_cast<int>(M.N1())-1, b,bmax=static_cast<int>(M.N2())-1;
                 for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ if(Equality((*this)(a,b),M(a,b))==false){ notequal=true;}}}
                }
           return notequal;
          }

#endif

