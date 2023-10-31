#if !defined(_CRVECTOR_MEMBERS)
#define _CRVECTOR_MEMBERS

template <class Type,std::size_t n> std::string CVECTOR<Type,n>::CLASS_NAME("CVECTOR - defined size");
template <class Type> std::string CVECTOR<Type,0>::CLASS_NAME("CVECTOR - undefined size");

// *******************************************************
// *********** PUBLIC CVECTOR MEMBERS ****************
// *******************************************************

template <typename Type,std::size_t n> CVECTOR<Type,n>::CVECTOR(Type CV0,...)
      	{ va_list val; va_start(val,CV0);
         (*this)[0]=CV0; 
         int a,amax=static_cast<int>(Size())-1;
         for(a=1;a<=amax;a++){ (*this)[a]=va_arg(val,Type);}
      	  va_end(val);
	 }

template <typename Type> CVECTOR<Type,0>::CVECTOR(const std::size_t N,Type CV0,...) : v(N)
      	{ va_list val; va_start(val,CV0);
         (*this)[0]=CV0; 
         int a,amax=static_cast<int>(Size())-1;
         for(a=1;a<=amax;a++){ (*this)[a]=va_arg(val,Type);}
      	  va_end(val);
	 }

// ******

template <typename Type,std::size_t n> CVECTOR<Type,n>::CVECTOR(std::vector<Type> const &V) : v(n,Zero<Type>())
         { if(n!=V.size()){ throw DIFFERENT_SIZES("CVECTOR(std::vector<Type> const&)");}
           v=V;
          }

// ******

template <typename Type,std::size_t n> CVECTOR<Type,n>::CVECTOR(CVECTOR<Type,0> const &CV)
         { if(n!=CV.Size()){ throw DIFFERENT_SIZES("CVECTOR(CVECTOR<Type,0> const&)");}
           v=CV.v;
          }

// ******

template <typename Type,std::size_t n> template <typename CVType> CVECTOR<Type,n>::CVECTOR(CVECTOREXPRESSION<Type,CVType> const &CV) : v(n,Zero<Type>())
        { if(n!=CV.Size()){ throw DIFFERENT_SIZES("CVECTOREXPRESSION<Type,CVType> const &");}
          int a,amax=static_cast<int>(Size())-1;
          for(a=0;a<=amax;a++){ (*this)[a]=CV[a];} 
         }

template <typename Type> template <typename CVType> CVECTOR<Type,0>::CVECTOR(CVECTOREXPRESSION<Type,CVType> const &CV) : v(CV.Size())
        { int a,amax=static_cast<int>(Size())-1;
          for(a=0;a<=amax;a++){ (*this)[a]=CV[a];} 
         }

// *******************************************************

template <typename Type,std::size_t n> CVECTOR<Type,n>& CVECTOR<Type,n>::operator=(std::vector<Type> const &V)
         { if(n!=V.size()){ throw DIFFERENT_SIZES("operator=(std::vector<Type> const&)");}
           v=V; 
           return (*this);
          }

template <typename Type> CVECTOR<Type,0>& CVECTOR<Type,0>::operator=(std::vector<Type> const& V)
         { v=V; 
           return (*this);
          }

// *******

template <typename Type,std::size_t n> CVECTOR<Type,n>& CVECTOR<Type,n>::operator=(CVECTOR<Type,n> const &CV)
         { v=CV.v; 
           return (*this);
          }

template <typename Type,std::size_t n> CVECTOR<Type,n>& CVECTOR<Type,n>::operator=(CVECTOR<Type,0> const &CV)
         { if(n!=CV.Size()){ throw DIFFERENT_SIZES("operator=(CVECTOR<Type,0> const&)");}
           v=CV.v; 
           return (*this);
          }

template <typename Type> CVECTOR<Type,0>& CVECTOR<Type,0>::operator=(CVECTOR<Type,0> const& CV)
         { v=CV.v; 
           return (*this);
          }

template <typename Type> template <std::size_t n> CVECTOR<Type,0>& CVECTOR<Type,0>::operator=(CVECTOR<Type,n> const& CV)
         { v=CV.v; 
           return (*this);
          }

// *******

template <typename Type,std::size_t n> template <typename CVType> CVECTOR<Type,n>& CVECTOR<Type,n>::operator=(CVECTOREXPRESSION<Type,CVType> const &CV)
         { if(n!=CV.Size()){ throw DIFFERENT_SIZES("operator=(CVECTOREXPRESSION<Type,CVType> const&)");}
           v=std::vector<Type>(CV.Size()); 
           int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]=CV[a];} 
           return (*this);
          }

template <typename Type> template <typename CVType> CVECTOR<Type,0>& CVECTOR<Type,0>::operator=(CVECTOREXPRESSION<Type,CVType> const &CV)
         { v=std::vector<Type>(CV.Size()); 
           int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]=CV[a];} 
           return (*this);
          }

// *******

template <typename Type> CVECTOR<Type,0>& CVECTOR<Type,0>::operator=(Type const& T)
         { v=std::vector<Type>(1,T); return (*this);}

// *******

template <typename Type,std::size_t n> CVECTOR<Type,n>& CVECTOR<Type,n>::operator*=(Type const& T)
         { int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]*=T;}
           return (*this);
          }

template <typename Type> CVECTOR<Type,0>& CVECTOR<Type,0>::operator*=(Type const& T)
         { int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]*=T;}
           return (*this);
          }

// *******

template <typename Type,std::size_t n> CVECTOR<Type,n>& CVECTOR<Type,n>::operator/=(Type const& T)
         {  
           #ifdef _LAERRORS
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]/=T;}
           return (*this);
          }

template <typename Type> CVECTOR<Type,0>& CVECTOR<Type,0>::operator/=(Type const& T)
         {  
           #ifdef _LAERRORS
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]/=T;}
           return (*this);
          }

// *******

template <typename Type,std::size_t n> template <typename CVType> CVECTOR<Type,n>& CVECTOR<Type,n>::operator+=(CVECTOREXPRESSION<Type,CVType> const &CV)
      	 {  
           #ifdef _LAERRORS
           if(Size()!=CV.Size()){ throw DIFFERENT_SIZES("operator+=(CVECTOREXPRESSION)",CLASS_NAME);}
           #endif
           int a;
           for(a=0;a<=static_cast<int>(Size())-1;a++){ (*this)[a]+=CV[a];}
           return (*this);
          }

template <typename Type> template <typename CVType> CVECTOR<Type,0>& CVECTOR<Type,0>::operator+=(CVECTOREXPRESSION<Type,CVType> const &CV)
      	 {  
           #ifdef _LAERRORS
           if(Size()!=CV.Size()){ throw DIFFERENT_SIZES("operator+=(CVECTOREXPRESSION)",CLASS_NAME);}
           #endif
           int a;
           for(a=0;a<=static_cast<int>(Size())-1;a++){ (*this)[a]+=CV[a];}
           return (*this);
          }

// *******

template <typename Type,std::size_t n> template <typename CVType> CVECTOR<Type,n>& CVECTOR<Type,n>::operator-=(CVECTOREXPRESSION<Type,CVType> const &CV)
      	 {  
           #ifdef _LAERRORS
           if(Size()!=CV.Size()){ throw DIFFERENT_SIZES("operator-=(CVECTOREXPRESSION)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]-=CV[a];}
           return (*this);
          }

template <typename Type> template <typename CVType> CVECTOR<Type,0>& CVECTOR<Type,0>::operator-=(CVECTOREXPRESSION<Type,CVType> const &CV)
      	 {  
           #ifdef _LAERRORS
           if(Size()!=CV.Size()){ throw DIFFERENT_SIZES("operator-=(CVECTOREXPRESSION)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;
           for(a=0;a<=amax;a++){ (*this)[a]-=CV[a];}
           return (*this);
          }

//*************************************************************************

template <typename Type,typename MType,typename CVType> Type CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION<Type,MType,CVType>::operator[](int i) const
         { Type T=Zero<Type>();
           int j,jmax=static_cast<int>(Size())-1;
           for(j=0;j<=jmax;j++){ T+=m(i,j)*cv[j];}
           return T;
          }

//*************************************************************************
//*************************************************************************
//*************************************************************************

template <class Type,std::size_t n> std::string RVECTOR<Type,n>::CLASS_NAME("RVECTOR - defined size");
template <class Type> std::string RVECTOR<Type,0>::CLASS_NAME("RVECTOR - undefined size");

// *******************************************************
// ************* PUBLIC RVECTOR MEMBERS ******************
// *******************************************************

template <typename Type,std::size_t n> RVECTOR<Type,n>::RVECTOR(Type RV0,...) : v(n)
         { va_list val; va_start(val,RV0);
           (*this)[0]=RV0;
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]=va_arg(val,Type);}
           va_end(val);
	  }

template <typename Type> RVECTOR<Type,0>::RVECTOR(const std::size_t N,Type RV0,...) : v(N)
         { va_list val; va_start(val,RV0);
           (*this)[0]=RV0;
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]=va_arg(val,Type);}
           va_end(val);
	  }

// ******

template <typename Type,std::size_t n> RVECTOR<Type,n>::RVECTOR(std::vector<Type> const &V) : v(n,Zero<Type>())
         { if(n!=V.size()){ throw DIFFERENT_SIZES("RVECTOR(std::vector<Type> const &)");}
           v=V;
          }

// ******

template <typename Type,std::size_t n> RVECTOR<Type,n>::RVECTOR(RVECTOR<Type,0> const &RV)
         { if(n!=RV.Size()){ throw DIFFERENT_SIZES("RVECTOR(RVECTOR<Type,0> const &)");}
           v=RV.v;
          }

// ******

template <typename Type,std::size_t n> template <typename RVType> RVECTOR<Type,n>::RVECTOR(RVECTOREXPRESSION<Type,RVType> const &RV) : v(n,Zero<Type>())
      	{ if(n!=RV.Size()){ throw DIFFERENT_SIZES("RVECTOR(RVECTOREXPRESSION<Type,RVType> const &)");}           
          int a,amax=static_cast<int>(Size())-1;  
          for(a=1;a<=amax;a++){ (*this)[a]=RV[a];} 
         }

template <typename Type> template <typename RVType> RVECTOR<Type,0>::RVECTOR(RVECTOREXPRESSION<Type,RVType> const&RV) : v(RV.Size())
      	{ int a,amax=static_cast<int>(Size())-1;  
          for(a=1;a<=amax;a++){ (*this)[a]=RV[a];} 
         }

// *******************************************************

template <typename Type,std::size_t n> RVECTOR<Type,n>& RVECTOR<Type,n>::operator=(std::vector<Type> const &V)
         { if(n!=V.size()){ throw DIFFERENT_SIZES("operator=(std::vector<Type> const &)");}
           v=V;
           return (*this);
          }

template <typename Type> RVECTOR<Type,0>& RVECTOR<Type,0>::operator=(std::vector<Type> const &V)
         { v=V;
           return (*this);
          }

// *******

template <typename Type,std::size_t n> RVECTOR<Type,n>& RVECTOR<Type,n>::operator=(RVECTOR<Type,n> const &RV)
         { v=RV.v;
           return (*this);
          }

template <typename Type,std::size_t n> RVECTOR<Type,n>& RVECTOR<Type,n>::operator=(RVECTOR<Type,0> const &RV)
         { if(n!=RV.Size()){ throw DIFFERENT_SIZES("operator=(RVECTOR<Type,0> const &)");}
           v=RV.v;
           return (*this);
          }

template <typename Type> template <std::size_t n> RVECTOR<Type,0>& RVECTOR<Type,0>::operator=(RVECTOR<Type,n> const &RV)
         { v=RV.v;
           return (*this);
          }

template <typename Type> RVECTOR<Type,0>& RVECTOR<Type,0>::operator=(RVECTOR<Type,0> const &RV)
         { v=RV.v;
           return (*this);
          }

// *******

template <typename Type,std::size_t n> template <typename RVType> RVECTOR<Type,n>& RVECTOR<Type,n>::operator=(RVECTOREXPRESSION<Type,RVType> const &RV)
         { if(n!=RV.Size()){ throw DIFFERENT_SIZES("operator=(RVECTOREXPRESSION<Type,RVType> const &)");}
           v=std::vector<Type>(RV.Size()); 
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]=RV[a];} 
           return (*this);
          }

template <typename Type> template <typename RVType> RVECTOR<Type,0>& RVECTOR<Type,0>::operator=(RVECTOREXPRESSION<Type,RVType> const &RV)
         { v=std::vector<Type>(RV.Size()); 
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]=RV[a];} 
           return (*this);
          }

// *******

template <typename Type> RVECTOR<Type,0>& RVECTOR<Type,0>::operator=(Type const&T)
         { v=std::vector<Type>(1,T); return (*this);}

// *******

template <typename Type,std::size_t n> RVECTOR<Type,n>& RVECTOR<Type,n>::operator*=(Type const&T)
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw Empty("operator*=(Type)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]*=T;}
           return (*this);
          }

template <typename Type> RVECTOR<Type,0>& RVECTOR<Type,0>::operator*=(Type const&T)
         { 
           #ifdef _LAERRORS
           if(Empty()==true){ throw Empty("operator*=(Type)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]*=T;}
           return (*this);
          }

// *******

template <typename Type,std::size_t n> RVECTOR<Type,n>& RVECTOR<Type,n>::operator/=(Type const&T)
         {  
           #ifdef _LAERRORS
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]/=T;}
           return (*this);
          }

template <typename Type> RVECTOR<Type,0>& RVECTOR<Type,0>::operator/=(Type const&T)
         {  
           #ifdef _LAERRORS
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/=(Type)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]/=T;}
           return (*this);
          }

// *******

template <typename Type> template <typename MType> RVECTOR<Type>& RVECTOR<Type>::operator*=(MATRIXEXPRESSION<Type,MType> const&MM)
         {  
           #ifdef _LAERRORS
           if(Size()!=MM.N1()){ throw DIFFERENT_SIZES("operator*=(MATRIX)",CLASS_NAME);}
           #endif
           RVECTOR<Type> MRV2(MM.N2());
           int a,amax=static_cast<int>(Size())-1,b,bmax=MM.N2()-1;
           for(a=0;a<=amax;a++){ for(b=0;b<=bmax;b++){ MRV2[b]+=(*this)[a]*MM(a,b);} }
	   return (*this)=MRV2;
          }

template <typename Type> template <typename MType> RVECTOR<Type>& RVECTOR<Type>::operator/=(MATRIXEXPRESSION<Type,MType> const&M)
      	 { try{ return (*this)*=Invert(M);}
           catch(SINGULAR S){ S.ChangeFunction("operator/=(MATRIX)"); throw S;}
           catch(DIFFERENT_SIZES DS){ DS.ChangeFunction("operator/=(MATRIX)"); throw DS;}
	  }

template <typename Type> template <typename RVType> RVECTOR<Type>& RVECTOR<Type>::operator+=(RVECTOREXPRESSION<Type,RVType> const&RV)
         {  
           #ifdef _LAERRORS
           if(Size()!=RV.Size()){ throw DIFFERENT_SIZES("operator+=(RVECTOR)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]+=RV[a];}
           return (*this);
          }

template <typename Type> template <typename RVType> RVECTOR<Type>& RVECTOR<Type>::operator-=(RVECTOREXPRESSION<Type,RVType> const&RV)
      	 {  
           #ifdef _LAERRORS
           if(Size()!=RV.Size()){ throw DIFFERENT_SIZES("operator-=(RVECTOR)",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1;  
           for(a=1;a<=amax;a++){ (*this)[a]-=RV[a];}
           return (*this);
          }

//*************************************************************************

template <typename Type,typename RVType,typename MType> Type RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION<Type,RVType,MType>::operator[](int i) const
         { Type T=Zero<Type>();
           int j,jmax=static_cast<int>(Size())-1;
           for(j=0;j<=jmax;j++){ T+=rv[j]*m(j,i);}
           return T;
          }

#endif



