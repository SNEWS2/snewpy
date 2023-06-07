#if !defined(_MVECTOR_MEMBERS)
#define _MVECTOR_MEMBERS

template <typename Type,std::size_t n> std::string MVECTOR<Type,n>::CLASS_NAME("MVECTOR - defined size");
template <typename Type> std::string MVECTOR<Type,0>::CLASS_NAME("MVECTOR - defined size");

// ******** MVECTOR OVERLOADED OPERATORS *****************

template <typename Type,std::size_t n> 
MVECTOR<Type,n>& MVECTOR<Type,n>::operator=(std::vector<Type> const &V)
  	      { if(Size()!=V.size()){ throw DIFFERENT_LENGTHS("operator=vector",CLASS_NAME);}
                v=V;
	        return (*this);
	       }

template <typename Type> 
MVECTOR<Type,0>& MVECTOR<Type,0>::operator=(std::vector<Type> const &V)
  	      { v=V;
	        return (*this);
	       }

// ******** 

template <typename Type,std::size_t n> 
MVECTOR<Type,n>& MVECTOR<Type,n>::operator=(MVECTOR<Type,n> const &V)
  	      { v=V.v;
	        return (*this);
	       }

template <typename Type,std::size_t n> 
MVECTOR<Type,n>& MVECTOR<Type,n>::operator=(MVECTOR<Type,0> const &V)
  	      { if(Size()!=V.Size()){ throw DIFFERENT_LENGTHS("operator=VECTOR",CLASS_NAME);}
                v=V.v;
	        return (*this);
	       }

template <typename Type> template <std::size_t n>
MVECTOR<Type,0>& MVECTOR<Type,0>::operator=(MVECTOR<Type,n> const &V)
  	      { v=V.v;
	        return (*this);
	       }

template <typename Type> 
MVECTOR<Type,0>& MVECTOR<Type,0>::operator=(MVECTOR<Type,0> const &V)
  	      { v=V.v;
	        return (*this);
	       }

// ******** 

template <typename Type,std::size_t n> template <typename VEType>
MVECTOR<Type,n>& MVECTOR<Type,n>::operator=(MVECTOREXPRESSION<Type,VEType> const &VE)
  	      { if(Size()!=VE.Size()){ throw DIFFERENT_LENGTHS("operator=VECTOREXPRESSION",CLASS_NAME);}
                v=std::vector<Type>(VE.Size());
	        int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]=VE[a];}
	        return (*this);
	       }

template <typename Type> template <typename VEType>
MVECTOR<Type,0>& MVECTOR<Type,0>::operator=(MVECTOREXPRESSION<Type,VEType> const &VE)
  	      { v=std::vector<Type>(VE.Size());
	        int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]=VE[a];}
	        return (*this);
	       }

// *******************************************************

template <typename Type,std::size_t n> template <typename VEType> MVECTOR<Type,n>& MVECTOR<Type,n>::operator+=(MVECTOREXPRESSION<Type,VEType> const &VE)
  	      { if(Size()!=VE.Size()){ throw DIFFERENT_LENGTHS("operator+=(MVECTOREXPRESSION<Type,VEType> const &)",CLASS_NAME);}
	        int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]+=VE[a];}
	        return (*this);
	       }

template <typename Type> template <typename VEType> MVECTOR<Type,0>& MVECTOR<Type,0>::operator+=(MVECTOREXPRESSION<Type,VEType> const &VE)
  	      { if(Size()!=VE.Size()){ throw DIFFERENT_LENGTHS("operator+=(MVECTOREXPRESSION<Type,VEType> const &)",CLASS_NAME);}
	        int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]+=VE[a];}
	        return (*this);
	       }

template <typename Type,std::size_t n> template <typename VEType> MVECTOR<Type,n>& MVECTOR<Type,n>::operator-=(MVECTOREXPRESSION<Type,VEType> const &VE)
  	      { if(Size()!=VE.Size()){ throw DIFFERENT_LENGTHS("operator-=(MVECTOREXPRESSION<Type,VEType> const &)",CLASS_NAME);}
                int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]-=VE[a];}
	        return (*this);
	       }

template <typename Type> template <typename VEType> MVECTOR<Type,0>& MVECTOR<Type,0>::operator-=(MVECTOREXPRESSION<Type,VEType> const &VE)
  	      { if(Size()!=VE.Size()){ throw DIFFERENT_LENGTHS("operator-=(MVECTOREXPRESSION<Type,VEType> const &)",CLASS_NAME);}
                int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]-=VE[a];}
	        return (*this);
	       }

// *******************************************************

template <typename Type,std::size_t n> MVECTOR<Type,n>& MVECTOR<Type,n>::operator*=(Type const &T)
 	      { int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]*=T;} 
                return (*this);
               }

template <typename Type> MVECTOR<Type,0>& MVECTOR<Type,0>::operator*=(Type const&T)
 	      { int a,amax=static_cast<int>(Size())-1; 
                for(a=0;a<=amax;a++){ (*this)[a]*=T;} 
                return (*this);
               }

template <typename Type,std::size_t n> MVECTOR<Type,n>& MVECTOR<Type,n>::operator/=(Type const &T)
         { 
           #ifdef _LAERRORS
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/=",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1; 
           for(a=0;a<=amax;a++){ (*this)[a]/=T;}
           return (*this);
          }

template <typename Type> MVECTOR<Type,0>& MVECTOR<Type,0>::operator/=(Type const&T)
         { 
           #ifdef _LAERRORS
           if(T==Zero<Type>()){ throw DIVISION_BY_ZERO("operator/=",CLASS_NAME);}
           #endif
           int a,amax=static_cast<int>(Size())-1; 
           for(a=0;a<=amax;a++){ (*this)[a]/=T;}
           return (*this);
          }

// ************** PUBLIC MVECTOR MEMBERS *****************

template <typename Type,std::size_t n> MVECTOR<Type,n>::MVECTOR(Type V0,...) : v(n)
    	   { va_list val; va_start(val,V0);
             (*this)[0]=V0; 
             int a,amax=static_cast<int>(Size())-1; 
             for(a=0;a<=amax;a++){ (*this)[a]=va_arg(val,Type);}
	     va_end(val);
	    }

template <typename Type> MVECTOR<Type,0>::MVECTOR(std::size_t N,Type V0,...) : v(N)
    	   { va_list val; va_start(val,V0);
             (*this)[0]=V0; 
             int a,amax=static_cast<int>(Size())-1; 
             for(a=0;a<=amax;a++){ (*this)[a]=va_arg(val,Type);}
	     va_end(val);
	    }

template <typename Type,std::size_t n> template <typename VType>
MVECTOR<Type,n>::MVECTOR(MVECTOREXPRESSION<Type,VType> const &V) : v(n)
    	   { if(Size()!=V.Size()){ throw DIFFERENT_LENGTHS("MVECTOR(MVECTOREXPRESSION<Type,VType> const &)",CLASS_NAME);}
             int a,amax=static_cast<int>(Size())-1; 
             for(a=0;a<=amax;a++){ (*this)[a]=V[a];} 
            }

template <typename Type> template <typename VType>
MVECTOR<Type,0>::MVECTOR(MVECTOREXPRESSION<Type,VType> const &V) : v(V.Size())
    	   { int a,amax=static_cast<int>(Size())-1; 
             for(a=0;a<=amax;a++){ (*this)[a]=V[a];} 
            }

#endif

