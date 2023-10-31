#if !defined(_MATH2_FUNCTIONS)
#define _MATH2_FUNCTIONS

template <class Type> Type Theta(const Type &T){ if(T<=Zero<Type>()){ return Zero<Type>();} else{ return One<Type>();} }

//*****************************************************************************************
//*****************************************************************************************
//*****************************************************************************************

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(Type A,Type B,unsigned int N)
        { int i=0; 
          if(A>B){ while(i<static_cast<int>(N)){ B=nextafter(B,+INFINITY); i++;}; return B >= A;}
          else{ while(i<static_cast<int>(N)){ A=nextafter(A,+INFINITY); i++;}; return A >= B;}
         }

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(Type A,std::complex<Type> B,unsigned int N)
         { return Equality(A,real(B),N) & Equality(Zero(A),imag(B),N);}

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(std::complex<Type> A,Type B,unsigned int N)
         { return Equality(real(A),B,N) & Equality(imag(A),Zero(B),N);}

template<class Type> typename std::enable_if<!std::numeric_limits<Type>::is_integer, bool>::type Equality(std::complex<Type> A,std::complex<Type> B,unsigned int N)
         { return Equality(real(A),real(B),N) & Equality(imag(A),imag(B),N);}

#endif



