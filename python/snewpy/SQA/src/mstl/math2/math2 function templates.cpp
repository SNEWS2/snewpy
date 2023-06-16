#if !defined(_MATH2_FUNCTIONS)
#define _MATH2_FUNCTIONS

template <class Type> Type Theta(const Type &T){ if(T<=Zero<Type>()){ return Zero<Type>();} else{ return One<Type>();} }

//*****************************************************************************************
//*****************************************************************************************
//*****************************************************************************************

template <class Type> bool Equality(Type A,Type B,int N)
         { if(std::isinf(A)==true || std::isinf(B)==true){ return A == B;}
           if(std::isnan(A)==true || std::isnan(B)==true){ return A == B;}
           if(std::numeric_limits<Type>::is_integer==true){ return A == B;}

           //convert to bit patterns 
           int i; 
           //std::bitset<sizeof(A)*CHAR_BIT> a(*reinterpret_cast<unsigned long*>(&A));
           //i=1; while(i*sizeof(unsigned long)<sizeof(A)){ a|=(std::bitset<sizeof(A)*CHAR_BIT>(*(reinterpret_cast<unsigned long*>(&A)+i))<<sizeof(unsigned long)*CHAR_BIT); i++;};
           //std::bitset<sizeof(B)*CHAR_BIT> b(*reinterpret_cast<unsigned long*>(&B));
           //i=1; while(i*sizeof(unsigned long)<sizeof(B)){ b|=(std::bitset<sizeof(B)*CHAR_BIT>(*(reinterpret_cast<unsigned long*>(&B)+i))<<sizeof(unsigned long)*CHAR_BIT); i++;};
           std::bitset<sizeof(A)*CHAR_BIT> a(*reinterpret_cast<char*>(&A));
           i=1; while(i*sizeof(unsigned long)<sizeof(A)){ a|=(std::bitset<sizeof(A)*CHAR_BIT>(*(reinterpret_cast<char*>(&A)+i))<<sizeof(char)*CHAR_BIT); i++;};
           std::bitset<sizeof(B)*CHAR_BIT> b(*reinterpret_cast<char*>(&B));
           i=1; while(i*sizeof(unsigned long)<sizeof(B)){ b|=(std::bitset<sizeof(B)*CHAR_BIT>(*(reinterpret_cast<char*>(&B)+i))<<sizeof(char)*CHAR_BIT); i++;};


           if(std::numeric_limits<Type>::is_signed==true) //flip from twos complement if negative
             { if(a[a.size()-1]==true){ a.flip()++;} 
               if(b[b.size()-1]==true){ b.flip()++;}
              }

           std::bitset<sizeof(B)*CHAR_BIT> d(a^b); // the difference between the two bit patterns
           for(i=0;i<=N-1;i++){ d.reset(i);}       // reset first N elements to zero
           if(d.any()==false){ return true;}       // if all other bits are zero then the two bit patterns are equal
           return false;
          }

template <class Type> bool Equality(std::complex<Type> A,std::complex<Type> B,int N)
         { return Equality(real(A),real(B),N) & Equality(imag(A),imag(B),N);}

#endif



