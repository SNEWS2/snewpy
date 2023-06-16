
#include <cstddef>
#include <bitset>

#include "mstl.h"

#if !defined(_BITSET2)
#define _BITSET2

template<std::size_t N> std::bitset<N> operator++(std::bitset<N> &,int);

template<std::size_t N> std::bitset<N> operator++(std::bitset<N> &bs,int)
        { std::bitset<N> bs0(bs); 
          bool carry=1;
          std::size_t n=0;
          do{ carry=(bs[n]|carry)^(bs[n]^carry); bs[n]=!carry; n++;} while(carry==1 && n<=N);
          return bs0;
         }

#endif
