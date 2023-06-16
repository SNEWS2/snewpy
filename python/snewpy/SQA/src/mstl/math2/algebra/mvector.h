
#include <algorithm>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>

#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_MVECTOR)
#define _MVECTOR

// *******************************************************************

template <typename Type,typename VType> class MVECTOREXPRESSION;
template <typename Type,std::size_t N> class MVECTOR;

template <typename Type,typename VType> class MVECTOREXPRESSION_NEGATE;

template <typename Type,typename VType1,typename VType2> class MVECTOREXPRESSION_VECTORADDITION;
template <typename Type,typename VType1,typename VType2> class MVECTOREXPRESSION_VECTORSUBTRACTION;

template <typename Type,typename VType,typename SType> class MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION;
template <typename Type,typename SType,typename VType> class MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION;
template <typename Type,typename VType,typename SType> class MVECTOREXPRESSION_VECTORSCALARDIVISION;

template <typename Type,typename VType1,typename VType2> class MVECTOREXPRESSION_CROSSPRODUCT;

// *******************************************************************
// *******************************************************************
// *******************************************************************

// operator prototypes

template <typename Type,typename VType> 
MVECTOREXPRESSION_NEGATE<Type,VType> operator-(MVECTOREXPRESSION<Type,VType> const&);

// ************************
 
template <typename Type,typename VType1,typename VType2>
MVECTOREXPRESSION_VECTORADDITION<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> > operator+(MVECTOREXPRESSION<Type,VType1> const&,MVECTOREXPRESSION<Type,VType2> const&);

template <typename Type,typename VType1,typename VType2>
MVECTOREXPRESSION_VECTORSUBTRACTION<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> > operator-(MVECTOREXPRESSION<Type,VType1> const&,MVECTOREXPRESSION<Type,VType2> const&);

template <typename Type,typename SType,typename VType>
MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION<Type,MVECTOREXPRESSION<Type,VType>,SType> operator*(MVECTOREXPRESSION<Type,VType> const&,SType const&);

template <typename Type,typename SType,typename VType>
MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION<Type,SType,MVECTOREXPRESSION<Type,VType> > operator*(SType const&,MVECTOREXPRESSION<Type,VType> const&);

template <typename Type,typename SType,typename VType>
MVECTOREXPRESSION_VECTORSCALARDIVISION<Type,MVECTOREXPRESSION<Type,VType>,SType> operator/(MVECTOREXPRESSION<Type,VType> const&,SType const&);

// cross product
template <typename Type,typename VType1,typename VType2>
MVECTOREXPRESSION_CROSSPRODUCT<Type,MVECTOREXPRESSION<Type,VType1>,MVECTOREXPRESSION<Type,VType2> > operator^(MVECTOREXPRESSION<Type,VType1> const&,MVECTOREXPRESSION<Type,VType2> const&);

// ******************************************

//dot product
template <typename Type,typename VType1,typename VType2> 
Type operator*(MVECTOREXPRESSION<Type,VType1> const&,MVECTOREXPRESSION<Type,VType2> const&);

template <typename Type,typename VType>
std::ostream& operator<<(std::ostream&,MVECTOREXPRESSION<Type,VType> const&);

// ***************************************************************
// ***************************************************************
// ***************************************************************

template <typename Type,typename VType>
class MVECTOREXPRESSION
      { public : virtual ~MVECTOREXPRESSION(void) {;}

                 Type operator[](int i) const { return reinterpret_cast<VType const&>(*this)[i];}

                 std::size_t Size(void) const { return reinterpret_cast<VType const&>(*this).size();}
                 std::size_t N(void) const { return Size();}

                 operator VType const&() const { return reinterpret_cast<VType const&>(*this);}
       };

// *******************************************************************

template <typename Type,std::size_t n=0> 
class MVECTOR : public MVECTOREXPRESSION<Type,MVECTOR<Type,n> >
	 { private : static std::string CLASS_NAME;

           protected : std::array<Type,n> v;

	   public : MVECTOR(void){ v.fill(Zero<Type>());}
                    explicit MVECTOR(Type,...);
	            explicit MVECTOR(std::vector<Type> const &V);
                    explicit MVECTOR(MVECTOR<Type,n> const &V) : v(V.v) {;}
                    explicit MVECTOR(MVECTOR<Type,0> const &V);
                    template <typename VType> MVECTOR(MVECTOREXPRESSION<Type,VType> const&);

                    bool Empty(void) const { return false;}

                    MVECTOR<Type,n>& operator=(std::vector<Type> const&);
                    MVECTOR<Type,n>& operator=(MVECTOR<Type,n> const&);
                    MVECTOR<Type,n>& operator=(MVECTOR<Type,0> const&);
                    template <typename VType> MVECTOR<Type,n>& operator=(MVECTOREXPRESSION<Type,VType> const&);

                    template <typename VType> MVECTOR<Type,n>& operator+=(MVECTOREXPRESSION<Type,VType> const&);
                    template <typename VType> MVECTOR<Type,n>& operator-=(MVECTOREXPRESSION<Type,VType> const&);

                    MVECTOR<Type,n>& operator*=(Type const&T);
                    MVECTOR<Type,n>& operator/=(Type const&T);

                    operator std::vector<Type>() { return v;}

                    Type operator[](int i) const { return v[i];}
                    Type& operator[](int i){ return v[i];}

                    std::size_t Size(void) const { return v.size();}
                    std::size_t N(void) const { return Size();}
	  };

template <typename Type> 
class MVECTOR<Type,0> : public MVECTOREXPRESSION<Type,MVECTOR<Type,0> >
	 { private : static std::string CLASS_NAME;

           protected : std::vector<Type> v;

	   public : MVECTOR(void) : v() {;}
                    explicit MVECTOR(std::size_t N) : v(N,Zero<Type>()) {;}
                    explicit MVECTOR(std::size_t,Type,...);
	            explicit MVECTOR(std::vector<Type> const &V) : v(V) {;} 
                    explicit MVECTOR(MVECTOR<Type,0> const &V) : v(V.v) {;}
                    template <std::size_t n> explicit MVECTOR(MVECTOR<Type,n> const &V) : v(V.v) {;}  
                    template <typename VType> MVECTOR(MVECTOREXPRESSION<Type,VType> const&);

                    bool Empty(void) const { return v.empty();}

                    MVECTOR<Type,0>& operator=(std::vector<Type> const&);
                    MVECTOR<Type,0>& operator=(MVECTOR<Type,0> const&);
                    template <std::size_t n> MVECTOR<Type,0>& operator=(MVECTOR<Type,n> const&);
                    template <typename VType> MVECTOR<Type>& operator=(MVECTOREXPRESSION<Type,VType> const&);

                    template <typename VType> MVECTOR<Type,0>& operator+=(MVECTOREXPRESSION<Type,VType> const&);
                    template <typename VType> MVECTOR<Type,0>& operator-=(MVECTOREXPRESSION<Type,VType> const&);

                    MVECTOR<Type,0>& operator*=(Type const&T);
                    MVECTOR<Type,0>& operator/=(Type const&T);

                    operator std::vector<Type>() { return v;}

                    Type operator[](int i) const { return v[i];}
                    Type& operator[](int i){ return v[i];}

                    std::size_t Size(void) const { return v.size();}
                    std::size_t N(void) const { return Size();}
	  };

typedef MVECTOR<int,0> IVECTOR;
typedef MVECTOR<double,0> DVECTOR;

// *******************************************************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename VType> 
class MVECTOREXPRESSION_NEGATE : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_NEGATE<Type,VType> >
      { private : VType const&v;
        public : explicit MVECTOREXPRESSION_NEGATE(VType const&V) : v(V) {;}
                 std::size_t Size() const { return v.size();}
                 Type operator[](int i) const { return -v[i];} 
       };

template <typename Type,typename VType1,typename VType2> 
class MVECTOREXPRESSION_VECTORADDITION : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_VECTORADDITION<Type,VType1,VType2> >
      { private : VType1 const&v1; VType2 const&v2;
        public : explicit MVECTOREXPRESSION_VECTORADDITION(VType1 const&V1,VType2 const&V2) : v1(V1), v2(V2) { assert(v1.size()==v2.size());}
                 std::size_t Size() const { return v1.size();}
                 Type operator[](int i) const { return v1[i] + v2[i];} 
       };

template <typename Type,typename VType1,typename VType2> 
class MVECTOREXPRESSION_VECTORSUBTRACTION : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_VECTORSUBTRACTION<Type,VType1,VType2> >
      { private : VType1 const&v1; VType2 const&v2;
        public : explicit MVECTOREXPRESSION_VECTORSUBTRACTION(VType1 const&V1,VType2 const&V2) : v1(V1), v2(V2) { assert(v1.size()==v2.size());}
                 std::size_t Size() const { return v1.size(); }
                 Type operator[](int i) const { return v1[i] - v2[i];} 
       };

template <typename Type,typename VType,typename SType> 
class MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION<Type,VType,SType> >
      { private : VType const&v; SType const&s;
        public : explicit MVECTOREXPRESSION_VECTORSCALARMULTIPLICATION(VType const&V,SType const&S) : v(V), s(S) {;}
                 std::size_t Size() const { return v.size(); }
                 Type operator[](int i) const { return v[i]*s;} 
       };

template <typename Type,typename SType,typename VType> 
class MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION<Type,SType,VType> >
      { private : VType const&v; SType const&s;
        public : explicit MVECTOREXPRESSION_SCALARVECTORMULTIPLICATION(SType const&S,VType const&V) : v(V), s(S) {;}
                 std::size_t Size() const { return v.size(); }
                 Type operator[](int i) const { return s*v[i];} 
       };

template <typename Type,typename VType,typename SType> 
class MVECTOREXPRESSION_VECTORSCALARDIVISION : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_VECTORSCALARDIVISION<Type,VType,SType> >
      { private : VType const&v; SType const&s;
        public : explicit MVECTOREXPRESSION_VECTORSCALARDIVISION(VType const&V,SType const&S) : v(V), s(S) {;}
                 std::size_t Size() const { return v.size(); }
                 Type operator[](int i) const { return v[i]/s;} 
       };

template <typename Type,typename VType1,typename VType2> 
class MVECTOREXPRESSION_CROSSPRODUCT : public MVECTOREXPRESSION<Type,MVECTOREXPRESSION_CROSSPRODUCT<Type,VType1,VType2> >
      { private : VType1 const&v1; VType2 const&v2;
        public : explicit MVECTOREXPRESSION_CROSSPRODUCT(VType1 const&V1,VType2 const&V2) : v1(V1), v2(V2) { assert(v1.size()==3); assert(v2.size()==3);}
                 std::size_t Size() const { return v1.size(); }
                 Type operator[](int i) const { return v1[(i+1)%3]*v2[(i+2)%3] - v1[(i+2)%3]*v2[(i+1)%3];}
       };

// *******************************************************************
// *******************************************************************
// *******************************************************************

#endif

