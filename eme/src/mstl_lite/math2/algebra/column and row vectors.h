
#include <cassert>
#include <algorithm>
#include <ostream>
#include <array>
#include <vector>

#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_CRVECTOR)
#define _CRVECTOR

// *******************************************************************

template <typename Type,typename CVType> class CVECTOREXPRESSION;

template <typename Type,std::size_t n> class CVECTOR;

typedef CVECTOR<int,0> ICVECTOR;
typedef CVECTOR<double,0> DCVECTOR;

//*************************************************************************

template <typename Type,typename CVType> class CVECTOREXPRESSION_NEGATE;

template <typename Type,typename RVType> class CVECTOREXPRESSION_RVECTORTRANSPOSE;

template <typename Type,typename RVType> class CVECTOREXPRESSION_RVECTORADJOINT;

template <typename Type,typename CVType> class CVECTOREXPRESSION_CONJUGATE;

template <typename Type,typename CVType1,typename CVType2> class CVECTOREXPRESSION_CVECTORADDITION;

template <typename Type,typename CVType1,typename CVType2> class CVECTOREXPRESSION_CVECTORSUBTRACTION;

template <typename Type,typename CVType> class CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION;

template <typename Type,typename CVType> class CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION;

template <typename Type,typename CVType> class CVECTOREXPRESSION_CVECTORSCALARDIVISION;

template <typename Type,typename MType,typename CVType> class CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION;

//*************************************************************************
//*************************************************************************
//************************* RVECTOR ***************************************
//*************************************************************************
//*************************************************************************

template <typename Type,typename RVType> class RVECTOREXPRESSION;

template <typename Type,std::size_t n> class RVECTOR;

typedef RVECTOR<int,0> IRVECTOR;
typedef RVECTOR<double,0> DRVECTOR;

//*************************************************************************

template <typename Type,typename RVType> class RVECTOREXPRESSION_NEGATE;

template <typename Type,typename CVType> class RVECTOREXPRESSION_CVECTORTRANSPOSE;

template <typename Type,typename CVType> class RVECTOREXPRESSION_CVECTORADJOINT;

template <typename Type,typename RVType> class RVECTOREXPRESSION_CONJUGATE;

template <typename Type,typename RVType1,typename RVType2> class RVECTOREXPRESSION_RVECTORADDITION;

template <typename Type,typename RVType1,typename RVType2> class RVECTOREXPRESSION_RVECTORSUBTRACTION;

template <typename Type,typename RVType> class RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION;

template <typename Type,typename RVType> class RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION;

template <typename Type,typename RVType> class RVECTOREXPRESSION_RVECTORSCALARDIVISION;

template <typename Type,typename RVType,typename MType> class RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION;

// *******************************************************************

// forward declarations
template <typename Type,typename MType> class MATRIXEXPRESSION;
template <typename Type,std::size_t n1,std::size_t n2> class MATRIX;

//*********************************************************************
//***************************** CVECTOR *******************************
//*********************************************************************

template <typename Type,typename CVType>
class CVECTOREXPRESSION 
	 { public : virtual ~CVECTOREXPRESSION(void) {;}

                    Type operator[](int i) const { return reinterpret_cast<CVType const&>(*this)[i];}

                    std::size_t Size(void) const { return reinterpret_cast<CVType const&>(*this).Size();}
                    std::size_t N(void) const { return Size();}

                    operator CVType const&() const { return reinterpret_cast<CVType const&>(*this); }
          };

//*************************************************************************

template <typename Type,std::size_t n=0> class CVECTOR : public CVECTOREXPRESSION<Type,CVECTOR<Type,n> >
         { private : static std::string CLASS_NAME;

           protected : std::array<Type,n> v;	   

	   public :

           CVECTOR(void){ v.fill(Zero<Type>());}
   	   explicit CVECTOR(Type,...);
           explicit CVECTOR(std::vector<Type> const &V);
           explicit CVECTOR(CVECTOR<Type,n> const &CV) : v(CV.v) {;}
           explicit CVECTOR(CVECTOR<Type,0> const &CV);
           template <typename CVType> CVECTOR(CVECTOREXPRESSION<Type,CVType> const&);

           ~CVECTOR(void) {;}

           bool Empty(void) const { return false;}

           CVECTOR<Type,n>& operator=(std::vector<Type> const&);
           CVECTOR<Type,n>& operator=(CVECTOR<Type,n> const&);
           CVECTOR<Type,n>& operator=(CVECTOR<Type,0> const&);
           template <typename CVType> CVECTOR<Type,n>& operator=(CVECTOREXPRESSION<Type,CVType> const&);

           Type operator[](int i) const { return v[i];}
           Type& operator[](int i){ return v[i];}

           std::size_t Size(void) const { return n;}
           std::size_t N(void) const { return Size();}

           CVECTOR<Type,n>& operator*=(Type const&);
           CVECTOR<Type,n>& operator/=(Type const&);

           template <typename CVType> CVECTOR<Type,n>& operator+=(CVECTOREXPRESSION<Type,CVType> const&);
           template <typename CVType> CVECTOR<Type,n>& operator-=(CVECTOREXPRESSION<Type,CVType> const&);
          };

template <typename Type> class CVECTOR<Type,0> : public CVECTOREXPRESSION<Type,CVECTOR<Type,0> >
         { private : static std::string CLASS_NAME;

           protected : std::vector<Type> v;	   

	   public :

           CVECTOR(void) : v() {;}
	   explicit CVECTOR(const std::size_t N) : v(N,Zero<Type>()) {;}
   	   explicit CVECTOR(const std::size_t,Type,...);
           explicit CVECTOR(std::vector<Type> const &V) : v(V) {;}
           template <std::size_t m> explicit CVECTOR(CVECTOR<Type,m> const &CV) : v(CV.v) {;}
           explicit CVECTOR(CVECTOR<Type,0> const &CV) : v(CV.v) {;}
           template <typename CVType> CVECTOR(CVECTOREXPRESSION<Type,CVType> const&);

           ~CVECTOR(void) {;}

           bool Empty(void) const { return v.empty();}

           CVECTOR<Type,0>& operator=(std::vector<Type> const&);
           CVECTOR<Type,0>& operator=(CVECTOR<Type,0> const&);
           template <std::size_t m> CVECTOR<Type,0>& operator=(CVECTOR<Type,m> const&);
           template <typename CVType> CVECTOR<Type,0>& operator=(CVECTOREXPRESSION<Type,CVType> const&);
   	   CVECTOR<Type,0>& operator=(Type const&T);

           Type operator[](int i) const { return v[i];}
           Type& operator[](int i){ return v[i];}

           std::size_t Size(void) const { return v.size();}
           std::size_t N(void) const { return Size();}

           CVECTOR<Type,0>& operator*=(Type const&);
           CVECTOR<Type,0>& operator/=(Type const&);

           template <typename CVType> CVECTOR<Type,0>& operator+=(CVECTOREXPRESSION<Type,CVType> const&);
           template <typename CVType> CVECTOR<Type,0>& operator-=(CVECTOREXPRESSION<Type,CVType> const&);
          };

//*************************************************************************

template <typename Type,typename CVType> class CVECTOREXPRESSION_NEGATE : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_NEGATE<Type,CVType> >
      { private : CVType const&cv;
        public : explicit CVECTOREXPRESSION_NEGATE(CVType const&CV) : cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return -cv[i];} 
       };

template <typename Type,typename RVType> class CVECTOREXPRESSION_RVECTORTRANSPOSE : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_RVECTORTRANSPOSE<Type,RVType> >
      { private : RVType const&rv;
        public : explicit CVECTOREXPRESSION_RVECTORTRANSPOSE(RVType const&RV) : rv(RV) {;}
                 std::size_t Size() const { return rv.Size();}
                 Type operator[](int i) const { return rv[i];} 
       };

template <typename Type,typename RVType> class CVECTOREXPRESSION_RVECTORADJOINT : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_RVECTORADJOINT<Type,RVType> >
      { private : RVType const&rv;
        public : explicit CVECTOREXPRESSION_RVECTORADJOINT(RVType const&RV) : rv(RV) {;}
                 std::size_t Size() const { return rv.sSize();}
                 Type operator[](int i) const { return conj(rv[i]);} 
       };

template <typename Type,typename CVType> class CVECTOREXPRESSION_CONJUGATE : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_CONJUGATE<Type,CVType> >
      { private : CVType const&cv;
        public : explicit CVECTOREXPRESSION_CONJUGATE(CVType const&CV) : cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return conj(cv[i]);} 
       };

template <typename Type,typename CVType1,typename CVType2> class CVECTOREXPRESSION_CVECTORADDITION : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_CVECTORADDITION<Type,CVType1,CVType2> >
      { private : CVType1 const&cv1; CVType2 const&cv2;
        public : explicit CVECTOREXPRESSION_CVECTORADDITION(CVType1 const&CV1,CVType2 const&CV2) : cv1(CV1), cv2(CV2) { assert(cv1.Size()==cv2.Size());}
                 std::size_t Size() const { return cv1.Size();}
                 Type operator[](int i) const { return cv1[i]+cv2[i];} 
       };

template <typename Type,typename CVType1,typename CVType2> class CVECTOREXPRESSION_CVECTORSUBTRACTION : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_CVECTORSUBTRACTION<Type,CVType1,CVType2> >
      { private : CVType1 const&cv1; CVType2 const&cv2;
        public : explicit CVECTOREXPRESSION_CVECTORSUBTRACTION(CVType1 const&CV1,CVType2 const&CV2) : cv1(CV1), cv2(CV2) { assert(cv1.Size()==cv2.Size());}
                 std::size_t Size() const { return cv1.Size();}
                 Type operator[](int i) const { return cv1[i]-cv2[i];} 
       };

template <typename Type,typename CVType> class CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION<Type,CVType> >
      { private : Type const&s; CVType const&cv;
        public : explicit CVECTOREXPRESSION_CVECTORSCALARMULTIPLICATION(CVType const&CV,Type const&S) : s(S), cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return cv[i]*s;} 
       };

template <typename Type,typename CVType> class CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION<Type,CVType> >
      { private : Type const&s; CVType const&cv;
        public : CVECTOREXPRESSION_SCALARCVECTORMULTIPLICATION(Type const&S,CVType const&CV) : s(S), cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return s*cv[i];} 
       };

//*************************************************************************

template <typename Type,typename CVType> class CVECTOREXPRESSION_CVECTORSCALARDIVISION : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_CVECTORSCALARDIVISION<Type,CVType> >
      { private : Type const&s; CVType const&cv;
        public : CVECTOREXPRESSION_CVECTORSCALARDIVISION(CVType const&CV,Type const&S) : s(S), cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return cv[i]/s;} 
       };

//*************************************************************************

template <typename Type,typename MType,typename CVType> class CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION : public CVECTOREXPRESSION<Type,CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION<Type,MType,CVType> >
      { private : MType const&m; CVType const&cv;
        public : explicit CVECTOREXPRESSION_MATRIXCVECTORMULTIPLICATION(MType const&M,CVType const&CV) : m(M), cv(CV) { assert(m.N2()==cv.N());}
                 std::size_t Size() const { return m.N1();}
                 Type operator[](int i) const; 
       };

//*************************************************************************
//*************************************************************************
//************************* RVECTOR ***************************************
//*************************************************************************
//*************************************************************************

template <typename Type,typename RVType>
class RVECTOREXPRESSION 
	 { public : virtual ~RVECTOREXPRESSION(void) {;}

                    Type operator[](int i) const { return reinterpret_cast<RVType const&>(*this)[i];}

                    std::size_t N(void) const { return Size();}
                    std::size_t Size(void) const { return reinterpret_cast<RVType const&>(*this).Size();}

                    operator RVType const&() const { return reinterpret_cast<RVType const&>(*this); }
          };

//*************************************************************************

template <typename Type,std::size_t n=0> class RVECTOR : public RVECTOREXPRESSION<Type,RVECTOR<Type,n> >
         { private : static std::string CLASS_NAME;

           protected : std::array<Type,n> v;	   

	   public :

           RVECTOR(void){ v.fill(Zero<Type>());}
   	   explicit RVECTOR(Type,...);
           explicit RVECTOR(std::vector<Type> const &V);
           explicit RVECTOR(RVECTOR<Type,n> const &RV) : v(RV.v) {;}
           explicit RVECTOR(RVECTOR<Type,0> const &RV);
           template <typename RVType> RVECTOR(RVECTOREXPRESSION<Type,RVType> const&);

           ~RVECTOR(void) {;}

           bool Empty(void) const { return false;}

           RVECTOR<Type,n>& operator=(std::vector<Type> const&);
           RVECTOR<Type,n>& operator=(RVECTOR<Type,n> const&);
           RVECTOR<Type,n>& operator=(RVECTOR<Type,0> const&);
           template <typename RVType> RVECTOR<Type,n>& operator=(RVECTOREXPRESSION<Type,RVType> const&);

           Type operator[](int i) const { return v[i];}
           Type& operator[](int i){ return v[i];}

           std::size_t Size(void) const { return v.size();}
           std::size_t N(void) const { return Size();}

           RVECTOR<Type,n>& operator*=(Type const&);
           RVECTOR<Type,n>& operator/=(Type const&);

           template <typename MType> RVECTOR<Type,n>& operator*=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> RVECTOR<Type,n>& operator/=(MATRIXEXPRESSION<Type,MType> const&);

           template <typename RVType> RVECTOR<Type,n>& operator+=(RVECTOREXPRESSION<Type,RVType> const&);
           template <typename RVType> RVECTOR<Type,n>& operator-=(RVECTOREXPRESSION<Type,RVType> const&);
          };

template <typename Type> class RVECTOR<Type,0> : public RVECTOREXPRESSION<Type,RVECTOR<Type,0> >
         { private : static std::string CLASS_NAME;

           protected : std::vector<Type> v;	   

	   public :

           RVECTOR(void) : v() {;}
	   explicit RVECTOR(const std::size_t N) : v(N,Zero<Type>()) {;}
   	   explicit RVECTOR(const std::size_t,Type,...);
           explicit RVECTOR(std::vector<Type> const &V) : v(V) {;}
           template<std::size_t m> explicit RVECTOR(RVECTOR<Type,m> const &RV) : v(RV.v) {;}
           template <typename RVType> RVECTOR(RVECTOREXPRESSION<Type,RVType> const&);

           ~RVECTOR(void) {;}

           bool Empty(void) const { return v.empty();}

           RVECTOR<Type,0>& operator=(std::vector<Type> const&);
           RVECTOR<Type,0>& operator=(RVECTOR<Type,0> const&);
           template<std::size_t m> RVECTOR<Type,0>& operator=(RVECTOR<Type,m> const&);
           template <typename RVType> RVECTOR<Type,0>& operator=(RVECTOREXPRESSION<Type,RVType> const&);
   	   RVECTOR<Type,0>& operator=(Type const&T);

           Type operator[](int i) const { return v[i];}
           Type& operator[](int i){ return v[i];}

           std::size_t Size(void) const { return v.size();}
           std::size_t N(void) const { return Size();}

           RVECTOR<Type,0>& operator*=(Type const&);
           RVECTOR<Type,0>& operator/=(Type const&);

           template <typename MType> RVECTOR<Type,0>& operator*=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> RVECTOR<Type,0>& operator/=(MATRIXEXPRESSION<Type,MType> const&);

           template <typename RVType> RVECTOR<Type,0>& operator+=(RVECTOREXPRESSION<Type,RVType> const&);
           template <typename RVType> RVECTOR<Type,0>& operator-=(RVECTOREXPRESSION<Type,RVType> const&);
          };

//*************************************************************************

template <typename Type,typename RVType> class RVECTOREXPRESSION_NEGATE : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_NEGATE<Type,RVType> >
      { private : RVType const&rv;
        public : explicit RVECTOREXPRESSION_NEGATE(RVType const&RV) : rv(RV) {;}
                 std::size_t Size() const { return rv.Size();}
                 Type operator[](int i) const { return -rv[i];} 
       };

template <typename Type,typename CVType> class RVECTOREXPRESSION_CVECTORTRANSPOSE : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_CVECTORTRANSPOSE<Type,CVType> >
      { private : CVType const&cv;
        public : explicit RVECTOREXPRESSION_CVECTORTRANSPOSE(CVType const&CV) : cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return cv[i];} 
       };

template <typename Type,typename CVType> class RVECTOREXPRESSION_CVECTORADJOINT : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_CVECTORADJOINT<Type,CVType> >
      { private : CVType const&cv;
        public : explicit RVECTOREXPRESSION_CVECTORADJOINT(CVType const&CV) : cv(CV) {;}
                 std::size_t Size() const { return cv.Size();}
                 Type operator[](int i) const { return conj(cv[i]);} 
       };

template <typename Type,typename RVType> class RVECTOREXPRESSION_CONJUGATE : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_CONJUGATE<Type,RVType> >
      { private : RVType const&rv;
        public : explicit RVECTOREXPRESSION_CONJUGATE(RVType const&RV) : rv(RV) {;}
                 std::size_t Size() const { return rv.Size();}
                 Type operator[](int i) const { return conj(rv[i]);} 
       };

template <typename Type,typename RVType1,typename RVType2> class RVECTOREXPRESSION_RVECTORADDITION : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_RVECTORADDITION<Type,RVType1,RVType2> >
      { private : RVType1 const&rv1; RVType2 const&rv2;
        public : explicit RVECTOREXPRESSION_RVECTORADDITION(RVType1 const&RV1,RVType2 const&RV2) : rv1(RV1), rv2(RV2) { assert(rv1.Size()==rv2.Size());}
                 std::size_t Size() const { return rv1.Size();}
                 Type operator[](int i) const { return rv1[i]+rv2[i];} 
       };

template <typename Type,typename RVType1,typename RVType2> class RVECTOREXPRESSION_RVECTORSUBTRACTION : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_RVECTORSUBTRACTION<Type,RVType1,RVType2> >
      { private : RVType1 const&rv1; RVType2 const&rv2;
        public : explicit RVECTOREXPRESSION_RVECTORSUBTRACTION(RVType1 const&RV1,RVType2 const&RV2) : rv1(RV1), rv2(RV2) { assert(rv1.Size()==rv2.Size());}
                 std::size_t Size() const { return rv1.Size();}
                 Type operator[](int i) const { return rv1[i]-rv2[i];} 
       };

template <typename Type,typename RVType> class RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION<Type,RVType> >
      { private : Type const&s; RVType const&rv;
        public : explicit RVECTOREXPRESSION_RVECTORSCALARMULTIPLICATION(RVType const&RV,Type const&S) : s(S), rv(RV) {;}
                 std::size_t Size() const { return rv.Size();}
                 Type operator[](int i) const { return rv[i]*s;} 
       };

template <typename Type,typename RVType> class RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION<Type,RVType> >
      { private : Type const&s; RVType const&rv;
        public : explicit RVECTOREXPRESSION_SCALARRVECTORMULTIPLICATION(Type const&S,RVType const&RV) : s(S), rv(RV) {;}
                 std::size_t Size() const { return rv.Size();}
                 Type operator[](int i) const { return s*rv[i];} 
       };

template <typename Type,typename RVType> class RVECTOREXPRESSION_RVECTORSCALARDIVISION : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_RVECTORSCALARDIVISION<Type,RVType> >
      { private : Type const&s; RVType const&rv;
        public : explicit RVECTOREXPRESSION_RVECTORSCALARDIVISION(RVType const&RV,Type const&S) : s(S), rv(RV) {;}
                 std::size_t Size() const { return rv.Size();}
                 Type operator[](int i) const { return rv[i]/s;} 
       };

template <typename Type,typename RVType,typename MType> class RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION : public RVECTOREXPRESSION<Type,RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION<Type,RVType,MType> >
      { private : RVType const&rv; MType const&m;
        public : explicit RVECTOREXPRESSION_RVECTORMATRIXMULTIPLICATION(RVType const&RV,MType const&M) : rv(RV), m(M) { assert(rv.N()==m.N1());}
                 std::size_t Size() const { return m.N2();}
                 Type operator[](int i) const; 
       };

#endif
