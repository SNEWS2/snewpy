
#include<cassert>
#include<cstddef>
#include <iostream>
#include <algorithm> 
#include <array>
#include <vector>

#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_MMATRIX)
#define _MMATRIX

// *******************************************************************
// *******************************************************************
// ******************************* MATRIX ****************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename MType> class MATRIXEXPRESSION;

// *******************************************************************

template <typename Type,std::size_t n1,std::size_t n2> class MATRIX;

typedef MATRIX<int,0,0> IMATRIX;
typedef MATRIX<double,0,0> DMATRIX;

// *******************************************************************

template <typename Type,typename MType> class MATRIXEXPRESSION_NEGATE;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXTRANSPOSE;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXANTITRANSPOSE;

// ********************
 
template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXADJOINT;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_CONJUGATE; 

// ********************

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXADDITION; 

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXSUBTRACTION; 

// ********************

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXMULTIPLICATION;

//template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXDIVISION;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARADDITION;

template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXADDITION;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION;

template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION;

template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION;

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARDIVISION;

//template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXDIVISION;

// ********************

template <typename Type,typename CVType,typename RVType> class MATRIXEXPRESSION_VECTORPRODUCT;

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_SUBMATRIXADDITION;

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_SUBMATRIXSUBTRACTION;

// *******************************************************************
// *******************************************************************
// *******************************************************************

// forward declarations
template <typename Type,typename CVType> class CVECTOREXPRESSION;
template <typename Type,std::size_t n> class CVECTOR;

template <typename Type,typename RVType> class RVECTOREXPRESSION;
template <typename Type,std::size_t n> class RVECTOR;

// *******************************************************************
// *******************************************************************
// ******************************* MATRIX ****************************
// *******************************************************************
// *******************************************************************

template <typename Type,typename MType>
class MATRIXEXPRESSION
      { public : virtual ~MATRIXEXPRESSION(void) {;}

                 Type operator()(int i,int j) const { return reinterpret_cast<MType const&>(*this)(i,j);}

                 std::size_t N1(void) const { return reinterpret_cast<MType const&>(*this).N1();}
                 std::size_t N2(void) const { return reinterpret_cast<MType const&>(*this).N2();}
 
                 operator MType const&() const { return reinterpret_cast<MType const&>(*this);}
       };

// *******************************************************************

template <typename Type,std::size_t n1=0,std::size_t n2=n1> class MATRIX : public MATRIXEXPRESSION<Type,MATRIX<Type,n1,n2> >
         { private :
           static std::string CLASS_NAME;

           protected :
           std::array<std::array<Type,n2>,n1> m;

           void Copy(std::vector<std::vector<Type> > const&);
           template <typename MType> void Copy(MATRIXEXPRESSION<Type,MType> const&);

           void Initialize(void);

           public :
           MATRIX(void){ Initialize();}
           explicit MATRIX(Type,...);
           explicit MATRIX(int N1,int N2){ Initialize();}
           explicit MATRIX(std::vector<std::vector<Type> > const&);
           explicit MATRIX(MATRIX<Type,n1,n2> const&);
           explicit MATRIX(MATRIX<Type,0,0> const&);
           template <typename MType> MATRIX(MATRIXEXPRESSION<Type,MType> const&);

           ~MATRIX(void) {;}

           bool Empty(void) const { return false;}
           std::size_t N1(void) const { return n1;} // Number of rows
           std::size_t N2(void) const { return n2;} // Number of cols
           std::size_t N(void);

           MATRIX<Type,n1,n2>& Conjugate(void);

           CVECTOR<Type,n1> Column(int i) const;
           RVECTOR<Type,n2> Row(int i) const;
           MATRIX<Type,n1,n2>& SwapRows(int i,int j);
           MATRIX<Type,n1,n2>& SwapColumns(int i,int j);

           template <typename MType> void SetSubMatrix(MATRIXEXPRESSION<Type,MType> const&,int,int);
           template <typename RVType> void SetSubRow(RVECTOREXPRESSION<Type,RVType> const&,int,int);
           template <typename CVType> void SetSubColumn(CVECTOREXPRESSION<Type,CVType> const&,int,int);

           std::array<Type,n2>& operator[](int i);       // Set or Return an element
           std::array<Type,n2> operator[](int i) const;  // Return an element

           Type& operator()(int i,int j);       // Set or Return an element
           Type operator()(int i,int j) const;  // Return an element

           MATRIX<Type,n1,n2>& operator=(MATRIX<Type,n1,n2> const &);
           MATRIX<Type,n1,n2>& operator=(MATRIX<Type,0,0> const &);

           template <typename MType> MATRIX<Type,n1,n2>& operator=(MATRIXEXPRESSION<Type,MType> const &);
           template <typename MType> MATRIX<Type,n1,n2>& operator+=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,n1,n2>& operator-=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,n1,n2>& operator*=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,n1,n2>& operator/=(MATRIXEXPRESSION<Type,MType> const&);

           MATRIX<Type,n1,n2>& operator*=(Type const&);
           MATRIX<Type,n1,n2>& operator/=(Type const&);

           template <typename MType> MATRIX<Type,n1,n2>& SubMatrixAddition(MATRIXEXPRESSION<Type,MType> const&,int,int);
           template <typename MType> MATRIX<Type,n1,n2>& SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType> const&,int,int);

           template <typename MType> bool operator==(MATRIXEXPRESSION<Type,MType> const&) const;
           template <typename MType> bool operator!=(MATRIXEXPRESSION<Type,MType> const&) const;
       	 };

template <typename Type,std::size_t n> class MATRIX<Type,n,n> : public MATRIXEXPRESSION<Type,MATRIX<Type,n,n> >
         { private :
           static std::string CLASS_NAME;

           protected :
           std::array<std::array<Type,n>,n> m;

           void Copy(std::vector<std::vector<Type> > const&);
           template <typename MType> void Copy(MATRIXEXPRESSION<Type,MType> const&);

           void Initialize(void);

           public :
           MATRIX(void){ Initialize();}
           explicit MATRIX(Type,...);
           explicit MATRIX(int N){ Initialize();}
           explicit MATRIX(std::vector<std::vector<Type> > const&);
           explicit MATRIX(MATRIX<Type,n,n> const&);
           explicit MATRIX(MATRIX<Type,0,0> const&);
           template <typename MType> MATRIX(MATRIXEXPRESSION<Type,MType> const&);

           ~MATRIX(void) {;}

           bool Empty(void) const { return false;}
           std::size_t N1(void) const { return n;} // Number of rows
           std::size_t N2(void) const { return n;} // Number of cols
           std::size_t N(void) const { return n;};

           MATRIX<Type,n,n>& Conjugate(void);
           MATRIX<Type,n,n>& Transpose(void);
           MATRIX<Type,n,n>& AntiTranspose(void);
           MATRIX<Type,n,n>& Adjoint(void);

           MATRIX<Type,n,n>& Inverse(void) { return Invert();}
           MATRIX<Type,n,n>& Invert(void);
           MATRIX<Type,n,n>& DInvert(void);
           MATRIX<Type,n,n>& GJInvert(void);
           MATRIX<Type,n,n>& LInvert(bool message=true); // the bool is to whether the warning message is printed
           MATRIX<Type,n,n>& LTInvert(void);
           MATRIX<Type,n,n>& UTInvert(void);
           MATRIX<Type,n,n>& LUInvert(void);
           MATRIX<Type,n,n>& MPInvert(void);
           MATRIX<Type,n,n>& QRInvert(void);

           CVECTOR<Type,n> Column(int i) const;
           RVECTOR<Type,n> Row(int i) const;
           MATRIX<Type,n>& SwapRows(int i,int j);
           MATRIX<Type,n>& SwapColumns(int i,int j);

           template <typename MType> void SetSubMatrix(MATRIXEXPRESSION<Type,MType> const&,int,int);
           template <typename RVType> void SetSubRow(RVECTOREXPRESSION<Type,RVType> const&,int,int);
           template <typename CVType> void SetSubColumn(CVECTOREXPRESSION<Type,CVType> const&,int,int);

           std::array<Type,n>& operator[](int i);      // Set or Return an element
           std::array<Type,n> operator[](int i) const; // Return an element

           Type& operator()(int i,int j);      // Set or Return an element
           Type operator()(int i,int j) const; // Return an element

           MATRIX<Type,n,n>& operator=(MATRIX<Type,n,n> const &);
           MATRIX<Type,n,n>& operator=(MATRIX<Type,0,0> const &);

           template <typename MType> MATRIX<Type,n>& operator=(MATRIXEXPRESSION<Type,MType> const &);
           template <typename MType> MATRIX<Type,n>& operator+=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,n>& operator-=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,n>& operator*=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,n>& operator/=(MATRIXEXPRESSION<Type,MType> const&);

           MATRIX<Type,n,n>& operator=(Type const&);
           MATRIX<Type,n,n>& operator+=(Type const&);
           MATRIX<Type,n,n>& operator-=(Type const&);
           MATRIX<Type,n,n>& operator*=(Type const&);
           MATRIX<Type,n,n>& operator/=(Type const&);

           template <typename MType> MATRIX<Type,n,n>& SubMatrixAddition(MATRIXEXPRESSION<Type,MType> const&,int,int);
           template <typename MType> MATRIX<Type,n,n>& SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType> const&,int,int);

           template <typename MType> bool operator==(MATRIXEXPRESSION<Type,MType> const&) const;
           template <typename MType> bool operator!=(MATRIXEXPRESSION<Type,MType> const&) const;
       	 };

template <typename Type> class MATRIX<Type,0,0> : public MATRIXEXPRESSION<Type,MATRIX<Type,0,0> >
         { private :
           static std::string CLASS_NAME;

           protected :
           bool empty;
           std::size_t n1, n2;  // n1 the number of rows, n2 the number of cols
           std::vector<std::vector<Type> > m;

           void Create(const std::size_t N1,const std::size_t N2);

           void Copy(std::vector<std::vector<Type> > const&);
           template <typename MType> void Copy(MATRIXEXPRESSION<Type,MType> const&);

           void Destroy(void);

           void Initialize(void);

           public :
           MATRIX(void){ Create(0,0);}
           explicit MATRIX(const std::size_t N1,const std::size_t N2){ Create(N1,N2); if(Empty()==false){ Initialize();} }
           explicit MATRIX(const std::size_t N){ Create(N,N); if(Empty()==false){ Initialize();} } // a square matrix_tag

           explicit MATRIX(const std::size_t,const std::size_t,Type,...);
           explicit MATRIX(std::vector<std::vector<Type> > const&);
           template <std::size_t m1,std::size_t m2> explicit MATRIX(MATRIX<Type,m1,m2> const&);
           explicit MATRIX(MATRIX<Type,0,0> const&);
           template <typename MType> MATRIX(MATRIXEXPRESSION<Type,MType> const&);

           ~MATRIX(void) { Destroy();}

           bool Empty(void) const { return empty;}
           std::size_t N1(void) const { return n1;} // Number of rows
           std::size_t N2(void) const { return n2;} // Number of cols
           std::size_t N(void);

           MATRIX<Type,0,0>& Transpose(void);
           MATRIX<Type,0,0>& AntiTranspose(void);
           MATRIX<Type,0,0>& Adjoint(void);
           MATRIX<Type,0,0>& Conjugate(void);

           MATRIX<Type,0,0>& Inverse(void) { return Invert();}
           MATRIX<Type,0,0>& Invert(void);
           MATRIX<Type,0,0>& DInvert(void);
           MATRIX<Type,0,0>& GJInvert(void);
           MATRIX<Type,0,0>& LInvert(bool message=true); // the bool is to whether the warning message is printed
           MATRIX<Type,0,0>& LTInvert(void);
           MATRIX<Type,0,0>& UTInvert(void);
           MATRIX<Type,0,0>& LUInvert(void);
           MATRIX<Type,0,0>& MPInvert(void);
           MATRIX<Type,0,0>& QRInvert(void);

           CVECTOR<Type,0> Column(int i) const;
           RVECTOR<Type,0> Row(int i) const;
           MATRIX<Type,0,0>& SwapRows(int i,int j);
           MATRIX<Type,0,0>& SwapColumns(int i,int j);

           template <typename MType> void SetSubMatrix(MATRIXEXPRESSION<Type,MType> const&,int,int);
           template <typename RVType> void SetSubRow(RVECTOREXPRESSION<Type,RVType> const&,int,int);
           template <typename CVType> void SetSubColumn(CVECTOREXPRESSION<Type,CVType> const&,int,int);

           //Type*& operator[](int i);      // Set or Return an element
           //Type* operator[](int i) const; // Return an element
           std::vector<Type>& operator[](int i);      // Set or Return an element
           std::vector<Type> operator[](int i) const; // Return an element

           Type& operator()(int i,int j);      // Set or Return an element
           Type operator()(int i,int j) const; // Return an element

           MATRIX<Type,0,0>& operator=(MATRIX<Type,0,0> const &);
           template <std::size_t m1,std::size_t m2> MATRIX<Type,0,0>& operator=(MATRIX<Type,m1,m2> const&);

           template <typename MType> MATRIX<Type,0,0>& operator=(MATRIXEXPRESSION<Type,MType> const &);
           template <typename MType> MATRIX<Type,0,0>& operator+=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,0,0>& operator-=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,0,0>& operator*=(MATRIXEXPRESSION<Type,MType> const&);
           template <typename MType> MATRIX<Type,0,0>& operator/=(MATRIXEXPRESSION<Type,MType> const&);

           MATRIX<Type,0,0>& operator=(Type const&);
           MATRIX<Type,0,0>& operator+=(Type const&);
           MATRIX<Type,0,0>& operator-=(Type const&);
           MATRIX<Type,0,0>& operator*=(Type const&);
           MATRIX<Type,0,0>& operator/=(Type const&);

           template <typename MType> MATRIX<Type,0,0>& SubMatrixAddition(MATRIXEXPRESSION<Type,MType> const&,int,int);
           template <typename MType> MATRIX<Type,0,0>& SubMatrixSubtraction(MATRIXEXPRESSION<Type,MType> const&,int,int);

           template <typename MType> bool operator==(MATRIXEXPRESSION<Type,MType> const&) const;
           template <typename MType> bool operator!=(MATRIXEXPRESSION<Type,MType> const&) const;
       	 };

//*************************************************************************
//*************************************************************************
//*************************************************************************

template <typename Type,typename MType> class MATRIXEXPRESSION_NEGATE : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_NEGATE<Type,MType> >
      { private : MType const &m;
        public : MATRIXEXPRESSION_NEGATE(MType const &M) : m(M) {;}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { return -m(i,j);} 
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXTRANSPOSE : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXTRANSPOSE<Type,MType> >
      { private : MType const &m;
        public : MATRIXEXPRESSION_MATRIXTRANSPOSE(MType const &M) : m(M) {;}
                 std::size_t N1() const { return m.N2();}
                 std::size_t N2() const { return m.N1();}
                 Type operator()(int i,int j) const { return m(j,i);} 
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXANTITRANSPOSE : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXANTITRANSPOSE<Type,MType> >
      { private : MType const &m;
        public : MATRIXEXPRESSION_MATRIXANTITRANSPOSE(MType const &M) : m(M) {;}
                 std::size_t N1() const { return m.N2();}
                 std::size_t N2() const { return m.N1();}
                 Type operator()(int i,int j) const { return m((int)N1()-1-i,(int)N2()-1-j);} 
       };

// ********************
 
template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXADJOINT : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXADJOINT<Type,MType> >
      { private : MType const &m;
        public : MATRIXEXPRESSION_MATRIXADJOINT(MType const &M) : m(M) {;}
                 std::size_t N1() const { return m.N2();}
                 std::size_t N2() const { return m.N1();}
                 Type operator()(int i,int j) const { return conj(m(j,i));} 
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_CONJUGATE : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_CONJUGATE<Type,MType> >
      { private : MType const &m;
        public : MATRIXEXPRESSION_CONJUGATE(MType const &M) : m(M) {;}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { return conj(m(i,j));} 
       }; 

// ********************

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXADDITION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXADDITION<Type,MType1,MType2> >
      { private : MType1 const &m1; MType2 const &m2;
        public : MATRIXEXPRESSION_MATRIXADDITION(MType1 const &M1,MType2 const &M2) : m1(M1), m2(M2) { assert(m1.N1()==m2.N1()); assert(m1.N2()==m2.N2());}
                 std::size_t N1() const { return m1.N1();}
                 std::size_t N2() const { return m1.N2();}
                 Type operator()(int i,int j) const { return m1(i,j)+m2(i,j);} 
       }; 

// ********************

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXSUBTRACTION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXSUBTRACTION<Type,MType1,MType2> >
      { private : MType1 const &m1; MType2 const &m2;
        public : MATRIXEXPRESSION_MATRIXSUBTRACTION(MType1 const &M1,MType2 const &M2) : m1(M1), m2(M2) { assert(m1.N1()==m2.N1()); assert(m1.N2()==m2.N2());}
                 std::size_t N1() const { return m1.N1();}
                 std::size_t N2() const { return m1.N2();}
                 Type operator()(int i,int j) const { return m1(i,j)-m2(i,j);} 
       }; 

// ********************

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXMULTIPLICATION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXMULTIPLICATION<Type,MType1,MType2> >
      { private : MType1 const &m1; MType2 const &m2;
        public : MATRIXEXPRESSION_MATRIXMULTIPLICATION(MType1 const &M1,MType2 const &M2) : m1(M1), m2(M2) { assert(m1.N2()==m2.N1());}
                 std::size_t N1() const { return m1.N1();}
                 std::size_t N2() const { return m2.N2();}
                 Type operator()(int i,int j) const; 
       };

template <typename Type,typename MType1,typename MType2> inline Type MATRIXEXPRESSION_MATRIXMULTIPLICATION<Type,MType1,MType2>::operator()(int i,int j) const
         { Type T=m1(i,0)*m2(0,j);
           int k,kmax=static_cast<int>(m1.N2())-1;
           for(k=1;k<=kmax;k++){ T+=m1(i,k)*m2(k,j);}
           return T;
          }

// ********************

/*template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_MATRIXDIVISION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXDIVISION<Type,MType1,MType2> >
      { private : MType1 const &m1; MType2 const &m2;
        public : MATRIXEXPRESSION_MATRIXDIVISION(MType1 const &M1,MType2 const &M2) : m1(M1), m2(Inverse(M2)) { assert(m1.N2()==m2.N1()); N=m1.N2();}
                 std::size_t N1() const { return m1.N1();}
                 std::size_t N2() const { return m2.N2();}
                 Type operator()(int i,int j) const; 
       };

template <typename Type,typename MType1,typename MType2> Type MATRIXEXPRESSION_MATRIXDIVISION<Type,MType1,MType2>::operator()(int i,int j) const
         { Type T=Zero(m1(i,0)*m2(0,j));
           int k,kmax=static_cast<int>(m1.N2())-1;
           for(k=0;k<=kmax;k++){ T+=m1(i,k)*m2(k,j);} 
           return T;
          }*/

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARADDITION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXSCALARADDITION<Type,MType> >
      { private : MType const &m; Type const &s;
        public : MATRIXEXPRESSION_MATRIXSCALARADDITION(MType const &M,Type const&S) : m(M), s(S) { assert(m.N1()==m.N2());}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { if(i==j){ return m(i,j)+s;} else{ return m(i,j);} } 
       };

template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXADDITION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_SCALARMATRIXADDITION<Type,MType> >
      { private : Type const &s; MType const &m;
        public : MATRIXEXPRESSION_SCALARMATRIXADDITION(Type const&S,MType const &M) : s(S), m(M) { assert(m.N1()==m.N2());}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { if(i==j){ return s+m(i,j);} else{ return m(i,j);} }
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION<Type,MType> >
      { private : MType const &m; Type const &s;
        public : MATRIXEXPRESSION_MATRIXSCALARSUBTRACTION(MType const &M,Type const&S) : m(M), s(S) { assert(m.N1()==m.N2());}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { if(i==j){ return m(i,j)-s;} else{ return m(i,j);} } 
       };

template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION<Type,MType> >
      { private : Type const &s; MType const &m;
        public : MATRIXEXPRESSION_SCALARMATRIXSUBTRACTION(Type const&S,MType const &M) : s(S), m(M) { assert(m.N1()==m.N2());}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { if(i==j){ return s-m(i,j);} else{ return m(i,j);} }
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION<Type,MType> >
      { private : MType const &m; Type const &s;
        public : MATRIXEXPRESSION_MATRIXSCALARMULTIPLICATION(MType const &M,Type const&S) : m(M), s(S) {;}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { return m(i,j)*s;} 
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION<Type,MType> >
      { private : Type const &s; MType const &m;
        public : MATRIXEXPRESSION_SCALARMATRIXMULTIPLICATION(Type const&S,MType const &M) : s(S), m(M) {;}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { return s*m(i,j);} 
       };

// ********************

template <typename Type,typename MType> class MATRIXEXPRESSION_MATRIXSCALARDIVISION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_MATRIXSCALARDIVISION<Type,MType> >
      { private : MType const &m; Type const &s; 
        public : MATRIXEXPRESSION_MATRIXSCALARDIVISION(MType const &M,Type const&S) : m(M), s(S) {;}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { return m(i,j)/s;} 
       };

/*template <typename Type,typename MType> class MATRIXEXPRESSION_SCALARMATRIXDIVISION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_SCALARMATRIXDIVISION<Type,MType> >
      { private : Type const&s; MType const &m;
        public : MATRIXEXPRESSION_SCALARMATRIXDIVISION(Type const&S,MType const &M) : s(S), m(Inverse(M)) {;}
                 std::size_t N1() const { return m.N1();}
                 std::size_t N2() const { return m.N2();}
                 Type operator()(int i,int j) const { return s*m(i,j);} 
       };*/

// ********************

template <typename Type,typename CVType,typename RVType> class MATRIXEXPRESSION_VECTORPRODUCT : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_VECTORPRODUCT<Type,CVType,RVType> >
      { private : CVType const &cv; RVType const &rv;
        public : MATRIXEXPRESSION_VECTORPRODUCT(CVType const&CV,RVType const&RV) : cv(CV), rv(RV) {;}
                 std::size_t N1() const { return cv.N();}
                 std::size_t N2() const { return rv.N();}
                 Type operator()(int i,int j) const { return cv[i]*rv[j];} 
       };

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_SUBMATRIXADDITION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_SUBMATRIXADDITION<Type,MType1,MType2> >
      { private : MType1 const &m1; MType2 const &m2; const int i,j;
        public : MATRIXEXPRESSION_SUBMATRIXADDITION(MType1 const &M1,MType2 const &M2,int I,int J) : m1(M1), m2(M2), i(I), j(J) { assert((int)m1.N1()==i+(int)m2.N1()); assert(m1.N2()==j+m2.N2());}
                 std::size_t N1() const { return m1.N1();}
                 std::size_t N2() const { return m1.N2();}
                 Type operator()(int a,int b) const { if(a>=i && a<=(int)m1.N1()-(int)m2.N1() && b>=j && b<=(int)m1.N2()-(int)m2.N2()){ return m1(a,b)+m2(a-i,b-j);} else{ return m1(a,b);} }  
       };

template <typename Type,typename MType1,typename MType2> class MATRIXEXPRESSION_SUBMATRIXSUBTRACTION : public MATRIXEXPRESSION<Type,MATRIXEXPRESSION_SUBMATRIXSUBTRACTION<Type,MType1,MType2> >
      { private : MType1 const &m1; MType2 const &m2; const int i,j;
        public : MATRIXEXPRESSION_SUBMATRIXSUBTRACTION(MType1 const &M1,MType2 const &M2,int I,int J) : m1(M1), m2(M2), i(I), j(J) { assert((int)m1.N1()==i+(int)m2.N1()); assert((int)m1.N2()==j+(int)m2.N2());}
                 std::size_t N1() const { return m1.N1();}
                 std::size_t N2() const { return m1.N2();}
                 Type operator()(int a,int b) const { if(a>=i && b>=j && b>=j && b<=(int)m1.N2()-(int)m2.N2()){ return m1(a,b)-m2(a-i,b-j);} else{ return m1(a,b);} }  
       };

//*************************************************************************
//*************************************************************************
//*************************************************************************

#endif

