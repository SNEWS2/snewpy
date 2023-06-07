
#include <algorithm>
#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mstl.h"

// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************

#if !defined(_INTERPOLATION_BASE)
#define _INTERPOLATION_BASE

namespace interpolation{

template <typename IType> class EXPRESSION;

// *************************************************************************************************

template <typename IType>
class EXPRESSION 
      { public :         

        virtual ~EXPRESSION(void) {;}

        double operator()(double T) const { return reinterpret_cast<IType const&>(*this)(T);}
        double Interpolate(double T) const { return reinterpret_cast<IType const&>(*this).Interpolate(T);}

        int N(void) const { return reinterpret_cast<IType const&>(*this).N();}

        double X(int i) const { return reinterpret_cast<IType const&>(*this).X(i);}
        std::vector<double> X(void) const { return reinterpret_cast<IType const&>(*this).X();}

        double XMin(void) const { return reinterpret_cast<IType const&>(*this).XMin();}
        double XMax(void) const { return reinterpret_cast<IType const&>(*this).XMax();}

        int XInterval(double X) const { return reinterpret_cast<IType const&>(*this).XInterval(X);}

        operator IType const& () const { return reinterpret_cast<IType const&>(*this);}
       };

// **************************

template <typename IType> class EXPRESSION_NEGATEEXPRESSION;

template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_ADDITION;
template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION;
template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION;
template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_DIVISION;

template <typename IType> class EXPRESSION_SCALAREXPRESSION_ADDITION;
template <typename IType> class EXPRESSION_SCALAREXPRESSION_SUBTRACTION;
template <typename IType> class EXPRESSION_SCALAREXPRESSION_MULTIPLICATION;

// **************************

template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION;
template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION;
template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION;
template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION;

template <typename IType> class EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION;
template <typename IType> class EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION;
template <typename IType> class EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION;

// **************************

template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION;
template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION;
template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION;
template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION;

template <typename IType> class EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION;
template <typename IType> class EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION;
template <typename IType> class EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION;

// **************************

template <typename IType1,typename IType2> class EXPRESSION_EXPRESSIONEXPRESSION_ADDITION;
template <typename IType1,typename IType2> class EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION;
template <typename IType1,typename IType2> class EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION;
template <typename IType1,typename IType2> class EXPRESSION_EXPRESSIONEXPRESSION_DIVISION;

// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************

template <typename IType> EXPRESSION_NEGATEEXPRESSION<EXPRESSION<IType> > operator-(const EXPRESSION<IType>&);

template <typename IType> EXPRESSION_EXPRESSIONSCALAR_ADDITION<EXPRESSION<IType> > operator+(const EXPRESSION<IType>&,const double&);
template <typename IType> EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION<EXPRESSION<IType> > operator-(const EXPRESSION<IType>&,const double&);
template <typename IType> EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION<EXPRESSION<IType> > operator*(const EXPRESSION<IType>&,const double&);
template <typename IType> EXPRESSION_EXPRESSIONSCALAR_DIVISION<EXPRESSION<IType> > operator/(const EXPRESSION<IType>&,const double&);

template <typename IType> EXPRESSION_SCALAREXPRESSION_ADDITION<EXPRESSION<IType> > operator+(const double&,const EXPRESSION<IType>&);
template <typename IType> EXPRESSION_SCALAREXPRESSION_SUBTRACTION<EXPRESSION<IType> > operator-(const double&,const EXPRESSION<IType>&);
template <typename IType> EXPRESSION_SCALAREXPRESSION_MULTIPLICATION<EXPRESSION<IType> > operator*(const double&,const EXPRESSION<IType>&);

template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION<EXPRESSION<IType> > operator+(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double,double(*)(void)> &NF); // C(x) + NF()
template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION<EXPRESSION<IType> > operator-(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double,double(*)(void)> &NF); // C(x) - NF()
template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION<EXPRESSION<IType> > operator*(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double,double(*)(void)> &NF); // C(x) * NF()
template <typename IType> EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION<EXPRESSION<IType> > operator/(const EXPRESSION<IType> &C,const NULLARYFUNCTOR<double,double(*)(void)> &NF); // C(x) / NF()

template <typename IType> EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION<EXPRESSION<IType> > operator+(const NULLARYFUNCTOR<double,double(*)(void)> &NF,const EXPRESSION<IType> &C); // NF() + C(x)
template <typename IType> EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION<EXPRESSION<IType> > operator-(const NULLARYFUNCTOR<double,double(*)(void)> &NF,const EXPRESSION<IType> &C); // NF() - C(x)
template <typename IType> EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION<EXPRESSION<IType> > operator*(const NULLARYFUNCTOR<double,double(*)(void)> &NF,const EXPRESSION<IType> &C); // NF() * C(x)

template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION<EXPRESSION<IType> > operator+(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double,double,double(*)(double)> &UF); // C(x) + UF(x)
template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION<EXPRESSION<IType> > operator-(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double,double,double(*)(double)> &UF); // C(x) - UF(x)
template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION<EXPRESSION<IType> > operator*(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double,double,double(*)(double)> &UF); // C(x) * UF(x)
template <typename IType> EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION<EXPRESSION<IType> > operator/(const EXPRESSION<IType> &C,const UNARYFUNCTOR<double,double,double(*)(double)> &UF); // C(x) / UF(x)

template <typename IType> EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION<EXPRESSION<IType> > operator+(const UNARYFUNCTOR<double,double,double(*)(double)> &UF,const EXPRESSION<IType> &C); // UF(x) + C(x)
template <typename IType> EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION<EXPRESSION<IType> > operator-(const UNARYFUNCTOR<double,double,double(*)(double)> &UF,const EXPRESSION<IType> &C); // UF(x) - C(x)
template <typename IType> EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION<EXPRESSION<IType> > operator*(const UNARYFUNCTOR<double,double,double(*)(double)> &UF,const EXPRESSION<IType> &C); // UF(x) * C(x)

// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************

template <typename IType1,typename IType2>
std::vector<double> OverlapX(EXPRESSION<IType1> const &C1,EXPRESSION<IType1> const &C2);

template <typename IType1,typename IType2>
EXPRESSION_EXPRESSIONEXPRESSION_ADDITION<EXPRESSION<IType1>,EXPRESSION<IType2> > 
operator+(EXPRESSION<IType1> const&,EXPRESSION<IType2> const&);  // C1(x) + C2(x)

template <typename IType1,typename IType2>
EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION<EXPRESSION<IType1>,EXPRESSION<IType2> > 
operator-(EXPRESSION<IType1> const&,EXPRESSION<IType2> const&); // C1(x) - C2(x)

template <typename IType1,typename IType2>
EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION<EXPRESSION<IType1>,EXPRESSION<IType2> > 
operator*(EXPRESSION<IType1> const&,EXPRESSION<IType2> const&); // C1(x) * C2(x)

template <typename IType1,typename IType2>
EXPRESSION_EXPRESSIONEXPRESSION_DIVISION<EXPRESSION<IType1>,EXPRESSION<IType2> > 
operator/(EXPRESSION<IType1> const&,EXPRESSION<IType2> const&); // C1(x) / C2(x)

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

template <typename IType> class EXPRESSION_NEGATEEXPRESSION
      : public EXPRESSION<EXPRESSION_NEGATEEXPRESSION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i;
        public : explicit EXPRESSION_NEGATEEXPRESSION(std::vector<double> X,IType const &I) : XDATA_SINGLESET(X), i(I) {;}

                 double operator()(double T) const { return -i(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_ADDITION
      : public EXPRESSION<EXPRESSION_EXPRESSIONSCALAR_ADDITION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_EXPRESSIONSCALAR_ADDITION(std::vector<double> X,IType const &I,double const &D) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return i(T) + d;}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION
      : public EXPRESSION<EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_EXPRESSIONSCALAR_SUBTRACTION(std::vector<double> X,IType const &I,double const &D) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return i(T) - d;}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION
      : public EXPRESSION<EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_EXPRESSIONSCALAR_MULTIPLICATION(std::vector<double> X,IType const &I,double const &D) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return i(T) * d;} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_EXPRESSIONSCALAR_DIVISION
      : public EXPRESSION<EXPRESSION_EXPRESSIONSCALAR_DIVISION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_EXPRESSIONSCALAR_DIVISION(std::vector<double> X,IType const &I,double const &D) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return i(T) / d;} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

// **********************************************************

template <typename IType> class EXPRESSION_SCALAREXPRESSION_ADDITION
      : public EXPRESSION<EXPRESSION_SCALAREXPRESSION_ADDITION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_SCALAREXPRESSION_ADDITION(std::vector<double> X,double const &D,IType const &I) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return d + i(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_SCALAREXPRESSION_SUBTRACTION
      : public EXPRESSION<EXPRESSION_SCALAREXPRESSION_SUBTRACTION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_SCALAREXPRESSION_SUBTRACTION(std::vector<double> X,double const &D,IType const &I) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return d - i(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_SCALAREXPRESSION_MULTIPLICATION
      : public EXPRESSION<EXPRESSION_SCALAREXPRESSION_MULTIPLICATION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; double const &d;
        public : explicit EXPRESSION_SCALAREXPRESSION_MULTIPLICATION(std::vector<double> X,double const &D,IType const &I) : XDATA_SINGLESET(X), i(I), d(D) {;}

                 double operator()(double T) const { return d * i(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

// **************************

template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION
      : public EXPRESSION<EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_EXPRESSIONNULLARYFUNCTOR_ADDITION(std::vector<double> X,IType const &I,NULLARYFUNCTOR<double,double(*)(void)> const &NF) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return i(T) + nf();} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION
      : public EXPRESSION<EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_EXPRESSIONNULLARYFUNCTOR_SUBTRACTION(std::vector<double> X,IType const &I,NULLARYFUNCTOR<double,double(*)(void)> const &NF) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return i(T) - nf();}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION
      : public EXPRESSION<EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_EXPRESSIONNULLARYFUNCTOR_MULTIPLICATION(std::vector<double> X,IType const &I,NULLARYFUNCTOR<double,double(*)(void)> const &NF) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return i(T) * nf();} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION
      : public EXPRESSION<EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_EXPRESSIONNULLARYFUNCTOR_DIVISION(std::vector<double> X,IType const &I,NULLARYFUNCTOR<double,double(*)(void)> const &NF) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return i(T) / nf();}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION
      : public EXPRESSION<EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_NULLARYFUNCTOREXPRESSION_ADDITION(std::vector<double> X,NULLARYFUNCTOR<double,double(*)(void)> const &NF,IType const &I) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return nf() + i(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION
      : public EXPRESSION<EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_NULLARYFUNCTOREXPRESSION_SUBTRACTION(std::vector<double> X,NULLARYFUNCTOR<double,double(*)(void)> const &NF,IType const &I) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return nf() - i(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION
      : public EXPRESSION<EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; NULLARYFUNCTOR<double,double(*)(void)> const &nf;
        public : explicit EXPRESSION_NULLARYFUNCTOREXPRESSION_MULTIPLICATION(std::vector<double> X,NULLARYFUNCTOR<double,double(*)(void)> const &NF,IType const &I) : XDATA_SINGLESET(X), i(I), nf(NF) {;}

                 double operator()(double T) const { return nf() * i(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

// **************************

template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION
      : public EXPRESSION<EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_EXPRESSIONUNARYFUNCTOR_ADDITION(std::vector<double> X,IType const &I,UNARYFUNCTOR<double,double,double(*)(double)> const &UF) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return i(T) + uf(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION
      : public EXPRESSION<EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_EXPRESSIONUNARYFUNCTOR_SUBTRACTION(std::vector<double> X,IType const &I,UNARYFUNCTOR<double,double,double(*)(double)> const &UF) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return i(T) - uf(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION
      : public EXPRESSION<EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_EXPRESSIONUNARYFUNCTOR_MULTIPLICATION(std::vector<double> X,IType const &I,UNARYFUNCTOR<double,double,double(*)(double)> const &UF) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return i(T) * uf(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType> class EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION
      : public EXPRESSION<EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_EXPRESSIONUNARYFUNCTOR_DIVISION(std::vector<double> X,IType const &I,UNARYFUNCTOR<double,double,double(*)(double)> const &UF) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return i(T) / uf(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION
      : public EXPRESSION<EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION<IType> >
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_UNARYFUNCTOREXPRESSION_ADDITION(std::vector<double> X,UNARYFUNCTOR<double,double,double(*)(double)> const &UF,IType const &I) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return uf(T) + i(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION
      : public EXPRESSION<EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_UNARYFUNCTOREXPRESSION_SUBTRACTION(std::vector<double> X,UNARYFUNCTOR<double,double,double(*)(double)> const &UF,IType const &I) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return uf(T) - i(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);} 
       };

template <typename IType> class EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION
      : public EXPRESSION<EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION<IType> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType const &i; UNARYFUNCTOR<double,double,double(*)(double)> const &uf;
        public : explicit EXPRESSION_UNARYFUNCTOREXPRESSION_MULTIPLICATION(std::vector<double> X,UNARYFUNCTOR<double,double,double(*)(double)> const &UF,IType const &I) : XDATA_SINGLESET(X), i(I), uf(UF) {;}

                 double operator()(double T) const { return uf(T) * i(T);}  

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

template <typename IType1,typename IType2> 
class EXPRESSION_EXPRESSIONEXPRESSION_ADDITION 
      : public EXPRESSION<EXPRESSION_EXPRESSIONEXPRESSION_ADDITION<IType1,IType2> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType1 const &i1; IType2 const &i2;
        public : explicit EXPRESSION_EXPRESSIONEXPRESSION_ADDITION(std::vector<double> X,IType1 const &I1,IType2 const &I2) : XDATA_SINGLESET(X), i1(I1), i2(I2) {;}

                 double operator()(double T) const { return i1(T) + i2(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType1,typename IType2> 
class EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION 
      : public EXPRESSION<EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION<IType1,IType2> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType1 const &i1; IType2 const &i2;
        public : explicit EXPRESSION_EXPRESSIONEXPRESSION_SUBTRACTION(std::vector<double> X,IType1 const &I1,IType2 const &I2) : XDATA_SINGLESET(X), i1(I1), i2(I2) {;}

                 double operator()(double T) const { return i1(T) - i2(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType1,typename IType2> 
class EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION 
      : public EXPRESSION<EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION<IType1,IType2> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType1 const &i1; IType2 const &i2;
        public : explicit EXPRESSION_EXPRESSIONEXPRESSION_MULTIPLICATION(std::vector<double> X,IType1 const &I1,IType2 const &I2) : XDATA_SINGLESET(X), i1(I1), i2(I2) {;}

                 double operator()(double T) const { return i1(T) * i2(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

template <typename IType1,typename IType2> 
class EXPRESSION_EXPRESSIONEXPRESSION_DIVISION 
      : public EXPRESSION<EXPRESSION_EXPRESSIONEXPRESSION_DIVISION<IType1,IType2> >, virtual public XDATA_SINGLESET, virtual public XLIMITS_SINGLESET, virtual public XINTERVAL_SINGLESET
      { private : IType1 const &i1; IType2 const &i2;
        public : explicit EXPRESSION_EXPRESSIONEXPRESSION_DIVISION(std::vector<double> X,IType1 const &I1,IType2 const &I2) : XDATA_SINGLESET(X), i1(I1), i2(I2) {;}

                 double operator()(double T) const { return i1(T) / i2(T);} 

                 int N(void) const { return XDATA_SINGLESET::NX();}

                 double X(int i) const { return XDATA_SINGLESET::X(i);}
                 std::vector<double> X(void) const { return XDATA_SINGLESET::X();}

                 double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
                 double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

                 int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}
       };

} // end of namespace interpolation

#endif

