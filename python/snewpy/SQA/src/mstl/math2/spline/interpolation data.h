
#include <algorithm>
#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mstl.h"

// *******************************************************************
// *******************************************************************
// *******************************************************************

#if !defined(_INTERPOLATION_DATA)
#define _INTERPOLATION_DATA

namespace interpolation{

class XDATA_SINGLESET;

template<int ND> class YDATA_SINGLESET; // the data is ND dimensional

class DELTAX_SINGLESET;

class DELTAY_SINGLESET; 

class GRADEBASE_SINGLEXSET_SINGLEYSET;
class LINEGRADE_SINGLEXSET_SINGLEYSET;
class LOCALGRADE_SINGLEXSET_SINGLEYSET;

class SPLINE_SINGLEXSET_SINGLEYSET;

class XLIMITS_SINGLESET;

class YLIMITS_SINGLESET;
      
class XYLIMITS_SINGLEXSET_SINGLEYSET;

class SORT_SINGLEXSET_SINGLEYSET;

class OPENWRITE_SINGLEXSET_SINGLEYSET;

class XINTERVAL_SINGLESET;
class YINTERVAL_SINGLESET;

class DISTINGUISH_SINGLEYSET;

class OPERATOR;

// ****************************************************************
// ****************************************************************
// ****************************************************************

class XDATA_SINGLESET
           { private :

             static std::string CLASS_NAME;

             mutable std::vector<double> x;  // x is the data

             // *******************
             // *******************
             // *******************

             protected :

             int NX(void) const { return x.size();}
             
             double _pX(int i) const { return x[i-1];}
             std::vector<double> _pX(void) const { return x;}

             void _pSetX(int i,double X) const { x[i-1]=X;}
             void _pSetX(const std::vector<double> &X) const { x=X;}

             // *******************

             void CreateX(int NX) const { x=std::vector<double>(NX);}
             void CreateX(void) const {;}

             void CopyX(const XDATA_SINGLESET &C) const { x=C.x;}

             void DestroyX(void) const { x.clear();}

             // ********************
             // ********************
             // ********************

	     public :

             XDATA_SINGLESET(void) { CreateX();}
             XDATA_SINGLESET(int NX) { CreateX(NX);}             
             XDATA_SINGLESET(const std::vector<double> &X) { _pSetX(X);}
             XDATA_SINGLESET(int NX,const double *X);
             XDATA_SINGLESET(const XDATA_SINGLESET &C) { CopyX(C);}
             
             virtual ~XDATA_SINGLESET(void){ DestroyX();}

             // *************************************

             bool Empty(void) const { return x.empty();}

             // retrieve or set the data
             double X(int i) const;
             std::vector<double> X(void) const { return _pX();}

             void SetX(int i,double X) const;
             void SetX(const std::vector<double> &X) const { _pSetX(X);}

             // extend the data
             void AddX(double X);
             void AddX(const std::vector<double> &X);

             // transform the data 
             template<typename UNARYFUNCTOR> void TransformX(const UNARYFUNCTOR &UF);

             void ShiftX(const double&);
             void RescaleX(const double&);

             void AbsX(void);

             // operators
             XDATA_SINGLESET& operator=(const XDATA_SINGLESET &C){ CopyX(C); return *this;}
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<int ND> class YDATA_SINGLESET
           { private :

             static std::string CLASS_NAME;
                        
             mutable TARRAY<double,ND> y;  // y is the data

             // *******************

             template<typename UNARYFUNCTOR> void _pTransformY(std::vector<int> i,int j,const UNARYFUNCTOR &UF);
             void _pShiftY(std::vector<int> i,int j,const double&);
             void _pRescaleY(std::vector<int> i,int j,const double&);
             void _pAbsY(std::vector<int> i,int j);

             // *******************
             // *******************
             // *******************

             protected :

             std::vector<int> NY(void) const { return y.Dimensions();}
             int NY(int j) const { return y.Dimension(j);}
             
             double _pY(std::vector<int> i) const;
             TARRAY<double,ND> _pY(void) const { return y;}

             void _pSetY(std::vector<int> i,double Y) const;
             void _pSetY(const TARRAY<double,ND> &Y) const { y=Y;}

             // *******************

             void CreateY(std::vector<int> NY) const { y=TARRAY<double,ND>(NY);}     
        
             void CopyY(const YDATA_SINGLESET &C) const { y=C.y;}

             void DestroyY(void) const {;} 

             // ********************
             // ********************
             // ********************

	     public :

	     // void constructor
             YDATA_SINGLESET(void){;}
             
             // set n DATA
             YDATA_SINGLESET(std::vector<int> NY){ CreateY(NY);}             
                          
             // construct from data in containers
             YDATA_SINGLESET(const TARRAY<double,ND> &Y) { _pSetY(Y);}

             YDATA_SINGLESET(const YDATA_SINGLESET &C) { CopyY(C);}
             
             virtual ~YDATA_SINGLESET(void){ DestroyY();}

             // *************************************

             bool Empty(void) const { return y.Empty();}

             // retrieve or set the data
             double Y(std::vector<int> i) const;
             TARRAY<double,ND> Y(void) const { return _pY();}

             void SetY(std::vector<int> i,double Y) const;
             void SetY(const TARRAY<double,ND> &Y) const { _pSetY(Y);}

             // transform the data 
             template<typename UNARYFUNCTOR> void TransformY(const UNARYFUNCTOR &UF);
             void ShiftY(const double&);
             void RescaleY(const double&);
             void AbsY(void);

             // operators
             YDATA_SINGLESET& operator=(const YDATA_SINGLESET &C){ CopyY(C); return *this;}
            };

//********************************************************************************
//********************************************************************************
//********************************************************************************

template<> class YDATA_SINGLESET<1>
           { private :

             static std::string CLASS_NAME;
                        
             mutable std::vector<double> y;  // y is the data

             // *******************
             // *******************
             // *******************

             protected :

             int NY(void) const { return y.size();}
             
             double _pY(int i) const { return y[i-1];}
             std::vector<double> _pY(void) const { return y;}

             void _pSetY(int i,double Y) const { y[i-1]=Y;}
             void _pSetY(const std::vector<double> &Y) const { y=Y;}

             // *******************

             void CreateY(int NY) const { y=std::vector<double>(NY);}

             void CopyY(const YDATA_SINGLESET &C) const { y=C.y;}

             void DestroyY(void) const { y.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

	     // void constructor
             YDATA_SINGLESET(void){ CreateY(0);}
             
             // set n DATA
             YDATA_SINGLESET(int NY){ CreateY(NY);}             
                          
             // construct from data in containers
             YDATA_SINGLESET(const std::vector<double> &Y) { _pSetY(Y);}
             YDATA_SINGLESET(int NY,const double *Y);

             YDATA_SINGLESET(const YDATA_SINGLESET &C) { CopyY(C);}
             
             virtual ~YDATA_SINGLESET(void){ DestroyY();}

             // *************************************

             bool Empty(void) const { return y.empty();}

             // retrieve or set the data
             double Y(int i) const;
             std::vector<double> Y(void) const { return _pY();}

             void SetY(int i,double Y) const;
             void SetY(const std::vector<double> &Y) const { _pSetY(Y);}

             // extend the data
             void AddY(double Y);
             void AddY(const std::vector<double> &Y);

             // transform the data 
             template<typename UNARYFUNCTOR> void TransformY(const UNARYFUNCTOR &UF);

             void ShiftY(const double&);
             void RescaleY(const double&);

             void AbsY(void);

             // operators
             YDATA_SINGLESET& operator=(const YDATA_SINGLESET &C){ CopyY(C); return *this;}
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class DELTAX_SINGLESET 
      : virtual public XDATA_SINGLESET
           { private :

             mutable bool xdifferenced;
             mutable std::vector<double> deltax;  // deltax are the differences between data points

             // *******************
             // *******************
             // *******************

             protected :

             int NDeltaX(void) const { return deltax.size();}
             
             double _pDeltaX(int i) const { return deltax[i-1];}
             std::vector<double> _pDeltaX(void) const { return deltax;}

             void _pSetDeltaX(int i,double DELTAX) const { deltax[i-1]=DELTAX;}
             void _pSetDeltaX(const std::vector<double> &DELTAX) const { deltax=DELTAX;}

             bool XDifferenced(void) const { return xdifferenced;}
             void SetXDifferenced(bool XDIFFERENCED) const { xdifferenced=XDIFFERENCED;}

             // *******************

             void CreateDeltaX(int NDeltaX) const { xdifferenced=false; deltax=std::vector<double>(NDeltaX);}
             void CreateDeltaX(void) const { xdifferenced=false;}

             void CopyDeltaX(const DELTAX_SINGLESET &C) const { xdifferenced=C.xdifferenced; deltax=C.deltax;}

             void DestroyDeltaX(void) const { xdifferenced=false; deltax.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

             DELTAX_SINGLESET(void) { CreateDeltaX();}
             DELTAX_SINGLESET(int NDeltaX) { CreateDeltaX(NDeltaX);}
             DELTAX_SINGLESET(const DELTAX_SINGLESET &C) { CopyDeltaX(C);}
             
             virtual ~DELTAX_SINGLESET(void) { DestroyDeltaX();}

             void XDifference(void) const;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class DELTAY_SINGLESET
      : virtual public YDATA_SINGLESET<1>
           { private :

             mutable bool ydifferenced;
             mutable std::vector<double> deltay;  // deltay are the differences between data points

             // *******************
             // *******************
             // *******************

             protected :

             int NDeltaY(void) const { return deltay.size();}
             
             double _pDeltaY(int i) const { return deltay[i-1];}
             std::vector<double> _pDeltaY(void) const { return deltay;}

             void _pSetDeltaY(int i,double DELTAY) const { deltay[i-1]=DELTAY;}
             void _pSetDeltaY(const std::vector<double> &DELTAY) const { deltay=DELTAY;}

             bool YDifferenced(void) const { return ydifferenced;}
             void SetYDifferenced(bool YDIFFERENCED) const { ydifferenced=YDIFFERENCED;}

             // *******************

             void CreateDeltaY(int NDeltaY) const { ydifferenced=false; deltay=std::vector<double>(NDeltaY);}

             void CopyDeltaY(const DELTAY_SINGLESET &C) const { ydifferenced=C.ydifferenced; deltay=C.deltay;}

             void DestroyDeltaY(void) const { ydifferenced=false; deltay.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

             DELTAY_SINGLESET(void){ CreateDeltaY(0);}
             DELTAY_SINGLESET(int NDeltaY){ CreateDeltaY(NDeltaY);}
             DELTAY_SINGLESET(const DELTAY_SINGLESET &C) { CopyDeltaY(C);}
             
             virtual ~DELTAY_SINGLESET(void) { DestroyDeltaY();}

             void YDifference(void) const;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class GRADEBASE_SINGLEXSET_SINGLEYSET
           { private :

             static std::string CLASS_NAME;
                        
             mutable bool graded;
             mutable std::vector<double> g;  // g are the gradients

             // *******************
             // *******************
             // *******************

             protected :

             int NG(void) const { return g.size();}
             
             double _pG(int i) const { return g[i-1];}
             std::vector<double> _pG(void) const { return g;}

             void _pSetG(int i,double G) const { g[i-1]=G;}
             void _pSetG(const std::vector<double> &G) const { g=G;}

             bool Graded(void) const { return graded;}
             void SetGraded(bool GRADED) const { graded=GRADED;}

             // *******************

             void CreateG(int NG) const { graded=false; g=std::vector<double>(NG);}

             void CopyG(const GRADEBASE_SINGLEXSET_SINGLEYSET &C) const { graded=C.graded; g=C.g;}

             void DestroyG(void) const { graded=false; g.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

	     // void constructor
             GRADEBASE_SINGLEXSET_SINGLEYSET(void){ CreateG(0);}
             GRADEBASE_SINGLEXSET_SINGLEYSET(int NG){ CreateG(NG);}
             GRADEBASE_SINGLEXSET_SINGLEYSET(const GRADEBASE_SINGLEXSET_SINGLEYSET &C) { CopyG(C);}
             
             virtual ~GRADEBASE_SINGLEXSET_SINGLEYSET(void) { DestroyG();}

             // *************************************

             virtual void Grade(void) const =0;

             // retrieve or set the data
             double G(int i) const;
             std::vector<double> G(void) const { return g;}

             void SetG(int i,double G) const;
             void SetG(const std::vector<double> &G) const { _pSetG(G);}
            };

// *************************************

class LINEGRADE_SINGLEXSET_SINGLEYSET 
      : virtual public GRADEBASE_SINGLEXSET_SINGLEYSET, 
        virtual public DELTAX_SINGLESET, virtual public DELTAY_SINGLESET
      { public :

        LINEGRADE_SINGLEXSET_SINGLEYSET(void) : GRADEBASE_SINGLEXSET_SINGLEYSET() { ;}
        LINEGRADE_SINGLEXSET_SINGLEYSET(int NG) : GRADEBASE_SINGLEXSET_SINGLEYSET(NG) {;}
        LINEGRADE_SINGLEXSET_SINGLEYSET(const LINEGRADE_SINGLEXSET_SINGLEYSET &C) : GRADEBASE_SINGLEXSET_SINGLEYSET(C) {;}

        void Grade(void) const;
       };

// *************************************

class LOCALGRADE_SINGLEXSET_SINGLEYSET 
      : virtual public GRADEBASE_SINGLEXSET_SINGLEYSET, 
        virtual public XDATA_SINGLESET, virtual public DELTAX_SINGLESET, 
        virtual public YDATA_SINGLESET<1>, virtual public DELTAY_SINGLESET
      { public :

        LOCALGRADE_SINGLEXSET_SINGLEYSET(void) : GRADEBASE_SINGLEXSET_SINGLEYSET() { ;}
        LOCALGRADE_SINGLEXSET_SINGLEYSET(int NG) : GRADEBASE_SINGLEXSET_SINGLEYSET(NG) {;}
        LOCALGRADE_SINGLEXSET_SINGLEYSET(const LOCALGRADE_SINGLEXSET_SINGLEYSET &C) : GRADEBASE_SINGLEXSET_SINGLEYSET(C) {;}

        virtual ~LOCALGRADE_SINGLEXSET_SINGLEYSET (void) {;}

        void Grade(void) const;
       };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class SPLINE_SINGLEXSET_SINGLEYSET
           { private :

             mutable bool fitted;
             // note that the second index is *not* offset by unity
             mutable std::vector<std::vector<double> > a; // a are the spline parameters

             // *******************
             // *******************
             // *******************

             protected :

             int NA(void) const { return a.size();}
             int NParameters(int i) const { return a[i-1].size();}

             double A(int i,int j) const { return a[i-1][j];}
             std::vector<double> A(int i) const { return a[i-1];}
             std::vector<std::vector<double> > A(void) const { return a;}

             void SetA(int i,int j,double A) const { a[i-1][j]=A;}
             void SetA(int i,std::vector<double> A) const { a[i-1]=A;}
             void SetA(std::vector<std::vector<double> > A) const { a=A;}

             bool Fitted(void) const { return fitted;}
             void SetFitted(bool FITTED) const { fitted=FITTED;}

             // *******************

             void CreateA(int NA,int NParameters) const { fitted=false; a=std::vector<std::vector<double> >(NA,std::vector<double>(NParameters));}
             void CreateA(int NA) const { fitted=false; a=std::vector<std::vector<double> >(NA);}
             void CreateA(void) const { fitted=false;}

             void CreateASet(int i,int NParameters) const { fitted=false; a[i-1]=std::vector<double>(NParameters);}

             void CopyA(const SPLINE_SINGLEXSET_SINGLEYSET &C) const { fitted=C.fitted; a=C.a;}

             void DestroyA(void) const;

             // ********************

             double Interpolate(int i,double T) const;
             double Derivative(int i,double T) const;
             double SecondDerivative(int i,double T) const;
             double Integral(int i,double Tmin,double Tmax) const;

             // ********************
             // ********************

	     public :

             SPLINE_SINGLEXSET_SINGLEYSET(void){ CreateA();}
             SPLINE_SINGLEXSET_SINGLEYSET(int NA){ CreateA(NA);}                          
             SPLINE_SINGLEXSET_SINGLEYSET(int NA,int NParameters){ CreateA(NA,NParameters);}
             SPLINE_SINGLEXSET_SINGLEYSET(const std::vector<std::vector<double> > &A) { SetA(A);}
             SPLINE_SINGLEXSET_SINGLEYSET(const SPLINE_SINGLEXSET_SINGLEYSET &C) { CopyA(C);}
             
             virtual ~SPLINE_SINGLEXSET_SINGLEYSET(void) { DestroyA();}

             virtual void Fit(void) const =0;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class XLIMITS_SINGLESET
      : virtual public XDATA_SINGLESET
           { private :
             
             mutable double xmin,xmax;
             mutable bool xlimited;

             // ********************
             // ********************
             // ********************

             protected :

             double _pXMin(void) const { return xmin;}
             double _pXMax(void) const { return xmax;}

             void SetXMin(double XMIN) const { xmin=XMIN;}
             void SetXMax(double XMAX) const { xmax=XMAX;}

             bool XLimited(void) const { return xlimited;}
             void SetXLimited(bool XLIMITED) const { xlimited=XLIMITED;}

             // *******************

             void CreateXLimits(void) const { xlimited=false;}

             void CopyXLimits(const XLIMITS_SINGLESET &C) const { xlimited=C.xlimited; xmin=C.xmin; xmax=C.xmax;}

             void DestroyXLimits(void) const { xlimited=false;} 

             virtual void XLimit(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             XLIMITS_SINGLESET(void) { CreateXLimits();}
             XLIMITS_SINGLESET(const XLIMITS_SINGLESET &C) { CopyXLimits(C);}
             
             virtual ~XLIMITS_SINGLESET(void) { DestroyXLimits();}

             double XMin(void) const;
             double XMax(void) const;
            };

// *************************************
// *************************************
// *************************************

class YLIMITS_SINGLESET
      : virtual public YDATA_SINGLESET<1>
           { private :
             
             mutable double ymin,ymax;
             mutable bool ylimited;

             // ********************
             // ********************
             // ********************

             protected :

             double _pYMin(void) const { return ymin;}
             double _pYMax(void) const { return ymax;}

             void SetYMin(double YMIN) const { ymin=YMIN;}
             void SetYMax(double YMAX) const { ymax=YMAX;}

             bool YLimited(void) const { return ylimited;}
             void SetYLimited(bool YLIMITED) const { ylimited=YLIMITED;}

             // *******************

             void CreateYLimits(void) const { ylimited=false;}

             void CopyYLimits(const YLIMITS_SINGLESET &C) const { ylimited=C.ylimited; ymin=C.ymin; ymax=C.ymax;}

             void DestroyYLimits(void) const { ylimited=false;} 

             virtual void YLimit(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             YLIMITS_SINGLESET(void) { CreateYLimits();}
             YLIMITS_SINGLESET(const YLIMITS_SINGLESET &C) { CopyYLimits(C);}
             
             virtual ~YLIMITS_SINGLESET(void) { DestroyYLimits();}

             double YMin(void) const;
             double YMax(void) const;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class XYLIMITS_SINGLEXSET_SINGLEYSET
      : virtual public XDATA_SINGLESET,
        virtual public YDATA_SINGLESET<1>,
        virtual public XLIMITS_SINGLESET,
        virtual public YLIMITS_SINGLESET
           { private :
             
             mutable double yatxmin,yatxmax;
             mutable double xatymin,xatymax;
             mutable bool xylimited;

             // ********************
             // ********************
             // ********************

             protected :

             double _pYAtXMin(void) const { return yatxmin;}
             double _pYAtXMax(void) const { return yatxmax;}
             double _pXAtYMin(void) const { return xatymin;}
             double _pXAtYMax(void) const { return xatymax;}

             void SetYAtXMin(double YATXMIN) const { yatxmin=YATXMIN;}
             void SetYAtXMax(double YATXMAX) const { yatxmax=YATXMAX;}
             void SetXAtYMin(double XATYMIN) const { xatymin=XATYMIN;}
             void SetXAtYMax(double XATYMAX) const { xatymax=XATYMAX;}

             bool XYLimited(void) const { return xylimited;}
             void SetXYLimited(bool XYLIMITED) const { xylimited=XYLIMITED;}

             // *******************

             void CreateXYLimits(void) const { xylimited=false;}

             void CopyXYLimits(const XYLIMITS_SINGLEXSET_SINGLEYSET &C) const { xylimited=C.xylimited; yatxmin=C.yatxmin; yatxmax=C.yatxmax; xatymin=C.xatymin; xatymax=C.xatymax;}

             void DestroyXYLimits(void) const { xylimited=false;} 

             virtual void XYLimit(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             XYLIMITS_SINGLEXSET_SINGLEYSET(void){ CreateXYLimits();}
             XYLIMITS_SINGLEXSET_SINGLEYSET(const XYLIMITS_SINGLEXSET_SINGLEYSET &C) { CopyXYLimits(C);}
             
             virtual ~XYLIMITS_SINGLEXSET_SINGLEYSET(void) { DestroyXYLimits();}

             double YAtXMin(void) const;
             double YAtXMax(void) const;
             double XAtYMin(void) const;
             double XAtYMax(void) const;
            };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class SORT_SINGLEXSET_SINGLEYSET
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>
      { private :
        mutable bool sorted;

        protected :
        bool Sorted(void) const { return sorted;}
        void SetSorted(bool SORTED) { sorted=SORTED;}

        public :
        SORT_SINGLEXSET_SINGLEYSET(void) { sorted=false;}
        virtual ~SORT_SINGLEXSET_SINGLEYSET(void) { sorted=false;}

        void Sort(void);
       };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class OPENWRITE_SINGLEXSET_SINGLEYSET
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>
      { public :
        OPENWRITE_SINGLEXSET_SINGLEYSET(void) {;}
        OPENWRITE_SINGLEXSET_SINGLEYSET(std::string filename) { Open(filename);}
        OPENWRITE_SINGLEXSET_SINGLEYSET(std::string filename,char ignore) { Open(filename,ignore);}

        virtual ~OPENWRITE_SINGLEXSET_SINGLEYSET(void) {;}

        void Open(std::string filename);
        void Open(std::string filename,char ignore); // all lines beginning with the char 'ignore' will be ignored
        void Write(std::string filename);

        std::ostream& operator<<(std::ostream&);
        std::istream& operator>>(std::istream&);
       };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class XINTERVAL_SINGLESET
      : virtual public XDATA_SINGLESET,
        virtual public XLIMITS_SINGLESET
      { private :
        mutable int interval; // this integer stores previous returns from calls to XInterval so as to make the routine faster

        protected :
        void CreateXInterval(void) const { interval=1;}

        public :
        XINTERVAL_SINGLESET(void) { interval=1;}
        virtual ~XINTERVAL_SINGLESET(void) { interval=1;}

        int XInterval(double X) const;
       };

// *******************************************************************

class YINTERVAL_SINGLESET
      : virtual public YDATA_SINGLESET<1>,
        virtual public YLIMITS_SINGLESET
      { private :
        mutable int interval; // this integer stores previous returns from calls to YInterval so as to make the routine faster

        protected :
        void CreateYInterval(void) const { interval=1;}

        public :
        YINTERVAL_SINGLESET(void) { interval=1;}
        virtual ~YINTERVAL_SINGLESET(void) { interval=1;}

        int YInterval(double Y) const;
       };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class DISTINGUISH_SINGLEYSET
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>, virtual public SORT_SINGLEXSET_SINGLEYSET
      { private :
        mutable bool distinguished;

        protected :
        bool Distinguished(void) const { return distinguished;}
        void SetDistinguished(bool DISTINGUISHED) const { distinguished=DISTINGUISHED;}

        public :
        DISTINGUISH_SINGLEYSET(void) { distinguished=false;}
        virtual ~DISTINGUISH_SINGLEYSET(void) { distinguished=false;}

        virtual int N(void){ return NX();}

        void Distinguish(void);
       };

// ****************************************************************
// ****************************************************************
// ****************************************************************

class OPERATOR
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>
      { public :
        OPERATOR& operator-(void);

        OPERATOR& operator+=(const double &A);
        OPERATOR& operator-=(const double &A);
        OPERATOR& operator*=(const double &A);
        OPERATOR& operator/=(const double &A);

        OPERATOR& operator=(const NULLARYFUNCTOR<double,double(*)(void)> &NF);
        OPERATOR& operator+=(const NULLARYFUNCTOR<double,double(*)(void)> &NF);
        OPERATOR& operator-=(const NULLARYFUNCTOR<double,double(*)(void)> &NF);
        OPERATOR& operator*=(const NULLARYFUNCTOR<double,double(*)(void)> &NF);
        OPERATOR& operator/=(const NULLARYFUNCTOR<double,double(*)(void)> &NF);

        OPERATOR& operator=(const UNARYFUNCTOR<double,double,double(*)(double)> &UF);
        OPERATOR& operator+=(const UNARYFUNCTOR<double,double,double(*)(double)> &UF);
        OPERATOR& operator-=(const UNARYFUNCTOR<double,double,double(*)(double)> &UF);
        OPERATOR& operator*=(const UNARYFUNCTOR<double,double,double(*)(double)> &UF);
        OPERATOR& operator/=(const UNARYFUNCTOR<double,double,double(*)(double)> &UF);

        OPERATOR& operator+=(const BINARYFUNCTOR<double,double,double,double(*)(double,double)> &BF);
        OPERATOR& operator-=(const BINARYFUNCTOR<double,double,double,double(*)(double,double)> &BF);
        OPERATOR& operator*=(const BINARYFUNCTOR<double,double,double,double(*)(double,double)> &BF);
        OPERATOR& operator/=(const BINARYFUNCTOR<double,double,double,double(*)(double,double)> &BF);

//        template<typename UMFType> OPERATOR& operator+=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF); // C(x) + UMF(x)
//        template<typename UMFType> OPERATOR& operator-=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF);
//        template<typename UMFType> OPERATOR& operator*=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF);
//        template<typename UMFType> OPERATOR& operator/=(const UNARYMEMBERFUNCTION<UMFType,double> &UMF);        

//        template<typename BMFType> OPERATOR& operator+=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF); // C(x) + BMF(x,C(x))
//        template<typename BMFType> OPERATOR& operator-=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF);
//        template<typename BMFType> OPERATOR& operator*=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF);
//        template<typename BMFType> OPERATOR& operator/=(const BINARYMEMBERFUNCTION<BMFType,double> &BMF);  
       };

// ****************************************************************
// ****************************************************************
// ****************************************************************

} // end of namespace interpolation

#endif
