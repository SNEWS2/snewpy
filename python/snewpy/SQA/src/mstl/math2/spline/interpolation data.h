
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
class XDATA_MULTIPLESETS;

template<size_t ND> class YDATA_SINGLESET; // the data is ND dimensional
class YDATA_MULTIPLESETS;

//class DATA_SINGLEXSET_SINGLEYSET;
//class DATA_SINGLEXSET_MULTIPLEYSETS;
//template<size_t ND> class DATA_MULTIPLEXSETS_SINGLEYSET; // the y data is ND dimensional

class DELTAX_SINGLESET;
class DELTAX_MULTIPLESETS;

class DELTAY_SINGLESET; 
class DELTAY_MULTIPLESETS;

class GRADEBASE_SINGLEXSET_SINGLEYSET;
//class LINEGRADE_SINGLEXSET_SINGLEYSET;
//class LOCALGRADE_SINGLEXSET_SINGLEYSET;

class GRADEBASE_SINGLEXSET_MULTIPLEYSETS;
//class LINEGRADE_SINGLEXSET_MULTIPLEYSETS;
//class LOCALGRADE_SINGLEXSET_MULTIPLEYSETS;

class SPLINE_SINGLEXSET_SINGLEYSET;
class SPLINE_SINGLEXSET_MULTIPLEYSETS;
template<size_t ND> class SPLINE_MULTIPLEXSETS_SINGLEYSET; // the y data is ND dimensional

class XLIMITS_SINGLESET;
class XLIMITS_MULTIPLESETS;

class YLIMITS_SINGLESET;
class YLIMITS_MULTIPLESETS;
      
class XYLIMITS_SINGLEXSET_SINGLEYSET;
class XYLIMITS_SINGLEXSET_MULTIPLEYSETS;

class SORT_SINGLEXSET_SINGLEYSET;
class SORT_SINGLEXSET_MULTIPLEYSETS;

class READWRITE_SINGLEXSET_SINGLEYSET;
class READWRITE_SINGLEXSET_MULTIPLEYSETS;
template<size_t ND> class READWRITE_MULTIPLEXSETS_SINGLEYSET; // the y data is ND dimensional

class XINTERVAL_SINGLESET;
class XINTERVAL_MULTIPLESETS;
class YINTERVAL_SINGLESET;
class YINTERVAL_MULTIPLESETS;

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
           
             double _pX(int i) const { return x[i-1];}
             std::vector<double> _pX(void) const { return x;}

             void _pSetX(int i,double X) const { x[i-1]=X;}
             void _pSetX(const std::vector<double> &X) const { x=X;}

             // *******************

             void CreateX(size_t NX) const { x=std::vector<double>(NX);}
             void CreateX(void) const {;}

             void CopyX(const XDATA_SINGLESET &C) const { x=C.X();}

             void DestroyX(void) const { x.clear();}

             // ********************
             // ********************
             // ********************

	     public :

             XDATA_SINGLESET(void) { CreateX();}
             XDATA_SINGLESET(size_t NX) { CreateX(NX);}             
             XDATA_SINGLESET(const std::vector<double> &X) { _pSetX(X);}
             XDATA_SINGLESET(size_t NX,const double *X);
             XDATA_SINGLESET(const XDATA_SINGLESET &C) { CopyX(C);}
             
             virtual ~XDATA_SINGLESET(void){ DestroyX();}

             // *************************************

             bool Empty(void) const { return x.empty();}

             size_t NX(void) const { return x.size();}

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

// *******************
// *******************
// *******************

class XDATA_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;
                        
             mutable std::vector<std::vector<double> > x;  // x is the data

             // *******************
             // *******************
             // *******************

             protected :
             
             double _pX(int i,int j) const { return x[i-1][j-1];}
             std::vector<double> _pX(int i) const { return x[i-1];}
             std::vector<std::vector<double> > _pX(void) const { return x;}

             void _pSetX(int i,int j,double X) const { x[i-1][j-1]=X;}
             void _pSetX(int i,const std::vector<double> &X) const { x[i-1]=X;}
             void _pSetX(const std::vector<std::vector<double> > &X) const { x=X;}

             // *******************

             void CreateX(size_t NXSets) const { x=std::vector<std::vector<double> >(NXSets);}
             void CreateX(size_t NXSets,size_t NX) const { x=std::vector<std::vector<double> >(NXSets,std::vector<double>(NX));}
             void CreateX(size_t NXSets,std::vector<size_t> NX) const;
             void CreateX(void) const {;}

             void CreateXSet(int i,size_t NX) const { x[i-1]=std::vector<double>(NX);}

             void CopyX(const XDATA_MULTIPLESETS &C) const { x=C.X();}

             void DestroyX(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             XDATA_MULTIPLESETS(void) { CreateX();}
             XDATA_MULTIPLESETS(size_t NXSets) { CreateX(NXSets);}             
             XDATA_MULTIPLESETS(size_t NXSets,size_t NX) { CreateX(NXSets,NX);}               // each x set has the same size           
             XDATA_MULTIPLESETS(size_t NXSets,std::vector<size_t> NX) { CreateX(NXSets,NX);}  // each x has the size given by the elements of NX
             XDATA_MULTIPLESETS(const std::vector<std::vector<double> > &X) { _pSetX(X);}
             XDATA_MULTIPLESETS(size_t NXSets,size_t NX,const double **X);
             XDATA_MULTIPLESETS(const XDATA_MULTIPLESETS &C) { CopyX(C);}
             
             virtual ~XDATA_MULTIPLESETS(void){ DestroyX();}

             // *************************************

             bool Empty(void) const { return x.empty();}
             bool Empty(int i) const { return x[i].empty();}

             size_t NXSets(void) const { return x.size();}
             size_t NX(int i) const { return x[i-1].size();}
             std::vector<size_t> NX(void) const; 

             // retrieve or set the data
             double X(int i,int j) const;
             std::vector<double> X(int i) const;
             std::vector<double> X(std::vector<int> i) const;
             std::vector<std::vector<double> > X(void) const { return _pX();}

             void SetX(int i,int j,double X) const;
             void SetX(int i,const std::vector<double> &X) const;
             void SetX(const std::vector<std::vector<double> > &X) const { _pSetX(X);}

             // extend the data
             void AddX(int i,const std::vector<double> &X);
             void AddX(const std::vector<std::vector<double> > &X);

             // transform the data 
             template<typename UNARYFUNCTOR> void TransformX(int i,const UNARYFUNCTOR &UF);
             template<typename UNARYFUNCTOR> void TransformX(const UNARYFUNCTOR &UF);

             void ShiftX(int i,const double&);
             void ShiftX(const double&);
             void RescaleX(int i,const double&);        
             void RescaleX(const double&);

             void AbsX(int i);
             void AbsX(void);

             // operators
             XDATA_MULTIPLESETS& operator=(const XDATA_MULTIPLESETS &C){ CopyX(C); return *this;}
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<size_t ND> class YDATA_SINGLESET
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
            
             double _pY(std::vector<int> i) const;
             TARRAY<double,ND> _pY(void) const { return y;}

             void _pSetY(std::vector<int> i,double Y) const;
             void _pSetY(const TARRAY<double,ND> &Y) const { y=Y;}

             // *******************

             void CreateY(std::vector<size_t> NY) const { y=TARRAY<double,ND>(NY);}     
        
             void CopyY(const YDATA_SINGLESET &C) const { y=C.Y();}

             void DestroyY(void) const {;} 

             // ********************
             // ********************
             // ********************

	     public :

	     // void constructor
             YDATA_SINGLESET(void){;}
             
             // set n DATA
             YDATA_SINGLESET(std::vector<size_t> NY){ CreateY(NY);}             
                          
             // construct from data in containers
             YDATA_SINGLESET(const TARRAY<double,ND> &Y) { _pSetY(Y);}

             YDATA_SINGLESET(const YDATA_SINGLESET &C) { CopyY(C);}
             
             virtual ~YDATA_SINGLESET(void){ DestroyY();}

             // *************************************

             bool Empty(void) const { return y.Empty();}

             std::vector<size_t> NY(void) const { return y.Dimensions();}
             size_t NY(int j) const { return y.Dimension(j);}

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
             
             double _pY(int i) const { return y[i-1];}
             std::vector<double> _pY(void) const { return y;}

             void _pSetY(int i,double Y) const { y[i-1]=Y;}
             void _pSetY(const std::vector<double> &Y) const { y=Y;}

             // *******************

             void CreateY(size_t NY) const { y=std::vector<double>(NY);}

             void CopyY(const YDATA_SINGLESET &C) const { y=C.Y();}

             void DestroyY(void) const { y.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

	     // void constructor
             YDATA_SINGLESET(void){ CreateY(0);}
             
             // set n DATA
             YDATA_SINGLESET(size_t NY){ CreateY(NY);}             
                          
             // construct from data in containers
             YDATA_SINGLESET(const std::vector<double> &Y) { _pSetY(Y);}
             YDATA_SINGLESET(size_t NY,const double *Y);

             YDATA_SINGLESET(const YDATA_SINGLESET &C) { CopyY(C);}
             
             virtual ~YDATA_SINGLESET(void){ DestroyY();}

             // *************************************

             bool Empty(void) const { return y.empty();}

             size_t NY(void) const { return y.size();}

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

class YDATA_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;
                        
             mutable std::vector<std::vector<double> > y;  // y is the data

             // *******************
             // *******************
             // *******************

             protected :

             double _pY(int i,int j) const { return y[i-1][j-1];}
             std::vector<double> _pY(int i) const { return y[i-1];}
             std::vector<std::vector<double> > _pY(void) const { return y;}

             void _pSetY(int i,int j,double Y) const { y[i-1][j-1]=Y;}
             void _pSetY(int i,const std::vector<double> &Y) const { y[i-1]=Y;}
             void _pSetY(const std::vector<std::vector<double> > &Y) const { y=Y;}

             // *******************

             void CreateY(size_t NYSets) const { y=std::vector<std::vector<double> >(NYSets);}
             void CreateY(size_t NYSets,size_t NY) const { y=std::vector<std::vector<double> >(NYSets,std::vector<double>(NY));}
             void CreateY(size_t NYSets,std::vector<size_t> NY) const;
             void CreateY(void) const {;}

             void CreateYSet(int i,size_t NY) const { y[i-1]=std::vector<double>(NY);}

             void CopyY(const YDATA_MULTIPLESETS &C) const { y=C.Y();}

             void DestroyY(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             YDATA_MULTIPLESETS(void){ CreateY();}
             YDATA_MULTIPLESETS(size_t NYSets,size_t NY){ CreateY(NYSets,NY);}             
             YDATA_MULTIPLESETS(size_t NYSets,std::vector<size_t> NY){ CreateY(NYSets,NY);}             
             YDATA_MULTIPLESETS(const std::vector<std::vector<double> > &Y) { _pSetY(Y);}
             YDATA_MULTIPLESETS(size_t NYSets,size_t NY,const double **Y);
             YDATA_MULTIPLESETS(const YDATA_MULTIPLESETS &C) { CopyY(C);}
             
             virtual ~YDATA_MULTIPLESETS(void){ DestroyY();}

             // *************************************

             bool Empty(void) const { return y.empty();}
             bool Empty(int i) const { return y[i].empty();}

             size_t NYSets(void) const { return y.size();}
             size_t NY(int i) const { return y[i-1].size();}
             std::vector<size_t> NY(void) const; 

             // retrieve or set the data
             double Y(int i,int j) const;
             std::vector<double> Y(int i) const;
             std::vector<double> Y(std::vector<int> i) const;
             std::vector<std::vector<double> > Y(void) const { return _pY();}

             void SetY(int i,int j,double Y) const;
             void SetY(int i,const std::vector<double> &Y) const;
             void SetY(const std::vector<std::vector<double> > &Y) const { _pSetY(Y);}

             // extend the data
             void AddY(int i,const std::vector<double> &Y);
             void AddY(const std::vector<std::vector<double> > &Y);

             // transform the data 
             template<typename UNARYFUNCTOR> void TransformY(const UNARYFUNCTOR &UF);
             template<typename UNARYFUNCTOR> void TransformY(int i,const UNARYFUNCTOR &UF);

             void ShiftY(const double&);
             void ShiftY(int i,const double&);
             void RescaleY(const double&);
             void RescaleY(int i,const double&);        

             void AbsY(void);
             void AbsY(int i);

             // operators
             YDATA_MULTIPLESETS& operator=(const YDATA_MULTIPLESETS &C){ CopyY(C); return *this;}
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

//class DATA_SINGLEXSET_SINGLEYSET : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1> {};

//class DATA_SINGLEXSET_MULTIPLEYSETS : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1> {};

//template<size_t ND> class DATA_MULTIPLEXSETS_SINGLEYSET : virtual public XDATA_MULTIPLESETS, virtual public YDATA_SINGLESET<ND> {};

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class DELTAX_SINGLESET 
      : virtual public XDATA_SINGLESET
           { private :

             static std::string CLASS_NAME;

             mutable bool xdifferenced;
             mutable std::vector<double> deltax;  // deltax are the differences between data points

             // *******************
             // *******************
             // *******************

             protected :
             
             double _pDeltaX(int i) const { return deltax[i-1];}
             std::vector<double> _pDeltaX(void) const { return deltax;}

             void _pSetDeltaX(int i,double DELTAX) const { deltax[i-1]=DELTAX;}
             void _pSetDeltaX(const std::vector<double> &DELTAX) const { deltax=DELTAX;}

             void SetXDifferenced(bool XDIFFERENCED) const { xdifferenced=XDIFFERENCED;}

             // *******************

             void CreateDeltaX(size_t NDeltaX) const { xdifferenced=false; deltax=std::vector<double>(NDeltaX);}
             void CreateDeltaX(void) const { xdifferenced=false;}

             void CopyDeltaX(const DELTAX_SINGLESET &C) const { xdifferenced=C.XDifferenced(); deltax=C.DeltaX();}

             void DestroyDeltaX(void) const { xdifferenced=false; deltax.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

             DELTAX_SINGLESET(void) { CreateDeltaX();}
             DELTAX_SINGLESET(size_t NDeltaX) { CreateDeltaX(NDeltaX);}
             DELTAX_SINGLESET(const DELTAX_SINGLESET &C) { CopyDeltaX(C);}

             virtual ~DELTAX_SINGLESET(void) { DestroyDeltaX();}

             size_t NDeltaX(void) const { return deltax.size();}

             double DeltaX(int i) const;
             std::vector<double> DeltaX(void) const { return _pDeltaX();}          

             bool XDifferenced(void) const { return xdifferenced;}
             void XDifference(void) const;
            };

// *************************************
// *************************************
// *************************************

class DELTAX_MULTIPLESETS 
      : virtual public XDATA_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;

             mutable bool xdifferenced;
             mutable std::vector<std::vector<double> > deltax;  // deltax are the differences between data points

             // *******************
             // *******************
             // *******************

             protected :
             
             double _pDeltaX(int i,int j) const { return deltax[i-1][j-1];}
             std::vector<double> _pDeltaX(int i) const { return deltax[i-1];}
             std::vector<std::vector<double> > _pDeltaX(void) const { return deltax;}

             void _pSetDeltaX(int i,int j,double DELTAX) const { deltax[i-1][j-1]=DELTAX;}
             void _pSetDeltaX(int i,const std::vector<double> &DELTAX) const { deltax[i-1]=DELTAX;}
             void _pSetDeltaX(const std::vector<std::vector<double> > &DELTAX) const { deltax=DELTAX;}

             void SetXDifferenced(bool XDIFFERENCED) const { xdifferenced=XDIFFERENCED;}

             // *******************

             void CreateDeltaX(size_t NDeltaXSets) const { xdifferenced=false; deltax=std::vector<std::vector<double> >(NDeltaXSets);}
             void CreateDeltaX(size_t NDeltaXSets,size_t NDeltaX) const { xdifferenced=false; deltax=std::vector<std::vector<double> >(NDeltaXSets,std::vector<double>(NDeltaX));}
             void CreateDeltaX(size_t NDeltaXSets,std::vector<size_t> NDeltaX) const;
             void CreateDeltaX(void) const { xdifferenced=false;}

             void CreateDeltaXSet(int i,size_t NDeltaX) const { xdifferenced=false; deltax[i-1]=std::vector<double>(NDeltaX);}

             void CopyDeltaX(const DELTAX_MULTIPLESETS &C) const { xdifferenced=C.XDifferenced(); deltax=C.DeltaX();}

             void DestroyDeltaX(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             DELTAX_MULTIPLESETS(void) { CreateDeltaX();}
             DELTAX_MULTIPLESETS(size_t NDeltaXSets) { CreateDeltaX(NDeltaXSets);}
             DELTAX_MULTIPLESETS(size_t NDeltaXSets,size_t NDeltaX) { CreateDeltaX(NDeltaXSets,NDeltaX);}
             DELTAX_MULTIPLESETS(size_t NDeltaXSets,std::vector<size_t> NDeltaX) { CreateDeltaX(NDeltaXSets,NDeltaX);}
             DELTAX_MULTIPLESETS(const DELTAX_MULTIPLESETS &C) { CopyDeltaX(C);}
             
             virtual ~DELTAX_MULTIPLESETS(void) { DestroyDeltaX();}

             size_t NDeltaXSets(void) const { return deltax.size();}
             size_t NDeltaX(int i) const { return deltax[i-1].size();}
             std::vector<size_t> NDeltaX(void) const; 

             double DeltaX(int i,int j) const;
             std::vector<double> DeltaX(int i) const;
             std::vector<double> DeltaX(std::vector<int> i) const;
             std::vector<std::vector<double> > DeltaX(void) const { return _pDeltaX();}

             bool XDifferenced(void) const { return xdifferenced;}
             void XDifference(void) const;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class DELTAY_SINGLESET
      : virtual public YDATA_SINGLESET<1>
           { private :

             static std::string CLASS_NAME;

             mutable bool ydifferenced;
             mutable std::vector<double> deltay;  // deltay are the differences between data points

             // *******************
             // *******************
             // *******************

             protected :
             
             double _pDeltaY(int i) const { return deltay[i-1];}
             std::vector<double> _pDeltaY(void) const { return deltay;}

             void _pSetDeltaY(int i,double DELTAY) const { deltay[i-1]=DELTAY;}
             void _pSetDeltaY(const std::vector<double> &DELTAY) const { deltay=DELTAY;}

             void SetYDifferenced(bool YDIFFERENCED) const { ydifferenced=YDIFFERENCED;}

             // *******************

             void CreateDeltaY(size_t NDeltaY) const { ydifferenced=false; deltay=std::vector<double>(NDeltaY);}

             void CopyDeltaY(const DELTAY_SINGLESET &C) const { ydifferenced=C.YDifferenced(); deltay=C.DeltaY();}

             void DestroyDeltaY(void) const { ydifferenced=false; deltay.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

             DELTAY_SINGLESET(void){ CreateDeltaY(0);}
             DELTAY_SINGLESET(size_t NDeltaY){ CreateDeltaY(NDeltaY);}
             DELTAY_SINGLESET(const DELTAY_SINGLESET &C) { CopyDeltaY(C);}
             
             virtual ~DELTAY_SINGLESET(void) { DestroyDeltaY();}

             size_t NDeltaY(void) const { return deltay.size();}

             double DeltaY(int i) const;
             std::vector<double> DeltaY(void) const { return _pDeltaY();}

             bool YDifferenced(void) const { return ydifferenced;}
             void YDifference(void) const;
            };

// *************************************
// *************************************
// *************************************

class DELTAY_MULTIPLESETS 
      : virtual public YDATA_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;

             mutable bool ydifferenced;
             mutable std::vector<std::vector<double> > deltay;  // deltax are the differences between data points

             // *******************
             // *******************
             // *******************

             protected :

             double _pDeltaY(int i,int j) const { return deltay[i-1][j-1];}
             std::vector<double> _pDeltaY(int i) const { return deltay[i-1];}
             std::vector<std::vector<double> > _pDeltaY(void) const { return deltay;}

             void _pSetDeltaY(int i,int j,double DELTAY) const { deltay[i-1][j-1]=DELTAY;}
             void _pSetDeltaY(int i,const std::vector<double> &DELTAY) const { deltay[i-1]=DELTAY;}
             void _pSetDeltaY(const std::vector<std::vector<double> > &DELTAY) const { deltay=DELTAY;}

             void SetYDifferenced(bool YDIFFERENCED) const { ydifferenced=YDIFFERENCED;}

             // *******************

             void CreateDeltaY(size_t NDeltaYSets,size_t NDeltaY) const { ydifferenced=false; deltay=std::vector<std::vector<double> >(NDeltaYSets,std::vector<double>(NDeltaY));}
             void CreateDeltaY(size_t NDeltaYSets) const { ydifferenced=false; deltay=std::vector<std::vector<double> >(NDeltaYSets);}
             void CreateDeltaY(size_t NDeltaYSets,std::vector<size_t> NDeltaY) const;
             void CreateDeltaY(void) const { ydifferenced=false;}

             void CreateDeltaYSet(int i,size_t NDeltaY) const { ydifferenced=false; deltay[i-1]=std::vector<double>(NDeltaY);}

             void CopyDeltaY(const DELTAY_MULTIPLESETS &C) const { ydifferenced=C.YDifferenced(); deltay=C.DeltaY();}

             void DestroyDeltaY(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             DELTAY_MULTIPLESETS(void) { CreateDeltaY();}
             DELTAY_MULTIPLESETS(size_t NDeltaYSets) { CreateDeltaY(NDeltaYSets);}
             DELTAY_MULTIPLESETS(size_t NDeltaYSets,size_t NDeltaY) { CreateDeltaY(NDeltaYSets,NDeltaY);}
             DELTAY_MULTIPLESETS(size_t NDeltaYSets,std::vector<size_t> NDeltaY) { CreateDeltaY(NDeltaYSets,NDeltaY);}
             DELTAY_MULTIPLESETS(const DELTAY_MULTIPLESETS &C) { CopyDeltaY(C);}
             
             virtual ~DELTAY_MULTIPLESETS(void) { DestroyDeltaY();}

             size_t NDeltaYSets(void) const { return deltay.size();}
             size_t NDeltaY(int i) const { return deltay[i-1].size();}
             std::vector<size_t> NDeltaY(void) const; 

             double DeltaY(int i,int j) const;
             std::vector<double> DeltaY(int i) const;
             std::vector<double> DeltaY(std::vector<int> i) const;
             std::vector<std::vector<double> > DeltaY(void) const { return _pDeltaY();}

             bool YDifferenced(void) const { return ydifferenced;}
             void YDifference(void) const;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class GRADEBASE_SINGLEXSET_SINGLEYSET
      : virtual public DELTAX_SINGLESET, 
        virtual public DELTAY_SINGLESET
           { private :

             static std::string CLASS_NAME;
                        
             mutable bool graded;
             mutable std::vector<double> gL,gR;  // gL and gR are the gradients at the left and right side gL[i+1] is the same x point as gR[i] 

             // *******************
             // *******************
             // *******************

             protected :
            
             double _pGL(int i) const { return gL[i-1];}
             double _pGR(int i) const { return gR[i-1];}
             std::vector<double> _pGL(void) const { return gL;}
             std::vector<double> _pGR(void) const { return gR;}

             void _pSetGL(int i,double GL) const { gL[i-1]=GL;}
             void _pSetGL(const std::vector<double> &GL) const { gL=GL;}
             void _pSetGR(int i,double GR) const { gR[i-1]=GR;}
             void _pSetGR(const std::vector<double> &GR) const { gR=GR;}

             void SetGraded(bool GRADED) const { graded=GRADED;}

             // *******************

             void CreateG(size_t NG) const { graded=false; gL=gR=std::vector<double>(NG-1);}
             void CreateG(void) const { graded=false;}

             void CopyG(const GRADEBASE_SINGLEXSET_SINGLEYSET &C) const { graded=C.Graded(); gL=C.GL(); gR=C.GR();}

             void DestroyG(void) const { graded=false; gL.clear(); gR.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

	     // void constructor
             GRADEBASE_SINGLEXSET_SINGLEYSET(void) { CreateG();}
             GRADEBASE_SINGLEXSET_SINGLEYSET(size_t NG) { CreateG(NG);}
             GRADEBASE_SINGLEXSET_SINGLEYSET(const GRADEBASE_SINGLEXSET_SINGLEYSET &C) { CopyG(C);}
             
             virtual ~GRADEBASE_SINGLEXSET_SINGLEYSET(void) { DestroyG();}

             // *************************************

             virtual void Grade(void) const =0;
             void LineGrade(void) const;
             void CatmullRomGrade(void) const { KochanekBartelsGrade(0.,0.,0.);}
             void KochanekBartelsGrade(double t,double b,double c) const;

             bool Graded(void) const { return graded;}

             size_t NGL(void) const { return gL.size();}
             size_t NGR(void) const { return gR.size();}
             size_t NG(void) const { return gL.size()+1;}

             // retrieve or set the data
             double GL(int i) const;
             std::vector<double> GL(void) const;

             double GR(int i) const;
             std::vector<double> GR(void) const;

             void SetGL(int i,double GL) const;
             void SetGL(const std::vector<double> &GL) const;

             void SetGR(int i,double GR) const;
             void SetGR(const std::vector<double> &GR) const;

             double G(int i) const;
             void SetG(int i,double G) const;
            };

// *************************************

/*class LINEGRADE_SINGLEXSET_SINGLEYSET 
      : virtual public GRADEBASE_SINGLEXSET_SINGLEYSET, 
        virtual public DELTAX_SINGLESET, virtual public DELTAY_SINGLESET
      { public :

        LINEGRADE_SINGLEXSET_SINGLEYSET(void) : GRADEBASE_SINGLEXSET_SINGLEYSET() { ;}
        LINEGRADE_SINGLEXSET_SINGLEYSET(size_t NG) : GRADEBASE_SINGLEXSET_SINGLEYSET(NG) {;}
        LINEGRADE_SINGLEXSET_SINGLEYSET(const LINEGRADE_SINGLEXSET_SINGLEYSET &C) : GRADEBASE_SINGLEXSET_SINGLEYSET(C) {;}

        void Grade(void) const;
       };*/

// *************************************

/*class LOCALGRADE_SINGLEXSET_SINGLEYSET 
      : virtual public GRADEBASE_SINGLEXSET_SINGLEYSET, 
        virtual public XDATA_SINGLESET, virtual public DELTAX_SINGLESET, 
        virtual public YDATA_SINGLESET<1>, virtual public DELTAY_SINGLESET
      { public :

        LOCALGRADE_SINGLEXSET_SINGLEYSET(void) : GRADEBASE_SINGLEXSET_SINGLEYSET() { ;}
        LOCALGRADE_SINGLEXSET_SINGLEYSET(size_t NG) : GRADEBASE_SINGLEXSET_SINGLEYSET(NG) {;}
        LOCALGRADE_SINGLEXSET_SINGLEYSET(const LOCALGRADE_SINGLEXSET_SINGLEYSET &C) : GRADEBASE_SINGLEXSET_SINGLEYSET(C) {;}

        virtual ~LOCALGRADE_SINGLEXSET_SINGLEYSET (void) {;}

        void Grade(void) const;
       };*/

// *************************************
// *************************************
// *************************************

class GRADEBASE_SINGLEXSET_MULTIPLEYSETS
      : virtual public DELTAX_SINGLESET, 
        virtual public DELTAY_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;

             mutable bool graded;                        
             mutable std::vector<std::vector<double> > gL,gR;  // g are the gradients

             // *******************
             // *******************
             // *******************

             protected :

             double _pGL(int i,int j) const { return gL[i-1][j-1];}
             std::vector<double> _pGL(int i) const { return gL[i-1];}
             std::vector<std::vector<double> > _pGL(void) const { return gL;}

             double _pGR(int i,int j) const { return gR[i-1][j-1];}
             std::vector<double> _pGR(int i) const { return gR[i-1];}
             std::vector<std::vector<double> > _pGR(void) const { return gR;}

             void _pSetGL(int i,int j,double GL) const { gL[i-1][j-1]=GL;}
             void _pSetGL(int i,const std::vector<double> &GL) const { gL[i-1]=GL;}

             void _pSetGR(int i,int j,double GR) const { gR[i-1][j-1]=GR;}
             void _pSetGR(int i,const std::vector<double> &GR) const { gR[i-1]=GR;}

             void SetGraded(bool GRADED) const { graded=GRADED;}

             // *******************

             void CreateG(size_t NGSets) const { graded=false; gL=gR=std::vector<std::vector<double> >(NGSets);}
             void CreateG(size_t NGSets,size_t NG) const { graded=false; gL=gR=std::vector<std::vector<double> >(NGSets,std::vector<double>(NG-1));}
             void CreateG(void) const { graded=false;}

             void CopyG(const GRADEBASE_SINGLEXSET_MULTIPLEYSETS &C) const { graded=C.Graded(); gL=C.GL(); gR=C.GR();}

             void DestroyG(void) const;

             // ********************
             // ********************
             // ********************

	     public :

             GRADEBASE_SINGLEXSET_MULTIPLEYSETS(void) { CreateG();}
             GRADEBASE_SINGLEXSET_MULTIPLEYSETS(size_t NGSets) { CreateG(NGSets);}
             GRADEBASE_SINGLEXSET_MULTIPLEYSETS(size_t NGSets,size_t NG) { CreateG(NGSets,NG);}
             GRADEBASE_SINGLEXSET_MULTIPLEYSETS(const GRADEBASE_SINGLEXSET_MULTIPLEYSETS &C) { CopyG(C);}
             
             virtual ~GRADEBASE_SINGLEXSET_MULTIPLEYSETS(void) { DestroyG();}

             // *************************************

             virtual void Grade(void) const =0;
             void LineGrade(void) const;
             void CatmullRomGrade(void) const { KochanekBartelsGrade(0.,0.,0.);}
             void KochanekBartelsGrade(double t,double b,double c) const;

             bool Graded(void) const { return graded;}

             size_t NGSets(void) const { return gL.size();}
             size_t NGL(int i) const { return gL[i-1].size();}
             size_t NGR(int i) const { return gR[i-1].size();}
             size_t NG(int i) const { return gL[i-1].size()+1;}

             // retrieve or set the data
             double GL(int i,int j) const;
             std::vector<double> GL(int i) const;
             std::vector<std::vector<double> > GL(void) const;

             double GR(int i,int j) const;
             std::vector<double> GR(int i) const;
             std::vector<std::vector<double> > GR(void) const;

             void SetGL(int i,int j,double GL) const;
             void SetGL(int i,const std::vector<double> &GL) const;

             void SetGR(int i,int j,double GR) const;
             void SetGR(int i,const std::vector<double> &GR) const;

             double G(int i,int j) const;
             void SetG(int i,int j,double G) const;
            };

// *******************************************************************

/*class LINEGRADE_SINGLEXSET_MULTIPLEYSETS 
      : virtual public GRADEBASE_SINGLEXSET_MULTIPLEYSETS, 
        virtual public DELTAX_SINGLESET, virtual public DELTAY_MULTIPLESETS
      { public :

        LINEGRADE_SINGLEXSET_MULTIPLEYSETS(void) : GRADEBASE_SINGLEXSET_MULTIPLEYSETS() { ;}
        LINEGRADE_SINGLEXSET_MULTIPLEYSETS(size_t NGSets,size_t NG) : GRADEBASE_SINGLEXSET_MULTIPLEYSETS(NGSets,NG) {;}
        LINEGRADE_SINGLEXSET_MULTIPLEYSETS(const LINEGRADE_SINGLEXSET_MULTIPLEYSETS &C) : GRADEBASE_SINGLEXSET_MULTIPLEYSETS(C) {;}

        virtual ~LINEGRADE_SINGLEXSET_MULTIPLEYSETS(void) {;}

        void Grade(void) const;
       };*/

/*class LOCALGRADE_SINGLEXSET_MULTIPLEYSETS 
      : virtual public GRADEBASE_SINGLEXSET_MULTIPLEYSETS, 
        virtual public XDATA_SINGLESET, virtual public DELTAX_SINGLESET, 
        virtual public YDATA_MULTIPLESETS, virtual public DELTAY_MULTIPLESETS
      { public :

        LOCALGRADE_SINGLEXSET_MULTIPLEYSETS(void) : GRADEBASE_SINGLEXSET_MULTIPLEYSETS() { ;}
        LOCALGRADE_SINGLEXSET_MULTIPLEYSETS(size_t NGSets,size_t NG) : GRADEBASE_SINGLEXSET_MULTIPLEYSETS(NGSets,NG) {;}
        LOCALGRADE_SINGLEXSET_MULTIPLEYSETS(const LOCALGRADE_SINGLEXSET_MULTIPLEYSETS &C) : GRADEBASE_SINGLEXSET_MULTIPLEYSETS(C) {;}

        virtual ~LOCALGRADE_SINGLEXSET_MULTIPLEYSETS(void) {;}

        void Grade(void) const;
       };*/

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class SPLINE_SINGLEXSET_SINGLEYSET
           { private :

             static std::string CLASS_NAME;

             mutable bool fitted;
             // note that the second index is *not* offset by unity
             mutable std::vector<std::vector<double> > a; // a are the spline parameters

             // *******************
             // *******************
             // *******************

             protected :

             double _pA(int i,int j) const { return a[i-1][j];}
             std::vector<double> _pA(int i) const { return a[i-1];}
             std::vector<std::vector<double> > _pA(void) const { return a;}

             void SetA(int i,int j,double A) const { a[i-1][j]=A;}
             void SetA(int i,std::vector<double> A) const { a[i-1]=A;}
             void SetA(std::vector<std::vector<double> > A) const { a=A;}

             void SetFitted(bool FITTED) const { fitted=FITTED;}

             // *******************

             void CreateA(size_t NA,size_t NParameters) const { fitted=false; a=std::vector<std::vector<double> >(NA,std::vector<double>(NParameters));}
             void CreateA(size_t NA) const { fitted=false; a=std::vector<std::vector<double> >(NA);}
             void CreateA(void) const { fitted=false;}

             void CreateASet(int i,size_t NParameters) const { fitted=false; a[i-1]=std::vector<double>(NParameters);}

             void CopyA(const SPLINE_SINGLEXSET_SINGLEYSET &C) const { fitted=C.Fitted(); a=C.a;}

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
             SPLINE_SINGLEXSET_SINGLEYSET(size_t NA){ CreateA(NA);}                          
             SPLINE_SINGLEXSET_SINGLEYSET(size_t NA,size_t NParameters){ CreateA(NA,NParameters);}
             SPLINE_SINGLEXSET_SINGLEYSET(const std::vector<std::vector<double> > &A) { SetA(A);}
             SPLINE_SINGLEXSET_SINGLEYSET(const SPLINE_SINGLEXSET_SINGLEYSET &C) { CopyA(C);}
             
             virtual ~SPLINE_SINGLEXSET_SINGLEYSET(void) { DestroyA();}

             size_t NA(void) const { return a.size();}
             size_t NParameters(int i) const { return a[i-1].size();}

             double A(int i,int j) const;
             std::vector<double> A(int i) const;
             std::vector<std::vector<double> > A(void) const;

             virtual void Fit(void) const =0;
             bool Fitted(void) const { return fitted;}
            };

// *************************************
// *************************************
// *************************************

class SPLINE_SINGLEXSET_MULTIPLEYSETS
           { private :

             static std::string CLASS_NAME;

             mutable bool fitted;
             mutable std::vector<std::vector<std::vector<double> > > a; // a are the spline parameters, note that the third index is *not* offset by unity i.e. it goes from 0,1,...

             // *******************
             // *******************
             // *******************

             protected :

             double _pA(int i,int j,int k) const { return a[i-1][j-1][k];}
             std::vector<double> _pA(int i,int j) const { return a[i-1][j-1];}
             std::vector<std::vector<double> > _pA(int i) const { return a[i-1];}
             std::vector<std::vector<std::vector<double> > > _pA(void) const { return a;}
             
             void SetA(int i,int j,int k,double A) const { a[i-1][j-1][k]=A;}
             void SetA(int i,int j,std::vector<double> A) const { a[i-1][j-1]=A;}
             void SetA(int i,std::vector<std::vector<double> > A) const { a[i-1]=A;}
             void SetA(std::vector<std::vector<std::vector<double> > > A) const { a=A;}

             void SetFitted(bool FITTED) const { fitted=FITTED;}

             // *******************

             void CreateA(size_t NASets,size_t NA,size_t NParameters) const { fitted=false; a=std::vector<std::vector<std::vector<double> > >(NASets,std::vector<std::vector<double> >(NA,std::vector<double>(NParameters)));}
             void CreateA(size_t NASets,size_t NA) const { fitted=false; a=std::vector<std::vector<std::vector<double> > >(NASets,std::vector<std::vector<double> >(NA));}
             void CreateA(size_t NASets) const { fitted=false; a=std::vector<std::vector<std::vector<double> > >(NASets);}
             void CreateA(void) const { fitted=false;}

             void CreateASet(int i,int j,size_t NParameters) const { fitted=false; a[i-1][j-1]=std::vector<double>(NParameters);}

             void CopyA(const SPLINE_SINGLEXSET_MULTIPLEYSETS &C) const { fitted=C.Fitted(); a=C.A();}

             void DestroyA(void) const;

             // ********************

             double Interpolate(int j,int i,double T) const;
             std::vector<double> Interpolate(int i,double T) const;

             double Derivative(int j,int i,double T) const;
             std::vector<double> Derivative(int i,double T) const;

             double SecondDerivative(int j,int i,double T) const;
             std::vector<double> SecondDerivative(int i,double T) const;

             double Integral(int j,int i,double Tmin,double Tmax) const;
             std::vector<double> Integral(int i,double Tmin,double Tmax) const;

             //std::vector<double> FindRoots(int j,int i,double Y) const;
             //std::vector<std::vector<double> > FindRoots(int i,std::vector<double> Y) const;

             // ********************
             // ********************

	     public :

             SPLINE_SINGLEXSET_MULTIPLEYSETS(void){ CreateA();}             
             SPLINE_SINGLEXSET_MULTIPLEYSETS(size_t NASets){ CreateA(NASets);}             
             SPLINE_SINGLEXSET_MULTIPLEYSETS(size_t NASets,size_t NA){ CreateA(NASets,NA);}             
             SPLINE_SINGLEXSET_MULTIPLEYSETS(size_t NASets,size_t NA,size_t NParameters){ CreateA(NASets,NA,NParameters);}             
             SPLINE_SINGLEXSET_MULTIPLEYSETS(const std::vector<std::vector<std::vector<double> > > &A) { SetA(A);}
             SPLINE_SINGLEXSET_MULTIPLEYSETS(const SPLINE_SINGLEXSET_MULTIPLEYSETS &C) { CopyA(C);}
             
             virtual ~SPLINE_SINGLEXSET_MULTIPLEYSETS(void) { DestroyA();}

             size_t NASets(void) const { return a.size();}
             size_t NA(int i) const { return a[i-1].size();}
             size_t NParameters(int i,int j) const { return a[i-1][j-1].size();}

             double A(int i,int j,int k) const;
             std::vector<double> A(int i,int j) const;
             std::vector<std::vector<double> > A(int i) const;
             std::vector<std::vector<std::vector<double> > > A(void) const;

             bool Fitted(void) const { return fitted;}
             virtual void Fit(void) const =0;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

template<size_t ND> 
class SPLINE_MULTIPLEXSETS_SINGLEYSET // the y data is ND dimensional
           { private :

             static std::string CLASS_NAME;

             mutable bool fitted;
             mutable TARRAY<TARRAY<double,ND>,ND> a; // a are the spline parameters. note that the second set of indicii are *not* offset by unity

             // *******************
             // *******************
             // *******************

             protected :

             double _pA(std::vector<int> i,std::vector<int> j) const;
             TARRAY<double,ND> _pA(std::vector<int> i) const;
             TARRAY<TARRAY<double,ND>,ND> _pA(void) const { return a;}

             void _pSetA(std::vector<int> i,std::vector<int> j,double D) const;
             void _pSetA(std::vector<int> i,TARRAY<double,ND> D) const;

             void SetA(std::vector<int> i,std::vector<int> j,double A) const;
             void SetA(std::vector<int> i,TARRAY<double,ND> A) const;
             void SetA(TARRAY<TARRAY<double,ND>,ND> A) const { a=A; fitted=false;}

             void SetFitted(bool FITTED) const { fitted=FITTED;}

             // *******************

             void CreateA(std::vector<size_t> NA,std::vector<size_t> NP) const { fitted=false; a=TARRAY<TARRAY<double,ND>,ND>(NA,TARRAY<double,ND>(NP));}
             void CreateA(std::vector<size_t> NA) const { fitted=false; a=TARRAY<TARRAY<double,ND>,ND>(NA);}
             void CreateA(void) const { fitted=false;}

             void CreateASet(std::vector<int> i,std::vector<size_t> NP);

             void CopyA(const SPLINE_SINGLEXSET_SINGLEYSET &C) const { fitted=C.Fitted(); a=C.A();}

             void DestroyA(void) const;

             // ********************

             double Interpolate(std::vector<int> i,std::vector<double> T) const;

             // ********************
             // ********************

	     public :

             SPLINE_MULTIPLEXSETS_SINGLEYSET(void){ CreateA();}
             SPLINE_MULTIPLEXSETS_SINGLEYSET(std::vector<size_t> NA){ CreateA(NA);}                          
             SPLINE_MULTIPLEXSETS_SINGLEYSET(std::vector<size_t> NA,std::vector<size_t> NParameters){ CreateA(NA,NParameters);}
             SPLINE_MULTIPLEXSETS_SINGLEYSET(const TARRAY<TARRAY<double,ND>,ND> &A) { SetA(A);}
             SPLINE_MULTIPLEXSETS_SINGLEYSET(const SPLINE_SINGLEXSET_SINGLEYSET &C) { CopyA(C);}
             
             virtual ~SPLINE_MULTIPLEXSETS_SINGLEYSET(void) { DestroyA();}

             virtual void Fit(void) const =0;
             bool Fitted(void) const { return fitted;}

             std::vector<size_t> NA(void) const { return a.Dimensions();}
             size_t NA(int i) const { return a.Dimension(i);}
             std::vector<size_t> NParameters(std::vector<int> i) const;

             double A(std::vector<int> i,std::vector<int> j) const;
             TARRAY<double,ND> A(std::vector<int> i) const;
             TARRAY<TARRAY<double,ND>,ND> A(void) const;
            };

// ***********************************************************************************************
// ***********************************************************************************************
// ***********************************************************************************************

class XLIMITS_SINGLESET
      : virtual public XDATA_SINGLESET
           { private :

             static std::string CLASS_NAME;
             
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

             void SetXLimited(bool XLIMITED) const { xlimited=XLIMITED;}

             // *******************

             void CreateXLimits(void) const { xlimited=false;}

             void CopyXLimits(const XLIMITS_SINGLESET &C) const { xlimited=C.XLimited(); xmin=C.XMin(); xmax=C.XMax();}

             void DestroyXLimits(void) const { xlimited=false;} 

             // ********************
             // ********************
             // ********************

	     public :

             XLIMITS_SINGLESET(void) { CreateXLimits();}
             XLIMITS_SINGLESET(const XLIMITS_SINGLESET &C) { CopyXLimits(C);}
             
             virtual ~XLIMITS_SINGLESET(void) { DestroyXLimits();}

             virtual void XLimit(void) const;
             bool XLimited(void) const { return xlimited;}

             double XMin(void) const;
             double XMax(void) const;
            };

// *************************************
// *************************************
// *************************************

class XLIMITS_MULTIPLESETS
      : virtual public XDATA_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;
             
             mutable std::vector<double> xmin,xmax;
             mutable bool xlimited;

             // ********************
             // ********************
             // ********************

             protected :

             double _pXMin(int i) const { return xmin[i-1];}
             double _pXMax(int i) const { return xmax[i-1];}

             std::vector<double> _pXMin(void) const { return xmin;}
             std::vector<double> _pXMax(void) const { return xmax;}

             void SetXMin(int i,double XMIN) const { xmin[i-1]=XMIN;}
             void SetXMax(int i,double XMAX) const { xmax[i-1]=XMAX;}

             void SetXMin(std::vector<double> XMIN) const { xmin=XMIN;}
             void SetXMax(std::vector<double> XMAX) const { xmax=XMAX;}

             void SetXLimited(bool XLIMITED) const { xlimited=XLIMITED;}

             // *******************

             void CreateXLimits(void) const { xlimited=false; xmin=xmax=std::vector<double>(NXSets());}

             void CopyXLimits(const XLIMITS_MULTIPLESETS &C) const { xlimited=C.XLimited(); xmin=C.XMin(); xmax=C.XMax();}

             void DestroyXLimits(void) const { xlimited=false; xmin.clear(); xmax.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

             XLIMITS_MULTIPLESETS(void){ CreateXLimits();}
             XLIMITS_MULTIPLESETS(const XLIMITS_MULTIPLESETS &C) { CopyXLimits(C);}
             
             virtual ~XLIMITS_MULTIPLESETS(void) { DestroyXLimits();}

             virtual void XLimit(void) const;
             bool XLimited(void) const { return xlimited;}

             double XMin(int i) const;
             double XMax(int i) const;

             std::vector<double> XMin(void) const;
             std::vector<double> XMax(void) const;
            };

// *************************************
// *************************************
// *************************************

class YLIMITS_SINGLESET
      : virtual public YDATA_SINGLESET<1>
           { private :

             static std::string CLASS_NAME;
             
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

             void SetYLimited(bool YLIMITED) const { ylimited=YLIMITED;}

             // *******************

             void CreateYLimits(void) const { ylimited=false;}

             void CopyYLimits(const YLIMITS_SINGLESET &C) const { ylimited=C.YLimited(); ymin=C.YMin(); ymax=C.YMax();}

             void DestroyYLimits(void) const { ylimited=false;} 

             // ********************
             // ********************
             // ********************

	     public :

             YLIMITS_SINGLESET(void) { CreateYLimits();}
             YLIMITS_SINGLESET(const YLIMITS_SINGLESET &C) { CopyYLimits(C);}
             
             virtual ~YLIMITS_SINGLESET(void) { DestroyYLimits();}

             virtual void YLimit(void) const;
             bool YLimited(void) const { return ylimited;}

             double YMin(void) const;
             double YMax(void) const;
            };

// *************************************
// *************************************
// *************************************

class YLIMITS_MULTIPLESETS
      : virtual public YDATA_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;
             
             mutable std::vector<double> ymin,ymax;
             mutable bool ylimited;

             // ********************
             // ********************
             // ********************

             protected :

             double _pYMin(int i) const { return ymin[i-1];}
             double _pYMax(int i) const { return ymax[i-1];}

             std::vector<double> _pYMin(void) const { return ymin;}
             std::vector<double> _pYMax(void) const { return ymax;}

             void SetYMin(int i,double YMIN) const { ymin[i-1]=YMIN;}
             void SetYMax(int i,double YMAX) const { ymax[i-1]=YMAX;}

             void SetYMin(std::vector<double> YMIN) const { ymin=YMIN;}
             void SetYMax(std::vector<double> YMAX) const { ymax=YMAX;}

             void SetYLimited(bool YLIMITED) const { ylimited=YLIMITED;}

             // *******************

             void CreateYLimits(void) const { ylimited=false; ymin=ymax=std::vector<double>(NYSets());}

             void CopyYLimits(const YLIMITS_MULTIPLESETS &C) const { ylimited=C.YLimited(); ymin=C.YMin(); ymax=C.YMax();}

             void DestroyYLimits(void) const { ylimited=false; ymin.clear(); ymax.clear();} 

             // ********************
             // ********************
             // ********************

	     public :

             YLIMITS_MULTIPLESETS(void) { CreateYLimits();}
             YLIMITS_MULTIPLESETS(const YLIMITS_MULTIPLESETS &C) { CopyYLimits(C);}
             
             virtual ~YLIMITS_MULTIPLESETS(void) { DestroyYLimits();}

             virtual void YLimit(void) const;
             bool YLimited(void) const { return ylimited;}

             double YMin(int i) const;
             double YMax(int i) const;

             std::vector<double> YMin(void) const;
             std::vector<double> YMax(void) const;
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

             static std::string CLASS_NAME;
             
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

             void SetXYLimited(bool XYLIMITED) const { xylimited=XYLIMITED;}

             // *******************

             void CreateXYLimits(void) const { xylimited=false;}

             void CopyXYLimits(const XYLIMITS_SINGLEXSET_SINGLEYSET &C) const { xylimited=C.XYLimited(); yatxmin=C.YAtXMin(); yatxmax=C.YAtXMax(); xatymin=C.XAtYMin(); xatymax=C.XAtYMax();}

             void DestroyXYLimits(void) const { xylimited=false;} 

             // ********************
             // ********************
             // ********************

	     public :

             XYLIMITS_SINGLEXSET_SINGLEYSET(void){ CreateXYLimits();}
             XYLIMITS_SINGLEXSET_SINGLEYSET(const XYLIMITS_SINGLEXSET_SINGLEYSET &C) { CopyXYLimits(C);}
             
             virtual ~XYLIMITS_SINGLEXSET_SINGLEYSET(void) { DestroyXYLimits();}

             virtual void XYLimit(void) const;
             bool XYLimited(void) const { return xylimited;}

             double YAtXMin(void) const;
             double YAtXMax(void) const;
             double XAtYMin(void) const;
             double XAtYMax(void) const;
            };

// ***********************************************************************************************

class XYLIMITS_SINGLEXSET_MULTIPLEYSETS
      : virtual public XDATA_SINGLESET,
        virtual public YDATA_MULTIPLESETS,
        virtual public XLIMITS_SINGLESET,
        virtual public YLIMITS_MULTIPLESETS
           { private :

             static std::string CLASS_NAME;
             
             mutable std::vector<double> yatxmin,yatxmax;
             mutable bool xylimited;

             // ********************
             // ********************
             // ********************

             protected :

             double _pYAtXMin(int i) const { return yatxmin[i-1];}
             std::vector<double> _pYAtXMin(void) const { return yatxmin;}
             double _pYAtXMax(int i) const { return yatxmax[i-1];}
             std::vector<double> _pYAtXMax(void) const { return yatxmax;}

             void SetYAtXMin(std::vector<double> YATXMIN) const { yatxmin=YATXMIN;}
             void SetYAtXMax(std::vector<double> YATXMAX) const { yatxmax=YATXMAX;}

             void SetXYLimited(bool XYLIMITED) const { xylimited=XYLIMITED;}

             // *******************

             void CreateXYLimits(void) const { xylimited=false;}
             void CreateXYLimits(size_t NYSets) const { xylimited=false; yatxmin=std::vector<double>(NYSets); yatxmax=std::vector<double>(NYSets);}

             void CopyXYLimits(const XYLIMITS_SINGLEXSET_MULTIPLEYSETS &C) const { xylimited=C.XYLimited(); yatxmin=C.YAtXMin(); yatxmax=C.YAtXMax();}

             void DestroyXYLimits(void) const { yatxmin.clear(); yatxmax.clear(); xylimited=false;} 

             // ********************
             // ********************
             // ********************

	     public :

             XYLIMITS_SINGLEXSET_MULTIPLEYSETS(void){ CreateXYLimits();}
             XYLIMITS_SINGLEXSET_MULTIPLEYSETS(size_t NYSets){ CreateXYLimits(NYSets);}
             XYLIMITS_SINGLEXSET_MULTIPLEYSETS(const XYLIMITS_SINGLEXSET_MULTIPLEYSETS &C) { CopyXYLimits(C);}
             
             virtual ~XYLIMITS_SINGLEXSET_MULTIPLEYSETS(void) { DestroyXYLimits();}

             virtual void XYLimit(void) const;
             bool XYLimited(void) const { return xylimited;}

             double YAtXMin(int i) const;
             std::vector<double> YAtXMin(void) const;

             double YAtXMax(int i) const;
             std::vector<double> YAtXMax(void) const;
            };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class SORT_SINGLEXSET_SINGLEYSET
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>
      { private :
        mutable bool sorted;

        protected :
        void SetSorted(bool SORTED) { sorted=SORTED;}

        public :
        SORT_SINGLEXSET_SINGLEYSET(void) { sorted=false;}
        virtual ~SORT_SINGLEXSET_SINGLEYSET(void) { sorted=false;}

        bool Sorted(void) const { return sorted;}
        void Sort(void);
       };

class SORT_SINGLEXSET_MULTIPLEYSETS
      : virtual public XDATA_SINGLESET, virtual public YDATA_MULTIPLESETS
      { private :
        mutable bool sorted;

        protected :
        void SetSorted(bool SORTED) { sorted=SORTED;}

        public :
        SORT_SINGLEXSET_MULTIPLEYSETS(void) { sorted=false;}
        virtual ~SORT_SINGLEXSET_MULTIPLEYSETS(void) { sorted=false;}

        bool Sorted(void) const { return sorted;}
        void Sort(void);
       };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class READWRITE_SINGLEXSET_SINGLEYSET
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>
      { public :
        READWRITE_SINGLEXSET_SINGLEYSET(void) {;}
        READWRITE_SINGLEXSET_SINGLEYSET(std::string filename,size_t N=0) { Read(filename,N);}
        explicit READWRITE_SINGLEXSET_SINGLEYSET(std::string filename,char ignore) { Read(filename,ignore);}

        virtual ~READWRITE_SINGLEXSET_SINGLEYSET(void) {;}

        void Read(std::string filename,size_t N=0);  // ignores the first N lines
        void Read(std::string filename,char ignore); // all lines beginning with the char 'ignore' will be ignored

        void Write(std::string filename);

        std::ostream& operator<<(std::ostream&);
        std::istream& operator>>(std::istream&);
       };

class READWRITE_SINGLEXSET_MULTIPLEYSETS
      : virtual public XDATA_SINGLESET, virtual public YDATA_MULTIPLESETS
      { public :
        READWRITE_SINGLEXSET_MULTIPLEYSETS(void) {;}
        READWRITE_SINGLEXSET_MULTIPLEYSETS(std::string filename,size_t N=0,bool read_size=false) { Read(filename,N,read_size);}
        explicit READWRITE_SINGLEXSET_MULTIPLEYSETS(std::string filename,char ignore,bool read_size=false) { Read(filename,ignore,read_size);}

        virtual ~READWRITE_SINGLEXSET_MULTIPLEYSETS(void) {;}

        void Read(std::string filename,size_t N=0,bool read_size=false);  // ignores the first N lines, size of date table is not read from file
        void Read(std::string filename,char ignore,bool read_size=false); // all lines beginning with the char 'ignore' will be ignored, size of date table is not read from file
        void Write(std::string filename);

        std::ostream& operator<<(std::ostream&);
        //std::istream& operator>>(std::istream&);
       };

template<size_t ND> 
class READWRITE_MULTIPLEXSETS_SINGLEYSET
      : virtual public XDATA_MULTIPLESETS, virtual public YDATA_SINGLESET<ND>
      { public :
        READWRITE_MULTIPLEXSETS_SINGLEYSET(void) {;}
        READWRITE_MULTIPLEXSETS_SINGLEYSET(std::string filename,size_t N=0,bool read_size=false) { Read(filename,N,read_size);}
        explicit READWRITE_MULTIPLEXSETS_SINGLEYSET(std::string filename,char ignore,bool read_size=false) { Read(filename,ignore,read_size);}

        virtual ~READWRITE_MULTIPLEXSETS_SINGLEYSET(void) {;}

        void Read(std::string filename,size_t N=0,bool read_size=false) {;}  // ignores the first N lines, size of date table is not read from file
        void Read(std::string filename,char ignore,bool read_size=false) {;} // all lines beginning with the char 'ignore' will be ignored, size of date table is not read from file
        void Write(std::string filename) {;}

        std::ostream& operator<<(std::ostream &os) { return os;}
        std::istream& operator>>(std::istream &is) { return is;}
       };

template<>
class READWRITE_MULTIPLEXSETS_SINGLEYSET<2>
      : virtual public XDATA_MULTIPLESETS, virtual public YDATA_SINGLESET<2>
      { public :
        READWRITE_MULTIPLEXSETS_SINGLEYSET(void) {;}
        READWRITE_MULTIPLEXSETS_SINGLEYSET(std::string filename,size_t N=0,bool read_size=false) { Read(filename,N,read_size);}
        explicit READWRITE_MULTIPLEXSETS_SINGLEYSET(std::string filename,char ignore,bool read_size=false) { Read(filename,ignore,read_size);}

        virtual ~READWRITE_MULTIPLEXSETS_SINGLEYSET(void) {;}

        void Read(std::string filename,size_t N=0,bool read_size=false);  // ignores the first N lines, size of date table is not read from file
        void Read(std::string filename,char ignore,bool read_size=false); // all lines beginning with the char 'ignore' will be ignored, size of date table is not read from file
        void Write(std::string filename);

        std::ostream& operator<<(std::ostream &os);
        std::istream& operator>>(std::istream &is);
       };

template<>
class READWRITE_MULTIPLEXSETS_SINGLEYSET<3>
      : virtual public XDATA_MULTIPLESETS, virtual public YDATA_SINGLESET<3>
      { public :
        READWRITE_MULTIPLEXSETS_SINGLEYSET(void) {;}
        READWRITE_MULTIPLEXSETS_SINGLEYSET(std::string filename,size_t N=0,bool read_size=false) { Read(filename,N,read_size);}
        explicit READWRITE_MULTIPLEXSETS_SINGLEYSET(std::string filename,char ignore,bool read_size=false) { Read(filename,ignore,read_size);}

        virtual ~READWRITE_MULTIPLEXSETS_SINGLEYSET(void) {;}

        void Read(std::string filename,size_t N=0,bool read_size=false);  // ignores the first N lines
        void Read(std::string filename,char ignore,bool read_size=false); // all lines beginning with the char 'ignore' will be ignored
        void Write(std::string filename);

        std::ostream& operator<<(std::ostream &os);
        std::istream& operator>>(std::istream &is);
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

class XINTERVAL_MULTIPLESETS
      : virtual public XDATA_MULTIPLESETS,
        virtual public XLIMITS_MULTIPLESETS
      { private :
        mutable std::vector<int> interval; // this vector stores previous returns from calls to XInterval so as to make the routine faster

        protected :
        void CreateXInterval(void) const { interval=std::vector<int>(NXSets(),1);}

        public :
        XINTERVAL_MULTIPLESETS(void) { interval=std::vector<int>(NXSets(),1);}
        virtual ~XINTERVAL_MULTIPLESETS(void) { interval=std::vector<int>(NXSets(),1);}

        std::vector<int> XInterval(std::vector<double> X) const;
        int XInterval(int,double X) const;
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

class YINTERVAL_MULTIPLESETS
      : virtual public YDATA_MULTIPLESETS,
        virtual public YLIMITS_MULTIPLESETS
      { private :
        mutable std::vector<int> interval; // this vector stores previous returns from calls to YInterval so as to make the routine faster

        protected :
        void CreateYInterval(void) const { interval=std::vector<int>(NYSets(),1);}

        public :
        YINTERVAL_MULTIPLESETS(void) { interval=std::vector<int>(NYSets(),1);}
        virtual ~YINTERVAL_MULTIPLESETS(void) { interval=std::vector<int>(NYSets(),1);}

        std::vector<int> YInterval(std::vector<double> Y) const;
        int YInterval(int,double Y) const;
       };

// *******************************************************************
// *******************************************************************
// *******************************************************************

class DISTINGUISH_SINGLEYSET
      : virtual public XDATA_SINGLESET, virtual public YDATA_SINGLESET<1>, virtual public SORT_SINGLEXSET_SINGLEYSET
      { private :
        mutable bool distinguished;

        protected :
        void SetDistinguished(bool DISTINGUISHED) const { distinguished=DISTINGUISHED;}

        public :
        DISTINGUISH_SINGLEYSET(void) { distinguished=false;}
        virtual ~DISTINGUISH_SINGLEYSET(void) { distinguished=false;}

        virtual size_t N(void){ return NX();}

        bool Distinguished(void) const { return distinguished;}
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
