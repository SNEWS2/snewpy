
#include <utility>
#include <vector>
#include <algorithm>
#include <functional>

#include "mstl.h"

// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************

#if !defined(_DISCONTINUOUS)
#define _DISCONTINUOUS

namespace interpolation{

class DISCONTINUOUS;

// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************

template <class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXModifyY(DISCONTINUOUS const &C,const double &A,Operator O);
template <class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXModifyY(const double &A,DISCONTINUOUS const &C,Operator O);

// *************************************

DISCONTINUOUS operator+(const DISCONTINUOUS&,const double&);
DISCONTINUOUS operator-(const DISCONTINUOUS&,const double&);
DISCONTINUOUS operator*(const DISCONTINUOUS&,const double&);
DISCONTINUOUS operator/(const DISCONTINUOUS&,const double&);

DISCONTINUOUS operator+(const double&,const DISCONTINUOUS&);
DISCONTINUOUS operator-(const double&,const DISCONTINUOUS&);
DISCONTINUOUS operator*(const double&,const DISCONTINUOUS&);

// *************************************

template <typename IType,class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXCombineY(DISCONTINUOUS const &C1,EXPRESSION<IType> const &C2,Operator O);
template <typename IType,class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXCombineY(EXPRESSION<IType> const &C1,DISCONTINUOUS const &C2,Operator O);
template <class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXCombineY(DISCONTINUOUS const &C1,DISCONTINUOUS const &C2,Operator O);

template <typename IType> DISCONTINUOUS operator+(DISCONTINUOUS const&,EXPRESSION<IType> const&);  // C1(x) + C2(x)
template <typename IType> DISCONTINUOUS operator+(EXPRESSION<IType> const&,DISCONTINUOUS const&);  // C1(x) + C2(x)
DISCONTINUOUS operator+(DISCONTINUOUS const&,DISCONTINUOUS const&);      // C1(x) + C2(x)

template <typename IType> DISCONTINUOUS operator-(DISCONTINUOUS const&,EXPRESSION<IType> const&);  // C1(x) - C2(x)
template <typename IType> DISCONTINUOUS operator-(EXPRESSION<IType> const&,DISCONTINUOUS const&);  // C1(x) - C2(x)
DISCONTINUOUS operator-(DISCONTINUOUS const&,DISCONTINUOUS const&);      // C1(x) - C2(x)

template <typename IType> DISCONTINUOUS operator*(DISCONTINUOUS const&,EXPRESSION<IType> const&);  // C1(x) * C2(x)
template <typename IType> DISCONTINUOUS operator*(EXPRESSION<IType> const&,DISCONTINUOUS const&);  // C1(x) * C2(x)
DISCONTINUOUS operator*(DISCONTINUOUS const&,DISCONTINUOUS const&);      // C1(x) * C2(x)

template <typename IType> DISCONTINUOUS operator/(DISCONTINUOUS const&,EXPRESSION<IType> const&);  // C1(x) / C2(x)
template <typename IType> DISCONTINUOUS operator/(EXPRESSION<IType> const&,DISCONTINUOUS const&);  // C1(x) / C2(x)
DISCONTINUOUS operator/(DISCONTINUOUS const&,DISCONTINUOUS const&);      // C1(x) / C2(x)

// *********************************************************************************
// *********************************************************************************
// *********************************************************************************

class DISCONTINUOUS 
      : virtual public XDATA_SINGLESET, 
        virtual public YDATA_SINGLESET<1>,
        virtual public XLIMITS_SINGLESET, 
        virtual public YLIMITS_SINGLESET, 
        virtual public XYLIMITS_SINGLEXSET_SINGLEYSET,
        virtual public XINTERVAL_SINGLESET, 
        virtual protected DELTAX_SINGLESET, 
        virtual protected DELTAY_SINGLESET,
        virtual protected GRADEBASE_SINGLEXSET_SINGLEYSET,
        virtual public SPLINE_SINGLEXSET_SINGLEYSET,
        virtual public SORT_SINGLEXSET_SINGLEYSET, 
        virtual public OPENWRITE_SINGLEXSET_SINGLEYSET,
        public EXPRESSION<DISCONTINUOUS>,
        virtual public OPERATOR
      { private :

        mutable bool founddomains;
        mutable std::vector<std::pair<int,int> > domains;

        protected :

        bool FoundDomains(void) const { return founddomains;}
        void SetFoundDomains(bool FOUNDDOMAINS) const { founddomains=FOUNDDOMAINS;}

        public : 
        explicit DISCONTINUOUS(void) { SetFoundDomains(false);}
        explicit DISCONTINUOUS(int N) : XDATA_SINGLESET(N), YDATA_SINGLESET<1>(N), DELTAX_SINGLESET(N-1), DELTAY_SINGLESET(N-1), GRADEBASE_SINGLEXSET_SINGLEYSET(N), SPLINE_SINGLEXSET_SINGLEYSET(N-1) { SetFoundDomains(false);}
        explicit DISCONTINUOUS(std::string filename){ Open(filename); SetFoundDomains(false);}
        explicit DISCONTINUOUS(std::string filename,char ignore){ Open(filename,ignore); SetFoundDomains(false);}
        explicit DISCONTINUOUS(std::pair<std::vector<double>,std::vector<double> > XY) : XDATA_SINGLESET(XY.first), YDATA_SINGLESET<1>(XY.second) { SetFoundDomains(false);}
        explicit DISCONTINUOUS(std::vector<double> X,std::vector<double> Y) : XDATA_SINGLESET(X), YDATA_SINGLESET<1>(Y) { SetFoundDomains(false);} 
        explicit DISCONTINUOUS(int N,const double* X,const double* Y) : XDATA_SINGLESET(N,X), YDATA_SINGLESET<1>(N,Y) { SetFoundDomains(false);}
        template <typename IType> DISCONTINUOUS(EXPRESSION<IType> const&);

        ~DISCONTINUOUS(void) { SetFoundDomains(false); domains.clear();}

        // *********************
        
        int N(void) const { return NX();}

        double X(int i) const { return XDATA_SINGLESET::X(i);}
        double Y(int i) const { return YDATA_SINGLESET<1>::Y(i);}

        std::vector<double> X(void) const { return XDATA_SINGLESET::X();}
        std::vector<double> Y(void) const { return YDATA_SINGLESET<1>::Y();}

        double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
        double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

        int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}

        void SetX(int i,double X){ XDATA_SINGLESET::SetX(i,X); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void SetX(const std::vector<double> &X){ XDATA_SINGLESET::SetX(X); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}

        void SetY(int i,double Y){ YDATA_SINGLESET<1>::SetY(i,Y); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false); }
        void SetY(const std::vector<double> &Y){ YDATA_SINGLESET<1>::SetY(Y); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        template<typename UNARYFUNCTION> void TransformX(const UNARYFUNCTION &UF){ XDATA_SINGLESET::TransformX(UF); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        template<typename UNARYFUNCTION> void TransformY(const UNARYFUNCTION &UF){ YDATA_SINGLESET<1>::TransformY(UF); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void ShiftX(const double &D){ XDATA_SINGLESET::ShiftX(D); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void ShiftY(const double &D){ YDATA_SINGLESET<1>::ShiftY(D); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void RescaleX(const double &D){ XDATA_SINGLESET::RescaleX(D); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void RescaleY(const double &D){ YDATA_SINGLESET<1>::RescaleY(D); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void AbsX(void){ XDATA_SINGLESET::AbsX(); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void AbsY(void){ YDATA_SINGLESET<1>::AbsY(); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void Open(std::string filename){ OPENWRITE_SINGLEXSET_SINGLEYSET::Open(filename); SetXDifferenced(false); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetYLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void Open(std::string filename,char ignore){ OPENWRITE_SINGLEXSET_SINGLEYSET::Open(filename,ignore); SetXDifferenced(false); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetYLimited(false); SetXYLimited(false); SetFoundDomains(false);}

        // *********************

        void FindDomains(void) const;

        int NDomains(void) const { if(FoundDomains()==false){ FindDomains();} return domains.size();}
        std::pair<double,double> Domain(int i) const { if(FoundDomains()==false){ FindDomains();} return std::pair<double,double>(_pX(domains[i-1].first),_pX(domains[i-1].second));}

        int NDiscontinuities(void) const { if(FoundDomains()==false){ FindDomains();} return domains.size()-1;}
        double Discontinuity(int i) const { if(FoundDomains()==false){ FindDomains();} return _pX(domains[i-1].second);}

        void Grade(void) const; 

        double Interpolate(double X) const;
        double Derivative(double X) const;

        double operator()(double X) const { return Interpolate(X);}
 
        void Fit(void) const;
        void LineFit(void) const;

        template <typename IType> DISCONTINUOUS& operator=(const EXPRESSION<IType> &I);
       };

} //end of namespace interpolation

#endif
