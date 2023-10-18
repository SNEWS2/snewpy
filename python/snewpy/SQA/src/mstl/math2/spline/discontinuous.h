
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
class DISCONTINUOUS_MULTIPLESETS;

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
        virtual public READWRITE_SINGLEXSET_SINGLEYSET
      { private :

        mutable bool founddomains;
        mutable std::vector<std::pair<int,int> > domains;

        protected :

        bool FoundDomains(void) const { return founddomains;}
        void SetFoundDomains(bool FOUNDDOMAINS) const { founddomains=FOUNDDOMAINS;}

        public : 
        explicit DISCONTINUOUS(void) { SetFoundDomains(false);}
        explicit DISCONTINUOUS(size_t N) : XDATA_SINGLESET(N), YDATA_SINGLESET<1>(N), DELTAX_SINGLESET(N-1), DELTAY_SINGLESET(N-1), GRADEBASE_SINGLEXSET_SINGLEYSET(N), SPLINE_SINGLEXSET_SINGLEYSET(N-1) { SetFoundDomains(false);}
        explicit DISCONTINUOUS(std::string filename,size_t Nignore=0){ Read(filename,Nignore);}
        explicit DISCONTINUOUS(std::string filename,char ignore){ Read(filename,ignore);}
        explicit DISCONTINUOUS(std::pair<std::vector<double>,std::vector<double> > XY) : XDATA_SINGLESET(XY.first), YDATA_SINGLESET<1>(XY.second) { SetFoundDomains(false);}
        explicit DISCONTINUOUS(std::vector<double> X,std::vector<double> Y) : XDATA_SINGLESET(X), YDATA_SINGLESET<1>(Y) { SetFoundDomains(false);} 
        explicit DISCONTINUOUS(size_t N,const double* X,const double* Y) : XDATA_SINGLESET(N,X), YDATA_SINGLESET<1>(N,Y) { SetFoundDomains(false);}

        ~DISCONTINUOUS(void) { SetFoundDomains(false); domains.clear();}

        // *********************
        
        size_t N(void) const { return NX();}

        double X(int i) const { return XDATA_SINGLESET::X(i);}
        double Y(int i) const { return YDATA_SINGLESET<1>::Y(i);}

        std::vector<double> X(void) const { return XDATA_SINGLESET::X();}
        std::vector<double> Y(void) const { return YDATA_SINGLESET<1>::Y();}

        double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
        double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

        double YMin(void) const { return YLIMITS_SINGLESET::YMin();}
        double YMax(void) const { return YLIMITS_SINGLESET::YMax();}

        double XAtYMin(void) const { return XYLIMITS_SINGLEXSET_SINGLEYSET::XAtYMin();}
        double XAtYMax(void) const { return XYLIMITS_SINGLEXSET_SINGLEYSET::XAtYMax();}

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

        void Read(std::string filename,size_t Nignore=0){ READWRITE_SINGLEXSET_SINGLEYSET::Read(filename,Nignore); SetXDifferenced(false); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetYLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void Read(std::string filename,char ignore){ READWRITE_SINGLEXSET_SINGLEYSET::Read(filename,ignore); SetXDifferenced(false); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetYLimited(false); SetXYLimited(false); SetFoundDomains(false);}

        // these are deprecated in favour of Read
        void Open(std::string filename,size_t Nignore=0){ Read(filename,Nignore);}
        void Open(std::string filename,char ignore){ Read(filename,ignore);}

        void Write(std::string filename){ READWRITE_SINGLEXSET_SINGLEYSET::Write(filename);}

        // *********************

        void FindDomains(void) const;

        size_t NDomains(void) const { if(FoundDomains()==false){ FindDomains();} return domains.size();}
        std::pair<double,double> Domain(int i) const { if(FoundDomains()==false){ FindDomains();} return std::pair<double,double>(_pX(domains[i-1].first),_pX(domains[i-1].second));}

        size_t NDiscontinuities(void) const { if(FoundDomains()==false){ FindDomains();} return domains.size()-1;}
        double Discontinuity(int i) const { if(FoundDomains()==false){ FindDomains();} return _pX(domains[i-1].second);}

        void Grade(void) const; 
        //void AkimaGrade(void) const; 
        //void CatmullRomGrade(void) const; 

        double Interpolate(double X) const;
        double Derivative(double X) const;

        double operator()(double X) const { return Interpolate(X);}
 
        void Fit(void) const;
        void LineFit(void) const;

        DISCONTINUOUS Derivative(void) const;
       };

// *********************************************************************************
// *********************************************************************************
// *********************************************************************************

class DISCONTINUOUS_MULTIPLESETS
      : virtual public XDATA_SINGLESET, 
        virtual public YDATA_MULTIPLESETS,
        virtual public XLIMITS_SINGLESET, 
        virtual public YLIMITS_MULTIPLESETS, 
        virtual public XYLIMITS_SINGLEXSET_MULTIPLEYSETS,
        virtual public XINTERVAL_SINGLESET, 
        virtual protected DELTAX_SINGLESET, 
        virtual protected DELTAY_MULTIPLESETS,
        virtual protected GRADEBASE_SINGLEXSET_MULTIPLEYSETS,
        virtual public SPLINE_SINGLEXSET_MULTIPLEYSETS,
        virtual public SORT_SINGLEXSET_MULTIPLEYSETS, 
        virtual public READWRITE_SINGLEXSET_MULTIPLEYSETS
      { private :

        mutable bool founddomains;
        mutable std::vector<std::pair<int,int> > domains;

        protected :

        bool FoundDomains(void) const { return founddomains;}
        void SetFoundDomains(bool FOUNDDOMAINS) const { founddomains=FOUNDDOMAINS;}

        public : 
        explicit DISCONTINUOUS_MULTIPLESETS(void) { SetFoundDomains(false);}
        explicit DISCONTINUOUS_MULTIPLESETS(size_t NSets,size_t N) : XDATA_SINGLESET(N), YDATA_MULTIPLESETS(NSets,N), XYLIMITS_SINGLEXSET_MULTIPLEYSETS(NSets), DELTAX_SINGLESET(N-1), DELTAY_MULTIPLESETS(NSets,N-1), GRADEBASE_SINGLEXSET_MULTIPLEYSETS(NSets,N), SPLINE_SINGLEXSET_MULTIPLEYSETS(NSets,N-1) 
                        { SetFoundDomains(false);}
        explicit DISCONTINUOUS_MULTIPLESETS(std::string filename){ Read(filename); SetFoundDomains(false);}
        explicit DISCONTINUOUS_MULTIPLESETS(std::string filename,char ignore){ Read(filename,ignore); SetFoundDomains(false);}
        explicit DISCONTINUOUS_MULTIPLESETS(std::vector<std::vector<double> > XY);
        explicit DISCONTINUOUS_MULTIPLESETS(std::vector<double> X, std::vector<std::vector<double> > Y);

        ~DISCONTINUOUS_MULTIPLESETS(void) { SetFoundDomains(false); domains.clear();}

        // *********************
        
        size_t N(void) const { return NX();}
        size_t NSets(void) const { return NYSets();}

        double X(int i) const { return XDATA_SINGLESET::X(i);}
        double Y(int i,int j) const { return YDATA_MULTIPLESETS::Y(i,j);}

        std::vector<double> X(void) const { return XDATA_SINGLESET::X();}
        std::vector<double> Y(int i) const { return YDATA_MULTIPLESETS::Y(i);}

        double XMin(void) const { return XLIMITS_SINGLESET::XMin();}
        double XMax(void) const { return XLIMITS_SINGLESET::XMax();}

        int XInterval(double X) const { return XINTERVAL_SINGLESET::XInterval(X);}

        void SetX(int i,double X){ XDATA_SINGLESET::SetX(i,X); SetXDifferenced(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetGraded(false); SetFitted(false); SetFoundDomains(false);}
        void SetX(const std::vector<double> &X){ XDATA_SINGLESET::SetX(X); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}

        void SetY(int i,int j,double Y){ YDATA_MULTIPLESETS::SetY(i,j,Y); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}
        void SetY(int i,const std::vector<double> &Y){ YDATA_MULTIPLESETS::SetY(i,Y); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        template<typename UNARYFUNCTION> void TransformX(const UNARYFUNCTION &UF){ XDATA_SINGLESET::TransformX(UF); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        template<typename UNARYFUNCTION> void TransformY(const UNARYFUNCTION &UF){ YDATA_MULTIPLESETS::TransformY(UF); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void ShiftX(const double &D){ XDATA_SINGLESET::ShiftX(D); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void ShiftY(const double &D){ YDATA_MULTIPLESETS::ShiftY(D); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void RescaleX(const double &D){ XDATA_SINGLESET::RescaleX(D); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void RescaleY(const double &D){ YDATA_MULTIPLESETS::RescaleY(D); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void AbsX(void){ XDATA_SINGLESET::AbsX(); SetXDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void AbsY(void){ YDATA_MULTIPLESETS::AbsY(); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetYLimited(false); SetXYLimited(false);}

        void Read(std::string filename){ READWRITE_SINGLEXSET_MULTIPLEYSETS::Read(filename); SetXDifferenced(false); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetYLimited(false); SetXYLimited(false); SetFoundDomains(false);}
        void Read(std::string filename,char ignore){ READWRITE_SINGLEXSET_MULTIPLEYSETS::Read(filename,ignore); SetXDifferenced(false); SetYDifferenced(false); SetGraded(false); SetFitted(false); SetSorted(false); SetXLimited(false); SetYLimited(false); SetXYLimited(false); SetFoundDomains(false);}

        // *********************

        void FindDomains(void) const;

        size_t NDomains(void) const { if(FoundDomains()==false){ FindDomains();} return domains.size();}
        std::pair<double,double> Domain(int i) const { if(FoundDomains()==false){ FindDomains();} return std::pair<double,double>(_pX(domains[i-1].first),_pX(domains[i-1].second));}

        size_t NDiscontinuities(void) const { if(FoundDomains()==false){ FindDomains();} return domains.size()-1;}
        double Discontinuity(int i) const { if(FoundDomains()==false){ FindDomains();} return _pX(domains[i-1].second);}

        void Grade(void) const; 
        //void AkimaGrade(void) const; 
        //void CatmullRomGrade(void) const; 

        double Interpolate(int j,double X) const;
        std::vector<double> Interpolate(double X) const;

        double Derivative(int j,double X) const;
        std::vector<double> Derivative(double X) const;

        double operator()(int j,double X) const { return Interpolate(j,X);}
        std::vector<double> operator()(double X) const { return Interpolate(X);}
 
        void Fit(void) const;
        void LineFit(void) const;
       };

} //end of namespace interpolation

#endif
