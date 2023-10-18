#if !defined(_DISCONTINUOUS)
#include "discontinuous.h"
#endif

#if !defined(_DISCONTINUOUS_FUNCTIONS)
#define _DISCONTINUOUS_FUNCTIONS

namespace interpolation{
/*
template <class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXModifyY(DISCONTINUOUS const &C,const double &A,Operator O)
     { std::vector<double> xpoints, ypoints;  
       for(int i=1;i<=static_cast<int>(C.N());i++){ xpoints.push_back(C.X(i)); ypoints.push_back(O(C.Y(i),A));}
       return std::pair<std::vector<double>,std::vector<double> >(xpoints,ypoints);
      }

template <class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXModifyY(const double &A,DISCONTINUOUS const &C,Operator O)
     { std::vector<double> xpoints, ypoints;  
       for(int i=1;i<=static_cast<int>(C.N());i++){ xpoints.push_back(C.X(i)); ypoints.push_back(O(A,C.Y(i)));}
       return std::pair<std::vector<double>,std::vector<double> >(xpoints,ypoints);
      }

// **************************

template <typename IType,class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXCombineY(DISCONTINUOUS const &C1,EXPRESSION<IType> const &C2,Operator O)
     { double xmin=std::max(C1.XMin(),C2.XMin()), xmax=std::min(C1.XMax(),C2.XMax()); 
       if(xmin>xmax){ throw NO_OVERLAP("OverlapXCombineY");}
       int a1min=C1.XInterval(xmin), a2min=C2.XInterval(xmin);

       std::vector<double> xpoints(1,xmin);
       std::vector<double> ypoints(1,O(C1(xmin),C2(xmin)));  

       int a1=a1min+1, a2=a2min+1;
       double x;

       while(C1.X(a1)<xmax && C2.X(a2)<xmax)
            { if(C1.X(a1)==C2.X(a2)){ x=C1.X(a1); xpoints.push_back(x); ypoints.push_back(O(C1.Y(a1),C2(x))); a1++; a2++;}
              else{ if(C1.X(a1)<C2.X(a2)){ x=C1.X(a1); xpoints.push_back(x); ypoints.push_back(O(C1.Y(a1),C2(x))); a1++;}
                    else{ if(C1.X(a1)>C2.X(a2)){ x=C2.X(a2); xpoints.push_back(x); ypoints.push_back(O(C1(x),C2(x))); a2++;} }
                   }
             };
       xpoints.push_back(xmax); ypoints.push_back(O(C1(xmax),C2(xmax)));

       return std::pair<std::vector<double>,std::vector<double> >(xpoints,ypoints);
      }

// **************************

template <typename IType,class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXCombineY(EXPRESSION<IType> const &C1,DISCONTINUOUS const &C2,Operator O)
     { double xmin=std::max(C1.XMin(),C2.XMin()), xmax=std::min(C1.XMax(),C2.XMax()); 
       if(xmin>xmax){ throw NO_OVERLAP("OverlapXCombineY");}
       int a1min=C1.XInterval(xmin), a2min=C2.XInterval(xmin);

       std::vector<double> xpoints(1,xmin);
       std::vector<double> ypoints(1,O(C1(xmin),C2(xmin)));  

       int a1=a1min+1, a2=a2min+1;
       double x;

       while(C1.X(a1)<xmax && C2.X(a2)<xmax)
            { if(C1.X(a1)==C2.X(a2)){ x=C1.X(a1); xpoints.push_back(x); ypoints.push_back(O(C1(x),C2.Y(a2))); a1++; a2++;}
              else{ if(C1.X(a1)<C2.X(a2)){ x=C1.X(a1); xpoints.push_back(x); ypoints.push_back(O(C1(x),C2(x))); a1++;}
                    else{ if(C1.X(a1)>C2.X(a2)){ x=C2.X(a2); xpoints.push_back(x); ypoints.push_back(O(C1(x),C2.Y(a2))); a2++;} }
                   }
             };
       xpoints.push_back(xmax); ypoints.push_back(O(C1(xmax),C2(xmax)));

       return std::pair<std::vector<double>,std::vector<double> >(xpoints,ypoints);
      }
// **************************

template <class Operator> 
std::pair<std::vector<double>,std::vector<double> > OverlapXCombineY(DISCONTINUOUS const &C1,DISCONTINUOUS const &C2,Operator O)
     { double xmin=std::max(C1.XMin(),C2.XMin()), xmax=std::min(C1.XMax(),C2.XMax()); 
       if(xmin>xmax){ throw NO_OVERLAP("OverlapXCombineY");}
       int a1min=C1.XInterval(xmin), a2min=C2.XInterval(xmin);

       std::vector<double> xpoints(1,xmin);
       std::vector<double> ypoints(1,O(C1(xmin),C2(xmin)));  

       int a1=a1min+1, a2=a2min+1;
       double x;

       while(C1.X(a1)<xmax && C2.X(a2)<xmax)
                  { if(C1.X(a1)==C2.X(a2)){ x=C1.X(a1); xpoints.push_back(x); ypoints.push_back(O(C1.Y(a1),C2.Y(a2))); a1++; a2++;}
                     else{ if(C1.X(a1)<C2.X(a2)){ x=C1.X(a1); xpoints.push_back(x); ypoints.push_back(O(C1.Y(a1),C2(x))); a1++;}
                                 else{ if(C1.X(a1)>C2.X(a2)){ x=C2.X(a2); xpoints.push_back(x); ypoints.push_back(O(C1(x),C2.Y(a2))); a2++;} }
                                            }
                    };
       xpoints.push_back(xmax); ypoints.push_back(O(C1(xmax),C2(xmax)));

       return std::pair<std::vector<double>,std::vector<double> >(xpoints,ypoints);
      }

// **************************
// **************************
// **************************

template <typename IType> DISCONTINUOUS operator+(DISCONTINUOUS const &D,EXPRESSION<IType> const &E)
         { return DISCONTINUOUS(OverlapXCombineY(D,E,std::plus<double>()));}

template <typename IType> DISCONTINUOUS operator+(EXPRESSION<IType> const &E,DISCONTINUOUS const &D)  
         { return DISCONTINUOUS(OverlapXCombineY(E,D,std::plus<double>()));}

// **************************

template <typename IType> DISCONTINUOUS operator-(DISCONTINUOUS const &D,EXPRESSION<IType> const &E) 
         { return DISCONTINUOUS(OverlapXCombineY(D,E,std::minus<double>()));}
 
template <typename IType> DISCONTINUOUS operator-(EXPRESSION<IType> const &E,DISCONTINUOUS const &D)  
         { return DISCONTINUOUS(OverlapXCombineY(E,D,std::minus<double>()));}

// **************************

template <typename IType> DISCONTINUOUS operator*(DISCONTINUOUS const &D,EXPRESSION<IType> const &E) 
         { return DISCONTINUOUS(OverlapXCombineY(D,E,std::multiplies<double>()));}

template <typename IType> DISCONTINUOUS operator*(EXPRESSION<IType> const &E,DISCONTINUOUS const &D) 
         { return DISCONTINUOUS(OverlapXCombineY(E,D,std::multiplies<double>()));}

// **************************

template <typename IType> DISCONTINUOUS operator/(DISCONTINUOUS const &D,EXPRESSION<IType> const &E)
         { return DISCONTINUOUS(OverlapXCombineY(D,E,std::divides<double>()));}

template <typename IType> DISCONTINUOUS operator/(EXPRESSION<IType> const &E,DISCONTINUOUS const &D) 
         { return DISCONTINUOUS(OverlapXCombineY(E,D,std::divides<double>()));}
*/
} //end of namespace interpolation

#endif


