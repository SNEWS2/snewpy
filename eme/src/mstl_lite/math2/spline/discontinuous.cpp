#include "discontinuous.h"

namespace interpolation{

using std::vector;
using std::pair;

void DISCONTINUOUS::FindDomains(void) const
     { domains.clear();
       int domainlowerlimit=1, domainupperlimit;
       for(int i=1;i<=NX()-1;i++)
          { if(Equality(_pX(i),_pX(i+1))==true)
              { domainupperlimit=i; 
                domains.push_back(pair<int,int>(domainlowerlimit,domainupperlimit));
                domainlowerlimit=domainupperlimit+1; 
               }
           }
       domainupperlimit=NX(); 
       domains.push_back(pair<int,int>(domainlowerlimit,domainupperlimit));
       SetFoundDomains(true);
      }

void DISCONTINUOUS::Grade(void) const
     { if(XLimited()==false){ XLimit();}
       if(FoundDomains()==false){ FindDomains();}

       CreateG(NX());
       for(int n=0;n<=NDomains()-1;n++)
          { int domainsize=domains[n].second-domains[n].first+1;
            if(domainsize>=5)
              { vector<double> XX(5), YY(5), D;
                int i,xoffset;
                for(i=domains[n].first;i<=domains[n].second;i++)
                      { if(i==domains[n].first){ xoffset=0;} 
                         else{ if(i==domains[n].first+1){ xoffset=-1;}  
                                     else{ if(i==domains[n].second-1){ xoffset=-3;} 
                                                 else{ if(i==domains[n].second){ xoffset=-4;} 
                                                             else{ xoffset=-2;} 
                                                           } 
                                                  } 
                                      }

                        for(int a=0;a<=4;a++){ XX[a]=_pX(i+a+xoffset)-_pX(i); YY[a]=_pY(i+a+xoffset)-_pY(i);}
                        D=FiniteDifference1D(0.,XX,YY);
                        _pSetG(i,D[1]);
                       }
               }
            else{ // don't find the gradient
                 }
           }
       SetGraded(true);
      }

void DISCONTINUOUS::Fit(void) const
     { if(XLimited()==false){ XLimit();}
       if(XDifferenced()==false){ XDifference();}
       if(YDifferenced()==false){ YDifference();}
       if(FoundDomains()==false){ FindDomains();}
       if(Graded()==false){ Grade();}
       CreateA(NX()-1);

       for(int n=0;n<=NDomains()-1;n++)
          { int domainsize=domains[n].second-domains[n].first+1;
            if(domainsize<=4)
              { for(int i=domains[n].first;i<=domains[n].second-1;i++)
   	           { CreateASet(i,2);
                     SetA(i,0,_pY(i)); 
                     SetA(i,1,_pDeltaY(i));
                    }
               }
            else{ int i; 
                  for(i=domains[n].first;i<=domains[n].second-1;i++)
           	              { CreateASet(i,4);
                                 SetA(i,0,_pY(i)); 
                                 SetA(i,1,_pDeltaX(i)*_pG(i)); 
                                 SetA(i,2,3.*_pDeltaY(i) -2.*_pDeltaX(i)*_pG(i) -_pDeltaX(i)*_pG(i+1)); 
                                 SetA(i,3,-2.*_pDeltaY(i) +_pDeltaX(i)*_pG(i) +_pDeltaX(i)*_pG(i+1));
                              }
                       }
           }

       SetFitted(true);
      }

void DISCONTINUOUS::LineFit(void) const
     { if(XLimited()==false){ XLimit();}
       if(XDifferenced()==false){ XDifference();}
       if(YDifferenced()==false){ YDifference();}
       if(FoundDomains()==false){ FindDomains();}
       CreateA(NX()-1);

       for(int n=0;n<=NDomains()-1;n++)
          { int i; 
            for(i=domains[n].first;i<=domains[n].second-1;i++)
    	           { CreateASet(i,2);
                      SetA(i,0,_pY(i)); 
                      SetA(i,1,_pDeltaY(i));
                    }
           }

       SetFitted(true);
      }


double DISCONTINUOUS::Interpolate(double X) const
                 { if(XLimited()==false){ XLimit();}
                   if(XDifferenced()==false){ XDifference();}
                   if(FoundDomains()==false){ FindDomains();}
                   if(Graded()==false){ Grade();}
                   if(Fitted()==false){ Fit();}
                   int i=XInterval(X);

                   return SPLINE_SINGLEXSET_SINGLEYSET::Interpolate(i,(X-_pX(i))/_pDeltaX(i));
                  }

double DISCONTINUOUS::Derivative(double X) const
                 { if(XLimited()==false){ XLimit();}
                   if(XDifferenced()==false){ XDifference();}
                   if(FoundDomains()==false){ FindDomains();}
                   if(Fitted()==false){ Fit();}
                   int i=XInterval(X);

                   return SPLINE_SINGLEXSET_SINGLEYSET::Derivative(i,(X-_pX(i))/_pDeltaX(i))/_pDeltaX(i);
                  }

// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************

DISCONTINUOUS operator+(const DISCONTINUOUS &D,const double &A){ return DISCONTINUOUS(OverlapXModifyY(D,A,std::plus<double>()));}
DISCONTINUOUS operator-(const DISCONTINUOUS &D,const double &A){ return DISCONTINUOUS(OverlapXModifyY(D,A,std::minus<double>()));}
DISCONTINUOUS operator*(const DISCONTINUOUS &D,const double &A){ return DISCONTINUOUS(OverlapXModifyY(D,A,std::multiplies<double>()));}
DISCONTINUOUS operator/(const DISCONTINUOUS &D,const double &A){ return DISCONTINUOUS(OverlapXModifyY(D,A,std::divides<double>()));}

DISCONTINUOUS operator+(const double &A,const DISCONTINUOUS &D){ return DISCONTINUOUS(OverlapXModifyY(A,D,std::plus<double>()));}
DISCONTINUOUS operator-(const double &A,const DISCONTINUOUS &D){ return DISCONTINUOUS(OverlapXModifyY(A,D,std::minus<double>()));}
DISCONTINUOUS operator*(const double &A,const DISCONTINUOUS &D){ return DISCONTINUOUS(OverlapXModifyY(A,D,std::multiplies<double>()));}

// ********************************************************************************************

DISCONTINUOUS operator+(DISCONTINUOUS const &D1,DISCONTINUOUS const &D2)    
         { return DISCONTINUOUS(OverlapXCombineY(D1,D2,std::plus<double>()));}

DISCONTINUOUS operator-(DISCONTINUOUS const &D1,DISCONTINUOUS const &D2)    
         { return DISCONTINUOUS(OverlapXCombineY(D1,D2,std::minus<double>()));}

DISCONTINUOUS operator*(DISCONTINUOUS const &D1,DISCONTINUOUS const &D2)  
         { return DISCONTINUOUS(OverlapXCombineY(D1,D2,std::multiplies<double>()));}

DISCONTINUOUS operator/(DISCONTINUOUS const &D1,DISCONTINUOUS const &D2)  
         { return DISCONTINUOUS(OverlapXCombineY(D1,D2,std::divides<double>()));}

} //end of namespace

