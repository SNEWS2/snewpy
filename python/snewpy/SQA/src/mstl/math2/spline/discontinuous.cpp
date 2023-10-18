#include "discontinuous.h"

namespace interpolation{

using std::vector;
using std::pair;

void DISCONTINUOUS::FindDomains(void) const
     { domains.clear();
       int domainlowerlimit=1, domainupperlimit;
       for(int i=1;i<=static_cast<int>(NX())-1;i++)
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
       for(int n=0;n<=static_cast<int>(NDomains())-1;n++)
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
                        SetG(i,D[1]);
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

       for(int n=0;n<=static_cast<int>(NDomains())-1;n++)
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
                                 SetA(i,1,_pDeltaX(i)*G(i)); 
                                 SetA(i,2,3.*_pDeltaY(i) -2.*_pDeltaX(i)*G(i) -_pDeltaX(i)*G(i+1)); 
                                 SetA(i,3,-2.*_pDeltaY(i) +_pDeltaX(i)*G(i) +_pDeltaX(i)*G(i+1));
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

       for(int n=0;n<=static_cast<int>(NDomains())-1;n++)
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

DISCONTINUOUS DISCONTINUOUS::Derivative(void) const
                 { if(Fitted()==false){ Fit();}

                   DISCONTINUOUS dfdx(N());
                   for(int i=1;i<=static_cast<int>(N());i++){ dfdx.SetX(i,X(i)); dfdx.SetY(i,Derivative(X(i)));}                       

                   return dfdx;
                  }

// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************

DISCONTINUOUS_MULTIPLESETS::DISCONTINUOUS_MULTIPLESETS(vector<vector<double> > XY)
         { size_t NSets=(XY).size()-1;
           size_t N=(XY[0]).size();
           CreateX(N); CreateY(N);
           _pSetX(XY[0]); 
           int i; 
           for(i=1;i<=static_cast<int>(NSets);i++){ _pSetY(i,XY[i]);}
           SetXDifferenced(false); 
           SetGraded(false); 
           SetFitted(false); 
           SetSorted(false); 
           SetXLimited(false);
           SetYLimited(false);  
           SetXYLimited(false); 
           SetFoundDomains(false); 
          }

DISCONTINUOUS_MULTIPLESETS::DISCONTINUOUS_MULTIPLESETS(vector<double> X,vector<vector<double> > Y)
         { size_t N=X.size();
           size_t NSets=Y.size();
           CreateX(N); CreateY(NSets,N);
           _pSetX(X); 
           int i; 
           for(i=1;i<=static_cast<int>(NSets);i++){ _pSetY(i,Y[i-1]);}
           SetXDifferenced(false); 
           SetGraded(false); 
           SetFitted(false); 
           SetSorted(false); 
           SetXLimited(false);
           SetYLimited(false);  
           SetXYLimited(false); 
           SetFoundDomains(false); 
          }

void DISCONTINUOUS_MULTIPLESETS::FindDomains(void) const
     { domains.clear();
       int domainlowerlimit=1, domainupperlimit;
       for(int i=1;i<=static_cast<int>(NX())-1;i++)
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

void DISCONTINUOUS_MULTIPLESETS::Grade(void) const
     { if(XLimited()==false){ XLimit();}
       if(FoundDomains()==false){ FindDomains();}

       CreateG(NSets(),NX());
       for(int n=0;n<=static_cast<int>(NDomains())-1;n++)
          { int domainsize=domains[n].second-domains[n].first+1;
            if(domainsize>=5)
              { vector<double> XX(5), YY(5), D;
                 int i,j,a,xoffset;
                 for(i=domains[n].first;i<=domains[n].second;i++)
                       { if(i==domains[n].first){ xoffset=0;} 
                          else{ if(i==domains[n].first+1){ xoffset=-1;}  
                                      else{ if(i==domains[n].second-1){ xoffset=-3;} 
                                                 else{ if(i==domains[n].second){ xoffset=-4;} 
                                                             else{ xoffset=-2;} 
                                                            } 
                                                  } 
                                     }

                          for(j=1;j<=static_cast<int>(NSets());j++)
                                { for(a=0;a<=4;a++){ XX[a]=_pX(i+a+xoffset)-_pX(i); YY[a]=_pY(j,i+a+xoffset)-_pY(j,i);}
                                   D=FiniteDifference1D(0.,XX,YY);
                                   SetG(j,i,D[1]);
                                 }
                         }
                }
            else{ // don't find the gradient
                     }
           }
       SetGraded(true);
      }

void DISCONTINUOUS_MULTIPLESETS::Fit(void) const
     { if(XLimited()==false){ XLimit();}
       if(XDifferenced()==false){ XDifference();}
       if(YDifferenced()==false){ YDifference();}
       if(FoundDomains()==false){ FindDomains();}
       if(Graded()==false){ Grade();}
       CreateA(NSets(),NX()-1);

       for(int n=0;n<=static_cast<int>(NDomains())-1;n++)
          { int domainsize=domains[n].second-domains[n].first+1;
            if(domainsize<=4)
              { for(int i=domains[n].first;i<=domains[n].second-1;i++)
   	           { for(int j=1;j<=static_cast<int>(NSets());j++)
                           { CreateASet(j,i,2);
                              SetA(j,i,0,_pY(j,i)); 
                              SetA(j,i,1,_pDeltaY(j,i));
                            }
                    }
               }
            else{ int i;
                  for(i=domains[n].first;i<=domains[n].second-1;i++)
         	              { for(int j=1;j<=static_cast<int>(NSets());j++)
                                       { CreateASet(j,i,4);
                                          SetA(j,i,0,_pY(j,i)); 
                                          SetA(j,i,1,_pDeltaX(i)*G(j,i)); 
                                          SetA(j,i,2,3.*_pDeltaY(j,i) -2.*_pDeltaX(i)*G(j,i) -_pDeltaX(i)*G(j,i+1)); 
                                          SetA(j,i,3,-2.*_pDeltaY(j,i) +_pDeltaX(i)*G(j,i) +_pDeltaX(i)*G(j,i+1));
                                        }
                             }
                     }
           }

       SetFitted(true);
      }

void DISCONTINUOUS_MULTIPLESETS::LineFit(void) const
     { if(XLimited()==false){ XLimit();}
       if(XDifferenced()==false){ XDifference();}
       if(YDifferenced()==false){ YDifference();}
       if(FoundDomains()==false){ FindDomains();}
       CreateA(NSets(),N()-1);

       for(int n=0;n<=static_cast<int>(NDomains())-1;n++)
          { int i;
            for(i=domains[n].first;i<=domains[n].second-1;i++)
                   { for(int j=1;j<=static_cast<int>(NSets());j++)
           	           { CreateASet(j,i,2);
                              SetA(j,i,0,_pY(j,i)); 
                              SetA(j,i,1,_pDeltaY(j,i));
                             }
                }
           }

       SetFitted(true);
      }

// ********************************

double DISCONTINUOUS_MULTIPLESETS::Interpolate(int j,double X) const
       { if(XLimited()==false){ XLimit();}
         if(XDifferenced()==false){ XDifference();}
         if(FoundDomains()==false){ FindDomains();}
         if(Graded()==false){ Grade();}
         if(Fitted()==false){ Fit();}
         int i=XInterval(X);

         return SPLINE_SINGLEXSET_MULTIPLEYSETS::Interpolate(j,i,(X-_pX(i))/_pDeltaX(i));
        }

vector<double> DISCONTINUOUS_MULTIPLESETS::Interpolate(double X) const
       { if(XLimited()==false){ XLimit();}
         if(XDifferenced()==false){ XDifference();}
         if(FoundDomains()==false){ FindDomains();}
         if(Graded()==false){ Grade();}
         if(Fitted()==false){ Fit();}
         int i=XInterval(X);

         return SPLINE_SINGLEXSET_MULTIPLEYSETS::Interpolate(i,(X-_pX(i))/_pDeltaX(i));
        }

double DISCONTINUOUS_MULTIPLESETS::Derivative(int j,double X) const
       { if(XLimited()==false){ XLimit();}
         if(XDifferenced()==false){ XDifference();}
         if(FoundDomains()==false){ FindDomains();}
         if(Fitted()==false){ Fit();}
         int i=XInterval(X);

         return SPLINE_SINGLEXSET_MULTIPLEYSETS::Derivative( j,i,(X-_pX(i))/_pDeltaX(i) )/_pDeltaX(i);
        }

vector<double> DISCONTINUOUS_MULTIPLESETS::Derivative(double X) const
       { if(XLimited()==false){ XLimit();}
         if(XDifferenced()==false){ XDifference();}
         if(FoundDomains()==false){ FindDomains();}
         if(Fitted()==false){ Fit();}
         int i=XInterval(X),j;

         vector<double> D=SPLINE_SINGLEXSET_MULTIPLEYSETS::Derivative(i,(X-_pX(i))/_pDeltaX(i) );
         for(j=0;j<=static_cast<int>(NSets())-1;j++){ D[j]/=_pDeltaX(i);}

         return D;
        }

} //end of namespace

