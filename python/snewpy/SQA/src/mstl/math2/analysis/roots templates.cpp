#if !defined(_ROOTS_TEMPLATES)
#define _ROOTS_TEMPLATES

#include "mstl.h"

//one dimension
template <class Functor>
double Root(double xlower,double xupper,const Functor &F,double C,double accuracy)
       { double a=xlower, b=(xlower+xupper)/2., c=xupper;
         double fa=F(a)-C;
         double fb=F(b)-C;
         double fc=F(c)-C;

         if(fa*fc>0.){ throw NO_SOLUTION("Root");}

         double P,dx,interval,priorendpoint;
         double Q,R,S,T;

         bool end; int counter=0;

         do{ end=false;
             if(fa==0.){ return a;}
             if(fb==0.){ return b;}
             if(fc==0.){ return c;}

             R=fb/fc, S=fb/fa, T=fa/fc;
             P=S*( T*(R-T)*(c-b) - (1.-R)*(b-a) );
             Q=(R-1.)*(S-1.)*(T-1.);

             if(Q!=0.){ dx=P/Q;
                        if(dx>0.)
                          { interval=(c-b);
                            if(std::fabs(dx)>interval){ dx=interval/2.;} 
                           }
                        else{ interval=std::fabs(b-a);
                              if(std::fabs(dx)>interval){ dx=-interval/2.;} 
                             }
                       }
             else{ if(fa*fb<0.){ interval=std::fabs(b-a); dx=-interval/2.;} 
                   else{ interval=std::fabs(c-b); dx=interval/2.;} 
                  }

             if(dx>0.){ priorendpoint=a; a=b; fa=fb;} else{ priorendpoint=c; c=b; fc=fb;}
             if(priorendpoint==b){ end=true;}

             b+=dx; fb=F(b)-C;
             ++counter;

             if(std::fabs(dx)<accuracy*std::fabs(interval)){ end=true;}
             else{ if(counter>=100){ throw DID_NOT_CONVERGE<double>(b,46,"Root");} }
            }while(end==false);

         return b;
        }

#endif
