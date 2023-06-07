#if !defined(_POLYNOMIAL_TEMPLATES)
#define _POLYNOMIAL_TEMPLATES

#include "mstl.h"

template <typename Type>
Type PolynomialInterpolation(Type X,std::vector<Type> Xpoints,std::vector<Type> Ypoints,double ORDER)
       { if(Xpoints.size()!=Ypoints.size()){ throw DIFFERENT_LENGTHS("PolynomialInterpolation");}
         if(Xpoints.size()==1){ return Ypoints[0];}

         int i,j;
         std::vector<Type> RaisedXpoints(Xpoints.size());
         for(i=0;i<=(int)Xpoints.size()-1;i++){ RaisedXpoints[i]=pow(Xpoints[i],ORDER);}

         Sort(RaisedXpoints,Ypoints,descending);

         for(i=0;i<=(int)Xpoints.size()-2;i++)
            { j=1; Type mean=Ypoints[i];
              while(RaisedXpoints[i]-RaisedXpoints[i+1]==Zero<Type>())
                   { RaisedXpoints.erase(RaisedXpoints.begin()+i+1);
                     mean=( mean*j + Ypoints[i+1] )/(One<Type>()+j); j++;
                     Ypoints.erase(Ypoints.begin()+i+1);
                     if(i==(int)RaisedXpoints.size()-1){ break;}
                    }
              Ypoints[i]=mean;

              if(i==(int)RaisedXpoints.size()-1){ break;}
             }
         if(RaisedXpoints.size()==1){ return Ypoints[0];}

         std::vector<std::vector<Type> > P(RaisedXpoints.size()); for(i=0;i<=(int)P.size()-1;i++){ P[i]=std::vector<Type>(RaisedXpoints.size());}
         Type RaisedX=pow(X,ORDER);

         for(i=0;i<=(int)RaisedXpoints.size()-1;i++){ P[i][0]=Ypoints[i];}
         for(j=1;j<=(int)RaisedXpoints.size()-1;j++)
            { for(i=0;i<=(int)RaisedXpoints.size()-1-j;i++)
                 { P[i][j]=( (RaisedX-RaisedXpoints[i+j])*P[i][j-1] + (RaisedXpoints[i]-RaisedX)*P[i+1][j-1] )/(RaisedXpoints[i]-RaisedXpoints[i+j]);}
             }
         return P[0].back();
        }

template <typename Type>
Type PolynomialInterpolation(Type X,int N,Type *Xpoints,Type *Ypoints,double ORDER)
       { return PolynomialInterpolation(X,std::vector<Type>(Xpoints,Xpoints+N),std::vector<Type>(Ypoints,Ypoints+N),ORDER);}

#endif


