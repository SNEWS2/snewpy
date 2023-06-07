#if !defined(_INTERPOLATION_TEMPLATES)
#define _INTERPOLATION_TEMPLATES

#include "mstl.h"

template <typename XType,typename YType>
YType PolynomialInterpolation(XType X,std::vector<XType> Xpoints,std::vector<YType> Ypoints,double ORDER)
       { try{ return RationalInterpolation(X,Xpoints,Ypoints,0,ORDER);}
         catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("YType PolynomialInterpolation(YType,vector<XType>,vector<YType>"); throw DL;}
        }

template <typename XType,typename YType>
YType PolynomialInterpolation(XType X,int N,XType *Xpoints,YType *Ypoints,double ORDER)
       { return PolynomialInterpolation(X,std::vector<XType>(Xpoints,Xpoints+N),std::vector<YType>(Ypoints,Ypoints+N),ORDER);}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

template <typename XType,typename YType>
YType RationalInterpolation(XType X,std::vector<XType> Xpoints,std::vector<YType> Ypoints,int DENEXP,double XPOWER)
       { if(Xpoints.size()!=Ypoints.size()){ throw DIFFERENT_LENGTHS("RationalInterpolation");}
         if(Xpoints.size()==1){ return Ypoints[0];}

         int i,j,k;
         std::vector<XType> RaisedXpoints(Xpoints.size());
         for(i=0;i<=(int)Xpoints.size()-1;i++){ RaisedXpoints[i]=Xpoints[i]; for(int j=2;j<=XPOWER;j++){ RaisedXpoints[i]*=Xpoints[i];} }

         /*Sort(RaisedXpoints,Ypoints,descending);

         // go through the points and remove duplicates in order to
         // avoid singularities in the later part of the algorithm
         // when duplicates of x are found average the y points.
         // actually if the y points are not identical then this is an error
         // and the polynomial is not one in x^ORDER
         for(i=0;i<=(int)Xpoints.size()-2;i++)
            { YType j=One(j); YType mean=Ypoints[i];
              while(Equality(RaisedXpoints[i]-RaisedXpoints[i+1],Zero<XType>())==true)
                   { RaisedXpoints.erase(RaisedXpoints.begin()+i+1);
                     mean=( mean*j + Ypoints[i+1] )/(One<YType>()+j); j+=One(j);
                     Ypoints.erase(Ypoints.begin()+i+1);
                     if(i==(int)RaisedXpoints.size()-1){ break;}
                    }
              Ypoints[i]=mean;

              if(i==(int)RaisedXpoints.size()-1){ break;}
             }*/

         int N=RaisedXpoints.size(), NUMEXP=N-1-DENEXP;
         if(N==1){ return Ypoints[0];}

         std::vector<std::vector<std::vector<YType> > > P(NUMEXP+1,std::vector<std::vector<YType> >(DENEXP+1));
         for(i=0;i<=NUMEXP;i++){ for(j=0;j<=DENEXP;j++){ P[i][j]=std::vector<YType>(N-i-j);} }

         XType RaisedX=pow(X,XPOWER);
         std::vector<XType> Alpha(N);
         for(k=0;k<=N-1;k++){ Alpha[k]=RaisedX-RaisedXpoints[k];}

         // initialize the first column of the P array
         for(k=0;k<=N-1;k++){ P[0][0][k]=Ypoints[k];}

         i=j=0;
         // walk up the power in the numerator
         if(NUMEXP>DENEXP)
           { while(i<NUMEXP-DENEXP)
                  { i++;
                    for(k=0;k<=(int)P[i][j].size()-1;k++)
                       { P[i][j][k]=( Alpha[k]*P[i-1][j][k+1] - Alpha[i+k]*P[i-1][j][k] )
                                   /( Alpha[k]-Alpha[i+k] );
                        }
                   };
            }
         // or walk up the power in the denominator
         if(DENEXP>NUMEXP)
           { while(j<DENEXP-NUMEXP)
                  { j++;
                    for(k=0;k<=(int)P[i][j].size()-1;k++)
                       { if(Equality(Alpha[k]*P[i][j-1][k] - Alpha[j+k]*P[i][j-1][k],Zero<YType>())==true)
                           { throw POLE<XType>(X,"RationalInterpolation");}
                         P[i][j][k]=(Alpha[k]-Alpha[j+k]) * P[i][j-1][k+1]*P[i][j-1][k]
                                   /( Alpha[k]*P[i][j-1][k] - Alpha[j+k]*P[i][j-1][k+1] );
                        }
                   };
            }

         // now step up power in denominator followed by step in numerator
         while(i+j<NUMEXP+DENEXP)
            { j++;
              if(i==0){ for(k=0;k<=(int)P[i][j].size()-1;k++)
                           { if(Equality(Alpha[k]*P[i][j-1][k] - Alpha[j+k]*P[i][j-1][k+1],Zero<YType>())==true)
                               { throw POLE<XType>(X,"RationalInterpolation");}
                             P[i][j][k]=(Alpha[k]-Alpha[j+k]) * P[i][j-1][k+1]*P[i][j-1][k]
                                       /( Alpha[k]*P[i][j-1][k] - Alpha[j+k]*P[i][j-1][k+1] );
                            }
                       }
              else{ for(k=0;k<=(int)P[i][j].size()-1;k++)
                       { if(Equality(Alpha[k]*(P[i][j-1][k]-P[i-1][j-1][k+1]) - Alpha[i+j+k]*(P[i][j-1][k+1]-P[i-1][j-1][k+1]),Zero<YType>())==true)
                           { throw POLE<XType>(X,"RationalInterpolation");}
                         P[i][j][k]=( Alpha[k]*P[i][j-1][k+1]*(P[i][j-1][k]-P[i-1][j-1][k+1]) - Alpha[i+j+k]*P[i][j-1][k]*(P[i][j-1][k+1]-P[i-1][j-1][k+1]) )
                                   /( Alpha[k]*(P[i][j-1][k]-P[i-1][j-1][k+1]) - Alpha[i+j+k]*(P[i][j-1][k+1]-P[i-1][j-1][k+1]) );
                        }
                   }

              i++;
              for(k=0;k<=(int)P[i][j].size()-1;k++)
                 { if(Equality(Alpha[k]*(P[i-1][j][k]-P[i-1][j-1][k+1]) - Alpha[i+j+k]*(P[i-1][j][k+1]-P[i-1][j-1][k+1]),Zero<YType>())==true)
                     { throw POLE<XType>(X,"RationalInterpolation");}
                       P[i][j][k]=( Alpha[k]*P[i-1][j][k+1]*(P[i-1][j][k]-P[i-1][j-1][k+1]) - Alpha[i+j+k]*P[i-1][j][k]*(P[i-1][j][k+1]-P[i-1][j-1][k+1]) )
                                 /( Alpha[k]*(P[i-1][j][k]-P[i-1][j-1][k+1]) - Alpha[i+j+k]*(P[i-1][j][k+1]-P[i-1][j-1][k+1]) );
                  }
             };
         return P[NUMEXP][DENEXP][0];
        }

template <typename XType,typename YType>
YType RationalInterpolation(XType X,int N,XType *Xpoints,YType *Ypoints,int DENORDER,double ORDER)
       { return RationalInterpolation(X,std::vector<XType>(Xpoints,Xpoints+N),std::vector<YType>(Ypoints,Ypoints+N),DENORDER,ORDER);}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// The ORDER parameter controls the powers of x in the polynomial
template <typename XType,typename YType>
YType DiagonalRationalInterpolation(XType X,std::vector<XType> Xpoints,std::vector<YType> Ypoints,double XPOWER)
       { try{ int DENEXP;
              if(Odd(Xpoints.size())==true){ DENEXP=(Xpoints.size()-1)/2;} else{ DENEXP=Xpoints.size()/2;}
              return RationalInterpolation(X,Xpoints,Ypoints,DENEXP,XPOWER);
             }
         catch(DIFFERENT_LENGTHS &DL){ DL.ChangeFunction("DiagonalRationalInterpolation(Type,vector<Type>,vector<Type>"); throw DL;}
         catch(POLE<XType> &P){ P.ChangeFunction("DiagonalRationalInterpolation(Type,vector<Type>,vector<Type>"); throw P;}
        }

template <typename XType,typename YType>
YType DiagonalRationalInterpolation(XType X,int N,XType *Xpoints,YType *Ypoints,double XPOWER)
       { return DiagonalRationalInterpolation(X,std::vector<XType>(Xpoints,Xpoints+N),std::vector<YType>(Ypoints,Ypoints+N),XPOWER);}

#endif


