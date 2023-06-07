#if !defined(_RUNGE_KUTTA_TEMPLATES)
#define _RUNGE_KUTTA_TEMPLATES

#include "mstl.h"

template <typename XType,typename YType,typename Functor>
YType RungeKuttaCashKarp(const Functor &dydx,YType &error,XType x0,YType y0,XType dx)
       { int NRK,NOrder;
         const double *a,**b,*c,*d;
         RungeKuttaCashKarpParameters(NRK,NOrder,a,b,c,d);

         return RungeKutta(dydx,error,x0,y0,dx,NRK,a,b,c,d);
        }

// *******************************************************
// *******************************************************
// *******************************************************

template <typename XType,typename YType,typename Functor>
YType RungeKutta(const Functor &dydx,YType &error,XType x0,YType y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D)
          { int m,n;
            XType x;
            YType y, yincrement;
            std::vector<YType> K(N,y0); 

            try{ for(m=0;m<=N-1;m++)
                    { x=x0+A[m]*dx;
                      y=y0;
                      for(n=0;n<=m-1;n++){ y+=B[m][n]*K[n];}
                      K[m]=dx*dydx(x,y);
                     }

                 yincrement=error=Zero(y);
                 for(m=0;m<=N-1;m++){ yincrement+=C[m]*K[m]; error+=(C[m]-D[m])*K[m];}

                }catch(OUT_OF_RANGE<int> R){ R.ChangeWhich("RungeKutta"); throw R;}

            return yincrement;
           }

// ******************

template <typename XType,typename YType,typename Functor>
std::vector<YType> RungeKutta(std::vector<Functor> &dydx,std::vector<YType> &error,XType x0,std::vector<YType> y0,XType dx,int N,const double *A,const double **B,const double *C,const double *D)
          { int Ny=y0.size();
            XType x;
            std::vector<YType> y(Ny), yincrement(Ny);
            std::vector<std::vector<YType> > K(N,y0); 

            if(Ny!=(int)dydx.size() || Ny!=(int)error.size()){ throw DIFFERENT_LENGTHS("RungeKutta");}

            int l,m,n;

            try{ for(l=0;l<=N-1;l++)
                    { x=x0+A[l]*dx;
                      y=y0;
                      for(m=0;m<=Ny-1;m++){ for(n=0;n<=l-1;n++){ y[m]+=B[l][n]*K[n][m];} }
                      for(m=0;m<=Ny-1;m++){ K[l][m]=dx*dydx[m](x,y);}
                     }

                 yincrement=error=Zero(y);
                 for(m=0;m<=Ny-1;m++){ for(l=0;l<=N-1;l++){ yincrement[m]+=C[l]*K[l][m]; error[m]+=(C[l]-D[l])*K[l][m];} }

                }catch(OUT_OF_RANGE<int> R){ R.ChangeWhich("RungeKutta"); throw R;}

            return yincrement;
           }

// ******************

template <typename XType,typename YType,typename Functor>
std::vector<YType> RungeKutta(const Functor &,std::vector<YType> &error,XType x0,std::vector<YType> y0,XType dx,
                              int N,const double *A,const double **B,const double *C,const double *D)
       { int Ny=y0.size();
         XType x;
         std::vector<YType> y(Ny), yincrement(Ny), fvalues;
         std::vector<std::vector<YType> > K(N,y0); 

         if(Ny!=(int)error.size() || Ny!=(int)y0.size()){ throw DIFFERENT_LENGTHS("RungeKutta");}

         int l,m,n;
         try{ for(l=0;l<=N-1;l++)
                 { x=x0+A[l]*dx;
                   y=y0;
                   for(m=0;m<=Ny-1;m++){ for(n=0;n<=l-1;n++){ y[m]+=B[l][n]*K[n][m];} }
                   fvalues=dydx(x,y);
                   for(m=0;m<=Ny-1;m++){ K[l][m]=dx*fvalues[m];}
                  }

              yincrement=error=Zero(y);
              for(m=0;m<=Ny-1;m++){ for(l=0;l<=N-1;l++){ yincrement[m]+=C[l]*K[l][m]; error[m]+=(C[l]-D[l])*K[l][m];} }

             }catch(OUT_OF_RANGE<int> R){ R.ChangeWhich("RungeKutta"); throw R;}

         return yincrement;
        }

// ******************

// for the case where a function returns a multicolumn list rather than a list of functions
template <typename XType,typename YType,typename Functor>
std::vector<std::vector<YType> > RungeKutta(const Functor &dydx,std::vector<std::vector<YType> > &error,XType x0,std::vector<std::vector<YType> > y0,XType dx,
                                            int N,const double *A,const double **B,const double *C,const double *D)
       { int k,l,m,n;

         int Ny1=(int)y0.size();
         std::vector<int> Ny2(Ny1);                  
         for(l=0;l<=Ny1-1;l++){ Ny2[l]=(int)y0[l].size();}

         //if(Ny1!=(int)error.size() || Ny1!=(int)y0.size()){ throw DIFFERENT_LENGTHS("RungeKutta");}
         //for(l=0;l<=Ny1-1;l++){ if(Ny2[l]!=(int)error[l].size() || Ny2[l]!=(int)y0[l].size()){ throw DIFFERENT_LENGTHS("RungeKutta");} }

         XType x;
         std::vector<std::vector<YType> > y(y0), yincrement(Zero(y0));
         std::vector<std::vector<std::vector<YType> > > K(N,y0); 
         std::vector<std::vector<YType> > fvalues;
         
         try{ for(k=0;k<=N-1;k++)
                 { x=x0+A[k]*dx;
                   y=y0; 
                   for(l=0;l<=Ny1-1;l++){ for(m=0;m<=Ny2[l]-1;m++){ for(n=0;n<=k-1;n++){ y[l][m]+=B[k][n]*K[n][l][m];} } }
                   fvalues=dydx(x,y); 

                   for(l=0;l<=Ny1-1;l++){ for(m=0;m<=Ny2[l]-1;m++){ K[k][l][m]=dx*fvalues[l][m];} }
                  }

              error=Zero(y);
              for(l=0;l<=Ny1-1;l++){ for(m=0;m<=Ny2[l]-1;m++){ for(k=0;k<=N-1;k++){ yincrement[l][m]+=C[k]*K[k][l][m]; error[l][m]+=(C[k]-D[k])*K[k][l][m];} } }

             }catch(OUT_OF_RANGE<int> R){ R.ChangeWhich("RungeKutta"); throw R;}

         return yincrement;
        }

// ******************

// for the case where a function returns a multicolumn list rather than a list of functions
template <typename XType,typename YType,typename Functor>
std::vector<std::vector<std::vector<YType> > > RungeKutta(const Functor &dydx,std::vector<std::vector<std::vector<YType> > > &error,XType x0,std::vector<std::vector<std::vector<YType> > > y0,XType dx,
                                            int NRK,const double *A,const double **B,const double *C,const double *D)
       { int j,k,l,m,n;

         XType x;
         std::vector<std::vector<std::vector<YType> > > y(y0), dy(Zero(y0));         
         std::vector<std::vector<std::vector<std::vector<YType> > > > K(NRK,y0); 
         std::vector<std::vector<std::vector<YType> > > fvalues;

         try{ for(k=0;k<=NRK-1;k++)
                 { x=x0+A[k]*dx;
                   y=y0; 
                   for(l=0;l<=(int)y.size()-1;l++){ for(m=0;m<=(int)y[l].size()-1;m++){ for(j=0;j<=(int)y[l][m].size()-1;j++){ for(n=0;n<=k-1;n++){ y[l][m][j]+=B[k][n]*K[n][l][m][j];} } } }
                   fvalues=dydx(x,y); 

                   for(l=0;l<=(int)y.size()-1;l++){ for(m=0;m<=(int)y[l].size()-1;m++){ for(j=0;j<=(int)y[l][m].size()-1;j++){ K[k][l][m][j]=dx*fvalues[l][m][j];} } }
                  }

              for(l=0;l<=(int)y.size()-1;l++){ for(m=0;m<=(int)y[l].size()-1;m++){ for(j=0;j<=(int)y[l][m].size()-1;j++){ error[l][m][j]=Zero<YType>();} } }
              for(l=0;l<=(int)y.size()-1;l++){ for(m=0;m<=(int)y[l].size()-1;m++){ for(j=0;j<=(int)y[l][m].size()-1;j++){ for(k=0;k<=NRK-1;k++){ dy[l][m][j]+=C[k]*K[k][l][m][j]; error[l][m][j]+=(C[k]-D[k])*K[k][l][m][j];} } } }

             }catch(OUT_OF_RANGE<int> R){ R.ChangeWhich("RungeKutta"); throw R;}

         return dy;
        }

#endif
