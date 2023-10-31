#include "algorithm3.h"

using std::vector;
using std::abs;

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

double LineSearch(vector<double> x0,vector<double> deltax,vector<double(*)(vector<double>)> functions)
       { double lambda,lambda0,lambda1,delta,F0,dF;
         DCVECTOR A(2), deltaF(2);
         DMATRIX inverselambda(2,2);
         bool calculate0, calculate1;
	 
	 vector<double> x(x0.size());

         int NX=(int)x0.size(),NF=(int)functions.size();

         delta=0.;
         for(int i=0;i<=(int)x0.size()-1;i++){ delta+=abs( (x0[i]+deltax[i])-x0[i]);}
         if(Equality(delta,0.)==true){ return 0.;}

         // determine lambda by using the minimimum of the quadratic fit
         // to the sum of the squared function values at lambda = 1 and 1/2
         // i.e. deltaF=A[0].lambda^2 + A[1].lambda, deltaF[0] is for +1, deltaF[1] for 1/2
         F0=0.;
         for(int i=0;i<=NX-1;i++){ x[i]=x0[i];}
         for(int i=0;i<=NF-1;i++){ F0+=pow((*functions[i])(x),2.);}

         calculate0=true; lambda0=-1.;
         calculate1=true; lambda1=1.;

         int counter=0,countermax=5;
         bool finish=false;

         do{ if(calculate0==true){ deltaF[0]=-F0;
                                   for(int i=0;i<=NX-1;i++){ x[i]=x0[i]+lambda0*deltax[i];}
                                   for(int i=0;i<=NF-1;i++){ deltaF[0]+=pow((*functions[i])(x),2.);}
                                   calculate0=false;
                                  }
             if(calculate1==true){ deltaF[1]=-F0;
                                   for(int i=0;i<=NX-1;i++){ x[i]=x0[i]+lambda1*deltax[i]/2.;}
                                   for(int i=0;i<=NF-1;i++){ deltaF[1]+=pow((*functions[i])(x),2.);}
                                   calculate1=false;
                                  }

             if(lambda1<lambda0){ std::swap(lambda0,lambda1); std::swap(deltaF[0],deltaF[1]);} // order the lambdas so that lambda0<lambda1
             if(Equality(lambda0,lambda1)==true){ throw EQUAL_VALUES<double>(lambda0,"LineSearch");}

             inverselambda[0][0]=1./lambda0/(lambda0-lambda1);       inverselambda[0][1]=-1./lambda1/(lambda0-lambda1);
             inverselambda[1][0]=-lambda1/lambda0/(lambda0-lambda1); inverselambda[1][1]=lambda0/lambda1/(lambda0-lambda1);

             // find value of lambda that minimizes deltaF
             A=inverselambda*deltaF;

             lambda=0.;
             if(A[0]>0.){ lambda=-A[1]/2./A[0];}
             if(Equality(A[0],0.)==true && Equality(A[1],0.)==false){ if(deltaF[1]<0.){ lambda=lambda1/2.;} else{ lambda=lambda0/3.;} }
             if(Equality(A[0],0.)==true && Equality(A[1],0.)==true){ lambda=0.;}
             if(A[0]<0.){ if(deltaF[1]<deltaF[0]){ lambda=lambda1;} else{ lambda=lambda0;} }

             if(lambda>lambda1){ lambda=lambda1;}
             if(lambda<lambda0){ lambda=lambda0;}

             delta=0.; dF=0.;
             if(Equality(lambda,0.)==false){ for(int i=0;i<=NX-1;i++){ delta+=abs( (x0[i]+lambda*deltax[i])-x0[i]);} }
             if(Equality(lambda,0.)==true || (Equality(lambda,lambda0)==true && deltaF[0]<0.) || (Equality(lambda,lambda1)==true && deltaF[1]<0.) )
               { finish=true;}
             else{ if(Equality(lambda,lambda0)==false && Equality(lambda,lambda1)==false)
                     { if(Equality(delta,0.)==false)
                         { dF=-F0;
                           for(int i=0;i<=NX-1;i++){ x[i]=x0[i]+lambda*deltax[i];}
                           for(int i=0;i<=NF-1;i++){ dF+=pow((*functions[i])(x),2.);}
                           if(dF<0.){ finish=true;}
                           else{ counter++; finish=false;
                                 if(lambda<0.){ lambda0=lambda; deltaF[0]=dF;} else{ lambda1=lambda; deltaF[1]=dF;}
                                }
                          }
                       else{ finish=true;}
                      }
                   else{ counter++; finish=false;
                         if(Equality(lambda,lambda0)==true){ lambda0+=lambda0; calculate0=true;}
                         if(Equality(lambda,lambda1)==true){ lambda1+=lambda1; calculate1=true;}
                        }
                   if(counter>=countermax){ finish=true;}
                  }
             if(counter>=countermax){ finish=true;}
            }while(finish==false);

        if(counter>=countermax){ throw NO_SOLUTION("LineSearch");}

        return lambda;
       }

//*****************************************************************************************************

double LineSearch(vector<double> x0,vector<double> deltax,vector<double>(*functions)(vector<double>))
       { double lambda,lambda0,lambda1,delta,F0,dF;
         DCVECTOR A(2), deltaF(2);
         DMATRIX inverselambda(2,2);
         bool calculate0, calculate1;

         int NX=(int)x0.size(),NF;
	 
	 vector<double> x(NX), ff;

         delta=0.;
         for(int i=0;i<=NX-1;i++){ delta+=abs( (x0[i]+deltax[i])-x0[i]);}
         if(Equality(delta,0.)==true){ return 0.;}

         // determine lambda by using the minimimum of the quadratic fit
         // to the sum of the squared function values at lambda = 1 and 1/2
         // i.e. deltaF=A[0].lambda^2 + A[1].lambda, deltaF[0] is for +1, deltaF[1] for 1/2
         F0=0.;
         x=x0;
         ff=(*functions)(x); 
         NF=(int)ff.size(); 
         for(int i=0;i<=NF-1;i++){ F0+=pow(ff[i],2.);}

         calculate0=true; lambda0=-1.;
         calculate1=true; lambda1=1.;

         int counter=0,countermax=5;
         bool finish=false;

         do{ if(calculate0==true){ deltaF[0]=-F0; 
                                   for(int i=0;i<=NX-1;i++){ x[i]=x0[i]+lambda0*deltax[i];}
                                   ff=(*functions)(x); for(int i=0;i<=NF-1;i++){ deltaF[0]+=pow(ff[i],2.);}
                                   calculate0=false;
                                  }
             if(calculate1==true){ deltaF[1]=-F0;
                                   for(int i=0;i<=NX-1;i++){ x[i]=x0[i]+lambda1*deltax[i]/2.;}
                                   ff=(*functions)(x); for(int i=0;i<=NF-1;i++){ deltaF[1]+=pow(ff[i],2.);}
                                   calculate1=false;
                                  }

             if(lambda1<lambda0){ std::swap(lambda0,lambda1); std::swap(deltaF[0],deltaF[1]);} // order the lambdas so that lambda0<lambda1
             if(Equality(lambda0,lambda1)==true){ throw EQUAL_VALUES<double>(lambda0,"LineSearch");}

             inverselambda[0][0]=1./lambda0/(lambda0-lambda1);       inverselambda[0][1]=-1./lambda1/(lambda0-lambda1);
             inverselambda[1][0]=-lambda1/lambda0/(lambda0-lambda1); inverselambda[1][1]=lambda0/lambda1/(lambda0-lambda1);

             // find value of lambda that minimizes deltaF
             A=inverselambda*deltaF;

             lambda=0.;
             if(A[0]>0.){ lambda=-A[1]/2./A[0];}
             if(Equality(A[0],0.)==true && Equality(A[1],0.)==false){ if(deltaF[1]<0.){ lambda=lambda1/2.;} else{ lambda=lambda0/3.;} }
             if(Equality(A[0],0.)==true && Equality(A[1],0.)==true){ lambda=0.;}
             if(A[0]<0.){ if(deltaF[1]<deltaF[0]){ lambda=lambda1;} else{ lambda=lambda0;} }

             if(lambda>lambda1){ lambda=lambda1;}
             if(lambda<lambda0){ lambda=lambda0;}

             delta=0.; dF=0.;
             if(Equality(lambda,0.)==false){ for(int i=0;i<=NX-1;i++){ delta+=abs( (x0[i]+lambda*deltax[i])-x0[i]);} }
             if(Equality(lambda,0.)==true || (Equality(lambda,lambda0)==true && deltaF[0]<0.) || (Equality(lambda,lambda1)==true && deltaF[1]<0.) )
               { finish=true;}
             else{ if(Equality(lambda,lambda0)==false && Equality(lambda,lambda1)==false)
                     { if(Equality(delta,0.)==false)
                         { dF=-F0;
                           for(int i=0;i<=NX-1;i++){ x[i]=x0[i]+lambda*deltax[i];}
                           ff=(*functions)(x); for(int i=0;i<=NF-1;i++){ dF+=pow(ff[i],2.);}
                           if(dF<0.){ finish=true;}
                           else{ counter++; finish=false;
                                 if(lambda<0.){ lambda0=lambda; deltaF[0]=dF;} else{ lambda1=lambda; deltaF[1]=dF;}
                                }
                          }
                       else{ finish=true;}
                      }
                   else{ counter++; finish=false;
                         if(Equality(lambda,lambda0)==true){ lambda0+=lambda0; calculate0=true;}
                         if(Equality(lambda,lambda1)==true){ lambda1+=lambda1; calculate1=true;}
                        }
                   if(counter>=countermax){ finish=true;}
                  }
             if(counter>=countermax){ finish=true;}
            }while(finish==false);

        if(counter>=countermax){ throw NO_SOLUTION("LineSearch");}

        return lambda;
       }

