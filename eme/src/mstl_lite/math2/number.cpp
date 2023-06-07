#include "mstl.h"

//******************************************************************************
//******************************************************************************
//******************************************************************************

/*
unsigned int Fibonacci(unsigned int N)
         { if(N==1 || N==2){ return 1;}
           unsigned int n=2, F, Fminus1=Fminus2=1;
           do{ F=Fminus1+Fminus2; n++; Fminus2=Fminus1; Fminus1=F;}
           while(n<N);
           return F;
          }

double Fibonacci(double N){ return ( pow((1.+sqrt(5.))/2.,N) - pow((1.-sqrt(5.))/2.,N) )/sqrt(5.);}
*/
unsigned int Graycode(unsigned int n) { return n ^ (n >> 1);}
