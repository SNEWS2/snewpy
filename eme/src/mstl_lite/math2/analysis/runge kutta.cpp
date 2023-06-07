#include "runge kutta.h"

using std::vector;

//*******************************************************
//*******************************************************
//*******************************************************

void RungeKuttaCashKarpParameters(int &NRK,int &NOrder,const double* &A,const double** &B,const double* &C,const double* &D)
     { NRK=6; NOrder=5;

       static const double a[]={ 0., 1./5., 3./10., 3./5., 1., 7./8. };
       static const double b0[]={};
       static const double b1[]={ 1./5. };
       static const double b2[]={ 3./40.,9./40. };
       static const double b3[]={ 3./10.,-9./10.,6./5. };
       static const double b4[]={ -11./54.,5./2.,-70./27.,35./27. };
       static const double b5[]={ 1631./55296.,175./512.,575./13824.,44275./110592.,253./4096. };
       static const double* b[]={ b0,b1,b2,b3,b4,b5 };
       static const double c[]={ 37./378.,0.,250./621.,125./594.,0.,512./1771. };
       static const double d[]={ 2825./27648.,0.,18575./48384.,13525./55296.,277./14336.,1./4. };

       A=a; B=b; C=c; D=d;
      } 

