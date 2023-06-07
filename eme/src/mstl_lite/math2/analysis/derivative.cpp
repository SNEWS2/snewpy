#include "derivative.h"

using std::vector;

//*******************************************************************************
//********************* Finite Difference Formulae ******************************
//*******************************************************************************

vector<double> FiniteDifference1D(double x0,vector<double> x,const vector<double> &y)
       { for(int i=0;i<=(int)x.size()-1;i++){ x[i]-=x0;}
         DCVECTOR FD(FactorialMatrix(x.size())*VandermondeMatrixInverse(x)*DCVECTOR(y));
         vector<double> fd(x.size());
         for(int i=0;i<=(int)x.size()-1;i++){ fd[i]=FD[i];}
         return fd;
        }

vector<vector<double> > FiniteDifference2D(double x0,double y0,vector<double> x,vector<double> y,const vector<vector<double> > &z)
       { for(int i=0;i<=(int)x.size()-1;i++){ x[i]-=x0;}
         for(int j=0;j<=(int)y.size()-1;j++){ y[j]-=y0;}

         DMATRIX FD(FactorialMatrix(x.size())*VandermondeMatrixInverse(x)*DMATRIX(z)*Transpose(VandermondeMatrixInverse(y))*FactorialMatrix(y.size()));
         vector<vector<double> > fd(x.size(),vector<double>(y.size()));
         for(int i=0;i<=(int)x.size()-1;i++){ for(int j=0;j<=(int)y.size()-1;j++){ fd[i][j]=FD[i][j];} }
         return fd;
        }




