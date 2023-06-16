
#include "jacobians.h"

// *********************************************************************

using std::array;

// ***************************************************************

MATRIX<double,10,8> J(array<double,NY> Y)
       { MATRIX<double,10,8> j;

         double cPsi1=cos(Y[0]),sPsi1=sin(Y[0]), cPsi2=cos(Y[1]),sPsi2=sin(Y[1]);
         double cPsi3=cos(Y[2]),sPsi3=sin(Y[2]), cPsi4=cos(Y[3]),sPsi4=sin(Y[3]), cPsi5=cos(Y[4]),sPsi5=sin(Y[4]);
         double cPsi6=cos(Y[5]),sPsi6=sin(Y[5]), cPsi7=cos(Y[6]),sPsi7=sin(Y[6]), cPsi8=cos(Y[7]),sPsi8=sin(Y[7]);

         j[0][0]=-sPsi1;

         j[1][0]=cPsi1*cPsi2; j[1][1]=-sPsi1*sPsi2;

         j[2][0]=cPsi1*sPsi2*cPsi3; j[2][1]=sPsi1*cPsi2*cPsi3; j[2][2]=-sPsi1*sPsi2*sPsi3; 

         j[3][0]=cPsi1*sPsi2*sPsi3*cPsi4; j[3][1]=sPsi1*cPsi2*sPsi3*cPsi4; j[3][2]=sPsi1*sPsi2*cPsi3*cPsi4; j[3][3]=-sPsi1*sPsi2*sPsi3*sPsi4;

         j[4][0]=cPsi1*sPsi2*sPsi3*sPsi4*cPsi5; j[4][1]=sPsi1*cPsi2*sPsi3*sPsi4*cPsi5;
         j[4][2]=sPsi1*sPsi2*cPsi3*sPsi4*cPsi5; j[4][3]=sPsi1*sPsi2*sPsi3*cPsi4*cPsi5; j[4][4]=-sPsi1*sPsi2*sPsi3*sPsi4*sPsi5;

         j[5][0]=cPsi1*sPsi2*sPsi3*sPsi4*sPsi5; j[5][1]=sPsi1*cPsi2*sPsi3*sPsi4*sPsi5;
         j[5][2]=sPsi1*sPsi2*cPsi3*sPsi4*sPsi5; j[5][3]=sPsi1*sPsi2*sPsi3*cPsi4*sPsi5; j[5][4]=sPsi1*sPsi2*sPsi3*sPsi4*cPsi5;

         j[6][0]=cPsi1*sPsi2*cPsi6; j[6][1]=sPsi1*cPsi2*cPsi6; j[6][5]=-sPsi1*sPsi2*sPsi6;

         j[7][0]=cPsi1*sPsi2*sPsi6*cPsi7; j[7][1]=sPsi1*cPsi2*sPsi6*cPsi7; j[7][5]=sPsi1*sPsi2*cPsi6*cPsi7; j[7][6]=-sPsi1*sPsi2*sPsi6*sPsi7;

         j[8][0]=cPsi1*sPsi2*sPsi6*sPsi7*cPsi8; j[8][1]=sPsi1*cPsi2*sPsi6*sPsi7*cPsi8;
         j[8][5]=sPsi1*sPsi2*cPsi6*sPsi7*cPsi8; j[8][6]=sPsi1*sPsi2*sPsi6*cPsi7*cPsi8; j[8][7]=-sPsi1*sPsi2*sPsi6*sPsi7*sPsi8;

         j[9][0]=cPsi1*sPsi2*sPsi6*sPsi7*sPsi8; j[9][1]=sPsi1*cPsi2*sPsi6*sPsi7*sPsi8;
         j[9][5]=sPsi1*sPsi2*cPsi6*sPsi7*sPsi8; j[9][6]=sPsi1*sPsi2*sPsi6*cPsi7*sPsi8; j[9][7]=sPsi1*sPsi2*sPsi6*sPsi7*cPsi8;

         return j;
        }

void J(array<double,NY> Y,MATRIX<double,10,8> &j)
       { double cPsi1=cos(Y[0]),sPsi1=sin(Y[0]), cPsi2=cos(Y[1]),sPsi2=sin(Y[1]);
         double cPsi3=cos(Y[2]),sPsi3=sin(Y[2]), cPsi4=cos(Y[3]),sPsi4=sin(Y[3]), cPsi5=cos(Y[4]),sPsi5=sin(Y[4]);
         double cPsi6=cos(Y[5]),sPsi6=sin(Y[5]), cPsi7=cos(Y[6]),sPsi7=sin(Y[6]), cPsi8=cos(Y[7]),sPsi8=sin(Y[7]);

         j[0][0]=-sPsi1;

         j[1][0]=cPsi1*cPsi2; j[1][1]=-sPsi1*sPsi2;

         j[2][0]=cPsi1*sPsi2*cPsi3; j[2][1]=sPsi1*cPsi2*cPsi3; j[2][2]=-sPsi1*sPsi2*sPsi3; 

         j[3][0]=cPsi1*sPsi2*sPsi3*cPsi4; j[3][1]=sPsi1*cPsi2*sPsi3*cPsi4; j[3][2]=sPsi1*sPsi2*cPsi3*cPsi4; j[3][3]=-sPsi1*sPsi2*sPsi3*sPsi4;

         j[4][0]=cPsi1*sPsi2*sPsi3*sPsi4*cPsi5; j[4][1]=sPsi1*cPsi2*sPsi3*sPsi4*cPsi5;
         j[4][2]=sPsi1*sPsi2*cPsi3*sPsi4*cPsi5; j[4][3]=sPsi1*sPsi2*sPsi3*cPsi4*cPsi5; j[4][4]=-sPsi1*sPsi2*sPsi3*sPsi4*sPsi5;

         j[5][0]=cPsi1*sPsi2*sPsi3*sPsi4*sPsi5; j[5][1]=sPsi1*cPsi2*sPsi3*sPsi4*sPsi5;
         j[5][2]=sPsi1*sPsi2*cPsi3*sPsi4*sPsi5; j[5][3]=sPsi1*sPsi2*sPsi3*cPsi4*sPsi5; j[5][4]=sPsi1*sPsi2*sPsi3*sPsi4*cPsi5;

         j[6][0]=cPsi1*sPsi2*cPsi6; j[6][1]=sPsi1*cPsi2*cPsi6; j[6][5]=-sPsi1*sPsi2*sPsi6;

         j[7][0]=cPsi1*sPsi2*sPsi6*cPsi7; j[7][1]=sPsi1*cPsi2*sPsi6*cPsi7; j[7][5]=sPsi1*sPsi2*cPsi6*cPsi7; j[7][6]=-sPsi1*sPsi2*sPsi6*sPsi7;

         j[8][0]=cPsi1*sPsi2*sPsi6*sPsi7*cPsi8; j[8][1]=sPsi1*cPsi2*sPsi6*sPsi7*cPsi8;
         j[8][5]=sPsi1*sPsi2*cPsi6*sPsi7*cPsi8; j[8][6]=sPsi1*sPsi2*sPsi6*cPsi7*cPsi8; j[8][7]=-sPsi1*sPsi2*sPsi6*sPsi7*sPsi8;

         j[9][0]=cPsi1*sPsi2*sPsi6*sPsi7*sPsi8; j[9][1]=sPsi1*cPsi2*sPsi6*sPsi7*sPsi8;
         j[9][5]=sPsi1*sPsi2*cPsi6*sPsi7*sPsi8; j[9][6]=sPsi1*sPsi2*sPsi6*cPsi7*sPsi8; j[9][7]=sPsi1*sPsi2*sPsi6*sPsi7*cPsi8;
        }

MATRIX<double,8,10> JInverse(array<double,NY> Y)
       { MATRIX<double,8,10> j;

         double cPsi1=cos(Y[0]),sPsi1=sin(Y[0]), cPsi2=cos(Y[1]),sPsi2=sin(Y[1]);
         double cPsi3=cos(Y[2]),sPsi3=sin(Y[2]), cPsi4=cos(Y[3]),sPsi4=sin(Y[3]), cPsi5=cos(Y[4]),sPsi5=sin(Y[4]);
         double cPsi6=cos(Y[5]),sPsi6=sin(Y[5]), cPsi7=cos(Y[6]),sPsi7=sin(Y[6]), cPsi8=cos(Y[7]),sPsi8=sin(Y[7]);

         j[0][0]=-sPsi1; j[0][1]=cPsi1*cPsi2; j[0][2]=cPsi1*sPsi2*cPsi3; j[0][3]=cPsi1*sPsi2*sPsi3*cPsi4; j[0][4]=cPsi1*sPsi2*sPsi3*sPsi4*cPsi5; j[0][5]=cPsi1*sPsi2*sPsi3*sPsi4*sPsi5;

         j[1][1]=-sPsi2/sPsi1; j[1][2]=cPsi2*cPsi3/sPsi1; j[1][3]=cPsi2*sPsi3*cPsi4/sPsi1; j[1][4]=cPsi2*sPsi3*sPsi4*cPsi5/sPsi1; j[1][5]=cPsi2*sPsi3*sPsi4*sPsi5/sPsi1;

         sPsi1*=sPsi2;
         j[2][2]=-sPsi3/sPsi1;           j[2][3]=cPsi3*cPsi4/sPsi1;       j[2][4]=cPsi3*sPsi4*cPsi5/sPsi1; j[2][5]=cPsi3*sPsi4*sPsi5/sPsi1;
         j[3][3]=-sPsi4/sPsi1/sPsi3;     j[3][4]=cPsi4*cPsi5/sPsi1/sPsi3; j[3][5]=cPsi4*sPsi5/sPsi1/sPsi3;
         sPsi3*=sPsi4;
         j[4][4]=-sPsi5/sPsi1/sPsi3;     j[4][5]=cPsi5/sPsi1/sPsi3;

         j[5][6]=-sPsi6/sPsi1;           j[5][7]=cPsi6*cPsi7/sPsi1; j[5][8]=cPsi6*sPsi7*cPsi8/sPsi1;    j[5][9]=cPsi6*sPsi7*sPsi8/sPsi1;
         sPsi1*=sPsi6;
         j[6][7]=-sPsi7/sPsi1;           j[6][8]=cPsi7*cPsi8/sPsi1; j[6][9]=cPsi7*sPsi8/sPsi1;
         sPsi1*=sPsi7;
         j[7][8]=-sPsi8/sPsi1;           j[7][9]=cPsi8/sPsi1;

         return j;
        }

void JInverse(array<double,NY> Y,MATRIX<double,8,10> &j)
       { double cPsi1=cos(Y[0]),sPsi1=sin(Y[0]), cPsi2=cos(Y[1]),sPsi2=sin(Y[1]);
         double cPsi3=cos(Y[2]),sPsi3=sin(Y[2]), cPsi4=cos(Y[3]),sPsi4=sin(Y[3]), cPsi5=cos(Y[4]),sPsi5=sin(Y[4]);
         double cPsi6=cos(Y[5]),sPsi6=sin(Y[5]), cPsi7=cos(Y[6]),sPsi7=sin(Y[6]), cPsi8=cos(Y[7]),sPsi8=sin(Y[7]);

         j[0][0]=-sPsi1; j[0][1]=cPsi1*cPsi2; j[0][2]=cPsi1*sPsi2*cPsi3; j[0][3]=cPsi1*sPsi2*sPsi3*cPsi4; j[0][4]=cPsi1*sPsi2*sPsi3*sPsi4*cPsi5; j[0][5]=cPsi1*sPsi2*sPsi3*sPsi4*sPsi5;

         j[1][1]=-sPsi2/sPsi1; j[1][2]=cPsi2*cPsi3/sPsi1; j[1][3]=cPsi2*sPsi3*cPsi4/sPsi1; j[1][4]=cPsi2*sPsi3*sPsi4*cPsi5/sPsi1; j[1][5]=cPsi2*sPsi3*sPsi4*sPsi5/sPsi1;

         sPsi1*=sPsi2;
         j[2][2]=-sPsi3/sPsi1;           j[2][3]=cPsi3*cPsi4/sPsi1;       j[2][4]=cPsi3*sPsi4*cPsi5/sPsi1; j[2][5]=cPsi3*sPsi4*sPsi5/sPsi1;
         j[3][3]=-sPsi4/sPsi1/sPsi3;     j[3][4]=cPsi4*cPsi5/sPsi1/sPsi3; j[3][5]=cPsi4*sPsi5/sPsi1/sPsi3;
         sPsi3*=sPsi4;
         j[4][4]=-sPsi5/sPsi1/sPsi3;     j[4][5]=cPsi5/sPsi1/sPsi3;

         j[5][6]=-sPsi6/sPsi1;           j[5][7]=cPsi6*cPsi7/sPsi1; j[5][8]=cPsi6*sPsi7*cPsi8/sPsi1;    j[5][9]=cPsi6*sPsi7*sPsi8/sPsi1;
         sPsi1*=sPsi6;
         j[6][7]=-sPsi7/sPsi1;           j[6][8]=cPsi7*cPsi8/sPsi1; j[6][9]=cPsi7*sPsi8/sPsi1;
         sPsi1*=sPsi7;
         j[7][8]=-sPsi8/sPsi1;           j[7][9]=cPsi8/sPsi1;
        }

