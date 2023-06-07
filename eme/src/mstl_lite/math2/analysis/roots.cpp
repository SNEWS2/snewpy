#include "roots.h"

using std::vector;
using std::abs;

// ************************************************************************************************
// ******************************* Newton's method ************************************************
// ************************************************************************************************

//two and higher dimensions

//each function must take the same arguments even if all those areguments are not used
vector<double> NewtonRoots(vector<double> x0,vector<double(*)(vector<double>)> functions,vector<vector<double(*)(vector<double>)> > Jacobian,double accuracy)
      { int NX=static_cast<int>(x0.size()), NF=static_cast<int>(functions.size());
        if(NX!=NF){ throw DIFFERENT_LENGTHS("NewtonRoot");}

        vector<double> x(NX), deltax(NX);
        DCVECTOR DX(NX); 
        DRVECTOR DXT;

        DCVECTOR F(NF), deltaF(NF);
        DMATRIX J(NF,NF);

        double maxerror=0.;
        for(int i=0;i<=NX-1;i++){ x[i]=x0[i];}

        for(int i=0;i<=NF-1;i++)
           { F[i]=(*functions[i])(x); 
             if(abs(F[i])>maxerror){ maxerror=abs(F[i]);} 
             for(int j=0;j<=NF-1;j++){ J[i][j]=(*Jacobian[i][j])(x);} 
            }
        if(maxerror<accuracy){ return x;} // early return

        int counter=0,countermax=5;
        double lambda;
        bool illconditionedJ=false,localminimum=false, repeat;
        do{ if(illconditionedJ==true || localminimum==true)
              { J=ZeroMatrix<double>(NF);
                for(int i=0;i<=NF-1;i++){ J[i][i]=1.;}
                illconditionedJ=localminimum=false;
               }
            DX=-LUInverse(J)*F;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}

            try{ lambda=LineSearch(x,deltax,functions);}
            catch(NO_SOLUTION &NS){ lambda=0.;}
            DX*=lambda;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}
            DXT=Transpose(DX);

            if(Equality(lambda,0.)==false){ // update x, F, and J
                            for(int i=0;i<=NX-1;i++){ x[i]+=deltax[i];}
                            for(int i=0;i<=NF-1;i++){ deltaF[i]=(*functions[i])(x)-F[i];}
                            F+=deltaF; 
                            for(int i=0;i<=NF-1;i++){ for(int j=0;j<=NF-1;j++){ J[i][j]=(*Jacobian[i][j])(x);} }
                            for(int i=0;i<=NX-1;i++){ deltax[i]=abs(x[i]-deltax[i])-abs(x[i]);}
                           }

            maxerror=0.;
            for(int i=0;i<=NF-1;i++){ if(abs(F[i])>maxerror){ maxerror=abs(F[i]);} }

            if(maxerror>accuracy){ repeat=true; if(Equality(DXT*DX,0.)==true){ counter++; localminimum=true;} }
            else{ repeat=false;}

            if(counter>=countermax){ repeat=false;}
           }while(repeat==true);

        if(counter>=countermax){ throw NO_SOLUTION("NewtonRoots");}

        return x;
       }

//each function must take the same arguments even if all those areguments are not used
vector<double> NewtonRoots(vector<double> x0,vector<double>(*functions)(vector<double>),vector<vector<double> >(*Jacobian)(vector<double>),double accuracy)
      { int NX=static_cast<int>(x0.size()), NF;

        vector<double> x(NX), deltax(NX);

        vector<double> ff;
        vector<vector<double> > JJ;

        DCVECTOR DX(NX); 
        DRVECTOR DXT;

        double maxerror=0.;
        for(int i=0;i<=NX-1;i++){ x[i]=x0[i];}

        ff=(*functions)(x);
        NF=static_cast<int>(ff.size()); 
        JJ=Jacobian(x);     

        DCVECTOR F(NF), deltaF(NF);
        DMATRIX J(NF,NF);

        for(int i=0;i<=NF-1;i++)
           { F[i]=ff[i]; 
             if(abs(F[i])>maxerror){ maxerror=abs(F[i]);} 
             for(int j=0;j<=NF-1;j++){ J[i][j]=JJ[i][j];} 
            }
        if(maxerror<accuracy){ return x;} // early return

        int counter=0,countermax=5;
        double lambda;
        bool illconditionedJ=false,localminimum=false, repeat;
        do{ if(illconditionedJ==true || localminimum==true)
              { J=ZeroMatrix<double>(NF);
                for(int i=0;i<=NF-1;i++){ J[i][i]=1.;}
                illconditionedJ=localminimum=false;
               }
            DX=-LUInverse(J)*F;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}

            try{ lambda=LineSearch(x,deltax,functions);}
            catch(NO_SOLUTION &NS){ lambda=0.;}
            DX*=lambda;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}
            DXT=Transpose(DX);

            if(Equality(lambda,0.)==false){ // update x, F, and J
                            for(int i=0;i<=NX-1;i++){ x[i]+=deltax[i];}
                            ff=(*functions)(x); for(int i=0;i<=NF-1;i++){ deltaF[i]=ff[i]-F[i];}
                            F+=deltaF; 
                            JJ=(*Jacobian)(x);  for(int i=0;i<=NF-1;i++){ for(int j=0;j<=NF-1;j++){ J[i][j]=JJ[i][j];} }
                            for(int i=0;i<=NX-1;i++){ deltax[i]=abs(x[i]-deltax[i])-abs(x[i]);}
                           }

            maxerror=0.;
            for(int i=0;i<=NF-1;i++){ if(abs(F[i])>maxerror){ maxerror=abs(F[i]);} }

            if(maxerror>accuracy){ repeat=true; if(Equality(DXT*DX,0.)==true){ counter++; localminimum=true;} }
            else{ repeat=false;}

            if(counter>=countermax){ repeat=false;}
           }while(repeat==true);

        if(counter>=countermax){ throw NO_SOLUTION("NewtonRoots");}

        return x;
       }

// ******************************* Broyden's method ************************************************

vector<double> BroydenRoots(double x01,double x02,double(*f1)(vector<double>),double(*f2)(vector<double>),double accuracy)
       { vector<double> x0s(2);                    x0s[0]=x01; x0s[1]=x02;
         vector<double(*)(vector<double>)> fs(2);  fs[0]=f1;   fs[1]=f2;
         return BroydenRoots(x0s,fs,accuracy);
        }

vector<double> BroydenRoots(vector<double> x0,vector<double(*)(vector<double>)> functions,double accuracy)
      { int NX=static_cast<int>(x0.size()), NF=static_cast<int>(functions.size());
        if(NX!=NF){ throw DIFFERENT_LENGTHS("BroydenRoot");}

        vector<double> x(NX), deltax(NX);
        DCVECTOR DX(NX); DRVECTOR DXT;

        DCVECTOR F(NF), deltaF(NF); 
        DMATRIX J(NF,NF);

        double maxerror=0., maxchange=0.;

        x=x0;
        for(int i=0;i<=NF-1;i++)
           { F[i]=(*functions[i])(x);
             if(abs(F[i])>maxerror){ maxerror=abs(F[i]);}
             J[i][i]=1.;
            }
        if(maxerror<accuracy){ return x;} // early return

        int counter=0,countermax=5;
        double lambda;
        bool illconditionedJ=false, localminimum=false, repeat;

        do{ if(illconditionedJ==true || localminimum==true)
              { J=ZeroMatrix<double>(NF);
                for(int i=0;i<=NF-1;i++){ J[i][i]=1.;}
                illconditionedJ=localminimum=false;
               }

            DX=-LUInverse(J)*F;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];} 

            try{ lambda=LineSearch(x,deltax,functions);}
            catch(NO_SOLUTION &NS){ lambda=0.;}
            DX*=lambda;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}
            DXT=Transpose(DX);

            if(Equality(lambda,0.)==false)
	      { for(int i=0;i<=NX-1;i++){ x[i]+=deltax[i];}
                for(int i=0;i<=NF-1;i++){ deltaF[i]=(*functions[i])(x)-F[i];}

                F+=deltaF; 
                J+=DMATRIX( (deltaF-J*DX)*DXT/(DXT*DX) );

                if(abs(Determinant(J))<1e-3){ illconditionedJ=true;}

                for(int i=0;i<=NX-1;i++){ deltax[i]=x[i]-(x[i]-deltax[i]);}
                if(Equality(DXT*DX,0.)==true){ counter++; localminimum=true;}

                //stopping criteria
                maxerror=maxchange=0.;
                for(int i=0;i<=NF-1;i++)
                   { if(abs(F[i])>maxerror){ maxerror=abs(F[i]);}
                     if(abs(deltaF[i])>maxchange){ maxchange=abs(deltaF[i]);}
                    }
                if(maxerror>accuracy && maxchange>0. && counter<countermax){ repeat=true;}
                else{ repeat=false;}
               }
            else{ if(counter<countermax){ counter++; localminimum=true; repeat=true;}
                  else{ repeat=false;}
                 }
           }while(repeat==true);

        if(counter>=countermax || (maxerror>accuracy && Equality(maxchange,0.)==true))
          { throw NO_SOLUTION("BroydenRoots");}

        return x;
       }

// *****************************************************************************************************

vector<double> BroydenRoots(vector<double> x0,vector<double(*)(vector<double>)> functions,vector<vector<double(*)(vector<double>)> > Jacobian,double accuracy)
      { int NX=static_cast<int>(x0.size()), NF=static_cast<int>(functions.size());
        if(NX!=NF){ throw DIFFERENT_LENGTHS("BroydenRoot");}

        vector<double> x(NX), deltax(NX);
        DCVECTOR DX(NX); DRVECTOR DXT;

        DCVECTOR F(NF), deltaF(NF); 
        DMATRIX J(NF,NF);

        double maxerror=0., maxchange=0.;

        x=x0;
        for(int i=0;i<=NF-1;i++)
           { F[i]=(*functions[i])(x);
             if(abs(F[i])>maxerror){ maxerror=abs(F[i]);}
             for(int j=0;j<=NF-1;j++){ J[i][j]=Jacobian[i][j](x);}
            }
        if(maxerror<accuracy){ return x;} // early return

        int counter=0,countermax=5;
        double lambda;
        bool illconditionedJ=false, localminimum=false, repeat;

        do{ if(illconditionedJ==true || localminimum==true)
              { J=ZeroMatrix<double>(NF);
                for(int i=0;i<=NF-1;i++){ J[i][i]=1.;}
                illconditionedJ=localminimum=false;
               }

            DX=-LUInverse(J)*F;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];} 

            try{ lambda=LineSearch(x,deltax,functions);}
            catch(NO_SOLUTION &NS){ lambda=0.;}
            DX*=lambda;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}
            DXT=Transpose(DX);

            if(Equality(lambda,0.)==false)
	      { for(int i=0;i<=NX-1;i++){ x[i]+=deltax[i];}
                for(int i=0;i<=NF-1;i++){ deltaF[i]=(*functions[i])(x)-F[i];}

                F+=deltaF; 
                J+=DMATRIX( (deltaF-J*DX)*DXT/(DXT*DX) );

                if(abs(Determinant(J))<1e-3){ illconditionedJ=true;}

                for(int i=0;i<=NX-1;i++){ deltax[i]=x[i]-(x[i]-deltax[i]);}
                if(Equality(DXT*DX,0.)==true){ counter++; localminimum=true;}

                //stopping criteria
                maxerror=maxchange=0.;
                for(int i=0;i<=NF-1;i++)
                   { if(abs(F[i])>maxerror){ maxerror=abs(F[i]);}
                     if(abs(deltaF[i])>maxchange){ maxchange=abs(deltaF[i]);}
                    }
                if(maxerror>accuracy && maxchange>0. && counter<countermax){ repeat=true;}
                else{ repeat=false;}
               }
            else{ if(counter<countermax){ counter++; localminimum=true; repeat=true;}
                  else{ repeat=false;}
                 }
           }while(repeat==true);

        if(counter>=countermax || (maxerror>accuracy && Equality(maxchange,0.)==true))
          { throw NO_SOLUTION("BroydenRoots");}

        return x;
       }

// *****************************************************************************************************

vector<double> BroydenRoots(vector<double> x0,vector<double>(*functions)(vector<double>),vector<vector<double> >(*Jacobian)(std::vector<double>),double accuracy)
      { int NX=static_cast<int>(x0.size()), NF;

        vector<double> x(NX), deltax(NX);
        DCVECTOR DX(x0.size()); DRVECTOR DXT;

        vector<double> ff;
        vector<vector<double> > JJ;

        double maxerror=0., maxchange=0.;

        x=x0;
        ff=(*functions)(x); 
        JJ=(*Jacobian)(x);
        NF=static_cast<int>(ff.size());

        DCVECTOR F(NF), deltaF(NF); 
        DMATRIX J(NF,NF);

        if(NX!=NF){ throw DIFFERENT_LENGTHS("BroydenRoot");}

        for(int i=0;i<=NF-1;i++)
           { F[i]=ff[i];
             if(abs(F[i])>maxerror){ maxerror=abs(F[i]);}
             for(int j=0;j<=NF-1;j++){ J[i][j]=JJ[i][j];}
            }
        if(maxerror<accuracy){ return x;} // early return

        int counter=0,countermax=5;
        double lambda;
        bool illconditionedJ=false, localminimum=false, repeat;

        do{ if(illconditionedJ==true || localminimum==true)
              { J=ZeroMatrix<double>(NF);
                for(int i=0;i<=NF-1;i++){ J[i][i]=1.;}
                illconditionedJ=localminimum=false;
               }

            DX=-LUInverse(J)*F;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];} 

            try{ lambda=LineSearch(x,deltax,functions);}
            catch(NO_SOLUTION &NS){ lambda=0.;} 
            DX*=lambda;
            for(int i=0;i<=NX-1;i++){ deltax[i]=DX[i];}
            DXT=Transpose(DX);

            if(Equality(lambda,0.)==false)
	      { for(int i=0;i<=NX-1;i++){ x[i]+=deltax[i];}
                ff=(*functions)(x); for(int i=0;i<=NF-1;i++){ deltaF[i]=ff[i]-F[i];}

                F+=deltaF; 
                J+=DMATRIX( (deltaF-J*DX)*DXT/(DXT*DX) );

                if(abs(Determinant(J))<1e-3){ illconditionedJ=true;}

                for(int i=0;i<=NX-1;i++){ deltax[i]=x[i]-(x[i]-deltax[i]);}
                if(Equality(DXT*DX,0.)==true){ counter++; localminimum=true;}

                //stopping criteria
                maxerror=maxchange=0.;
                for(int i=0;i<=NF-1;i++)
                   { if(abs(F[i])>maxerror){ maxerror=abs(F[i]);}
                     if(abs(deltaF[i])>maxchange){ maxchange=abs(deltaF[i]);}
                    }
                if(maxerror>accuracy && maxchange>0. && counter<countermax){ repeat=true;}
                else{ repeat=false;}
               }
            else{ if(counter<countermax){ counter++; localminimum=true; repeat=true;}
                  else{ repeat=false;}
                 }
           }while(repeat==true);

        if(counter>=countermax || (maxerror>accuracy && Equality(maxchange,0.)==true))
          { throw NO_SOLUTION("BroydenRoots");}

        return x;
       }
