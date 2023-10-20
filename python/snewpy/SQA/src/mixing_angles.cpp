
#include "mixing_angles.h"

// *********************************************************************

using std::complex;
using std::array;

// ********************************************************************************************************

MATRIX<complex<double>,NF,NF> MixingMatrix(MATRIX<complex<double>,NF,NF> Hf,array<double,NF> k,array<double,NF> dk,array<array<double,NF>,NF> A);

// *********************************************************************
// *********************************************************************
// *********************************************************************

MATRIX<complex<double>,NF,NF> MixingMatrix(MATRIX<complex<double>,NF,NF> Hf,array<double,NF> k,array<double,NF> dk,array<array<double,NF>,NF> A)
       { MATRIX<complex<double>,NF,NF> u;

         double d;
         array<double,NF> r2;

         for(int j=0;j<=NF-1;j++)
            { if(j==0){ d=dk[0]*dk[1];}  // first column 
              if(j==1){ d=-dk[0]*dk[2];} // second column
              if(j==2){ d=dk[1]*dk[2];}  // third column

              r2[e]=real(C<e,e>(Hf,k[j]))*d;
              r2[mu]=real(C<mu,mu>(Hf,k[j]))*d;
              r2[tau]=real(C<tau,tau>(Hf,k[j]))*d;

              if(r2[e]>=r2[mu] && r2[e]>=r2[tau]){ u[e][j]=A[j][e]*C<e,e>(Hf,k[j])/sqrt(r2[e]);         u[mu][j]=A[j][e]*C<e,mu>(Hf,k[j])/sqrt(r2[e]);       u[tau][j]=A[j][e]*C<e,tau>(Hf,k[j])/sqrt(r2[e]);}
              if(r2[mu]>=r2[e] && r2[mu]>=r2[tau]){ u[e][j]=A[j][mu]*C<mu,e>(Hf,k[j])/sqrt(r2[mu]);     u[mu][j]=A[j][mu]*C<mu,mu>(Hf,k[j])/sqrt(r2[mu]);    u[tau][j]=A[j][mu]*C<mu,tau>(Hf,k[j])/sqrt(r2[mu]);}
              if(r2[tau]>=r2[e] && r2[tau]>=r2[mu]){ u[e][j]=A[j][tau]*C<tau,e>(Hf,k[j])/sqrt(r2[tau]); u[mu][j]=A[j][tau]*C<tau,mu>(Hf,k[j])/sqrt(r2[tau]); u[tau][j]=A[j][tau]*C<tau,tau>(Hf,k[j])/sqrt(r2[tau]);}
              // set the element in the e (top) row to be pure real for the first two columns
              //if(j==0 || j==1){ u[e][j]=A[j][e]*C<e,e>(Hf,k[j])/sqrt(r2[e]); u[mu][j]=A[j][e]*C<e,mu>(Hf,k[j])/sqrt(r2[e]); u[tau][j]=A[j][e]*C<e,tau>(Hf,k[j])/sqrt(r2[e]);}
              // set the element in the tau (bottom) row to be pure real for the 3rd column
              //if(j==2){ u[e][j]=A[j][tau]*C<tau,e>(Hf,k[j])/sqrt(r2[tau]); u[mu][j]=A[j][tau]*C<tau,mu>(Hf,k[j])/sqrt(r2[tau]); u[tau][j]=A[j][tau]*C<tau,tau>(Hf,k[j])/sqrt(r2[tau]);}
            }
       
        return u;
       }

MATRIX<complex<double>,NF,NF> MixingMatrix(array<double,NF> dk,array<MATRIX<complex<double>,NF,NF>,NF> &C,array<array<double,NF>,NF> A)
       { MATRIX<complex<double>,NF,NF> u;

         double d;

         array<double,NF> r2;

         for(int j=0;j<=NF-1;j++)
            { if(j==0){ d=dk[0]*dk[1];}  // first column 
              if(j==1){ d=-dk[0]*dk[2];} // second column
              if(j==2){ d=dk[1]*dk[2];}  // third column

              r2[e]=real(C[j][e][e])*d;
              r2[mu]=real(C[j][mu][mu])*d;
              r2[tau]=real(C[j][tau][tau])*d;

              if(r2[e]>=r2[mu] && r2[e]>=r2[tau]){ u[e][j]=A[j][e]*C[j][e][e]/sqrt(r2[e]);         u[mu][j]=A[j][e]*C[j][e][mu]/sqrt(r2[e]);       u[tau][j]=A[j][e]*C[j][e][tau]/sqrt(r2[e]);}
              if(r2[mu]>=r2[e] && r2[mu]>=r2[tau]){ u[e][j]=A[j][mu]*C[j][mu][e]/sqrt(r2[mu]);     u[mu][j]=A[j][mu]*C[j][mu][mu]/sqrt(r2[mu]);    u[tau][j]=A[j][mu]*C[j][mu][tau]/sqrt(r2[mu]);}
              if(r2[tau]>=r2[e] && r2[tau]>=r2[mu]){ u[e][j]=A[j][tau]*C[j][tau][e]/sqrt(r2[tau]); u[mu][j]=A[j][tau]*C[j][tau][mu]/sqrt(r2[tau]); u[tau][j]=A[j][tau]*C[j][tau][tau]/sqrt(r2[tau]);}
              // set the element in the e (top) row to be pure real for the first two columns
              //if(j==0 || j==1){ u[e][j]=A[j][e]*C[j][e][e]/sqrt(r2[e]); u[mu][j]=A[j][e]*C[j][e][mu]/sqrt(r2[e]); u[tau][j]=A[j][e]*C[j][e][tau]/sqrt(r2[e]);}
              // set the element in the tau (bottom) row to be pure real for the 3rd column
              //if(j==2){ u[e][j]=A[j][tau]*C[j][tau][e]/sqrt(r2[tau]); u[mu][j]=A[j][tau]*C[j][tau][mu]/sqrt(r2[tau]); u[tau][j]=A[j][tau]*C[j][tau][tau]/sqrt(r2[tau]);}
            }
       
        return u;
       }

void Evaluate_UV(void) 
       { UV[nu][0][0] = c12V*c13V*exp(I*etaV[0]);
         UV[nu][0][1] = s12V*c13V*exp(I*etaV[1]);
         UV[nu][0][2] = s13V*exp(-I*deltaV);

         UV[nu][1][0] = -(s12V*c23V + c12V*s13V*s23V*exp(I*deltaV)) * exp(I*etaV[0]);
         UV[nu][1][1] =  (c12V*c23V - s12V*s13V*s23V*exp(I*deltaV)) * exp(I*etaV[1]);
         UV[nu][1][2] =  c13V*s23V;

         UV[nu][2][0] =  (s12V*s23V - c12V*s13V*c23V*exp(I*deltaV)) * exp(I*etaV[0]);
         UV[nu][2][1] = -(c12V*s23V + s12V*s13V*c23V*exp(I*deltaV)) * exp(I*etaV[1]);
         UV[nu][2][2] =  c13V*c23V;

         UV[antinu]=Conjugate(UV[nu]);
        }

// ********************************************************************************

void Evaluate_CV(void)
     { for(int i=0;i<=NE-1;i++)
          { CV[nu][i][0][e][mu]=C<e,mu>(HfV[nu][i],kV[i][0]);   CV[nu][i][1][e][mu]=C<e,mu>(HfV[nu][i],kV[i][1]);   CV[nu][i][2][e][mu]=C<e,mu>(HfV[nu][i],kV[i][2]);
            CV[nu][i][0][e][tau]=C<e,tau>(HfV[nu][i],kV[i][0]); CV[nu][i][1][e][tau]=C<e,tau>(HfV[nu][i],kV[i][1]); CV[nu][i][2][e][tau]=C<e,tau>(HfV[nu][i],kV[i][2]);

            CV[nu][i][0][mu][e]=C<mu,e>(HfV[nu][i],kV[i][0]);     CV[nu][i][1][mu][e]=C<mu,e>(HfV[nu][i],kV[i][1]);     CV[nu][i][2][mu][e]=C<mu,e>(HfV[nu][i],kV[i][2]);
            CV[nu][i][0][mu][tau]=C<mu,tau>(HfV[nu][i],kV[i][0]); CV[nu][i][1][mu][tau]=C<mu,tau>(HfV[nu][i],kV[i][1]); CV[nu][i][2][mu][tau]=C<mu,tau>(HfV[nu][i],kV[i][2]);

            CV[nu][i][0][tau][e]=C<tau,e>(HfV[nu][i],kV[i][0]);   CV[nu][i][1][tau][e]=C<tau,e>(HfV[nu][i],kV[i][1]);   CV[nu][i][2][tau][e]=C<tau,e>(HfV[nu][i],kV[i][2]);
            CV[nu][i][0][tau][mu]=C<tau,mu>(HfV[nu][i],kV[i][0]); CV[nu][i][1][tau][mu]=C<tau,mu>(HfV[nu][i],kV[i][1]); CV[nu][i][2][tau][mu]=C<tau,mu>(HfV[nu][i],kV[i][2]);

            CV[antinu][i][0][e][mu]=C<e,mu>(HfV[antinu][i],kV[i][0]);   CV[antinu][i][1][e][mu]=C<e,mu>(HfV[antinu][i],kV[i][1]);   CV[antinu][i][2][e][mu]=C<e,mu>(HfV[antinu][i],kV[i][2]);
            CV[antinu][i][0][e][tau]=C<e,tau>(HfV[antinu][i],kV[i][0]); CV[antinu][i][1][e][tau]=C<e,tau>(HfV[antinu][i],kV[i][1]); CV[antinu][i][2][e][tau]=C<e,tau>(HfV[antinu][i],kV[i][2]);

            CV[antinu][i][0][mu][e]=C<mu,e>(HfV[antinu][i],kV[i][0]);     CV[antinu][i][1][mu][e]=C<mu,e>(HfV[antinu][i],kV[i][1]);     CV[antinu][i][2][mu][e]=C<mu,e>(HfV[antinu][i],kV[i][2]);
            CV[antinu][i][0][mu][tau]=C<mu,tau>(HfV[antinu][i],kV[i][0]); CV[antinu][i][1][mu][tau]=C<mu,tau>(HfV[antinu][i],kV[i][1]); CV[antinu][i][2][mu][tau]=C<mu,tau>(HfV[antinu][i],kV[i][2]);

            CV[antinu][i][0][tau][e]=C<tau,e>(HfV[antinu][i],kV[i][0]);   CV[antinu][i][1][tau][e]=C<tau,e>(HfV[antinu][i],kV[i][1]);   CV[antinu][i][2][tau][e]=C<tau,e>(HfV[antinu][i],kV[i][2]);
            CV[antinu][i][0][tau][mu]=C<tau,mu>(HfV[antinu][i],kV[i][0]); CV[antinu][i][1][tau][mu]=C<tau,mu>(HfV[antinu][i],kV[i][1]); CV[antinu][i][2][tau][mu]=C<tau,mu>(HfV[antinu][i],kV[i][2]);
           } 
      }

// ********************************************************************************

void Evaluate_AV(void)
     { double Delta;
       for(int i=0;i<=NE-1;i++)
          { for(int j=0;j<=NF-1;j++)
               { if(j==0){ Delta=(kV[i][1]-kV[i][0])*(kV[i][2]-kV[i][0]);}
                 if(j==1){ Delta=(kV[i][0]-kV[i][1])*(kV[i][2]-kV[i][1]);} 
                 if(j==2){ Delta=(kV[i][0]-kV[i][2])*(kV[i][1]-kV[i][2]);}  

                 double re2=Delta*real(C<e,e>(HfV[nu][i],kV[i][j]));
                 double rmu2=Delta*real(C<mu,mu>(HfV[nu][i],kV[i][j]));
                 double rtau2=Delta*real(C<tau,tau>(HfV[nu][i],kV[i][j]));

                 if(norm(UV[nu][e][j])>norm(UV[nu][mu][j]) && norm(UV[nu][e][j])>norm(UV[nu][tau][j]) )
                   { AV[nu][i][j][e]=real( UV[nu][e][j]*sqrt(re2) / C<e,e>(HfV[nu][i],kV[i][j]) );
                     AV[nu][i][j][mu]=real( UV[nu][e][j]*sqrt(rmu2) / C<mu,e>(HfV[nu][i],kV[i][j]) );
                     AV[nu][i][j][tau]=real( UV[nu][e][j]*sqrt(rtau2) / C<tau,e>(HfV[nu][i],kV[i][j]) );
                    }
                 if(norm(UV[nu][mu][j])>norm(UV[nu][e][j]) && norm(UV[nu][mu][j])>norm(UV[nu][tau][j]) )
                   { AV[nu][i][j][e]=real( UV[nu][mu][j]*sqrt(re2) / C<e,mu>(HfV[nu][i],kV[i][j]) );
                     AV[nu][i][j][mu]=real( UV[nu][mu][j]*sqrt(rmu2) / C<mu,mu>(HfV[nu][i],kV[i][j]) );
                     AV[nu][i][j][tau]=real( UV[nu][mu][j]*sqrt(rtau2) / C<tau,mu>(HfV[nu][i],kV[i][j]) );
                    }
                 if(norm(UV[nu][tau][j])>norm(UV[nu][e][j]) && norm(UV[nu][tau][j])>norm(UV[nu][mu][j]) )
                   { AV[nu][i][j][e]=real( UV[nu][tau][j]*sqrt(re2) / C<e,tau>(HfV[nu][i],kV[i][j]) );
                     AV[nu][i][j][mu]=real( UV[nu][tau][j]*sqrt(rmu2) / C<mu,tau>(HfV[nu][i],kV[i][j]) );
                     AV[nu][i][j][tau]=real( UV[nu][tau][j]*sqrt(rtau2) / C<tau,tau>(HfV[nu][i],kV[i][j]) );
                    }

                 re2=Delta*real(C<e,e>(HfV[antinu][i],kV[i][j]));
                 rmu2=Delta*real(C<mu,mu>(HfV[antinu][i],kV[i][j]));
                 rtau2=Delta*real(C<tau,tau>(HfV[antinu][i],kV[i][j]));

                 if(norm(UV[antinu][e][j])>norm(UV[antinu][mu][j]) && norm(UV[antinu][e][j])>norm(UV[antinu][tau][j]) )
                   { AV[antinu][i][j][e]=real( UV[antinu][e][j]*sqrt(re2) / C<e,e>(HfV[antinu][i],kV[i][j]) );
                     AV[antinu][i][j][mu]=real( UV[antinu][e][j]*sqrt(rmu2) / C<mu,e>(HfV[antinu][i],kV[i][j]) );
                     AV[antinu][i][j][tau]=real( UV[antinu][e][j]*sqrt(rtau2) / C<tau,e>(HfV[antinu][i],kV[i][j]) );
                    }
                 if(norm(UV[antinu][mu][j])>norm(UV[antinu][e][j]) && norm(UV[antinu][mu][j])>norm(UV[antinu][tau][j]) )
                   { AV[antinu][i][j][e]=real( UV[antinu][mu][j]*sqrt(re2) / C<e,mu>(HfV[antinu][i],kV[i][j]) );
                     AV[antinu][i][j][mu]=real( UV[antinu][mu][j]*sqrt(rmu2) / C<mu,mu>(HfV[antinu][i],kV[i][j]) );
                     AV[antinu][i][j][tau]=real( UV[antinu][mu][j]*sqrt(rtau2) / C<tau,mu>(HfV[antinu][i],kV[i][j]) );
                    }
                 if(norm(UV[antinu][tau][j])>norm(UV[antinu][e][j]) && norm(UV[antinu][tau][j])>norm(UV[antinu][mu][j]) )
                   { AV[antinu][i][j][e]=real( UV[antinu][tau][j]*sqrt(re2) / C<e,tau>(HfV[antinu][i],kV[i][j]) );
                     AV[antinu][i][j][mu]=real( UV[antinu][tau][j]*sqrt(rmu2) / C<mu,tau>(HfV[antinu][i],kV[i][j]) );
                     AV[antinu][i][j][tau]=real( UV[antinu][tau][j]*sqrt(rtau2) / C<tau,tau>(HfV[antinu][i],kV[i][j]) );
                    }
                }
           } 
      }

// ********************************************************************************

array<array<double,NF>,NF> MixingMatrixFactors(array<MATRIX<complex<double>,NF,NF>,NF> &C,array<MATRIX<complex<double>,NF,NF>,NF> &C0,array<array<double,NF>,NF> A0)
       { array<array<double,NF>,NF> A(A0);

         for(int j=0;j<=NF-1;j++)
            { if(real(C[j][e][mu]*C0[j][e][mu])<0. && real(C[j][e][tau]*C0[j][e][tau])<0.){ A[j][e]*=-1.;}
              if(real(C[j][mu][e]*C0[j][mu][e])<0. && real(C[j][mu][tau]*C0[j][mu][tau])<0.){ A[j][mu]*=-1.;}  
              if(real(C[j][tau][e]*C0[j][tau][e])<0. && real(C[j][tau][mu]*C0[j][tau][mu])<0.){ A[j][tau]*=-1.;}  
             }

         return A;
        }

// ********************************************************************************

array<MATRIX<complex<double>,NF,NF>,NF> CofactorMatrices(MATRIX<complex<double>,NF,NF> H,array<double,NF> k)
      { array<MATRIX<complex<double>,NF,NF>,NF> CC;
         
        for(int j=0;j<=NF-1;j++)
           { CC[j][e][e] = (H[mu][mu]-k[j])*(H[tau][tau]-k[j]) - norm(H[mu][tau]);
             CC[j][e][mu] = H[mu][tau]*H[tau][e]-H[mu][e]*(H[tau][tau]-k[j]);
             CC[j][e][tau] = H[tau][mu]*H[mu][e]-H[tau][e]*(H[mu][mu]-k[j]);

             CC[j][mu][e] = H[e][tau]*H[tau][mu]-H[e][mu]*(H[tau][tau]-k[j]);
             CC[j][mu][mu] = (H[e][e]-k[j])*(H[tau][tau]-k[j]) - norm(H[e][tau]);
             CC[j][mu][tau] = H[tau][e]*H[e][mu]-H[tau][mu]*(H[e][e]-k[j]);

             CC[j][tau][e] = H[e][mu]*H[mu][tau]-H[e][tau]*(H[mu][mu]-k[j]);
             CC[j][tau][mu] = H[mu][e]*H[e][tau]-H[mu][tau]*(H[e][e]-k[j]);
             CC[j][tau][tau] = (H[e][e]-k[j])*(H[mu][mu]-k[j]) - norm(H[e][mu]);
            }

        return CC;
       }

void CofactorMatrices(MATRIX<complex<double>,NF,NF> H,array<double,NF> k,array<MATRIX<complex<double>,NF,NF>,NF> &CC)
      { for(int j=0;j<=NF-1;j++)
           { CC[j][e][e] = (H[mu][mu]-k[j])*(H[tau][tau]-k[j]) - norm(H[mu][tau]);
             CC[j][e][mu] = H[mu][tau]*H[tau][e]-H[mu][e]*(H[tau][tau]-k[j]);
             CC[j][e][tau] = H[tau][mu]*H[mu][e]-H[tau][e]*(H[mu][mu]-k[j]);

             CC[j][mu][e] = H[e][tau]*H[tau][mu]-H[e][mu]*(H[tau][tau]-k[j]);
             CC[j][mu][mu] = (H[e][e]-k[j])*(H[tau][tau]-k[j]) - norm(H[e][tau]);
             CC[j][mu][tau] = H[tau][e]*H[e][mu]-H[tau][mu]*(H[e][e]-k[j]);

             CC[j][tau][e] = H[e][mu]*H[mu][tau]-H[e][tau]*(H[mu][mu]-k[j]);
             CC[j][tau][mu] = H[mu][e]*H[e][tau]-H[mu][tau]*(H[e][e]-k[j]);
             CC[j][tau][tau] = (H[e][e]-k[j])*(H[mu][mu]-k[j]) - norm(H[e][mu]);
            }
       }


