
#include "Sqa3.h"

// *************************************************************************

//#include <complex>
using std::complex;
using std::polar;
using std::abs;
using std::arg;
using std::real;
using std::imag;
using std::norm;

//#include <cstdarg>
using std::va_list;

//#include<iostream>
using::std::cout;

//#include<ostream>
using std::ostream;
using std::endl;
using std::flush;

//#include<fstream>
using std::ifstream;
using std::ofstream;

//#include<sstream>
using std::stringstream;

//#include<algorithm>
using std::min;
using std::max;
using std::sort;
using std::swap;
using std::lower_bound;
using std::upper_bound;

//#include<string>
using std::string;

//#include <utility>
using std::pair;

//#include<limits>
using std::numeric_limits;

//#include<vector>
using std::vector;

//#include<array>
using std::array;

//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
using namespace pybind11;

//#include <boost/python.hpp>
//using namespace boost::python;

// ************************

//#include "mstl.h"
using namespace prefixes;
using interpolation::DISCONTINUOUS;

// ************************ Neutrino Potentials **************************

// DISCONTINUOUS is a cubic spline interpolator based on Akima's algorithm but it can handle discontinuities 
interpolation::DISCONTINUOUS rho, Ye;

// ********************************************************************** 

PYBIND11_MODULE(Sqa3Earth, m)
{   class_<InputDataSqa3Earth>(m, "InputDataSqa3Earth")
        .def(init<>())
        .def_readwrite("altitude", &InputDataSqa3Earth::altitude)
        .def_readwrite("azimuth", &InputDataSqa3Earth::azimuth)
        .def_readwrite("outputfilenamestem", &InputDataSqa3Earth::outputfilenamestem)
        .def_readwrite("densityprofile", &InputDataSqa3Earth::densityprofile)
        .def_readwrite("electronfraction", &InputDataSqa3Earth::electronfraction)
        .def_readwrite("Emin", &InputDataSqa3Earth::Emin)
        .def_readwrite("Emax", &InputDataSqa3Earth::Emax)
        .def_readwrite("deltam_21", &InputDataSqa3Earth::deltam_21)
        .def_readwrite("deltam_32", &InputDataSqa3Earth::deltam_32)
        .def_readwrite("theta12", &InputDataSqa3Earth::theta12)
        .def_readwrite("theta13", &InputDataSqa3Earth::theta13)
        .def_readwrite("theta23", &InputDataSqa3Earth::theta23)
        .def_readwrite("deltaCP", &InputDataSqa3Earth::deltaCP)
        .def_readwrite("accuracy", &InputDataSqa3Earth::accuracy)
        .def_readwrite("stepcounterlimit", &InputDataSqa3Earth::stepcounterlimit)
        .def_readwrite("NE", &InputDataSqa3Earth::NE)
        .def_readwrite("outputflag", &InputDataSqa3Earth::outputflag)
        ;

    /*class_<vector<vector<vector<vector<double> > > > >("stl_vectorx4")
        .def(vector_indexing_suite<vector<vector<vector<vector<double> > > > >());
    class_<vector<vector<vector<double> > > >("stl_vectorx3")
        .def(vector_indexing_suite<vector<vector<vector<double> > > >());
    class_<vector<vector<double> > >("stl_vectorx2")
        .def(vector_indexing_suite<vector<vector<double> > >());
    class_<vector<double> >("stl_vector_double")
        .def(vector_indexing_suite<vector<double> >());*/


    m.def("RunSqa3Earth", &RunSqa3Earth);
}

// ********************************* MAIN ****************

vector<vector<vector<vector<double> > > > RunSqa3Earth(InputDataSqa3Earth ID)
    { vector<vector<vector<vector<double> > > > PPfm;

      try{ string outputfilenamestem; 

           // **********************************************

           outputfilenamestem  = ID.outputfilenamestem;

           // *****

           Profile_loader(ID,outputfilenamestem);  
           Neutrino_loader(ID,outputfilenamestem);

           cout<<"\n\noutput filename stem\t"<<outputfilenamestem; cout.flush();  

           // ******************************************************
           // integration domians

           // ND is the number of domains, rs the rmin,rmax and the discontinutities
           double lambda,lambda0,dlambda,dlambdamin, r;           
           int ND;
           vector<double> lambdas;

           rho.FindDomains();

           lambdas.push_back(lambdamin+0.0*cgs::units::cm); 
           for(int d=1;d<=static_cast<int>(rho.NDiscontinuities());d++)
               { r=rho.Discontinuity(d); 
                 if(r>RE*cos(altitude))
                   { lambdas.push_back( -RE*sin(altitude)-sqrt(r-RE*cos(altitude))*sqrt(r+RE*cos(altitude)) );
                     lambdas.push_back( -RE*sin(altitude)+sqrt(r-RE*cos(altitude))*sqrt(r+RE*cos(altitude)) );
                    } 
                }
           lambdas.push_back(lambdamax-0.0*cgs::units::cm); 

           sort(lambdas.begin(),lambdas.end());
           ND=lambdas.size()-1;

           cout<<"\n\nNumber of domains\t"<<ND;

           // ******************************************************

           ofstream fPvslambda;           
           ofstream fPvsE;
           ofstream fHvslambda;               // a file for anything else
           ofstream fSvsE;               

           if(ID.outputflag==true){ Initialize_Output(outputfilenamestem,fPvslambda,fHvslambda);}

           // *****************************************************
           // *****************************************************
           // *****************************************************

           int i;

           E = vector<double>(NE);
           kV = vector<array<double,NF> >(NE);
           HfV=vector<vector<MATRIX<complex<double>,NF,NF> > >(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
           CV=vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > >(NM,vector<array<MATRIX<complex<double>,NF,NF>,NF> >(NE));
           AV=vector<vector<array<array<double,NF>,NF> > >(NM,vector<array<array<double,NF>,NF> >(NE));

	   // vectors of energies at infinity and vacuum eigenvalues at infinity
           for(i=0;i<=NE-1;i++)
              { E[i] = ((NE-1.-i)*Emin + i*Emax) / (NE-1.);
                kV[i][0] = m1*m1 * cgs::constants::c4 /2./E[i];
                kV[i][1] = (m1*m1 + dm21) * cgs::constants::c4 /2./E[i];
                kV[i][2] = (m1*m1 + dm21 + dm32) * cgs::constants::c4 /2./E[i];
               }

           // determine eigenvalue ordering
           if(dm32>0.){ ordering[0]=0; ordering[1]=1; ordering[2]=2;}
           else{ ordering[0]=2; ordering[1]=0; ordering[2]=1;}

           Evaluate_omega1();     Evaluate_omega2();     Evaluate_omega3();
           Evaluate_comega1();    Evaluate_comega2();    Evaluate_comega3();
           Evaluate_somega1();    Evaluate_somega2();    Evaluate_somega3();
           Evaluate_comega1p60(); Evaluate_comega2p60(); Evaluate_comega3p60();
           Evaluate_somega1p60(); Evaluate_somega2p60(); Evaluate_somega3p60();

           Evaluate_somega12();        Evaluate_somega13();        Evaluate_somega23();
           Evaluate_comega12_2();      Evaluate_comega13_2();      Evaluate_comega23_2();
           Evaluate_somega12_2();      Evaluate_somega13_2();      Evaluate_somega23_2();
           Evaluate_comega12p120_2();  Evaluate_comega13p120_2();  Evaluate_comega23p120_2();
           Evaluate_somega12p120_2();  Evaluate_somega13p120_2();  Evaluate_somega23p120_2();

           // vaccum mixing matrices and Hamiltonians at infinity
           Evaluate_UV(); 
	   Evaluate_HfV(); 
           // cofactor matrices in vacuum
           Evaluate_CV();
           // mixing matrix element prefactors in vacuum
           Evaluate_AV();

           // *****************************************************
           // quantities evaluated at inital point

           // MSW potential matrix
           double rrho = rho(lambdas.front());
           double YYe=Ye(lambdas.front());

           MATRIX<complex<double>,NF,NF> VfMSW0, VfMSWbar0, Hf0,Hfbar0;
           array<double,NF> k0, kbar0, deltak0, deltakbar0;

           VfMSW0[e][e]=Ve(rrho,YYe); 
           VfMSW0[mu][mu]=Vmu(rrho,YYe); 
           VfMSW0[tau][tau]=Vtau(rrho,YYe); 
           VfMSWbar0=-Conjugate(VfMSW0);

           // cofactor matrices at initial point - will be recycled as cofactor matrices at beginning of every step           
           vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > C0(NM,vector<array<MATRIX<complex<double>,NF,NF>,NF> >(NE));
           // mixing matrix element prefactors at initial point - will be recycled like C0         
           vector<vector<array<array<double,NF>,NF> > > A0(NM,vector<array<array<double,NF>,NF> >(NE));
           // mixing matrices at initial point, not recycled
           U0 = vector<vector<MATRIX<complex<double>,NF,NF> > >(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)); 

           // mixing angles to MSW basis at initial point and assign A0
	   for(i=0;i<=NE-1;i++)
              { Hf0=HfV[nu][i]+VfMSW0;
                k0=k(Hf0);
                deltak0=deltak(k0);
                C0[nu][i]=CofactorMatrices(Hf0,k0);
                for(int j=0;j<=NF-1;j++)
                   { if( real(C0[nu][i][j][mu][e]*CV[nu][i][j][mu][e])<0. ){ A0[nu][i][j][e]=-AV[nu][i][j][e];} else{ A0[nu][i][j][e]=AV[nu][i][j][e];}
                     A0[nu][i][j][mu]=AV[nu][i][j][mu];
                     if( real(C0[nu][i][j][mu][tau]*CV[nu][i][j][mu][tau])<0. ){ A0[nu][i][j][tau]=-AV[nu][i][j][tau];} else{ A0[nu][i][j][tau]=AV[nu][i][j][tau];}
                    }
                U0[nu][i]=U(deltak0,C0[nu][i],A0[nu][i]);

                Hfbar0=HfV[antinu][i]-VfMSW0;
                kbar0=kbar(Hfbar0);
                deltakbar0=deltakbar(kbar0);
                C0[antinu][i]=CofactorMatrices(Hfbar0,kbar0);
                for(int j=0;j<=NF-1;j++)
                   { if( real(C0[antinu][i][j][mu][e]*CV[antinu][i][j][mu][e])<0. ){ A0[antinu][i][j][e]=-AV[antinu][i][j][e];} else{ A0[antinu][i][j][e]=AV[antinu][i][j][e];}
                     A0[antinu][i][j][mu]=AV[antinu][i][j][mu];
                     if( real(C0[antinu][i][j][mu][tau]*CV[antinu][i][j][mu][tau])<0. ){ A0[antinu][i][j][tau]=-AV[antinu][i][j][tau];} else{ A0[antinu][i][j][tau]=AV[antinu][i][j][tau];}
                    }
                U0[antinu][i]=Conjugate(U(deltakbar0,C0[antinu][i],A0[antinu][i]));
               }
 
           // ******************************************************
           // ******************************************************
           // ******************************************************
           // quantities needed for the calculation
           
           double maxerror,increase=3.,accuracy;
           bool repeat, finish, resetflag, output;
           int counterout,step;

           // *************

           accuracy = ID.accuracy;
           step = ID.stepcounterlimit;

           // *************
           // variables followed as a function of lambda
           // Y are the parameters for the S matrices, Y0 the initial values for each RK step, Yerror the RK errors on the parameters

           vector<vector<array<double,NY> > > Y(NM,vector<array<double,NY> >(NE));
           vector<vector<array<double,NY> > > Y0(NM,vector<array<double,NY> >(NE));
           vector<vector<array<double,NY> > > Yerror(NM,vector<array<double,NY> >(NE));

           // cofactor matrices
           vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > C=C0;

           // mixing matrix prefactors
           vector<vector<array<array<double,NF>,NF> > > A=A0;

           // accumulated S matrices from prior integration domains
           vector<vector<MATRIX<complex<double>,NF,NF> > > Scumulative(NM,vector<MATRIX<complex<double>,NF,NF> >(NE,UnitMatrix<complex<double> >(NF)));
           PPfm=vector<vector<vector<vector<double> > > >(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF,0.))));

           // *************

           // Runge-Kutta quantities
           int NRK,NRKOrder;
           const double *AA=NULL,**BB=NULL,*CC=NULL,*DD=NULL;
           RungeKuttaCashKarpParameters(NRK,NRKOrder,AA,BB,CC,DD);

           // RK intermediate results
           vector<vector<vector<array<double,NY> > > > Ks(NRK,vector<vector<array<double,NY> > >(NM,vector<array<double,NY> >(NE)));

           // temporaries 
           MATRIX<complex<double>,NF,NF> SS;

           // ******************************************************
           // ******************************************************
           // ******************************************************

           // start of calculation

           // loop through the domains
           for(int d=0;d<=ND-1;d++)
              { if(d==0){ lambdamin=lambdas[d];} else{ lambdamin=lambdas[d]+1.*cgs::units::cm;}
                if(d==ND-1){ lambdamax=lambdas[d+1];} else{ lambdamax=lambdas[d+1]-1.*cgs::units::cm;}

                cout<<"\nDomain\t"<<d<<":\t"<<lambdamin<<"\tto\t"<<lambdamax; cout.flush();

                // **********************************************************    

                // initialize at beginning of every domain
                lambda=lambdamin; dlambda=1e-3*cgs::units::cm; dlambdamin=4.*lambda*numeric_limits<double>::epsilon();

		#pragma omp parallel for schedule(static)
                for(i=0;i<=NE-1;i++){
                    for(state m=nu;m<=antinu;m++){ 
                        Y[m][i][0]=M_PI/2.;
	                Y[m][i][1]=M_PI/2.;
		        Y[m][i][2]=M_PI/2.;
                	Y[m][i][3]=M_PI/2.;
                        Y[m][i][4]=0.;
                        Y[m][i][5]=M_PI/2.;
                	Y[m][i][6]=M_PI/2.;
                        Y[m][i][7]=0;

                        Y[m][i][8]=1.; // The determinant of the S matrix

                        Y[m][i][9]=0.;
  	            	Y[m][i][10]=0.;
        	        Y[m][i][11]=0.;
		       }
		   }

                finish=false;
                counterout=1;

                // *********

                if(ID.outputflag==true){ output=true;}
                if(output==true){ 
                    Output_Pvslambda(fPvslambda,lambda,Y,C,A,Scumulative);
                    Output_Hvslambda(fHvslambda,lambda,Y,C,A,Scumulative);
                   } 

                // **********************************************************    

                // within each domain integrate over r
                do{ if(lambda+dlambda>lambdamax){ dlambda=lambdamax-lambda; finish=true; if(ID.outputflag==true){ output=true;} else{ output=false;};}

                    lambda0=lambda;
                    Y0=Y;
                    C0=C;
                    A0=A;

                    // beginning of RK section
                    do{ repeat=false;                         
                        // first step: assumes derivatives are evaluated at lambda
                        K(lambda,dlambda,Y,C,A,Ks[0]); 

                        // second step
                        lambda=lambda0+AA[1]*dlambda;
                        #pragma omp parallel for schedule(static)
                        for(i=0;i<=NE-1;i++){
                            for(state m=nu;m<=antinu;m++){
                                for(int j=0;j<=NY-1;j++){ Y[m][i][j] += BB[1][0] * Ks[0][m][i][j];}
		               } 
			   }
                        K(lambda,dlambda,Y,C,A,Ks[1]);

                        // remaining steps
                        for(int k=2;k<=NRK-1;k++){
                            lambda=lambda0+AA[k]*dlambda;
                            Y=Y0; 
                            #pragma omp parallel for schedule(static)
                            for(i=0;i<=NE-1;i++){
                                for(state m = nu; m <= antinu; m++){
 		                    for(int j=0;j<=NY-1;j++){
				        for(int l=0;l<=k-1;l++){ Y[m][i][j] += BB[k][l] * Ks[l][m][i][j];}
				       }
				   } 
			       } 
                            K(lambda,dlambda,Y,C,A,Ks[k]);
                           }

                        // increment all quantities and update C and A arrays
                        lambda=lambda0+dlambda;             
                        #pragma omp parallel for schedule(static)
                        for(i=0;i<=NE-1;i++){
                            for(state m=nu;m<=antinu;m++){                            
                                for(int j=0;j<=NY-1;j++){
                                    Y[m][i][j]=Y0[m][i][j]; 
     			            Yerror[m][i][j]=0.;
                        	    for(int k=0;k<=NRK-1;k++){
                                        Y[m][i][j]+=CC[k]*Ks[k][m][i][j]; 
				        Yerror[m][i][j]+=(CC[k]-DD[k])*Ks[k][m][i][j];
				       }
                                   }
                               } 
                           }

                        C=UpdateC(lambda);
                        A=UpdateA(C,C0,A0);

                        // find largest error
                        maxerror=0.; 
                        for(state m=nu;m<=antinu;m++){
                            for(i=0;i<=NE-1;i++){
                                for(int j=0;j<=NY-1;j++){ maxerror = max( maxerror, fabs(Yerror[m][i][j]) );}
			       }
                           }

                       // decide whether to accept step, if not adjust step size
                       if(maxerror>accuracy){
                          dlambda*=0.9*pow(accuracy/maxerror,1./(NRKOrder-1.));
                          if(dlambda>dlambdamin){ repeat=true;}
                         }

                        // reset integration variables to those at beginning of step
                        if(repeat==true){ lambda=lambda0; Y=Y0; C=C0; A=A0; finish=output=false;} 

                       }while(repeat==true);
                    // end of RK section

                    // check S matrices are diagonal dominated, if not then accumulate S and reset variables
                    #pragma omp parallel for schedule(static) private(SS,resetflag)
                    for(i=0;i<=NE-1;i++){
                        for(state m=nu;m<=antinu;m++){                         
                            SS=W(Y[m][i])*B(Y[m][i]); 

                            resetflag=false;

                            // test that the S matrix is close to diagonal
          	            if( norm(SS[0][0])+0.1<norm(SS[0][1]) || norm(SS[0][0])+0.1<norm(SS[0][2]) ){ resetflag=true;}
                	    if( norm(SS[0][1])+0.1<norm(SS[0][2]) ){ resetflag=true;}
	                    if( norm(SS[2][2])+0.1<norm(SS[1][2]) ){ resetflag=true;}

	                    if(resetflag!=false)
                              { // reset the S matrices
                                Scumulative[m][i]=MATRIX<complex<double>,NF,NF>( SS*Scumulative[m][i] );

                                Y[m][i][0]=Y[m][i][1]=Y[m][i][2]=Y[m][i][3]=M_PI/2.; Y[m][i][4]=0.;
                                Y[m][i][5]=Y[m][i][6]=M_PI/2.;                       Y[m][i][7]=0.;
                	        Y[m][i][8]=1.;	
                                Y[m][i][9]=Y[m][i][10]=Y[m][i][11]=0.;
                               }
                            else{ // take modulo 2 pi of phase angles
	                          Y[m][i][4]=fmod(Y[m][i][4],M_2PI);
	                          Y[m][i][7]=fmod(Y[m][i][7],M_2PI);

	                          double ipart;
	                          Y[m][i][9]=modf(Y[m][i][9],&ipart);
	                          Y[m][i][10]=modf(Y[m][i][10],&ipart); 
	                          Y[m][i][11]=modf(Y[m][i][11],&ipart);
	                         }
			   }
                       }

                    if(counterout==step){ if(ID.outputflag==true){ output=true;} else{ output=false;}; counterout=1;} else{ counterout++;}

                    if(output==true)
                      { cout<<"\nOutput at\t"<<lambda<<flush;
                        Output_Pvslambda(fPvslambda,lambda,Y,C,A,Scumulative);
                        Output_Hvslambda(fHvslambda,lambda,Y,C,A,Scumulative);
                        //Output_PvsE(fPvsE,outputfilenamestem,lambda,Y,C,A,Scumulative);
                        output=false;
                       }

                    // adjust step size based on RK error - could be moved up to RK section but better left here in case adjustments are necessary based on new S matrices
                    dlambda=min(dlambda*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dlambda);
                    dlambdamin=4.*lambda*numeric_limits<double>::epsilon();
                    dlambda=max(dlambda,dlambdamin); 

                   }while(finish==false);

                // if this is not the last domain then carry the S matrix across the domain boundaries
                if(d<=ND-2)
                  { double lambdaminus=lambdamax;
                    double lambdaplus=lambdamax+2.*cgs::units::cm;
                    Scumulative=UpdateSm(lambdaminus,lambdaplus,Y,C,A,Scumulative);

                    C0=C;
                    C=UpdateC(lambdaplus);
                    A=UpdateA(C,C0,A0);
                   } 
                else{ // output at the end of the code
                      if(ID.outputflag==true){ output=true;}
                      if(output==true){ 
                          Output_PvsE(fPvsE,outputfilenamestem,lambdamax,Y,C,A,Scumulative);
                          output=false;
                         }
                     } 

               }// end of domain loop

           Pfm(lambdas.back(),Y,C,A,Scumulative,PPfm);

           // ********************************

           if(ID.outputflag==true){ Close_Output(fHvslambda);}

          }catch(OUT_OF_RANGE<unsigned int> &OOR){ OOR.Message();}
           catch(OUT_OF_RANGE<int> &OOR){ OOR.Message();}
           catch(OUT_OF_RANGE<double> &OOR){ OOR.Message();}
           catch(EMPTY &E){ E.Message();}
           catch(FUNCTION_ERROR &FE){ FE.Message();}
           catch(BASIC_ERROR &BE){ BE.Message();}
           catch(...){ UNKNOWN_ERROR("main");}

      cout<<"\nFinished\n\a"; cout.flush();

      return PPfm;
     }


