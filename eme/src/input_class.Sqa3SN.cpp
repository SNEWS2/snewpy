
#include "input_class.h"

// ****************************************************

using std::string;
using std::stringstream;

using std::ifstream;
using std::cout;

using std::min;
using std::max;

using std::vector;

using std::unary_function;

using std::log;

using namespace prefixes;
using interpolation::DISCONTINUOUS;

// **********************************************************
// **********************************************************
// **********************************************************

void Profile_loader(InputDataSqa3SN ID,string &outputfilenamestem)
     { string rhofilename, Yefilename, vfilename, Mfilename;
       
       rmin = ID.rmin;
       rmax = ID.rmax;
 
       rhofilename = ID.densityprofile;
       Yefilename = ID.electronfraction;

       //******************************************************

       cout<<"\n\n*********************************************************\n";
       cout<<"\nrho\t"<<rhofilename<<"\nYe\t"<<Yefilename;
       cout.flush();

       // *********************
       // load rho and Ye data

       rho.Open(rhofilename,'#'); 
       Ye.Open(Yefilename,'#');

       lnrho=rho;
       lnrho.TransformX(static_cast<double(*)(double)>(log));
       lnrho.TransformY(static_cast<double(*)(double)>(log));

       // adjust rmin and rmax given the input profiles
       rmin=max(rmin,rho.XMin());
       rmax=min(rmax,rho.XMax());

       //******************************************************
       cout<<"\nrmin,rmax\t"<<rmin<<"\t"<<rmax; 
       cout.flush(); 

       // *********************

       stringstream filename;
       string comma(","), colon(":"), dash("-");

       filename.str("");
       filename<<outputfilenamestem;
       filename<<colon<<rmin<<comma<<rmax;  
      }

// **********************************************************
// **********************************************************
// **********************************************************

void Neutrino_loader(InputDataSqa3SN ID,std::string &outputfilenamestem)
     { NE = ID.NE;
       EminMeV = ID.Emin;
       EmaxMeV = ID.Emax; // in MeV

       m1 = 0.;
       dm21 = ID.deltam_21;
       dm32 = ID.deltam_32; // in eV^2

       theta12V = ID.theta12; 
       theta13V = ID.theta13;
       theta23V = ID.theta23; 

       alphaV[0] = 0.;
       alphaV[1] = 0.;
       alphaV[2] = 0.;

       betaV[0] = 0.;
       betaV[1] = 0.;
       betaV[2] = 0.;

       epsilonV = ID.deltaCP;

       // *********************

       cout << "\n\nm1\t" << m1 << "\tdm21^2\t" << dm21 << "\tdm32\t" << dm32;
       cout << "\ntheta12V\t" << theta12V << "\ttheta13V\t" << theta13V << "\ttheta23V\t" << theta23V;
       cout << "\nalpha1V\t" << alphaV[0] << "\talpha2V\t" << alphaV[1] << "\talpha3V\t" << alphaV[2];
       cout << "\nbeta1V\t" << betaV[0] << "\tbeta2V\t" << betaV[1] << "\tbeta3V\t" << betaV[2];

       cout<<"\n\nNE\t"<<NE<<"\tEmin\t"<<EminMeV<<"\tEmax\t"<<EmaxMeV;
       cout.flush();

       // *********************

       stringstream filename;
       string comma(","), colon(":");

       filename<<outputfilenamestem;
       filename<<colon<<dm21<<comma<<dm32<<comma<<theta12V<<comma<<theta13V<<comma<<theta23V<<comma<<epsilonV;
       filename<<colon<<NE<<comma<<EminMeV<<comma<<EmaxMeV;

       outputfilenamestem=filename.str();

       // ******************************************************

       // unit conversion to cgs
       Emin=EminMeV*mega*cgs::units::eV;
       Emax=EmaxMeV*mega*cgs::units::eV;

       m1*=1.*cgs::units::eV/cgs::constants::c2;
       dm21*=1.*cgs::units::eV*cgs::units::eV/cgs::constants::c4;
       dm32*=1.*cgs::units::eV*cgs::units::eV/cgs::constants::c4;

       theta12V*=M_PI/180.; c12V=cos(theta12V); s12V=sin(theta12V);
       theta13V*=M_PI/180.; c13V=cos(theta13V); s13V=sin(theta13V);
       theta23V*=M_PI/180.; c23V=cos(theta23V); s23V=sin(theta23V);

       epsilonV*=M_PI/180.; cepsilonV=cos(epsilonV); sepsilonV=sin(epsilonV);
      }

       
