
#include "input_class.h"

// ****************************************************

using std::string;
using std::stringstream;

using std::ifstream;
using std::cout;

using std::min;
using std::max;

using std::vector;

using namespace prefixes;
using interpolation::DISCONTINUOUS;

// **********************************************************
// **********************************************************
// **********************************************************

void Profile_loader(InputDataSqa3Earth ID,string &outputfilenamestem)
     { string rhofilename, Yefilename, vfilename, Mfilename;
       
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

       RE=rho.XMax();

       // adjust rmin and rmax given the input profiles
       altitude = ID.altitude;
       azimuth = ID.azimuth;

       // *********************

       stringstream filename;
       string comma(","), colon(":"), dash("-");

       filename.str("");
       filename<<outputfilenamestem;
       filename<<colon<<altitude<<comma<<azimuth;  

       outputfilenamestem=filename.str();

       altitude *= M_PI/180.;
       azimuth *= M_PI/180.;

       // *********************

       lambdamin = 0.;
       if(altitude>=0.){ lambdamax=0.;}
       else{ lambdamax=2.*RE*sin(-altitude);}
      }

// **********************************************************
// **********************************************************
// **********************************************************

void Neutrino_loader(InputDataSqa3Earth ID,std::string &outputfilenamestem)
     { NE = ID.NE;
       EminMeV = ID.Emin;
       EmaxMeV = ID.Emax; // in MeV

       m1 = 0.;
       dm21 = ID.deltam_21;
       dm32 = ID.deltam_32; // in eV^2

       theta12V = ID.theta12; 
       theta13V = ID.theta13;
       theta23V = ID.theta23; 

       etaV[0] = 0.;
       etaV[1] = 0.;

       deltaV = ID.deltaCP;

       // *********************

       cout << "\n\nm1\t" << m1 << "\tdm21^2\t" << dm21 << "\tdm32\t" << dm32;
       cout << "\ntheta12V\t" << theta12V << "\ttheta13V\t" << theta13V << "\ttheta23V\t" << theta23V<< "\tdeltaCP\t" << deltaV;

       cout<<"\n\nNE\t"<<NE<<"\tEmin\t"<<EminMeV<<"\tEmax\t"<<EmaxMeV;
       cout.flush();

       // *********************

       stringstream filename;
       string comma(","), colon(":");

       filename<<outputfilenamestem;
       filename<<colon<<dm21<<comma<<dm32<<comma<<theta12V<<comma<<theta13V<<comma<<theta23V<<comma<<deltaV;
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

       deltaV*=M_PI/180.; cdeltaV=cos(deltaV); sdeltaV=sin(deltaV);
      }

       
