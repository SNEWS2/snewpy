
#include "EMEWS.h"

#ifndef input_class_H
#define input_class_H

// **************************************

struct InputDataEMEWS;

void Profile_loader(InputDataEMEWS ID,std::string &outputfilenamestem);
void Neutrino_loader(InputDataEMEWS ID,std::string &outputfilenamestem);

// **********************************************************
// **********************************************************
// **********************************************************

struct InputDataEMEWS 
       { double altitude, azimuth; // in decimal degrees, altitude is negative for angles below the horizon
         std::string outputfilenamestem;
         std::string densityprofile;
         std::string electronfraction;

         int NE; // number of energies
         double Emin, Emax; // in MeV
         double deltam_21, deltam_32; // in eV^2	
         double theta12, theta13, theta23, deltaCP; // all in degrees
         double accuracy;
         double stepcounterlimit; // how often it spits out data
         bool outputflag; // whether the code outputs data as it does the integgration

         InputDataEMEWS(void) {;}  
        };

#endif
