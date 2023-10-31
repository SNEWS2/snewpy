
#include "EMEWS.h"

#ifndef parameters_H
#define parameters_H

constexpr int NM=2;
enum state { nu, antinu};
inline state operator++(state &n,int){ state tmp=n; n=(state)( (int)n+1 ); return tmp;};

constexpr int NF=3;
enum flavour { e, mu, tau };
inline flavour operator++(flavour &n,int){ flavour tmp=n; n=(flavour)( (int)n+1 ); return tmp;};

// number of parametrs needed to describe neutrino S matrix
constexpr int NY=12; 

// number of energy bins 
extern int NE;

// min and max energy, the first pair are in erg, the second in MeV
extern double Emin,Emax, EminMeV,EmaxMeV;

// array of neutrino energies
extern std::vector<double> E; 

// *******************************************************

// mass of mass state1, delta m^2 differences
extern double m1,dm21,dm32;

extern double theta12V, theta13V, theta23V, deltaV;
extern std::array<double,NF-1> etaV;
extern double c12V,s12V, c13V,s13V, c23V,s23V, cdeltaV,sdeltaV;

// vacuum eigenvalues
extern std::vector<std::array<double,NF> > kV;
extern std::array<int,NF> ordering;

// angles used in computing the eigenvalues plus a bunch of them in various combinations
extern double omega1,omega2,omega3;
extern double comega1,comega2,comega3, somega1,somega2,somega3;
extern double comega1p60,comega2p60,comega3p60, somega1p60,somega2p60,somega3p60;
extern double somega12,somega13,somega23; 
extern double comega12_2,comega13_2,comega23_2, somega12_2,somega13_2,somega23_2;
extern double comega12p120_2,comega13p120_2,comega23p120_2, somega12p120_2,somega13p120_2,somega23p120_2;

// vacuum Hamiltonian and mixing matrices
extern std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > HfV;
extern std::array<MATRIX<std::complex<double>,NF,NF>,NM> UV;

// vacuum values of the off-diagonal elements of the cofactor matrices
extern std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > CV;

// mixing matrix element prefactors
extern std::vector<std::vector<std::array<std::array<double,NF>,NF> > > AV;

// *******************************************************

// miniumum and maximum radius for calculation
extern double lambdamin, lambdamax;
extern double rmin, rmax;
extern double RE;
extern double altitude, azimuth;

// *******************************************************

// initial mixing matrices, needs to be global
extern std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > U0; 

#endif






