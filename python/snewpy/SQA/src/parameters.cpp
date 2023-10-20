
#include "parameters.h"

using std::complex;
using std::array;
using std::vector;

// *******************************************************

// number of energy bins 
int NE;

// min and max energy, the first pair are in erg, the second in MeV
double Emin,Emax, EminMeV,EmaxMeV;

// array of neutrino energies
//array<double,NE> E; 
vector<double> E; 

// *******************************************************

// mass of mass state1, delta m^2 differences
double m1,dm21,dm32;

double theta12V, theta13V, theta23V, deltaV;
array<double,NF-1> etaV;
double c12V,s12V, c13V,s13V, c23V,s23V, cdeltaV,sdeltaV;

// vacuum eigenvalues
vector<array<double,NF> > kV;

array<int,NF> ordering;

// angles used in computing the eigenvalues plus a bunch of them in various combinations
double omega1,omega2,omega3;
double comega1,comega2,comega3, somega1,somega2,somega3;
double comega1p60,comega2p60,comega3p60, somega1p60,somega2p60,somega3p60;
double somega12,somega13,somega23; 
double comega12_2,comega13_2,comega23_2, somega12_2,somega13_2,somega23_2;
double comega12p120_2,comega13p120_2,comega23p120_2, somega12p120_2,somega13p120_2,somega23p120_2;

// vacuum Hamiltonian and mixing matrices
vector<vector<MATRIX<complex<double>,NF,NF> > > HfV;
array<MATRIX<complex<double>,NF,NF>,NM> UV;

// vacuum values of the off-diagonal elements of the cofactor matrices
vector<vector<array<MATRIX<complex<double>,NF,NF>,NF> > > CV;

// mixing matrix element prefactors
vector<vector<array<array<double,NF>,NF> > > AV;

// *******************************************************

// miniumum and maximum radius for calculation
double lambdamin, lambdamax;
double rmin, rmax;
double RE;
double altitude, azimuth;

// *******************************************************

// initial mixing matrices, needs to be global
vector<vector<MATRIX<complex<double>,NF,NF> > > U0; 
















