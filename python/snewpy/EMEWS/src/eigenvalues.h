
#include "EMEWS.h"

#ifndef eigenvalues_H
#define eigenvalues_H

double Q(MATRIX<std::complex<double>,NF,NF> Hf);
double R(MATRIX<std::complex<double>,NF,NF> Hf);

double omega(double Q,double sqrtminusQ,double R);
double piminusomegabar(double Qbar,double sqrtminusQbar,double Rbar);

// ********************************************************************

double k1(double T,double sqrtminusQ,double comega_3,double somega_3);
double k2(double T,double sqrtminusQ,double comega_3,double somega_3);
double k3(double T,double sqrtminusQ,double comega_3,double somega_3);

double k1bar(double Tbar,double sqrtminusQbar,double cpiminusomegabar_3,double spiminusomegabar_3);
double k2bar(double Tbar,double sqrtminusQbar,double cpiminusomegabar_3,double spiminusomegabar_3);
double k3bar(double Tbar,double sqrtminusQbar,double cpiminusomegabar_3,double spiminusomegabar_3);

std::array<double,NF> k(MATRIX<std::complex<double>,NF,NF> Hf);
std::array<double,NF> kbar(MATRIX<std::complex<double>,NF,NF> Hfbar);

// *********************************************************************

std::array<double,NF> deltak(MATRIX<std::complex<double>,NF,NF> Hf);
std::array<double,NF> deltakbar(MATRIX<std::complex<double>,NF,NF> Hfbar);

std::array<double,NF> deltak(std::array<double,NF> k);
std::array<double,NF> deltakbar(std::array<double,NF> kbar);

// *********************************************************************
// *********************************************************************
// *********************************************************************

double Evaluate_C(MATRIX<std::complex<double>,NF,NF> Hf,double comegai,double somegai);
double Evaluate_Cbar(MATRIX<std::complex<double>,NF,NF> Hf,double comegaip60,double somegaip60);
double Evaluate_a(MATRIX<std::complex<double>,NF,NF> Hf,double comegai,double somegai);

void Evaluate_omega1(void);
void Evaluate_omega2(void);
void Evaluate_omega3(void);

void Evaluate_comega1(void); void Evaluate_somega1(void);
void Evaluate_comega2(void); void Evaluate_somega2(void);
void Evaluate_comega3(void); void Evaluate_somega3(void);

void Evaluate_comega1p60(void); void Evaluate_somega1p60(void);
void Evaluate_comega2p60(void); void Evaluate_somega2p60(void);
void Evaluate_comega3p60(void); void Evaluate_somega3p60(void);

void Evaluate_somega12(void);
void Evaluate_somega13(void);
void Evaluate_somega23(void);

void Evaluate_comega12_2(void);
void Evaluate_comega13_2(void);
void Evaluate_comega23_2(void);
void Evaluate_somega12_2(void);
void Evaluate_somega13_2(void);
void Evaluate_somega23_2(void);

void Evaluate_comega12p120_2(void);  
void Evaluate_comega13p120_2(void);  
void Evaluate_comega23p120_2(void);
void Evaluate_somega12p120_2(void);  
void Evaluate_somega13p120_2(void);  
void Evaluate_somega23p120_2(void);

#endif
