
#include "EMEWS.h"

#ifndef output_H
#define output_H

extern std::vector<std::string> fPvsrfilename;

void Initialize_Output(std::string outputfilenamestem,std::ofstream &fPvslambda,std::ofstream &fHvslambda);

void Close_Output(std::ofstream &fHvslambda);

void Output_Pvslambda(std::ofstream &fPvslambda,double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_PvsE(std::ofstream &fPvsE,std::string outputfilenamestem,double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_Hvslambda(std::ofstream &fHvslambda,double lambda,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_Pvsr(std::ofstream &fPvsr,double r,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_PvsE(std::ofstream &fPvsE,std::string outputfilenamestem,double r,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_Hvsr(std::ofstream &fHvsr,double r,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

void Output_PvsEat10kpc(std::ofstream &fFvsE,std::string outputfilenamestem,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Scumulative);

#endif
