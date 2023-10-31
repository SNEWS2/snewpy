
#include "EMEWS.h"

#ifndef update_H
#define update_H

std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > UpdateSm(double lambdaminus,double lambdaplus,std::vector<std::vector<std::array<double,NY> > > &Y,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0,std::vector<std::vector<MATRIX<std::complex<double>,NF,NF> > > &Smprior);

std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > UpdateC(double lambda);

std::vector<std::vector<std::array<std::array<double,NF>,NF> > > UpdateA(std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C,std::vector<std::vector<std::array<MATRIX<std::complex<double>,NF,NF>,NF> > > &C0,std::vector<std::vector<std::array<std::array<double,NF>,NF> > > A0);

#endif
