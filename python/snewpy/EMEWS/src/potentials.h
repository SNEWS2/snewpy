
#include "EMEWS.h"

#ifndef potentials_H
#define potentials_H

double Ve(double rho, double Ye);
double dVedr(double rho, double drhodr, double Ye, double dYedr);

double Vmu(double rho, double Ye);
double dVmudr(double rho, double drhodr, double Ye, double dYedr);

double Vtau(double rho, double Ye);
double dVtaudr(double rho, double drhodr, double Ye, double dYedr);

#endif
