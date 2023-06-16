
#include "potentials.h"

// ***************************************

double Ve(double rho, double Ye){ return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;}

double dVedr(double rho, double drhodr, double Ye, double dYedr){ return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp) * (drhodr*Ye + rho*dYedr );}

double Vmu(double rho, double Ye){ return 0.;}

double dVmudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}

double Vtau(double rho, double Ye){ return 0.;}

double dVtaudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}



