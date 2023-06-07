
#include <cmath>

#include "mstl.h"

// ******************************************************************************
// ******************************************************************************
// ******************************************************************************

#if !defined(_UNITS_AND_CONSTANTS)
#define _UNITS_AND_CONSTANTS

namespace units{ 

struct base_unit;
struct derived_unit;
template<int x> struct per;

struct base_unit { 
                   //static const double value=1.;                   
                   static constexpr double value=1.;

                   base_unit(void){;}
                  };
                  
struct derived_unit { const double value;     
                      explicit derived_unit(const double &VALUE) : value(VALUE) {;}                      
                      virtual ~derived_unit(void){;}
                     };
                     
template<int x> struct per : public derived_unit
                { per(void) : derived_unit(1./x) {;}
                  per(const per &P) : derived_unit(P) {;}
                 };
                  
derived_unit operator*(const base_unit&,const base_unit&);
derived_unit operator*(const base_unit&,const derived_unit&);
derived_unit operator*(const derived_unit&,const base_unit&);
derived_unit operator*(const derived_unit&,const derived_unit&);

derived_unit operator/(const base_unit&,const base_unit&);
derived_unit operator/(const base_unit&,const derived_unit&);
derived_unit operator/(const derived_unit&,const base_unit&);
derived_unit operator/(const derived_unit&,const derived_unit&);

double operator*(const base_unit&,const double&);
double operator*(const double&,const base_unit&);
double operator*(const derived_unit&,const double&);
double operator*(const double&,const derived_unit&);

double operator/(const base_unit&,const double&);
double operator/(const double&,const base_unit&);
double operator/(const derived_unit&,const double&);
double operator/(const double&,const derived_unit&);

} // end of namespace units 

// ******************************************************************************
// ******************************************************************************
// ******************************************************************************

// prefixes
namespace prefixes{ 

extern const double yotta, zetta, exa, peta, tera;
extern const double giga, mega, kilo, hecto, deca;
extern const double deci, centi, milli, micro, nano;
extern const double pico, femto, atto, zepto, yocto;

} // end of namespace prefixes

// ******************************************************************************
// ******************************************************************************
// ******************************************************************************

namespace units{ 

extern const derived_unit rad, sr;
extern const derived_unit radian, steradian;
                     
extern const derived_unit degree;

extern const per<100> percent;
extern const per<1000> permille;

} // end of namespace units

// ******************************************************************************
// ******************************************************************************
// ******************************************************************************

namespace SI{ 
namespace units{ 

extern const ::units::base_unit m,kg,s,A,K,mol,cd;
extern const ::units::base_unit metre, kilogram, second, ampere, kelvin, mole, candela; 

extern const ::units::derived_unit Hz, J, N, Pa, W;
extern const ::units::derived_unit hertz, joule, newton, pascal, watt;
              
extern const ::units::derived_unit C, V,ohm, S, F, Wb, H, T;
extern const ::units::derived_unit coulomb, volt, siemen, farad, weber, henry, tesla;
              
extern const ::units::derived_unit lm, lx; 

extern const ::units::derived_unit Bq, Gy, Sv;

extern const ::units::derived_unit kat;

// non-SI units that are accepted for use with SI
extern const ::units::derived_unit min, h, d, y;
extern const ::units::derived_unit minute, hour, day, year;

extern const ::units::derived_unit m2, m3, L, ha, Angstrom, b, t;

// non-SI units that are accepted for use with SI but must be experimentally measured
extern const ::units::derived_unit eV, u;
extern const ::units::derived_unit au, pc, ly;

// units of other systems in terms of SI
extern const ::units::derived_unit cm, g, erg, dyn, esu, statV;

} // end of namespace units
} // end of namespace SI

// ************************************************

namespace cgs{ 
namespace units{ 

extern const ::units::base_unit cm,g,s,K,mol,esu;
extern const ::units::base_unit centimetre,gram,second,kelvin,mole;

extern const ::units::derived_unit dyn, erg;

extern const ::units::derived_unit statV;

extern const ::units::derived_unit min, h, d, y;
extern const ::units::derived_unit minute, hour, day, year;

extern const ::units::derived_unit cm2, cm3;

// non-cgs units that must be experimentally measured
extern const ::units::derived_unit eV, u;
extern const ::units::derived_unit au, pc, ly;

// units of other systems in terms of cgs
extern const ::units::derived_unit m, kg, J, N, C, V;
extern const ::units::derived_unit b;

} // end of namespace units
} // end of namespace cgs

// ************************************************

namespace natural{ 
namespace units{ 

extern const ::units::base_unit eV;

extern const ::units::derived_unit u;

// units of other systems in terms of natural
extern const ::units::derived_unit m, kg, s, C, K;
extern const ::units::derived_unit J, N, W, V;

extern const ::units::derived_unit min, h, d, y;
extern const ::units::derived_unit minute, hour, day, year;

} // end of namespace units
} // end of namespace natural

// ************************************************

namespace Planck{ 
namespace units{ 

extern const ::units::base_unit tPl, lPl, mPl, qPl, TPl;

extern const ::units::derived_unit omegaPl, EPl, FPl, pPl, PPl; 
extern const ::units::derived_unit rhoPl; 
extern const ::units::derived_unit IPl, VPl, ZPl;

// units of other systems in terms of Planck units
extern const ::units::derived_unit m, kg, s, C, K;
extern const ::units::derived_unit J, N, W, V;

extern const ::units::derived_unit min, h, d, y;
extern const ::units::derived_unit minute, hour, day, year;

extern const ::units::derived_unit eV, u;

} // end of namespace units
} // end of namespace Planck
                  
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************
// ******************************************************************************

namespace constants{ // scientific quantities

extern const double NA;
extern const double alphaEM, alphaS;
extern const double sin2ThetaW, gA;
extern const double Vud;

extern const double h0; 
extern const double OmegaG0, OmegaNu0, OmegaR0, OmegaB0, OmegaM0, OmegaK0, OmegaL0, Omega0;
extern const double sigma8, ns;

} // end of namespace constants

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

namespace SI{ 
namespace constants{ 

extern const double c, c2, c4, hbar, hbarc, kB;
extern const double e, GN, GF;
extern const double mu0, epsilon0;

extern const double Me, Mmu, Mtau;
extern const double Mu, Md,  Ms;
extern const double Mc, Mb,  Mt;
extern const double MW, MZ;

extern const double Mpi0, MK0;

extern const double Mp, DeltaMnp, Mn;
extern const double BD, MD;
              
extern const double Ry, muB,muN, re, sigmaT;

extern const double Mearth, Msun;
extern const double Rearth, Rsun;
extern const double RSearth, RSsun;
extern const double Lsun;

extern const double H0, rhoC0;
extern const double TCMB0, Tnu0;

} // end of namespace constants
} // end namespace SI

// ***************************************************

namespace cgs{ 
namespace constants{ 

extern const double c, c2, c4, hbar, hbarc,kB;
extern const double e, GN, GF;

extern const double Me, Mmu, Mtau;
extern const double Mu, Md,  Ms;
extern const double Mc, Mb,  Mt;
extern const double MW, MZ;

extern const double Mpi0, MK0;

extern const double Mp, DeltaMnp, Mn;
extern const double BD, MD;
              
extern const double Ry, muB,muN, re, sigmaT;

extern const double Mearth, Msun;
extern const double Rearth, Rsun;
extern const double RSearth, RSsun;
extern const double Lsun;

extern const double H0, rhoC0;
extern const double TCMB0, Tnu0;

} // end of namespace constants
} // end namespace cgs

// ***************************************************
 
namespace natural{ 
namespace constants{ 

extern const double c, c2, c4, hbar, hbarc, kB;
extern const double e, GN, GF;

extern const double Me, Mmu, Mtau;
extern const double Mu, Md,  Ms;
extern const double Mc, Mb,  Mt;
extern const double MW, MZ;

extern const double Mpi0, MK0;

extern const double Mp, DeltaMnp, Mn;
extern const double BD, MD;
              
extern const double Ry, muB,muN, re, sigmaT; 

extern const double au, pc, ly;
extern const double Mearth, Msun;
extern const double Rearth, Rsun;
extern const double RSearth, RSsun;
extern const double Lsun;

extern const double H0, rhoC0;
extern const double TCMB0, Tnu0; 

} // end of namespace constants
} // end of namespace natural

// ***************************************************
 
namespace Planck{ 
namespace constants{ 

extern const double c, c2, c4, hbar, hbarc, kB;
extern const double e, GN, GF;

extern const double Me, Mmu, Mtau;
extern const double Mu, Md,  Ms;
extern const double Mc, Mb,  Mt;
extern const double MW, MZ;

extern const double Mpi0, MK0;

extern const double Mp, DeltaMnp, Mn;
extern const double BD, MD;
              
extern const double Ry, muB,muN, re, sigmaT;  

extern const double au, pc, ly;
extern const double Mearth, Msun;
extern const double Rearth, Rsun;
extern const double RSearth, RSsun;
extern const double Lsun;

extern const double H0, rhoC0;
extern const double TCMB0, Tnu0;

} // end of namespace constants
} // end of namespace Planck

#endif

