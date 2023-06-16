
#include "units and constants.h"

// *****************************************************************************

namespace units{ 
             
// **************************

derived_unit operator*(const base_unit &b1,const base_unit &b2){ return derived_unit(b1.value*b2.value);}
derived_unit operator*(const base_unit &b,const derived_unit &d){ return derived_unit(b.value*d.value);}             
derived_unit operator*(const derived_unit &d,const base_unit &b){ return derived_unit(d.value*b.value);}             
derived_unit operator*(const derived_unit &d1,const derived_unit &d2){ return derived_unit(d1.value*d2.value);}
             
// **************************

derived_unit operator/(const base_unit &b1,const base_unit &b2){ return derived_unit(b1.value/b2.value);}             
derived_unit operator/(const base_unit &b,const derived_unit &d){ return derived_unit(b.value/d.value);}             
derived_unit operator/(const derived_unit &d,const base_unit &b){ return derived_unit(d.value/b.value);}             
derived_unit operator/(const derived_unit &d1,const derived_unit &d2){ return derived_unit(d1.value/d2.value);}

// **************************

double operator*(const base_unit &b,const double &D){ return b.value*D;}             
double operator*(const double &D,const base_unit &b){ return (D*b.value);}             
double operator*(const derived_unit &d,const double &D){ return d.value*D;}             
double operator*(const double &D,const derived_unit &d){ return (D*d.value);}
             
// **************************

double operator/(const base_unit &b,const double &D){ return b.value/D;}             
double operator/(const double &D,const base_unit &b){ return D/b.value;}             
double operator/(const derived_unit &d,const double &D){ return d.value/D;}             
double operator/(const double &D,const derived_unit &d){ return D/d.value;}
             
} // end of namespace units              

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// prefixes
namespace prefixes{ 

const double yotta=1e24, zetta=1e21, exa=1e18, peta=1e15, tera=1e12;
const double giga=1e9, mega=1e6, kilo=1e3, hecto=1e2, deca=1e1;
const double deci=1e-1, centi=1e-2, milli=1e-3, micro=1e-6, nano=1e-9;
const double pico=1e-12, femto=1e-15, atto=1e-18, zepto=1e-21, yocto=1e-24;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

namespace units{ 

const derived_unit rad(1.), sr(1.);
const derived_unit radian(rad), steradian(sr);
                     
const derived_unit degree(rad*180./M_PI);

const per<100> percent;
const per<1000> permille;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

namespace SI{ 
namespace units{ 
using namespace ::units;
using namespace prefixes;

const base_unit m, kg, s, A, K, mol, cd; 
const base_unit metre(m), kilogram(kg), second(s), ampere(A), kelvin(K), mole(mol), candela(cd);

const derived_unit Hz(1./s), J(kg*m*m/s/s), N(J/m), Pa(N/m/m), W(J/s);
const derived_unit hertz(Hz), joule(J), newton(N), pascal(Pa), watt(W);

const derived_unit C(A*s), V(W/A), ohm(V/A), S(A/V), F(C/V), Wb(m*m*kg/s/C), H(Wb/A), T(kg/s/C);
const derived_unit coulomb(C), volt(V), siemens(S), farad(F), weber(Wb), henry(H), tesla(T);
              
const derived_unit lm(cd*sr), lx(cd/m/m); 

const derived_unit Bq(1./s), Gy(J/kg), Sv(J/kg);

const derived_unit kat(mol/s);

// non-SI units that are accepted for use with SI
const derived_unit min(60.*s), h(60.*min), d(24.*h), y(31556925.2*s);             
const derived_unit minute(min), hour(h), day(d), year(y);
             
const derived_unit m2(m*m), m3(m*m*m), L(milli*m3), ha(1e4*m2), Angstrom(1e-10*m), b(1e-28*m2), t(kilo*kg);

const derived_unit au(149597870660.*m), AU(au), pc(3.0856775807e16*m), ly(0.9461e16*m);

} // end of namespace units
} // end of namespace SI

// *********************************************

namespace cgs{ 
namespace units{ 
using namespace ::units;
using namespace prefixes;

const base_unit cm, g, s, K, mol, esu;
const base_unit centimetre(cm), gram(g), second(s), kelvin(K), mole(mol);

const derived_unit dyn(g*cm/s/s), erg(dyn*cm);
const derived_unit statV(erg/esu);

const derived_unit min(60.*s), h(60.*min), d(24.*h), y(31556925.2*s);             
const derived_unit minute(min), hour(h), day(d), year(y);

const derived_unit cm2(cm*cm), cm3(cm*cm*cm);
} // end of namespace units
} // end namespace cgs

// *********************************************

namespace natural{ 
namespace units{ 
using namespace ::units;
using namespace prefixes;

const base_unit eV;

const derived_unit u(931.494043*mega*eV);
} // end of namespace units
} // end of namespace natural

// *********************************************

namespace Planck{ 
namespace units{ 
using namespace ::units;
using namespace prefixes;

const base_unit tPl, lPl, mPl, qPl, TPl;

const derived_unit omegaPl(1./tPl), EPl(1.*mPl), FPl(EPl/lPl), pPl(FPl/lPl/lPl), PPl(EPl/tPl); 

const derived_unit rhoPl(mPl/lPl/lPl/lPl); 

const derived_unit IPl(qPl/tPl), VPl(EPl/qPl), ZPl(VPl/IPl); 
} // end of namespace units
} // end of namespace Planck

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

namespace constants{ 

const double NA=6.02214076e23;
const double alphaEM=1./137.03599911, alphaS=0.1172;
const double sin2ThetaW=0.23120, gA=-1.2720;
const double Vud=0.9734;

const double h0=0.73;
const double OmegaG0=4.6e-5, OmegaNu0=3.*7./8.*OmegaG0*pow(4./11.,4./3.), OmegaR0=OmegaG0+OmegaNu0; 
const double OmegaB0=0.042, OmegaM0=0.24, OmegaK0=0., OmegaL0=0.76;
const double Omega0=/*OmegaR0+*/OmegaM0+OmegaK0+OmegaL0;
const double sigma8=0.74, ns=0.951;

} // end of namespace constants

//******************************************************************************
//******************************************************************************
//******************************************************************************

namespace SI{ 
namespace constants{
using namespace prefixes;
using namespace units;
using namespace ::constants;

const double c=299792458.*m/s, c2=c*c, c4=c2*c2, hbar=1.054571817e-34*J*s, hbarc=hbar*c, kB=1.380649e-23*J/K, e=1.602176634e-19*C;
const double GN=6.6742e-11*m3/kg/s/s;
const double mu0=4.*M_PI*1e-7*N/A/A, epsilon0=1./mu0/c2;

const double Mearth=5.974e24*kg,      Msun=1.9889e30*kg;
const double Rearth=6.37814*mega*m,   Rsun=6.961e8*m;
const double RSearth=2.*GN*Mearth/c2, RSsun=2.*GN*Msun/c2;
const double Lsun=3.846e26*J/s;

const double H0=h0*100.*kilo*m/s/(mega*pc), rhoC0=3.*H0*H0/(8.*M_PI*GN);
const double TCMB0=2.725*K, Tnu0=TCMB0*cbrt(4./11.);
} // end of namespace constants

} //end of namespace SI

//************************************
//************************************
//************************************

namespace natural{ 
namespace constants{
using namespace prefixes;
using namespace units;

const double c=1., c2=c*c, c4=c2*c2, hbar=1., hbarc=hbar*c, kB=1.;
const double e=1., GF=1.16637e-23/eV/eV;

const double Me=0.510998918*mega*eV, Mmu=105.6583568*mega*eV, Mtau=1776.99*mega*eV;
const double Mu=3.*mega*eV, Md=6.25*mega*eV, Ms=117.5*mega*eV, Mc=1.2*giga*eV, Mb=4.25*giga*eV, Mt=174.3*giga*eV;
const double MW=80.425*giga*eV, MZ=91.1876*giga*eV;

const double Mpi0=134.9766*mega*eV, MK0=497.648*mega*eV;

const double Mp=938.272029*mega*eV, DeltaMnp=1.2933318*mega*eV, Mn=Mp+DeltaMnp;
const double BD=2223.938*kilo*eV,  MD=Mp+Mn-BD;
                   
const double Ry=13.6056923*eV;
const double muB=e*hbar/2./Me/eV, muN=e*hbar/2./Mp/eV;
} // end of namespace constants

} // end of namespace natural

//************************************

namespace Planck{ 
namespace constants{
using namespace prefixes;
using namespace units;

const double c=1.*lPl/tPl, c2=c*c, c4=c2*c2, hbar=1.*EPl*tPl, hbarc=hbar*c, kB=1.*EPl/TPl;
const double GN=1.*lPl*lPl*lPl/mPl/tPl/tPl;
} // end of namespace constants
} // end of namespace Planck

//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************

// Convert defined constants and units in one unit system into other unit systems

namespace SI{ 

namespace units{
using namespace ::units;
using namespace constants;

// units of other systems in terms of SI
const derived_unit cm(centi*m), g(milli*kg), erg(g*cm*cm/s/s), dyn(g*cm/s/s);

const derived_unit eV(e*V), u(natural::units::u*eV/c2); 

const derived_unit esu(sqrt(dyn/N*mu0/4./M_PI*cm/m*c)*C), statV(erg/esu);
} //end of namespace units

namespace constants{
using namespace prefixes;
using namespace units;
using namespace ::constants;

const double GF=natural::constants::GF*pow(hbarc,3.)/eV/eV;

const double Me=natural::constants::Me*eV/c2, Mmu=natural::constants::Mmu*eV/c2, Mtau=natural::constants::Mtau*eV/c2;
const double Mu=natural::constants::Mu*eV/c2, Md=natural::constants::Md*eV/c2,   Ms=natural::constants::Ms*eV/c2;
const double Mc=natural::constants::Mc*eV/c2, Mb=natural::constants::Mb*eV/c2,   Mt=natural::constants::Mt*eV/c2;
const double MW=natural::constants::MW*eV/c2, MZ=natural::constants::MZ*eV/c2;

const double Mpi0=natural::constants::Mpi0*eV/c2, MK0=natural::constants::MK0*eV/c2;

const double Mp=natural::constants::Mp*eV/c2, DeltaMnp=natural::constants::DeltaMnp*eV/c2, Mn=natural::constants::Mn*eV/c2;
const double BD=natural::constants::BD*eV/c2, MD=natural::constants::MD*eV/c2;

const double Ry=natural::constants::Ry*eV;
const double muB=natural::constants::muB*eV/T, muN=natural::constants::muN*eV/T;

const double re=e*e/(4.*M_PI*epsilon0*Me*c2)*m, sigmaT=8.*M_PI/3.*re*re*m2;
} // end of namespace constants

} // end of namespace SI

//************************************

namespace cgs{ 

namespace units{
using namespace ::units;
using namespace constants;

// units of other systems in terms of cgs
const derived_unit m(hecto*cm), kg(kilo*g), J(kg*m*m/s/s), N(J/m);
const derived_unit b(1e-28*m*m);

const derived_unit C(sqrt(N/dyn*4.*M_PI/SI::constants::mu0/SI::constants::c)*esu), V(J/C), T(kg/s/C);

const derived_unit eV(SI::units::eV*J), u(SI::units::u*kg);
const derived_unit au(SI::units::au*m), AU(SI::units::AU*m), pc(SI::units::pc*m), ly(SI::units::ly*m);
} // end of namespace units

namespace constants{
using namespace prefixes;
using namespace units; 
using namespace ::constants;

const double c=SI::constants::c*m/s, c2=c*c, c4=c2*c2;
const double hbar=SI::constants::hbar*J*s, hbarc=hbar*c, kB=SI::constants::kB*J/K;
const double e=SI::constants::e*C, GN=SI::constants::GN*m*m*m/kg/s/s, GF=natural::constants::GF*pow(hbarc,3.)/eV/eV;

const double Me=SI::constants::Me*kg, Mmu=SI::constants::Mmu*kg, Mtau=SI::constants::Mtau*kg;
const double Mu=SI::constants::Mu*kg, Md=SI::constants::Md*kg,   Ms=SI::constants::Ms*kg;
const double Mc=SI::constants::Mc*kg, Mb=SI::constants::Mb*kg,   Mt=SI::constants::Mt*kg;
const double MW=SI::constants::MW*kg, MZ=SI::constants::MZ*kg;

const double Mpi0=SI::constants::Mpi0*kg, MK0=SI::constants::MK0*kg;

const double Mp=SI::constants::Mp*kg, DeltaMnp=SI::constants::DeltaMnp*kg, Mn=SI::constants::Mn*kg;
const double BD=SI::constants::BD*kg, MD=SI::constants::MD*kg;
               
const double Ry=SI::constants::Ry*J, re=SI::constants::re*m, sigmaT=SI::constants::sigmaT*m*m;
const double muB=SI::constants::muB*J/T, muN=SI::constants::muN*J/T;

const double Mearth=SI::constants::Mearth*kg, Rearth=SI::constants::Rearth*m;
const double Msun=SI::constants::Msun*kg, Rsun=SI::constants::Rsun*m;
const double RSearth=2.*GN*Mearth/c2, RSsun=2.*GN*Msun/c2;
const double Lsun=SI::constants::Lsun*J/s;

const double H0=SI::constants::H0/s, rhoC0=SI::constants::rhoC0*kg/m/m/m;
const double TCMB0=SI::constants::TCMB0*K, Tnu0=SI::constants::Tnu0*K;
} // end of namespace constants

} // end of namespace cgs

//************************************

namespace natural{ 

namespace units{
using namespace ::units;
using namespace constants;

// units of other systems in terms of natural
const derived_unit m(SI::units::eV/SI::constants::hbarc);
const derived_unit kg(SI::constants::c2/SI::units::eV);
const derived_unit s(SI::units::eV/SI::constants::hbar);
const derived_unit C(1./SI::constants::e); 
const derived_unit K(SI::constants::kB/SI::units::eV); 

const derived_unit J(kg*m*m/s/s), N(J/m), W(J/s), V(J/C); 

const derived_unit min(60.*s), h(60.*min), d(24.*h), y(31556925.2*s);             
const derived_unit minute(min), hour(h), day(d), year(y);
} // end of namespace units

namespace constants{
using namespace prefixes;
using namespace units;
using namespace ::constants;

const double GN=SI::constants::GN*m*m*m/kg/s/s;
                  
const double re=SI::constants::re*m, sigmaT=SI::constants::sigmaT*m*m;

const double au=1.*SI::constants::au*m, pc=1.*SI::constants::pc*m, ly=1.*SI::constants::ly*m;

const double Mearth=SI::constants::Mearth*kg, Rearth=SI::constants::Rearth*m;
const double Msun=SI::constants::Msun*kg, Rsun=SI::constants::Rsun*m;
const double RSearth=2.*GN*Mearth/c2, RSsun=2.*GN*Msun/c2;
const double Lsun=SI::constants::Lsun*J/s;

const double H0=SI::constants::H0/s, rhoC0=SI::constants::rhoC0*kg/m/m/m;
const double TCMB0=SI::constants::TCMB0*K, Tnu0=SI::constants::Tnu0*K;
} // end of namespace constants

} // end of namespace natural

//************************************

namespace Planck{ 

namespace units{
using namespace prefixes;
using namespace ::units;

// units of other systems in terms of Planck units
const derived_unit m(sqrt(SI::constants::c4/SI::constants::hbarc/SI::constants::GN) *lPl);
const derived_unit kg(sqrt(SI::constants::GN/SI::constants::hbarc) *mPl);
const derived_unit s(sqrt(SI::constants::c4*SI::constants::c/SI::constants::hbar/SI::constants::GN) *tPl);
const derived_unit C(1./sqrt(4.*M_PI*SI::constants::epsilon0*SI::constants::hbarc) *qPl); 
const derived_unit K(sqrt(SI::constants::GN*SI::constants::kB*SI::constants::kB/SI::constants::hbarc/SI::constants::c4) *TPl); 

const derived_unit J(kg*m*m/s/s), N(J/m), W(J/s), V(J/C); 

const derived_unit min(60.*s), h(60.*min), d(24.*h), y(31556925.2*s);             
const derived_unit minute(min), hour(h), day(d), year(y);                

const derived_unit eV(SI::units::eV*J), u(SI::units::u*kg);
} // end of namespace units

namespace constants{
using namespace prefixes;
using namespace units;
using namespace ::constants;

const double e=sqrt(::constants::alphaEM)*qPl, GF=natural::constants::GF/eV/eV;

const double Me=SI::constants::Me*kg, Mmu=SI::constants::Mmu*kg, Mtau=SI::constants::Mtau*kg;
const double Mu=SI::constants::Mu*kg, Md=SI::constants::Md*kg,   Ms=SI::constants::Ms*kg;
const double Mc=SI::constants::Mc*kg, Mb=SI::constants::Mb*kg,   Mt=SI::constants::Mt*kg;
const double MW=SI::constants::MW*kg, MZ=SI::constants::MZ*kg;

const double Mpi0=SI::constants::Mpi0*kg, MK0=SI::constants::MK0*kg;

const double Mp=SI::constants::Mp*kg, DeltaMnp=SI::constants::DeltaMnp*kg, Mn=SI::constants::Mn*kg;
const double BD=SI::constants::BD*kg, MD=SI::constants::MD*kg;
               
const double Ry=SI::constants::Ry*J, re=SI::constants::re*m, sigmaT=SI::constants::sigmaT*m*m;

const double au=1.*SI::constants::au*m, pc=1.*SI::constants::pc*m, ly=1.*SI::constants::ly*m;

const double Mearth=SI::constants::Mearth*kg, Rearth=SI::constants::Rearth*m;
const double Msun=SI::constants::Msun*kg, Rsun=SI::constants::Rsun*m;
const double RSearth=2.*GN*Mearth/c2, RSsun=2.*GN*Msun/c2;
const double Lsun=SI::constants::Lsun*J/s;

const double H0=SI::constants::H0/s, rhoC0=SI::constants::rhoC0*kg/m/m/m;
const double TCMB0=SI::constants::TCMB0*K, Tnu0=SI::constants::Tnu0*K;
} // end of namespace constants

} // end of namespace Planck



