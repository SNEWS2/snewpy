# -*- coding: utf-8 -*-
"""Supernova oscillation physics 

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import AltAz

from .neutrino import MassHierarchy
from .flavor  import TwoFlavor,ThreeFlavor, FourFlavor, FlavorMatrix
from .neutrino import MixingParameters, ThreeFlavorMixingParameters, FourFlavorMixingParameters
###############################################################################

class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""
    def __str__(self):
        return self.__class__.__name__
    
    @abstractmethod
    def P_ff(self, t, E)->FlavorMatrix:
        pass
    
    def apply(self, flux):
        M = self.P_ff(flux.time, flux.energy)
        M = (flux.flavor_scheme<<M<<flux.flavor_scheme)
        return M@flux
###############################################################################

class ThreeFlavorTransformation:
    "A class to ensure the mixing parameters are ThreeFlavor"
    @property
    def mix_pars(self)->ThreeFlavorMixingParameters:
        return self._mix_pars
    @mix_pars.setter
    def mix_pars(self, value:ThreeFlavorMixingParameters):
        assert isinstance(value, ThreeFlavorMixingParameters)
        self._mix_pars = value

class FourFlavorTransformation:
    "A class to ensure the mixing parameters are FourFlavor"
    @property
    def mix_pars(self)->FourFlavorMixingParameters:
        return self._mix_pars
    @mix_pars.setter
    def mix_pars(self, value:FourFlavorMixingParameters):
        if not isinstance(value, FourFlavorMixingParameters):
            raise TypeError(f"{self.__class__.__name__} requires {m_type}, but received {type(value)}")
        self._mix_pars = value
       
###############################################################################
class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def P_ff(self, t, E):
        p = FlavorMatrix.eye(ThreeFlavor)
        return p

    def apply(self, flux):
        """This transformation returns the object without transform"""
        return flux
###############################################################################
class CompleteExchange(FlavorTransformation):
    """Survival probabilities for the case when the electron flavors 
       are half exchanged with the mu flavors and the half with the tau flavors.
    """
    def P_ff(self, t, E):
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1,f2):
            return (f1.is_neutrino==f2.is_neutrino)*(f1!=f2)*0.5

        return P

###############################################################################
class ThreeFlavorDecoherence(FlavorTransformation):
    """Equal mixing of all threen eutrino matter states and antineutrino matter states"""

    def __init__(self):
        """Initialize ThreeFlavorTransformation to default case"""
        super().__init__() 

    def __str__(self):
        return f'ThreeFlavorDecoherence'

    def P_ff(self, t, E): 
        """Equal mixing so Earth matter has no effect"""
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1,f2):
            return (f1.is_neutrino==f2.is_neutrino)*1/3.
        return P

###############################################################################
# SN Transformations
###############################################################################
class SNTransformation(ABC):
    def __init__(self, mixing_params):
        self.mix_pars = mixing_params

    @abstractmethod
    def P_mf(self, t, E)->FlavorMatrix:
        pass
###############################################################################
class AdiabaticMSW(SNTransformation):
    """Adiabatic MSW effect."""

    def P_mf(self, t, E):
        return self.mix_pars.Pmf_HighDensityLimit()
###############################################################################
class NonAdiabaticMSWH(SNTransformation):
    """Nonadiabatic MSW H resonance. 
    The NonAdiabaticMSWH transformation assumes that the H resonance mixing is nonadiabatic.
    This case is relevant when a shock is present at the H resonance densities (Schirato & Fuller 2002).
    
    For the NMO the H resonance occurs in the neutrinos (Kneller & McLaughlin 2009) between ‘matter’
    states ν2 and ν3.
    In the IMO the H resonance mixes the antineutrino matter states ν̄1 and ν̄3.
    """
    def P_mf(self, t, E):
        Pmf = self.mix_pars.Pmf_HighDensityLimit()
        if self.mix_pars.mass_order == MassHierarchy.NORMAL:
            Pmf[['NU_2','NU_3'],:] = Pmf[['NU_3','NU_2'],:]
        else:
            Pmf[['NU_1_BAR','NU_3_BAR'],:] = Pmf[['NU_3_BAR','NU_1_BAR'],:]
        return Pmf
###############################################################################
class TwoFlavorDecoherence(SNTransformation):
    """equal mixing of whatever two matter states form the MSW H resonance.

    The TwoFlavorDecoherence transformation is relevant when the size of the density fluctuations
    is ≲ 10% for densities around the H resonance density —see Kneller (2010); Kneller & Mauney (2013). 
    This prescription models the effect of the turbulence as leading to 50% mixing between the
    matter states which participate in the H resonance. 
    
    In the NMO this is ν2 and ν3,
    For the IMO, the H resonance occurs in the antineutrinos between antineutrino matter states ν̄1
and ν̄3. 
    """
    
    def P_mf(self, t, E):
        Pmf = self.mix_pars.Pmf_HighDensityLimit()
        if self.mix_pars.mass_order == MassHierarchy.NORMAL:
            Pmf['NU_2']=Pmf['NU_3']=0.5*(Pmf['NU_2']+Pmf['NU_3'])
        else:
            Pmf['NU_1_BAR']=Pmf['NU_3_BAR']=0.5*(Pmf['NU_1_BAR']+Pmf['NU_3_BAR'])

        return Pmf

###############################################################################

try:
    import SNOSHEWS
except:
    SNOSHEWS = None

class SNprofile:
    """A placeholder class for the density and electron fraction profiles. Currently SNOSHEWS 
    reads the profiles from files but that might change in the future.
    """

    def __init__(self, rhofilename, Yefilename):         
        self.rhofilename = rhofilename
        self.Yefilename = Yefilename
        
class MSWEffect(SNTransformation):
    """The MSW effect using a density profile and electron 
       fraction provided by the user. Uses the SNOSHEWS module.
    """

    def __init__(self, SNprofile, rmin = 0, rmax = 1e99):
        """Initialize flavor transformation
        
        Parameters
        ----------
        SNprofile : instance of profile class
        mix_params : ThreeFlavorMixingParameters instance or None
        rmin : starting radius for calculation. rmin will be corrected by SNOSHEWS to the minimum radius of the profile if rmin is less than that value        
        rmax : the ending radius of the calculation. rmax will be corrected by SNOSHEWS to the maximum radius of the profile if rmax is greater than that value
        """
        self.SNprofile = SNprofile

    def P_mf(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        Pmf : array of 6 x 6 matrices
        """
                     
        if SNOSHEWS == None:
            print("The SNOSHEWS module cannot be found. Returning results using AdiabaticMSW prescription")
            return AdiabticMSW(self.mix_pars).get_probabilities(t,E)
        
        #input data object for SNOSHEWS
        ID = SNOSHEWS.InputDataSNOSHEWS()

        ID.outputfilenamestem = "./out/SNOSHEWS"      # stem of output filenames 

        ID.rmin = self.rmin
        ID.rmax = self.rmax

        ID.densityprofile = self.SNprofile.rhofilename     # mass density profile    
        ID.electronfraction = self.SNprofile.Yefilename    # electron fraction profile     
        
        ID.NE = len(E)         # number of energy bins
        E = E.to_value('MeV')
        ID.Emin = E[0]         # in MeV
        ID.Emax = E[-1]        # in MeV

        #MixingParameters
        ID.deltam_21 = self.mix_pars.dm21_2.value   # in eV^2
        ID.deltam_32 = self.mix_pars.dm32_2.value   # in eV^2
        ID.theta12 = self.mix_pars.theta12.value    # in degrees
        ID.theta13 = self.mix_pars.theta13.value    # in degrees
        ID.theta23 = self.mix_pars.theta23.value    # in degrees
        ID.deltaCP = self.mix_pars.deltaCP.value    # in degrees

        ID.accuracy = 1.01E-009      # controls accuracy of integrator: smaller is more accurate
        ID.stepcounterlimit = 10000  # output frequency if outputflag = True: larger is less frequent
        ID.outputflag = False        # set to True if output is desired
 
        # Do the calculation. The return is a four dimensional 
        # array of transition probabilities nu_alpha -> nu_i: 
        # the index order is matter/antimatter, energy, i, alpha
        pSN = SNOSHEWS.Run(ID)

        # restructure the results
        Pmf = FlavorMatrix(np.zeros((6,6,ID.NE)))

        for m in range(ID.NE):
            Pmf[0,ThreeFlavor.NU_E,m] = pSN[0][m][0][0] 
            Pmf[1,ThreeFlavor.NU_E,m] = pSN[0][m][1][0]
            Pmf[2,ThreeFlavor.NU_E,m] = pSN[0][m][2][0]
            Pmf[0,ThreeFlavor.NU_MU,m] = pSN[0][m][0][1] 
            Pmf[1,ThreeFlavor.NU_MU,m] = pSN[0][m][1][1]
            Pmf[2,ThreeFlavor.NU_MU,m] = pSN[0][m][2][1]
            Pmf[0,ThreeFlavor.NU_TAU,m] = pSN[0][m][0][2] 
            Pmf[1,ThreeFlavor.NU_TAU,m] = pSN[0][m][1][2]
            Pmf[2,ThreeFlavor.NU_TAU,m] = pSN[0][m][2][2]

            Pmf[3,ThreeFlavor.NU_E_BAR,m] = pSN[1][m][0][0] 
            Pmf[4,ThreeFlavor.NU_E_BAR,m] = pSN[1][m][1][0]
            Pmf[5,ThreeFlavor.NU_E_BAR,m] = pSN[1][m][2][0]
            Pmf[3,ThreeFlavor.NU_MU_BAR,m] = pSN[1][m][0][1] 
            Pmf[4,ThreeFlavor.NU_MU_BAR,m] = pSN[1][m][1][1]
            Pmf[5,ThreeFlavor.NU_MU_BAR,m] = pSN[1][m][2][1]
            Pmf[3,ThreeFlavor.NU_TAU_BAR,m] = pSN[1][m][0][2] 
            Pmf[4,ThreeFlavor.NU_TAU_BAR,m] = pSN[1][m][1][2]
            Pmf[5,ThreeFlavor.NU_TAU_BAR,m] = pSN[1][m][2][2]
            
        return Pmf
        
###############################################################################

class AdiabaticMSWes(SNTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large 
       that the electron-sterile resonances are inside the neutrinosphere;
    2. the outer electron-sterile MSW resonance is adiabatic;
    3. the inner electron-sterile MSW resonance (where the electron 
       fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, 
        Phys. Rev. D 90, 033013 (2014).
    """
    def P_mf(self, t, E):
        Pmf = self.mix_pars.Pmf_HighDensityLimit()
        return Pmf
        
###############################################################################

class NonAdiabaticMSWes(SNTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large 
       that the electron-sterile resonances are inside the neutrinosphere;
    2. the outer electron-sterile MSW resonance is non-adiabatic;
    3. the inner electron-sterile MSW resonance (where the electron 
       fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, 
        Phys. Rev. D 90, 033013 (2014).
    """
    def P_mf(self, t, E): 
        Pmf = self.mix_pars.Pmf_HighDensityLimit()
        if self.mix_pars.mass_order == MassHierarchy.NORMAL:
            Pmf[['NU_3','NU_4'],:] = Pmf[['NU_4','NU_3'],:]
            Pmf[['NU_2_BAR','NU_3_BAR','NU_4_BAR'],:] = Pmf[['NU_3_BAR','NU_4_BAR','NU_2_BAR'],:]
        else:
            Pmf[['NU_2','NU_4'],:] = Pmf[['NU_4','NU_2'],:]
            Pmf[['NU_1_BAR','NU_2_BAR','NU_4_BAR'],:] = Pmf[['NU_2_BAR','NU_4_BAR','NU_1_BAR'],:]

        return Pmf

###############################################################################
# Vacuum transformations
###############################################################################
class VacuumTransformation(ABC):
    @abstractmethod
    def P_mm(self, t, E)->FlavorMatrix:
        pass
###############################################################################
class NoVacuumTransformation(VacuumTransformation):
    def P_mm(self, t, E)->FlavorMatrix:
        return FlavorMatrix.eye(self.mix_pars.basis_mass)
        
class NeutrinoDecay(VacuumTransformation):
    """Decay effect, where the heaviest neutrino decays to the lightest
    neutrino. For a description and typical parameters, see A. de Gouvêa et al.,
    PRD 101:043013, 2020, arXiv:1910.01127.
    """
    def __init__(self, mass=1*u.eV/c.c**2, tau=1*u.day, dist=10*u.kpc):
        """Initialize transformation matrix.

        Parameters
        ----------
        mass : astropy.units.quantity.Quantity
            Mass of the heaviest neutrino; expect in eV/c^2.
        tau : astropy.units.quantity.Quantity
            Lifetime of the heaviest neutrino.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        """
        self.m = mass
        self.tau = tau
        self.d = dist

    def gamma(self, E):
        """Decay width of the heaviest neutrino mass.

        Parameters
        ----------
        E : float
            Energy of the nu3.

        Returns
        -------
        Gamma : float
            Decay width of the neutrino mass, in units of 1/length.

        :meta private:
        """
        return self.m*c.c / (E*self.tau)

    def P_mm(self, t, E)->FlavorMatrix:
        
        decay_factor = np.exp(-self.gamma(E)*self.d)
        PND = FlavorMatrix.eye(self.mix_pars.basis_mass, extra_dims=[len(E)])

        if self.mix_pars.mass_order == MassHierarchy.NORMAL:
            PND['NU_1','NU_3'] = 1 - decay_factor
            PND['NU_3','NU_3'] = decay_factor

        else:
            PND['NU_2','NU_2'] = decay_factor
            PND['NU_3','NU_2'] = 1 - decay_factor

            PND['NU_BAR','NU_BAR'] = PND['NU','NU']
        return PND

###############################################################################

class QuantumDecoherence(VacuumTransformation):
    """Quantum Decoherence, where propagation in vacuum leads to equipartition
    of states. For a description and typical parameters, see M. V. dos Santos et al.,
    2023, arXiv:2306.17591.
    """
    def __init__(self, Gamma3=1e-27*u.eV, Gamma8=1e-27*u.eV, dist=10*u.kpc, n=0, E0=10*u.MeV):
        """Initialize transformation matrix.

        Parameters
        ----------
        Gamma3 : astropy.units.quantity.Quantity
            Quantum decoherence parameter; expect in eV.
        Gamma8 : astropy.units.quantity.Quantity
            Quantum decoherence parameter; expect in eV.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        n : float
            Exponent of power law for energy dependent quantum decoherence parameters,
            i.e. Gamma = Gamma0*(E/E0)**n. If not specified, it is taken as zero.
        E0 : astropy.units.quantity.Quantity
            Reference energy in the power law Gamma = Gamma0*(E/E0)**n. If not specified, 
            it is taken as 10 MeV. Note that if n = 0, quantum decoherence parameters are independent
            of E0.
        """
        self.Gamma3 = (Gamma3 / (c.hbar.to('eV s') * c.c)).to('1/kpc')
        self.Gamma8 = (Gamma8 / (c.hbar.to('eV s') * c.c)).to('1/kpc')
        self.d = dist
        self.n = n
        self.E0 = E0

    def P_mm(self, t, E)->FlavorMatrix: 
        PQD = FlavorMatrix.zeros(self.mix_pars.basis_mass, extra_dims=(len(E)))

        PQD[1,1] = 1/3 + 1/2 * np.exp(-(self.Gamma3 * (E/self.E0)**self.n + self.Gamma8 * (E/self.E0)**self.n / 3) * self.d) \
                  + 1/6 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        PQD[1,2] = 1/3 - 1/2 * np.exp(-(self.Gamma3 * (E/self.E0)**self.n + self.Gamma8 * (E/self.E0)**self.n / 3) * self.d) \
                  + 1/6 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        PQD[1,3] = 1/3 - 1/3 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        PQD[2,[1,2,3]] = PQD[1,[2,1,3]]
        PQD[2,2] = PQD[1,1]
        PQD[2,3] = PQD[1,3]

        PQD[3,1] = PQD[1,3]
        PQD[3,2] = PQD[2,3]
    
        PQD[3,3] = 1/3 + 2/3 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)
        PQD['NU_BAR','NU_BAR'] = PQD['NU','NU']
        return PQD

###############################################################################
# Earth transformations
###############################################################################

class EarthTransformation(ABC):
    @abstractmethod
    def P_fm(self, t, E)->FlavorMatrix:
        pass
###############################################################################

class NoEarthMatter(EarthTransformation):
    def P_fm(self, t, E)->FlavorMatrix:
        D = self.mix_pars.VacuumMixingMatrix().abs2()
        return D
###############################################################################

try:
    import EMEWS
except:
    EMEWS = None
        
class EarthMatter(EarthTransformation):

    def __init__(self, SNAltAz):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        SNAltAz : astropy AltAz object
        """
           
        self.SNAltAz = SNAltAz

        self.prior_E = None # used to store energy array from previous calls to get_probabilities
        self.prior_D = None

    def P_fm(self, t, E):
        """the D matrix for the case of Earth matter effects

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        D : an array of length of the E array with each element being a 6 x 6 matrix
        """
        if EMEWS == None:
            print("The EMEWS module cannot be found. Results do not include the Earth-matter effect.")
            return ThreeFlavorNoEarthMatter(self.mix_pars).get_probabilities(t,E)

        if self.prior_E != None:
            if u.allclose(self.prior_E, E) == True:
                return self.prior_D

        self.prior_D = np.zeros((6,6,len(E))) 
        self.prior_E = E
        
        #input data object for EMEWS
        ID = EMEWS.InputDataEMEWS()

        ID.altitude = self.SNAltAz.alt.deg
        ID.azimuth = self.SNAltAz.az.deg

        ID.outputfilenamestem = "./out/EMEWS:PREM"   # stem of output filenames 
        ID.densityprofile = "./PREM.rho.dat"         # PREM density profile    
        ID.electronfraction = "./PREM.Ye.dat"        # Electron fraction of density profile     
        
        ID.NE = len(E)         # number of energy bins
        E = E.to_value('MeV')
        ID.Emin = E[0]         # in MeV
        ID.Emax = E[-1]        # in MeV

        #MixingParameters
        ID.deltam_21 = self.mix_pars.dm21_2.value   # in eV^2
        ID.deltam_32 = self.mix_pars.dm32_2.value   # in eV^2
        ID.theta12 = self.mix_pars.theta12.value    # in degrees
        ID.theta13 = self.mix_pars.theta13.value    # in degrees
        ID.theta23 = self.mix_pars.theta23.value    # in degrees
        ID.deltaCP = self.mix_pars.deltaCP.value    # in degrees

        ID.accuracy = 1.01E-007      # controls accuracy of integrtaor: smaller is more accurate
        ID.stepcounterlimit = 1      # output frequency if outputflag = True: larger is less frequent
        ID.outputflag = False        # set to True if output is desired
 
        #matrix from EMEWS needs to be rearranged to match SNEWPY flavor indicii ordering
        Pfm = EMEWS.Run(ID)

        for m in range(len(E)):
            self.prior_D[ThreeFlavor.NU_E,0,m] = Pfm[0][m][0][0] 
            self.prior_D[ThreeFlavor.NU_E,1,m] = Pfm[0][m][0][1]
            self.prior_D[ThreeFlavor.NU_E,2,m] = Pfm[0][m][0][2]
            self.prior_D[ThreeFlavor.NU_MU,0,m] = Pfm[0][m][1][0] 
            self.prior_D[ThreeFlavor.NU_MU,1,m] = Pfm[0][m][1][1]
            self.prior_D[ThreeFlavor.NU_MU,2,m] = Pfm[0][m][1][2]
            self.prior_D[ThreeFlavor.NU_TAU,0,m] = Pfm[0][m][2][0] 
            self.prior_D[ThreeFlavor.NU_TAU,1,m] = Pfm[0][m][2][1]
            self.prior_D[ThreeFlavor.NU_TAU,2,m] = Pfm[0][m][2][2]

            self.prior_D[ThreeFlavor.NU_E_BAR,3,m] = Pfm[1][m][0][0] 
            self.prior_D[ThreeFlavor.NU_E_BAR,4,m] = Pfm[1][m][0][1]
            self.prior_D[ThreeFlavor.NU_E_BAR,5,m] = Pfm[1][m][0][2]
            self.prior_D[ThreeFlavor.NU_MU_BAR,3,m] = Pfm[1][m][1][0] 
            self.prior_D[ThreeFlavor.NU_MU_BAR,4,m] = Pfm[1][m][1][1]
            self.prior_D[ThreeFlavor.NU_MU_BAR,5,m] = Pfm[1][m][1][2]
            self.prior_D[ThreeFlavor.NU_TAU_BAR,3,m] = Pfm[1][m][2][0] 
            self.prior_D[ThreeFlavor.NU_TAU_BAR,4,m] = Pfm[1][m][2][1]
            self.prior_D[ThreeFlavor.NU_TAU_BAR,5,m] = Pfm[1][m][2][2]

        return self.prior_D

#####################################################################
class TransformationChain(FlavorTransformation):
    def __init__(self, 
                 in_sn: SNTransformation,
                 in_vacuum: VacuumTransformation=NoVacuumTransformation(),
                 in_earth: EarthTransformation=NoEarthMatter(),
                 *, mix_pars=MixingParameters()):
        if in_sn is None:
            return NoTransformation()
        self.in_sn = in_sn
        self.in_vacuum = in_vacuum
        self.in_earth = in_earth
        self.transforms = [in_sn, in_vacuum, in_earth]

        self.mix_pars = mix_pars
        #set the mixing parameters to all the inner classes
        for t in self.transforms:
            t.mix_pars = mix_pars
        
    def P_ff(self, t, E)->FlavorMatrix:
        in_sn, in_vacuum, in_earth = self.transforms
        return in_earth.P_fm(t,E) @ in_vacuum.P_mm(t,E) @ in_sn.P_mf(t,E)
    
    def __str__(self):
        s = '+'.join([t.__class__.__name__ for t in self.chain])+'_'+self.mix_pars.mass_order.name
        return s
