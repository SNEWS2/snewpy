from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c

from snewpy.flavor  import FlavorMatrix, ThreeFlavor
from snewpy.neutrino import MixingParameters, MassHierarchy

from .base import ThreeFlavorTransformation, FourFlavorTransformation
try:
    import SNOSHEWS
except:
    SNOSHEWS = None


class SNTransformation(ABC):
    """Base class for all transformations in SN"""
    @abstractmethod
    def P_mf(self, t, E)->FlavorMatrix:
        pass

        

###############################################################################
class AdiabaticMSW(SNTransformation, ThreeFlavorTransformation):
    """Adiabatic MSW effect."""

    def P_mf(self, t, E):
        return self.mixing_params.Pmf_HighDensityLimit()
###############################################################################
class NonAdiabaticMSWH(SNTransformation, ThreeFlavorTransformation):
    """Nonadiabatic MSW H resonance. 
    The NonAdiabaticMSWH transformation assumes that the H resonance mixing is nonadiabatic.
    This case is relevant when a shock is present at the H resonance densities (Schirato & Fuller 2002).
    
    For the NMO the H resonance occurs in the neutrinos (Kneller & McLaughlin 2009) between ‘matter’
    states ν2 and ν3.
    In the IMO the H resonance mixes the antineutrino matter states ν̄1 and ν̄3.
    """
    def P_mf(self, t, E):
        Pmf = self.mixing_params.Pmf_HighDensityLimit()
        if self.mixing_params.mass_order == MassHierarchy.NORMAL:
            Pmf[['NU_2','NU_3'],:] = Pmf[['NU_3','NU_2'],:]
        else:
            Pmf[['NU_1_BAR','NU_3_BAR'],:] = Pmf[['NU_3_BAR','NU_1_BAR'],:]
        return Pmf
###############################################################################
class TwoFlavorDecoherence(SNTransformation, ThreeFlavorTransformation):
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
        Pmf = self.mixing_params.Pmf_HighDensityLimit()
        if self.mixing_params.mass_order == MassHierarchy.NORMAL:
            Pmf['NU_2']=Pmf['NU_3']=0.5*(Pmf['NU_2']+Pmf['NU_3'])
        else:
            Pmf['NU_1_BAR']=Pmf['NU_3_BAR']=0.5*(Pmf['NU_1_BAR']+Pmf['NU_3_BAR'])

        return Pmf

###############################################################################

class SNprofile:
    """A placeholder class for the density and electron fraction profiles. Currently SNOSHEWS 
    reads the profiles from files but that might change in the future.
    """

    def __init__(self, rhofilename, Yefilename):         
        self.rhofilename = rhofilename
        self.Yefilename = Yefilename
        
class MSWEffect(SNTransformation, ThreeFlavorTransformation):
    """The MSW effect using a density profile and electron 
       fraction provided by the user. Uses the SNOSHEWS module.
    """

    def __init__(self, SNprofile, rmin = 0, rmax = 1e99):
        """Initialize flavor transformation
        
        Parameters
        ----------
        SNprofile : instance of profile class
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
            return AdiabticMSW(self.mixing_params).get_probabilities(t,E)
        
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
        ID.deltam_21 = self.mixing_params.dm21_2.value   # in eV^2
        ID.deltam_32 = self.mixing_params.dm32_2.value   # in eV^2
        ID.theta12 = self.mixing_params.theta12.value    # in degrees
        ID.theta13 = self.mixing_params.theta13.value    # in degrees
        ID.theta23 = self.mixing_params.theta23.value    # in degrees
        ID.deltaCP = self.mixing_params.deltaCP.value    # in degrees

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

class AdiabaticMSWes(SNTransformation, FourFlavorTransformation):
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
        Pmf = self.mixing_params.Pmf_HighDensityLimit()
        return Pmf
        
###############################################################################

class NonAdiabaticMSWes(SNTransformation, FourFlavorTransformation):
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
        Pmf = self.mixing_params.Pmf_HighDensityLimit()
        if self.mixing_params.mass_order == MassHierarchy.NORMAL:
            Pmf[['NU_3','NU_4'],:] = Pmf[['NU_4','NU_3'],:]
            Pmf[['NU_2_BAR','NU_3_BAR','NU_4_BAR'],:] = Pmf[['NU_3_BAR','NU_4_BAR','NU_2_BAR'],:]
        else:
            Pmf[['NU_2','NU_4'],:] = Pmf[['NU_4','NU_2'],:]
            Pmf[['NU_1_BAR','NU_2_BAR','NU_4_BAR'],:] = Pmf[['NU_2_BAR','NU_4_BAR','NU_1_BAR'],:]

        return Pmf