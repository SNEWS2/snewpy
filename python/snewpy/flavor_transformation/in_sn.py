r"""
Transformations in supernova
=============================
Transition of neutrino flavor states to mass states :math:`\nu_\alpha\to\nu_i` inside supernova
"""
from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c

from snewpy.flavor  import FlavorMatrix, ThreeFlavor
from snewpy.neutrino import MixingParameters, MassHierarchy
from dataclasses import dataclass

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
    """Equal mixing of whatever two matter states form the MSW H resonance.
    
    The TwoFlavorDecoherence transformation is relevant when the size of the density fluctuations
    is ≲ 10% for densities around the H resonance density —see Kneller (2010); Kneller & Mauney (2013). 
    This prescription models the effect of the turbulence as leading to 50% mixing between the
    matter states which participate in the H resonance. 
    
    In the NMO this is ν2 and ν3.
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
class MSWEffect(SNTransformation, ThreeFlavorTransformation):
    """The MSW effect using a density profile and electron 
       fraction provided by the user. Uses the SNOSHEWS module.
    """
    @dataclass
    class SNprofile:
        """A placeholder class for the density and electron fraction profiles. Currently SNOSHEWS 
        reads the profiles from files but that might change in the future.
        """
        rhofilename:str
        Yefilename:str
        radius_min:float = 0
        radius_max:float = 1e99
        
    def __init__(self, SNprofile:SNprofile):
        """Initialize flavor transformation
        
        Parameters
        ----------
        SNprofile : instance of profile class
        """
        if SNOSHEWS == None:
            raise ModuleNotFoundError("The SNOSHEWS module not be found. Please make sure SNOSHEWS is installed to use MSWEffect transformation")
        #input data object for SNOSHEWS
        self.settings = SNOSHEWS.InputDataSNOSHEWS()

        self.settings.outputfilenamestem = "./out/SNOSHEWS"      # stem of output filenames 

        self.settings.rmin = SNprofile.radius_min
        self.settings.rmax = SNprofile.radius_max

        self.settings.densityprofile = SNprofile.rhofilename     # mass density profile    
        self.settings.electronfraction = SNprofile.Yefilename    # electron fraction profile     
        
        self.settings.accuracy = 1.01E-009      # controls accuracy of integrator: smaller is more accurate
        self.settings.stepcounterlimit = 10000  # output frequency if outputflag = True: larger is less frequent
        self.settings.outputflag = False        # set to True if output is desired
        self._update_settings()
        
    def _update_settings(self):
        """Put the values from mixing_parameters into self.settings"""
        self.settings.deltam_21 = self.mixing_params.dm21_2.to_value('eV**2')
        self.settings.deltam_32 = self.mixing_params.dm32_2.to_value('eV**2')
        self.settings.theta12 = self.mixing_params.theta12.to_value('deg')
        self.settings.theta13 = self.mixing_params.theta13.to_value('deg')
        self.settings.theta23 = self.mixing_params.theta23.to_value('deg')
        self.settings.deltaCP = self.mixing_params.deltaCP.to_value('deg')
        
    def P_mf(self, t, E): 
        #update the settings - in case mixing_params were changed
        self._update_settings()
        #- Set the input energy bins
        E = E.to_value('MeV')
        self.settings.NE = len(E)         # number of energy bins
        self.settings.Emin = E[0]         # in MeV
        self.settings.Emax = E[-1]        # in MeV
        #run the calculation
        pSN = SNOSHEWS.Run(self.settings)
        #matrix from SNOSHEWS needs to be rearranged to match SNEWPY flavor indicii ordering
        #pSN contains P(nu_i -> nu_alpha) index order is (nu/nubar, energy, i, alpha)
        #We convert the array dimensions: 
        pSN = np.swapaxes(pSN, 1,3) #(nu/nubar, alpha, i, energy)
        # restructure the results
        P = FlavorMatrix.zeros(
            flavor=self.mixing_params.basis_flavor,
            from_flavor=self.mixing_params.basis_mass,
            extra_dims=E.shape)
        
        P["NU","NU"] = pSN[0]
        P["NU_BAR","NU_BAR"] = pSN[1]
        return P
        
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