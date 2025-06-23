r"""
Transformations in Earth
=============================
Transition of neutrino mass states to flavor states :math:`\nu_i\to\nu_\alpha` in the Earth matter
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u

from snewpy.flavor  import FlavorMatrix
from .base import ThreeFlavorTransformation, FourFlavorTransformation
from snewpy.neutrino import MixingParameters

from importlib.resources import files

try:
    import BEMEWS
    import BEMEWS.data
except ModuleNotFoundError:
    BEMEWS = None


###############################################################################
# Earth transformations
###############################################################################
class EarthTransformation(ABC):
    @abstractmethod
    def P_fm(self, t, E)->FlavorMatrix:
        pass
###############################################################################

class NoEarthMatter(EarthTransformation):
    def __init__(self, mixing_params=None):
         self.mixing_params = mixing_params or MixingParameters('NORMAL')
         
    def P_fm(self, t, E)->FlavorMatrix:
        D = self.mixing_params.VacuumMixingMatrix().abs2()
        return D
###############################################################################

class EarthMatter(ThreeFlavorTransformation, EarthTransformation):

    def __init__(self, SNAltAz, mixing_params=None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mixing_params : ThreeFlavorMixingParameters instance or None
        SNAltAz : astropy AltAz object
        """       
        if BEMEWS is None:
            raise ModuleNotFoundError('BEMEWS module is not found. Please make sure BEMEWS is installed to use EarthMatter transformation')
        if(mixing_params):
            self.mixing_params=mixing_params
        self.SNAltAz = SNAltAz

        self.prior_E = None # used to store energy array from previous calls to get_probabilities
        self.prior_D = None

        
        # Initialize BEMEWS input data object
        self.settings = BEMEWS.InputDataBEMEWS()

        self.settings.altitude = self.SNAltAz.alt.deg
        self.settings.azimuth = self.SNAltAz.az.deg

        self.settings.densityprofile = str(files(BEMEWS.data).joinpath('PREM.rho.dat'))
        self.settings.electronfraction = str(files(BEMEWS.data).joinpath('PREM.Ye.dat'))
        self.settings.accuracy = 1.01e-9
        self.settings.outputflag = False
        self.settings.stepcounterlimit = False
        super().__init__()

    def _update_settings(self):
        """Put the values from mixing_parameters into self.settings"""
        self.settings.deltam_21 = self.mixing_params.dm21_2.to_value('eV**2')
        self.settings.deltam_32 = self.mixing_params.dm32_2.to_value('eV**2')
        self.settings.theta12 = self.mixing_params.theta12.to_value('deg')
        self.settings.theta13 = self.mixing_params.theta13.to_value('deg')
        self.settings.theta23 = self.mixing_params.theta23.to_value('deg')
        self.settings.deltaCP = self.mixing_params.deltaCP.to_value('deg')
            
    def P_fm(self, t, E):
        #update the settings - in case mixing_params were changed
        self._update_settings()
        
        if self.prior_E != None:
            # Use cached result if possible
            if u.allclose(self.prior_E, E) == True:
                return self.prior_D

        self.prior_E = E
        
        #- Set the input energy bins
        E = E.to_value('MeV')
        self.settings.NE = len(E)
        self.settings.Emin = E[0]
        self.settings.Emax = E[-1]

        #run the calculation
        Pfm = np.asarray(BEMEWS.Run(self.settings))
        #matrix from BEMEWS needs to be rearranged to match SNEWPY flavor indicii ordering
        #Pfm contains P(nu_alpha -> nu_i) index order is (nu/nubar, energy, alpha, i)
        #We convert the array dimensions: 
        Pfm = np.swapaxes(Pfm, 1,3) #(nu/nubar, i, alpha, energy)
    
        P = FlavorMatrix.zeros(
            flavor=self.mixing_params.basis_mass,
            from_flavor=self.mixing_params.basis_flavor,
            extra_dims=E.shape)

        P["NU","NU"] = Pfm[0]
        P["NU_BAR","NU_BAR"] = Pfm[1]
        self.prior_D = P
        return P
