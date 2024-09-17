from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import AltAz

from snewpy.flavor  import FlavorMatrix, ThreeFlavor

try:
    import BEMEWS
    import BEMEWS.data
except ImportError as e:
    BEMEWS = None


###############################################################################
# Earth transformations
###############################################################################
class EarthTransformation(ABC):
    def __init__(self, mixing_params):
        self.mixing_params = mixing_params

    @abstractmethod
    def P_fm(self, t, E)->FlavorMatrix:
        pass
###############################################################################

class NoEarthMatter(EarthTransformation):
    def __init__(self):
        super().__init__(None)

    def P_fm(self, t, E)->FlavorMatrix:
        D = self.mixing_params.VacuumMixingMatrix().abs2()
        return D
###############################################################################

class EarthMatter(EarthTransformation):

    def __init__(self, mixing_params, SNAltAz):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mixing_params : ThreeFlavorMixingParameters instance or None
        SNAltAz : astropy AltAz object
        """

        super().__init__(mixing_params)

        self.SNAltAz = SNAltAz

        self.prior_E = None # used to store energy array from previous calls to get_probabilities
        self.prior_D = None

        # Initialize BEMEWS input data object
        self.settings = None

        if BEMEWS is not None:
            self.settings = BEMEWS.InputDataBEMEWS()

            self.settings.altitude = self.SNAltAz.alt.deg
            self.settings.azimuth = self.SNAltAz.az.deg

            self.settings.densityprofile = str(files(BEMEWS.data).joinpath('PREM.rho.dat'))
            self.settings.electronfraction = str(files(BEMEWS.data).joinpath('PREM.Ye.dat'))

            self.settings.deltam_21 = self.mixing_params.dm21_2.to_value('eV**2')
            self.settings.deltam_32 = self.mixing_params.dm32_2.to_value('eV**2')
            self.settings.theta12 = self.mixing_params.theta12.to_value('deg')
            self.settings.theta13 = self.mixing_params.theta13.to_value('deg')
            self.settings.theta23 = self.mixing_params.theta23.to_value('deg')
            self.settings.deltaCP = self.mixing_params.deltaCP.to_value('deg')

            self.settings.accuracy = 1.01e-9
            self.settings.outputflag = False
            self.settings.stepcounterlimit = False

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
        if BEMEWS is None:
            logger.warn('BEMEWS is not found. Not computing Earth-matter effect.')
            return NoEarthMatter(self.mixing_params).P_fm(t, E)

        if self.prior_E != None:
            # Use cached result if possible
            if u.allclose(self.prior_E, E) == True:
                return self.prior_D

        self.prior_D = np.zeros((6,6,len(E))) 
        self.prior_E = E
        
        #- Set the input energy bins
        E = E.to_value('MeV')
        self.settings.NE = len(E)
        self.settings.Emin = E[0]
        self.settings.Emax = E[-1]

        #matrix from EMEWS needs to be rearranged to match SNEWPY flavor indicii ordering
        Pfm = BEMEWS.Run(self.settings)

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
