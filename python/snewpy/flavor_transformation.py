# -*- coding: utf-8 -*-
"""Supernova oscillation physics for flavors e, X, e-bar, X-bar.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import AltAz

from .neutrino import MassHierarchy, ThreeFlavor, FourFlavor

###############################################################################

class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    @abstractmethod
    def get_probabilities(self, t, E):
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : a N x N array, or an array of N x N arrays 
            where N is either 6 or 8
        """    
        pass

###############################################################################

class ThreeFlavorTransformation(FlavorTransformation):
    """Base class defining common data and methods for all three flavor transformations"""

    def __init__(self, mix_params):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        """
        if mix_params == None:
            self.mix_params = ThreeFlavorMixingParameters(MassHierarchy.NORMAL)
        else:
            self.mix_params = mix_params


    def Pmf_HighDensityLimit(self):
        """ The probability that a given flavor state is 'detected' as a particular matter state in the 
        infinite density limit.  

        Returns
        -------
        6 x 6 matrix
        """
        PmfHDL = np.zeros((6,6))

        NU_E, NU_MU, NU_TAU, NU_E_BAR, NU_MU_BAR, NU_TAU_BAR = \
             ThreeFlavor.NU_E, ThreeFlavor.NU_MU, ThreeFlavor.NU_TAU, \
             ThreeFlavor.NU_E_BAR, ThreeFlavor.NU_MU_BAR, ThreeFlavor.NU_TAU_BAR

        M2 = np.zeros((6,6)) 
        M2[1,1] = self.mix_params.dm21_2.value 
        M2[2,2] = self.mix_params.dm31_2.value 
        M2[1+3,1+3] = self.mix_params.dm21_2.value 
        M2[2+3,2+3] = self.mix_params.dm31_2.value 

        U = self.mix_params.VacuumMixingMatrix()

        HV = U @ M2 @ np.conjugate(np.transpose(U))

        T = np.real( ( HV[NU_MU,NU_MU] + HV[NU_TAU,NU_TAU] ) / 2 )
        D = np.abs( HV[NU_MU,NU_TAU] )**2 - np.abs( HV[NU_MU,NU_MU]-HV[NU_TAU,NU_TAU] )**2 / 4 

        Tbar = np.real( ( HV[NU_MU_BAR,NU_MU_BAR] + HV[NU_TAU_BAR,NU_TAU_BAR] ) / 2 )
        Dbar = np.abs( HV[NU_MU_BAR,NU_TAU_BAR] )**2 \
              -np.abs( HV[NU_MU_BAR,NU_MU_BAR]-HV[NU_TAU_BAR,NU_TAU_BAR] )**2 / 4 

        """ The NMO case. Matter state 3 is the electron neutrino, matter state 1bar is the electron 
        antineutrino. How mu and tau flavors are assigned to matter states 1 and 2 (2bar and 3bar 
        for the antineutrinos) is given by 
        """
        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            k1 = T - np.sqrt(D)
            k2 = T + np.sqrt(D)
            k2bar = Tbar - np.sqrt(Dbar)
            k3bar = Tbar + np.sqrt(Dbar)

            PmfHDL[2,NU_E] = 1 
            PmfHDL[0,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k1 ) / ( k2 - k1 )
            PmfHDL[1,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k2 ) / ( k1 - k2 )
            PmfHDL[0,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k1 ) / ( k2 - k1 )
            PmfHDL[1,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k2 ) / ( k1 - k2 )

            PmfHDL[3,NU_E_BAR] = 1
            PmfHDL[4,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k2bar ) / ( k3bar - k2bar )
            PmfHDL[5,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k3bar ) / ( k2bar - k3bar )
            PmfHDL[4,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k2bar ) / ( k3bar - k2bar )
            PmfHDL[5,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k3bar ) / ( k2bar - k3bar )

        """ The IMO case. Matter state 2 is the electron neutrino, matter state 3bar is the electron 
        antineutrino. How mu and tau flavors are assigned to matter states 1 and 3 (1bar and 2bar 
        for the antineutrinos) is given by 
        """
        if self.mix_params.mass_order == MassHierarchy.INVERTED:
            k1 = T + np.sqrt(D)
            k3 = T - np.sqrt(D)
            k1bar = Tbar - np.sqrt(Dbar)
            k2bar = Tbar + np.sqrt(Dbar)

            PmfHDL[1,NU_E] = 1 
            PmfHDL[0,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k1 ) / ( k3 - k1 )
            PmfHDL[2,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k3 ) / ( k1 - k3 )
            PmfHDL[0,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k1 ) / ( k3 - k1 )
            PmfHDL[2,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k3 ) / ( k1 - k3 )

            PmfHDL[5,NU_E_BAR] = 1
            PmfHDL[3,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k1bar ) / ( k2bar - k1bar )
            PmfHDL[4,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k2bar ) / ( k1bar - k2bar )
            PmfHDL[3,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k1bar ) / ( k2bar - k1bar )
            PmfHDL[4,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k2bar ) / ( k1bar - k2bar )

        return PmfHDL

###############################################################################

class FourFlavorTransformation:
    """Base class defining common data and method for all four flavor transformations"""

    def __init__(self, mix_params):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : FourFlavorMixingParameters instance
        """
        if mix_params == None:
            self.mix_params = FourFlavorMixingParameters(MassHierarchy.NORMAL)
        else:
            self.mix_params = mix_params

    def Pmf_HighDensityLimit(self):
        """ The probability that a given flavor state is 'detected' as a particular matter state in the 
        infinite density limit. This method assumes the 4th matter state is the 
        heaviest and that the electron fraction remains larger than 1/3

        Returns
        -------
        8 x 8 matrix
        """
        PmfHDL = np.zeros((8,8))

        NU_E, NU_MU, NU_TAU, NU_S, NU_E_BAR, NU_MU_BAR, NU_TAU_BAR, NU_S_BAR = \
             FourFlavor.NU_E, FourFlavor.NU_MU, FourFlavor.NU_TAU, FourFlavor.NU_S, \
             FourFlavor.NU_E_BAR, FourFlavor.NU_MU_BAR, FourFlavor.NU_TAU_BAR, FourFlavor.NU_S_BAR

        M2 = np.zeros((8,8)) 
        M2[1,1] = self.mix_params.dm21_2.value 
        M2[2,2] = self.mix_params.dm31_2.value 
        M2[3,3] = self.mix_params.dm41_2.value 
        M2[1+4,1+4] = self.mix_params.dm21_2.value 
        M2[2+4,2+4] = self.mix_params.dm31_2.value 
        M2[3+4,3+4] = self.mix_params.dm41_2.value 

        U = VacuumMixingMatrix()

        HV = np.zeros((6,6,len(E)))        
        for m in range(len(E)):
            HV[:,:,m] = U @ M2 @ np.conjugate(np.transpose(U)) / ( 2 * E[m] )

        T = np.real( (HV[NU_MU,NU_MU] + HV[NU_TAU,NU_TAU]) / 2 )
        D = np.abs( HV[NU_MU,NU_TAU] )**2 - np.abs( HV[NU_MU,NU_MU]-HV[NU_TAU,NU_TAU] )**2 / 4 

        Tbar = np.real( ( HV[NU_MU_BAR,NU_MU_BAR] + HV[NU_TAU_BAR,NU_TAU_BAR] ) / 2 )
        Dbar = np.abs( HV[NU_MU_BAR,NU_TAU_BAR] )**2 \
              -np.abs( HV[NU_MU_BAR,NU_MU_BAR]-HV[NU_TAU_BAR,NU_TAU_BAR] )**2 / 4 

        U = Vacuum_Mixing_Matrix()

        """ The NMO case. Matter state 4 is the electron neutrino, matter state 1bar is the electron 
        antineutrino. Matter state 3 is the sterile neutrino, matter state 2bar is the sterile 
        antineutrino. How mu and tau flavors are assigned to matter states 1 and 2 (3bar and 4bar 
        for the antineutrinos) is given by the transpose of the vacuum mixing matrix
        """
        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            k1 = T - np.sqrt(D)
            k2 = T + np.sqrt(D)
            k3bar = Tbar - np.sqrt(Dbar)
            k4bar = Tbar + np.sqrt(Dbar)

            PmfHDL[3,NU_E] = 1 
            PmfHDL[0,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k1 ) / ( k2 - k1 )
            PmfHDL[1,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k2 ) / ( k1 - k2 )
            PmfHDL[0,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k1 ) / ( k2 - k1 )
            PmfHDL[1,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k2 ) / ( k1 - k2 )
            PmfHDL[2,NU_S] = 1

            PmfHDL[4,NU_E_BAR] = 1
            PmfHDL[6,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k3bar ) / ( k4bar - k3bar )
            PmfHDL[7,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k4bar ) / ( k3bar - k4bar )
            PmfHDL[6,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k3bar ) / ( k4bar - k3bar )
            PmfHDL[7,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k4bar ) / ( k3bar - k4bar )
            PmfHDL[5,NU_S_BAR] = 1

        """ The IMO case. Matter state 4 is the electron neutrino, matter state 3bar is the electron 
        antineutrino. Matter state 2 is the sterile neutrino, matter state 1bar is the sterile 
        antineutrino. How mu and tau flavors are assigned to matter states 1 and 3 (2bar and 4bar 
        for the antineutrinos) is given by the transpose of the vacuum mixing matrix
        """
        if self.mix_params.mass_order == MassHierarchy.INVERTED:
            k1 = T + np.sqrt(D)
            k3 = T - np.sqrt(D)
            k2bar = Tbar - np.sqrt(Dbar)
            k4bar = Tbar + np.sqrt(Dbar)

            PmfHDL[3,NU_E] = 1 
            PmfHDL[0,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k1 ) / ( k3 - k1 )
            PmfHDL[2,NU_MU] = ( np.real(HV[NU_TAU,NU_TAU]) - k3 ) / ( k1 - k3 )
            PmfHDL[0,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k1 ) / ( k3 - k1 )
            PmfHDL[2,NU_TAU] = ( np.real(HV[NU_MU,NU_MU]) - k3 ) / ( k1 - k3 )
            PmfHDL[1,NU_S] = 1 

            PmfHDL[6,NU_E_BAR] = 1
            PmfHDL[5,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k2bar ) / ( k4bar - k2bar )
            PmfHDL[7,NU_MU_BAR] = ( np.real(HV[NU_TAU_BAR,NU_TAU_BAR]) - k4bar ) / ( k2bar - k4bar )
            PmfHDL[5,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k2bar ) / ( k4bar - k2bar )
            PmfHDL[7,NU_TAU_BAR] = ( np.real(HV[NU_MU_BAR,NU_MU_BAR]) - k4bar ) / ( k2bar - k4bar )
            PmfHDL[4,NU_S_BAR] = 1


        return PmfHDL

###############################################################################

class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def __init__(self):
        pass

    def __str__(self):
        return f'NoTransformation'

    def get_probabilities(self, t, E):
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    
        p = np.zeros((6,6,len(E)))
        for f in ThreeFlavor:
            p[f,f] = np.ones(len(E))

        return p

###############################################################################

class CompleteExchange(FlavorTransformation):
    """Survival probabilities for the case when the electron flavors 
       are half exchanged with the mu flavors and the half with the tau flavors.
    """

    def __init__(self):
        pass

    def __str__(self):
        return f'CompleteExchange'

    def get_probabilities(self, t, E):
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    
        p = np.zeros((6,6,len(E))) 

        p[ThreeFlavor.NU_E,ThreeFlavor.NU_MU] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_E,ThreeFlavor.NU_TAU] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_MU,ThreeFlavor.NU_E] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_MU,ThreeFlavor.NU_MU] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_TAU,ThreeFlavor.NU_E] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_TAU,ThreeFlavor.NU_TAU] = np.ones(len(E)) / 2

        p[ThreeFlavor.NU_E_BAR,ThreeFlavor.NU_MU_BAR] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_E_BAR,ThreeFlavor.NU_TAU_BAR] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_MU_BAR,ThreeFlavor.NU_E_BAR] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_MU_BAR,ThreeFlavor.NU_MU_BAR] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_TAU_BAR,ThreeFlavor.NU_E_BAR] = np.ones(len(E)) / 2
        p[ThreeFlavor.NU_TAU_BAR,ThreeFlavor.NU_TAU_BAR] = np.ones(len(E)) / 2

        return p        
        
###############################################################################

class AdiabaticMSW(ThreeFlavorTransformation):
    """Adiabatic MSW effect."""

    def __init__(self, mix_params = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)        

    def __str__(self):
        return f'AdiabaticMSW_' + str(self.mix_params.mass_order)

    def get_SNprobabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities in the supernova

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        Pmf : 6 x 6 array
        """    
        Pmf = self.Pmf_HighDensityLimit()

        return Pmf

    def get_probabilities(self, t, E):         
        """neutrino and antineutrino transition probabilities from the 
        neutrinosphere to Earth

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    
        Pmf = self.get_SNprobabilities(t,E)
        D = ThreeFlavorNoEarthMatter(self.mix_params).get_probabilities(t,E)

        p = np.empty((6,6,len(E))) 
        for m in range(len(E)):
            p[:,:,m] = D[:,:,m] @ Pmf[:,:]

        return p
        
###############################################################################

class NonAdiabaticMSWH(ThreeFlavorTransformation):
    """Nonadiabatic MSW effect."""

    def __init__(self, mix_params = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)  

    def __str__(self):
        return f'NonAdiabaticMSWH_' + str(self.mix_params.mass_order)

    def get_SNprobabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities in the supernova

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        Pmf : 6 x 6 array
        """            
        PmfHDL = self.Pmf_HighDensityLimit()

        Pmf = PmfHDL
        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            for f in ThreeFlavor:
                Pmf[1,f], Pmf[2,f] = Pmf[2,f], Pmf[1,f]
        if self.mix_params.mass_order == MassHierarchy.INVERTED:
            for f in ThreeFlavor:
                Pmf[3,f], Pmf[5,f] = Pmf[5,f], Pmf[3,f]

        return Pmf

    def get_probabilities(self, t, E):         
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.
            
        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    

        Pmf = self.get_SNprobabilities(t,E)      
        D = ThreeFlavorNoEarthMatter(self.mix_params).get_probabilities(t,E)        
        
        p = np.empty((6,6,len(E))) 
        for m in range(len(E)):
            p[:,:,m] = D[:,:,m] @ Pmf[:,:]
        
        return p

###############################################################################

class TwoFlavorDecoherence(ThreeFlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def __init__(self, mix_params = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)    

    def __str__(self):
        return f'TwoFlavorDecoherence_' + str(self.mix_params.mass_order)

    def get_SNprobabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    
        PmfHDL = self.Pmf_HighDensityLimit(E)  

        Pmf = PmfHDL
        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            for f in ThreeFlavor:
                Pmf[1,f], Pmf[2,f] = 0.5 * ( Pmf[1,f] + Pmf[2,f] ), 0.5 * ( Pmf[1,f] + Pmf[2,f] )
        if self.mix_params.mass_order == MassHierarchy.INVERTED:
            for f in ThreeFlavor:
                Pmf[3,f], Pmf[5,f] = 0.5 * ( Pmf[3,f] + Pmf[5,f] ), 0.5 * ( Pmf[3,f] + Pmf[5,f] )

        return Pmf

    def get_probabilities(self, t, E):         
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.
            
        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    

        Pmf = self.get_SNprobabilities(t,E)      
        D = ThreeFlavorNoEarthMatter(self.mix_params).get_probabilities(t,E)        
        
        p = np.empty((6,6,len(E))) 
        for m in range(len(E)):
            p[:,:,m] = D[:,:,m] @ Pmf[:,:]
        
        return p

###############################################################################

class ThreeFlavorDecoherence(ThreeFlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def __init__(self):
        """Initialize ThreeFlavorTransformation to default case"""
        super().__init__(None) 

    def __str__(self):
        return f'ThreeFlavorDecoherence'

    def get_SNprobabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """    
        p = np.ones((6,6,len(E))) / 3

        return p

    def get_probabilities(self, t, E): 
        return self.get_SNprobabilities(t,E) 

###############################################################################

try:
    import SNOSHEWS
except:
    SNOSHEWS = None

class SNprofile:
    """Simple container for density and electron fraction profiles"""

    def __init__(self, rhofilename, Yefilename):         
        self.rhofilename = rhofilename
        self.Yefilename = Yefilename
        
class MSWEffect(ThreeFlavorTransformation):
    """The MSW effect using a density profile and electron 
       fraction provided by the user. Uses the SNOSHEWS module.
    """

    def __init__(self, SNprofile, mix_params = None, rmin = None, rmax = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        SNprofile : instance of profile class
        mix_params : ThreeFlavorMixingParameters instance or None
        """
        self.SNprofile = SNprofile
     
        super().__init__(mix_params) 

        # rmmin will be corrected by SNOSHEWS to the minimum radius of the profile if rmin is less than that value        
        if rmin == None:
            rmin = 0 
        else: 
            self.rmin = rmin

        # rmax will be corrected by SNOSHEWS to the maximum radius of the profile if rmax is greater than that value
        if rmax == None:
            rmax = 1e99 
        else: 
            self.rmax = rmax

        super().__init__(mix_params)        

    def __str__(self):
        return f'MSW_' + str(self.mix_params.mass_order)

    def get_SNprobabilities(self, t, E): 
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
            return AdiabticMSW(self.mix_params).get_probabilities(t,E)
        
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
        ID.deltam_21 = self.mix_params.dm21_2.value   # in eV^2
        ID.deltam_32 = self.mix_params.dm32_2.value   # in eV^2
        ID.theta12 = self.mix_params.theta12.value    # in degrees
        ID.theta13 = self.mix_params.theta13.value    # in degrees
        ID.theta23 = self.mix_params.theta23.value    # in degrees
        ID.deltaCP = self.mix_params.deltaCP.value    # in degrees

        ID.accuracy = 1.01E-009      # controls accuracy of integrator: smaller is more accurate
        ID.stepcounterlimit = 10000  # output frequency if outputflag = True: larger is less frequent
        ID.outputflag = False        # set to True if output is desired
 
        # Do the calculation. The return is a four dimensional 
        # array of transition probabilities nu_alpha -> nu_i: 
        # the index order is matter/antimatter, energy, i, alpha
        pSN = SNOSHEWS.Run(ID)

        # restructure the results
        Pmf = np.zeros((6,6,ID.NE))
        m = 0
        while m<ID.NE:
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
            m += 1
            
    return Pmf

    def get_probabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """
        p = np.empty((6,6,len(E))

        Pmf = self.get_SNprobabilities(t,E)
        D = ThreeFlavorNoEarthMatter(self.mix_params).get_probabilities(t,E)                                     

        # multiply the D matrix and the Pmf matrix together
        for m in range(len(E)):
            p[:,:,m] = self.D[:,:,m] @ Pmf[:,:,m]

        return p
        
###############################################################################

class AdiabaticMSWes(FourFlavorTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large 
       that the electron-sterile resonances are inside the neutrinosphere;
    2. the outer electron-sterile MSW resonance is adiabatic;
    3. the inner electron-sterile MSW resonance (where the electron 
       fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, 
        Phys. Rev. D 90, 033013 (2014).
    """
    def __init__(self, mix_params = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : FourFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)   

    def __str__(self):
        return f'AdiabaticMSWes_' + str(self.mix_params.mass_order)
    
    def get_SNprobabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        Pmf : 8 x 8 matrix
        """   
        Pmf = self.Pmf_HighDensityLimit(E)

        return Pmf

    def get_probabilities(self, t, E):         
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 matrix
        """   
        p = np.empty((8,8,len(E))
        
        D = FourFlavorNoEarthMatter(self.mix_params).get_probabilities(t,E)                                     
        Pmf = self.get_SNprobabilities(t,E)
        
        # multiply the D matrix and the Pmf matrix together
        for m in range(len(E)):
            p[:,:,m] = D[:,:,m] @ Pmf[:,:]

        # remove sterile rows/columns: A is a 6 x 8 matrix
        A = Remove_Steriles()
        p = A @ p @ np.transpose(A)

        return p 
        
###############################################################################

class NonAdiabaticMSWes(FourFlavorTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large 
       that the electron-sterile resonances are inside the neutrinosphere;
    2. the outer electron-sterile MSW resonance is non-adiabatic;
    3. the inner electron-sterile MSW resonance (where the electron 
       fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, 
        Phys. Rev. D 90, 033013 (2014).
    """
    def __init__(self, mix_params = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : FourFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)   

    def __str__(self):
        return f'NonAdiabaticMSWes_' + str(self.mix_params.mass_order)
    
    def get_SNprobabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        Pmf : an 8 x 8 matrix
        """                
        PmfHDL = self.Pmf_HighDensityLimit()
        Pmf = PmfHDL

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            for f in ThreeFlavor:
                Pmf[2,f], Pmf[3,f] = Pmf[3,f], Pmf[2,f]
                Pmf[5,f], Pmf[6,f], Pmf[7,f] = Pmf[6,f], Pmf[7,f], Pmf[5,f]

        if self.mix_params.mass_order == MassHierarchy.INVERTED:
            for f in ThreeFlavor:
                Pmf[1,f], Pmf[3,f] = Pmf[3,f], Pmf[1,f]
                Pmf[4,f], Pmf[5,f], Pmf[7,f] = Pmf[5,f], Pmf[7,f], Pmf[4,f]

        return Pmf

    def get_probabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays       
        """

        D = FourFlavorNoEarthMatter(self.mix_params).get_probabilities(t,E)                                     
        Pmf = self.get_SNprobabilities(t,E)

        p = np.empty((8,8,len(E))) 
        for m in range(len(E)):
            p[:,:,m] = D[:,:,m] @ Pmf[:,:]

        # remove sterile rows/columns
        A = Remove_Steriles()
        p = A @ p @ np.transpose(A)
        
        return p        

###############################################################################

""" The modifier transformtions. These cannot be used by themselves but can be 
    combined with a flavor transformtion class from the first group using the Catenate class
"""

class NeutrinoDecay(ThreeFlavorTransformation):
    """Decay effect, where the heaviest neutrino decays to the lightest
    neutrino. For a description and typical parameters, see A. de Gouvêa et al.,
    PRD 101:043013, 2020, arXiv:1910.01127.
    """
    def __init__(self, mix_params=None, mass=1*u.eV/c.c**2, tau=1*u.day, dist=10*u.kpc):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        mass : astropy.units.quantity.Quantity
            Mass of the heaviest neutrino; expect in eV/c^2.
        tau : astropy.units.quantity.Quantity
            Lifetime of the heaviest neutrino.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        """
        super().__init__(mix_params)        
            
        self.m = mass
        self.tau = tau
        self.d = dist

    def __str__(self):
        return f'NeutrinoDecay_' + str(self.mix_params.mass_order)

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

    def get_probabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        PND : an array of 6 x 6 matrices
        """        
        decay_factor = np.exp(-self.gamma(E)*self.d) 
        PND = np.zeros((6,6,len(E)))  

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            PND[0,0] = 1
            PND[0,2] = 1 - decay_factor
            PND[1,1] = 1
            PND[2,2] = decay_factor

        if self.mix_params.mass_order == MassHierarchy.INVERTED:
            PND[0,0] = 1
            PND[1,1] = decay_factor
            PND[2,1] = 1 - decay_factor
            PND[2,2] = 1

        for i in range(3):
            for j in range(3):
                PND[i+3,j+3] = PND[i,j]        

        return PND

###############################################################################

class QuantumDecoherence(ThreeFlavorTransformation):
    """Quantum Decoherence, where propagation in vacuum leads to equipartition
    of states. For a description and typical parameters, see M. V. dos Santos et al.,
    2023, arXiv:2306.17591.
    """
    def __init__(self, mix_params=None, Gamma3=1e-27*u.eV, Gamma8=1e-27*u.eV, dist=10*u.kpc, n=0, E0=10*u.MeV):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None

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
        super().__init__(mix_params) 

        self.Gamma3 = (Gamma3 / (c.hbar.to('eV s') * c.c)).to('1/kpc')
        self.Gamma8 = (Gamma8 / (c.hbar.to('eV s') * c.c)).to('1/kpc')
        self.d = dist
        self.n = n
        self.E0 = E0

    def __str__(self):
        return f'QuantumDecoherence_' + str(self.mix_params.mass_order)

    def get_probabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        PQD : an array of 6 x 6 matrices
        """        
        PQD = np.zeros((6,6,len(E)))

        PQD[1,1] = 1/3 + 1/2 * np.exp(-(self.Gamma3 * (E/self.E0)**self.n + self.Gamma8 * (E/self.E0)**self.n / 3) * self.d) \ 
                  + 1/6 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        PQD[1,2] = 1/3 - 1/2 * np.exp(-(self.Gamma3 * (E/self.E0)**self.n + self.Gamma8 * (E/self.E0)**self.n / 3) * self.d) \
                  + 1/6 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        PQD[1,3] = 1/3 - 1/3 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        PQD[2,1] = PQD[1,2]
        PQD[2,2] = PQD[1,1]
        PQD[2,3] = PQD[1,3]

        PQD[3,1] = PQD[1,3]
        PQD[3,2] = PQD[2,3]
    
        PQD[3,3] = 1/3 + 2/3 * np.exp(-self.Gamma8 * (E/self.E0)**self.n * self.d)

        for i in range(3):
            for j in range(3):
                PQD[i+3,j+3] = PQD[i,j]

        return PQD

###############################################################################

class ThreeFlavorNoEarthMatter(ThreeFlavorTransformation):

    def __init__(self, mix_params = None ):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)  

    def __str__(self):
        pass

    def get_probabilities(self, t, E):
        """the D matrix for the case of no Earth matter effects and three neutrino flavors

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        an array of length of the E array with each element being the 
        D matrix after computing Earth-matter effects
        """
        U = self.mix_params.VacuumMixingMatrix()

        D = np.zeros((6,6,len(E))) # note the first index is a flavor, the second is a mass state
        for f in ThreeFlavor:
            for i in range(6):
                for m in range(len(E)): 
                    D[f,i,m] = float(np.abs(U[f,i])**2)    
                    
        return D

###############################################################################

class FourFlavorNoEarthMatter(FourFlavorTransformation):

    def __init__(self, mix_params = None ):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : FourFlavorMixingParameters instance or None
        """
        super().__init__(mix_params)  

    def __str__(self):
        pass

    def get_probabilities(self, t, E):
        """the D matrix for the case of no Earth matter effects and four neutrino flavors

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        an array of 8 x 8 matrices of length equal to the length of the E array
        """
        U = self.mix_params.VacuumMixingMatrix()

        D = np.zeros((8,8,len(E))) # note the first index is a flavor, the second is a mass state
        for f in FourFlavor:
            for i in range(8):
                for m in range(len(E)): 
                    D[f,i,m] = float(np.abs(U[f,i])**2)    
                    
        return D

###############################################################################

try:
    import EMEWS
except:
    EMEWS = None
        
class EarthMatter(ThreeFlavorTransformation):

    def __init__(self, SNAltAz, mix_params = None ):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : ThreeFlavorMixingParameters instance or None
        SNAltAz : astropy AltAz object
        """
        super().__init__(mix_params)  
            
        self.SNAltAz = SNAltAz

        self.prior_E = None # used to store energy array from previous calls to get_probabilities

    def __str__(self):
        return f'EarthMatter_' + str(self.mix_params.mass_order)

    def get_probabilities(self, t, E):
        """the D matrix for the case of Earth matter effects

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        an array of length of the E array with each element being the 
        D matrix after computing Earth-matter effects
        """
        if EMEWS == None:
            print("The EMEWS module cannot be found. Results do not include the Earth-matter effect.")
            return self.D

        if self.prior_E != None:
            if u.allclose(self.prior_E, E) == True:
                return self.D

        self.D = np.zeros((6,6,len(E))) 
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
        ID.deltam_21 = self.mix_params.dm21_2.value   # in eV^2
        ID.deltam_32 = self.mix_params.dm32_2.value   # in eV^2
        ID.theta12 = self.mix_params.theta12.value    # in degrees
        ID.theta13 = self.mix_params.theta13.value    # in degrees
        ID.theta23 = self.mix_params.theta23.value    # in degrees
        ID.deltaCP = self.mix_params.deltaCP.value    # in degrees

        ID.accuracy = 1.01E-007      # controls accuracy of integrtaor: smaller is more accurate
        ID.stepcounterlimit = 1      # output frequency if outputflag = True: larger is less frequent
        ID.outputflag = False        # set to True if output is desired
 
        #matrix from EMEWS needs to be rearranged to match SNEWPY flavor indicii ordering
        Pfm = EMEWS.Run(ID)

        m = 0
        while m<len(E):
            self.D[ThreeFlavor.NU_E,0,m] = Pfm[0][m][0][0] 
            self.D[ThreeFlavor.NU_E,1,m] = Pfm[0][m][0][1]
            self.D[ThreeFlavor.NU_E,2,m] = Pfm[0][m][0][2]
            self.D[ThreeFlavor.NU_MU,0,m] = Pfm[0][m][1][0] 
            self.D[ThreeFlavor.NU_MU,1,m] = Pfm[0][m][1][1]
            self.D[ThreeFlavor.NU_MU,2,m] = Pfm[0][m][1][2]
            self.D[ThreeFlavor.NU_TAU,0,m] = Pfm[0][m][2][0] 
            self.D[ThreeFlavor.NU_TAU,1,m] = Pfm[0][m][2][1]
            self.D[ThreeFlavor.NU_TAU,2,m] = Pfm[0][m][2][2]

            self.D[ThreeFlavor.NU_E_BAR,3,m] = Pfm[1][m][0][0] 
            self.D[ThreeFlavor.NU_E_BAR,4,m] = Pfm[1][m][0][1]
            self.D[ThreeFlavor.NU_E_BAR,5,m] = Pfm[1][m][0][2]
            self.D[ThreeFlavor.NU_MU_BAR,3,m] = Pfm[1][m][1][0] 
            self.D[ThreeFlavor.NU_MU_BAR,4,m] = Pfm[1][m][1][1]
            self.D[ThreeFlavor.NU_MU_BAR,5,m] = Pfm[1][m][1][2]
            self.D[ThreeFlavor.NU_TAU_BAR,3,m] = Pfm[1][m][2][0] 
            self.D[ThreeFlavor.NU_TAU_BAR,4,m] = Pfm[1][m][2][1]
            self.D[ThreeFlavor.NU_TAU_BAR,5,m] = Pfm[1][m][2][2]
            m += 1

        return self.D


class Catenate:
    """Catenate flavor transformation effects together."""

    def __init__(self, SNTransformation, InVacuumTransformation = None, AtEarthTransformation = None):
        """        
        Parameters
        ----------
        SNTransformation : the transformation that occurs in the supernova
        InVacuumTransformation : the transforamtion that occurs in the vacuum
        AtEarthTransformation : the transformation that occurs at Earth
        The order is that SNTransformation is applied first, then InVacuumTransformation, 
        then AtEarthTransformation
        """  
        self.transform1 = SNTransformation
        
        if  InVacuumTransformation == None:
            self.transform2 = NoTransformation(self.transform1.mx_params)
        else:
            self.transform2 = InVacuumTransformation

        if  AtEarthTransformation == None:
            self.transform3 = ThreeFlavorNoEarthMatter(self.transform1.mx_params)
        else:
            self.transform3 = AtEarthTransformation


    def __str__(self):
        return str(self.transform1) + '+' + str(self.transform2) + '+' + str(self.transform3)

    def get_probabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        p : 6 x 6 array or array of 6 x 6 arrays 
        """ 
        p1 = self.transform1.get_SNprobabilities(t,E)
        p2 = self.transform2.get_probabilities(t,E)
        p3 = self.transform2.get_probabilities(t,E)        

        p = np.empty((6,6,len(E)))
        for m in range(len(E)):
            p[:,:,m] = p3[:,:,m] @ p2[:,:,m] @ p1[:,:,m]

        return p   
