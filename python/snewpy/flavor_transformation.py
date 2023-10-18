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

from .neutrino import MassHierarchy, Flavor

import Sqa3Earth

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
        p : 4x4 array or array of 4 x 4 arrays 
        """    
        pass

    def prob_ee(self, t, E):
        """Electron neutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_E,Flavor.NU_E]

    def prob_ex(self, t, E):
        """x -> e neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_E,Flavor.NU_X]

    def prob_xx(self, t, E):
        """Flavor x neutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_X,Flavor.NU_X]

    def prob_xe(self, t, E):
        """e -> x neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_X,Flavor.NU_E]

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_E_BAR,Flavor.NU_E_BAR]

    def prob_exbar(self, t, E):
        """x -> e antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_E_BAR,Flavor.NU_X_BAR]

    def prob_xxbar(self, t, E):
        """x -> x antineutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_X_BAR,Flavor.NU_X_BAR]

    def prob_xebar(self, t, E):
        """e -> x antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        float or ndarray
            Transition probability.
        """
        p = get_probabilities(t, E)
        return p[Flavor.NU_X_BAR,Flavor.NU_E_BAR]


class ThreeFlavorTransformation(FlavorTransformation):
    """Base class defining common data and methods for all three flavor transformations"""

    def __init__(self, mix_params = None, AltAz = None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance or None
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done
        """
        if mix_params == None:
            self.mix_params = MixingParameters3Flavor(MassHierarchy.NORMAL)
        else:
            self.mix_params = mix_params
            
        self.AltAz = AltAz

        self.DV_e1 = float((np.cos(mix_params.theta12) * np.cos(mix_params.theta13))**2)
        self.DV_e2 = float((np.sin(mix_params.theta12) * np.cos(mix_params.theta13))**2)
        self.DV_e3 = float((np.sin(mix_params.theta13))**2) 

        self.D_e1 = None
        self.D_e2 = None
        self.D_e3 = None
        self.Dbar_e1 = None
        self.Dbar_e2 = None
        self.Dbar_e3 = None

        self.prior_E = None


    def vacuum(self, E):
        """the first rows of the D and Dbar matrices for the case of no Earth matter effect

        Parameters
        ----------
        E : float or ndarray
            List of energies.

        Returns
        -------
        the six floats or six ndarrays of the first rows of the D and Dbar matrices
        """
        if self.prior_E != None:
            if u.allclose(self.prior_E, E) == True:
                return self.D_e1, self.D_e2, self.D_e3, self.Dbar_e1, self.Dbar_e2, self.Dbar_e3

        if self.prior_E == None:
            self.D_e1 = np.zeros(len(E))
            self.D_e2 = np.zeros(len(E))
            self.D_e3 = np.zeros(len(E))
            self.Dbar_e1 = np.zeros(len(E))
            self.Dbar_e2 = np.zeros(len(E))
            self.Dbar_e3 = np.zeros(len(E))

        self.prior_E = E            

        self.D_e1 = np.full(len(E),self.DV_e1)
        self.D_e2 = np.full(len(E),self.DV_e2)
        self.D_e3 = np.full(len(E),self.DV_e3)

        self.Dbar_e1 = np.full(len(E),self.DV_e1)
        self.Dbar_e2 = np.full(len(E),self.DV_e2)
        self.Dbar_e3 = np.full(len(E),self.DV_e3)

        return self.D_e1, self.D_e2, self.D_e3, self.Dbar_e1, self.Dbar_e2, self.Dbar_e3
        

    def Earth_matter(self, E):
        """the first rows of the D and Dbar matrices for the case of Earth matter effects

        Parameters
        ----------
        E : float or ndarray
            List of energies.

        Returns
        -------
        the six floats or six ndarrays of the first rows of the D and Dbar matrices
        """
        if self.prior_E != None:
            if u.allclose(self.prior_E, E) == True:
                return self.D_e1, self.D_e2, self.D_e3, self.Dbar_e1, self.Dbar_e2, self.Dbar_e3

        if self.prior_E == None:
            self.D_e1 = np.zeros(len(E))
            self.D_e2 = np.zeros(len(E))
            self.D_e3 = np.zeros(len(E))
            self.Dbar_e1 = np.zeros(len(E))
            self.Dbar_e2 = np.zeros(len(E))
            self.Dbar_e3 = np.zeros(len(E))

        self.prior_E = E
         
        E = E.to_value('MeV')

        #input data object for Sqa
        ID = Sqa3Earth.InputDataSqa3Earth()

        ID.altitude = self.AltAz.alt.deg
        ID.azimuth = self.AltAz.az.deg

        ID.outputfilenamestem = "./out/Sqa3Earth:PREM"   # stem of output filenames 
        ID.densityprofile = "./PREM.rho.dat"             # PREM density profile    
        ID.electronfraction = "./PREM.Ye.dat"            # Electron fraction of density profile     
        
        ID.NE = len(E)   # number of energy bins
        ID.Emin = E[0]   # in MeV
        ID.Emax = E[-1]  # in MeV

        #MixingParameters
        ID.deltam_21 = self.mix_params.dm21_2.value   # in eV^2
        ID.deltam_32 = self.mix_params.dm32_2.value  # in eV^2
        ID.theta12 = self.mix_params.theta12.value   # in degrees
        ID.theta13 = self.mix_params.theta13.value    # in degrees
        ID.theta23 = self.mix_params.theta23.value    # in degrees
        ID.deltaCP = self.mix_params.deltaCP.value    # in degrees

        ID.accuracy = 1.01E-007       # controls accuracy of integrtaor: smaller is more accurate
        ID.stepcounterlimit = 10    # output frequency if outputflag = True: larger is less frequent
        ID.outputflag = True         # set to True if output is desired
 
        #matrix from Sqa3Earth needs to be rearranged to match SNEWPY indici
        Pfm = Sqa3Earth.RunSqa3Earth(ID)

        m = 0
        while m<len(E):
            self.D_e1[m] = Pfm[0][m][0][0] 
            self.D_e2[m] = Pfm[0][m][0][1]
            self.D_e3[m] = Pfm[0][m][0][2]

            self.Dbar_e1[m] = Pfm[1][m][0][0]
            self.Dbar_e2[m] = Pfm[1][m][0][1]
            self.Dbar_e3[m] = Pfm[1][m][0][2]
            m += 1

        return self.D_e1, self.D_e2, self.D_e3, self.Dbar_e1, self.Dbar_e2, self.Dbar_e3


class FourFlavorTransformation:
    """Base class defining common data and method for all four flavor transformations"""

    def __init__(self, mix_params):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters4Flavor instance
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done
        """
        if mix_params == None:
            self.mix_params = MixingParameters4Flavor(MassHierarchy.NORMAL)
        else:
            self.mix_params = mix_params

        self.DV_e1 = float((np.cos(mix_params.theta12) * np.cos(mix_params.theta13) * np.cos(mix_params.theta14))**2)
        self.DV_e2 = float((np.sin(mix_params.theta12) * np.cos(mix_params.theta13) * np.cos(mix_params.theta14))**2)
        self.DV_e3 = float((np.sin(mix_params.theta13) * np.cos(mix_params.theta14))**2)
        self.DV_e4 = float((np.sin(mix_params.theta14))**2)  

        self.DV_s1 = float((np.cos(mix_params.theta12) * np.cos(mix_params.theta13) * np.sin(mix_params.theta14))**2)
        self.DV_s2 = float((np.sin(mix_params.theta12) * np.cos(mix_params.theta13) * np.sin(mix_params.theta14))**2)
        self.DV_s3 = float((np.sin(mix_params.theta13) * np.sin(mix_params.theta14))**2)
        self.DV_s4 = float((np.cos(mix_params.theta14))**2)  


    def vacuum(self, E):
        """the first and last rows of the D and Dbar matrices for the case of no Earth matter effect

        Parameters
        ----------
        E : float or ndarray
            List of energies.

        Returns
        -------
        the 16 floats or 16 arrays of the first and last rows of the D and Dbar matrices
        """
        D_e1 = np.full(len(E),self.DV_e1)
        D_e2 = np.full(len(E),self.DV_e2)
        D_e3 = np.full(len(E),self.DV_e3)
        D_e4 = np.full(len(E),self.DV_e4)

        D_s1 = np.full(len(E),self.DV_s1)
        D_s2 = np.full(len(E),self.DV_s2)
        D_s3 = np.full(len(E),self.DV_s3)
        D_s4 = np.full(len(E),self.DV_s4)

        Dbar_e1 = np.full(len(E),self.DV_e1)
        Dbar_e2 = np.full(len(E),self.DV_e2)
        Dbar_e3 = np.full(len(E),self.DV_e3)
        Dbar_e4 = np.full(len(E),self.DV_e4)

        Dbar_s1 = np.full(len(E),self.DV_s1)
        Dbar_s2 = np.full(len(E),self.DV_s2)
        Dbar_s3 = np.full(len(E),self.DV_s3)
        Dbar_s4 = np.full(len(E),self.DV_s4)

        return D_e1, D_e2, D_e3, D_e4, D_s1, D_s2, D_s3, D_s4, Dbar_e1, Dbar_e2, Dbar_e3, Dbar_e4, Dbar_s1, Dbar_s2, Dbar_s3, Dbar_s4


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
        p : 4x4 array or array of 4 x 4 arrays 
        """    
        p_ee = np.ones(len(E))    #p_ee is the probability the an electron neutrino will stay an electron neutrino
        pbar_ee = np.ones(len(E)) #pee_bar is the probability that an electron antineutrino will stay an electron antineutrino

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = 1 - p_ee
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = (1 + p_ee)/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = 1 - pbar_ee
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = (1 + pbar_ee)/2

        return p

    
class CompleteExchange(FlavorTransformation):
    """Survival probabilities for the case when the electron flavors are completely exchanged with the x flavor."""

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
        p : 4x4 array or array of 4 x 4 arrays 
        """    
        p_ee = np.zeros(len(E))    #p_ee is the probability the an electron neutrino will stay an electron neutrino
        pbar_ee = np.zeros(len(E)) #pee_bar is the probability that an electron antineutrino will stay an electron antineutrino

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = 1 - p_ee
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = (1 + p_ee)/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = 1 - pbar_ee
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = (1 + pbar_ee)/2

        return p        
        

class AdiabaticMSW(ThreeFlavorTransformation):
    """Adiabatic MSW effect."""

    def __init__(self, mix_params, AltAz=None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done
        """
        super().__init__(mix_params, AltAz)        

    def __str__(self):
        return f'AdiabaticMSW'

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
        p : 4x4 array or array of 4 x 4 arrays 
        """    
        if self.AltAz == None or self.AltAz.alt > 0:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.vacuum(E) 
        else:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.Earth_matter(E) 

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            p_ee = D_e3
            pbar_ee = Dbar_e1
        else:
            p_ee = D_e2
            pbar_ee = Dbar_e3

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee
        p[Flavor.NU_E,Flavor.NU_X] = 1 - p_ee
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = (1 + p_ee)/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = 1 - pbar_ee
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = (1 + pbar_ee)/2

        return p
        

class NonAdiabaticMSWH(ThreeFlavorTransformation):
    """Nonadiabatic MSW effect."""

    def __init__(self, mix_params, AltAz=None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done
        """
        super().__init__(mix_params, AltAz)   

    def __str__(self):
        return f'NonAdiabaticMSWH'

    def get_probabilities(self, t, E): 
        """neutrino and antineutrino transition probabilities.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.
        """
        if self.AltAz == None or self.AltAz.alt > 0:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.vacuum(E) 
        else:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.Earth_matter(E) 

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            p_ee = D_e2
            pbar_ee = Dbar_e1
        else:
            p_ee = D_e2
            pbar_ee = Dbar_e1
 
        p = np.zeros((4,4,len(E)))

        T[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = 1 - p_ee
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = (1 + p_ee)/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = 1 - pbar_ee
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = (1 + pbar_ee)/2

        return p

class TwoFlavorDecoherence(ThreeFlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def __init__(self, mix_params, AltAz=None):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done
        """
        super().__init__(mix_params, AltAz)       

    def __str__(self):
        return f'TwoFlavorDecoherence'                

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
        p : 4x4 array or array of 4 x 4 arrays 
        """        
        if self.AltAz == None or self.AltAz.alt > 0:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.vacuum(E) 
        else:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.Earth_matter(E) 

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            p_ee = (D_e2 + D_e3)/2
            pbar_ee = Dbar_e1
        else:
            p_ee = D_e2
            pbar_ee = (Dbar_e1 + Dbar_e3)/2

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = 1 - p_ee
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = (1 + p_ee)/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = 1 - pbar_ee
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = (1 + pbar_ee)/2

        return p


class ThreeFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def __init__(self):
        pass

    def __str__(self):
        return f'ThreeFlavorDecoherence'

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
        p : 4x4 array or array of 4 x 4 arrays 
        """    
        p_ee.full(len(E),1/3)
        pbar_ee.full(len(E),1/3)

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = 1 - p_ee
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = (1 + p_ee)/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = 1 - pbar_ee
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = (1 + pbar_ee)/2

        return p


class NeutrinoDecay(ThreeFlavorTransformation):
    """Decay effect, where the heaviest neutrino decays to the lightest
    neutrino. For a description and typical parameters, see A. de Gouvêa et al.,
    PRD 101:043013, 2020, arXiv:1910.01127.
    """
    def __init__(self, mix_params, AltAz=None, mass=1*u.eV/c.c**2, tau=1*u.day, dist=10*u.kpc):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done            
        mass : astropy.units.quantity.Quantity
            Mass of the heaviest neutrino; expect in eV/c^2.
        tau : astropy.units.quantity.Quantity
            Lifetime of the heaviest neutrino.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        AltAz : astropy AltAz object 
            If not None the Earth-matter effect calculation will be done
        """
        super().__init__(mix_params, AltAz)        
            
        self.m = mass
        self.tau = tau
        self.d = dist

    def __str__(self):
        return f'NeutrinoDecay'

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
        p : 4x4 array or array of 4 x 4 arrays 
        """        
        if self.AltAz == None or self.AltAz.alt > 0:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.vacuum(E) 
        else:
            D_e1, D_e2, D_e3, Dbar_e1, Dbar_e2, Dbar_e3 = self.Earth_matter(E) 

        decay_factor = math.e**(-self.gamma(E)*self.dist) 
            
        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            p_ee = D_e1 * (1 - decay_factor) + D_e3 * decay_factor
            p_ex = D_e1 + D_e2
            pbar_ee = Dbar_e1
            pbar_ex = Dbar_e1 * (1 - decay_factor) + Dbar_e2 + Dbar_e3 * decay_factor
        else:
            p_ee = D_e2 * decay_factor + D_e3 * (1 - decay_factor)
            p_ex = D_e1 + D_e3
            pbar_ee = Dbar_e3
            pbar_ex = Dbar_e1 + Dbar_e2 * decay_factor + Dbar_e3 * (1 - decay_factor)

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = p_ex
        p[Flavor.NU_X,Flavor.NU_E] = (1 - p_ee)/2
        p[Flavor.NU_X,Flavor.NU_X] = 1 - p_ex/2

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = pbar_ex
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = (1 - pbar_ee)/2
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = 1 - pbar_ex/2

        return p


class AdiabaticMSWes(FourFlavorTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large that the electron-sterile resonances
       are inside the neutrinosphere;
    2. the “outer” or H' electron-sterile MSW resonance is adiabatic;
    3. the “inner” or H'' electron-sterile MSW resonance (where the electron fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, Phys. Rev. D 90, 033013 (2014).
    """
    def __init__(self, mix_params):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance
        """
        super().__init__(mix_params)   

    def __str__(self):
        return f'AdiabaticMSWes'
    
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
        p : 4x4 array or array of 4 x 4 arrays 
        """    
        D_e1, D_e2, D_e3, D_e4, D_s1, D_s2, D_s3, D_s4, Dbar_e1, Dbar_e2, Dbar_e3, Dbar_e4, Dbar_s1, Dbar_s2, Dbar_s3, Dbar_s4 = self.vacuum(E) 

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            p_ee = D_e4
            p_ex = D_e1 + D_e2
            p_xe = (1 - D_e4 - D_s4) / 2
            p_xx = (2 - D_e1 - D_e2 - D_s1 - D_s2) / 2

            pbar_ee = Dbar_e1
            pbar_ex = Dbar_e3 + Dbar_e4
            pbar_xe = (1 - Dbar_e1 - Dbar_s1) / 2
            pbar_xx = (2 - Dbar_e3 - Dbar_e4 - Dbar_s3 - Dbar_s4) / 2
        else:
            p_ee = D_e4
            p_ex = D_e1 + D_e3
            p_xe = (1 - D_e4 - D_s4) / 2
            p_xx = (2 - D_e1 - D_e3 - D_s1 - D_s3) / 2

            pbar_ee = Dbar_e4
            pbar_ex = Dbar_e1 + Dbar_e3
            pbar_xe = (1 - Dbar_e4 - Dbar_s4) / 2
            pbar_xx = (2 - Dbar_e1 - Dbar_e3 - Dbar_s1 - Dbar_s3) / 2

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = p_ex
        p[Flavor.NU_X,Flavor.NU_E] = p_xe
        p[Flavor.NU_X,Flavor.NU_X] = p_xx

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = pbar_ex
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = pbar_xe
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = pbar_xx

        return p 
        

class NonAdiabaticMSWes(FourFlavorTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large that the electron-sterile resonances
       are inside the neutrinosphere;
    2. the “outer” or H' electron-sterile MSW resonance is non-adiabatic;
    3. the “inner” or H'' electron-sterile MSW resonance (where the electron fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, Phys. Rev. D 90, 033013 (2014).
    """
    def __init__(self, mix_params, mh=MassHierarchy.NORMAL):
        """Initialize flavor transformation
        
        Parameters
        ----------
        mix_params : MixingParameters3Flavor instance
        """
        super().__init__(mix_params)   

    def __str__(self):
        return f'NonAdiabaticMSWes'
    
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
        p : 4x4 array or array of 4 x 4 arrays 
        """                
        D_e1, D_e2, D_e3, D_e4, D_s1, D_s2, D_s3, D_s4, Dbar_e1, Dbar_e2, Dbar_e3, Dbar_e4, Dbar_s1, Dbar_s2, Dbar_s3, Dbar_s4 = self.vacuum(E) 

        if self.mix_params.mass_order == MassHierarchy.NORMAL:
            p_ee = D_e3
            p_ex = D_e1 + D_e2
            p_xe = (1 - D_e3) / 2
            p_xx = (2 - D_e1 - D_e2 - D_s1 - D_s2) / 2

            pbar_ee = Dbar_e1
            pbar_ex = Dbar_e2 + Dbar_e3
            pbar_xe = (1 - Dbar_e1) / 2
            pbar_xx = (2 - Dbar_e2 - Dbar_e3 - Dbar_s2 - Dbar_s3) / 2
        else:
            p_ee = D_e2
            p_ex = D_e1 + D_e3
            p_xe = (1 - D_e2) / 2
            p_xx = (2 - D_e1 - D_e3 - D_s1 - D_s3) / 2

            pbar_ee = Dbar_e3
            pbar_ex = Dbar_e1 + Dbar_e2
            pbar_xe = (1 - Dbar_e3) / 2
            pbar_xx = (2 - Dbar_e1 - Dbar_e2 - Dbar_s1 - Dbar_s2) / 2

        p = np.zeros((4,4,len(E)))

        p[Flavor.NU_E,Flavor.NU_E] = p_ee 
        p[Flavor.NU_E,Flavor.NU_X] = p_ex
        p[Flavor.NU_X,Flavor.NU_E] = p_xe
        p[Flavor.NU_X,Flavor.NU_X] = p_xx

        p[Flavor.NU_E_BAR,Flavor.NU_E_BAR] = pbar_ee
        p[Flavor.NU_E_BAR,Flavor.NU_X_BAR] = pbar_ex
        p[Flavor.NU_X_BAR,Flavor.NU_E_BAR] = pbar_xe
        p[Flavor.NU_X_BAR,Flavor.NU_X_BAR] = pbar_xx

        return p   
        
   
