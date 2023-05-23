# -*- coding: utf-8 -*-
"""Supernova oscillation physics for flavors e, X, e-bar, X-bar.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c

from .neutrino import MassHierarchy, MixingParameters, Flavor
from functools import lru_cache

class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""
        
    @abstractmethod
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
        pass

    @abstractmethod
    def prob_ex(self, t, E):
        """X -> e neutrino transition probability.

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
        pass

    @abstractmethod
    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability.

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
        pass

    @abstractmethod
    def prob_xe(self, t, E):
        """e -> X neutrino transition probability.

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
        return self.matrix(t,E)[Flavor.NU_X, Flavor.NU_E]

    @abstractmethod
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
        pass

    @abstractmethod
    def prob_exbar(self, t, E):
        """X -> e antineutrino transition probability.

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
        pass

    @abstractmethod
    def prob_xxbar(self, t, E):
        """X -> X antineutrino survival probability.

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
        pass

    @abstractmethod
    def prob_xebar(self, t, E):
        """e -> X antineutrino survival probability.

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
        pass

   

class MatrixFlavorTransformation(FlavorTransformation):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    @abstractmethod
    def matrix(self, t, E):
        pass
        
    def prob_ee(self, t, E):
        return self.matrix(t,E)[Flavor.NU_E, Flavor.NU_E]
        
    def prob_ex(self, t, E):
        return self.matrix(t,E)[Flavor.NU_E, Flavor.NU_X]
        
    def prob_xx(self, t, E):
        return self.matrix(t,E)[Flavor.NU_X, Flavor.NU_X]
        
    def prob_xe(self, t, E):
        return self.matrix(t,E)[Flavor.NU_X, Flavor.NU_E]

    def prob_eebar(self, t, E):
        return self.matrix(t,E)[Flavor.NU_E_BAR, Flavor.NU_E_BAR]

    def prob_exbar(self, t, E):
        return self.matrix(t,E)[Flavor.NU_E_BAR, Flavor.NU_X_BAR]

    def prob_xxbar(self, t, E):
        return self.matrix(t,E)[Flavor.NU_X_BAR, Flavor.NU_X_BAR]

    def prob_xebar(self, t, E):
        return self.matrix(t,E)[Flavor.NU_X_BAR, Flavor.NU_E_BAR]



class ConstMatrixTransformation(MatrixFlavorTransformation):
    """Transition matrix that doesn't depend on energy and time"""
    def __init__(self, matrix:np.ndarray):
        assert matrix.shape==(len(Flavor),len(Flavor)), f'Wrong matrix shape {matrix.shape}'
        self._matrix_base = matrix

    def matrix(self, t, E):
        return self.__broadcast_matrix(self._matrix_base, t,E)

    @staticmethod
    def __broadcast_matrix(M, *arrays):
        #convert additional axes to arrays
        #arrays = [np.array(a, copy=False) for a in axes]
        extended_shape  = [*M.shape, *(1 for a in arrays if not np.isscalar(a))]
        broadcast_shape = [*M.shape, *(len(a) for a in arrays if not np.isscalar(a))]
        M1 = M.reshape(extended_shape)
        return np.broadcast_to(M1, broadcast_shape)
        
    @staticmethod
    def from_Pee_Peebar(prob_ee:float, prob_eebar:float):
        def _m2(p):
            return np.array([ [p, 1-p],
                              [1-p, p]])
        z = np.zeros((2,2))
        M = np.block([[_m2(prob_ee), z],
                        [z, _m2(prob_eebar)]])
        return ConstMatrixTransformation(M)


def NoTransformation():
    return ConstMatrixTransformation.from_Pee_Peebar(prob_ee = 1, prob_eebar = 1)
def CompleteExchange():
    return ConstMatrixTransformation.from_Pee_Peebar(prob_ee = 0, prob_eebar = 0)

def AdiabaticMSW(pars:MixingParameters):
    """Adiabatic MSW effect."""
    De1 = float((np.cos(pars.theta12) * np.cos(pars.theta13))**2)
    De2 = float((np.sin(pars.theta12) * np.cos(pars.theta13))**2)
    De3 = float(np.sin(pars.theta13)**2)
    if pars.mass_order == MassHierarchy.NORMAL:
        Pee, Pee_bar = De3, De1
    else: 
        Pee, Pee_bar = De2, De3
    return ConstMatrixTransformation.from_Pee_Peebar(Pee, Pee_bar)
    
def NonAdiabaticMSWH(pars:MixingParameters):
    """Nonadiabatic MSW effect."""
    De1 = float((np.cos(pars.theta12) * np.cos(pars.theta13))**2)
    De2 = float((np.sin(pars.theta12) * np.cos(pars.theta13))**2)
    De3 = float(np.sin(pars.theta13)**2)
    Pee, Pee_bar = (De2, De1)
    return ConstMatrixTransformation.from_Pee_Peebar(Pee, Pee_bar)

def TwoFlavorDecoherence(pars:MixingParameters):
    """Star-earth transit survival probability: two flavor case."""
    De1 = float((np.cos(pars.theta12) * np.cos(pars.theta13))**2)
    De2 = float((np.sin(pars.theta12) * np.cos(pars.theta13))**2)
    De3 = float(np.sin(pars.theta13)**2)
    if pars.mass_order == MassHierarchy.NORMAL:
        Pee, Pee_bar = (De2+De3)/2, De1
    else:
        Pee, Pee_bar = De1, (De1+De3)/2
    return ConstMatrixTransformation.from_Pee_Peebar(Pee, Pee_bar)

def ThreeFlavorDecoherence(pars:MixingParameters):
    """Star-earth transit survival probability: three flavor case."""
    return ConstMatrixTransformation.from_Pee_Peebar(1/3, 1/3)

class NeutrinoDecay(FlavorTransformation):
    """Decay effect, where the heaviest neutrino decays to the lightest
    neutrino. For a description and typical parameters, see A. de Gouvêa et al.,
    PRD 101:043013, 2020, arXiv:1910.01127.
    """
    def __init__(self, mix_angles=None, mass=1*u.eV/c.c**2, tau=1*u.day, dist=10*u.kpc, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple or None
            If not None, override default mixing angles using tuple (theta12, theta13, theta23).
        mass : astropy.units.quantity.Quantity
            Mass of the heaviest neutrino; expect in eV/c^2.
        tau : astropy.units.quantity.Quantity
            Lifetime of the heaviest neutrino.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        mh : MassHierarchy
            MassHierarchy.NORMAL or MassHierarchy.INVERTED.
        """
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

        if mix_angles is not None:
            theta12, theta13, theta23 = mix_angles
        else:
            pars = MixingParameters(mh)
            theta12, theta13, theta23 = pars.get_mixing_angles()

        self.De1 = float((np.cos(theta12) * np.cos(theta13))**2)
        self.De2 = float((np.sin(theta12) * np.cos(theta13))**2)
        self.De3 = float(np.sin(theta13)**2)

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
        prob : float or ndarray
            Transition probability.
        """
        # NMO case.
        if self.mass_order == MassHierarchy.NORMAL:
            pe_array = self.De1*(1-np.exp(-self.gamma(E)*self.d)) + \
                       self.De3*np.exp(-self.gamma(E)*self.d)
        # IMO case.
        else:
            pe_array = self.De2*np.exp(-self.gamma(E)*self.d) + \
                       self.De3*(1-np.exp(-self.gamma(E)*self.d))
        return pe_array

    def prob_ex(self, t, E):
        """X -> e neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        # NMO case.
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1 + self.De3
        # IMO case.
        else:
            return self.De1 + self.De2

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        return 1. - self.prob_ex(t,E) / 2.

    def prob_xe(self, t, E):
        """e -> X neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        return (1. - self.prob_ee(t,E)) / 2.

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
        prob : float or ndarray
            Transition probability.
        """
        return self.De3  

    def prob_exbar(self, t, E):
        """X -> e antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        # NMO case.
        if self.mass_order == MassHierarchy.NORMAL:
            pxbar_array = self.De1*(1-np.exp(-self.gamma(E)*self.d)) + \
                          self.De2 + self.De3*np.exp(-self.gamma(E)*self.d)
        # IMO case.
        else:
            pxbar_array = self.De1 + self.De2*np.exp(-self.gamma(E)*self.d) + \
                          self.De3*(1-np.exp(-self.gamma(E)*self.d))
        return pxbar_array

    def prob_xxbar(self, t, E):
        """X -> X antineutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        return 1. - self.prob_exbar(t,E) / 2.

    def prob_xebar(self, t, E):  
        """e -> X antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        return (1. - self.prob_eebar(t,E)) / 2.       


class AdiabaticMSWes(FlavorTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large that the electron-sterile resonances
       are inside the neutrinosphere;
    2. the “outer” or H' electron-sterile MSW resonance is adiabatic;
    3. the “inner” or H'' electron-sterile MSW resonance (where the electron fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, Phys. Rev. D 90, 033013 (2014).
    """
    def __init__(self, mix_angles, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple
            Values for mixing angles (theta12, theta13, theta23, theta14).
        mh : MassHierarchy
            MassHierarchy.NORMAL or MassHierarchy.INVERTED.
        """
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

        theta12, theta13, theta23, theta14 = mix_angles

        self.De1 = float((np.cos(theta12) * np.cos(theta13) * np.cos(theta14))**2)
        self.De2 = float((np.sin(theta12) * np.cos(theta13) * np.cos(theta14))**2)
        self.De3 = float((np.sin(theta13) * np.cos(theta14))**2)
        self.De4 = float((np.sin(theta14))**2)
        self.Ds1 = float((np.cos(theta12) * np.cos(theta13) * np.sin(theta14))**2)
        self.Ds2 = float((np.sin(theta12) * np.cos(theta13) * np.sin(theta14))**2)
        self.Ds3 = float((np.sin(theta13) * np.sin(theta14))**2)
        self.Ds4 = float((np.cos(theta14))**2)
    
    def prob_ee(self, t, E):
        """e -> e neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """
        return self.De4       
        
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
        prob : float or ndarray
            Transition probability.
        """       
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1 + self.De2
        else:
            return self.De1 + self.De3     
        
    def prob_xx(self, t, E):
        """x -> x neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return ( 2 - self.De1 - self.De2 - self.Ds1 - self.Ds2 ) / 2 
        else:
            return ( 2 - self.De1 - self.De3 - self.Ds1 - self.Ds3 ) / 2          
        
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
        prob : float or ndarray
            Transition probability.
        """        
        return ( 1 - self.De4 - self.Ds4 )/2         
    
    def prob_eebar(self, t, E):
        """e -> e antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1
        else:
            return self.De3
        
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
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De3 + self.De4
        else:
            return self.De2 + self.De4   
        
    def prob_xxbar(self, t, E):
        """x -> x antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return ( 2 - self.De3 - self.De4 - self.Ds3 - self.Ds4 ) / 2 
        else:
            return ( 2 - self.De2 - self.De4 - self.Ds2 - self.Ds4 ) / 2     
        
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
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return ( 1 - self.De1 - self.Ds1 ) / 2
        else:
            return ( 1 - self.De3 - self.Ds3 ) / 2     


class NonAdiabaticMSWes(FlavorTransformation):
    """A four-neutrino mixing prescription. The assumptions used are that:

    1. the fourth neutrino mass is the heaviest but not so large that the electron-sterile resonances
       are inside the neutrinosphere;
    2. the “outer” or H' electron-sterile MSW resonance is non-adiabatic;
    3. the “inner” or H'' electron-sterile MSW resonance (where the electron fraction = 1/3) is non-adiabatic.

    For further insight see, for example, Esmaili, Peres, and Serpico, Phys. Rev. D 90, 033013 (2014).
    """
    def __init__(self, mix_angles, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.
        
        Parameters
        ----------
        mix_angles : tuple
            Values for mixing angles (theta12, theta13, theta23, theta14).
        mh : MassHierarchy
            MassHierarchy.NORMAL or MassHierarchy.INVERTED.
        """
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

        theta12, theta13, theta23, theta14 = mix_angles

        self.De1 = float((np.cos(theta12) * np.cos(theta13) * np.cos(theta14))**2)
        self.De2 = float((np.sin(theta12) * np.cos(theta13) * np.cos(theta14))**2)
        self.De3 = float((np.sin(theta13) * np.cos(theta14))**2)
        self.De4 = float((np.sin(theta14))**2)
        self.Ds1 = float((np.cos(theta12) * np.cos(theta13) * np.sin(theta14))**2)
        self.Ds2 = float((np.sin(theta12) * np.cos(theta13) * np.sin(theta14))**2)
        self.Ds3 = float((np.sin(theta13) * np.sin(theta14))**2)
        self.Ds4 = float((np.cos(theta14))**2)
    
    def prob_ee(self, t, E):
        """e -> e neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """                
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De3
        else:
            return self.De2        
        
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
        prob : float or ndarray
            Transition probability.
        """                
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1 + self.De2
        else:
            return self.De1 + self.De3        
        
    def prob_xx(self, t, E):
        """x -> x neutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """               
        if self.mass_order == MassHierarchy.NORMAL:
            return ( 2 - self.De1 - self.De2 - self.Ds1 - self.Ds2 ) / 2 
        else:
            return ( 2 - self.De1 - self.De3 - self.Ds1 - self.Ds3 ) / 2         
        
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
        prob : float or ndarray
            Transition probability.
        """                
        if self.mass_order == MassHierarchy.NORMAL:
            return (1 - self.De3 - self.Ds3)/2
        else:
            return (1 - self.De2 - self.Ds2) / 2        
    
    def prob_eebar(self, t, E):
        """e -> e antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """               
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1
        else:
            return self.De3    
        
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
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De2 + self.De3
        else:
            return self.De1 + self.De2        
        
    def prob_xxbar(self, t, E):
        """x -> x antineutrino transition probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.

        Returns
        -------
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return ( 2 - self.De2 - self.De3 - self.Ds2 - self.Ds3 ) / 2 
        else:
            return ( 2 - self.De1 - self.De2 - self.Ds1 - self.Ds2 ) / 2        
        
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
        prob : float or ndarray
            Transition probability.
        """        
        if self.mass_order == MassHierarchy.NORMAL:
            return ( 1 - self.De1 - self.Ds1 ) / 2
        else:
            return ( 1 - self.De3 - self.Ds3 ) / 2        
    
