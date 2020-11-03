# -*- coding: utf-8 -*-
"""Simple supernova oscillation physics.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC
from enum import Enum

import numpy as np
from astropy import units as u
from astropy import constants as c


class MassHierarchy(Enum):
    """Neutrino mass ordering: normal or inverted.
    """
    NORMAL = 1
    INVERTED = 2


class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    @abstractmethod
    def prob_ee(self):
        """Electron neutrino survival probability.
        """
        pass

    @abstractmethod
    def prob_ex(self):
        """Electron -> flavor X neutrino transition probability.
        """
        pass

    @abstractmethod
    def prob_xx(self):
        """Flavor X neutrino survival probability.
        """
        pass

    @abstractmethod
    def prob_xe(self):
        """Flavor X -> electron neutrino transition probability.
        """
        pass    

    @abstractmethod
    def prob_eebar(self):
        """Electron antineutrino survival probability.
        """
        pass    

    @abstractmethod
    def prob_exbar(self):
        """Electron -> flavor X antineutrino transition probability.
        """
        pass   

    @abstractmethod
    def prob_xxbar(self):
        """Flavor X -> flavor X antineutrino survival probability.
        """
        pass

    @abstractmethod
    def prob_xebar(self):
        """Flavor X -> electron antineutrino transition probability.
        """
        pass    


class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def __init__(self):
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
        prob : float or ndarray
            Transition probability.
        """
        return 1.    

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability.

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
        return 1. - self.prob_ee(t,E)

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
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability.

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
        return 1.

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability.

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
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability.

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
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability.

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


class AdiabaticMSW(FlavorTransformation):
    """Adiabatic MSW effect."""

    def __init__(self, theta12, theta13, theta23, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        mh : MassHierarchy
            MassHierarchy.NORMAL or MassHierarchy.INVERTED.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2
        
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

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
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De3
        else:
            return self.De2

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability.

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
        return 1. - self.prob_ee(t,E)

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
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability.

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
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1  
        else:
            return self.De3

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability.

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
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability.

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
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability.

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


class NonAdiabaticMSW(FlavorTransformation):
    """Nonadiabatic MSW effect."""

    def __init__(self, theta12, theta13, theta23, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        mh : MassHierarchy
            MassHierarchy.NORMAL or MassHierarchy.INVERTED.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2
        
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

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
        return self.De2

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability.

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
        return 1. - self.prob_ee(t,E)

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
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability.

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
        return self.De1  

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability.

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
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability.

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
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability.

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


class TwoFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def __init__(self):
        pass

    def prob_ee(self, t, E):
        """Electron neutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.
        """
        return 0.5    

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability.

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
        return 1. - self.prob_ee(t,E)

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
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability.

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
        return 0.5

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability.

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
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability.

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
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability.

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


class ThreeFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def __init__(self):
        pass

    def prob_ee(self, t, E):
        """Electron neutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.
        """
        return 1./3.    

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability.

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
        return 1. - self.prob_ee(t,E)

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
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability.

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
        return 1./3.

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability.

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
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability.

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
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability.

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


class NeutrinoDecay(FlavorTransformation):
    """Decay effect, where the heaviest neutrino of mass m3 decays to the
    lightest neutrino of mass m1. For a description and typical parameters,
    see A. de Gouvêa et al., PRD 101:043013, 2020, arXiv:1910.01127.
    """

    def __init__(self, theta12, theta13, theta23, mass, tau, dist, mh):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        mass : astropy.units.quantity.Quantity
            Mass m3 of the heaviest neutrino; expect in eV/c^2.
        tau : astropy.units.quantity.Quantity
            Lifetime of the heaviest neutrino.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        mh : MassHierarchy
            MassHierarchy.NORMAL or MassHierarchy.INVERTED.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2
        self.m = mass
        self.tau = tau
        self.d = dist
        
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

    def gamma(self, E):
        """Decay width of the m3 in the nu3 restframe.

        Parameters
        ----------
        E : float
            Energy of the nu3.

        Returns
        -------
        Gamma : float
            Decay width of the m3, in units of 1/length.
        """
        return m*c.c / (E*tau)

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
        pe_array = []

        # NMO case.
        if self.mass_order == MassHierarchy.NORMAL:
            for energy in E:
                pe_array.append(
                    self.De1*(1-np.exp(-self.gamma(energy)*self.d)) +
                    self.De3*np.exp(-self.gamma(energy)*self.d))
            pe_array = np.array(pe_array)
        # IMO case.
        else:
            for energy in E:
                pe_array.append(
                    self.De2*np.exp(-self.gamma(energ)*self.d) +
                    self.De3*(1-np.exp(-self.gamma(energy)*self.d)))
            pe_array = np.array(pe_array)

        return pe_array

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability.

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
        """Flavor X -> electron neutrino transition probability.

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
        """Electron -> flavor X antineutrino transition probability.

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
        pxbar_array = []

        # NMO case.
        if self.mass_order == MassHierarchy.Normal:
            for energy in E:
                pxbar_array.append(
                    self.De1*(1-np.exp(-self.gamma(energy)*self.d)) +
                    self.De2 + self.De3*np.exp(-self.gamma(energy)*self.d))
            pxbar_array = np.array(pxbar_array)
        # IMO case.
        else:
            for energy in E:
                pxbar_array.append(
                    self.De1 + self.De2*np.exp(-self.gamma(energy)*self.d) +
                    self.De3*(1-np.exp(-self.gamma(energy)*self.d)))
            pxbar_array = np.array(pxbar_array)

        return pxbar_array

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability.

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
        """Flavor X -> electron antineutrino transition probability.

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
        
