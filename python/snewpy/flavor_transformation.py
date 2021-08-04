# -*- coding: utf-8 -*-
"""Supernova oscillation physics for flavors e, X, e-bar, X-bar.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c

from .neutrino import MassHierarchy, MixingParameters


class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    @abstractmethod
    def prob_ee(self):
        """Electron neutrino survival probability.
        """
        pass

    @abstractmethod
    def prob_ex(self):
        """X -> e neutrino transition probability.
        """
        pass

    @abstractmethod
    def prob_xx(self):
        """Flavor X neutrino survival probability.
        """
        pass

    @abstractmethod
    def prob_xe(self):
        """e -> X neutrino transition probability.
        """
        pass    

    @abstractmethod
    def prob_eebar(self):
        """Electron antineutrino survival probability.
        """
        pass    

    @abstractmethod
    def prob_exbar(self):
        """X -> e antineutrino transition probability.
        """
        pass   

    @abstractmethod
    def prob_xxbar(self):
        """X -> X antineutrino survival probability.
        """
        pass

    @abstractmethod
    def prob_xebar(self):
        """e -> X antineutrino transition probability.
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
        return 1.

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
        return 1. - self.prob_eebar(t,E)       

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
        return (1. + self.prob_eebar(t,E)) / 2.

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

    
class CompleteExchange(FlavorTransformation):
    """Survival probabilities for the case when the electron flavors are completely exchanged with the x flavor."""

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
        return 0.    

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
        return 0.

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
        return 1. - self.prob_eebar(t,E)       

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
        return (1. + self.prob_eebar(t,E)) / 2.

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


class AdiabaticMSW(FlavorTransformation):
    """Adiabatic MSW effect."""

    def __init__(self, mix_angles=None, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple or None
            If not None, override default mixing angles using tuple (theta12, theta13, theta23).
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
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1  
        else:
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
        return 1. - self.prob_eebar(t,E)       

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
        return (1. + self.prob_eebar(t,E)) / 2.

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


class NonAdiabaticMSWH(FlavorTransformation):
    """Nonadiabatic MSW effect."""

    def __init__(self, mix_angles=None, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple or None
            If not None, override default mixing angles using tuple (theta12, theta13, theta23).
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
        return self.De1  

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
        return 1. - self.prob_eebar(t,E)       

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
        return (1. + self.prob_eebar(t,E)) / 2.

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


class TwoFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def __init__(self, mix_angles=None, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple or None
            If not None, override default mixing angles using tuple (theta12, theta13, theta23).
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

    def prob_ee(self, t, E):
        """Electron neutrino survival probability.

        Parameters
        ----------
        t : float or ndarray
            List of times.
        E : float or ndarray
            List of energies.
        """
        if self.mass_order == MassHierarchy.NORMAL:
            return (self.De2+self.De3)/2.    
        else:
            return self.De2

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
        if self.mass_order == MassHierarchy.NORMAL:
            return self.De1
        else:
            return (self.De1+self.De3)/2

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
        return 1. - self.prob_eebar(t,E)       

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
        return (1. + self.prob_eebar(t,E)) / 2.

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
        return 1./3.

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
        return 1. - self.prob_eebar(t,E)       

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
        return (1. + self.prob_eebar(t,E)) / 2.

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
                    self.De2*np.exp(-self.gamma(energy)*self.d) +
                    self.De3*(1-np.exp(-self.gamma(energy)*self.d)))
            pe_array = np.array(pe_array)

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
        pxbar_array = []

        # NMO case.
        if self.mass_order == MassHierarchy.NORMAL:
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
    
    def __init__(self, mix_angles, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple
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
    
    def __init__(self, mix_angles, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple
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
    
