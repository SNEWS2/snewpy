# -*- coding: utf-8 -*-
"""Supernova oscillation physics for flavors e, X, e-bar, X-bar.

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


class MixingParameters:
    """Best-fit parameters of the PMNS matrix and mass differences, assuming
    three neutrino flavors. See www.nu-fit.org for current global fits.
    """
    def __init__(self, mh=MassHierarchy.NORMAL):
        """Initialize the neutrino mixing parameters.

        Parameters
        ----------
        mh : MassHierarchy
            Desired mass ordering: NORMAL or INVERTED.
        """
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

        # Values from JHEP 09 (2020) 178 [arXiv:2007.14792] and www.nu-fit.org.
        # The reported precision is not significant given current uncertainties
        # on these parameters, but is useful for comparing to the table of
        # parameters presented on nu-fit.org.
        if self.mass_order == MassHierarchy.NORMAL:
            self.theta12 = 33.44 * u.deg
            self.theta13 =  8.57 * u.deg
            self.theta23 = 49.20 * u.deg
            self.deltaCP = 197 * u.deg
            self.dm21_2  = 7.42e-5 * u.eV**2
            self.dm31_2  = 2.517e-3 * u.eV**2
        else:
            self.theta12 = 33.45 * u.deg
            self.theta13 =  8.60 * u.deg
            self.theta23 = 49.30 * u.deg
            self.deltaCP = 282 * u.deg
            self.dm21_2  = 7.42e-5 * u.eV**2
            self.dm32_2  = -2.498e-3 * u.eV**2

    def get_mixing_angles(self):
        """Mixing angles of the PMNS matrix.
        
        Returns
        -------
        angles : tuple
            Angles theta12, theta13, theta23.
        """
        return (self.theta12, self.theta13, self.theta23)


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

        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

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

        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

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
    def __init__(self, mix_angles=None, mass=1*u.eV/c.c**2, tau=1*u.day, dist=10*u.kpc, mh=MassHierarchy.NORMAL):
        """Initialize transformation matrix.

        Parameters
        ----------
        mix_angles : tuple or None
            If not None, override default mixing angles using tuple (theta12, theta13, theta23).
        mass : astropy.units.quantity.Quantity
            Mass m3 of the heaviest neutrino; expect in eV/c^2.
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

        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

        self.m = mass
        self.tau = tau
        self.d = dist

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
        
