# -*- coding: utf-8 -*-
"""Simple supernova oscillation physics.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u


# Mixing angles, in degrees. Data source:
# - JHEP 09 (2020) 178 [arXiv:2007.14792]
# - NuFIT 5.0 (2020), www.nu-fit.org
theta12 = 33.44 * u.degree
theta13 =  8.57 * u.degree
theta23 = 49.0  * u.degree


class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival
    probability."""

    @abstractmethod
    def p(self):
        """Neutrino survival probability."""
        pass
    
    @abstractmethod
    def pbar(self):
        """Antineutrino survival probability."""
        pass    


class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def p(self):
        """Neutrino survival probability."""
        return 1.
    
    def pbar(self):
        """Antineutrino survival probability."""
        return 1.


class AdiabaticMSW_NMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming normal mass ordering."""

    def p(self):
        """Neutrino survival probability."""
        return np.sin(theta13)**2
    
    def pbar(self):
        """Antineutrino survival probability."""
        return (np.cos(theta12)*np.cos(theta13))**2    


class AdiabaticMSW_IMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming inverted mass ordering."""

    def p(self):
        """Neutrino survival probability."""
        return (np.sin(theta12)*np.cos(theta13))**2

    def pbar(self):
        """Antineutrino survival probability."""
        return np.sin(theta13)**2


class TwoFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def p(self):
        """Neutrino survival probability."""
        return 0.5
    
    def pbar(self):
        """Antineutrino survival probability."""
        return 0.5


class ThreeFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def p(self):
        """Neutrino survival probability."""
        return 1./3.
    
    def pbar(self):
        """Antineutrino survival probability."""
        return 1./3.
