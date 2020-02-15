# -*- coding: utf-8 -*-
"""Simple supernova oscillation physics.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np


theta12 = np.deg2rad(33.) 
theta13 = np.deg2rad(9.) 
theta23 = np.deg2rad(45.)


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
        return 1.
    
    def pbar(self):
        return 1.


class AdiabaticMSW_NMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming normal mass ordering."""

    def p(self):
        return pow(np.sin(theta13),2.)
    
    def pbar(self):
        return pow(np.cos(theta12)*np.cos(theta13),2.)    


class AdiabaticMSW_IMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming inverted mass ordering."""

    def p(self):
        return pow(np.sin(theta12)*np.cos(theta13),2.)

    def pbar(self):
        return pow(np.sin(theta13),2.)


class TwoFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def p(self):
        return 0.5
    
    def pbar(self):
        return 0.5


class ThreeFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def p(self):
        return 1./3.
    
    def pbar(self):
        return 1./3.
