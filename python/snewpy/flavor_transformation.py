# -*- coding: utf-8 -*-
"""Simple supernova oscillation physics.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c

# Used to convert model parameters in eV and kpc to ergs and cm.
c = 29979245800 # cm/s
eV = 1.60218e-12 # ergs
kpc = 1000. * 3.086e+18 # cm


class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    @abstractmethod
    def prob_ee(self):
        """Electron neutrino survival probability."""
        pass

    @abstractmethod
    def prob_ex(self):
        """Electron -> flavor X neutrino transition probability."""
        pass

    @abstractmethod
    def prob_xx(self):
        """Flavor X neutrino survival probability."""
        pass

    @abstractmethod
    def prob_xe(self):
        """Flavor X -> electron neutrino transition probability."""
        pass    

    @abstractmethod
    def prob_eebar(self):
        """Electron antineutrino survival probability."""
        pass    

    @abstractmethod
    def prob_exbar(self):
        """Electron -> flavor X antineutrino transition probability."""
        pass   

    @abstractmethod
    def prob_xxbar(self):
        """Flavor X -> flavor X antineutrino survival probability."""
        pass

    @abstractmethod
    def prob_xebar(self):
        """Flavor X -> electron antineutrino transition probability."""
        pass    


class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def __init__(self,parameters):
        pass

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return 1.    

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return 1.

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.


class AdiabaticMSW_NMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming normal mass ordering."""

    def __init__(self, theta12, theta13, theta23):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return self.De3

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return self.De1  

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.


class AdiabaticMSW_IMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming inverted mass ordering."""

    def __init__(self, theta12, theta13, theta23):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return self.De2

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return self.De3  

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.
 

class NonAdiabaticMSWH_NMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming normal mass ordering."""

    def __init__(self, theta12, theta13, theta23):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return self.De2

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return self.De1  

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.


class NonAdiabaticMSWH_IMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming inverted mass ordering."""

    def __init__(self, theta12, theta13, theta23):
        """Initialize transformation matrix.

        Parameters
        ----------
        theta12 : astropy.units.quantity.Quantity
            Mixing angle 1->2 in PMNS matrix.
        theta13 : astropy.units.quantity.Quantity
            Mixing angle 1->3 in PMNS matrix.
        theta23 : astropy.units.quantity.Quantity
            Mixing angle 2->3 in PMNS matrix.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return self.De2

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return self.De1  

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.
       

class TwoFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def __init__(self):
        pass

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return 0.5    

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return 0.5

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.


class ThreeFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def __init__(self):
        pass

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        return 1./3.    

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return 1. - self.prob_ee(t,E)

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return (1. + self.prob_ee(t,E)) / 2.   

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2.  

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return 1./3.

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        return 1. - self.prob_eebar(t,E)       

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return (1. + self.prob_eebar(t,E)) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.


class NeutrinoDecay_IMO(FlavorTransformation):
    """Decay effect, assuming inverted mass ordering."""

    def __init__(self, theta12, theta13, theta23, mass, tau, dist):
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
            Neutrino mass.
        tau : astropy.units.quantity.Quantity
            Neutrino decay constant.
        dist : astropy.units.quantity.Quantity
            Source distance.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2
        self.m = mass * eV / (c**2)
        self.tau = tau
        self.d = dist * kpc

    def getgamma(self, E):
        return m*c/(E*tau)

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        pe_array = []
        for energy in E:
            pe_array.append( self.De2*np.exp(-self.getgamma(energ)*d) + self.De3*(1-np.exp(-self.getgamma(energy)*d)) )
        pe_array = np.array(pe_array)
        return pe_array

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return self.De1 + self.De2

    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return 1. - self.prob_ex(t,E) / 2.

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) /2.

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return self.De3  

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        pxbar_array = []
        for energy in E:
            pxbar_array.append( self.De1 + self.De2*np.exp(-self.getgamma(energy)*d) + self.De3*(1-np.exp(-self.getgamma(energy)*d)) )
        pxbar_array = np.array(pxbar_array)
        return pxbar_array

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return 1. - self.prob_exbar(t,E) / 2.

    def prob_xebar(self, t, E):  
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.       
        

class NeutrinoDecay_NMO(FlavorTransformation):
    """Decay effect, assuming normal mass ordering.""" 

    def __init__(self, theta12, theta13, theta23, mass, tau, dist):
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
            Neutrino mass.
        tau : astropy.units.quantity.Quantity
            Neutrino decay constant.
        dist : astropy.units.quantity.Quantity
            Source distance.
        """
        self.De1 = (np.cos(theta12) * np.cos(theta13))**2
        self.De2 = (np.sin(theta12) * np.cos(theta13))**2
        self.De3 = np.sin(theta13)**2
        self.m = mass * eV / (c**2)
        self.tau = tau
        self.d = dist * kpc

    def getgamma(self, E):
        return m*c/(E*tau)

    def prob_ee(self, t, E):
        """Electron neutrino survival probability."""
        pe_array = []
        for energy in E:
            pe_array.append( self.De1*(1-np.exp(-self.getgamma(energy)*d)) + self.De3*np.exp(-self.getgamma(energy)*d) ) 
        pe_array = np.array(pe_array)
        return pe_array

    def prob_ex(self, t, E):
        """Electron -> flavor X neutrino transition probability."""
        return self.De1 + self.De3
 
    def prob_xx(self, t, E):
        """Flavor X neutrino survival probability."""
        return 1. - self.prob_ex(t,E) / 2.

    def prob_xe(self, t, E):
        """Flavor X -> electron neutrino transition probability."""
        return (1. - self.prob_ee(t,E)) / 2. 

    def prob_eebar(self, t, E):
        """Electron antineutrino survival probability."""
        return self.De3

    def prob_exbar(self, t, E):
        """Electron -> flavor X antineutrino transition probability."""
        pxbar_array = []
        for energy in E:
            pxbar_array.append( self.De1*(1-np.exp(-self.getgamma(energy)*d)) + self.De2 + self.De3*np.exp(-self.getgamma(energy)*d) )
        pxbar_array = np.array(pxbar_array)
        return pxbar_array

    def prob_xxbar(self, t, E):
        """Flavor X -> flavor X antineutrino survival probability."""
        return 1. - self.prob_exbar(t,E) / 2.

    def prob_xebar(self, t, E):
        """Flavor X -> electron antineutrino transition probability."""
        return (1. - self.prob_eebar(t,E)) / 2.

