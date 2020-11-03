# -*- coding: utf-8 -*-
"""Simple supernova oscillation physics.

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

from abc import abstractmethod, ABC

import numpy as np
from astropy import units as u
from astropy import constants as c

c = 29979245800 # cm/s
eV = 1.60218e-12 # ergs
kpc = 1000. * 3.086e+18 # cm

class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    @abstractmethod
    def pee(self):
        """Electron Neutrino survival probability."""
        pass

    @abstractmethod
    def pex(self):
        pass
    
    @abstractmethod
    def pxx(self):
        pass
    
    @abstractmethod
    def pxe(self):
        pass    
    
    @abstractmethod
    def peebar(self):
        """Electron Antineutrino survival probability."""
        pass    
    
    @abstractmethod
    def pexbar(self):
        pass   
    
    @abstractmethod
    def pxxbar(self):
        pass
    
    @abstractmethod
    def pxebar(self):
        pass    

class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def __init__(self,parameters):
        pass

    def pee(self, t, E):
        return 1.    
    
    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.
    
    def peebar(self, t, E):
        return 1.
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.

class AdiabaticMSW_NMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming normal mass ordering."""

    def __init__(self,parameters):
        self.De1 = (np.cos(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De2 = (np.sin(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De3 = np.sin(parameters[1]*np.pi/180.)**2
    
    def pee(self, t, E):
        return self.De3

    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.
    
    def peebar(self, t, E):
        return self.De1  
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.

class AdiabaticMSW_IMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming inverted mass ordering."""
    
    def __init__(self,parameters):
        self.De1 = (np.cos(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De2 = (np.sin(parameters[0])**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De3 = np.sin(parameters[1])**2
    
    def pee(self, t, E):
        return self.De2

    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.
    
    def peebar(self, t, E):
        return self.De3  
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.
    
class NonAdiabaticMSWH_NMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming normal mass ordering."""

    def __init__(self,parameters):
        self.De1 = (np.cos(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De2 = (np.sin(parameters[0])**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De3 = np.sin(parameters[1])**2
    
    def pee(self, t, E):
        return self.De2

    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.
    
    def peebar(self, t, E):
        return self.De1  
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.
    
class NonAdiabaticMSWH_IMO(FlavorTransformation):
    """Adiabatic MSW effect, assuming inverted mass ordering."""
    
    def __init__(self,parameters):
        self.De1 = (np.cos(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De2 = (np.sin(parameters[0])**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De3 = np.sin(parameters[1])**2
    
    def pee(self, t, E):
        return self.De2

    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.
    
    def peebar(self, t, E):
        return self.De1  
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.
       
class TwoFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: two flavor case."""

    def __init__(self,parameters):
        pass

    def pee(self, t, E):
        return 0.5    
    
    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.
    
    def peebar(self, t, E):
        return 0.5
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.

class ThreeFlavorDecoherence(FlavorTransformation):
    """Star-earth transit survival probability: three flavor case."""

    def __init__(self,parameters):
        pass

    def pee(self, t, E):
        return 1./3.    
    
    def pex(self, t, E):
        return 1. - self.pee(t,E)

    def pxx(self, t, E):
        return (1. + self.pee(t,E)) / 2.   
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2.  
    
    def peebar(self, t, E):
        return 1./3.
    
    def pexbar(self, t, E):
        return 1. - self.peebar(t,E)       
    
    def pxxbar(self, t, E):
        return (1. + self.peebar(t,E)) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.
    
class NeutrinoDecay_IMO(FlavorTransformation):
    """Decay effect, assuming inverted mass ordering."""
    
    def __init__(self,parameters):
        self.De1 = (np.cos(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De2 = (np.sin(parameters[0])**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De3 = np.sin(parameters[1])**2
        
        self.m=parameters[3] * eV / (c**2)
        self.tau=parameters[4]
        self.d=parameters[5] * kpc
    
    def getgamma(self, E):
        return m*c/(E*tau)

    def pee(self, t, E):
        pe_array = []
        for energy in E:
            pe_array.append( self.De2*np.exp(-self.getgamma(energ)*d) + self.De3*(1-np.exp(-self.getgamma(energy)*d)) )
        pe_array = np.array(pe_array)
        return pe_array
    
    def pex(self, t, E):
        return self.De1 + self.De2
    
    def pxx(self, t, E):
        return 1. - self.pex(t,E) / 2.
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) /2.
    
    def peebar(self, t, E):
        return self.De3  
   
    def pexbar(self, t, E):
        pxbar_array = []
        for energy in E:
            pxbar_array.append( self.De1 + self.De2*np.exp(-self.getgamma(energy)*d) + self.De3*(1-np.exp(-self.getgamma(energy)*d)) )
        pxbar_array = np.array(pxbar_array)
        return pxbar_array
    
    def pxxbar(self, t, E):
        return 1. - self.pexbar(t,E) / 2.
   
    def pxebar(self, t, E):  
        return (1. - self.peebar(t,E)) / 2.       
        
class NeutrinoDecay_NMO(FlavorTransformation):
    """decay effect, assuming normal mass ordering.""" 
    
    def __init__(self,parameters):
        self.De1 = (np.cos(parameters[0]*np.pi/180.)**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De2 = (np.sin(parameters[0])**2)*(np.cos(parameters[1]*np.pi/180.)**2)
        self.De3 = np.sin(parameters[1])**2
        
        self.m=parameters[3] * eV / (c**2)
        self.tau=parameters[4]
        self.d=parameters[5] * kpc
    
    def getgamma(self, E):
        return m*c/(E*tau)
    
    def pee(self, t, E):
        pe_array = []
        for energy in E:
            pe_array.append( self.De1*(1-np.exp(-self.getgamma(energy)*d)) + self.De3*np.exp(-self.getgamma(energy)*d) ) 
        pe_array = np.array(pe_array)
        return pe_array
    
    def pex(self, t, E):
        return self.De1 + self.De3
    
    def pxx(self, t, E):
        return 1. - self.pex(t,E) / 2.
    
    def pxe(self, t, E):
        return (1. - self.pee(t,E)) / 2. 
    
    def peebar(self, t, E):
        return self.De3
  
    def pexbar(self, t, E):
        pxbar_array = []
        for energy in E:
            pxbar_array.append( self.De1*(1-np.exp(-self.getgamma(energy)*d)) + self.De2 + self.De3*np.exp(-self.getgamma(energy)*d) )
        pxbar_array = np.array(pxbar_array)
        return pxbar_array
    
    def pxxbar(self, t, E):
        return 1. - self.pexbar(t,E) / 2.
    
    def pxebar(self, t, E):
        return (1. - self.peebar(t,E)) / 2.
    
