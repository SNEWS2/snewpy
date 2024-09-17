# -*- coding: utf-8 -*-
"""Supernova oscillation physics 

For measured mixing angles and latest global analysis results, visit
http://www.nu-fit.org/.
"""

import logging
logger = logging.getLogger()

from abc import abstractmethod, ABC

from importlib.resources import files

from snewpy.flavor  import FlavorMatrix, ThreeFlavor
from snewpy.neutrino import MassHierarchy
from snewpy.neutrino import MixingParameters
from .in_sn import SNTransformation
from .in_vacuum import VacuumTransformation, NoVacuumTransformation
from .in_earth import EarthTransformation, NoEarthMatter
###############################################################################

class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""
    def __str__(self):
        return self.__class__.__name__
    
    @abstractmethod
    def P_ff(self, t, E)->FlavorMatrix:
        pass
    
    def apply(self, flux):
        M = self.P_ff(flux.time, flux.energy)
        M = (flux.flavor_scheme<<M<<flux.flavor_scheme)
        return M@flux
     
###############################################################################
class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def P_ff(self, t, E):
        p = FlavorMatrix.eye(ThreeFlavor)
        return p

    def apply(self, flux):
        """This transformation returns the object without transform"""
        return flux
###############################################################################
class CompleteExchange(FlavorTransformation):
    """Survival probabilities for the case when the electron flavors 
       are half exchanged with the mu flavors and the half with the tau flavors.
    """
    def P_ff(self, t, E):
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1,f2):
            return (f1.is_neutrino==f2.is_neutrino)*(f1!=f2)*0.5

        return P
###############################################################################
class ThreeFlavorDecoherence(FlavorTransformation):
    """Equal mixing of all threen eutrino matter states and antineutrino matter states"""

    def __init__(self):
        """Initialize ThreeFlavorTransformation to default case"""
        super().__init__() 

    def P_ff(self, t, E): 
        """Equal mixing so Earth matter has no effect"""
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1,f2):
            return (f1.is_neutrino==f2.is_neutrino)*1/3.
        return P

#####################################################################
class TransformationChain(FlavorTransformation):
    def __init__(self, 
                 in_sn: SNTransformation,
                 in_vacuum: VacuumTransformation=NoVacuumTransformation(),
                 in_earth: EarthTransformation=NoEarthMatter(),
                 *,
                 mixing_params=MixingParameters()
                ):
        if in_sn is None:
            return NoTransformation()
        self.in_sn = in_sn
        self.in_vacuum = in_vacuum
        self.in_earth = in_earth
        self.transforms = [in_sn, in_vacuum, in_earth]
        self.set_mixing_params(mixing_params)

    def set_mixing_params(self, mixing_params):
        #set the mixing parameters to all the inner classes
        self.mixing_params = mixing_params
        for t in self.transforms:
            t.mixing_params = mixing_params
                    
    def P_ff(self, t, E)->FlavorMatrix:
        in_sn, in_vacuum, in_earth = self.transforms
        return in_earth.P_fm(t,E) @ in_vacuum.P_mm(t,E) @ in_sn.P_mf(t,E)
    
    def __str__(self):
        s = '+'.join([t.__class__.__name__ for t in self.transforms])+'_'+self.mixing_params.mass_order.name
        return s
