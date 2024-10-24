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
from snewpy.neutrino import MixingParameters, ThreeFlavorMixingParameters, FourFlavorMixingParameters
from .in_sn import SNTransformation
from .in_vacuum import VacuumTransformation, NoVacuumTransformation
from .in_earth import EarthTransformation, NoEarthMatter

from snewpy.flux import Container
###############################################################################

class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""
    def __str__(self):
        return self.__class__.__name__
    
    @abstractmethod
    def P_ff(self, t, E)->FlavorMatrix:
        r"""Transition probability matrix in flavor basis :math:`P_{\alpha\to\beta}`
        
        Parameters
        ----------
        """
        pass
    
    def apply(self, flux:Container)->Container:
        r"""Apply this transformation to the given flux, return transformaed flux"""
        M = self.P_ff(flux.time, flux.energy)
        M = (flux.flavor_scheme<<M<<flux.flavor_scheme)
        return M@flux
     
###############################################################################
class NoTransformation(FlavorTransformation):
    """Survival probabilities for no oscillation case."""

    def P_ff(self, t, E):
        r"""This transformation returns the object without transform, 
        so the transformation probability matrix is unit:
        
        .. math::
        
            P_{\alpha\beta} = I_{\alpha\beta}
        """
        p = FlavorMatrix.eye(ThreeFlavor)
        return p

    def apply(self, flux):
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

    def P_ff(self, t, E): 
        """Equal mixing so Earth matter has no effect"""
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1,f2):
            return (f1.is_neutrino==f2.is_neutrino)*1/3.
        return P

#####################################################################
class TransformationChain(FlavorTransformation):
    r"""This class calculates the probability matrix :math:`P_{\beta\alpha}` of the :math:`\nu_\alpha\to\nu_\beta` flavor transition as multiplication of :math:`P^{SN}` transformation in SN, :math:`P^{Vac}` transformation in vacuum and :math:`P^{Earth}` transformation in Earth:

    .. math:: 
    
        P_{\beta\alpha} = \sum\limits_{i,j}P^{Earth}_{\beta j} \cdot P^{Vac}_{ji} \cdot P^{SN}_{i\alpha}
    """
    def __init__(self, 
                 in_sn: SNTransformation,
                 in_vacuum: VacuumTransformation=NoVacuumTransformation(),
                 in_earth: EarthTransformation=NoEarthMatter(),
                 *,
                 mixing_params:ThreeFlavorMixingParameters|FourFlavorMixingParameters=MixingParameters()
                ):
        """
        Parameters
        ----------
        in_sn 
            Transformation in Supernova.
        in_vacuum
            Transformation in Vacuum.
            By default NoVacuumTransformation is applied
        in_earth
            Transformation in Earth.
            By default NoEarthTransformation is applied

        Keyword mixing_params
            Neutrino mixing parameters (to be applied to all individual transformations in chain)
            By default use standard `MixingParameters` with normal neutrino ordering
        """
        if in_sn is None:
            return NoTransformation()
        self.in_sn = in_sn
        self.in_vacuum = in_vacuum
        self.in_earth = in_earth
        self.transforms = [in_sn, in_vacuum, in_earth]
        self.set_mixing_params(mixing_params)
    
    def set_mixing_params(self, mixing_params):
        """Update the mixing parameters in all transformations in chain"""
        #set the mixing parameters to all the inner classes
        self.mixing_params = mixing_params
        for t in self.transforms:
            t.mixing_params = mixing_params
        
    def __call__(self, mixing_params):
        """Convenience method: update mixing parameters and return self"""
        self.set_mixing_params(mixing_params)
        return self
        
    def P_ff(self, t, E)->FlavorMatrix:
        in_sn, in_vacuum, in_earth = self.transforms
        return in_earth.P_fm(t,E) @ in_vacuum.P_mm(t,E) @ in_sn.P_mf(t,E)
    
    def __str__(self):
        s = '+'.join([t.__class__.__name__ for t in self.transforms])+'_'+self.mixing_params.mass_order.name
        return s
