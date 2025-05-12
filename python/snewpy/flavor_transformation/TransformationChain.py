# -*- coding: utf-8 -*-
"""TransformationChain

Class for calculating the combined probability matrix of the neutrino flavor transformation in SN, vacuum and Earth matter."""


from snewpy.flavor import FlavorMatrix
from snewpy.neutrino import MixingParameters, ThreeFlavorMixingParameters, FourFlavorMixingParameters
from .base import FlavorTransformation
from .in_sn import SNTransformation
from .in_vacuum import VacuumTransformation, NoVacuumTransformation
from .in_earth import EarthTransformation, NoEarthMatter
from collections import namedtuple

#a tuple to hold the transformations list
TransformationsTuple = namedtuple('TransformationsTuple',['in_sn','in_vacuum','in_earth'])

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
        self.transforms = TransformationsTuple(in_sn, in_vacuum, in_earth)
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
        return self.transforms.in_earth.P_fm(t,E) @ \
               self.transforms.in_vacuum.P_mm(t,E) @ \
               self.transforms.in_sn.P_mf(t,E)
    
    def __str__(self):
        s = '+'.join([t.__class__.__name__ for t in self.transforms])+'_'+self.mixing_params.mass_order.name
        return s
