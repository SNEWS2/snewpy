r""" This module implements flavor transformations that describe how neutrinos of
different flavors change into each other between production inside the
supernova and detection on Earth.

Base Class for Flavor Transformations
-------------------------------------
.. autoclass:: snewpy.flavor_transformation.FlavorTransformation
   :members:


Available Transformations
-------------------------
.. autoclass:: snewpy.flavor_transformation.NoTransformation
    :members:

.. autoclass:: snewpy.flavor_transformation.CompleteExchange
    :members:

.. autoclass:: snewpy.flavor_transformation.ThreeFlavorDecoherence
    :members:

.. autoclass:: snewpy.flavor_transformation.TransformationChain
    :members:

Concrete transformations
------------------------
.. automodule:: snewpy.flavor_transformation.in_sn
   :members:
   :member-order: bysource

.. automodule:: snewpy.flavor_transformation.in_vacuum
   :members:
   :member-order: bysource
   
.. automodule:: snewpy.flavor_transformation.in_earth
   :members:
   :member-order: bysource
"""
from snewpy.flavor import FlavorMatrix, ThreeFlavor
from . import in_sn, in_earth, in_vacuum
from .base import FlavorTransformation
from .TransformationChain import TransformationChain
from snewpy.neutrino import ThreeFlavorMixingParameters, FourFlavorMixingParameters

def construct_chain(*transformation_classes):
    """Create a function that constructs a chain from given transformation_steps"""
    def construct(mixing_params:ThreeFlavorMixingParameters|FourFlavorMixingParameters,
                 )->TransformationChain:
        """Construct a transformation chain with 
        Parameters
        ----------
        mixing_params: 
            neutrino mixing parameters to be passed to the internal transformation steps
        
        Returns
        -------
            A TransformationChain object with {step_names} transformations
        """
        #initialize the transformations
        transformations = [value() for value in transformation_classes]
        #create the chain
        return TransformationChain(*transformations, mixing_params=mixing_params)
    #update the docstring
    step_names=[cls.__qualname__ for cls in transformation_classes]
    construct.__doc__ = construct.__doc__.format(step_names=step_names)
    return construct
 
# define default values for backward compatibility
AdiabaticMSW = construct_chain(in_sn.AdiabaticMSW)
NonAdiabaticMSWH = construct_chain(in_sn.NonAdiabaticMSWH)
AdiabaticMSWes = construct_chain(in_sn.AdiabaticMSWes)
NonAdiabaticMSWes = construct_chain(in_sn.NonAdiabaticMSWes)
TwoFlavorDecoherence = construct_chain(in_sn.TwoFlavorDecoherence)
NeutrinoDecay = construct_chain(in_sn.AdiabaticMSW, 
                                in_vacuum.NeutrinoDecay)
QuantumDecoherence = construct_chain(in_sn.AdiabaticMSW, 
                                     in_vacuum.QuantumDecoherence)
EarthMatter = lambda mixing_params,AltAz: TransformationChain(
    in_sn.AdiabaticMSW(),
    in_earth=in_earth.EarthMatter(SNAltAz=AltAz),
    mixing_params=mixing_params
)


# Phenomenological transformations that cannot be represented as a TransformationChain
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

    def apply_to(self, flux):
        return flux


class CompleteExchange(FlavorTransformation):
    """Survival probabilities for the case when the electron flavors
       are half exchanged with the mu flavors and the half with the tau flavors.
    """

    def P_ff(self, t, E):
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1, f2):
            return (f1.is_neutrino == f2.is_neutrino)*(f1 != f2)*0.5

        return P


class ThreeFlavorDecoherence(FlavorTransformation):
    """Equal mixing of all threen eutrino matter states and antineutrino matter states"""

    def P_ff(self, t, E):
        """Equal mixing so Earth matter has no effect"""
        @FlavorMatrix.from_function(ThreeFlavor)
        def P(f1, f2):
            return (f1.is_neutrino == f2.is_neutrino)*1/3.
        return P
