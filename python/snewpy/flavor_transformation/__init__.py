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
from . import in_sn, in_earth, in_vacuum
from .transforms import NoTransformation, CompleteExchange, ThreeFlavorDecoherence, FlavorTransformation
from .transforms import TransformationChain

#define default values for backward compatibility
AdiabaticMSW = TransformationChain(in_sn.AdiabaticMSW())
NonAdiabaticMSWH = TransformationChain(in_sn.NonAdiabaticMSWH())
AdiabaticMSWes = TransformationChain(in_sn.AdiabaticMSWes())
NonAdiabaticMSWes = TransformationChain(in_sn.NonAdiabaticMSWes())
TwoFlavorDecoherence = TransformationChain(in_sn.TwoFlavorDecoherence())
NeutrinoDecay = TransformationChain(in_sn.AdiabaticMSW(), in_vacuum.NeutrinoDecay())
QuantumDecoherence = TransformationChain(in_sn.AdiabaticMSW(), in_vacuum.QuantumDecoherence())
EarthMatter = lambda mixing_params,AltAz: TransformationChain(
    in_sn.AdiabaticMSW(),
    in_earth=in_earth.EarthMatter(SNAltAz=AltAz), 
    mixing_params=mixing_params
)