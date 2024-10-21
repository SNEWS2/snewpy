from . import in_sn, in_earth, in_vacuum
from .transforms import NoTransformation, CompleteExchange, ThreeFlavorDecoherence
from .transforms import TransformationChain
from snewpy.neutrino import MassHierarchy, MixingParameters

#define default values for backward compatibility
AdiabaticMSW = TransformationChain(in_sn.AdiabaticMSW())
NonAdiabaticMSWH = TransformationChain(in_sn.NonAdiabaticMSWH())
AdiabaticMSWes = TransformationChain(in_sn.AdiabaticMSWes())
NonAdiabaticMSWes = TransformationChain(in_sn.NonAdiabaticMSWes())
TwoFlavorDecoherence = TransformationChain(in_sn.TwoFlavorDecoherence())
NeutrinoDecay = TransformationChain(in_sn.AdiabaticMSW(), in_vacuum.NeutrinoDecay())