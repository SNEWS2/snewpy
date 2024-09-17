from . import in_sn, in_earth, in_vacuum
from .transforms import NoTransformation, CompleteExchange, ThreeFlavorDecoherence
from .transforms import TransformationChain

#define default values
AdiabaticMSW = TransformationChain(in_sn.AdiabaticMSW())
NonAdiabaticMSWH = TransformationChain(in_sn.NonAdiabaticMSWH())
AdiabaticMSWes = TransformationChain(in_sn.AdiabaticMSWes())
NonAdiabaticMSWHes = TransformationChain(in_sn.NonAdiabaticMSWes())