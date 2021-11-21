from . import base,ccsn,presn
import importlib
from warnings import warn

def __getattr__(name):
    if name in dir(ccsn):
        warn(f"{__name__}.{name} is moved to {__name__}.ccsn.{name}")
        return getattr(ccsn,name)
    raise AttributeError(f"module {__name__} has no attribute {name}")
