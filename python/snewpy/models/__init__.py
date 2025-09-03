from warnings import warn

from . import ccsn

"""
def __getattr__(name):
    if name in dir(ccsn):
        warn(f"{__name__}.{name} is moved to {__name__}.ccsn.{name}", FutureWarning, stacklevel=2)
        return getattr(ccsn, name)
    raise AttributeError(f"module {__name__} has no attribute {name}")
"""
