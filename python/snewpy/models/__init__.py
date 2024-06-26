from warnings import warn

from . import ccsn, presn


def __getattr__(name):
    if name in dir(ccsn):
        warn(f"{__name__}.{name} is moved to {__name__}.ccsn.{name}", FutureWarning, stacklevel=2)
        return getattr(ccsn, name)
    raise AttributeError(f"module {__name__} has no attribute {name}")


def _init_model(model_name, **user_param):
    """Attempts to retrieve instantiated SNEWPY model using model class name and model parameters.

    Parameters
    ----------
    model_name : str
        Name of SNEWPY model to import, must exactly match the name of the corresponding model class
    user_param : varies
        User-requested model parameters used to initialize the model, if one is found.
        Error checking is performed during model initialization

    Raises
    ------
    ValueError
        If the requested model_name does not match any SNEWPY models

    See Also
    --------
    snewpy.models.ccsn
    snewpy.models.presn

    Example
    -------
    >>> from snewpy.models import _init_model; import astropy.units as u
    >>> _init_model('Nakazato_2013', progenitor_mass=13*u.Msun, metallicity=0.004, revival_time=0*u.s, eos='shen')
    Nakazato_2013 Model: nakazato-shen-BH-z0.004-s30.0.fits
    Progenitor mass  : 30.0 solMass
    EOS              : Shen
    Metallicity      : 0.004
    Revival time     : 0.0 ms

    :meta private:
    """
    if model_name in dir(ccsn):
        module = ccsn
    elif model_name in dir(presn):
        module = presn
    else:
        raise ValueError(f"Unable to find model with name '{model_name}' in snewpy.models.ccsn or snewpy.models.presn")

    return getattr(module, model_name)(**user_param)
