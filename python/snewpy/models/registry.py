# -*- coding: utf-8 -*-
"""
A submodule with classes used for accessing supernova model files stored on disk.
It assumes models are available under a directory `data_folder` that is specified when
SNEWPY is installed under /home/user/.astropy/cache via snewpy.__init__.py.
"""


from astropy.units.quantity import Quantity
from astropy.units import UnitTypeError, get_physical_type
from . import ccsn, presn


def get_model(model_name, **user_param):
    """Attempts to retrieve instantiated SNEWPY model using model class name and model parameters

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
    >>> from snewpy.models.registry import get_model; import astropy.units as u
    >>> get_model('Nakazato_2013', progenitor_mass=13*u.Msun, metallicity=0.004, revival_time=0*u.s, eos='shen')
    Nakazato_2013 Model: nakazato-shen-BH-z0.004-s30.0.fits
    Progenitor mass  : 30.0 solMass
    EOS              : Shen
    Metallicity      : 0.004
    Revival time     : 0.0 ms
    """
    if model_name in dir(ccsn):
        return getattr(ccsn, model_name)(**user_param)
    elif model_name in dir(presn):
        return getattr(presn, model_name)(**user_param)
    else:
        raise ValueError(f"Unable to find model with name '{model_name}'")


def check_param_values(model, **user_param):
    """Performs generic check for that the requested model parameters have valid values and units for the requested
    SNEWPY model. Model arguments MUST be provided as keyword arguments that match the model_class `param` class member.

    Parameters
    ----------
    model : snewpy.model.SupernovaModel
        Model class used to perform parameter check
    user_param : varies
        User-requested model parameters to be tested for validity. MUST be provided as keyword arguments that match the
        model `param` class member

    Raises
    ------
    ValueError
        If invalid model parameters are provided based on units, allowed values, etc.
    UnitTypeError
        If invalid units are provided for a model parameter

    See Also
    --------
    snewpy.models.ccsn
    snewpy.models.presn

    """
    model_param = model.param
    # Check that the appropriate number of params are provided
    if len(user_param) != len(model_param):
        raise ValueError(f"Invalid model parameters, expected {len(model_param)} "
                         f"but {len(user_param)} were given")

    # Check that user-requested params have valid units and values
    for (key, allowed_params), user_param in zip(model_param.items(), user_param.values()):
        # If both have units, check that the user param value is valid. If valid, continue. Else, error
        if type(user_param) == Quantity and type(allowed_params) == Quantity:
            if get_physical_type(user_param.unit) != get_physical_type(allowed_params.unit):
                raise UnitTypeError(f"Incorrect units {user_param.unit} provided for parameter {key}, "
                                    f"expected {allowed_params.unit}")
            elif user_param.to(allowed_params.unit).value in allowed_params.value:
                continue
            else:
                raise ValueError(f"Invalid value '{user_param}' provided for parameter {key}, "
                                 f"allowed value(s): {allowed_params}")

        # If one only one has units, then error
        elif (type(user_param) == Quantity) ^ (type(allowed_params) == Quantity):
            # User param has units, model param is unitless
            if type(user_param) == Quantity:
                raise ValueError(f"Invalid units {user_param.unit} for parameter {key} provided, expected None")
            else:
                raise ValueError(f"Missing units for parameter {key}, expected {allowed_params.unit}")

        # Check that unitless user param value is valid. If valid, continue. Else, Error
        elif user_param in allowed_params:
            continue
        else:
            raise ValueError(f"Invalid value '{user_param}' provided for parameter {key}, "
                             f"allowed value(s): {allowed_params}")
