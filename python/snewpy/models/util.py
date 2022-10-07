# -*- coding: utf-8 -*-
"""
A submodule with classes used for accessing supernova model files stored on disk.
It assumes models are available under a directory `data_folder` that is specified when
SNEWPY is installed under /home/user/.astropy/cache via snewpy.__init__.py.
"""


from astropy.units.quantity import Quantity
from astropy.units import UnitTypeError, get_physical_type
from snewpy import model_path, get_models
import numpy as np
import logging
from . import ccsn, presn
import itertools as it


def get_param_combinations(param, func_isvalid=None):
    """Returns all valid combinations of parameters for a given SNEWPY register model.

    Parameters
    ----------
    param : dict
        Dictionary of SNEWPY model parameter values.
    func_isvalid : callable or None
        Callable that acts upon argument param that returns True if a particular combinations of parameters is valid.
        If None is provided, all combinations are considered valid

    Returns
    -------
    valid_combinations: tuple[dict]
        A tuple of all valid parameter combinations stored as Dictionaries
    """
    # Input sanitization
    param = dict(param)
    for key, val in param.items():
        if not isinstance(val, (list, Quantity)):
            param[key] = [val]
        elif isinstance(val, Quantity) and val.size == 1:
            param[key] = [val]
    combos = tuple(dict(zip(param, combo)) for combo in it.product(*param.values()))
    return tuple(c for c in filter(func_isvalid, combos))


def init_model(model_name, download=True, download_dir=model_path, **user_param):
    """Attempts to retrieve instantiated SNEWPY model using model class name and model parameters.
    If a model name is valid, but is not found and `download`=True, this function will attempt to download the model

    Parameters
    ----------
    model_name : str
        Name of SNEWPY model to import, must exactly match the name of the corresponding model class
    download : bool
        Switch for attempting to download model data if the first load attempt failed due to a missing file.
    download_dir : str
        Local directory to download model files to.
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
    >>> from snewpy.models.util import init_model; import astropy.units as u
    >>> init_model('Nakazato_2013', progenitor_mass=13*u.Msun, metallicity=0.004, revival_time=0*u.s, eos='shen')
    Nakazato_2013 Model: nakazato-shen-BH-z0.004-s30.0.fits
    Progenitor mass  : 30.0 solMass
    EOS              : Shen
    Metallicity      : 0.004
    Revival time     : 0.0 ms
    """
    if model_name in dir(ccsn):
        module = ccsn
    elif model_name in dir(presn):
        module = presn
    else:
        raise ValueError(f"Unable to find model with name '{model_name}' in snewpy.models.ccsn or snewpy.models.presn")

    try:
        return getattr(module, model_name)(**user_param)
    except FileNotFoundError as e:
        logger = logging.getLogger()
        logger.warning(f"Unable to find model {model_name} in {download_dir}")
        if not download:
            raise e
        logger.warning(f"Attempting to download model...")
        get_models(model_name, download_dir)
        return getattr(module, model_name)(**user_param)


def check_valid_params(model, **user_params):
    """Checks that the model-specific values, units, names and conbinations of requested parameters are valid.

    Parameters
    ----------
    model : snewpy.model.SupernovaModel
        Model class used to perform parameter check
    user_params : varies
        User-requested model parameters to be tested for validity.
        NOTE: This must be provided as kwargs that match the keys of model.param

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

    # Check that the appropriate number of params are provided
    if not all(key in user_params for key in model.param.keys()):
        raise ValueError(f"Missing parameter! Expected {model.param.keys()} but was given {user_params.keys()}")

    # Check parameter units and values
    for (key, allowed_params), user_param in zip(model.param.items(), user_params.values()):

        # If both have units, check that the user param value is valid. If valid, continue. Else, error
        if type(user_param) == Quantity and type(allowed_params) == Quantity:
            if get_physical_type(user_param.unit) != get_physical_type(allowed_params.unit):
                raise UnitTypeError(f"Incorrect units {user_param.unit} provided for parameter {key}, "
                                    f"expected {allowed_params.unit}")

            elif np.isin(user_param.to(allowed_params.unit).value, allowed_params.value):
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

    # Check Combinations (Logic lives inside model subclasses under model.isvalid_param_combo)
    if user_params not in model.param_combinations:
        raise ValueError(
            f"Invalid parameter combination. See {model.__class__.__name__}.param_combinations for a "
            "list of allowed parameter combinations.")
