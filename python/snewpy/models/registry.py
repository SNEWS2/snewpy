# -*- coding: utf-8 -*-
"""
A submodule with classes used for accessing supernova model files stored on disk.
It assumes models are available under a directory `data_folder` that is specified when
SNEWPY is installed under /home/user/.astropy/cache via snewpy.__init__.py.
"""


from astropy.units.quantity import Quantity
from astropy.units import UnitTypeError, get_physical_type
from snewpy import model_path, get_models
import logging
from . import ccsn, presn
import itertools as it


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
    >>> from snewpy.models.registry import init_model; import astropy.units as u
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


def check_param_values(model, **user_param):
    """Performs generic check that the requested model parameters have valid values and units for the requested
    SNEWPY model.

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
    check_param_names(model, **user_param)

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


def check_param_combo(model, **user_param):
    """Performs generic check that the model-specific combination of requested parameters is valid according
      to model.isvalid_param_combo. Valid parameter combinations may be found in model.param_combinations.

    Parameters
    ----------
    model : snewpy.model.SupernovaModel
        Model class used to perform parameter combination check
    user_param : varies
        User-requested model parameters to be tested for validity.
        NOTEL These must be provided as kwargs that match the keys of model.param

    Raises
    ------
    ValueError
        If invalid combination of model parameters are provided.

    See Also
    --------
    snewpy.models.ccsn
    snewpy.models.presn
    """
    check_param_names(model, **user_param)
    if not model.isvalid_param_combo(**user_param):
        raise ValueError(f"Invalid parameter combination. See {model.__name__}.param_combinations for a list of "
                         "allowed parameter combinations.")


def check_param_names(model, **user_param):
    if not all(key in user_param for key in model.param.keys()):
        raise ValueError(f"Missing parameter! Expected {model.param.keys()} but was given {user_param.keys()}")


def check_valid_params(model, **user_param):
    """Test for valid model parameter combination validity
        See model.__init__ for description of "Other Parameters"
    """
    # Check parameters for valid values and units
    check_param_values(model, **user_param)
    # Check parameter combination for validity (model-specific)
    check_param_combo(model, **user_param)


def get_param_combinations(model, **user_param):
    check_param_names(model, **user_param)
    param = {}
    for key, allowed_val, user_val in zip(model.param.keys(), model.param.values(), user_param.values()):
        if user_val is not None and user_val in allowed_val:
            param.update({key: [user_val]})
        elif user_val is not None and user_val not in allowed_val:
            raise ValueError(f"Invalid value for parameter {key}. Given {user_val}, but expected one from {val}")
        else:
            param.update({key: allowed_val})
    return tuple(dict(zip(param, c)) for c in it.product(*param.values()) if model.isvalid_param_combo(*c))
