# -*- coding: utf-8 -*-
"""
A submodule with classes used for accessing supernova model files stored on disk.
It assumes models are available under a directory `data_folder` that is specified when
SNEWPY is installed under /home/user/.astropy/cache via snewpy.__init__.py.
"""


from astropy.units.quantity import Quantity
from astropy.units import UnitTypeError, get_physical_type


def check_param_values(model, **user_param):
    """Performs generic check for a model that has an existing entry in the registry that the use provided args
    and kwargs have valid values and units (if applicable) for the specified model.

    Parameters
    ----------
    model : snewpy.model.SupernovaModel
        Name of model entry in registry, it should match exactly with a SNEWPY model class name
    user_param : varies
        Model parameters provided with keywords, type varies based on model entry
        It is strongly recommended that the user provides keyword arguments to this function other than ``model``

    Raises
    ------
    ValueError
        If invalid model parameters are provided based on units, allowed values, etc.
    UnitTypeError
        If invalid units are provided for a model parameter

    See Also
    --------
    snewpy.models

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
