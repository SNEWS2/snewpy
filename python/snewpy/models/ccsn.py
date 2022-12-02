# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.ccsn`` contains models of core-collapse supernovae
derived from the :class:`SupernovaModel` base class.

Since SNEWPY v1.3, the prefered method is to initialise a model based on its physics parameters.
Initialisation from a file name is deprecated.

Use the ``param`` class property to view all physics parameters and their possible values:

>>> from snewpy.models.ccsn import Nakazato_2013
>>> Nakazato_2013.param
{'progenitor_mass': <Quantity [13., 20., 30., 50.] solMass>,
 'revival_time': <Quantity [  0., 100., 200., 300.] ms>,
 'metallicity': [0.02, 0.004],
 'eos': ['LS220', 'shen', 'togashi']}

For some models, not all combinations of parameters are valid. Use the ``get_param_combinations()``
class method to get a list of all valid combinations and filter it:

>>> list(params for params in Nakazato_2013.get_param_combinations() if params['eos'] != 'shen')
[{'progenitor_mass': <Quantity 30. solMass>, 'revival_time': <Quantity 0. ms>,
  'metallicity': 0.004, 'eos': 'LS220'},
 {'progenitor_mass': <Quantity 30. solMass>, 'revival_time': <Quantity 0. ms>,
  'metallicity': 0.004, 'eos': 'togashi'}]

.. _Garching Supernova Archive: https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/
"""
import logging
import os
import re
import tarfile

import numpy as np
from astropy import units as u
from astropy.table import Table

from snewpy.models import loaders
from .base import PinchedModel, _RegistryModel

from warnings import warn
from functools import wraps


def _warn_deprecated_filename_argument(func):
    """Decorator for model.__new__ methods, causing them to issue a deprecation warning
     if argument `filename` is used. Initialization by filename will be moved to
     snewpy.models.init_model, and model classes are initialized from physical parameters.
    """
    @wraps(func)  # Ensures docstrings are preserved
    def decorator(cls, *args, **kwargs):
        filename = args[0] if len(args) > 0 else None  # Assumes filename is first pos. arg if provided
        if filename is not None:
            msg = ''.join(['Initializing this model with a filename is deprecated. ',
                           f'Instead, use keyword arguments {list(cls.param.keys())}. ',
                           f'See `{cls.__name__}.param`, `{cls.__name__}.get_param_combinations()` for more info.'])
            warn(FutureWarning(msg), stacklevel=2)
        return func(cls, *args, **kwargs)
    return decorator


class Analytic3Species(PinchedModel):
    """An analytical model calculating spectra given total luminosity,
    average energy, and rms or pinch, for each species.
    """

    param = "There are no input files available for this class. Use `doc/scripts/Analytic.py` in the SNEWPY GitHub repo to create a custom input file."

    def get_param_combinations(cls):
        print(cls.param)
        return []

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data.
        """

        simtab = Table.read(filename,format='ascii')
        self.filename = filename
        super().__init__(simtab, metadata={})


class Nakazato_2013(_RegistryModel):
    """Model based on simulations from Nakazato et al., ApJ S 205:2
    (2013), ApJ 804:75 (2015), PASJ 73:639 (2021). See also http://asphwww.ph.noda.tus.ac.jp/snn/.
    """

    param = {'progenitor_mass': [13, 20, 30, 50] * u.Msun,
             'revival_time': [0, 100, 200, 300] * u.ms,
             'metallicity': [0.02, 0.004],
             'eos': ['LS220', 'shen', 'togashi']}

    _param_validator = lambda p: (p['revival_time'] == 0 * u.ms and p['progenitor_mass'] == 30 * u.Msun
                                  and p['metallicity'] == 0.004) or \
                                 (p['revival_time'] != 0 * u.ms and p['eos'] == 'shen'
                                  and not (p['progenitor_mass'] == 30 * u.Msun and p['metallicity'] == 0.004))

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, *, progenitor_mass=None, revival_time=None, metallicity=None, eos=None):
        """Model initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data. This argument is deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.
        revival_time: astropy.units.Quantity
            Time of shock revival in model in units ms. Valid values are {revival_time}.
            Selecting 0 ms will load a black hole formation model
        metallicity: float
            Progenitor metallicity. Valid values are {metallicity}.
        eos: str
            Equation of state. Valid values are {eos}.

        Raises
        ------
        ValueError
            If a combination of parameters is invalid when loading from parameters

        Examples
        --------
        >>> from snewpy.models.ccsn import Nakazato_2013; import astropy.units as u
        >>> Nakazato_2013(progenitor_mass=30*u.Msun, metallicity=0.004, revival_time=0*u.s, eos='togashi')
        Nakazato_2013 Model: nakazato-togashi-BH-z0.004-s30.0.fits
        Progenitor mass  : 30.0 solMass
        EOS              : togashi
        Metallicity      : 0.004
        Revival time     : 0.0 ms
        """
        # Attempt to load model from parameters
        if filename is not None:
            metadata = {'Progenitor mass': float(filename.split('-')[-1].strip('s%.fits')) * u.Msun}
            if 't_rev' in filename:
                metadata.update({
                    'EOS': filename.split('-')[-4].upper(),
                    'Metallicity': float(filename.split('-')[-3].strip('z%')),
                    'Revival time': float(filename.split('-')[-2].strip('t_rev%ms')) * u.ms
                })
            # No revival time because the explosion "failed" (BH formation).
            else:
                metadata.update({
                    'EOS': filename.split('-')[-4].upper(),
                    'Metallicity': float(filename.split('-')[-2].strip('z%')),
                    'Revival time': 0 * u.ms
                })
            return loaders.Nakazato_2013(os.path.abspath(filename), metadata)

        # Load from model parameters
        user_params = dict(zip(cls.param.keys(), (progenitor_mass, revival_time, metallicity, eos)))
        cls.check_valid_params(cls, **user_params)

        # Store model metadata.
        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos,
            'Metallicity': metallicity,
            'Revival time': revival_time
        }

        # Strip units for filename construction
        progenitor_mass = progenitor_mass.to(u.Msun).value
        revival_time = revival_time.to(u.ms).value

        if revival_time != 0:
            filename = f"nakazato-{eos}-z{metallicity}-t_rev{int(revival_time)}ms-s{progenitor_mass:3.1f}.fits"
        else:
            filename = f"nakazato-{eos}-BH-z{metallicity}-s{progenitor_mass:3.1f}.fits"

        return loaders.Nakazato_2013(filename, metadata)

    # Populate Docstring with param values
    __new__.__doc__ = __new__.__doc__.format(**param)


class Sukhbold_2015(_RegistryModel):
    """Model based on simulations from Sukhbold et al., ApJ 821:38,2016. Models were shared privately by email.
    """
    param = {'progenitor_mass': [27., 9.6] * u.Msun,
             'eos': ['LS220', 'SFHo']}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, *, progenitor_mass=None, eos=None):
        """Model Initialization

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data. This argument is deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.
        eos: str
            Equation of state. Valid values are {eos}.

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters
        """
        if filename is not None:
            metadata = {
                'Progenitor mass': float(filename.split('-')[-1].strip('z%.fits')) * u.Msun,
                'EOS': filename.split('-')[-2]
            }
            return loaders.Sukhbold_2015(os.path.abspath(filename), metadata)

        user_params = dict(zip(cls.param.keys(), (progenitor_mass, eos)))
        cls.check_valid_params(cls, **user_params)

        if progenitor_mass.value == 9.6:
            filename = f'sukhbold-{eos}-z{progenitor_mass.value:3.1f}.fits'
        else:
            filename = f'sukhbold-{eos}-s{progenitor_mass.value:3.1f}.fits'

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos
        }
        return loaders.Sukhbold_2015(filename, metadata)

    # Populate Docstring with param values
    __new__.__doc__ = __new__.__doc__.format(**param)


class Tamborra_2014(_RegistryModel):
    """Model based on 3D simulations from `Tamborra et al., PRD 90:045032, 2014 <https://arxiv.org/abs/1406.0006>`_.
    Data files are from the `Garching Supernova Archive`_.
    """

    param = {'progenitor_mass': [20., 27.] * u.Msun}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='LS220', *, progenitor_mass=None):
        if filename is not None:
            # Metadata creation is implemented in snewpy.models.base._GarchingArchiveModel
            return loaders.Tamborra_2014(os.path.abspath(filename))

        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        filename = f's{progenitor_mass.value:3.1f}c_3D_dir1'

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': 'LS220'
        }

        # Metadata is handled by __init__ in _GarchingArchiveModel
        return loaders.Tamborra_2014(filename=filename, metadata=metadata)

    # Populate Docstring with param values
    __new__.__doc__ = loaders.Tamborra_2014.__init__.__doc__.format(**param)


class Bollig_2016(_RegistryModel):
    """Model based on simulations from `Bollig et al. (2016) <https://arxiv.org/abs/1508.00785>`_.
    Models were taken, with permission, from the `Garching Supernova Archive`_.
    """

    param = {'progenitor_mass': [11.2, 27.] * u.Msun}

    def __new__(cls, filename=None,  eos='LS220', *, progenitor_mass=None):
        if filename is not None:
            # Metadata creation is implemented in snewpy.models.base._GarchingArchiveModel
            return loaders.Bollig_2016(os.path.abspath(filename))

        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        filename = f's{progenitor_mass.value:3.1f}c'

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': 'LS220'
        }

        return loaders.Bollig_2016(filename=filename, metadata=metadata)

    # Populate Docstring with param values (Docstring is inherited from base._GarchingArchiveModel.__init__)
    __new__.__doc__ = loaders.Bollig_2016.__init__.__doc__.format(**param)


class Walk_2018(_RegistryModel):
    """Model based on SASI-dominated simulations from `Walk et al.,
    PRD 98:123001, 2018 <https://arxiv.org/abs/1807.02366>`_. Data files are from
    the `Garching Supernova Archive`_.
    """

    param = {'progenitor_mass': 15. * u.Msun}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='LS220', *, progenitor_mass=None):
        if filename is not None:
            # Metadata creation is implemented in snewpy.models.base._GarchingArchiveModel
            return loaders.Walk_2018(os.path.abspath(filename))

        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        filename = f's{progenitor_mass.value:3.1f}c_3D_nonrot_dir1'

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': 'LS220'
        }

        return loaders.Walk_2018(filename=filename, metadata=metadata)

    # Populate Docstring with param values (Docstring is inherited from base._GarchingArchiveModel.__init__)
    __new__.__doc__ = loaders.Walk_2018.__init__.__doc__.format(**param)


class Walk_2019(_RegistryModel):
    """Model based on SASI-dominated simulations from `Walk et al.,
    PRD 101:123013, 2019 <https://arxiv.org/abs/1910.12971>`_. Data files are
    from the `Garching Supernova Archive`_.
    """

    param = {'progenitor_mass': 40 * u.Msun}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='LS220', *, progenitor_mass=None):
        if filename is not None:
            # Metadata creation is implemented in snewpy.models.base._GarchingArchiveModel
            return loaders.Walk_2019(os.path.abspath(filename))

        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        filename = f's{progenitor_mass.value:3.1f}c_3DBH_dir1'

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': 'LS220'
        }

        return loaders.Walk_2019(filename=filename, metadata=metadata)

    # Populate Docstring with param values (Docstring is inherited from base._GarchingArchiveModel.__init__)
    __new__.__doc__ = loaders.Walk_2019.__init__.__doc__.format(**param)


class OConnor_2013(_RegistryModel):
    """Model based on the black hole formation simulation in `O'Connor & Ott (2013) <https://arxiv.org/abs/1207.1100>`_.
    """

    param = {'progenitor_mass': (list(range(12, 34)) +
                                 list(range(35, 61, 5)) +
                                 [70, 80, 100, 120]) * u.Msun,
             'eos': ['HShen', 'LS220']}

    _param_abbrv = {'progenitor_mass': '[12..33, 35..5..60, 70, 80, 100, 120] solMass',
                    'eos': ['HShen', 'LS220']}

    @_warn_deprecated_filename_argument
    def __new__(cls, base=None, mass=None, eos='LS220', *, progenitor_mass=None):
        """Model Initialization.

        Parameters
        ----------
        base : str
            Absolute or relative path folder with model data. This argument is deprecated.
        mass: int
            Mass of model progenitor in units Msun. This argument is deprecated.
        eos: str
            Equation of state. Valid values are {eos}.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters

        """
        # TODO: (For v2.0) Change `base` to filename, move compressed model files to OConnor_2013 model folder
        if mass is not None:
            warn(f'Argument `mass` of type int is deprecated. To initialize this model, use keyword arguments '
                 f'{list(cls.param.keys())}. See {cls.__name__}.param, {cls.__name__}.get_param_combinations() for more info',
                 category=DeprecationWarning, stacklevel=2)
        else:
            mass = 15  # Default Value, this is handled this way for backwards compatibility -- TODO (For V2.0) Remove

        if base is not None:
            # If base is provided, do not attempt to load from param.
            if mass * u.Msun not in cls.param['progenitor_mass']:
                raise ValueError(f'Invalid value for argument `progenitor mass` or `mass`, see {cls.__name__}.param'
                                 f' for allowed values')
            metadata = {'Progenitor mass': progenitor_mass if progenitor_mass is not None else mass * u.Msun,
                        'EOS': eos}

            filename = os.path.join(base, f"{eos}_timeseries.tar.gz")
            return loaders.OConnor_2013(os.path.abspath(filename), metadata)

        # Load from Parameters
        cls.check_valid_params(cls, progenitor_mass=progenitor_mass, eos=eos)
        filename = f'{eos}_timeseries.tar.gz'

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos,
        }

        return loaders.OConnor_2013(filename=filename, metadata=metadata)

    # Populate Docstring with param values
    __new__.__doc__ = __new__.__doc__.format(**_param_abbrv)


class OConnor_2015(_RegistryModel):
    """Model based on the black hole formation simulation in `O'Connor (2015) <https://arxiv.org/abs/1411.7058>`_.
    """

    param = {'progenitor_mass': 40 * u.Msun}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='LS220', *, progenitor_mass=None):
        """Model Initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to tar.gz file with model data. This argument is deprecated.
        eos: str
            Equation of state. Valid value is 'LS220'. This argument is deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters
        """
        metadata = {
            'Progenitor mass': 40*u.Msun,
            'EOS': 'LS220',
        }

        if filename is not None:
            return loaders.OConnor_2015(os.path.abspath(filename), metadata)

        # Load from Parameters
        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        # Filename is currently the same regardless of parameters
        filename = 'M1_neutrinos.dat'

        return loaders.OConnor_2015(filename, metadata)

    # Populate Docstring with param values
    __new__.__doc__ = __new__.__doc__.format(**param)


class Zha_2021(_RegistryModel):
    """Model based on the hadron-quark phse transition models from `Zha et al. 2021 <https://arxiv.org/abs/2103.02268>`_.
    """

    param = {'progenitor_mass': (list(range(16, 27)) + [19.89, 22.39, 30, 33]) * u.Msun}

    _param_abbrv = {'progenitor_mass': '[16..26, 19.89, 22.39, 30, 33] solMass'}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='ST0S_B145', *, progenitor_mass=None):
        """Model Initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data. This argument is deprecated.
        eos : str
            Equation of state. Valid value is 'ST0S_B145'. This argument is deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters
        """
        if filename is not None:
            metadata = {'Progenitor mass': float(os.path.splitext(os.path.basename(filename))[0][1:]) * u.Msun,
                        'EOS': 'STOS_B145'}
            return loaders.Zha_2021(os.path.abspath(filename), metadata)

        # Load from Parameters
        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': 'STOS_B145',
        }

        filename = f's{progenitor_mass.value:g}.dat'

        return loaders.Zha_2021(filename, metadata)

    # Populate Docstring with abbreviated param values
    __new__.__doc__ = __new__.__doc__.format(**_param_abbrv)


class Warren_2020(_RegistryModel):
    """Model based on simulations from Warren et al., ApJ 898:139, 2020.
    Neutrino fluxes available at https://doi.org/10.5281/zenodo.3667908."""
    # TODO: (For v2.0) Resolve Zenodo issues (Missing files)
    # np.arange with decimal increments can produce floating point errors
    # Though it may be more intuitive to use np.arange, these fp-errors quickly become troublesome
    param = {'progenitor_mass': np.concatenate((np.linspace(9.0, 12.75, 16),
                                                np.linspace(13, 30., 171),
                                                np.linspace(31., 33., 3),
                                                np.linspace(35, 55, 5),
                                                np.linspace(60, 80, 3),
                                                np.linspace(100, 120, 2))) * u.Msun,
             'turbmixing_param': [1.23, 1.25, 1.27]}

    _param_abbrv = {'progenitor_mass': '[9..0.25..13, 13..0.1..30, 31..35, 35..5..60, 70..10..90, 100, 120] solMass',
                    'turbmixing_param': [1.23, 1.25, 1.27]}


    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='SFHo', *, progenitor_mass=None, turbmixing_param=None):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data. This argument is deprecated.
        eos : str
            Equation of state. Valid value is 'SFHo'. This argument is deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.
        turbmixing_param: float
            Turbulent mixing parameter alpha_lambda. Valid Values are {turbmixing_param}

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters
        """
        if filename is not None:
            _, _, turbmixing_param, progenitor_mass = os.path.splitext(os.path.basename(filename))[0].split('_')
            metadata = {'Progenitor mass': float(progenitor_mass[1:]) * u.Msun,
                        'Turb. mixing param.': float(turbmixing_param[1:]),
                        'EOS': 'SFHo'}

            return loaders.Warren_2020(os.path.abspath(filename), metadata)

        # Load from Parameters
        user_params = dict(zip(cls.param.keys(), (progenitor_mass, turbmixing_param)))
        cls.check_valid_params(cls, **user_params)

        if progenitor_mass.value.is_integer() and progenitor_mass.value <= 30.:
            fname = os.path.join(f'stir_a{turbmixing_param:3.2f}',
                                 f'stir_multimessenger_a{turbmixing_param:3.2f}_m{progenitor_mass.value:.1f}.h5')
        else:
            fname = os.path.join(f'stir_a{turbmixing_param:3.2f}',
                                 f'stir_multimessenger_a{turbmixing_param:3.2f}_m{progenitor_mass.value:g}.h5')

        # Set model metadata.
        metadata = {
            'Progenitor mass': progenitor_mass,
            'Turb. mixing param.': turbmixing_param,
            'EOS': 'SFHo',
        }

        return loaders.Warren_2020(fname, metadata)

    # Populate Docstring with abbreviated param values
    __new__.__doc__ = __new__.__doc__.format(**_param_abbrv)


class Kuroda_2020(_RegistryModel):
    """Model based on simulations from `Kuroda et al. (2020) <https://arxiv.org/abs/2009.07733>`_."""

    param = {'rotational_velocity': [0, 1] * u.rad / u.s,
             'magnetic_field_exponent': [0, 12, 13]}
    _param_validator = lambda p: (p['rotational_velocity'].value == 1 and p['magnetic_field_exponent'] in (12, 13)) or \
                               (p['rotational_velocity'].value == 0 and p['magnetic_field_exponent'] == 0)

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, eos='LS220', mass=20*u.Msun, *, rotational_velocity=None,
                magnetic_field_exponent=None):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data. This argument is deprecated.
        eos : str
            Equation of state. Valid value is 'LS220'. This argument is deprecated.
        mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid value is 20 * u.Msun. This argument is deprecated.

        Other Parameters
        ----------------
        rotational_velocity: astropy.units.Quantity
            Rotational velocity of progenitor. Valid values are {rotational_velocity}
        magnetic_field_exponent: int
            Exponent of magnetic field (See Eq. 46). Valid Values are {magnetic_field_exponent}

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters
        """
        if filename is not None:
            _, rotational_velocity, magnetic_field_exponent = re.split('R|B',
                                                                       os.path.splitext(os.path.basename(filename))[0])
            metadata = {
                'Progenitor mass': 20 * u.Msun,
                'Rotational Velocity': int(rotational_velocity[0]),
                'B_0 Exponent': int(magnetic_field_exponent),
                'EOS': 'LS220'
            }
            return loaders.Kuroda_2020(os.path.abspath(filename), metadata)

        # Load from Parameters
        cls.check_valid_params(cls, rotational_velocity=rotational_velocity,
                               magnetic_field_exponent=magnetic_field_exponent)
        filename = f'LnuR{int(rotational_velocity.value):1d}0B{int(magnetic_field_exponent):02d}.dat'

        metadata = {
            'Progenitor mass': 20 * u.Msun,
            'Rotational Velocity': rotational_velocity,
            'B_0 Exponent': magnetic_field_exponent,
            'EOS': 'LS220',
            }

        return loaders.Kuroda_2020(filename, metadata)

    __new__.__doc__ = __new__.__doc__.format(**param)


class Fornax_2019(_RegistryModel):
    """Model based on 3D simulations from D. Vartanyan, A. Burrows, D. Radice, M.  A. Skinner and J. Dolence, MNRAS 482(1):351, 2019. 
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.3d/
    """
    param = {'progenitor_mass': [9, 10, 12, 13, 14, 15, 16, 19, 25, 60] * u.Msun}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, cache_flux=False, *, progenitor_mass=None, ):
        """Model Initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data. This argument is deprecated.
        cache_flux : bool
            If true, pre-compute the flux on a fixed angular grid and store the values in a FITS file.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.
        """
        if filename is not None:
            progenitor_mass = os.path.splitext(os.path.basename(filename))[0].split('_')[2]
            metadata = {'Progenitor mass': int(progenitor_mass[:-1]) * u.Msun}
            return loaders.Fornax_2019(os.path.abspath(filename), metadata, cache_flux=cache_flux)

        # Load from Parameters
        metadata = {'Progenitor mass': progenitor_mass}

        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        if progenitor_mass.value == 16:
            filename = f'lum_spec_{int(progenitor_mass.value):d}M_r250.h5'
        else:
            filename = f'lum_spec_{int(progenitor_mass.value):d}M.h5'

        return loaders.Fornax_2019(filename, metadata, cache_flux=cache_flux)

    # Populate Docstring with abbreviated param values
    __new__.__doc__ = __new__.__doc__.format(**param)


class Fornax_2021(_RegistryModel):
    """Model based on 3D simulations from D. Vartanyan, A. Burrows, D. Radice, M.  A. Skinner and J. Dolence, MNRAS 482(1):351, 2019. 
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.3d/
        """
    param = {'progenitor_mass': (list(range(12, 24)) + [25, 26, 26.99]) * u.Msun}

    _param_abbrv = {'progenitor_mass': '[12..26, 26.99] solMass'}

    @_warn_deprecated_filename_argument
    def __new__(cls, filename=None, *, progenitor_mass=None):
        """Model Initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data. This argument is deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.
        """
        if filename is not None:
            progenitor_mass = os.path.splitext(os.path.basename(filename))[0].split('_')[2]
            metadata = {'Progenitor mass': float(progenitor_mass[:-1]) * u.Msun}
            return loaders.Fornax_2021(os.path.abspath(filename), metadata)

        # Load from Parameters
        cls.check_valid_params(cls, progenitor_mass=progenitor_mass)
        if progenitor_mass.value.is_integer():
            filename = f'lum_spec_{int(progenitor_mass.value):2d}M_r10000_dat.h5'
        else:
            filename = f'lum_spec_{progenitor_mass.value:.2f}M_r10000_dat.h5'

        metadata = {'Progenitor mass': progenitor_mass}

        return loaders.Fornax_2021(filename, metadata)

    # Populate Docstring with abbreviated param values
    __new__.__doc__ = __new__.__doc__.format(**_param_abbrv)

class SNOwGLoBES:
    """A model that does not inherit from SupernovaModel (yet) and imports a group of SNOwGLoBES files."""

    def __init__(self, tarfilename):
        """
        Parameters
        ----------
        tarfilename: str
            Absolute or relative path to tar archive with SNOwGLoBES files.
        """
        self.tfname = tarfilename
        tf = tarfile.open(self.tfname)

        # For now just pull out the "NoOsc" files.
        datafiles = sorted([f.name for f in tf if '.dat' in f.name])
        noosc = [df for df in datafiles if 'NoOsc' in df]
        noosc.sort(key=len)

        # Loop through the noosc files and pull out the number fluxes.
        self.time = []
        self.energy = None
        self.flux = {}
        self.fmin = 1e99
        self.fmax = -1e99

        for nooscfile in noosc:
            with tf.extractfile(nooscfile) as f:
                logging.debug('Reading {}'.format(nooscfile))
                meta = f.readline()
                metatext = meta.decode('utf-8')
                t = float(metatext.split('TBinMid=')[-1].split('sec')[0])
                dt = float(metatext.split('tBinWidth=')[-1].split('s')[0])
                dE = float(metatext.split('eBinWidth=')[-1].split('MeV')[0])

                data = Table.read(f, format='ascii.commented_header', header_start=-1)
                data.meta['t'] = t
                data.meta['dt'] = dt
                data.meta['dE'] = dE

                self.time.append(t)
                if self.energy is None:
                    self.energy = (data['E(GeV)'].data*1000).tolist()

            for flavor in ['NuE', 'NuMu', 'NuTau', 'aNuE', 'aNuMu', 'aNuTau']:
                if flavor in self.flux:
                    self.flux[flavor].append(data[flavor].data.tolist())
                else:
                    self.flux[flavor] = [data[flavor].data.tolist()]

        # We now have a table with rows=times and columns=energies. Transpose
        # so that rows=energy and cols=time.
        for k, v in self.flux.items():
            self.flux[k] = np.transpose(self.flux[k])
            self.fmin = np.minimum(self.fmin, np.min(self.flux[k]))
            self.fmax = np.maximum(self.fmax, np.max(self.flux[k]))

    def get_fluence(self, t):
        """Return the fluence at a given time t.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        fluence : dict
            A dictionary giving fluence at time t, keyed by flavor.
        """
        # Get index of closest element in the array
        idx = np.abs(np.asarray(self.time) - t).argmin()

        fluence = {}
        for k, fl in self.flux.items():
            fluence[k] = fl[:,idx]

        return fluence
