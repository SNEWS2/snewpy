# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.ccsn`` contains models of core-collapse supernovae
derived from the :class:`SupernovaModel` base class.

You can :ref:`download neutrino fluxes for each of these models <sec-download_models>`
using ``snewpy.get_models("<model_name>")``.

.. _Garching Supernova Archive: https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/
"""

import logging
import os
import sys
import tarfile

import h5py
import numpy as np
from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table
from scipy.special import gamma, lpmv
import re

try:
    import healpy as hp
except ImportError:
    pass

from snewpy import model_path
from snewpy.neutrino import Flavor
from snewpy.models import loaders
from .base import PinchedModel, SupernovaModel, _GarchingArchiveModel
from .registry import check_valid_params, get_param_combinations

class _RegistryModel():
    """TODO: empty base class for now?"""
    pass

class Analytic3Species(PinchedModel):
    """An analytical model calculating spectra given total luminosity,
    average energy, and rms or pinch, for each species.
    """

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
    # Todo: Add to docs, parameters names and allowed values (for readthedocs) See if add-to-docs it can be done auto
    param = {'progenitor_mass': [13, 20, 30, 50] * u.Msun,
             'revival_time': [0, 100, 200, 300] * u.ms,
             'metallicity': [0.02, 0.004],
             'eos': ['LS220', 'shen', 'togashi']}

    def __new__(cls, *, progenitor_mass=None, revival_time=None, metallicity=None, eos=None):
        """Model initialization.

        Parameters
        ----------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun
        revival_time: astropy.units.Quantity
            Time of shock revival in model in units ms
            Selecting 0 ms will load a black hole formation model
        metallicity: float
            Progenitor metallicity
        eos: str
            Equation of state

        Raises
        ------
        ValueError
            If a combination of parameters is invalid when loading form parameters

        See also
        --------
        snewpy._model_registry : Describes allowed values for parameters `progenitor_mass`, `revival_time`,
                                 `metallicity`, and `eos`
        Examples
        --------
        >>> from snewpy.models.ccsn import Nakazato_2013; import astropy.units as u
        >>> Nakazato_2013(progenitor_mass=13*u.Msun, metallicity=0.004, revival_time=0*u.s, eos='togashi')
        Nakazato_2013 Model: nakazato-togashi-BH-z0.004-s30.0.fits
        Progenitor mass  : 30.0 solMass
        EOS              : Togashi
        Metallicity      : 0.004
        Revival time     : 0.0 ms
        """
        # Store model metadata.
        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos,
            'Metallicity': metallicity,
            'Revival time': revival_time
        }

        # TODO: Check GitHub PR for error in this example
        # Attempt to load model from parameters

        # Build user params, check validity, construct filename, then load from filename
        user_params = dict(zip(cls.param.keys(), (progenitor_mass, revival_time, metallicity, eos)))
        check_valid_params(cls, **user_params)

        # Strip units for filename construction
        progenitor_mass = progenitor_mass.to(u.Msun).value
        revival_time = revival_time.to(u.ms).value

        if revival_time != 0:
            fname = f"nakazato-{eos}-z{metallicity}-t_rev{int(revival_time)}ms-s{progenitor_mass:3.1f}.fits"
        else:
            fname = f"nakazato-{eos}-BH-z{metallicity}-s{progenitor_mass:3.1f}.fits"

        filename = os.path.join(model_path, cls.__name__, fname)

        if not os.path.isfile(filename):
            # download file from GitHub/Zenodo
            raise NotImplementedError()

        return loaders.Nakazato_2013(filename, metadata)


    @classmethod
    def get_param_combinations(cls, *, progenitor_mass=None, revival_time=None, metallicity=None, eos=None):
        user_param = dict(zip(cls.param.keys(), (progenitor_mass, revival_time, metallicity, eos)))
        return get_param_combinations(cls, **user_param)

    @staticmethod
    def isvalid_param_combo(progenitor_mass, revival_time, metallicity, eos):
        """Returns True if the parameter combination is valid (Corresponds to an existing Nakazato_2013 model)
        See __init__ for a full parameter descriptions.
        """
        if revival_time.to(u.s).value == 0:
            return progenitor_mass.to(u.Msun).value == 30 and metallicity == 0.004
        else:
            return eos in ['shen'] and not (progenitor_mass.to(u.Msun).value == 30 and metallicity == 0.004)


class Sukhbold_2015(_RegistryModel):
    """Model based on simulations from Sukhbold et al., ApJ 821:38,2016. Models were shared privately by email.
    """
    pass

class Tamborra_2014(_RegistryModel):
    """Model based on 3D simulations from `Tamborra et al., PRD 90:045032, 2014 <https://arxiv.org/abs/1406.0006>`_.
    Data files are from the `Garching Supernova Archive`_.
    """
    pass

class Bollig_2016(_RegistryModel):
    """Model based on simulations from `Bollig et al. (2016) <https://arxiv.org/abs/1508.00785>`_. Models were taken, with permission, from the Garching Supernova Archive.
    """
    pass

class Walk_2018(_RegistryModel):
    """Model based on SASI-dominated simulations from `Walk et al.,
    PRD 98:123001, 2018 <https://arxiv.org/abs/1807.02366>`_. Data files are from
    the `Garching Supernova Archive`_.
    """
    pass

class Walk_2019(_RegistryModel):
    """Model based on SASI-dominated simulations from `Walk et al.,
    PRD 101:123013, 2019 <https://arxiv.org/abs/1910.12971>`_. Data files are
    from the `Garching Supernova Archive`_.
    """
    pass


class OConnor_2013(PinchedModel): # TODO: Requires changes to the model file to have one file per model instead of a single gzip archive!
    """Model based on the black hole formation simulation in `O'Connor & Ott (2013) <https://arxiv.org/abs/1207.1100>`_.
    """
    def __init__(self, base, mass=15, eos='LS220'):
        """
        Parameters
        ----------
        base : str
            Path of directory containing model files
        mass : int
            Progenitor mass
        eos : string
            Equation of state used in simulation
        """

        # Open luminosity file.
        tf = tarfile.open(base+'{}_timeseries.tar.gz'.format(eos))

        # Extract luminosity data.
        dataname = 's{:d}_{}_timeseries.dat'.format(mass, eos)
        datafile = tf.extractfile(dataname)
        simtab = ascii.read(datafile, names=['TIME', 'L_NU_E', 'L_NU_E_BAR', 'L_NU_X',
                                               'E_NU_E', 'E_NU_E_BAR', 'E_NU_X',
                                               'RMS_NU_E', 'RMS_NU_E_BAR', 'RMS_NU_X'])

        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2)/(simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2)/(simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2)/(simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        #note, here L_NU_X is already divided by 4

        self.filename = datafile
        self.EOS = eos
        self.progenitor_mass = mass * u.Msun

        metadata = {
            'Progenitor mass':self.progenitor_mass,
            'EOS':self.EOS,
        }
        super().__init__(simtab, metadata)


class OConnor_2015(_RegistryModel):
    """Model based on the black hole formation simulation in `O'Connor (2015) <https://arxiv.org/abs/1411.7058>`_.
    """
    def __new__(cls, progenitor_mass=None, eos='LS220'):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos,
        }

        filename = os.path.join(model_path, cls.__name__, 'M1_neutrinos.dat')

        return loaders.OConnor_2015(filename, metadata)


class Zha_2021(_RegistryModel):
    """Model based on the hadron-quark phse transition models from `Zha et al. 2021 <https://arxiv.org/abs/2103.02268>`_.
    """
    def __new__(cls, progenitor_mass=None, eos='LS220'):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """

        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos,
        }
        
        filename = os.path.join(model_path, cls.__name__, f's{progenitor_mass}.dat')

        return loaders.Zha_2021(filename, metadata)


class Warren_2020(_RegistryModel):
    """Model based on simulations from Warren et al., ApJ 898:139, 2020.
    Neutrino fluxes available at https://doi.org/10.5281/zenodo.3667908."""

    def __new__(cls, progenitor_mass=None, turbmixing_param=None, eos='SFHo'):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        # Set model metadata.
        metadata = {
            'Progenitor mass': progenitor_mass,
            'Turb. mixing param.': turbmixing_param,
            'EOS': eos,
        }

        filename = os.path.join(model_path, cls.__name__, f'stir_a{turbmixing_param}/stir_multimessenger_a{turbmixing_param}_m{progenitor_mass}.h5')

        return loaders.Warren_2020(filename, metadata)


class Kuroda_2020(_RegistryModel):
    """Model based on simulations from `Kuroda et al. (2020) <https://arxiv.org/abs/2009.07733>`_."""

    def __new__(cls, eos='LS220', progenitor_mass=20*u.Msun):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        metadata = {
            'Progenitor mass': progenitor_mass,
            'EOS': eos,
            }

        filename = os.path.join(model_path, cls.__name__, f'LnuR00B00.dat')  # TODO: replace hardcoded filename with one based on input parameters

        return loaders.Kuroda_2020(filename, metadata)


class Fornax_2019(_RegistryModel):
    """Model based on 3D simulations from D. Vartanyan, A. Burrows, D. Radice, M.  A. Skinner and J. Dolence, MNRAS 482(1):351, 2019. 
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.3d/
    """

    def __new__(cls, progenitor_mass=None, cache_flux=False):
        metadata = {
            'Progenitor mass': progenitor_mass,
            }

        filename = os.path.join(model_path, cls.__name__, f'lum_spec_{progenitor_mass}M.h5')

        return loaders.Fornax_2019(filename, metadata, cache_flux=cache_flux)


class Fornax_2021(_RegistryModel):
    """Model based on 3D simulations from D. Vartanyan, A. Burrows, D. Radice, M.  A. Skinner and J. Dolence, MNRAS 482(1):351, 2019. 
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.3d/
        """

    def __new__(cls, progenitor_mass=None):
        metadata = {
            'Progenitor mass': progenitor_mass,
            }

        filename = os.path.join(model_path, cls.__name__, f'lum_spec_{progenitor_mass}M_r10000_dat.h5')

        return loaders.Fornax_2021(filename, metadata)

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
