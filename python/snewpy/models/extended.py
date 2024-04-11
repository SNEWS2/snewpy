import itertools as it
import os
from abc import ABC, abstractmethod
from warnings import warn

import numpy as np
from astropy import units as u
from astropy.table import Table, join
from astropy.units import UnitTypeError, get_physical_type
from astropy.units.quantity import Quantity
from scipy.special import loggamma
from snewpy import _model_downloader

from snewpy.neutrino import Flavor
from snewpy.flavor_transformation import NoTransformation
from functools import wraps

from snewpy.flux import Flux
from snewpy.models.base import SupernovaModel


class ExtendedModel(SupernovaModel):
    """Class defining a supernova model with a cooling tail extension."""

    def __init__(self, *args):
        """Initialize extended supernova model class."""
        if isinstance(args[0],SupernovaModel):
            self.__dict__ = args[0].__dict__.copy()
            for method_name in dir(args[0]):
                if callable(getattr(args[0], method_name)) and method_name[0] != '_':
                    if method_name == 'get_initial_spectra':
                        self._get_initial_spectra = getattr(args[0], method_name)
                    else:
                        self.method_name  = getattr(args[0], method_name)
            self.t_final = self.time[-1]
            self.L_final = {Flavor.NU_E: self.luminosity[Flavor.NU_E][-1],
                            Flavor.NU_X: self.luminosity[Flavor.NU_X][-1],
                            Flavor.NU_E_BAR: self.luminosity[Flavor.NU_E_BAR][-1],
                            Flavor.NU_X_BAR: self.luminosity[Flavor.NU_X_BAR][-1]}
        else:
            raise TypeError("ExtendedModel.__init__ requires a SupernovaModel object")

    def get_initial_spectra(self, *args):
        """Get neutrino spectra/luminosity curves before oscillation"""
        return self._get_initial_spectra(*args)

    def get_extended_luminosity(self, t, k = -1., A = 1e51 * u.erg / u.s, tau_c = 36. * u.s, alpha = 2.66, flavor = Flavor.NU_E):
        """Get neutrino luminosity from supernova cooling tail luminosity model.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate luminosity.
        k : float
            Power law factor (default: -1)
        A : astropy.Quantity
            Normalization factor (default: 1e51 erg/s)
        tau_c : astropy.Quantity
            Exponential decay characteristic timescale (default: 36 s)
        alpha : float
            Exponential decay factor (default: 2.66)

        Returns
        -------
        astropy.Quantity
            Luminosity calculated from cooling tail model.
        """
        if t.value < 0.5:
            warn("Extended luminosity model not applicable to early times")
        if A is None:
            tf = self.t_final
            Lf = self.L_final[flavor]
            A = Lf/((tf.value**k) * np.exp(-1.*((tf/tau_c)**alpha)) * u.erg / u.s)
        return A * (t.value**k) * np.exp(-1.*((t/tau_c)**alpha)) * u.erg / u.s

    def extend(self, ts, k = -1., A = None, tau_c = 36. * u.s, alpha = 2.66):
        for t in ts:
            self.time = np.append(self.time, t)
            self.meanE = np.append(self.meanE, self.meanE[-1])
            for flavor in [Flavor.NU_E, Flavor.NU_X, Flavor.NU_E_BAR, Flavor.NU_X_BAR]:
                L_ext = self.get_extended_luminosity(t, k = k, A = A, tau_c = tau_c, alpha = alpha, flavor = flavor)
                self.luminosity[flavor] = np.append(self.luminosity[flavor], L_ext)
