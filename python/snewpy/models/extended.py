from warnings import warn

import numpy as np
from astropy import units as u

from snewpy.neutrino import Flavor
from snewpy.models.base import SupernovaModel


class ExtendedModel(SupernovaModel):
    """Class defining a supernova model with a cooling tail extension."""

    def __init__(self, base_model):
        """Initialize extended supernova model class."""
        if not isinstance(base_model, SupernovaModel):
            raise TypeError("ExtendedModel.__init__ requires a SupernovaModel object")

        self.__dict__ = base_model.__dict__.copy()
        for method_name in dir(base_model):
            if callable(getattr(base_model, method_name)) and method_name[0] != '_':
                if method_name == 'get_initial_spectra':
                    self._get_initial_spectra = getattr(base_model, method_name)
                else:
                    setattr(self, method_name, getattr(base_model, method_name))
        self.t_final = self.time[-1]
        self.L_final = {flv: self.luminosity[flv][-1] for flv in Flavor}

    def get_initial_spectra(self, *args, **kwargs):
        """Get neutrino spectra/luminosity curves before oscillation"""
        return self._get_initial_spectra(*args, **kwargs)

    def get_extended_luminosity(self, t, k=-1., A=None, tau_c=36. * u.s, alpha=2.66, flavor = Flavor.NU_E):
        """Get neutrino luminosity from supernova cooling tail luminosity model.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate luminosity.
        k : float
            Power law factor (default: -1)
        A : astropy.Quantity
            Normalization factor (default: None, automatically match original model data)
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
            A = Lf / (tf.value**k * np.exp(-(tf/tau_c)**alpha))
        return A * t.value**k * np.exp(-(t/tau_c)**alpha)

    def extend(self, ts, k=-1., A=None, tau_c=36. * u.s, alpha=2.66):
        """Extend supernova model to specific times.

        Parameters
        ----------
        ts : astropy.Quantity
            Times to add to supernova model.
        k : float
            Power law factor (default: -1)
        A : astropy.Quantity
            Normalization factor (default: None, automatically match original model data)
        tau_c : astropy.Quantity
            Exponential decay characteristic timescale (default: 36 s)
        alpha : float
            Exponential decay factor (default: 2.66)
        """
        # Select times after the end of the model
        select = ts > self.t_final

        for t in ts[select]:
            self.time = np.append(self.time, t)
            for flavor in Flavor:
                L_ext = self.get_extended_luminosity(t, k = k, A = A, tau_c = tau_c, alpha = alpha, flavor = flavor)
                self.luminosity[flavor] = np.append(self.luminosity[flavor], L_ext)
                self.meanE[flavor] = np.append(self.meanE[flavor], self.meanE[flavor][-1])
                self.pinch[flavor] = np.append(self.pinch[flavor], self.pinch[flavor][-1])
