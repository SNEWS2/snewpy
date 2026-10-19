"""
Diffuse Supernova Neutrino Background (DSNB) model for SNEWPY.

The DSNB is built from ordinary SNEWPY supernova models and exposes the same
interface as they do: ``get_flux(t, E, flavor_xform)`` returns a
:class:`snewpy.flux.Container`, ``flux.integrate('time')`` gives the fluence
over a detector exposure, and the fluence can be passed to ``RateCalculator``.

Only public methods of the supernova models are used
(``get_transformed_spectra`` and the ``time`` attribute), so any model class
works, pinched or not.

Physics follows Li, Vagins & Wurm (2022) [arXiv:2201.12920].
Core-collapse rate: Hopkins & Beacom (2006), ApJ 651, 142.
Cosmology: Planck 2018.

Example
-------
>>> import numpy as np
>>> import astropy.units as u
>>> from snewpy.models.ccsn import Nakazato_2013
>>> from snewpy.dsnb import DSNB
>>> from snewpy.neutrino import MixingParameters
>>> from snewpy.flavor_transformation import AdiabaticMSW
>>> from snewpy.rate_calculator import RateCalculator
>>>
>>> sn_models = [Nakazato_2013(progenitor_mass=m*u.Msun, revival_time=100*u.ms,
...                            metallicity=0.02, eos='shen') for m in (13, 20)]
>>> model    = DSNB(sn_models)
>>> xform    = AdiabaticMSW(MixingParameters('NORMAL'))
>>> times    = np.linspace(0, 1, 2) << u.yr
>>> energies = np.linspace(0, 50, 501) << u.MeV
>>> flux     = model.get_flux(t=times, E=energies, flavor_xform=xform)
>>> fluence  = flux.integrate('time')
>>> events   = RateCalculator().run(fluence, "wc100kt30prct")
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from snewpy.flavor import ThreeFlavor
from snewpy.flux import Container
from snewpy.flavor_transformation import FlavorTransformation, NoTransformation
from typing import Optional, Sequence

__all__ = ["CoreCollapseRate", "salpeter_imf", "DSNB"]

_MPC_TO_CM  = 3.0856775815e24
_YR_TO_S    = 3.15576e7
_PLANCK18   = FlatLambdaCDM(H0=67.4, Om0=0.315)


class CoreCollapseRate:
    """
    Redshift-dependent core-collapse supernova rate.

    Implements the parametrisation of Hopkins & Beacom (2006), ApJ 651, 142,
    which uses the functional form of Cole et al. (2001),

        rho(z) = (a + b z) h / (1 + (z/c)^d),

    with a=0.0170, b=0.13, c=3.3, d=5.3 and h=0.7.  The overall factor h
    cancels here, because the rate is normalised to ``r0`` at z=0.

    Parameters
    ----------
    r0 : astropy.units.Quantity
        Present-day rate.  Default: 1e-4 yr^{-1} Mpc^{-3}.

    Examples
    --------
    >>> import astropy.units as u
    >>> rcc = CoreCollapseRate()
    >>> rcc(0)
    <Quantity 1.e-4 1 / (Mpc3 yr)>
    >>> rcc(1) > rcc(0)
    True
    """

    _a, _b, _c, _d, _h = 0.0170, 0.13, 3.3, 5.3, 0.7

    def __init__(self, r0: u.Quantity = 1e-4 * u.yr**-1 * u.Mpc**-3):
        if not r0.unit.is_equivalent(u.yr**-1 * u.Mpc**-3):
            raise u.UnitsError(
                f"r0 must have units yr^-1 Mpc^-3, got {r0.unit}"
            )
        self.r0 = r0.to(u.yr**-1 * u.Mpc**-3)

    def __call__(self, z) -> u.Quantity:
        """
        Evaluate R_CC(z).

        Parameters
        ----------
        z : float or array_like
            Redshift.

        Returns
        -------
        astropy.units.Quantity
            Rate in yr^{-1} Mpc^{-3}.
        """
        z   = np.asarray(z, dtype=float)
        num = (self._a + self._b * z) * self._h
        den = self._a * self._h * (1.0 + (z / self._c)**self._d)
        return self.r0 * num / den


def salpeter_imf(masses: np.ndarray) -> np.ndarray:
    """
    Salpeter (1955) initial mass function weights.

    Returns weights proportional to M^{-2.35}.

    Parameters
    ----------
    masses : array_like
        Progenitor masses in solar masses.

    Returns
    -------
    numpy.ndarray
        Unnormalised weights.
    """
    return np.asarray(masses, dtype=float) ** (-2.35)


class DSNB:
    """
    Diffuse Supernova Neutrino Background model.

    The observed flux is the redshift integral of the time-integrated
    emission spectrum of a single supernova, weighted by the cosmic
    core-collapse rate (Eq. 1 of Li, Vagins & Wurm 2022):

        dPhi/dE = c/H0 * integral_0^z_max R_CC(z) / E(z)
                         * dN/dE'[(1+z) E]  dz

    where dN/dE' is the IMF-weighted, time-integrated, flavor-transformed
    NU_E_BAR spectrum emitted by one supernova.  That spectrum is taken from
    ``model.get_transformed_spectra(t, E, flavor_xform)`` followed by
    ``.integrate('time')``, for every supernova model in the list.

    Parameters
    ----------
    sn_models : list
        Supernova model instances (any SNEWPY ``SupernovaModel``).  The
        progenitor mass used for the IMF weight is read from
        ``model.metadata['Progenitor mass']``.  A ``(model, mass)`` tuple
        can be given instead for models without that metadata entry
        (``mass`` in solar masses, or a Quantity).
    weights : array_like, optional
        Explicit weight of each model (normalised internally).  When given,
        the IMF weighting is skipped, and models with equal progenitor
        masses are allowed.
    rcc : CoreCollapseRate, optional
        Default: Hopkins & Beacom (2006), R_CC(0) = 1e-4 yr^-1 Mpc^-3.
    cosmology : astropy cosmology instance, optional
        Default: Planck 2018 (H0=67.4, Om0=0.315).
    z_max : float, optional
        Upper redshift limit.  Default: 5.0.
    n_z : int, optional
        Number of redshift grid points.  Default: 800.
    """

    #: lowest emitted energy, MeV, and number of points of the emission grid
    _E_EMIT_MIN = 0.1
    _N_E_EMIT   = 600

    def __init__(
        self,
        sn_models : Sequence,
        weights   : Optional[Sequence[float]] = None,
        rcc       : Optional[CoreCollapseRate] = None,
        cosmology                              = None,
        z_max     : float                      = 5.0,
        n_z       : int                        = 800,
    ):
        models, masses = [], []
        for item in sn_models:
            if isinstance(item, tuple):
                model, mass = item
            else:
                model, mass = item, None
            models.append(model)
            masses.append(self._progenitor_mass(model, mass, need=weights is None))

        if weights is not None:
            w = np.asarray(weights, dtype=float)
            if w.shape != (len(models),) or np.any(w < 0) or w.sum() <= 0:
                raise ValueError(
                    "weights must be non-negative, with one entry per model."
                )
            self._models  = models
            self._masses  = np.array([np.nan if m is None else m for m in masses])
            self._weights = w / w.sum()
        else:
            if len(set(masses)) != len(masses):
                raise ValueError(
                    "Duplicate progenitor masses in sn_models. Give explicit "
                    "`weights`, or keep one model per mass."
                )
            order         = np.argsort(masses)
            self._models  = [models[i] for i in order]
            self._masses  = np.array([masses[i] for i in order], dtype=float)
            self._weights = self._imf_weights(self._masses)

        self.rcc       = rcc or CoreCollapseRate()
        self.cosmology = cosmology or _PLANCK18
        self.z_max     = float(z_max)
        self.n_z       = int(n_z)

    # ------------------------------------------------------------------
    @staticmethod
    def _progenitor_mass(model, mass, need: bool):
        """Progenitor mass in solar masses (float), or None if not needed."""
        if mass is None:
            try:
                mass = model.metadata['Progenitor mass']
            except (AttributeError, KeyError, TypeError):
                if need:
                    raise ValueError(
                        f"{type(model).__name__} has no 'Progenitor mass' in its "
                        "metadata. Pass (model, mass) tuples or explicit `weights`."
                    )
                return None
        return float(u.Quantity(mass, u.solMass).to_value(u.solMass))

    @staticmethod
    def _imf_weights(sorted_masses: np.ndarray) -> np.ndarray:
        """
        Salpeter IMF integral weights for unevenly spaced progenitor masses.

        Weight for progenitor i = integral of M^{-2.35} dM over its mass
        interval, where intervals are defined by midpoints between adjacent
        masses (8 and 100 Msun as outer limits).
        """
        n          = len(sorted_masses)
        edges      = np.empty(n + 1)
        edges[0]   = 8.0
        edges[-1]  = 100.0
        for i in range(1, n):
            edges[i] = 0.5 * (sorted_masses[i - 1] + sorted_masses[i])
        weights = np.array([
            quad(lambda m: m**(-2.35), edges[i], edges[i + 1])[0]
            for i in range(n)
        ])
        return weights / weights.sum()

    # ------------------------------------------------------------------
    def _dNdE_per_sn(
        self,
        E_emit       : u.Quantity,
        flavor_xform : FlavorTransformation = NoTransformation(),
    ) -> np.ndarray:
        """
        Weighted, time-integrated NU_E_BAR spectrum of one core collapse.

        Uses only the public interface of the supernova models:
        ``get_transformed_spectra`` gives d2N/dE/dt after the flavor
        transformation (any time and energy dependence of the transformation
        is handled by the model), and ``integrate('time')`` over the model's
        own time grid gives dN/dE per supernova.

        Parameters
        ----------
        E_emit : astropy.units.Quantity
            Emitted-energy grid.
        flavor_xform : FlavorTransformation
            Transformation inside the supernova.  Default:
            ``NoTransformation()``.

        Returns
        -------
        numpy.ndarray, shape (len(E_emit),)
            dN/dE in MeV^-1.
        """
        dNdE = np.zeros(E_emit.size)
        for model, weight in zip(self._models, self._weights):
            spectra  = model.get_transformed_spectra(model.time, E_emit, flavor_xform)
            per_sn   = spectra.integrate('time')          # (flavor, 1, energy), 1/MeV
            i_nuebar = list(per_sn.flavor).index(ThreeFlavor.NU_E_BAR)
            dNdE    += weight * per_sn.array[i_nuebar, 0, :].to_value(1 / u.MeV)
        return dNdE

    # ------------------------------------------------------------------
    def get_flux(
        self,
        t            : u.Quantity,
        E            : u.Quantity,
        flavor_xform : FlavorTransformation = NoTransformation(),
    ) -> Container:
        """
        DSNB flux, same call signature as the supernova models' get_flux
        (without the distance).

        The DSNB is steady-state: the returned flux has the same value at
        every time in ``t``.  ``flux.integrate('time')`` gives the fluence
        over the exposure ``[t[0], t[-1]]``.

        Parameters
        ----------
        t : astropy.units.Quantity
            Time sample points of the exposure window, e.g. ``[0, 1] * u.yr``.
            At least 2 points are needed to integrate over time.
        E : astropy.units.Quantity
            Observed neutrino energies.
        flavor_xform : FlavorTransformation
            Flavor transformation inside the source supernovae, e.g.
            ``AdiabaticMSW(MixingParameters('NORMAL'))``.  Default:
            ``NoTransformation()``, as for the supernova models.

        Returns
        -------
        snewpy.flux.Container
            Flux in cm^-2 s^-1 MeV^-1, shape (6, len(t), len(E)).  Only the
            NU_E_BAR component is non-zero.
        """
        E_obs = np.atleast_1d(E.to_value(u.MeV))
        z     = np.linspace(1e-4, self.z_max, self.n_z)

        # Emission grid has to reach E_max * (1 + z_max)
        E_hi   = max(300.0, E_obs.max() * (1.0 + self.z_max))
        E_grid = np.geomspace(self._E_EMIT_MIN, E_hi, self._N_E_EMIT) << u.MeV
        dNdE   = self._dNdE_per_sn(E_grid, flavor_xform)

        # Emitted energy at each (observed energy, redshift): (N_E, N_z)
        E_emit    = np.outer(E_obs, 1.0 + z)
        dNdE_grid = np.interp(
            E_emit.ravel(), E_grid.value, dNdE, left=0.0, right=0.0
        ).reshape(E_emit.shape)

        Hz  = self.cosmology.efunc(z)
        rcc = self.rcc(z).to_value(u.yr**-1 * u.Mpc**-3)

        integrand = dNdE_grid * (rcc / Hz)[np.newaxis, :]
        integral  = np.trapezoid(integrand, z, axis=1)

        H0_s      = self.cosmology.H(0).to_value(u.s**-1)
        c_cm_s    = const.c.to_value(u.cm / u.s)
        rate_conv = 1.0 / (_YR_TO_S * _MPC_TO_CM**3)
        phi       = (c_cm_s / H0_s) * integral * rate_conv      # cm^-2 s^-1 MeV^-1

        t       = u.Quantity(t, ndmin=1)
        flavors = list(ThreeFlavor)
        data    = np.zeros((len(flavors), t.size, E_obs.size))
        data[flavors.index(ThreeFlavor.NU_E_BAR), :, :] = phi
        return Container(
            data * u.cm**-2 * u.s**-1 * u.MeV**-1,
            flavor        = flavors,
            time          = t.to(u.s),
            energy        = u.Quantity(E, ndmin=1).to(u.MeV),
            flavor_scheme = ThreeFlavor,
        )
