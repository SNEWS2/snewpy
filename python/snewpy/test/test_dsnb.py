"""
Tests for snewpy.dsnb

The DSNB only uses the public interface of supernova models, so the tests
use a small non-pinched toy model (Gaussian spectrum) and need no downloads.

Run with:
    pytest python/snewpy/test/test_dsnb.py -v
"""

import numpy as np
import pytest
from astropy import units as u

from snewpy.dsnb import CoreCollapseRate, DSNB, salpeter_imf
from snewpy.flavor import ThreeFlavor
from snewpy.flavor_transformation import AdiabaticMSW, NoTransformation
from snewpy.models.base import SupernovaModel
from snewpy.neutrino import MixingParameters


class GaussianModel(SupernovaModel):
    """Toy model with a Gaussian (non-pinched) NU_E_BAR spectrum, 1 s long."""

    def __init__(self, n_nu=1e57, mean=15.0, sigma=4.0, mass=15.0, with_mass=True):
        meta = {'Progenitor mass': mass * u.solMass} if with_mass else {}
        super().__init__(np.linspace(0, 1, 11) << u.s, meta)
        self.n_nu, self.mean, self.sigma = n_nu, mean, sigma

    def _get_initial_spectra_dict(self, t, E, flavors=ThreeFlavor):
        t = u.Quantity(t, ndmin=1)
        E = u.Quantity(E, ndmin=1).to_value(u.MeV)
        g = np.exp(-0.5 * ((E - self.mean) / self.sigma) ** 2)
        g = g / np.trapezoid(g, E)                       # unit area, 1/MeV
        spec = np.tile(g, (t.size, 1)) * self.n_nu       # 1/(MeV s)
        unit = 1 / (u.MeV * u.s)
        return {f: (spec if f == ThreeFlavor.NU_E_BAR else 0 * spec) << unit
                for f in flavors}


class FlatModel(SupernovaModel):
    """Toy model: dN/dE = k (constant) for all E, emitted over 1 s."""

    def __init__(self, k=1e57):
        super().__init__(np.linspace(0, 1, 5) << u.s, {'Progenitor mass': 15 * u.solMass})
        self.k = k

    def _get_initial_spectra_dict(self, t, E, flavors=ThreeFlavor):
        t = u.Quantity(t, ndmin=1)
        E = u.Quantity(E, ndmin=1)
        spec = np.full((t.size, E.size), self.k) << (1 / (u.MeV * u.s))
        return {f: (spec if f == ThreeFlavor.NU_E_BAR else 0 * spec) for f in flavors}


class ConstantRate:
    """R_CC(z) = r0, independent of z."""

    def __init__(self, r0=1e-4 * u.yr**-1 * u.Mpc**-3):
        self.r0 = r0

    def __call__(self, z):
        return np.ones_like(np.asarray(z, dtype=float)) * self.r0


T = np.array([0.0, 1.0]) << u.yr
E = np.linspace(2, 50, 120) << u.MeV
I_NUEBAR = list(ThreeFlavor).index(ThreeFlavor.NU_E_BAR)


def nuebar(flux):
    return flux.array[I_NUEBAR, 0, :].to_value(1 / (u.cm**2 * u.s * u.MeV))


# ---------------------------------------------------------------------------
class TestCoreCollapseRate:

    def test_present_day_rate(self):
        r0 = 1e-4 * u.yr**-1 * u.Mpc**-3
        assert np.isclose(CoreCollapseRate(r0=r0)(0).value, r0.value)

    def test_rate_peaks_at_intermediate_redshift(self):
        z = np.linspace(0, 5, 500)
        assert 1.0 < z[np.argmax(CoreCollapseRate()(z).value)] < 3.0

    def test_rate_rises_by_factor_8_to_z1(self):
        """Cole et al. form with multiplicative h: R(1)/R(0) = 8.6."""
        rcc = CoreCollapseRate()
        assert np.isclose(rcc(1.0).value / rcc(0.0).value, 8.63, rtol=0.01)

    def test_rate_positive(self):
        assert np.all(CoreCollapseRate()(np.linspace(0, 10, 200)).value > 0)

    def test_scales_with_r0(self):
        z = np.array([0.0, 0.5, 1.0, 2.0])
        r1 = CoreCollapseRate(1e-4 * u.yr**-1 * u.Mpc**-3)(z).value
        r2 = CoreCollapseRate(2e-4 * u.yr**-1 * u.Mpc**-3)(z).value
        assert np.allclose(r2 / r1, 2.0)

    def test_wrong_units_raise(self):
        with pytest.raises(u.UnitsError):
            CoreCollapseRate(r0=1e-4 * u.s**-1)


class TestIMF:

    def test_salpeter_decreasing(self):
        w = salpeter_imf(np.array([10.0, 20.0, 40.0]))
        assert np.all(np.diff(w) < 0)

    def test_weights_normalised_and_ordered(self):
        w = DSNB._imf_weights(np.array([9.0, 15.0, 25.0, 40.0]))
        assert np.isclose(w.sum(), 1.0)
        assert np.all(np.diff(w) < 0)


# ---------------------------------------------------------------------------
class TestDSNB:

    @pytest.fixture
    def dsnb(self):
        return DSNB([GaussianModel()], n_z=300)

    def test_works_for_non_pinched_model(self, dsnb):
        assert nuebar(dsnb.get_flux(T, E)).max() > 0

    def test_container_shape_and_unit(self, dsnb):
        flux = dsnb.get_flux(T, E)
        assert flux.array.shape == (len(ThreeFlavor), len(T), len(E))
        assert flux.array.unit.is_equivalent(u.cm**-2 * u.s**-1 * u.MeV**-1)

    def test_only_nuebar_is_filled(self, dsnb):
        arr = dsnb.get_flux(T, E).array.value
        others = [i for i in range(arr.shape[0]) if i != I_NUEBAR]
        assert np.all(arr[others] == 0)

    def test_flux_is_steady(self, dsnb):
        t = np.linspace(0, 1, 5) << u.yr
        arr = dsnb.get_flux(t, E).array[I_NUEBAR].value
        assert np.allclose(arr, arr[0])

    def test_fluence_is_flux_times_exposure(self, dsnb):
        flux = dsnb.get_flux(T, E)
        fluence = flux.integrate('time')
        expected = nuebar(flux) * (1 * u.yr).to_value(u.s)
        got = fluence.array[I_NUEBAR, 0, :].to_value(1 / (u.cm**2 * u.MeV))
        assert np.allclose(got, expected, rtol=1e-6)

    def test_flux_non_negative(self, dsnb):
        assert np.all(nuebar(dsnb.get_flux(T, E)) >= 0)

    def test_linear_in_number_of_neutrinos(self):
        f1 = nuebar(DSNB([GaussianModel(n_nu=1e57)], n_z=200).get_flux(T, E))
        f2 = nuebar(DSNB([GaussianModel(n_nu=3e57)], n_z=200).get_flux(T, E))
        assert np.allclose(f2, 3 * f1)

    def test_linear_in_r0(self):
        m = [GaussianModel()]
        f1 = nuebar(DSNB(m, rcc=CoreCollapseRate(1e-4 * u.yr**-1 * u.Mpc**-3), n_z=200).get_flux(T, E))
        f2 = nuebar(DSNB(m, rcc=CoreCollapseRate(2e-4 * u.yr**-1 * u.Mpc**-3), n_z=200).get_flux(T, E))
        assert np.allclose(f2, 2 * f1)

    def test_redshift_weight_is_one_over_H(self):
        """
        For constant R_CC and flat dN/dE = k the flux is analytic,
            c/H0 * R * k * integral dz / E(z),
        i.e. there is no extra factor (1+z) in the integrand.
        """
        from astropy import constants as const
        from scipy.integrate import quad
        k, r0 = 1e57, 1e-4 * u.yr**-1 * u.Mpc**-3
        d = DSNB([FlatModel(k)], rcc=ConstantRate(r0), z_max=5.0, n_z=4000)
        got = nuebar(d.get_flux(T, np.array([10.0, 20.0]) << u.MeV))
        z_lo = 1e-4                       # first point of the redshift grid
        integral = quad(lambda z: 1.0 / d.cosmology.efunc(z), z_lo, 5.0)[0]
        expected = (const.c / d.cosmology.H(0) * r0 * k / u.MeV * integral)
        expected = expected.to_value(1 / (u.cm**2 * u.s * u.MeV))
        assert np.allclose(got, expected, rtol=2e-3)

    def test_matches_first_principles_geometry(self):
        """
        Independent derivation without the c/H0/E(z) formula: supernovae per
        observed time per dz = R/(1+z) * 4 pi dVc/dz (time dilation), number
        fluence per supernova = dN/dE'(E') (1+z) / (4 pi d_C^2).  The product,
        integrated over z, must equal DSNB.get_flux.
        """
        d = DSNB([GaussianModel()], z_max=5.0, n_z=1500)
        Eo = np.array([10.0, 20.0, 30.0]) << u.MeV
        got = nuebar(d.get_flux(T, Eo))
        z = np.linspace(1e-4, 5.0, 20001)
        cosmo = d.cosmology
        dC, dVc = cosmo.comoving_distance(z), cosmo.differential_comoving_volume(z)
        Rz = d.rcc(z)
        Eg = np.geomspace(0.1, 300, 600) << u.MeV
        dN = d._dNdE_per_sn(Eg)
        for k, En in enumerate(Eo):
            dNe = np.interp(En.value * (1 + z), Eg.value, dN, left=0, right=0) / u.MeV
            integrand = (Rz / (1 + z) * 4 * np.pi * u.sr * dVc) * (dNe * (1 + z) / (4 * np.pi * dC**2))
            expected = np.trapezoid(integrand.to_value(1 / (u.s * u.cm**2 * u.MeV)), z)
            assert np.isclose(got[k], expected, rtol=3e-3)

    def test_flavor_transformation_is_applied(self, dsnb):
        none = nuebar(dsnb.get_flux(T, E, NoTransformation()))
        nmo = nuebar(dsnb.get_flux(T, E, AdiabaticMSW(MixingParameters('NORMAL'))))
        assert not np.allclose(none, nmo)

    def test_default_is_no_transformation(self, dsnb):
        assert np.allclose(nuebar(dsnb.get_flux(T, E)),
                           nuebar(dsnb.get_flux(T, E, NoTransformation())))

    def test_high_energy_reaches_redshifted_spectrum(self):
        """Emission grid must extend to E_max * (1 + z_max)."""
        model = DSNB([GaussianModel(mean=150.0, sigma=20.0)], z_max=5.0, n_z=200)
        flux = nuebar(model.get_flux(T, np.array([40.0, 50.0]) << u.MeV))
        assert np.all(flux > 0)

    def test_single_energy(self, dsnb):
        assert dsnb.get_flux(T, 20 << u.MeV).array.shape == (len(ThreeFlavor), len(T), 1)

    def test_two_models_are_weighted(self):
        lo, hi = GaussianModel(mean=10.0, mass=10.0), GaussianModel(mean=25.0, mass=40.0)
        f_lo = nuebar(DSNB([lo], n_z=200).get_flux(T, E))
        f_hi = nuebar(DSNB([hi], n_z=200).get_flux(T, E))
        both = DSNB([lo, hi], n_z=200)
        w = both._weights
        assert np.isclose(w.sum(), 1.0)
        assert np.allclose(nuebar(both.get_flux(T, E)), w[0] * f_lo + w[1] * f_hi, rtol=1e-6)

    def test_explicit_weights_allow_equal_masses(self):
        models = [GaussianModel(mean=10.0), GaussianModel(mean=25.0)]
        d = DSNB(models, weights=[1.0, 3.0], n_z=100)
        assert np.allclose(d._weights, [0.25, 0.75])

    def test_duplicate_masses_raise(self):
        with pytest.raises(ValueError):
            DSNB([GaussianModel(), GaussianModel()])

    def test_missing_mass_metadata_raises(self):
        with pytest.raises(ValueError):
            DSNB([GaussianModel(with_mass=False)])

    def test_mass_tuple_is_accepted(self):
        d = DSNB([(GaussianModel(with_mass=False), 20.0)], n_z=100)
        assert np.isclose(d._masses[0], 20.0)

    def test_bad_weights_raise(self):
        with pytest.raises(ValueError):
            DSNB([GaussianModel()], weights=[1.0, 2.0])
