"""
Tests for snewpy_dsnb.dsnb

Run with:
    pytest snewpy_dsnb/tests/test_dsnb.py -v
"""

import numpy as np
import pytest
import astropy.units as u

from snewpy_dsnb.dsnb import CoreCollapseRate, salpeter_imf, DSNB
from snewpy.models.ccsn import Nakazato_2013
from snewpy.flavor_transformation import AdiabaticMSW, NoTransformation
from snewpy.neutrino import MassHierarchy, Flavor


# Shared fixtures
@pytest.fixture(scope="module")
def sn_models():
    """Two Nakazato progenitors, loaded once for the whole test session."""
    return [
        (Nakazato_2013(progenitor_mass=13*u.Msun, revival_time=100*u.ms,
                       metallicity=0.02, eos='shen'), 13),
        (Nakazato_2013(progenitor_mass=20*u.Msun, revival_time=100*u.ms,
                       metallicity=0.02, eos='shen'), 20),
    ]


@pytest.fixture(scope="module")
def model(sn_models):
    """DSNB instance with coarser z grid for speed."""
    return DSNB(sn_models, n_z=200)


@pytest.fixture(scope="module")
def E():
    return np.linspace(5, 50, 30) * u.MeV


@pytest.fixture(scope="module")
def t():
    return [0, 1] * u.yr

# CoreCollapseRate
class TestCoreCollapseRate:

    def test_present_day_rate(self):
        """R_CC(0) must equal the reference r0."""
        r0  = 1e-4 * u.yr**-1 * u.Mpc**-3
        rcc = CoreCollapseRate(r0=r0)
        assert np.isclose(
            rcc(0).to(u.yr**-1 * u.Mpc**-3).value, r0.value, rtol=1e-6
        )

    def test_rate_peaks_around_z1(self):
        """Star-formation history peaks near z ~ 1--2."""
        rcc    = CoreCollapseRate()
        z      = np.linspace(0, 5, 200)
        z_peak = z[np.argmax(rcc(z).value)]
        assert 1.0 < z_peak < 3.0, f"Peak at z={z_peak:.2f}, expected 1--3"

    def test_rate_positive_everywhere(self):
        """Rate must be positive for all z >= 0."""
        rcc = CoreCollapseRate()
        assert np.all(rcc(np.linspace(0, 10, 500)).value > 0)

    def test_custom_r0_scales_linearly(self):
        """Custom r0 scales the entire rate proportionally."""
        rcc1  = CoreCollapseRate(r0=1e-4 * u.yr**-1 * u.Mpc**-3)
        rcc2  = CoreCollapseRate(r0=2e-4 * u.yr**-1 * u.Mpc**-3)
        z     = np.array([0.0, 0.5, 1.0, 2.0])
        ratio = rcc2(z).value / rcc1(z).value
        assert np.allclose(ratio, 2.0, rtol=1e-10)

    def test_wrong_units_raises(self):
        """Passing r0 in wrong units should raise UnitsError."""
        with pytest.raises(u.UnitsError):
            CoreCollapseRate(r0=1e-4 * u.s**-1)

# salpeter_imf
class TestSalpeterIMF:

    def test_power_law(self):
        """Weights must follow M^{-2.35}."""
        masses  = np.array([10.0, 20.0, 40.0])
        weights = salpeter_imf(masses)
        ratios  = weights[:-1] / weights[1:]
        expected = (masses[:-1] / masses[1:])**(-2.35)
        assert np.allclose(ratios, expected, rtol=1e-10)

    def test_positive(self):
        assert np.all(salpeter_imf(np.array([8., 15., 30., 100.])) > 0)

# DSNB construction
class TestDSNBConstruction:

    def test_duplicate_mass_raises(self, sn_models):
        """Duplicate progenitor masses must raise ValueError."""
        doubled = sn_models + [(sn_models[0][0], 13)]
        with pytest.raises(ValueError, match="Duplicate"):
            DSNB(doubled)

    def test_imf_weights_sum_to_one(self, model):
        """IMF weights must sum to 1."""
        assert np.isclose(model._weights.sum(), 1.0, rtol=1e-10)

    def test_imf_weights_positive(self, model):
        """All IMF weights must be positive."""
        assert np.all(model._weights > 0)

    def test_lower_mass_gets_higher_weight(self, model):
        """Salpeter IMF gives more weight to lower-mass progenitors."""
        # masses are sorted; weight[0] is for the 13 Msun progenitor
        assert model._weights[0] > model._weights[1]



# DSNB._oscillation_probs
class TestOscillationProbs:

    def test_no_transformation_is_identity(self):
        p_ee, p_xe = DSNB._oscillation_probs(NoTransformation())
        assert np.isclose(p_ee, 1.0)
        assert np.isclose(p_xe, 0.0)

    def test_none_is_identity(self):
        p_ee, p_xe = DSNB._oscillation_probs(None)
        assert np.isclose(p_ee, 1.0)
        assert np.isclose(p_xe, 0.0)

    def test_nh_p_ee_near_cos2_theta12(self):
        """NH survival probability ~ cos^2(theta_12) ~ 0.68."""
        p_ee, _ = DSNB._oscillation_probs(AdiabaticMSW(mh=MassHierarchy.NORMAL))
        assert 0.60 < p_ee < 0.75, f"NH p_ee={p_ee:.3f}, expected 0.60--0.75"

    def test_ih_p_ee_near_sin2_theta13(self):
        """IH survival probability ~ sin^2(theta_13) ~ 0.02."""
        p_ee, _ = DSNB._oscillation_probs(AdiabaticMSW(mh=MassHierarchy.INVERTED))
        assert p_ee < 0.05, f"IH p_ee={p_ee:.3f}, expected < 0.05"

    def test_nh_ih_p_ee_differ(self):
        """NH and IH must give different survival probabilities."""
        p_nh, _ = DSNB._oscillation_probs(AdiabaticMSW(mh=MassHierarchy.NORMAL))
        p_ih, _ = DSNB._oscillation_probs(AdiabaticMSW(mh=MassHierarchy.INVERTED))
        assert not np.isclose(p_nh, p_ih)

    def test_probabilities_between_zero_and_one(self):
        for xf in [NoTransformation(),
                   AdiabaticMSW(mh=MassHierarchy.NORMAL),
                   AdiabaticMSW(mh=MassHierarchy.INVERTED)]:
            p_ee, p_xe = DSNB._oscillation_probs(xf)
            assert 0.0 <= p_ee <= 1.0
            assert 0.0 <= p_xe <= 1.0

    def test_invalid_xform_raises(self):
        """An object without prob_eebar must raise AttributeError."""
        class Fake:
            pass
        with pytest.raises(AttributeError):
            DSNB._oscillation_probs(Fake())


# DSNB.get_flux — shape and units
class TestGetFluxShapeAndUnits:

    def test_flux_shape(self, model, E, t):
        flux = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        # shape is (N_flavor, N_time, N_energy)
        assert flux.shape == (len(list(Flavor)), 2, len(E))

    def test_flux_units(self, model, E, t):
        flux = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        arr  = flux[Flavor.NU_E_BAR].array
        assert arr.unit.is_equivalent(u.cm**-2 * u.s**-1 * u.MeV**-1)

    def test_nuebar_nonzero(self, model, E, t):
        flux = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        phi  = flux[Flavor.NU_E_BAR].array.value
        assert np.any(phi > 0)

    def test_other_flavors_zero(self, model, E, t):
        """Only NU_E_BAR should carry flux; all other flavors must be zero."""
        flux = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        for flav in Flavor:
            if flav == Flavor.NU_E_BAR:
                continue
            assert np.allclose(flux[flav].array.value, 0.0), \
                f"Flavor {flav} is non-zero"

    def test_steady_state_same_at_all_times(self, model, E):
        """DSNB is steady-state: flux must be identical at all time samples."""
        t_multi = np.linspace(0, 1, 5) * u.yr
        flux    = model.get_flux(t=t_multi, E=E, flavor_xform=NoTransformation())
        phi     = flux[Flavor.NU_E_BAR].array.value   # (1, N_t, N_E)
        assert np.allclose(phi[0, 0, :], phi[0, -1, :])

    def test_flux_positive(self, model, E, t):
        flux = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        phi  = flux[Flavor.NU_E_BAR].array.value
        assert np.all(phi >= 0)

#DSNB.get_flux — oscillation physics
class TestOscillationPhysics:

    def test_nh_ih_fluxes_differ(self, model, E, t):
        """NH and IH must produce different NU_E_BAR fluxes."""
        flux_NH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.NORMAL))
        flux_IH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.INVERTED))
        phi_NH = flux_NH[Flavor.NU_E_BAR].array.value
        phi_IH = flux_IH[Flavor.NU_E_BAR].array.value
        assert not np.allclose(phi_NH, phi_IH)

    def test_nh_larger_than_ih(self, model, E, t):
        """
        NH NU_E_BAR flux must exceed IH: in NH the high-density MSW resonance
        leaves NU_E_BAR mostly intact (p_ee ~ 0.68), while in IH it nearly
        fully converts NU_E_BAR into NU_X_BAR (p_ee ~ 0.02).
        Lunardini & Tamborra (2012) find NH > IH for antineutrinos.
        """
        flux_NH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.NORMAL))
        flux_IH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.INVERTED))
        phi_NH = flux_NH[Flavor.NU_E_BAR].array.value.mean()
        phi_IH = flux_IH[Flavor.NU_E_BAR].array.value.mean()
        assert phi_NH > phi_IH, \
            f"NH mean flux {phi_NH:.2e} not greater than IH {phi_IH:.2e}"

    def test_nh_ih_ratio_consistent_with_probabilities(self, model, E, t):
        """
        NH/IH ratio must be consistent with the ratio of their p_ee values.
        Since dN/dE_osc = p_ee * spec_eb + p_xe * spec_x, and spec_x is
        subdominant, the flux ratio should be roughly p_ee(NH)/p_ee(IH).
        We check that it is at least > 2 (the actual ratio is ~10-30x).
        """
        flux_NH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.NORMAL))
        flux_IH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.INVERTED))
        phi_NH = flux_NH[Flavor.NU_E_BAR].array.value.mean()
        phi_IH = flux_IH[Flavor.NU_E_BAR].array.value.mean()
        assert phi_NH / phi_IH > 2.0

    def test_msw_changes_flux_significantly(self, model, E, t):
        """
        MSW oscillations must change the flux by > 10% relative to
        unoscillated (Lunardini & Tamborra 2012 find ~50-60% for NH).
        """
        flux_0  = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        flux_NH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.NORMAL))
        phi_0  = flux_0[Flavor.NU_E_BAR].array.value.mean()
        phi_NH = flux_NH[Flavor.NU_E_BAR].array.value.mean()
        change = abs(phi_NH - phi_0) / phi_0
        assert change > 0.10, \
            f"MSW change only {change*100:.1f}%, expected > 10%"

    def test_no_transformation_between_nh_and_ih(self, model, E, t):
        """
        Unoscillated flux must lie between NH and IH fluxes, since
        p_ee(NH) ~ 0.68 > 1.0 is not possible — actually unoscillated
        (p_ee=1) gives the maximum possible NU_E_BAR flux.
        """
        flux_0  = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        flux_NH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.NORMAL))
        flux_IH = model.get_flux(t=t, E=E,
                                 flavor_xform=AdiabaticMSW(mh=MassHierarchy.INVERTED))
        phi_0  = flux_0[Flavor.NU_E_BAR].array.value.mean()
        phi_NH = flux_NH[Flavor.NU_E_BAR].array.value.mean()
        phi_IH = flux_IH[Flavor.NU_E_BAR].array.value.mean()
        # Unoscillated has p_ee=1 so it must be the largest
        assert phi_0 > phi_NH > phi_IH

#DSNB.get_flux — cosmology scaling
class TestCosmologyScaling:

    def test_higher_r0_scales_flux_linearly(self, sn_models, E, t):
        """Flux scales linearly with R_CC(0)."""
        m1 = DSNB(sn_models,
                  rcc=CoreCollapseRate(1e-4 * u.yr**-1 * u.Mpc**-3), n_z=200)
        m2 = DSNB(sn_models,
                  rcc=CoreCollapseRate(2e-4 * u.yr**-1 * u.Mpc**-3), n_z=200)
        phi1 = m1.get_flux(t=t, E=E,
                           flavor_xform=NoTransformation())[Flavor.NU_E_BAR].array.value
        phi2 = m2.get_flux(t=t, E=E,
                           flavor_xform=NoTransformation())[Flavor.NU_E_BAR].array.value
        assert np.allclose(phi2 / phi1, 2.0, rtol=1e-3)

    def test_flux_ballpark_magnitude(self, model, E, t):
        """
        Peak DSNB NU_E_BAR flux should be in the range 1e-13 to 1e-10
        cm^-2 s^-1 MeV^-1, consistent with published predictions.
        """
        flux = model.get_flux(t=t, E=E, flavor_xform=NoTransformation())
        peak = flux[Flavor.NU_E_BAR].array.to(
            u.cm**-2 * u.s**-1 * u.MeV**-1).value.max()
        assert 1e-13 < peak < 1e-10, \
            f"Peak flux {peak:.2e} outside expected range 1e-13 to 1e-10"
