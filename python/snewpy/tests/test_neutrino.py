# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""

from snewpy.neutrino import Flavor, MassHierarchy, MixingParameters

from astropy import units as u


def test_flavor():
    nue  = Flavor.NU_E
    assert(nue.is_electron)
    assert(nue.is_neutrino)
    assert(not nue.is_antineutrino)

    nux  = Flavor.NU_X
    assert(not nux.is_electron)
    assert(nux.is_neutrino)
    assert(not nux.is_antineutrino)

    nueb = Flavor.NU_E_BAR
    assert(nueb.is_electron)
    assert(not nueb.is_neutrino)
    assert(nueb.is_antineutrino)

    nuxb = Flavor.NU_X_BAR
    assert(not nuxb.is_electron)
    assert(not nuxb.is_neutrino)
    assert(nuxb.is_antineutrino)


def test_mixing_nmo():
    # By default, return mixing parameters for NMO.
    mixpars = MixingParameters()
    assert(mixpars.theta12 == 33.44 * u.deg)
    assert(mixpars.theta13 ==  8.57 * u.deg)
    assert(mixpars.theta23 == 49.20 * u.deg)
    assert(mixpars.deltaCP == 197 * u.deg)
    assert(mixpars.dm21_2  == 7.42e-5 * u.eV**2)
    assert(mixpars.dm32_2  == 2.517e-3 * u.eV**2)


def test_mixing_imo():
    # By default, return mixing parameters for NMO.
    mixpars = MixingParameters(MassHierarchy.INVERTED)
    assert(mixpars.theta12 == 33.45 * u.deg)
    assert(mixpars.theta13 ==  8.60 * u.deg)
    assert(mixpars.theta23 == 49.30 * u.deg)
    assert(mixpars.deltaCP == 282 * u.deg)
    assert(mixpars.dm21_2  == 7.42e-5 * u.eV**2)
    assert(mixpars.dm31_2  == -2.498e-3 * u.eV**2)
