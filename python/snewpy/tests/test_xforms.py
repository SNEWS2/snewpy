# -*- coding: utf-8 -*-
from snewpy.flavor_transformation \
import MassHierarchy, NoTransformation, AdiabaticMSW, NonAdiabaticMSW, \
       TwoFlavorDecoherence, ThreeFlavorDecoherence, NeutrinoDecay

from astropy import units as u
from astropy import constants as c
import numpy as np
from numpy import sin, cos, abs

# Dummy time and energy arrays, with proper dimensions.
t = np.arange(10) * u.s
E = np.linspace(1,100,21) * u.MeV

# Dummy mixing angles.
theta12 = 33 * u.deg
theta13 =  9 * u.deg
theta23 = 49 * u.deg

# Dummy neutrino decay parameters; see arXiv:1910.01127.
mass3 = 0.5 * u.eV/c.c**2
lifetime = 1 * u.day
distance = 10 * u.kpc


def test_noxform():
    # No transformations.
    xform = NoTransformation()

    assert(xform.prob_ee(t, E) == 1)
    assert(xform.prob_ex(t, E) == 0)
    assert(xform.prob_xx(t, E) == 1)
    assert(xform.prob_xe(t, E) == 0)

    assert(xform.prob_eebar(t, E) == 1)
    assert(xform.prob_exbar(t, E) == 0)
    assert(xform.prob_xxbar(t, E) == 1)
    assert(xform.prob_xebar(t, E) == 0)


def test_adiabaticmsw_nmo():
    # Adiabatic MSW: normal mass ordering.
    xform = AdiabaticMSW(theta12, theta13, theta23, mh=MassHierarchy.NORMAL)

    assert(xform.prob_ee(t, E) == sin(theta13)**2)
    assert(xform.prob_ex(t, E) == 1. - sin(theta13)**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + sin(theta13)**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - sin(theta13)**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))


def test_adiabaticmsw_imo():
    # Adiabatic MSW: inverted mass ordering.
    xform = AdiabaticMSW(theta12, theta13, theta23, mh=MassHierarchy.INVERTED)

    assert(xform.prob_ee(t, E) == (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(theta12)*cos(theta13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(theta12)*cos(theta13))**2))

    assert(xform.prob_eebar(t, E) == sin(theta13)**2)
    assert(xform.prob_exbar(t, E) == 1. - sin(theta13)**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + sin(theta13)**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - sin(theta13)**2))


def test_nonadiabaticmsw_nmo():
    # Adiabatic MSW: normal mass ordering.
    xform = NonAdiabaticMSW(theta12, theta13, theta23, mh=MassHierarchy.NORMAL)

    assert(xform.prob_ee(t, E) == (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(theta12)*cos(theta13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(theta12)*cos(theta13))**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))


def test_nonadiabaticmsw_imo():
    # Adiabatic MSW: inverted mass ordering.
    xform = NonAdiabaticMSW(theta12, theta13, theta23, mh=MassHierarchy.NORMAL)

    assert(xform.prob_ee(t, E) == (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(theta12)*cos(theta13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(theta12)*cos(theta13))**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))


def test_2fd():
    # Two-flavor decoherence.
    xform = TwoFlavorDecoherence()

    assert(xform.prob_ee(t, E) == 0.5)
    assert(xform.prob_ex(t, E) == 0.5)
    assert(xform.prob_xx(t, E) == 0.75)
    assert(xform.prob_xe(t, E) == 0.25)

    assert(xform.prob_eebar(t, E) == 0.5)
    assert(xform.prob_exbar(t, E) == 0.5)
    assert(xform.prob_xxbar(t, E) == 0.75)
    assert(xform.prob_xebar(t, E) == 0.25)


def test_3fd():
    # Three-flavor decoherence.
    xform = ThreeFlavorDecoherence()

    assert(xform.prob_ee(t, E) == 1./3)
    assert(abs(xform.prob_ex(t, E) - 2./3) < 1e-12)
    assert(abs(xform.prob_xx(t, E) - 2./3) < 1e-12)
    assert(abs(xform.prob_xe(t, E) - 1./3) < 1e-12)

    assert(xform.prob_eebar(t, E) == 1./3)
    assert(abs(xform.prob_exbar(t, E) - 2./3) < 1e-12)
    assert(abs(xform.prob_xxbar(t, E) - 2./3) < 1e-12)
    assert(abs(xform.prob_xebar(t, E) - 1./3) < 1e-12)


def test_nudecay_nmo():
    # Neutrino decay.
    xform = NeutrinoDecay(theta12, theta13, theta23, mass3, lifetime, distance,
                          mh=MassHierarchy.NORMAL)

    E = 10*u.MeV
    assert(xform.gamma(E) == mass3*c.c / (E*lifetime))
