# -*- coding: utf-8 -*-
from snewpy.flavor_transformation \
import MassHierarchy, MixingParameters, \
       NoTransformation, AdiabaticMSW, NonAdiabaticMSWH, \
       TwoFlavorDecoherence, ThreeFlavorDecoherence, NeutrinoDecay

from astropy import units as u
from astropy import constants as c
import numpy as np
from numpy import sin, cos, exp, abs

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
    # Adiabatic MSW: normal mass ordering; override default mixing angles.
    xform = AdiabaticMSW(mix_angles=(theta12, theta13, theta23), mh=MassHierarchy.NORMAL)

    assert(xform.prob_ee(t, E) == sin(theta13)**2)
    assert(xform.prob_ex(t, E) == 1. - sin(theta13)**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + sin(theta13)**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - sin(theta13)**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))

    # Test interface using default mixing angles defined in the submodule.
    mixpars = MixingParameters(MassHierarchy.NORMAL)
    th12, th13, th23 = mixpars.get_mixing_angles()

    xform = AdiabaticMSW()

    assert(xform.prob_ee(t, E) == sin(th13)**2)
    assert(xform.prob_ex(t, E) == 1. - sin(th13)**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + sin(th13)**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - sin(th13)**2))

    assert(xform.prob_eebar(t, E) == (cos(th12)*cos(th13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(th12)*cos(th13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(th12)*cos(th13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(th12)*cos(th13))**2))


def test_adiabaticmsw_imo():
    # Adiabatic MSW: inverted mass ordering; override default mixing angles.
    xform = AdiabaticMSW(mix_angles=(theta12, theta13, theta23), mh=MassHierarchy.INVERTED)

    assert(xform.prob_ee(t, E) == (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(theta12)*cos(theta13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(theta12)*cos(theta13))**2))

    assert(xform.prob_eebar(t, E) == sin(theta13)**2)
    assert(xform.prob_exbar(t, E) == 1. - sin(theta13)**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + sin(theta13)**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - sin(theta13)**2))

    # Test interface using default mixing angles defined in the submodule.
    mixpars = MixingParameters(MassHierarchy.INVERTED)
    th12, th13, th23 = mixpars.get_mixing_angles()

    xform = AdiabaticMSW(mh=MassHierarchy.INVERTED)

    assert(xform.prob_ee(t, E) == (sin(th12)*cos(th13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(th12)*cos(th13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(th12)*cos(th13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(th12)*cos(th13))**2))

    assert(xform.prob_eebar(t, E) == sin(th13)**2)
    assert(xform.prob_exbar(t, E) == 1. - sin(th13)**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + sin(th13)**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - sin(th13)**2))


def test_nonadiabaticmswh_nmo():
    # Adiabatic MSW: normal mass ordering; override the default mixing angles.
    xform = NonAdiabaticMSWH(mix_angles=(theta12, theta13, theta23), mh=MassHierarchy.NORMAL)

    assert(xform.prob_ee(t, E) == (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(theta12)*cos(theta13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(theta12)*cos(theta13))**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))

    # Test interface using default mixing angles defined in the submodule.
    mixpars = MixingParameters(MassHierarchy.NORMAL)
    th12, th13, th23 = mixpars.get_mixing_angles()

    xform = NonAdiabaticMSWH()

    assert(xform.prob_ee(t, E) == (sin(th12)*cos(th13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(th12)*cos(th13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(th12)*cos(th13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(th12)*cos(th13))**2))

    assert(xform.prob_eebar(t, E) == (cos(th12)*cos(th13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(th12)*cos(th13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(th12)*cos(th13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(th12)*cos(th13))**2))


def test_nonadiabaticmswh_imo():
    # Adiabatic MSW: inverted mass ordering; override default mixing angles.
    xform = NonAdiabaticMSWH(mix_angles=(theta12, theta13, theta23), mh=MassHierarchy.INVERTED)

    assert(xform.prob_ee(t, E) == (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(theta12)*cos(theta13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(theta12)*cos(theta13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(theta12)*cos(theta13))**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))

    # Test interface using default mixing angles defined in the submodule.
    mixpars = MixingParameters(MassHierarchy.INVERTED)
    th12, th13, th23 = mixpars.get_mixing_angles()

    xform = NonAdiabaticMSWH(mh=MassHierarchy.INVERTED)

    assert(xform.prob_ee(t, E) == (sin(th12)*cos(th13))**2)
    assert(xform.prob_ex(t, E) == 1. - (sin(th12)*cos(th13))**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + (sin(th12)*cos(th13))**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - (sin(th12)*cos(th13))**2))

    assert(xform.prob_eebar(t, E) == (cos(th12)*cos(th13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(th12)*cos(th13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(th12)*cos(th13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(th12)*cos(th13))**2))


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
    # Neutrino decay with NMO, overriding the default mixing angles.
    xform = NeutrinoDecay(mix_angles=(theta12, theta13, theta23), mass=mass3, tau=lifetime, dist=distance, mh=MassHierarchy.NORMAL)

    # Test computation of the decay length.
    _E = 10*u.MeV
    assert(xform.gamma(_E) == mass3*c.c / (_E*lifetime))

    De1 = (cos(theta12) * cos(theta13))**2
    De2 = (sin(theta12) * cos(theta13))**2
    De3 = sin(theta13)**2

    # Check transition probabilities.
    prob_ee = np.asarray([De1*(1.-exp(-xform.gamma(_E)*distance)) + De3*exp(-xform.gamma(_E)*distance) for _E in E])

    assert(np.array_equal(xform.prob_ee(t, E), prob_ee))
    assert(xform.prob_ex(t, E) == De1 + De3)
    assert(xform.prob_xx(t, E) == 1 - 0.5*(De1 + De3))
    assert(np.array_equal(xform.prob_xe(t, E), 0.5*(1 - prob_ee)))

    prob_exbar = np.asarray([De1*(1.-exp(-xform.gamma(_E)*distance)) + De2 + De3*exp(-xform.gamma(_E)*distance) for _E in E])

    assert(xform.prob_eebar(t, E) == De3)
    assert(np.array_equal(xform.prob_exbar(t, E), prob_exbar))
    assert(np.array_equal(xform.prob_xxbar(t, E), 1. - 0.5*prob_exbar))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - De3))


def test_nudecay_nmo_default_mixing():
    # Test interface using default mixing angles defined in the submodule.
    xform = NeutrinoDecay(mass=mass3, tau=lifetime, dist=distance)

    # Check transition probabilities (normal hierarchy is default).
    mixpars = MixingParameters()
    th12, th13, th23 = mixpars.get_mixing_angles()

    De1 = (cos(th12) * cos(th13))**2
    De2 = (sin(th12) * cos(th13))**2
    De3 = sin(th13)**2

    prob_ee = np.asarray([De1*(1.-exp(-xform.gamma(_E)*distance)) + De3*exp(-xform.gamma(_E)*distance) for _E in E])

    assert(np.array_equal(xform.prob_ee(t, E), prob_ee))
    assert(xform.prob_ex(t, E) == De1 + De3)
    assert(xform.prob_xx(t, E) == 1 - 0.5*(De1 + De3))
    assert(np.array_equal(xform.prob_xe(t, E), 0.5*(1 - prob_ee)))

    prob_exbar = np.asarray([De1*(1.-exp(-xform.gamma(_E)*distance)) + De2 + De3*exp(-xform.gamma(_E)*distance) for _E in E])

    assert(xform.prob_eebar(t, E) == De3)
    assert(np.array_equal(xform.prob_exbar(t, E), prob_exbar))
    assert(np.array_equal(xform.prob_xxbar(t, E), 1. - 0.5*prob_exbar))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - De3))


def test_nudecay_imo():
    # Neutrino decay with IMO, overriding the default mixing angles.
    xform = NeutrinoDecay(mix_angles=(theta12, theta13, theta23), mass=mass3, tau=lifetime, dist=distance, mh=MassHierarchy.INVERTED)

    De1 = (cos(theta12) * cos(theta13))**2
    De2 = (sin(theta12) * cos(theta13))**2
    De3 = sin(theta13)**2

    # Check transition probabilities.
    prob_ee = np.asarray([De2*exp(-xform.gamma(_E)*distance) +
                          De3*(1.-exp(-xform.gamma(_E)*distance)) for _E in E])

    assert(np.array_equal(xform.prob_ee(t, E), prob_ee))
    assert(xform.prob_ex(t, E) == De1 + De2)
    assert(xform.prob_xx(t, E) == 1 - 0.5*(De1 + De2))
    assert(np.array_equal(xform.prob_xe(t, E), 0.5*(1 - prob_ee)))

    prob_exbar = np.asarray([De1 + De2*np.exp(-xform.gamma(_E)*distance) +
                             De3*(1-np.exp(-xform.gamma(_E)*distance)) for _E in E])

    assert(xform.prob_eebar(t, E) == De3)
    assert(np.array_equal(xform.prob_exbar(t, E), prob_exbar))
    assert(np.array_equal(xform.prob_xxbar(t, E), 1. - 0.5*prob_exbar))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - De3))


def test_nudecay_imo_default_mixing():
    # Test interface using default mixing angles defined in the submodule.
    xform = NeutrinoDecay(mass=mass3, tau=lifetime, dist=distance, mh=MassHierarchy.INVERTED)

    # Check transition probabilities.
    mixpars = MixingParameters(MassHierarchy.INVERTED)
    th12, th13, th23 = mixpars.get_mixing_angles()

    De1 = (cos(th12) * cos(th13))**2
    De2 = (sin(th12) * cos(th13))**2
    De3 = sin(th13)**2

    prob_ee = np.asarray([De2*exp(-xform.gamma(_E)*distance) +
                          De3*(1.-exp(-xform.gamma(_E)*distance)) for _E in E])

    assert(np.array_equal(xform.prob_ee(t, E), prob_ee))
    assert(xform.prob_ex(t, E) == De1 + De2)
    assert(xform.prob_xx(t, E) == 1 - 0.5*(De1 + De2))
    assert(np.array_equal(xform.prob_xe(t, E), 0.5*(1 - prob_ee)))

    prob_exbar = np.asarray([De1 + De2*np.exp(-xform.gamma(_E)*distance) +
                             De3*(1-np.exp(-xform.gamma(_E)*distance)) for _E in E])

    assert(xform.prob_eebar(t, E) == De3)
    assert(np.array_equal(xform.prob_exbar(t, E), prob_exbar))
    assert(np.array_equal(xform.prob_xxbar(t, E), 1. - 0.5*prob_exbar))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - De3))
