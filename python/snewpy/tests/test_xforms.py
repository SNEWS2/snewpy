# -*- coding: utf-8 -*-
#from snewpy.FlavorTransformation import *
from snewpy.flavor_transformation import NoTransformation, AdiabaticMSW_NMO

from astropy import units as u
import numpy as np
from numpy import sin, cos

# Dummy time and energy arrays, with proper dimensions.
t = np.arange(10) * u.s
E = np.linspace(1,100,21) * u.MeV

# Dummy mixing angles.
theta12 = 33 * u.deg
theta13 =  9 * u.deg
theta23 = 49 * u.deg

def test_noxform():
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
    xform = AdiabaticMSW_NMO(theta12, theta13, theta23)

    assert(xform.prob_ee(t, E) == sin(theta13)**2)
    assert(xform.prob_ex(t, E) == 1. - sin(theta13)**2)
    assert(xform.prob_xx(t, E) == 0.5*(1. + sin(theta13)**2))
    assert(xform.prob_xe(t, E) == 0.5*(1. - sin(theta13)**2))

    assert(xform.prob_eebar(t, E) == (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_exbar(t, E) == 1. - (cos(theta12)*cos(theta13))**2)
    assert(xform.prob_xxbar(t, E) == 0.5*(1. + (cos(theta12)*cos(theta13))**2))
    assert(xform.prob_xebar(t, E) == 0.5*(1. - (cos(theta12)*cos(theta13))**2))

#def test_mixing_angles():
#    assert(theta12 == 33.44 * u.degree)
#    assert(theta13 ==  8.57 * u.degree)
#    assert(theta23 == 49.0 * u.degree)
#
#def test_noxform():
#    # No transformation.
#    xform = NoTransformation()
#    assert(xform.p()==1 and xform.pbar()==1)
#
#def test_adiabaticmsw_nmo():
#    # Adiabatic MSW (normal ordering).
#    xform = AdiabaticMSW_NMO()
#    assert(xform.p()==sin(theta13)**2)
#    assert(xform.pbar()==(cos(theta12)*cos(theta13))**2)
#
#def test_adiabaticmsw_imo():
#    # Adiabatic MSW (inverted ordering).
#    xform = AdiabaticMSW_IMO()
#    assert(xform.p()==(sin(theta12)*cos(theta13))**2)
#    assert(xform.pbar()==sin(theta13)**2)
#
#def test_2fd():
#    # Two-flavor decoherence.
#    xform = TwoFlavorDecoherence()
#    assert(xform.p()==0.5)
#    assert(xform.pbar()==0.5)
#
#def test_3fd():
#    # Three-flavor decoherence.
#    xform = ThreeFlavorDecoherence()
#    assert(xform.p()==1./3)
#    assert(xform.pbar()==1./3)
