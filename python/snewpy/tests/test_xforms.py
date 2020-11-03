# -*- coding: utf-8 -*-
from snewpy.FlavorTransformation import *
  
import astropy.units as u

def test_mixing_angles():
    assert(theta12 == 33.44 * u.degree)
    assert(theta13 ==  8.57 * u.degree)
    assert(theta23 == 49.0 * u.degree)

def test_noxform():
    xform = NoTransformation()
    assert(xform.p()==1 and xform.pbar()==1)

def test_2fd():
    xform = TwoFlavorDecoherence()
    assert(xform.p()==0.5 and xform.pbar()==0.5)

def test_3fd():
    xform = ThreeFlavorDecoherence()
    assert(xform.p()==1./3 and xform.pbar()==1./3)
