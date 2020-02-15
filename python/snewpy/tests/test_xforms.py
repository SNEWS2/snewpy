# -*- coding: utf-8 -*-
from snewpy.FlavorTransformation import *
  
import astropy.units as u

import numpy as np

def test_mixing_angles():
    assert(theta12 == np.deg2rad(33.))
    assert(theta13 == np.deg2rad(9.))
    assert(theta23 == np.deg2rad(45.))

def test_noxform():
    xform = NoTransformation()
    assert(xform.p()==1 and xform.pbar()==1)

def test_2fd():
    xform = TwoFlavorDecoherence()
    assert(xform.p()==0.5 and xform.pbar()==0.5)

def test_3fd():
    xform = ThreeFlavorDecoherence()
    assert(xform.p()==1./3 and xform.pbar()==1./3)
