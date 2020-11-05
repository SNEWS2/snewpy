# -*- coding: utf-8 -*-
from snewpy.neutrino import Flavor, MassHierarchy, MixingParameters

from astropy import units as u


def test_mixing_nmo():
    # By default, return mixing parameters for NMO.
    mixpars = MixingParameters()
    assert(True)


def test_mixing_imo():
    # By default, return mixing parameters for NMO.
    mixpars = MixingParameters(MassHierarchy.INVERTED)
    assert(True)

