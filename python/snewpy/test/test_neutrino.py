# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""

import unittest

from snewpy.neutrino import Flavor, MassHierarchy, MixingParameters

from astropy import units as u

class TestNeutrino(unittest.TestCase):

    def test_flavor(self):
        nue  = Flavor.NU_E
        self.assertTrue(nue.is_electron)
        self.assertTrue(nue.is_neutrino)
        self.assertTrue(not nue.is_antineutrino)

        nux  = Flavor.NU_X
        self.assertTrue(not nux.is_electron)
        self.assertTrue(nux.is_neutrino)
        self.assertTrue(not nux.is_antineutrino)

        nueb = Flavor.NU_E_BAR
        self.assertTrue(nueb.is_electron)
        self.assertTrue(not nueb.is_neutrino)
        self.assertTrue(nueb.is_antineutrino)

        nuxb = Flavor.NU_X_BAR
        self.assertTrue(not nuxb.is_electron)
        self.assertTrue(not nuxb.is_neutrino)
        self.assertTrue(nuxb.is_antineutrino)


    def test_mixing_nmo(self):
        # By default, return mixing parameters for NMO.
        mixpars = MixingParameters()
        self.assertTrue(mixpars.theta12 == 33.44 * u.deg)
        self.assertTrue(mixpars.theta13 ==  8.57 * u.deg)
        self.assertTrue(mixpars.theta23 == 49.20 * u.deg)
        self.assertTrue(mixpars.deltaCP == 197 * u.deg)
        self.assertTrue(mixpars.dm21_2  == 7.42e-5 * u.eV**2)
        self.assertTrue(mixpars.dm32_2  == 2.517e-3 * u.eV**2)


    def test_mixing_imo(self):
        # By default, return mixing parameters for NMO.
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        self.assertTrue(mixpars.theta12 == 33.45 * u.deg)
        self.assertTrue(mixpars.theta13 ==  8.60 * u.deg)
        self.assertTrue(mixpars.theta23 == 49.30 * u.deg)
        self.assertTrue(mixpars.deltaCP == 282 * u.deg)
        self.assertTrue(mixpars.dm21_2  == 7.42e-5 * u.eV**2)
        self.assertTrue(mixpars.dm31_2  == -2.498e-3 * u.eV**2)
