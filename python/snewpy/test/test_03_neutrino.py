# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""

import unittest

from snewpy.neutrino import Flavor, MassHierarchy, MixingParameters

from astropy import units as u

class TestNeutrino(unittest.TestCase):

    def test_flavor(self):
        """
        Neutrino flavor types
        """
        nue  = Flavor.NU_E
        self.assertTrue(nue.is_electron)
        self.assertTrue(nue.is_neutrino)
        self.assertFalse(nue.is_antineutrino)

        nux  = Flavor.NU_X
        self.assertFalse(nux.is_electron)
        self.assertTrue(nux.is_neutrino)
        self.assertFalse(nux.is_antineutrino)

        nueb = Flavor.NU_E_BAR
        self.assertTrue(nueb.is_electron)
        self.assertFalse(nueb.is_neutrino)
        self.assertTrue(nueb.is_antineutrino)

        nuxb = Flavor.NU_X_BAR
        self.assertFalse(nuxb.is_electron)
        self.assertFalse(nuxb.is_neutrino)
        self.assertTrue(nuxb.is_antineutrino)


    def test_mixing_nmo(self):
        """
        Mixing parameter values; NMO
        """
        # By default, return mixing parameters for NMO.
        mixpars = MixingParameters()
        self.assertEqual(mixpars.theta12, 33.44 * u.deg)
        self.assertEqual(mixpars.theta13,  8.57 * u.deg)
        self.assertEqual(mixpars.theta23, 49.20 * u.deg)
        self.assertEqual(mixpars.deltaCP, 197 * u.deg)
        self.assertEqual(mixpars.dm21_2, 7.42e-5 * u.eV**2)
        self.assertEqual(mixpars.dm32_2, 2.517e-3 * u.eV**2)


    def test_mixing_imo(self):
        """
        Mixing parameter values; IMO
        """
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        self.assertEqual(mixpars.theta12, 33.45 * u.deg)
        self.assertEqual(mixpars.theta13,  8.60 * u.deg)
        self.assertEqual(mixpars.theta23, 49.30 * u.deg)
        self.assertEqual(mixpars.deltaCP, 282 * u.deg)
        self.assertEqual(mixpars.dm21_2, 7.42e-5 * u.eV**2)
        self.assertEqual(mixpars.dm31_2, -2.498e-3 * u.eV**2)
