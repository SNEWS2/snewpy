# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""

import unittest

from snewpy.neutrino import TwoFlavor, ThreeFlavor, FourFlavor, MassHierarchy, MixingParameters

from astropy import units as u

class TestNeutrino(unittest.TestCase):

    def test_twoflavor(self):
        """
        Neutrino flavor types
        """
        nue  = TwoFlavor.NU_E
        self.assertTrue(nue.is_electron)
        self.assertTrue(nue.is_neutrino)
        self.assertFalse(nue.is_antineutrino)

        nux  = TwoFlavor.NU_X
        self.assertFalse(nux.is_electron)
        self.assertTrue(nux.is_neutrino)
        self.assertFalse(nux.is_antineutrino)

        nueb = TwoFlavor.NU_E_BAR
        self.assertTrue(nueb.is_electron)
        self.assertFalse(nueb.is_neutrino)
        self.assertTrue(nueb.is_antineutrino)

        nuxb = TwoFlavor.NU_X_BAR
        self.assertFalse(nuxb.is_electron)
        self.assertFalse(nuxb.is_neutrino)
        self.assertTrue(nuxb.is_antineutrino)

    def test_threeflavor(self):
        """
        Neutrino flavor types
        """
        nue  = ThreeFlavor.NU_E
        self.assertTrue(nue.is_electron)
        self.assertTrue(nue.is_neutrino)
        self.assertFalse(nue.is_antineutrino)

        numu  = ThreeFlavor.NU_MU
        self.assertFalse(numu.is_electron)
        self.assertTrue(numu.is_neutrino)
        self.assertFalse(numu.is_antineutrino)

        nutau  = ThreeFlavor.NU_TAU
        self.assertFalse(nutau.is_electron)
        self.assertTrue(nutau.is_neutrino)
        self.assertFalse(nutau.is_antineutrino)        

        nueb = ThreeFlavor.NU_E_BAR
        self.assertTrue(nueb.is_electron)
        self.assertFalse(nueb.is_neutrino)
        self.assertTrue(nueb.is_antineutrino)

        numub  = ThreeFlavor.NU_MU_NAR
        self.assertFalse(numub.is_electron)
        self.assertTrue(numub.is_neutrino)
        self.assertFalse(numub.is_antineutrino)

        nutaub  = ThreeFlavor.NU_TAU_BAR
        self.assertFalse(nutaub.is_electron)
        self.assertTrue(nutaub.is_neutrino)
        self.assertFalse(nutaub.is_antineutrino)  

        def test_fourflavor(self):
        """
        Neutrino flavor types
        """
        nue  = FourFlavor.NU_E
        self.assertTrue(nue.is_electron)
        self.assertTrue(nue.is_neutrino)
        self.assertFalse(nue.is_antineutrino)

        numu  = FourFlavor.NU_MU
        self.assertFalse(numu.is_electron)
        self.assertTrue(numu.is_neutrino)
        self.assertFalse(numu.is_antineutrino)

        nutau  = FourFlavor.NU_TAU
        self.assertFalse(nutau.is_electron)
        self.assertTrue(nutau.is_neutrino)
        self.assertFalse(nutau.is_antineutrino)        

        nus  = FourFlavor.NU_S
        self.assertFalse(nus.is_electron)
        self.assertTrue(nus.is_neutrino)
        self.assertFalse(nus.is_antineutrino)        

        nueb = FourFlavor.NU_E_BAR
        self.assertTrue(nueb.is_electron)
        self.assertFalse(nueb.is_neutrino)
        self.assertTrue(nueb.is_antineutrino)

        numub  = FourFlavor.NU_MU_NAR
        self.assertFalse(numub.is_electron)
        self.assertTrue(numub.is_neutrino)
        self.assertFalse(numub.is_antineutrino)

        nutaub  = FourFlavor.NU_TAU_BAR
        self.assertFalse(nutaub.is_electron)
        self.assertTrue(nutaub.is_neutrino)
        self.assertFalse(nutaub.is_antineutrino)  
   
        nusb  = FourFlavor.NU_S_BAR
        self.assertFalse(nusb.is_electron)
        self.assertTrue(nusb.is_neutrino)
        self.assertFalse(nusb.is_antineutrino)  


    def test_mixing_nmo(self):
        """
        Mixing parameter values; NMO
        """
        # By default, return mixing parameters for NMO.
        mixpars = MixingParameters(version="NuFIT5.0")
        self.assertEqual(mixpars.theta12, 33.44 * u.deg)
        self.assertEqual(mixpars.theta13,  8.57 * u.deg)
        self.assertEqual(mixpars.theta23, 49.20 * u.deg)
        self.assertEqual(mixpars.deltaCP, 197 * u.deg)
        self.assertEqual(mixpars.dm21_2, 7.42e-5 * u.eV**2)
        self.assertEqual(mixpars.dm31_2, 2.517e-3 * u.eV**2)


    def test_mixing_imo(self):
        """
        Mixing parameter values; IMO
        """
        mixpars = MixingParameters(MassHierarchy.INVERTED, version="NuFIT5.0")
        self.assertEqual(mixpars.theta12, 33.45 * u.deg)
        self.assertEqual(mixpars.theta13,  8.60 * u.deg)
        self.assertEqual(mixpars.theta23, 49.30 * u.deg)
        self.assertEqual(mixpars.deltaCP, 282 * u.deg)
        self.assertEqual(mixpars.dm21_2, 7.42e-5 * u.eV**2)
        self.assertEqual(mixpars.dm32_2, -2.498e-3 * u.eV**2)
