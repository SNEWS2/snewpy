# -*- coding: utf-8 -*-
import unittest

from snewpy.flavor_transformation import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import NoTransformation, \
    CompleteExchange, AdiabaticMSW, NonAdiabaticMSWH, \
    TwoFlavorDecoherence, ThreeFlavorDecoherence, \
    NeutrinoDecay, AdiabaticMSWes, NonAdiabaticMSWes

from astropy import units as u
from astropy import constants as c
import numpy as np
from numpy import sin, cos, exp, abs

class TestFlavorTransformations(unittest.TestCase):

    def setUp(self):
        # Dummy time and energy arrays, with proper dimensions.
        self.t = np.arange(10) * u.s
        self.E = np.linspace(1,100,21) * u.MeV

        # Dummy mixing angles.
        self.theta12 = 33 * u.deg
        self.theta13 =  9 * u.deg
        self.theta23 = 49 * u.deg

        # Dummy neutrino decay parameters; see arXiv:1910.01127.
        self.mass3 = 0.5 * u.eV/c.c**2
        self.lifetime = 1 * u.day
        self.distance = 10 * u.kpc

    def test_noxform(self):
        """
        Survival probabilities for no oscillations
        """
        xform = NoTransformation()

        self.assertEqual(xform.prob_ee(self.t, self.E), 1)
        self.assertEqual(xform.prob_ex(self.t, self.E), 0)
        self.assertEqual(xform.prob_xx(self.t, self.E), 1)
        self.assertEqual(xform.prob_xe(self.t, self.E), 0)

        self.assertEqual(xform.prob_eebar(self.t, self.E), 1)
        self.assertEqual(xform.prob_exbar(self.t, self.E), 0)
        self.assertEqual(xform.prob_xxbar(self.t, self.E), 1)
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0)

    def test_fullex(self):
        """
        Survival probabilities for complete electron->X transformation
        """
        xform = CompleteExchange()

        self.assertEqual(xform.prob_ee(self.t, self.E), 0)
        self.assertEqual(xform.prob_ex(self.t, self.E), 1)
        self.assertEqual(xform.prob_xx(self.t, self.E), 0.5)
        self.assertEqual(xform.prob_xe(self.t, self.E), 0.5)

        self.assertEqual(xform.prob_eebar(self.t, self.E), 0)
        self.assertEqual(xform.prob_exbar(self.t, self.E), 1)
        self.assertEqual(xform.prob_xxbar(self.t, self.E), 0.5)
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5)

    def test_adiabaticmsw_nmo(self):
        """
        Adiabatic MSW with normal ordering
        """
        xform = AdiabaticMSW(mix_angles=(self.theta12, self.theta13, self.theta23), mh=MassHierarchy.NORMAL)

        self.assertEqual(xform.prob_ee(self.t, self.E), sin(self.theta13)**2)
        self.assertEqual(xform.prob_ex(self.t, self.E), 1. - sin(self.theta13)**2)
        self.assertEqual(xform.prob_xx(self.t, self.E), 0.5*(1. + sin(self.theta13)**2))
        self.assertEqual(xform.prob_xe(self.t, self.E), 0.5*(1. - sin(self.theta13)**2))

        self.assertEqual(xform.prob_eebar(self.t, self.E), (cos(self.theta12)*cos(self.theta13))**2)
        self.assertEqual(xform.prob_exbar(self.t, self.E), 1. - (cos(self.theta12)*cos(self.theta13))**2)
        self.assertEqual(xform.prob_xxbar(self.t, self.E), 0.5*(1. + (cos(self.theta12)*cos(self.theta13))**2))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - (cos(self.theta12)*cos(self.theta13))**2))

        # Test interface using default mixing angles defined in the submodule.
        mixpars = MixingParameters(MassHierarchy.NORMAL)
        th12, th13, th23 = mixpars.get_mixing_angles()

        xform = AdiabaticMSW()

        self.assertEqual(xform.prob_ee(self.t, self.E), sin(th13)**2)
        self.assertEqual(xform.prob_ex(self.t, self.E), 1. - sin(th13)**2)
        self.assertEqual(xform.prob_xx(self.t, self.E), 0.5*(1. + sin(th13)**2))
        self.assertEqual(xform.prob_xe(self.t, self.E), 0.5*(1. - sin(th13)**2))

        self.assertEqual(xform.prob_eebar(self.t, self.E), (cos(th12)*cos(th13))**2)
        self.assertEqual(xform.prob_exbar(self.t, self.E), 1. - (cos(th12)*cos(th13))**2)
        self.assertEqual(xform.prob_xxbar(self.t, self.E), 0.5*(1. + (cos(th12)*cos(th13))**2))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - (cos(th12)*cos(th13))**2))


    def test_adiabaticmsw_imo(self):
        """
        Adiabatic MSW with inverted ordering
        """
        # Adiabatic MSW: inverted mass ordering; override default mixing angles.
        xform = AdiabaticMSW(mix_angles=(self.theta12, self.theta13, self.theta23), mh=MassHierarchy.INVERTED)

        self.assertEqual(xform.prob_ee(self.t, self.E), (sin(self.theta12)*cos(self.theta13))**2)
        self.assertEqual(xform.prob_ex(self.t, self.E), 1. - (sin(self.theta12)*cos(self.theta13))**2)
        self.assertEqual(xform.prob_xx(self.t, self.E), 0.5*(1. + (sin(self.theta12)*cos(self.theta13))**2))
        self.assertEqual(xform.prob_xe(self.t, self.E), 0.5*(1. - (sin(self.theta12)*cos(self.theta13))**2))

        self.assertEqual(xform.prob_eebar(self.t, self.E), sin(self.theta13)**2)
        self.assertEqual(xform.prob_exbar(self.t, self.E), 1. - sin(self.theta13)**2)
        self.assertEqual(xform.prob_xxbar(self.t, self.E), 0.5*(1. + sin(self.theta13)**2))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - sin(self.theta13)**2))

        # Test interface using default mixing angles defined in the submodule.
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        th12, th13, th23 = mixpars.get_mixing_angles()

        xform = AdiabaticMSW(mh=MassHierarchy.INVERTED)

        self.assertEqual(xform.prob_ee(self.t, self.E), (sin(th12)*cos(th13))**2)
        self.assertEqual(xform.prob_ex(self.t, self.E), 1. - (sin(th12)*cos(th13))**2)
        self.assertEqual(xform.prob_xx(self.t, self.E), 0.5*(1. + (sin(th12)*cos(th13))**2))
        self.assertEqual(xform.prob_xe(self.t, self.E), 0.5*(1. - (sin(th12)*cos(th13))**2))

        self.assertEqual(xform.prob_eebar(self.t, self.E), sin(th13)**2)
        self.assertEqual(xform.prob_exbar(self.t, self.E), 1. - sin(th13)**2)
        self.assertEqual(xform.prob_xxbar(self.t, self.E), 0.5*(1. + sin(th13)**2))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - sin(th13)**2))


    #def test_nonadiabaticmsw_nmo():
    #    # Adiabatic MSW: normal mass ordering; override the default mixing angles.
    #    xform = NonAdiabaticMSW(mix_angles=(self.theta12, self.theta13, self.theta23), mh=MassHierarchy.NORMAL)
    #
    #    self.assertTrue(xform.prob_ee(self.t, self.E) == (sin(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_ex(self.t, self.E) == 1. - (sin(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_xx(self.t, self.E) == 0.5*(1. + (sin(self.theta12)*cos(self.theta13))**2))
    #    self.assertTrue(xform.prob_xe(self.t, self.E) == 0.5*(1. - (sin(self.theta12)*cos(self.theta13))**2))
    #
    #    self.assertTrue(xform.prob_eebar(self.t, self.E) == (cos(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_exbar(self.t, self.E) == 1. - (cos(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_xxbar(self.t, self.E) == 0.5*(1. + (cos(self.theta12)*cos(self.theta13))**2))
    #    self.assertTrue(xform.prob_xebar(self.t, self.E) == 0.5*(1. - (cos(self.theta12)*cos(self.theta13))**2))
    #
    #    # Test interface using default mixing angles defined in the submodule.
    #    mixpars = MixingParameters(MassHierarchy.NORMAL)
    #    th12, th13, th23 = mixpars.get_mixing_angles()
    #
    #    xform = NonAdiabaticMSW()
    #
    #    self.assertTrue(xform.prob_ee(self.t, self.E) == (sin(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_ex(self.t, self.E) == 1. - (sin(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_xx(self.t, self.E) == 0.5*(1. + (sin(th12)*cos(th13))**2))
    #    self.assertTrue(xform.prob_xe(self.t, self.E) == 0.5*(1. - (sin(th12)*cos(th13))**2))
    #
    #    self.assertTrue(xform.prob_eebar(self.t, self.E) == (cos(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_exbar(self.t, self.E) == 1. - (cos(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_xxbar(self.t, self.E) == 0.5*(1. + (cos(th12)*cos(th13))**2))
    #    self.assertTrue(xform.prob_xebar(self.t, self.E) == 0.5*(1. - (cos(th12)*cos(th13))**2))
    #
    #
    #def test_nonadiabaticmsw_imo():
    #    # Adiabatic MSW: inverted mass ordering; override default mixing angles.
    #    xform = NonAdiabaticMSW(mix_angles=(self.theta12, self.theta13, self.theta23), mh=MassHierarchy.NORMAL)
    #
    #    self.assertTrue(xform.prob_ee(self.t, self.E) == (sin(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_ex(self.t, self.E) == 1. - (sin(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_xx(self.t, self.E) == 0.5*(1. + (sin(self.theta12)*cos(self.theta13))**2))
    #    self.assertTrue(xform.prob_xe(self.t, self.E) == 0.5*(1. - (sin(self.theta12)*cos(self.theta13))**2))
    #
    #    self.assertTrue(xform.prob_eebar(self.t, self.E) == (cos(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_exbar(self.t, self.E) == 1. - (cos(self.theta12)*cos(self.theta13))**2)
    #    self.assertTrue(xform.prob_xxbar(self.t, self.E) == 0.5*(1. + (cos(self.theta12)*cos(self.theta13))**2))
    #    self.assertTrue(xform.prob_xebar(self.t, self.E) == 0.5*(1. - (cos(self.theta12)*cos(self.theta13))**2))
    #
    #    # Test interface using default mixing angles defined in the submodule.
    #    mixpars = MixingParameters(MassHierarchy.INVERTED)
    #    th12, th13, th23 = mixpars.get_mixing_angles()
    #
    #    xform = NonAdiabaticMSW(mh=MassHierarchy.INVERTED)
    #
    #    self.assertTrue(xform.prob_ee(self.t, self.E) == (sin(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_ex(self.t, self.E) == 1. - (sin(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_xx(self.t, self.E) == 0.5*(1. + (sin(th12)*cos(th13))**2))
    #    self.assertTrue(xform.prob_xe(self.t, self.E) == 0.5*(1. - (sin(th12)*cos(th13))**2))
    #
    #    self.assertTrue(xform.prob_eebar(self.t, self.E) == (cos(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_exbar(self.t, self.E) == 1. - (cos(th12)*cos(th13))**2)
    #    self.assertTrue(xform.prob_xxbar(self.t, self.E) == 0.5*(1. + (cos(th12)*cos(th13))**2))
    #    self.assertTrue(xform.prob_xebar(self.t, self.E) == 0.5*(1. - (cos(th12)*cos(th13))**2))


#    def test_2fd(self):
#        # Two-flavor decoherence.
#        xform = TwoFlavorDecoherence()
#
#        self.assertTrue(xform.prob_ee(self.t, self.E) == 0.5)
#        self.assertTrue(xform.prob_ex(self.t, self.E) == 0.5)
#        self.assertTrue(xform.prob_xx(self.t, self.E) == 0.75)
#        self.assertTrue(xform.prob_xe(self.t, self.E) == 0.25)
#
#        self.assertTrue(xform.prob_eebar(self.t, self.E) == 0.5)
#        self.assertTrue(xform.prob_exbar(self.t, self.E) == 0.5)
#        self.assertTrue(xform.prob_xxbar(self.t, self.E) == 0.75)
#        self.assertTrue(xform.prob_xebar(self.t, self.E) == 0.25)
#
#
    def test_3fd(self):
        # Three-flavor decoherence.
        xform = ThreeFlavorDecoherence()

        self.assertEqual(xform.prob_ee(self.t, self.E), 1./3)
        self.assertTrue(abs(xform.prob_ex(self.t, self.E) - 2./3) < 1e-12)
        self.assertTrue(abs(xform.prob_xx(self.t, self.E) - 2./3) < 1e-12)
        self.assertTrue(abs(xform.prob_xe(self.t, self.E) - 1./3) < 1e-12)

        self.assertEqual(xform.prob_eebar(self.t, self.E) , 1./3)
        self.assertTrue(abs(xform.prob_exbar(self.t, self.E) - 2./3) < 1e-12)
        self.assertTrue(abs(xform.prob_xxbar(self.t, self.E) - 2./3) < 1e-12)
        self.assertTrue(abs(xform.prob_xebar(self.t, self.E) - 1./3) < 1e-12)


    def test_nudecay_nmo(self):
        # Neutrino decay with NMO, overriding the default mixing angles.
        xform = NeutrinoDecay(mix_angles=(self.theta12, self.theta13, self.theta23), mass=self.mass3, tau=self.lifetime, dist=self.distance, mh=MassHierarchy.NORMAL)

        # Test computation of the decay length.
        _E = 10*u.MeV
        self.assertTrue(xform.gamma(_E) == self.mass3*c.c / (_E*self.lifetime))

        De1 = (cos(self.theta12) * cos(self.theta13))**2
        De2 = (sin(self.theta12) * cos(self.theta13))**2
        De3 = sin(self.theta13)**2

        # Check transition probabilities.
        prob_ee = np.asarray([De1*(1.-exp(-xform.gamma(_E)*self.distance)) + De3*exp(-xform.gamma(_E)*self.distance) for _E in self.E])

        self.assertTrue(np.array_equal(xform.prob_ee(self.t, self.E), prob_ee))
        self.assertEqual(xform.prob_ex(self.t, self.E), De1 + De3)
        self.assertEqual(xform.prob_xx(self.t, self.E), 1 - 0.5*(De1 + De3))
        self.assertTrue(np.array_equal(xform.prob_xe(self.t, self.E), 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1*(1.-exp(-xform.gamma(_E)*self.distance)) + De2 + De3*exp(-xform.gamma(_E)*self.distance) for _E in self.E])

        self.assertEqual(xform.prob_eebar(self.t, self.E), De3)
        self.assertTrue(np.array_equal(xform.prob_exbar(self.t, self.E), prob_exbar))
        self.assertTrue(np.array_equal(xform.prob_xxbar(self.t, self.E), 1. - 0.5*prob_exbar))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - De3))


    def test_nudecay_nmo_default_mixing(self):
        # Test interface using default mixing angles defined in the submodule.
        xform = NeutrinoDecay(mass=self.mass3, tau=self.lifetime, dist=self.distance)

        # Check transition probabilities (normal hierarchy is default).
        mixpars = MixingParameters()
        th12, th13, th23 = mixpars.get_mixing_angles()

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        prob_ee = np.asarray([De1*(1.-exp(-xform.gamma(_E)*self.distance)) + De3*exp(-xform.gamma(_E)*self.distance) for _E in self.E])

        self.assertTrue(np.array_equal(xform.prob_ee(self.t, self.E), prob_ee))
        self.assertEqual(xform.prob_ex(self.t, self.E), De1 + De3)
        self.assertEqual(xform.prob_xx(self.t, self.E), 1 - 0.5*(De1 + De3))
        self.assertTrue(np.array_equal(xform.prob_xe(self.t, self.E), 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1*(1.-exp(-xform.gamma(_E)*self.distance)) + De2 + De3*exp(-xform.gamma(_E)*self.distance) for _E in self.E])

        self.assertEqual(xform.prob_eebar(self.t, self.E), De3)
        self.assertTrue(np.array_equal(xform.prob_exbar(self.t, self.E), prob_exbar))
        self.assertTrue(np.array_equal(xform.prob_xxbar(self.t, self.E), 1. - 0.5*prob_exbar))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - De3))


    def test_nudecay_imo(self):
        # Neutrino decay with IMO, overriding the default mixing angles.
        xform = NeutrinoDecay(mix_angles=(self.theta12, self.theta13, self.theta23), mass=self.mass3, tau=self.lifetime, dist=self.distance, mh=MassHierarchy.INVERTED)

        De1 = (cos(self.theta12) * cos(self.theta13))**2
        De2 = (sin(self.theta12) * cos(self.theta13))**2
        De3 = sin(self.theta13)**2

        # Check transition probabilities.
        prob_ee = np.asarray([De2*exp(-xform.gamma(_E)*self.distance) + De3*(1.-exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        self.assertTrue(np.array_equal(xform.prob_ee(self.t, self.E), prob_ee))
        self.assertEqual(xform.prob_ex(self.t, self.E), De1 + De2)
        self.assertEqual(xform.prob_xx(self.t, self.E), 1 - 0.5*(De1 + De2))
        self.assertTrue(np.array_equal(xform.prob_xe(self.t, self.E), 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1 + De2*np.exp(-xform.gamma(_E)*self.distance) +
                                 De3*(1-np.exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        self.assertEqual(xform.prob_eebar(self.t, self.E), De3)
        self.assertTrue(np.array_equal(xform.prob_exbar(self.t, self.E), prob_exbar))
        self.assertTrue(np.array_equal(xform.prob_xxbar(self.t, self.E), 1. - 0.5*prob_exbar))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - De3))


    def test_nudecay_imo_default_mixing(self):
        # Test interface using default mixing angles defined in the submodule.
        xform = NeutrinoDecay(mass=self.mass3, tau=self.lifetime, dist=self.distance, mh=MassHierarchy.INVERTED)

        # Check transition probabilities.
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        th12, th13, th23 = mixpars.get_mixing_angles()

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        prob_ee = np.asarray([De2*exp(-xform.gamma(_E)*self.distance) +
                              De3*(1.-exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        self.assertTrue(np.array_equal(xform.prob_ee(self.t, self.E), prob_ee))
        self.assertEqual(xform.prob_ex(self.t, self.E), De1 + De2)
        self.assertEqual(xform.prob_xx(self.t, self.E), 1 - 0.5*(De1 + De2))
        self.assertTrue(np.array_equal(xform.prob_xe(self.t, self.E), 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1 + De2*np.exp(-xform.gamma(_E)*self.distance) +
                                 De3*(1-np.exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        self.assertEqual(xform.prob_eebar(self.t, self.E), De3)
        self.assertTrue(np.array_equal(xform.prob_exbar(self.t, self.E), prob_exbar))
        self.assertTrue(np.array_equal(xform.prob_xxbar(self.t, self.E), 1. - 0.5*prob_exbar))
        self.assertEqual(xform.prob_xebar(self.t, self.E), 0.5*(1. - De3))
