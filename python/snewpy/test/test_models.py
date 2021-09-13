# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""
import unittest

from snewpy.neutrino import Flavor
from snewpy.flavor_transformation import NoTransformation
from snewpy.models import Nakazato_2013, Tamborra_2014, OConnor_2015, \
                          Sukhbold_2015, Bollig_2016, Walk_2018

from astropy import units as u

import numpy as np

class TestModels(unittest.TestCase):

    def test_Nakazato_2013(self):
        """
        Instantiate a set of 'Nakazato 2013' models
        """
        for z in [0.004, 0.02]:
            for trev in [100, 200, 300]:
                for mass in [13., 20., 50.]:
                    mfile = 'models/Nakazato_2013/nakazato-shen-z{}-t_rev{}ms-s{:.1f}.fits'.format(z, trev, mass)
                    model = Nakazato_2013(mfile)

                    self.assertEqual(model.EOS, 'SHEN')
                    self.assertEqual(model.progenitor_mass, mass*u.Msun)
                    self.assertEqual(model.revival_time, trev*u.ms)
                    self.assertEqual(model.metallicity, z)

                    # Check that times are in proper units.
                    t = model.get_time()
                    self.assertTrue(t.unit, u.s)
                    self.assertEqual(model.time[0], -50*u.ms)

                    # Check that we can compute flux dictionaries.
                    f = model.get_initial_spectra(0*u.s, 10*u.MeV)
                    self.assertEqual(type(f), dict)
                    self.assertEqual(len(f), len(Flavor))
                    self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Tamborra_2014(self):
        """
        Instantiate a set of 'Tamborra 2014' models
        """
        for mass in [20., 27.]:
            mfile = 'models/Tamborra_2014/s{:.1f}c_3D_dir1'.format(mass)
            model = Tamborra_2014(mfile, eos='LS220')

            self.assertEqual(model.EOS, 'LS220')
            self.assertEqual(model.progenitor_mass, mass*u.Msun)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_OConnor_2015(self):
        """
        Instantiate a set of "O'Connor 2015" models
        """
        mfile = 'models/OConnor_2015/M1_neutrinos.dat'
        model = OConnor_2015(mfile, eos='LS220')

        self.assertEqual(model.EOS, 'LS220')
        self.assertEqual(model.progenitor_mass, 40*u.Msun)

        # Check that times are in proper units.
        t = model.get_time()
        self.assertTrue(t.unit, u.s)

        # Check that we can compute flux dictionaries.
        f = model.get_initial_spectra(0*u.s, 10*u.MeV)
        self.assertEqual(type(f), dict)
        self.assertEqual(len(f), len(Flavor))
        self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Sukhbold_2015(self):
        """
        Instantiate a set of 'Sukhbold 2015' models
        """
        for mass in ['z9.6', 's27.0']:
            for eos in ['LS220', 'SFHo']:
                mfile = 'models/Sukhbold_2015/sukhbold-{}-{}.fits'.format(eos, mass)
                massval = float(mass[1:]) * u.Msun
                model = Sukhbold_2015(mfile)

                self.assertEqual(model.EOS, eos)
                self.assertEqual(model.progenitor_mass, massval)

                # Check that times are in proper units.
                t = model.get_time()
                self.assertTrue(t.unit, u.s)

                # Check that we can compute flux dictionaries.
                f = model.get_initial_spectra(0*u.s, 10*u.MeV)
                self.assertEqual(type(f), dict)
                self.assertEqual(len(f), len(Flavor))
                self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Bollig_2016(self):
        """
        Instantiate a set of 'Bollig 2016' models
        """
        for mass in [11.2, 27.]:
            mfile = 'models/Bollig_2016/s{:.1f}c'.format(mass)
            model = Bollig_2016(mfile, eos='LS220')

            self.assertEqual(model.EOS, 'LS220')
            self.assertEqual(model.progenitor_mass, mass*u.Msun)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Walk_2018(self):
        """
        Instantiate a set of 'Walk 2018' models
        """
        mass = 15.
        mfile = 'models/Walk_2018/s{:.1f}c_3D_nonrot_dir1'.format(mass)
        model = Walk_2018(mfile, eos='LS220')

        self.assertEqual(model.EOS, 'LS220')
        self.assertEqual(model.progenitor_mass, mass*u.Msun)

        # Check that times are in proper units.
        t = model.get_time()
        self.assertTrue(t.unit, u.s)

        # Check that we can compute flux dictionaries.
        f = model.get_initial_spectra(0*u.s, 10*u.MeV)
        self.assertEqual(type(f), dict)
        self.assertEqual(len(f), len(Flavor))
        self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

