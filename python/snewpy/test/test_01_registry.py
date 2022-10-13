# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""
import unittest

from snewpy.neutrino import Flavor
from snewpy.flavor_transformation import NoTransformation
from snewpy.models.ccsn import Nakazato_2013, Tamborra_2014, OConnor_2013, OConnor_2015, \
                          Sukhbold_2015, Bollig_2016, Walk_2018, \
                          Walk_2019, Fornax_2019, Warren_2020, \
                          Kuroda_2020, Fornax_2021, Zha_2021

from astropy import units as u

import numpy as np


class TestModels(unittest.TestCase):

    def test_param_combinations(self):
        models = [Nakazato_2013,
                  Tamborra_2014,
                  OConnor_2013,
                  OConnor_2015,
                  Sukhbold_2015,
                  Bollig_2016,
                  Walk_2018,
                  Walk_2019,
                  Fornax_2019,
                  Warren_2020,
                  Kuroda_2020,
                  Fornax_2021,
                  Zha_2021]

        for model in models:
            for pc in model.param_combinations:
                model(**pc)  # Initialize

    def test_Nakazato_2013(self):
        """
        Instantiate a set of 'Nakazato 2013' models
        """
        for z in [0.004, 0.02]:
            for trev in [100, 200, 300]:
                for mass in [13., 20., 50.]:
                    model = Nakazato_2013(progenitor_mass=mass * u.Msun, metallicity=z, eos='shen',
                                          revival_time=trev * u.ms)

                    self.assertEqual(model.metadata['EOS'], 'shen')
                    self.assertEqual(model.metadata['Progenitor mass'], mass*u.Msun)
                    self.assertEqual(model.metadata['Revival time'], trev*u.ms)
                    self.assertEqual(model.metadata['Metallicity'], z)

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
            model = Tamborra_2014(progenitor_mass=mass*u.Msun, eos='LS220')

            self.assertEqual(model.metadata['EOS'], 'LS220')
            self.assertEqual(model.metadata['Progenitor mass'], mass*u.Msun)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_OConnor_2013(self):
        """
        Instantiate a set of "O'Connor 2015" models
        """
        for mass in list(range(12, 34)) + list(range(35, 61, 5)) + [70, 80, 100, 120]:
            for eos in ['LS220', 'HShen']:
                model = OConnor_2013(progenitor_mass=mass*u.Msun, eos=eos)

                self.assertEqual(model.metadata['EOS'], eos)
                self.assertEqual(model.metadata['Progenitor mass'].value, mass)

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
        model = OConnor_2015(progenitor_mass=40*u.Msun, eos='LS220')

        self.assertEqual(model.metadata['EOS'], 'LS220')
        self.assertEqual(model.metadata['Progenitor mass'], 40*u.Msun)

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
        for mass in [9.6, 27.]:
            for eos in ['LS220', 'SFHo']:
                model = Sukhbold_2015(progenitor_mass=mass * u.Msun, eos=eos)

                self.assertEqual(model.metadata['EOS'], eos)
                self.assertEqual(model.metadata['Progenitor mass'].value, mass)

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
            model = Bollig_2016(progenitor_mass=mass*u.Msun, eos='LS220')

            self.assertEqual(model.metadata['EOS'], 'LS220')
            self.assertEqual(model.metadata['Progenitor mass'], mass*u.Msun)

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
        model = Walk_2018(progenitor_mass=mass*u.Msun)

        self.assertEqual(model.metadata['EOS'], 'LS220')
        self.assertEqual(model.metadata['Progenitor mass'], mass*u.Msun)

        # Check that times are in proper units.
        t = model.get_time()
        self.assertTrue(t.unit, u.s)

        # Check that we can compute flux dictionaries.
        f = model.get_initial_spectra(0*u.s, 10*u.MeV)
        self.assertEqual(type(f), dict)
        self.assertEqual(len(f), len(Flavor))
        self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Walk_2019(self):
        """
        Instantiate a set of 'Walk 2019' models
        """
        model = Walk_2019(progenitor_mass=40*u.Msun, eos='LS220')

        self.assertEqual(model.metadata['EOS'], 'LS220')
        self.assertEqual(model.metadata['Progenitor mass'], 40*u.Msun)

        # Check that times are in proper units.
        t = model.get_time()
        self.assertTrue(t.unit, u.s)

        # Check that we can compute flux dictionaries.
        f = model.get_initial_spectra(0*u.s, 10*u.MeV)
        self.assertEqual(type(f), dict)
        self.assertEqual(len(f), len(Flavor))
        self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Fornax_2019(self):
        """
        Instantiate a set of 'Fornax 2019' models
        """
        for mass in [9, 10, 12, 13, 14, 15, 19, 25, 60]:
            model = Fornax_2019(progenitor_mass=mass*u.Msun)

            self.assertEqual(model.metadata['Progenitor mass'].value, mass)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV, theta=23*u.degree, phi=22*u.degree)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, u.erg/(u.MeV * u.s))

    def test_Warren_2020(self):
        """
        Instantiate a set of 'Warren 2020' models
        """
        masses = np.concatenate((np.linspace(9.0, 12.75, 16),
                                 np.linspace(13, 30., 171),
                                 np.linspace(31., 33., 3),
                                 np.linspace(35, 55, 5),
                                 np.linspace(60, 80, 3),
                                 np.linspace(100, 120, 2)))

        for mixing in [1.23, 1.25, 1.27]:
            for mass in masses:
                model = Warren_2020(progenitor_mass=mass * u.Msun, turbmixing_param=mixing)

                self.assertEqual(model.metadata['Progenitor mass'], float(mass) * u.Msun)
                self.assertEqual(model.metadata['Turb. mixing param.'], mixing)
                self.assertEqual(model.metadata['EOS'], 'SFHo')

                # Check that times are in proper units.
                t = model.get_time()
                self.assertTrue(t.unit, u.s)

                # Check that we can compute flux dictionaries.
                f = model.get_initial_spectra(0 * u.s, 10 * u.MeV)
                self.assertEqual(type(f), dict)
                self.assertEqual(len(f), len(Flavor))
                self.assertEqual(f[Flavor.NU_E].unit, 1. / (u.erg * u.s))

    def test_Kuroda_2020(self):
        """
        Instantiate a set of 'Kuroda 2020' models
        """

        for rot_vel, exponents in [(0, [0]),
                                   (1, [12, 13])]:
            for exponent in exponents:
                model = Kuroda_2020(rotational_velocity=rot_vel * u.rad / u.s,
                                    magnetic_field_exponent=exponent, eos='LS220')

                self.assertEqual(model.metadata['EOS'], 'LS220')
                self.assertEqual(model.metadata['Progenitor mass'], 20*u.Msun)

                # Check that times are in proper units.
                t = model.get_time()
                self.assertTrue(t.unit, u.s)

                # Check that we can compute flux dictionaries.
                f = model.get_initial_spectra(0*u.s, 10*u.MeV)
                self.assertEqual(type(f), dict)
                self.assertEqual(len(f), len(Flavor))
                self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Fornax_2021(self):
        """
        Instantiate a set of 'Fornax 2021' models
        """
        for mass in list(range(12, 24)) + [25, 26, 26.99]:
            model = Fornax_2021(progenitor_mass=mass*u.Msun)

            self.assertEqual(model.metadata['Progenitor mass'], float(mass)*u.Msun)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Zha_2021(self):
        """
        Instantiate a set of 'Zha 2021' models
        """

        for mass in list(range(16, 27)) + [19.89, 22.39, 30, 33]:
            model = Zha_2021(progenitor_mass=mass * u.Msun, eos='STOS_B145')

            self.assertEqual(model.metadata['Progenitor mass'], float(mass)*u.Msun)
            self.assertEqual(model.metadata['EOS'], 'STOS_B145')

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1./(u.erg * u.s))

