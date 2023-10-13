# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""
import unittest

from snewpy.neutrino import Flavor
from snewpy.flavor_transformation import NoTransformation
from snewpy.models.ccsn import Nakazato_2013, Tamborra_2014, OConnor_2013, OConnor_2015, \
                          Sukhbold_2015, Bollig_2016, Walk_2018, \
                          Walk_2019, Fornax_2019, Warren_2020, \
                          Kuroda_2020, Fornax_2021, Zha_2021, \
                          Fornax_2022

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
                  Zha_2021,
                  Fornax_2022,
                  ]

        for model in models:
            for pc in model.get_param_combinations():
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
        for mass in [11.2, 20., 27.]:
            for angles in [1,2,3]:
                if mass != 27. and angles > 1:
                    with self.assertRaises(ValueError):
                        model = Tamborra_2014(progenitor_mass=mass*u.Msun, eos='LS220', direction=angles)
                else:
                    model = Tamborra_2014(progenitor_mass=mass*u.Msun, eos='LS220',direction=angles)
                    
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
        for angles in [1,2,3]:
            for rot in ['fast','slow','non']:
                model = Walk_2018(progenitor_mass=mass*u.Msun, direction=angles, rotation=rot)

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

        for mass in [40.,75.]:
            for angles in [1,2]:
                model = Walk_2019(progenitor_mass=mass*u.Msun, direction=angles)

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

    def test_Fornax_2022(self):
        """
        Instantiate a set of 'Fornax 2022' models
        """
        progenitors = [
            '9.0',     '9.25',     '9.5',      '9.75',     '10.0',
            '10.25',    '10.5',     '10.75',    '11.0',     '11.25',
            '11.5',     '11.75',    '12.00.bh', '12.03.bh', '12.07.bh',
            '12.1.bh',  '12.13',    '12.15',    '12.18.bh', '12.20.bh',
            '12.25',    '12.33.bh', '12.40.bh', '12.45.bh', '12.50.bh',
            '12.54.bh', '12.60.bh', '12.63',    '12.70',    '12.72.bh',
            '12.75',    '12.80.bh', '12.85.bh', '12.90.bh', '12.93',
            '12.97.bh', '13.00.bh', '13.05.bh', '13.11',    '13.25.bh',
            '13.27.bh', '13.32.bh', '13.40.bh', '13.45',    '13.50.bh',
            '13.60.bh', '13.75',    '13.82.bh', '13.90.bh', '13.96',
            '14.01',    '14.13.bh', '14.25.bh', '14.40.bh', '14.41.bh',
            '14.43',    '14.44.bh', '14.70.bh', '14.87.bh', '15.00.bh',
            '15.01',    '15.04.bh', '15.05',    '15.38.bh', '16.43',
            '16.65',    '16.99',    '17.00',    '17.07',    '17.10',
            '17.40',    '17.48',    '17.50',    '17.51',    '17.83',
            '18.04',    '18.05',    '18.09',    '18.10',    '18.50',
            '19.02',    '19.56',    '19.83',    '19.99',    '20.08',
            '20.09',    '20.18',    '20.37',    '21.00',    '21.68',
            '22.00',    '22.30',    '22.82',    '23.00',    '23.04',
            '23.43',    '24.00',    '25.00',    '26.00',    '26.99' ]

        for progenitor in progenitors:
            model = Fornax_2022(progenitor=progenitor)

            mass = float(progenitor[:-3] if 'bh' in progenitor else progenitor) * u.Msun
            self.assertEqual(model.metadata['Progenitor'], progenitor)
            self.assertEqual(model.metadata['Progenitor mass'], mass)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

