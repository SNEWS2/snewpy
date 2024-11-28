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
                          Fornax_2022, Mori_2023, Bugli_2021, Fischer_2020

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
                  Mori_2023,
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

    def test_Bugli_2021(self):
        """
        Instantiate a set of 'Bugli 2021' models
        """
        for Bfield in ['hydro', 'L1', 'L2']:
            for direction in ['average','equator','north','south']:
                if Bfield == 'L2':
                    for grav in ['A','B']:
                        model = Bugli_2021(Bfield=Bfield, grav=grav, rotation=None, direction=direction)
                if Bfield == 'L1':
                    for rotation in [0,90]:
                        model = Bugli_2021(Bfield=Bfield, grav=None, rotation=rotation, direction=direction)
                if Bfield == 'hydro':
                    model = Bugli_2021(Bfield=Bfield, grav=None, rotation=None, direction=direction)

                self.assertEqual(model.metadata['EOS'], 'LS220')
                self.assertEqual(model.metadata['Progenitor mass'], 35*u.Msun)

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

    def test_Fischer_2020(self):
        """
        Instantiate a set of 'Fischer 2020' models
        """
        model = Fischer_2020()
        self.assertEqual(model.metadata['Progenitor mass'], [18 * u.Msun])
        self.assertEqual(model.metadata['EOS'], 'HS(DD2)')

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
        progenitors = Fornax_2022.param['progenitor_mass']
        bh_masses = [12.  , 12.03, 12.07, 12.1 , 12.18, 12.2 , 12.33, 12.4 , 12.45,
                      12.5 , 12.54, 12.6 , 12.72, 12.8 , 12.85, 12.9 , 12.97, 13.  ,
                      13.05, 13.25, 13.27, 13.32, 13.4 , 13.5 , 13.6 , 13.82, 13.9 ,
                      14.13, 14.25, 14.4 , 14.41, 14.44, 14.7 , 14.87, 15.  , 15.04,
                      15.38]<<u.Msun
        # PNS masses for each simulation, D. Vartanyan, private comm.
        pns_masses = [  1.35 , 1.38 , 1.4  , 1.45 , 1.5  ,
                        1.52 , 1.54 , 1.57 , 1.55 , 1.54 ,
                        1.55 , 1.57 , 1.69 , 1.73 , 1.72 ,
                        1.87 , 1.67 , 1.67 , 1.73 , 1.69 ,
                        1.68 , 1.73 , 1.75 , 1.75 , 1.73 ,
                        1.75 , 1.78 , 1.83 , 1.81 , 1.75 ,
                        1.84 , 1.74 , 1.8  , 1.74 , 1.86 ,
                        1.77 , 1.83 , 1.84 , 1.78 , 1.85 ,
                        1.79 , 1.87 , 1.87 , 1.84 , 1.75 ,
                        1.85 , 1.9  , 1.85 , 1.77 , 1.84 ,
                        1.87 , 1.79 , 1.77 , 1.87 , 1.82 ,
                        1.76 , 1.86 , 1.89 , 1.95 , 1.94 ,
                        1.69 , 1.94 , 1.74 , 1.99 , 1.97 ,
                        2.05 , 2.08 , 2.07 , 2.02 , 2.03 ,
                        1.72 , 1.68 , 2.12 , 1.81 , 1.83 ,
                        2.16 , 1.84 , 1.68 , 2.12 , 2.23 ,
                        1.85 , 2.39 , 2.34 , 2.31 , 1.99 ,
                        2.37 , 2.25 , 2.25 , 2.36 , 2.23 ,
                        2.1  , 2.16 , 1.98 , 2.01 , 2.1  ,
                        1.98 , 2.11 , 2.18 , 2.14 , 2.2 ]<<u.Msun

        for progenitor, m_pns in zip(progenitors, pns_masses):
            model = Fornax_2022(progenitor_mass=progenitor)

            self.assertEqual(model.metadata['Progenitor mass'], progenitor)
            self.assertEqual(model.metadata['PNS mass'], m_pns)
            self.assertEqual(model.metadata['Black hole'], progenitor in bh_masses)
            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1/(u.erg * u.s))

    def test_Mori_2023(self):
        """
        Instantiate a set of 'Mori 2023' models
        """
        axion_mass_coupling = [ (0,0),
            (100,2), (100,4), (100,10), (100,12), (100,14), (100,16), (100,20),
            (200,2), (200,4), (200,6), (200,8), (200,10), (200,20)]

        pns_mass = [1.78, 1.77, 1.76, 1.77, 1.77, 1.77, 1.77, 1.74, 1.77, 1.76, 1.75, 1.74, 1.73, 1.62] * u.Msun

        for ((am, ac), mpns) in zip(axion_mass_coupling, pns_mass):
            model = Mori_2023(axion_mass=am*u.MeV, axion_coupling=ac*1e-10/u.GeV)

            axion_mass = float(am) * u.MeV
            self.assertEqual(model.metadata['Axion mass'], axion_mass)

            axion_coupling = float(ac) * 1e-10/u.GeV
            axion_coupling = np.round(axion_coupling.to('1e-10/GeV'))
            self.assertEqual(model.metadata['Axion coupling'], axion_coupling)

            self.assertEqual(model.metadata['Progenitor mass'], 20*u.Msun)
            self.assertEqual(model.metadata['PNS mass'], mpns)

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1./(u.erg * u.s))
