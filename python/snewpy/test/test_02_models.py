# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""
import unittest

from snewpy.neutrino import Flavor
from snewpy.models.loaders import Nakazato_2013, Tamborra_2014, OConnor_2013, OConnor_2015, \
                                  Sukhbold_2015, Bollig_2016, Walk_2018, \
                                  Walk_2019, Fornax_2019, Warren_2020, \
                                  Kuroda_2020, Fornax_2021, Zha_2021
from snewpy._model_urls import model_urls
from astropy import units as u
from snewpy import model_path
import os


class TestModels(unittest.TestCase):

    def test_model_urls(self):
        """Test that snewpy._model_urls.model_urls is empty. This should be populated if snewpy is downloaded from PyPI.
        This serves as a guard against accidentally committing/merging a populated model_urls to main.
        """
        self.assertFalse(model_urls)

    def test_Nakazato_2013(self):
        """
        Instantiate a set of 'Nakazato 2013' models
        """
        for z in [0.004, 0.02]:
            for trev in [100, 200, 300]:
                for mass in [13., 20., 50.]:
                    mfile = 'Nakazato_2013/nakazato-shen-z{}-t_rev{}ms-s{:.1f}.fits'.format(z, trev, mass)
                    model = Nakazato_2013(os.path.join(model_path, mfile))

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
            mfile = 'Tamborra_2014/s{:.1f}c_3D_dir1'.format(mass)
            model = Tamborra_2014(os.path.join(model_path, mfile), eos='LS220')

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
                metadata = {
                    'Progenitor mass': mass * u.Msun,
                    'EOS': eos,
                }
                model = OConnor_2013(filename=f'{eos}_timeseries.tar.gz', metadata=metadata)

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
        mfile = 'M1_neutrinos.dat'
        model = OConnor_2015(mfile)

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
                mfile = 'Sukhbold_2015/sukhbold-{}-{}.fits'.format(eos, mass)
                model = Sukhbold_2015(os.path.join(model_path, mfile))

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
            mfile = 'Bollig_2016/s{:.1f}c'.format(mass)
            model = Bollig_2016(os.path.join(model_path, mfile), eos='LS220')

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
        mfile = 'Walk_2018/s{:.1f}c_3D_nonrot_dir1'.format(mass)
        model = Walk_2018(os.path.join(model_path, mfile), eos='LS220')

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
        mass = 40.
        mfile = 'Walk_2019/s{:.1f}c_3DBH_dir1'.format(mass)
        model = Walk_2019(os.path.join(model_path, mfile), eos='LS220')

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
            mfile = 'Fornax_2019/lum_spec_{}M.h5'.format(mass)
            model = Fornax_2019(os.path.join(model_path, mfile), metadata={'Progenitor mass': mass*u.Msun})

            self.assertEqual(model.metadata['Progenitor mass'], mass*u.Msun)

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
        masses = [
            '10.0', '10.25', '10.5', '10.75', '100', '11.0', '11.25', '11.5', '11.75',
            '12.0', '12.25', '12.5', '12.75', '120', '13.0', '13.1', '13.2', '13.3',
            '13.4', '13.5', '13.6', '13.7', '13.8', '13.9', '14.0', '14.1', '14.2',
            '14.3', '14.4', '14.5', '14.6', '14.7', '14.8', '14.9', '15.0', '15.1',
            '15.2', '15.3', '15.4', '15.5', '15.6', '15.7', '15.8', '15.9', '16.0',
            '16.1', '16.2', '16.3', '16.4', '16.5', '16.6', '16.7', '16.8', '16.9',
            '17.0', '17.1', '17.2', '17.3', '17.4', '17.5', '17.6', '17.7', '17.8',
            '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7',
            '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6',
            '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5',
            '20.6', '20.7', '20.8', '20.9', '21.0', '21.1', '21.2', '21.3', '21.4',
            '21.5', '21.6', '21.7', '21.8', '21.9', '22.0', '22.1', '22.2', '22.3',
            '22.4', '22.5', '22.6', '22.7', '22.8', '22.9', '23.0', '23.1', '23.2',
            '23.3', '23.4', '23.5', '23.6', '23.7', '23.8', '23.9', '24.0', '24.1',
            '24.2', '24.3', '24.4', '24.5', '24.6', '24.7', '24.8', '24.9', '25.0',
            '25.1', '25.2', '25.3', '25.4', '25.5', '25.6', '25.7', '25.8', '25.9',
            '26.0', '26.1', '26.2', '26.3', '26.4', '26.5', '26.6', '26.7', '26.8',
            '26.9', '27.0', '27.1', '27.2', '27.3', '27.4', '27.5', '27.6', '27.7',
            '27.8', '27.9', '28.0', '28.1', '28.2', '28.3', '28.4', '28.5', '28.6',
            '28.7', '28.8', '28.9', '29.0', '29.1', '29.2', '29.3', '29.4', '29.5',
            '29.6', '29.7', '29.8', '29.9', '30.0', '31', '32', '33', '35', '40', '45',
            '50', '55', '60', '70', '80', '9.0', '9.25', '9.5', '9.75']

        for mixing in [1.23, 1.25, 1.27]:
            for mass in masses:
                mfile = f'stir_a{mixing}/stir_multimessenger_a{mixing}_m{mass}.h5'
                metadata = {
                    'Progenitor mass': float(mass) * u.Msun,
                    'Turb. mixing param.': mixing,
                    'EOS': 'SFHo',
                }
                model = Warren_2020(mfile, metadata)

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
        for field in ['R00B00', 'R10B12', 'R10B13']:
            mfile = 'Kuroda_2020/Lnu{}.dat'.format(field)
            model = Kuroda_2020(os.path.join(model_path, mfile))

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
        for mass in ['12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '25', '26', '26.99']:
            mfile = 'Fornax_2021/lum_spec_{}M_r10000_dat.h5'.format(mass)
            model = Fornax_2021(os.path.join(model_path, mfile), metadata={'Progenitor mass': mass*u.Msun})

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
        for mass in ['16', '17', '18', '19', '19.89', '20', '21', '22.39', '23', '24', '25', '26', '30', '33']:
            mfile = 'Zha_2021/s{}.dat'.format(mass)
            model = Zha_2021(os.path.join(model_path, mfile))

            # Check that times are in proper units.
            t = model.get_time()
            self.assertTrue(t.unit, u.s)

            # Check that we can compute flux dictionaries.
            f = model.get_initial_spectra(0*u.s, 10*u.MeV)
            self.assertEqual(type(f), dict)
            self.assertEqual(len(f), len(Flavor))
            self.assertEqual(f[Flavor.NU_E].unit, 1./(u.erg * u.s))

