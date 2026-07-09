# -*- coding: utf-8 -*-
"""Integration test based on SNEWS2.0_rate_table_singleexample.py
"""
import unittest

from snewpy import snowglobes
from snewpy.models.ccsn import Bollig_2016

import astropy.units as u

class TestSimpleRate(unittest.TestCase):

    def test_simplerate(self):
        """Integration test based on SNEWS2.0_rate_table_singleexample.py
        """
        SNOwGLoBES_path = None

        distance = 10 * u.kpc # Supernova distance 
        detector = "wc100kt30prct" #SNOwGLoBES detector for water Cerenkov
        model = Bollig_2016(progenitor_mass=11.2<<u.Msun) # SN model
        transformation = 'AdiabaticMSW_NMO' # Desired flavor transformation
       
        # Now, do the main work:
        print("Generating fluence files ...")
        fluences = snowglobes.generate(model, transformation, distance)

        print("Simulating detector effects with SNOwGLoBES ...")
        collated_events = snowglobes.calculate(SNOwGLoBES_path, fluences, detector=detector, detector_effects=True)

        # Use results to print the number of events in different interaction channels
        total_events  = 0
        for channel in collated_events[detector]:
            n_events = collated_events[detector][channel].integrate_or_sum('energy').array.squeeze().value
            print(f"{detector}:{channel}: {n_events:.3f} events")
            total_events += n_events

        #Super-K has 32kT inner volume
        print("Total events in Super-K-like detector (with smearing):" , 0.32*total_events)

        # We do not use the SNOwGLoBES scaling factors but use other constants so we do not
        # expect the results to agree to 7 digits. Here sub-permille agreement is good enough.
        sk_expected = 4065.662374
        sk_computed = 0.32 * total_events
        discrepancy = abs(sk_computed - sk_expected)/sk_expected

        assert discrepancy < 0.001, f"Number of unsmeared events computed for SK is {sk_computed}, should be {sk_expected}"
