# -*- coding: utf-8 -*-
"""Integration test based on SNEWS2.0_rate_table_singleexample.py
"""
import unittest

from snewpy import snowglobes
from snewpy.models.ccsn import Bollig_2016
from snewpy.neutrino import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.rate_calculator import RateCalculator

import numpy as np
import astropy.units as u

class TestSimpleRate(unittest.TestCase):

    def test_simplerate(self):
        """Integration test based on SNEWS2.0_rate_table_singleexample.py
        """
        model = Bollig_2016(progenitor_mass=11.2<<u.Msun) # SN model
        transformation = AdiabaticMSW(MixingParameters('NORMAL')) # Desired flavor transformation
       
        # Now, do the main work:
        print("Generating fluence files ...")
        times    = model.get_time()
        energies = np.linspace(0,100,501)<<u.MeV
        distance = 10*u.kpc
        #get the flux from the model
        flux = model.get_flux(t=times, E=energies, distance=distance, flavor_xform=transformation)
        fluence = flux.integrate('time')

        print("Simulating detector effects ...")
        detector = "wc100kt30prct"
        rc = RateCalculator()
        events = rc.run(fluence, detector, detector_effects=True)
        
        # Compute number of events in all interaction channels
        total_events  = sum([chan.integrate_or_sum('energy').array.squeeze().value for chan in events.values()])

        #Super-K has 32kT inner volume
        print("Total events in Super-K-like detector (with smearing):" , 0.32*total_events)

        # We do not use the SNOwGLoBES scaling factors but use other constants so we do not
        # expect the results to agree to 7 digits. Here sub-permille agreement is good enough.
        sk_expected = 4065.662374
        sk_computed = 0.32 * total_events
        discrepancy = abs(sk_computed - sk_expected)/sk_expected

        assert discrepancy < 0.001, f"Number of unsmeared events computed for SK is {sk_computed}, should be {sk_expected}"

