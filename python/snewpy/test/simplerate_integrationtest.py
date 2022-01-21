# -*- coding: utf-8 -*-
"""Integration test based on SNEWS2.0_rate_table_singleexample.py
"""
import unittest
import os
from snewpy import snowglobes


class TestSimpleRate(unittest.TestCase):

    def test_simplerate(self):
        """Integration test based on SNEWS2.0_rate_table_singleexample.py
        """
        # Hardcoded paths on GitHub Action runner machines
        SNOwGLoBES_path = os.environ['SNOWGLOBES']
        SNEWPY_model_dir = "models/"

        distance = 10  # Supernova distance in kpc
        detector = "wc100kt30prct" #SNOwGLoBES detector for water Cerenkov
        modeltype = 'Bollig_2016' # Model type from snewpy.models
        model = 's11.2c' # Name of model
        transformation = 'AdiabaticMSW_NMO' # Desired flavor transformation

        # Construct file system path of model file and name of output file
        model_path = SNEWPY_model_dir + "/" + modeltype + "/" + model
        outfile = modeltype + "_" + model + "_" + transformation

        # Now, do the main work:
        print("Generating fluence files ...")
        tarredfile = snowglobes.generate_fluence(model_path, modeltype, transformation, distance, outfile)

        print("Simulating detector effects with SNOwGLoBES ...")
        snowglobes.simulate(SNOwGLoBES_path, tarredfile, detector_input=detector, detector_effects=False)

        print("Collating results ...")
        tables = snowglobes.collate(SNOwGLoBES_path, tarredfile, detector_input=detector, skip_plots=True, smearing=False)

        # Use results to print the number of events in different interaction channels
        key = f"Collated_{outfile}_{detector}_events_unsmeared_weighted.dat"
        total_events = 0
        for i, channel in enumerate(tables[key]['header'].split()):
            if i == 0:
                continue
            n_events = sum(tables[key]['data'][i])
            total_events += n_events
            print(f"{channel:10}: {n_events:.3f} events")

        #Super-K has 32kT inner volume
        print("Total events in Super-K-like detector:" , 0.32*total_events)

        # We do not use the SNOwGLoBES scaling factors but use other constants so we do not
        # expect the results to agree to 7 digits. Here sub-permille agreement is good enough.
        sk_expected = 4486.929197175579
        sk_computed = 0.32 * total_events
        discrepancy = abs(sk_computed - sk_expected)/sk_expected

        assert discrepancy < 0.001, f"Number of events computed for SK is {sk_computed}, should be {sk_expected}"
