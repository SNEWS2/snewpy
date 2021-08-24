#!/usr/bin/env python

from snewpy import snowglobes

SNOwGLoBES_path = "/path/to/snowglobes/"  # where snowglobes is located

# arguments for generate_time_series
model_file = "/path/to/snewpy/models/Nakazato_2013/nakazato-LS220-BH-z0.004-s30.0.fits"
modeltype = 'Nakazato_2013'
transformation = 'AdiabaticMSW_NMO'
d = 10  # Supernova distance in kpc

# Running the modules
outfile = snowglobes.generate_time_series(model_file, modeltype, transformation, d)
snowglobes.simulate(SNOwGLoBES_path, outfile, detector_input="icecube")
snowglobes.collate(SNOwGLoBES_path, outfile, detector_input="icecube")

# An additional, optional argument in both simulate() and collate() is the detector name, if one wants to only run 1 detector, rather than all of them.
# The value of this argument is the same for both functions, except in the case of wc100kt30prct and ar40kt.
# For wc100kt30prct: simulate() takes in "wc100kt30prct", while collate() takes in "wc100kt30prct_eve"
# For ar40kt: simulate() takes in "ar40kt", while collate() takes in "ar40kt_eve"
