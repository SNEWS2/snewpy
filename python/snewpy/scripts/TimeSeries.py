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
snowglobes.collate(SNOwGLoBES_path, outfile)

# An additional, optional argument in simulate() is the detector name, if one wants to only run 1 detector, rather than all of them.
