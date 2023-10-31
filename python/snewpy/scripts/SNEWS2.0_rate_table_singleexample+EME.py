#!/usr/bin/env python
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from snewpy import snowglobes
from snewpy.flavor_transformation import *
from snewpy.neutrino import *

SNOwGLoBES_path = None  # change to SNOwGLoBES directory if using a custom detector configuration
SNEWPY_model_dir = "/path/to/snewpy/models/"  # directory containing model input files

# skycoordinates of neutrino source
Betelgeuse = SkyCoord.from_name('Betelgeuse') 
    
# neutrino detector
SuperK = EarthLocation.of_site('SuperK')

# time when the supernova occured 
# the first time option means the neutrinos traveled through the Earth, the second means they did not
time = Time('2021-5-26 14:14:00')
#time = Time('2021-5-26 14:14:00') - 12*u.hour    

# altaz of supernovae at detector
SNaltaz = Betelgeuse.transform_to(AltAz(obstime=time,location=SuperK)) 

distance = 10  # Supernova distance in kpc
detector = "wc100kt30prct" #SNOwGLoBES detector for water Cerenkov
modeltype = 'Bollig_2016' # Model type from snewpy.models
model = 's11.2c' # Name of model

mix_params = MixingParameters(MassHierarchy.NORMAL)

transformation = AdiabaticMSW(mix_params,SNaltaz) # Desired flavor transformation

# Construct file system path of model file and name of output file
model_path = SNEWPY_model_dir + "/" + modeltype + "/" + model
outfile = modeltype + "_" + model + "_" + str(transformation)

# Now, do the main work:
print("Generating fluence files ...")
tarredfile = snowglobes.generate_fluence(model_path, modeltype, transformation, distance, outfile)

print("Simulating detector effects with SNOwGLoBES ...")
snowglobes.simulate(SNOwGLoBES_path, tarredfile, detector_input=detector)

print("Collating results ...")
tables = snowglobes.collate(SNOwGLoBES_path, tarredfile, skip_plots=True)


# Use results to print the number of events in different interaction channels
key = f"Collated_{outfile}_{detector}_events_smeared_weighted.dat"
total_events = 0
for i, channel in enumerate(tables[key]['header'].split()):
    if i == 0:
        continue
    n_events = sum(tables[key]['data'][i])
    total_events += n_events
    print(f"{channel:10}: {n_events:.3f} events")

#Super-K has 32kT inner volume
print("Total events in Super-K-like detector:",0.32*total_events)
