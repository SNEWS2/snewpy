#!/usr/bin/env python
from snewpy.models.ccsn import Bollig_2016
from snewpy.neutrino import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.rate_calculator import RateCalculator

import numpy as np
import astropy.units as u

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
