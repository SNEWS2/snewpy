#!/usr/bin/env python
from snewpy.rate_calculator import RateCalculator
from snewpy.models import ccsn, ccsn_loaders
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.neutrino import MixingParameters

import astropy.units as u
import numpy as np

path = "~/.astropy/cache/snewpy/models/PUSH/"

model = ccsn_loaders.PUSH("luminosity.d","mutau_luminosity.d") # SN model
transformation = AdiabaticMSW(MixingParameters('NORMAL')) # Desired flavor transformation
       
times    = model.get_time()
energies = np.linspace(0,100,501)<<u.MeV
distance = 10*u.kpc

flux = model.get_flux(t=times, E=energies, distance=distance, flavor_xform=transformation)
fluence = flux.integrate('time')

detector = "wc100kt30prct"
rc = RateCalculator()

events = rc.run(fluence, detector, detector_effects=False)        
events_smeared = rc.run(fluence, detector, detector_effects=True)
        
# Compute number of events in all interaction channels
total_events  = sum([chan.integrate_or_sum('energy').array.squeeze().value for chan in events.values()])        
total_events_smeared  = sum([chan.integrate_or_sum('energy').array.squeeze().value for chan in events_smeared.values()])

print("Total events in Super-K-like detector (with smearing):" , 0.32*total_events_smeared)

