#!/usr/bin/env python

from snewpy.models.ccsn import Nakazato_2013
from snewpy.neutrino import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.rate_calculator import RateCalculator

import numpy as np
import astropy.units as u

model = Nakazato_2013(progenitor_mass=30*u.solMass, revival_time=0*u.ms, metallicity=0.004, eos='LS220')

transformation = AdiabaticMSW(MixingParameters('NORMAL')) # Desired flavor transformation
       
times    = model.get_time()
energies = np.linspace(0,100,501)<<u.MeV
distance = 10*u.kpc

flux = model.get_flux(t=times, E=energies, distance=distance, flavor_xform=transformation)
fluence = flux.integrate('time')

detector = "wc100kt30prct"
rc = RateCalculator()
events = rc.run(fluence, detector, detector_effects=True)

with channel in events:
    filename = f"{model}.{transformation}.{times[0]:.3f}-{times[-1]:.3f},{energies[0]:.3f}-{energies[-1]:.3f},{distance:.3f}.{channel}.npz"
    events[channel].save(filename)

