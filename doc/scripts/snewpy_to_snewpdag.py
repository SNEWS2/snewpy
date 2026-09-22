
from snewpy.models.ccsn import Tamborra_2014
from snewpy.neutrino import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.rate_calculator import RateCalculator, center

import numpy as np
import astropy.units as u

model = Tamborra_2014(progenitor_mass=20*u.solMass, direction=1)
transformation = AdiabaticMSW(MixingParameters('NORMAL')) # Desired flavor transformation
       
times    = model.get_time()
energies = np.linspace(0,100,501)<<u.MeV
distance = 10*u.kpc

#Specify sequence of time intervals
window_tstart = 0.001
window_tend = 0.331
window_bins = 330
tbins =  np.linspace(window_tstart,window_tend,window_bins) * u.s
                     
#get the flux from the model
flux = model.get_flux(t=times, E=energies, distance=distance, flavor_xform=transformation)
fluence = flux.integrate('time', limits=tbins)

#set SNOwGLoBES detector to use
detector = "icecube"
rc = RateCalculator()
tables = rc.run(fluence, detector, detector_effects=True)

output_path = "/path/to/output/" #where the output files will be located
fout = open(output_path+"snewpy_output_"+detector+"_Tamborra2014_20solMass,direction=1_1msbin.txt", "a")

tmid = center(tbins)
nevents = np.zeros(len(tmid))
for i in range(len(tmid)):
    nevents[i] = sum([(chan.integrate_or_sum('energy').array.squeeze())[i].value for chan in tables.values()])
    print(i, "\t" , nevents[i], file=fout)
