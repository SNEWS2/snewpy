import re
import tarfile
import io
import json

import numpy as np
from astropy import units as u

from snewpy.models.ccsn import Bollig_2016
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.neutrino import MixingParameters
from snewpy.flavor import ThreeFlavor

def marley_name(flavor):
    marley_name = "v"
    if 'E' in flavor.name:
        marley_name += "e"
    elif 'MU' in flavor.name:
        marley_name += "u"
    else:
        marley_name += "t"
    if 'BAR' in flavor.name:
        marley_name += "bar"

    return marley_name  
    
def save_as_marley(fluence, output_filename, additional_inputs=None):
    """Save the contents of a fluence in MARLEY format.

       For each time bin, 6 Marley config files are generated and placed into a tarfile.

       Parameters
       ----------
       fluence : flux container
           Time-integrated neutrino flux from supernova model.
       additional_inputs : dict
           Additional dictionary that will be inserted into the Marley config files. 

       Returns
       -------
       str
    	   Path of tarfile with the Marley config files.
       """

    times=fluence.time
    ntbins = len(times)

    if additional_inputs != None:
        marley_data = additional_inputs
    else: 
        marley_data = {}
        
    if 'seed' not in marley_data:        
        marley_data['seed'] = 123456
       
    if 'direction' not in marley_data:                
        marley_data['direction'] = { 'x': 0.0, 'y': 0.0, 'z': 1.0 }
        
    if 'target' not in marley_data:                
        marley_data['target'] = { 'nuclides': [1000180400], 'atom_fractions': [1.0] }

    if 'reactions' not in marley_data:
        marley_data['reactions'] = [ "ve40ArCC_Bhattacharya2009.react", "ES.react" ] 
       
    marley_data["source"] = { 'type': "histogram", 'E_bin_lefts': fluence.energy[:-1].value.tolist(), 'Emax': fluence.energy[-1].value } 
    
    with tarfile.open(output_filename, 'w:bz2') as tf:
        for i in range(ntbins):
            for f in ThreeFlavor:
                marley_data['source'].update({"neutrino": marley_name(f), "weights": np.squeeze(fluence[f,i].array).value.tolist() })

                marley_output_name = marley_name(f) + f'.{times[i].value:.3f}' + ".hepevt"                
                marley_data['executable_settings'] = { 
                    'events': 1000, 
                    'output': [ { 'file': marley_output_name, 'format': "hepevt", 'mode': "overwrite" } ]
                    }                

                json_str = json.dumps(marley_data,indent=4)
                json_str = re.sub(r'"(\w+)":', r'\1:',json_str) # remove the quotes around the keys in the json_str                
                json_bytes = json_str.encode('ascii')                  

                name_in_tar = marley_name(f) + f'.{times[i]:.3f}' + ".js"
                tar_info = tarfile.TarInfo(name=name_in_tar)
                tar_info.size = len(json_bytes)           
               
                tf.addfile(tar_info, io.BytesIO(json_bytes))

    return output_filename


snmodel = Bollig_2016(progenitor_mass=27<<u.Msun)

times    = snmodel.time
energies = np.linspace(0,50,501)<<u.MeV

d = 10 * u.kpc  # distance to SN
    
mp_nmo = MixingParameters()
xform = AdiabaticMSW(mp_nmo)

flux = snmodel.get_flux(t=times, E=energies,  distance=d, flavor_xform=xform)
fluence = flux.integrate('time', limits = times).integrate('energy', limits = energies)

output_filename = "Bollig_2016.MARLEY.tar.bz2"
save_as_marley(fluence,output_filename)

