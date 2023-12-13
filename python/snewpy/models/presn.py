import snewpy.models.presn_loaders as loaders
from snewpy.models.registry_model import RegistryModel
import numpy as np
from scipy.interpolate import interp1d
from astropy import units as u

@RegistryModel(
    progenitor_mass = [15, 25]<<u.Msun
)
class Odrzywolek_2010(loaders.Odrzywolek_2010):
    def __init__(self, progenitor_mass:u.Quantity):
        filename=f"s{progenitor_mass.to_value('Msun'):.0f}_nuebar_data.txt"
        super().__init__(filename)
    
@RegistryModel(
    progenitor_mass = [15, 30]<<u.Msun
)
class Patton_2017(loaders.Patton_2017):
    def __init__(self, progenitor_mass:u.Quantity):
        filename=f"totalLuminosity_{progenitor_mass.to_value('Msun'):.0f}SolarMass.dat"
        super().__init__(filename)
        
@RegistryModel(
    progenitor_mass = [12, 15]<<u.Msun
)
class Kato_2017(loaders.Kato_2017):
    def __init__(self, progenitor_mass:u.Quantity):
        path=f"pre_collapse/m{progenitor_mass.to_value('Msun'):.0f}"
        super().__init__(path)
        
        