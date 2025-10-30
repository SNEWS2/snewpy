# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.pisn`` contains models of pair-instability supernova neutrino fluxes,
"""
from snewpy.models.base import SupernovaModel
from snewpy.models.registry_model import RegistryModel
from astropy import units as u

@RegistryModel(
    progenitor_mass = [150, 250] * u.Msun,
    eos = ['SFHo', 'Helm'],    
)
class Wright_2017(SupernovaModel):
    """PISN model described in the paper `Neutrino signal from pair-instability supernovae' by Wright et al.,  
    Phys. Rev. D96 (2017) 103008  <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.103008>`_
    """
    def _metadata_from_filename(self, filename:str):
        metadata = {
            'Progenitor mass': float(filename.split('_')[1].strip('Msun')) * u.Msun,
            'EOS': filename.split('_')[2].strip('EOS=')
        }
        return metadata    
        
    def __init__(self, progenitor_mass:u.Quantity, eos:str):
        filename=f"PISN_{progenitor_mass.to_value('Msun'):.0f}Msun_EOS={eos}_NeutrinoFlux.tar.bz2"
        super().__init__(filename, self.metadata)
