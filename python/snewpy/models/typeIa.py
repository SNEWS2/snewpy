# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.typeIa`` contains models of neutrino fluxes from Type Ia supernovae
"""
import logging
import os

from astropy import units as u

from snewpy.models import typeIa_loaders as loaders
from .base import SupernovaModel

from snewpy.models.registry_model import RegistryModel, Parameter
from snewpy.models.registry_model import all_models
        
@RegistryModel(
    mechanism = ['GCD', 'DDT'],    
)
class TypeIa(loaders.TypeIa):
    """The DDT model is described in 'Neutrinos from type Ia supernovae: The deflagration-to-detonation transition scenario', by Warren P. Wright et al.,
    [Phys. Rev. D94 (2016) 025026](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.025026), [arXiv:1605.01408](https://arxiv.org/abs/1605.01408).  
    There are 30 snapshots in time and the format of each data file is the SNOwGLoBES format.  

    The GCD model is described in 'Neutrinos from type Ia supernovae: The gravitationally confined detonation scenario' by Warren P. Wright et al.,
    [Phys. Rev. D95 (2017) 043006](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.043006), [arXiv:1609.07403](https://arxiv.org/abs/1609.07403).  
    There are 64 snapshots in time and the format of each data file is the SNOwGLoBES format.
    """
    def _metadata_from_filename(self, filename:str):
        metadata = {
            'mechanism': filename.split('_')[0]
        }
        return metadata    
        
    def __init__(self, mechanism:str):
        filename = f"{mechanism}_NeutrinoFlux.tar.bz2"
        super().__init__(filename)        
    

