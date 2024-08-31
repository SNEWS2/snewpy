# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.presn`` contains models of presupernova neutrino fluxes,
derived from the :class:`SupernovaModel` base class.
"""
import snewpy.models.presn_loaders as loaders
from snewpy.models.registry_model import RegistryModel
from astropy import units as u

@RegistryModel(
    progenitor_mass = [15, 25]<<u.Msun
)
class Odrzywolek_2010(loaders.Odrzywolek_2010):
    """Presupernova model based on
    `A.Odrzywolek, Alexander Heger, Acta Phys.Polon.B 41 (2010) <https://inspirehep.net/literature/870759>`_
    
    Dataset available on `Odrzywolek’s website <http://th.if.uj.edu.pl/~odrzywolek/psns/index.html>`_
    """
    def __init__(self, progenitor_mass:u.Quantity):
        filename=f"s{progenitor_mass.to_value('Msun'):.0f}_nuebar_data.txt"
        super().__init__(filename)
    
@RegistryModel(
    progenitor_mass = [15, 30]<<u.Msun
)
class Patton_2017(loaders.Patton_2017):
    """Presupernova model based on
    `Kelly M. Patton et al. 2017 ApJ 851 6 <https://iopscience.iop.org/article/10.3847/1538-4357/aa95c4>`_
    
    Dataset available on Zenodo (`DOI:10.5281/zenodo.2626645 <http://doi.org/10.5281/zenodo.2626645>`_)
    """
    def __init__(self, progenitor_mass:u.Quantity):
        filename=f"totalLuminosity_{progenitor_mass.to_value('Msun'):.0f}SolarMass.dat"
        super().__init__(filename)
        
@RegistryModel(
    progenitor_mass = [12, 15]<<u.Msun
)
class Kato_2017(loaders.Kato_2017):
    """Presupernova model based on
    `Chinami Kato et al. 2017 ApJ 848 48 <https://iopscience.iop.org/article/10.3847/1538-4357/aa8b72>`_
    
    Dataset available on `Zenodo <https://zenodo.org/records/3768052>`__
    """
    def __init__(self, progenitor_mass:u.Quantity):
        path=f"pre_collapse/m{progenitor_mass.to_value('Msun'):.0f}"
        super().__init__(path)


@RegistryModel(
    progenitor_mass = [12, 15, 20]<<u.Msun
)
class Yoshida_2016(loaders.Yoshida_2016):
    """Presupernova model based on
    `Yoshida et al. (2016), PRD 93, 123012 <https://doi.org/10.1103/PhysRevD.93.123012>`_
    
    Dataset available on `Zenodo <https://zenodo.org/records/3778014>`__
    """
    def __init__(self, progenitor_mass:u.Quantity):
        path=f"t_spc_m{progenitor_mass.to_value('Msun'):.0f}_1.2.txt"
        super().__init__(path)
