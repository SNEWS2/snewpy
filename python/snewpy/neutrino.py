# -*- coding: utf-8 -*-
"""This module implements basic neutrino properties that are used throughout SNEWPY."""

from enum import IntEnum
from astropy import units as u
from dataclasses import dataclass
from typing import Optional
import numpy as np
from collections.abc import Mapping

class MassHierarchy(IntEnum):
    """Neutrino mass ordering: ``NORMAL`` or ``INVERTED``."""
    NORMAL = 1
    INVERTED = 2
    
    @classmethod
    def derive_from_dm2(cls, dm12_2, dm32_2, dm31_2):
        """derive the mass hierarchy based on the given mass squared differences"""
        assert dm12_2>0,f'dm12_2(dm12_2) should be positive'
        assert (dm32_2*dm31_2>=0),f'dm32_2 ({dm32_2}) and dm31_2 ({dm31_2}) should be of the same sign'
        if(dm32_2>=0):
            return MassHierarchy.NORMAL
        else:
            return MassHierarchy.INVERTED
            
class Flavor(IntEnum):
    """Enumeration of CCSN Neutrino flavors."""
    NU_E = 0
    NU_X = 1
    NU_E_BAR = 2
    NU_X_BAR = 3
    
    def to_tex(self):
        """LaTeX-compatible string representations of flavor."""
        if '_BAR' in self.name:
            return r'$\overline{{\nu}}_{0}$'.format(self.name[3].lower())
        return r'$\{0}$'.format(self.name.lower())

    @property
    def is_electron(self):
        """Return ``True`` for ``Flavor.NU_E`` and ``Flavor.NU_E_BAR``."""
        return self.value in (Flavor.NU_E.value, Flavor.NU_E_BAR.value)

    @property
    def is_neutrino(self):
        """Return ``True`` for ``Flavor.NU_E`` and ``Flavor.NU_X``."""
        return self.value in (Flavor.NU_E.value, Flavor.NU_X.value)

    @property
    def is_antineutrino(self):
        return self.value in (Flavor.NU_E_BAR.value, Flavor.NU_X_BAR.value)
        
@dataclass
class MixingParameters3Flavor(Mapping):
    """Mixing angles and mass differences, assuming three neutrino flavors.
    This class contains the default values used throughout SNEWPY, which are
    based on `NuFIT 5.0 <http://www.nu-fit.org>`_ results from July 2020,
    published in `JHEP 09 (2020) 178 <https://dx.doi.org/10.1007/JHEP09(2020)178>`_
    [`arXiv:2007.14792 <https://arxiv.org/abs/2007.14792>`_].
    """
    #mixing angles
    theta12: u.Quantity[u.deg]
    theta13: u.Quantity[u.deg]
    theta23: u.Quantity[u.deg]
    #CP violation phase         
    deltaCP: u.Quantity[u.deg]
    #square mass difference
    dm21_2: u.Quantity[u.eV**2]
    dm32_2: Optional[u.Quantity] = None
    dm31_2: Optional[u.Quantity] = None
    #mass ordering
    mass_order: Optional[MassHierarchy] = None
    # Note: in IH, the mass splittings are: m3..............m1..m2.

    def __iter__(self):
        return iter(self.__dict__)
    def __getitem__(self, item):
        return self.__dict__[item]
    def __len__(self):
        return len(self.__dict__)
    def _repr_markdown_(self):
        s = [f'**{self.__class__.__name__}**']
        s+=  ['|Parameter|Value|',
              '|:--------|:----:|']
        for name, v in self.__dict__.items():
            try: 
                s += [f"|{name}|{v._repr_latex_()}"]
            except:
                try: 
                    s += [f"|{name}|{v.name}"]
                except:
                    s += [f"|{name}|{v}|"]
        return '\n'.join(s)
        
    def __post_init__(self):
        #calculate the missing dm2
        if self.dm31_2 is None: 
            self.dm31_2 = self.dm32_2+self.dm21_2
        if self.dm32_2 is None: 
            self.dm32_2 = self.dm31_2-self.dm21_2  
        #evaluate mass ordering
        if self.mass_order is None:
            self.mass_order = MassHierarchy.derive_from_dm2(*self.get_mass_squared_differences())
        #validate angles
        #angles_sum = sum(self.get_mixing_angles())
        #assert angles_sum==90<<u.deg, f'Mixing angles sum is {angles_sum}!=90 degrees'
        #check hierarchy
        if self.mass_order==MassHierarchy.NORMAL:
            assert self.dm32_2>0, 'dm32_2 should be positive for NH'
            assert self.dm31_2>0, 'dm31_2 should be positive for NH'
        else:
            assert self.dm32_2<0, 'dm32_2 should be negative for IH'
            assert self.dm31_2<0, 'dm31_2 should be negative for IH'
        #validate dm2
        dm2_sum = self.dm32_2+self.dm21_2-self.dm31_2
        assert np.isclose(dm2_sum,0), f'dm32_2+dm31_2-dm31_2 = {dm2_sum} !=0'

    def get_mixing_angles(self):
        """Mixing angles of the PMNS matrix.
        
        Returns
        -------
        tuple
            Angles theta12, theta13, theta23.
        """
        return (self.theta12, self.theta13, self.theta23)
        
    def get_mass_squared_differences(self):
        """Mass squared differences .
        
        Returns
        -------
        tuple
            dm21_2, dm31_2, dm32_2.
        """        
        return (self.dm21_2, self.dm31_2, self.dm32_2)

@dataclass
class MixingParameters4Flavor(MixingParameters3Flavor):
    """A class for four flavor neutrino mixing. 
    ..Note: it is an extension of :class:`MixingParameters3Flavor`, and can be constructed using it:
    
        >>> pars_3f = MixingParameters() #standard 3flavor mixing
        >>> pars_4f = MixingParameters4Flavor(**pars_3f, theta14=90<<u.deg, dm41_2=1<<u.GeV**2)
    """
    #sterile neutrino miging angles
    theta14: u.Quantity[u.deg] = 0<<u.deg
    theta24: u.Quantity[u.deg] = 0<<u.deg
    theta34: u.Quantity[u.deg] = 0<<u.deg
    #sterile neutrino mass squared differences
    dm41_2: u.Quantity[u.eV**2] = 0<<u.eV**2
    dm42_2: Optional[u.Quantity] = None
    dm43_2: Optional[u.Quantity] = None
    
    def __post_init__(self):
        super().__post_init__()
        self.dm42_2 = self.dm42_2 or self.dm41_2 - self.dm21_2
        self.dm43_2 = self.dm43_2 or self.dm41_2 - self.dm31_2

        dm2_sum = self.dm41_2 - self.dm42_2 - self.dm21_2
        assert np.isclose(dm2_sum,0), f'dm41_2 - dm42_2 - dm21_2 = {dm2_sum} !=0'
        dm2_sum = self.dm41_2 - self.dm43_2 - self.dm31_2
        assert np.isclose(dm2_sum,0), f'dm41_2 - dm43_2 - dm31_2 = {dm2_sum} !=0'

parameter_presets = {
    # Values from JHEP 09 (2020) 178 [arXiv:2007.14792] and www.nu-fit.org.
    'NuFIT5.0': {
        MassHierarchy.NORMAL:
        MixingParameters3Flavor(
            theta12 = 33.44<< u.deg,
            theta13 = 8.57 << u.deg,
            theta23 = 49.20 << u.deg,
            deltaCP = 197  << u.deg,
            dm21_2 =  7.42e-5  << u.eV**2,
            dm31_2 =  2.517e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        MixingParameters3Flavor(    
            theta12 = 33.45 << u.deg,
            theta13 = 8.60 << u.deg,
            theta23 = 49.30 << u.deg,
            deltaCP = 282 << u.deg,
            dm21_2 = 7.42e-5 << u.eV**2,
            dm32_2 = -2.498e-3 << u.eV**2
        )
    },
    'NuFIT5.2': {
        MassHierarchy.NORMAL:
        MixingParameters3Flavor(
            theta12 = 33.41 << u.deg,
            theta13 = 8.58 << u.deg,
            theta23 = 42.20 << u.deg,
            deltaCP = 232 << u.deg,
            dm21_2 = 7.41e-5 << u.eV**2,
            dm31_2 = 2.507e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        MixingParameters3Flavor(        
            theta12 = 33.41 << u.deg,
            theta13 = 8.57 << u.deg,
            theta23 = 49.00 << u.deg,
            deltaCP = 276 << u.deg,
            dm21_2 = 7.41e-5 << u.eV**2,
            dm32_2 = -2.486e-3 << u.eV**2
        )
    },
    'PDG2022':{
    # Values from https://pdg.lbl.gov
        MassHierarchy.NORMAL:
        MixingParameters3Flavor(
            theta12 = 33.65 << u.deg,
            theta13 = 8.53 << u.deg,
            theta23 = 47.64 << u.deg,
            deltaCP = 245 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = 2.453e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        MixingParameters3Flavor(
            theta12 = 33.65 << u.deg,
            theta13 = 8.53 << u.deg,
            theta23 = 47.24 << u.deg,
            deltaCP = 245 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = -2.536e-3 << u.eV**2
        )
    }
}
   

def MixingParameters(mh:MassHierarchy=MassHierarchy.NORMAL, version:str='NuFIT5.0'):
    return parameter_presets[version][mh]
