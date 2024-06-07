# -*- coding: utf-8 -*-
"""This module implements basic neutrino properties that are used throughout SNEWPY."""

from enum import IntEnum
from astropy import units as u
from dataclasses import dataclass
import numpy as np
from collections.abc import Mapping
from .flavor import ThreeFlavor as Flavor
from .flavor import ThreeFlavor, FourFlavor
from .flavor import TwoFlavor as Flavor

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

    def __str__(self):
        if self.value == MassHierarchy.NORMAL: 
           return f'NMO'
        if self.value == MassHierarchy.INVERTED: 
           return f'IMO'

        
@dataclass
class ThreeFlavorMixingParameters(Mapping):
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
    dm32_2: u.Quantity | None = None
    dm31_2: u.Quantity | None = None
    #mass ordering
    mass_order: MassHierarchy | None = None
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


    def VacuumMixingMatrix(self):
        """The vacuum mixing matrix given the mixing paramters
           N.B. This is a 6 x 6 matrix
        """

        U = self.ComplexRotationMatrix(1,2,self.theta23,0) \
            @ self.ComplexRotationMatrix(0,2,self.theta13,self.deltaCP) \
            @ self.ComplexRotationMatrix(0,1,self.theta12,0) 

        """Reorder rows to match SNEWPY's flavor ordering convention"""
        U[ThreeFlavor.NU_E], U[ThreeFlavor.NU_MU], U[ThreeFlavor.NU_TAU], \
            U[ThreeFlavor.NU_E_BAR], U[ThreeFlavor.NU_MU_BAR], U[ThreeFlavor.NU_TAU_BAR] = \
                U[0], U[1], U[2], U[3], U[4], U[5]

        return U


    def ComplexRotationMatrix(self,i,j,theta,phase):
        """A complex rotation matrix. N.B. the minus sign in the complex exponential matches PDG convention"""
        V = np.zeros((6,6),dtype = 'complex_')
        for k in range(6): 
            V[k,k] = 1

        V[i,i] = np.cos(theta) 
        V[j,j] = V[i,i]
        V[i,j] = np.sin(theta) * ( np.cos(phase) - 1j * np.sin(phase) )
        V[j,i] = - np.conjugate(V[i,j])

        V[i+3,i+3] = V[i,i]
        V[j+3,j+3] = V[j,j]
        V[i+3,j+3] = np.conjugate(V[i,j]) 
        V[j+3,i+3] = np.conjugate(V[j,i])

        return V


@dataclass
class FourFlavorMixingParameters(ThreeFlavorMixingParameters):
    """A class for four flavor neutrino mixing. 
    ..Note: it is an extension of :class:`ThreeFlavorMixingParameters`, and can be constructed using it:
    
        >>> pars_3f = ThreeFlavorMixingParameters() #standard 3flavor mixing
        >>> pars_4f = FpourFlavorMixingParameters(**pars_3f, theta14=90<<u.deg, dm41_2=1<<u.GeV**2)
    """
    #sterile neutrino mixing angles. 
    theta14: u.Quantity[u.deg] = 0<<u.deg
    theta24: Optional[u.Quantity[u.deg]] = 0<<u.deg
    theta34: Optional[u.Quantity[u.deg]] = 0<<u.deg

    #sterile CP violating phases
    delta12: Optional[u.Quantity[u.deg]] = 0<<u.deg
    #delta13: u.Quantity[u.deg] '''same as deltaCP in 3Flavor Mixing'''
    delta24: Optional[u.Quantity[u.deg]] = 0<<u.deg

    #sterile neutrino mass squared differences
    dm41_2: u.Quantity[u.eV**2] = 0<<u.eV**2
    dm42_2: u.Quantity | None = None
    dm43_2: u.Quantity | None = None
    
    def __post_init__(self):
        super().__post_init__()
        self.dm42_2 = self.dm42_2 or self.dm41_2 - self.dm21_2
        self.dm43_2 = self.dm43_2 or self.dm41_2 - self.dm31_2

        dm2_sum = self.dm41_2 - self.dm42_2 - self.dm21_2
        assert np.isclose(dm2_sum,0), f'dm41_2 - dm42_2 - dm21_2 = {dm2_sum} !=0'
        dm2_sum = self.dm41_2 - self.dm43_2 - self.dm31_2
        assert np.isclose(dm2_sum,0), f'dm41_2 - dm43_2 - dm31_2 = {dm2_sum} !=0'

    def get_mixing_angles(self):
        """Mixing angles of the PMNS matrix.
        
        Returns
        -------
        tuple
            Angles theta12, theta13, theta23, theta14, theta24, theta34. 
        """
        return (self.theta12, self.theta13, self.theta23, self.theta14, self.theta24, self.theta34)
        
    def get_mass_squared_differences(self):
        """Mass squared differences .
        
        Returns
        -------
        tuple
            dm21_2, dm31_2, dm32_2, dm41_2, dm42_2, dm43_2.
        """        
        return (self.dm21_2, self.dm31_2, self.dm32_2, self.dm41_2, self.dm42_2, self.dm43_2)

    def VacuumMixingMatrix(self):
        """The vacuum mixing matrix given the mixing paramters
           N.B. This is a 8 x 8 matrix
        """

        U = self.ComplexRotationMatrix(2,3,self.theta34,0) \
            @ self.ComplexRotationMatrix(1,3,self.theta24,self.delta24) \
            @ self.ComplexRotationMatrix(0,3,self.theta14,0) \
            @ self.ComplexRotationMatrix(1,2,self.theta23,0) \
            @ self.ComplexRotationMatrix(0,2,self.theta13,self.deltaCP) \
            @ self.ComplexRotationMatrix(0,1,self.theta12,self.delta12) 

        """Reorder rows to match SNEWPY's flavor ordering convention"""
        U[FourFlavor.NU_E], U[FourFlavor.NU_MU], U[FourFlavor.NU_TAU], U[FourFlavor.NU_S], \
            U[FourFlavor.NU_E_BAR], U[FourFlavor.NU_MU_BAR], U[FourFlavor.NU_TAU_BAR], U[FourFlavor.NU_S_BAR] = \
                U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7]

        return U


    def ComplexRotationMatrix(self,i,j,theta,phase):
        """A complex rotation matrix. N.B. the minus sign in the complex exponential matches PDG convention"""
        V = np.zeros((8,8),dtype = 'complex_')
        for k in range(8): 
            V[k,k] = 1

        V[i,i] = np.cos(theta) 
        V[j,j] = V[i,i]
        V[i,j] = np.sin(theta) * ( np.cos(phase) - 1j * np.sin(phase) )
        V[j,i] = - np.conjugate(V[i,j])

        V[i+4,i+4] = V[i,i]
        V[j+4,j+4] = V[j,j]
        V[i+4,j+4] = np.conjugate(V[i,j]) 
        V[j+4,i+4] = np.conjugate(V[j,i])

        return V


parameter_presets = {
    'NuFIT5.0': {
        # Values from http://www.nu-fit.org/?q=node/228; cite as JHEP 09 (2020) 178 [arXiv:2007.14792]
        MassHierarchy.NORMAL:
        ThreeFlavorMixingParameters(
            theta12 = 33.44<< u.deg,
            theta13 = 8.57 << u.deg,
            theta23 = 49.2 << u.deg,
            deltaCP = 197 << u.deg,
            dm21_2 = 7.42e-5 << u.eV**2,
            dm31_2 = 2.517e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        ThreeFlavorMixingParameters(    
            theta12 = 33.45 << u.deg,
            theta13 = 8.60 << u.deg,
            theta23 = 49.3 << u.deg,
            deltaCP = 282 << u.deg,
            dm21_2 = 7.42e-5 << u.eV**2,
            dm32_2 = -2.498e-3 << u.eV**2
        )
    },
    'NuFIT5.2': {
        # Values from http://www.nu-fit.org/?q=node/256; cite as JHEP 09 (2020) 178 [arXiv:2007.14792]
        MassHierarchy.NORMAL:
        ThreeFlavorMixingParameters(
            theta12 = 33.41 << u.deg,
            theta13 = 8.58 << u.deg,
            theta23 = 42.2 << u.deg,
            deltaCP = 232 << u.deg,
            dm21_2 = 7.41e-5 << u.eV**2,
            dm31_2 = 2.507e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        ThreeFlavorMixingParameters(        
            theta12 = 33.41 << u.deg,
            theta13 = 8.57 << u.deg,
            theta23 = 49.0 << u.deg,
            deltaCP = 276 << u.deg,
            dm21_2 = 7.41e-5 << u.eV**2,
            dm32_2 = -2.486e-3 << u.eV**2
        )
    },
    'NuFIT6.0': {
        # Values from http://www.nu-fit.org/?q=node/294; cite as arXiv:2410.05380
        MassHierarchy.NORMAL:
        ThreeFlavorMixingParameters(
            theta12 = 33.68 << u.deg,
            theta13 = 8.56 << u.deg,
            theta23 = 43.3 << u.deg,
            deltaCP = 212 << u.deg,
            dm21_2 = 7.49e-5 << u.eV**2,
            dm31_2 = 2.513e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        ThreeFlavorMixingParameters(
            theta12 = 33.68 << u.deg,
            theta13 = 8.59 << u.deg,
            theta23 = 47.9 << u.deg,
            deltaCP = 274 << u.deg,
            dm21_2 = 7.49e-5 << u.eV**2,
            dm32_2 = -2.484e-3 << u.eV**2
        )
    },
    'PDG2022':{
        # Cite as R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
        MassHierarchy.NORMAL:
        ThreeFlavorMixingParameters(
            theta12 = 33.65 << u.deg,
            theta13 = 8.53 << u.deg,
            theta23 = 47.64 << u.deg,
            deltaCP = 245 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = 2.453e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        ThreeFlavorMixingParameters(
            theta12 = 33.65 << u.deg,
            theta13 = 8.53 << u.deg,
            theta23 = 47.24 << u.deg,
            deltaCP = 245 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = -2.536e-3 << u.eV**2
        )
    },
    'PDG2024':{
        # Values from https://pdglive.lbl.gov/Particle.action?node=S067&init=0
        # Cite as S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024)
        MassHierarchy.NORMAL:
        MixingParameters3Flavor(
            theta12 = 33.65 << u.deg,
            theta13 = 8.51 << u.deg,
            theta23 = 48.33 << u.deg,
            deltaCP = 214 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = 2.455e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        MixingParameters3Flavor(
            theta12 = 33.65 << u.deg,
            theta13 = 8.51 << u.deg,
            theta23 = 48.04 << u.deg,
            deltaCP = 214 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = -2.529e-3 << u.eV**2
        )
    }
}
   

def MixingParameters(mass_order:MassHierarchy=MassHierarchy.NORMAL, version:str='NuFIT5.0'):
    return parameter_presets[version][mass_order]
