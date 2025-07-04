# -*- coding: utf-8 -*-
"""This module implements basic neutrino properties that are used throughout SNEWPY."""

from enum import IntEnum
from astropy import units as u
from dataclasses import dataclass
import numpy as np
from collections.abc import Mapping

from .flavor import ThreeFlavor as Flavor # unused import needed for backward compatibility (see example notebooks)
from .flavor import FlavorScheme, FlavorMatrix


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
    based on `NuFIT 6.0 <http://www.nu-fit.org>`_ results from September 2024,
    published in `JHEP 12 (2024) 216 <https://dx.doi.org/10.1007/JHEP12(2024)216>`_
    [`arXiv:2410.05380 <https://arxiv.org/abs/2410.05380>`_].
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
    
    #define the basis states
    basis_mass = FlavorScheme("ThreeFlavor_MassBasis", start=0,
                                names=['NU_1','NU_2','NU_3','NU_1_BAR','NU_2_BAR','NU_3_BAR'])
    basis_flavor = FlavorScheme("ThreeFlavor_FlavorBasis", start=0,
                                names=['NU_E','NU_MU','NU_TAU','NU_E_BAR','NU_MU_BAR','NU_TAU_BAR'])

    def update(self, **parameters):
        """Update the parameters which are provided in the input arguments. 
        Ignore the input parameters, which are not present in this class.
        """
        keys_matching = {key:val for key,val in parameters.items() if key in self.__dict__}
        self.__dict__.update(**keys_matching)
        
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
        return FlavorMatrix(U, flavor=self.basis_flavor, from_flavor=self.basis_mass)


    def ComplexRotationMatrix(self,i,j,theta,phase):
        """A complex rotation matrix. N.B. the minus sign in the complex exponential matches PDG convention"""
        theta = (theta<<u.radian).value
        phase = (phase<<u.radian).value
        V = np.eye(6,dtype = 'cdouble')
        V[j,j] = V[i,i] = np.cos(theta)
        V[i,j] = np.sin(theta) * np.exp(-1j*phase)
        V[j,i] = -np.conjugate(V[i,j])
        #making the block matrix
        V[3:,3:] = np.conjugate(V[:3,:3])
        return V

    def Pmf_HighDensityLimit(self):
        """ The probability that a given flavor state is a particular matter state in the 
        infinite density limit.  

        Returns
        -------
        6 x 6 matrix
        """
        M2 = FlavorMatrix.zeros(self.basis_mass)
        M2[1,1] = self.dm21_2.value
        M2[2,2] = self.dm31_2.value
        M2[3:,3:]=M2[:3,:3]
        
        U = self.VacuumMixingMatrix()
        HV = U @ M2 @ np.conjugate(U.T)
        
        T = np.real( ( HV['mu','mu'] + HV['tau','tau'] ) / 2 )
        D = np.abs( HV['mu','tau'] )**2 \
           -np.abs( HV['mu','mu']-HV['tau','tau'] )**2 / 4
        Tbar = np.real( ( HV['mu_bar','mu_bar'] + HV['tau_bar','tau_bar'] ) / 2 )
        Dbar = np.abs( HV['mu_bar','tau_bar'] )**2 \
              -np.abs( HV['mu_bar','mu_bar']-HV['tau_bar','tau_bar'] )**2 / 4
    
        PmfHDL = FlavorMatrix.zeros(self.basis_mass, self.basis_flavor)
        
        if self.mass_order == MassHierarchy.NORMAL:
            # The NMO case. Matter state 3 is the electron neutrino, matter state 1bar is the electron antineutrino. 
            k1 = T - np.sqrt(D)
            k2 = T + np.sqrt(D)
            k2bar = Tbar - np.sqrt(Dbar)
            k3bar = Tbar + np.sqrt(Dbar)
            PmfHDL['1','mu'] =  (HV['tau','tau'].real - k1)/(k2-k1)
            PmfHDL['1','tau'] = (HV['mu','mu'].real - k1)/(k2-k1)
            PmfHDL['2','mu'] =  (HV['tau','tau'].real - k2)/(k1-k2)
            PmfHDL['2','tau'] = (HV['mu','mu'].real - k2)/(k1 -k2)
            PmfHDL['3','e'] = 1
        
            PmfHDL['1_bar','e_bar'] = 1
            PmfHDL['2_bar','mu_bar'] =  (HV['tau_bar','tau_bar'].real - k2bar ) / ( k3bar - k2bar )
            PmfHDL['2_bar','tau_bar'] = (HV['mu_bar','mu_bar'].real - k2bar ) / ( k3bar - k2bar )
            PmfHDL['3_bar','mu_bar'] =  (HV['tau_bar','tau_bar'].real - k3bar ) / ( k2bar - k3bar )
            PmfHDL['3_bar','tau_bar'] = (HV['mu_bar','mu_bar'].real - k3bar ) / ( k2bar - k3bar )
        
        elif self.mass_order == MassHierarchy.INVERTED:
            # The IMO case. Matter state 2 is the electron neutrino, matter state 3bar is the electron antineutrino. 
            k1 = T + np.sqrt(D)
            k3 = T - np.sqrt(D)
            k1bar = Tbar - np.sqrt(Dbar)
            k2bar = Tbar + np.sqrt(Dbar)
        
            PmfHDL['1','mu'] = ( HV['tau','tau'].real - k1 ) / ( k3 - k1 )
            PmfHDL['1','tau'] = ( HV['mu','mu'].real - k1 ) / ( k3 - k1 )
            PmfHDL['2','e'] = 1
            PmfHDL['3','mu'] = ( HV['tau','tau'].real - k3) / ( k1 - k3 )
            PmfHDL['3','tau'] = ( HV['mu','mu'].real - k3 ) / ( k1 - k3 )
        
            PmfHDL['1_bar','mu_bar'] = ( HV['tau_bar','tau_bar'].real - k1bar ) / ( k2bar - k1bar )
            PmfHDL['1_bar','tau_bar'] = ( HV['mu_bar','mu_bar'].real - k1bar ) / ( k2bar - k1bar )
            PmfHDL['2_bar','mu_bar'] = ( HV['tau_bar','tau_bar'].real - k2bar ) / ( k1bar - k2bar )
            PmfHDL['2_bar','tau_bar'] = ( HV['mu_bar','mu_bar'].real - k2bar ) / ( k1bar - k2bar )
            PmfHDL['3_bar','e_bar'] = 1
        return PmfHDL

@dataclass
class FourFlavorMixingParameters(ThreeFlavorMixingParameters):
    """A class for four flavor neutrino mixing. 
    ..Note: it is an extension of :class:`ThreeFlavorMixingParameters`, and can be constructed using it:
    
        >>> pars_3f = ThreeFlavorMixingParameters() #standard 3flavor mixing
        >>> pars_4f = FourFlavorMixingParameters(**pars_3f, theta14=90<<u.deg, dm41_2=1<<u.GeV**2)
    """
    #sterile neutrino mixing angles. 
    theta14: u.Quantity[u.deg] = 0<<u.deg
    theta24: u.Quantity[u.deg]|None = 0<<u.deg
    theta34: u.Quantity[u.deg]|None = 0<<u.deg

    #sterile CP violating phases
    delta12: u.Quantity[u.deg]|None = 0<<u.deg
    #delta13: u.Quantity[u.deg] '''same as deltaCP in 3Flavor Mixing'''
    delta24: u.Quantity[u.deg]|None = 0<<u.deg

    #sterile neutrino mass squared differences
    dm41_2: u.Quantity[u.eV**2] = 0<<u.eV**2
    dm42_2: u.Quantity | None = None
    dm43_2: u.Quantity | None = None

    #define the basis states
    basis_mass = FlavorScheme("FourFlavor_MassBasis", start=0,
                                names=['NU_1','NU_2','NU_3','NU_4','NU_1_BAR','NU_2_BAR','NU_3_BAR','NU_4_BAR'])
    basis_flavor = FlavorScheme("FourFlavor_FlavorBasis", start=0,
                                names=['NU_E','NU_MU','NU_TAU','NU_S','NU_E_BAR','NU_MU_BAR','NU_TAU_BAR','NU_S_BAR'])
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
        """The vacuum mixing matrix given the mixing parameters
           N.B. This is a 8 x 8 matrix
        """

        U = self.ComplexRotationMatrix(2,3,self.theta34,0) \
            @ self.ComplexRotationMatrix(1,3,self.theta24,self.delta24) \
            @ self.ComplexRotationMatrix(0,3,self.theta14,0) \
            @ self.ComplexRotationMatrix(1,2,self.theta23,0) \
            @ self.ComplexRotationMatrix(0,2,self.theta13,self.deltaCP) \
            @ self.ComplexRotationMatrix(0,1,self.theta12,self.delta12) 

        return FlavorMatrix(U, flavor=self.basis_flavor, from_flavor=self.basis_mass)


    def ComplexRotationMatrix(self,i,j,theta,phase):
        """A complex rotation matrix. N.B. the minus sign in the complex exponential matches PDG convention"""
        theta = (theta<<u.radian).value
        phase = (phase<<u.radian).value
        V = np.eye(8,dtype = 'cdouble')
        V[j,j] = V[i,i] = np.cos(theta)
        V[i,j] = np.sin(theta) * np.exp(-1j*phase)
        V[j,i] = -np.conjugate(V[i,j])
        #making the block matrix
        V[4:,4:] = np.conjugate(V[:4,:4])
        return V

    def Pmf_HighDensityLimit(self):
        """ The probability that a given flavor state is a particular matter state in the 
        infinite density limit.  This method assumes the 4th matter state is the 
        heaviest and that the electron fraction remains larger than 1/3

        Returns
        -------
        8 x 8 matrix
        """
        M2 = FlavorMatrix.zeros(self.basis_mass)
        M2[1,1] = self.dm21_2.value
        M2[2,2] = self.dm31_2.value
        M2[3,3] = self.dm41_2.value
        M2[4:,4:]=M2[:4,:4]
        
        U = self.VacuumMixingMatrix()
        HV = U @ M2 @ np.conjugate(U.T)
        
        T = np.real( ( HV['mu','mu'] + HV['tau','tau'] ) / 2 )
        D = np.abs(HV['mu','tau'])**2-np.abs( HV['mu','mu']-HV['tau','tau'] )**2 / 4
        Tbar = np.real( ( HV['mu_bar','mu_bar'] + HV['tau_bar','tau_bar'] ) / 2 )
        Dbar = np.abs( HV['mu_bar','tau_bar'] )**2 \
              -np.abs( HV['mu_bar','mu_bar']-HV['tau_bar','tau_bar'] )**2 / 4
        PmfHDL = FlavorMatrix.zeros(self.basis_mass, self.basis_flavor)
        
        if self.mass_order == MassHierarchy.NORMAL:
            # The NMO case. Matter state 4 is the electron neutrino, matter state 1bar is the electron antineutrino.
            # Matter state 3 is the sterile neutrino, matter state 2bar is the sterile antineutrino. 
            k1 = T - np.sqrt(D)
            k2 = T + np.sqrt(D)
            k3bar = Tbar - np.sqrt(Dbar)
            k4bar = Tbar + np.sqrt(Dbar)
            
            PmfHDL['1','mu'] =  (HV['tau','tau'].real - k1)/(k2-k1)
            PmfHDL['1','tau'] = (HV['mu','mu'].real - k1)/(k2-k1)
            PmfHDL['2','mu'] =  (HV['tau','tau'].real - k2)/(k1-k2)
            PmfHDL['2','tau'] = (HV['mu','mu'].real - k2)/(k1-k2)
            PmfHDL['3','S'] = 1
            PmfHDL['4','e'] = 1
        
            PmfHDL['1_bar','e_bar'] = 1
            PmfHDL['2_bar','S_bar'] = 1
            PmfHDL['3_bar','mu_bar'] =  (HV['tau_bar','tau_bar'].real - k3bar ) / ( k4bar - k3bar )
            PmfHDL['3_bar','tau_bar'] = (HV['mu_bar','mu_bar'].real - k3bar ) / ( k4bar - k3bar )
            PmfHDL['4_bar','mu_bar'] =  (HV['tau_bar','tau_bar'].real - k4bar ) / ( k3bar - k4bar )
            PmfHDL['4_bar','tau_bar'] = (HV['mu_bar','mu_bar'].real - k4bar ) / ( k3bar - k4bar )
        
        elif self.mass_order == MassHierarchy.INVERTED:
            # The IMO case. Matter state 4 is the electron neutrino, matter state 3bar is the electron antineutrino.
            # Matter state 2 is the sterile neutrino, matter state 1bar is the sterile antineutrino.
            k1 = T + np.sqrt(D)
            k3 = T - np.sqrt(D)
            k2bar = Tbar - np.sqrt(Dbar)
            k4bar = Tbar + np.sqrt(Dbar)
        
            PmfHDL['1','mu'] = ( HV['tau','tau'].real - k1 ) / ( k3 - k1 )
            PmfHDL['1','tau'] = ( HV['mu','mu'].real - k1 ) / ( k3 - k1 )
            PmfHDL['2','S'] = 1
            PmfHDL['3','mu'] = ( HV['tau','tau'].real - k3) / ( k1 - k3 )
            PmfHDL['3','tau'] = ( HV['mu','mu'].real - k3 ) / ( k1 - k3 )
            PmfHDL['4','e'] = 1
        
            PmfHDL['1_bar','S_bar'] = 1
            PmfHDL['2_bar','mu_bar'] = ( HV['tau_bar','tau_bar'].real - k2bar ) / ( k4bar - k2bar )
            PmfHDL['2_bar','tau_bar'] = ( HV['mu_bar','mu_bar'].real - k2bar ) / ( k4bar - k2bar )
            PmfHDL['3_bar','e_bar'] = 1
            PmfHDL['4_bar','mu_bar'] = ( HV['tau_bar','tau_bar'].real - k4bar ) / ( k2bar - k4bar )
            PmfHDL['4_bar','tau_bar'] = ( HV['mu_bar','mu_bar'].real - k4bar ) / ( k2bar - k4bar )
        return PmfHDL

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
        ThreeFlavorMixingParameters(
            theta12 = 33.65 << u.deg,
            theta13 = 8.51 << u.deg,
            theta23 = 48.33 << u.deg,
            deltaCP = 214 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = 2.455e-3 << u.eV**2
        ),
        MassHierarchy.INVERTED:
        ThreeFlavorMixingParameters(
            theta12 = 33.65 << u.deg,
            theta13 = 8.51 << u.deg,
            theta23 = 48.04 << u.deg,
            deltaCP = 214 << u.deg,
            dm21_2 = 7.53e-5 << u.eV**2,
            dm32_2 = -2.529e-3 << u.eV**2
        )
    }
}
   

def MixingParameters(mass_order:MassHierarchy|str='NORMAL', version:str='NuFIT6.0'):
    if isinstance(mass_order,str):
        mass_order = MassHierarchy[mass_order]
    return parameter_presets[version][mass_order]
