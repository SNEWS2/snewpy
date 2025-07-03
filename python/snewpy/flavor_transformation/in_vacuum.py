r"""
Transformations in vacuum
=========================

Transitions of neutrino mass states :math:`\nu_i\to\nu_j` in vacuum
"""
from abc import abstractmethod, ABC
import numpy as np
from astropy import units as u
from astropy import constants as c

from snewpy.flavor  import FlavorMatrix
from .base import ThreeFlavorTransformation, FourFlavorTransformation
from snewpy.neutrino import MassHierarchy

###############################################################################
# Vacuum transformations
###############################################################################
class VacuumTransformation(ABC):
    @abstractmethod
    def P_mm(self, t, E)->FlavorMatrix:
        pass
###############################################################################
class NoVacuumTransformation(VacuumTransformation):
    def P_mm(self, t, E)->FlavorMatrix:
        return FlavorMatrix.eye(self.mixing_params.basis_mass)
        
class NeutrinoDecay(VacuumTransformation, ThreeFlavorTransformation):
    """Decay effect, where the heaviest neutrino decays to the lightest
    neutrino. For a description and typical parameters, see A. de Gouvêa et al.,
    PRD 101:043013, 2020, arXiv:1910.01127.
    """
    def __init__(self, mass=1*u.eV/c.c**2, tau=1*u.day, dist=10*u.kpc):
        """Initialize transformation matrix.

        Parameters
        ----------
        mass : astropy.units.quantity.Quantity
            Mass of the heaviest neutrino; expect in eV/c^2.
        tau : astropy.units.quantity.Quantity
            Lifetime of the heaviest neutrino.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        """
        super().__init__()
        self.m = mass
        self.tau = tau
        self.d = dist
        
    def gamma(self, E):
        """Decay width of the heaviest neutrino mass.

        Parameters
        ----------
        E : float
            Energy of the nu3.

        Returns
        -------
        Gamma : float
            Decay width of the neutrino mass, in units of 1/length.

        :meta private:
        """
        return self.m*c.c / (E*self.tau)

    def P_mm(self, t, E)->FlavorMatrix:
        decay_factor = np.exp(-self.gamma(E)*self.d)
        PND = FlavorMatrix.eye(self.mixing_params.basis_mass, extra_dims=E.shape)

        if self.mixing_params.mass_order == MassHierarchy.NORMAL:
            PND['NU_1','NU_3'] = 1 - decay_factor
            PND['NU_3','NU_3'] = decay_factor

        else:
            PND['NU_2','NU_2'] = decay_factor
            PND['NU_3','NU_2'] = 1 - decay_factor
            
        PND['NU_BAR','NU_BAR'] = PND['NU','NU']
        return PND

###############################################################################

class QuantumDecoherence(VacuumTransformation, ThreeFlavorTransformation):
    """Quantum Decoherence, where propagation in vacuum leads to equipartition
    of states. For a description and typical parameters, see M. V. dos Santos et al.,
    2023, arXiv:2306.17591.
    """
    def __init__(self, Gamma3=1e-27*u.eV, Gamma8=1e-27*u.eV, dist=10*u.kpc, n=0, E0=10*u.MeV):
        """Initialize transformation matrix.

        Parameters
        ----------
        Gamma3 : astropy.units.quantity.Quantity
            Quantum decoherence parameter; expect in eV.
        Gamma8 : astropy.units.quantity.Quantity
            Quantum decoherence parameter; expect in eV.
        dist : astropy.units.quantity.Quantity
            Distance to the supernova.
        n : float
            Exponent of power law for energy dependent quantum decoherence parameters,
            i.e. Gamma = Gamma0*(E/E0)**n. If not specified, it is taken as zero.
        E0 : astropy.units.quantity.Quantity
            Reference energy in the power law Gamma = Gamma0*(E/E0)**n. If not specified, 
            it is taken as 10 MeV. Note that if n = 0, quantum decoherence parameters are independent
            of E0.
        """
        super().__init__()
        self.Gamma3 = (Gamma3 / (c.hbar.to('eV s') * c.c)).to('1/kpc')
        self.Gamma8 = (Gamma8 / (c.hbar.to('eV s') * c.c)).to('1/kpc')
        self.d = dist
        self.n = n
        self.E0 = E0

    def P_mm(self, t, E)->FlavorMatrix: 
        PQD = FlavorMatrix.zeros(self.mixing_params.basis_mass, extra_dims=E.shape)
        x =  (E/self.E0)**self.n
        PQD['NU_1','NU_1'] = 1/3 + 1/2 * np.exp(-(self.Gamma3 + self.Gamma8 / 3) * x * self.d) \
                  + 1/6 * np.exp(-self.Gamma8 * x * self.d)

        PQD['NU_1','NU_2'] = 1/3 - 1/2 * np.exp(-(self.Gamma3 + self.Gamma8 / 3) * x * self.d) \
                  + 1/6 * np.exp(-self.Gamma8 * x * self.d)

        PQD['NU_1','NU_3'] = 1/3 - 1/3 * np.exp(-self.Gamma8 * x * self.d)

        PQD['NU_2',['NU_1','NU_2','NU_3']] = PQD['NU_1',['NU_2','NU_1','NU_3']]
        PQD['NU_2','NU_2'] = PQD['NU_1','NU_1']
        PQD['NU_2','NU_3'] = PQD['NU_1','NU_3']

        PQD['NU_3','NU_1'] = PQD['NU_1','NU_3']
        PQD['NU_3','NU_2'] = PQD['NU_2','NU_3']
    
        PQD['NU_3','NU_3'] = 1/3 + 2/3 * np.exp(-self.Gamma8 * x * self.d)
        PQD['NU_BAR','NU_BAR'] = PQD['NU','NU']
        return PQD
