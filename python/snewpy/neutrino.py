# -*- coding: utf-8 -*-
"""This module implements basic neutrino properties that are used throughout SNEWPY."""

from enum import IntEnum
from astropy import units as u
from dataclasses import dataclass

class MassHierarchy(IntEnum):
    """Neutrino mass ordering: ``NORMAL`` or ``INVERTED``."""
    NORMAL = 1
    INVERTED = 2


class Flavor(IntEnum):
    """Enumeration of CCSN Neutrino flavors."""
    NU_E = 2
    NU_E_BAR = 1
    NU_X = 3
    NU_X_BAR = 0
    
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
class MixingParameters:
    """Mixing angles and mass differences, assuming three neutrino flavors.
    This class contains the default values used throughout SNEWPY, which are
    based on `NuFIT 5.0 <http://www.nu-fit.org>`_ results from July 2020,
    published in `JHEP 09 (2020) 178 <https://dx.doi.org/10.1007/JHEP09(2020)178>`_
    [`arXiv:2007.14792 <https://arxiv.org/abs/2007.14792>`_].
    """
    mass_order: MassHierarchy = MassHierarchy.NORMAL,
    #mixing angles
    theta12: u.Quantity[u.deg] = 33.44<< u.deg
    theta13: u.Quantity[u.deg] = 8.57 << u.deg
    theta23: u.Quantity[u.deg] = 49.20 << u.deg
    #CP violation phase         
    deltaCP: u.Quantity[u.deg] = 197  << u.deg
    #square mass difference     
    dm21_2: u.Quantity[u.eV**2] = 7.42e-5  << u.eV**2
    dm32_2: u.Quantity[u.eV**2] = 2.517e-3 << u.eV**2
    # Note: in IH, the mass splittings are: m3..............m1..m2.

    def get_mixing_angles(self):
        """Mixing angles of the PMNS matrix.
        
        Returns
        -------
        tuple
            Angles theta12, theta13, theta23.
        """
        return (self.theta12, self.theta13, self.theta23)

    
# Values from JHEP 09 (2020) 178 [arXiv:2007.14792] and www.nu-fit.org.
NuFIT50_NH = MixingParameters(
        mass_order = MassHierarchy.NORMAL,
        theta12 = 33.44<< u.deg,
        theta13 = 8.57 << u.deg,
        theta23 = 49.20 << u.deg,
        deltaCP = 197  << u.deg,
        dm21_2 =  7.42e-5  << u.eV**2,
        dm32_2 =  2.517e-3 << u.eV**2
)
# Values from JHEP 09 (2020) 178 [arXiv:2007.14792] and www.nu-fit.org.
NuFIT50_IH = MixingParameters(
        mass_order = MassHierarchy.INVERTED,
        theta12 = 33.45 << u.deg,
        theta13 = 8.60 << u.deg,
        theta23 = 49.30 << u.deg,
        deltaCP = 282 << u.deg,
        dm21_2 = 7.42e-5 << u.eV**2,
        dm32_2 = -2.498e-3 << u.eV**2
)

NuFIT52_NH = MixingParameters(
        mass_order = MassHierarchy.NORMAL,
        theta12 = 33.41 << u.deg,
        theta13 = 8.58 << u.deg,
        theta23 = 42.20 << u.deg,
        deltaCP = 232 << u.deg,
        dm21_2 = 7.41e-5 << u.eV**2,
        dm32_2 = 2.507e-3 << u.eV**2
)
NuFIT52_IH = MixingParameters(
        mass_order = MassHierarchy.INVERTED,
        theta12 = 33.41 << u.deg,
        theta13 = 8.57 << u.deg,
        theta23 = 49.00 << u.deg,
        deltaCP = 276 << u.deg,
        dm21_2 = 7.41e-5 << u.eV**2,
        dm32_2 = -2.486e-3 << u.eV**2
)

# Values from https://pdg.lbl.gov
PDG2022_NH = MixingParameters(
        mass_order = MassHierarchy.NORMAL,
        theta12 = 33.65 << u.deg,
        theta13 = 8.53 << u.deg,
        theta23 = 47.64 << u.deg,
        deltaCP = 245 << u.deg,
        dm21_2 = 7.53e-5 << u.eV**2,
        dm32_2 = 2.453e-3 << u.eV**2
)
PDG2022_IH = MixingParameters(
        mass_order = MassHierarchy.INVERTED,
        theta12 = 33.65 << u.deg,
        theta13 = 8.53 << u.deg,
        theta23 = 47.24 << u.deg,
        deltaCP = 245 << u.deg,
        dm21_2 = 7.53e-5 << u.eV**2,
        dm32_2 = -2.536e-3 << u.eV**2
)
