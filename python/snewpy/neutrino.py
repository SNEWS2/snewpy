# -*- coding: utf-8 -*-
"""Basic neutrino properties.
"""

from enum import IntEnum
from astropy import units as u


class MassHierarchy(IntEnum):
    """Neutrino mass ordering: normal or inverted.
    """
    NORMAL = 1
    INVERTED = 2


class Flavor(IntEnum):
    """Enumeration of CCSN Neutrino flavors.
    """
    NU_E = 2
    NU_E_BAR = 1
    NU_X = 3
    NU_X_BAR = 0
    
    def to_tex(self):
        """LaTeX-compatible string representations of flavor.
        """
        
        if '_BAR' in self.name:
            return r'$\overline{{\nu}}_{0}$'.format(self.name[3].lower())
        return r'$\{0}$'.format(self.name.lower())

    @property
    def is_electron(self):
        return self.value in (Flavor.NU_E.value, Flavor.NU_E_BAR.value)

    @property
    def is_neutrino(self):
        return self.value in (Flavor.NU_E.value, Flavor.NU_X.value)

    @property
    def is_antineutrino(self):
        return self.value in (Flavor.NU_E_BAR.value, Flavor.NU_X_BAR.value)


class MixingParameters:
    """Best-fit parameters of the PMNS matrix and mass differences, assuming
    three neutrino flavors. See www.nu-fit.org for current global fits.
    """
    def __init__(self, mh=MassHierarchy.NORMAL):
        """Initialize the neutrino mixing parameters.

        Parameters
        ----------
        mh : MassHierarchy
            Desired mass ordering: NORMAL or INVERTED.
        """
        if type(mh) == MassHierarchy:
            self.mass_order = mh
        else:
            raise TypeError('mh must be of type MassHierarchy')

        # Values from JHEP 09 (2020) 178 [arXiv:2007.14792] and www.nu-fit.org.
        # The reported precision is not significant given current
        # uncertainties, but is useful for comparing to the table of
        # parameters presented on nu-fit.org.
        if self.mass_order == MassHierarchy.NORMAL:
            # Note: in NH, the mass splittings are: m1..m2..............m3.
            self.theta12 = 33.44 * u.deg
            self.theta13 =  8.57 * u.deg
            self.theta23 = 49.20 * u.deg
            self.deltaCP = 197 * u.deg
            self.dm21_2  = 7.42e-5 * u.eV**2
            self.dm32_2  = 2.517e-3 * u.eV**2
        else:
            # Note: in IH, the mass splittings are: m3..............m1..m2.
            self.theta12 = 33.45 * u.deg
            self.theta13 =  8.60 * u.deg
            self.theta23 = 49.30 * u.deg
            self.deltaCP = 282 * u.deg
            self.dm21_2  = 7.42e-5 * u.eV**2
            self.dm31_2  = -2.498e-3 * u.eV**2

    def get_mixing_angles(self):
        """Mixing angles of the PMNS matrix.
        
        Returns
        -------
        angles : tuple
            Angles theta12, theta13, theta23.
        """
        return (self.theta12, self.theta13, self.theta23)

