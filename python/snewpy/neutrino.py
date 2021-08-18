# -*- coding: utf-8 -*-
"""This module implements basic neutrino properties that are used throughout SNEWPY."""

from enum import IntEnum
from astropy import units as u


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


class MixingParameters:
    """Mixing angles and mass differences, assuming three neutrino flavors.
    
    This class contains the default values used throughout SNEWPY, which are
    based on `NuFIT 5.0 <http://www.nu-fit.org>`_ results from July 2020,
    published in `JHEP 09 (2020) 178 <https://dx.doi.org/10.1007/JHEP09(2020)178>`_
    [`arXiv:2007.14792 <https://arxiv.org/abs/2007.14792>`_].
    Note that the best fit values vary between normal and inverted mass hierarchy.
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
        tuple
            Angles theta12, theta13, theta23.
        """
        return (self.theta12, self.theta13, self.theta23)

