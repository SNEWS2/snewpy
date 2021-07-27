# -*- coding: utf-8 -*-
"""
A submodule with classes used for supernova model files stored on disk. It
assumes models are available in a format usable by the AstroPy unified table
reader; see https://docs.astropy.org/en/stable/index.html for details.

Based on the ASTERIA (https://github.com/IceCubeOpenSource/ASTERIA) models
developed by Navya Uberoi and Spencer Griswold.

Updated summer 2020 by Jim Kneller & Arkin Worlikar.
"""

from abc import abstractmethod, ABC
from enum import IntEnum

import astropy
from astropy.io import ascii
from astropy.table import Table, join
from astropy.units.quantity import Quantity

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
from scipy.special import loggamma

import os
import re

import logging
logging.basicConfig(level=logging.INFO)

# Conditional import of photospline while we determine whether to make it a
# strong dependency of snewpy (July 2021).
try:
    import photospline
    from photospline import glam_fit, ndsparse, bspline, SplineTable
    usephotospline = True
except ImportError as e:
    usephotospline = False
    logging.warning(e)

import tarfile
import h5py

from .neutrino import Flavor
from .flavor_transformation import *


def get_value(x):
    """If quantity x has is an astropy Quantity with units, return just the
    value.

    Parameters
    ----------
    x : Quantity, float, or ndarray
        Input quantity.

    Returns
    -------
    value : float or ndarray
    """
    if type(x) == Quantity:
        return x.value
    return x


def get_closest(arr, x):
    """Get index of closest element in an array to input value.

    Parameters
    ----------
    arr : list or ndarray
        Array of values.
    x : float or int or str
        Value to search.

    Returns
    -------
    idx : int
        Index of closest element in the array.
    """
    return np.abs(np.asarray(arr) - x).argmin()


class SupernovaModel(ABC):
    """Base class defining an interface to a supernova model."""
    
    def __init__(self):
        pass
    
    @abstractmethod
    def get_time(self):
        """Returns
        -------
            returns array of snapshot times from the simulation
        """
        pass

    @abstractmethod
    def get_initialspectra(self, t, E):
        """Get neutrino spectra at the source.

        Parameters
        ----------
        t : float
            Time to evaluate spectra.
        E : float or ndarray
            Energies to evaluate spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of neutrino spectra, keyed by neutrino flavor.
        """
        pass

    def get_oscillatedspectra(self, t, E, flavor_xform):
        """Get neutrino spectra after applying oscillation.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial and oscillated spectra.
        E : astropy.Quantity or ndarray
            Energies to evaluate the initial and oscillated spectra.
        flavor_xform : FlavorTransformation
            An instance from the flavor_transformation module.

        Returns
        -------
        oscillatedspectra : dict
            Dictionary of oscillated spectra, keyed by neutrino flavor.
        """
        initialspectra = self.get_initialspectra(t, E)
        oscillatedspectra = {}

        oscillatedspectra[Flavor.NU_E] = \
            flavor_xform.prob_ee(t, E) * initialspectra[Flavor.NU_E] + \
            flavor_xform.prob_ex(t, E) * initialspectra[Flavor.NU_X]

        oscillatedspectra[Flavor.NU_X] = \
            flavor_xform.prob_xe(t, E) * initialspectra[Flavor.NU_E] + \
            flavor_xform.prob_xx(t, E) * initialspectra[Flavor.NU_X] 

        oscillatedspectra[Flavor.NU_E_BAR] = \
            flavor_xform.prob_eebar(t, E) * initialspectra[Flavor.NU_E_BAR] + \
            flavor_xform.prob_exbar(t, E) * initialspectra[Flavor.NU_X_BAR]

        oscillatedspectra[Flavor.NU_X_BAR] = \
            flavor_xform.prob_xebar(t, E) * initialspectra[Flavor.NU_E_BAR] + \
            flavor_xform.prob_xxbar(t, E) * initialspectra[Flavor.NU_X_BAR] 

        return oscillatedspectra   
    
    
class Analytic3Species(SupernovaModel):
    """Allow to generate an analytic model given total luminosity,
    average energy, and rms or pinch, for each species.
    """

    def __init__(self, filename):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data.
        """

        simtab = Table.read(filename,format='ascii')
        self.filename = filename
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        # Get grid of model times.
        self.time = simtab['TIME'] * u.s

        for flavor in Flavor:
            # Note: file only contains NU_E, NU_E_BAR, and NU_X, so double up
            # the use of NU_X for NU_X_BAR.
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor

            self.luminosity[flavor] = simtab['L_{}'.format(_flav.name)] * u.erg/u.s
            self.meanE[flavor] = simtab['E_{}'.format(_flav.name)] * u.MeV
            self.pinch[flavor] = simtab['ALPHA_{}'.format(_flav.name)]

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time
    
    def get_initialspectra(self, t, E):
        """Get neutrino spectra/luminosity curves after oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra={}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Analytic Model: {}\n'.format(self.filename)
        return mod

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Analytic Model**: {}\n\n'.format(self.filename)
        return mod  


class Nakazato_2013(SupernovaModel):
    """Set up a model based on simulations from Nakazato et al., ApJ S 205:2,
    2013 and ApJ 804:75, 2015. See also http://asphwww.ph.noda.tus.ac.jp/snn/.
    """

    def __init__(self, filename):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """
        # Store model metadata.
        if 't_rev' in filename:
            self.progenitor_mass = float(filename.split('-')[-1].strip('s%.fits')) * u.Msun
            self.revival_time = float(filename.split('-')[-2].strip('t_rev%ms')) * u.ms
            self.metallicity = float(filename.split('-')[-3].strip('z%'))
            self.EOS = filename.split('-')[-4].upper()
        # No revival time because the explosion "failed" (BH formation).
        else:
            self.progenitor_mass = float(filename.split('-')[-1].strip('s%.fits')) * u.Msun
            self.metallicity = float(filename.split('-')[-2].strip('z%'))
            self.revival_time = 0 * u.ms
            self.EOS = filename.split('-')[-4].upper()

        # Read FITS table using the astropy reader.
        simtab = Table.read(filename)
        self.filename = os.path.basename(filename)

        # Get grid of model times.
        self.time = simtab['TIME'].to('s')

        # Set up dictionary of luminosity, mean energy and shape parameter
        # alpha, keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        for flavor in Flavor:
            # Note: file only contains NU_E, NU_E_BAR, and NU_X, so double up
            # the use of NU_X for NU_X_BAR.
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor

            self.luminosity[flavor] = simtab['L_{}'.format(_flav.name)].to('erg/s')
            self.meanE[flavor] = simtab['E_{}'.format(_flav.name)].to('MeV')
            self.pinch[flavor] = simtab['ALPHA_{}'.format(_flav.name)]

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time
    
    def get_initialspectra(self, t, E):
        """Get neutrino spectra/luminosity at the source.

        Parameters
        ----------
        t : float or astropy.Quantity
            Time to evaluate initial spectra.
        E : float or ndarray of astropy.Quantity
            Energies to evaluate the initial spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a) 
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Nakazato_2013 Model: {}\n'.format(self.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Metallicity     : {}'.format(self.metallicity),
             'Revival time    : {}'.format(self.revival_time),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)
        
    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Nakazato_2013 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
             '|Metallicity | ${:g}$|'.format(self.metallicity),
             '|Revival time | ${0.value:g}$ {0.unit:latex}|'.format(self.revival_time),
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)


class Sukhbold_2015(SupernovaModel):
    """Set up a model based on simulations from Sukhbold et al., ApJ 821:38,2016. Models were shared privately by email.
    """

    def __init__(self, filename):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """
        # Store model metadata.
        self.progenitor_mass = float(filename.split('-')[-1].strip('z%.fits')) * u.Msun
        self.EOS = filename.split('-')[-2]

        # Read FITS table using the astropy unified Table reader.
        simtab = Table.read(filename)
        self.filename = os.path.basename(filename)

        # Get grid of model times.
        self.time = simtab['TIME'].to('s')

        # Set up dictionary of luminosity, mean energy, and shape parameter,
        # keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        for flavor in Flavor:
            self.luminosity[flavor] = simtab['L_{}'.format(flavor.name)].to('erg/s')
            self.meanE[flavor] = simtab['E_{}'.format(flavor.name)].to('MeV')
            self.pinch[flavor] = simtab['ALPHA_{}'.format(flavor.name)]
            
    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time
    
    def get_initialspectra(self,t,E):
        """Get neutrino spectra/luminosity curves after oscillation.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity or ndarray
            Energies to evaluate the initial spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)>, and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the global time grid or
        # the interpolation will not work properly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a) 
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Sukhbold_2015 Model: {}\n'.format(self.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Sukhbold_2015 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)


class Bollig_2016(SupernovaModel):
    """Set up a model based on simulations from Bollig et al. (2016). Models were taken, with permission, from the Garching Supernova Archive.
    """
    def __init__(self, filename, eos='LS220'):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux.
        eos : string
            Equation of state used in simulation.
        """
        self.time = {}
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        # Store model metadata.
        self.filename = os.path.basename(filename)
        self.EOS = eos
        self.progenitor_mass = float(self.filename.strip('s%c')) * u.Msun

        # Read through the several ASCII files for the chosen simulation and
        # merge the data into one giant table.
        mergtab = None
        for flavor in Flavor:
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor
            _sfx = _flav.name.replace('_', '').lower()
            _filename = '{}_{}_{}'.format(filename, eos, _sfx)
            _lname  = 'L_{}'.format(flavor.name)
            _ename  = 'E_{}'.format(flavor.name)
            _e2name = 'E2_{}'.format(flavor.name)
            _aname  = 'ALPHA_{}'.format(flavor.name)

            simtab = Table.read(_filename,
                                names=['TIME', _lname, _ename, _e2name],
                                format='ascii')
            simtab['TIME'].unit = 's'
            simtab[_lname].unit = '1e51 erg/s'
            simtab[_aname] = (2*simtab[_ename]**2 - simtab[_e2name]) / (simtab[_e2name] - simtab[_ename]**2)
            simtab[_ename].unit = 'MeV'
            del simtab[_e2name]

            if mergtab is None:
                mergtab = simtab
            else:
                mergtab = join(mergtab, simtab, keys='TIME', join_type='left')
                mergtab[_lname].fill_value = 0.
                mergtab[_ename].fill_value = 0.
                mergtab[_aname].fill_value = 0.
        simtab = mergtab.filled()

        self.time = simtab['TIME'].to('s')

        for flavor in Flavor:
            # Set the dictionary of luminosity, mean energy, and shape
            # parameter keyed by NU_E, NU_X, NU_E_BAR, NU_X_BAR.
            _lname  = 'L_{}'.format(flavor.name)
            self.luminosity[flavor] = simtab[_lname].to('erg/s')

            _ename  = 'E_{}'.format(flavor.name)
            self.meanE[flavor] = simtab[_ename].to('MeV')

            _aname  = 'ALPHA_{}'.format(flavor.name)
            self.pinch[flavor] = simtab[_aname]

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time

    def get_initialspectra(self,t,E):
        """Get neutrino spectra/luminosity curves before oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial spectra.
        E : float or ndarray
            Energies to evaluate the initial spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            if L==0:
                initialspectra[flavor] = 0.0*E/(u.erg*u.s)
            else:
                # For numerical stability, evaluate log PDF and then exponentiate.
                initialspectra[flavor] = \
                  np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                        - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Bollig_2016 Model: {}\n'.format(self.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Bollig_2016 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)


class OConnor_2015(SupernovaModel):
    """Set up a model based on the black hole formation simulation in O'Connor (2015). 
    """
    def __init__(self, filename, eos='LS220'):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        simtab = Table.read(filename, 
                     names= ['TIME','L_NU_E','L_NU_E_BAR','L_NU_X',
                                    'E_NU_E','E_NU_E_BAR','E_NU_X',
                                    'RMS_NU_E','RMS_NU_E_BAR','RMS_NU_X'],
                     format='ascii')

        header = ascii.read(simtab.meta['comments'], delimiter='=',format='no_header', names=['key', 'val'])
        tbounce = float(header['val'][0])
        simtab['TIME'] -= tbounce
        
        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2)/(simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2)/(simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2)/(simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        # SYB: double-check on this factor of 4. Should be factor of 2?
        simtab['L_NU_X'] /= 4.0

        self.filename = 'OConnor2015_s40WH07_LS220'
        self.EOS = eos
        self.progenitor_mass = 40 * u.Msun

        # Get grid of model times.
        self.time = simtab['TIME'] * u.s

        # Set up dictionary of luminosity, mean energy and shape parameter
        # alpha, keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        for flavor in Flavor:
            # Note: file only contains NU_E, NU_E_BAR, and NU_X, so double up
            # the use of NU_X for NU_X_BAR.
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor

            self.luminosity[flavor] = simtab['L_{}'.format(_flav.name)] * u.erg/u.s
            self.meanE[flavor] = simtab['E_{}'.format(_flav.name)] * u.MeV
            self.pinch[flavor] = simtab['ALPHA_{}'.format(_flav.name)]

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time

    def get_initialspectra(self,t,E):
        """Get neutrino spectra/luminosity curves before oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'OConnor_2015 Model: {}\n'.format(self.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**OConnor_2015 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

class Zha_2021(SupernovaModel):
    """Set up a model based on the hadron-quark phse transition models from Zha et al. 2021. 
    """
    def __init__(self, filename, eos='STOS_B145'):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        simtab = Table.read(filename, 
                     names= ['TIME','L_NU_E','L_NU_E_BAR','L_NU_X',
                                    'E_NU_E','E_NU_E_BAR','E_NU_X',
                                    'RMS_NU_E','RMS_NU_E_BAR','RMS_NU_X'],
                     format='ascii')

        header = ascii.read(simtab.meta['comments'], delimiter='=',format='no_header', names=['key', 'val'])
        tbounce = float(header['val'][0])
        simtab['TIME'] -= tbounce
        
        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2)/(simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2)/(simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2)/(simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        # SYB: double-check on this factor of 4. Should be factor of 2?
        simtab['L_NU_X'] /= 4.0

        #prevent neagative lums
        simtab['L_NU_E'][simtab['L_NU_E'] < 0] = 1
        simtab['L_NU_E_BAR'][simtab['L_NU_E_BAR'] < 0] = 1
        simtab['L_NU_X'][simtab['L_NU_X'] < 0] = 1

        
        basename =os.path.basename(filename)[:-4]
        
        self.filename = 'Zha2021_'+basename
        self.EOS = eos
        self.progenitor_mass =  float(basename[1:])* u.Msun

        # Get grid of model times.
        self.time = simtab['TIME'] * u.s

        # Set up dictionary of luminosity, mean energy and shape parameter
        # alpha, keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        for flavor in Flavor:
            # Note: file only contains NU_E, NU_E_BAR, and NU_X, so double up
            # the use of NU_X for NU_X_BAR.
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor

            self.luminosity[flavor] = simtab['L_{}'.format(_flav.name)] * u.erg/u.s
            self.meanE[flavor] = simtab['E_{}'.format(_flav.name)] * u.MeV
            self.pinch[flavor] = simtab['ALPHA_{}'.format(_flav.name)]

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time

    def get_initialspectra(self,t,E):
        """Get neutrino spectra/luminosity curves before oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Zha_2021 Model: {}\n'.format(self.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Zha_2021 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)    


class Warren_2020(SupernovaModel):
    """Set up a model based on simulations from Warren et al. (2020)."""

    def __init__(self, filename, eos='LS220'):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        # Read data from HDF5 files, then store.
        f = h5py.File(filename, 'r')
        simtab = Table()

        for i in range(len(f['nue_data']['lum'])):
            if f['sim_data']['shock_radius'][i][1] > 0.00001:
                bounce = f['sim_data']['shock_radius'][i][0]
                break

        simtab['TIME'] = f['nue_data']['lum'][:, 0] - bounce
        simtab['L_NU_E'] = f['nue_data']['lum'][:, 1] * 1e51
        simtab['L_NU_E_BAR'] = f['nuae_data']['lum'][:, 1] * 1e51
        simtab['L_NU_X'] = f['nux_data']['lum'][:, 1] * 1e51 / 4.0
        simtab['E_NU_E'] = f['nue_data']['avg_energy'][:, 1]
        simtab['E_NU_E_BAR'] = f['nuae_data']['avg_energy'][:, 1]
        simtab['E_NU_X'] = f['nux_data']['avg_energy'][:, 1]
        simtab['RMS_NU_E'] = f['nue_data']['rms_energy'][:, 1]
        simtab['RMS_NU_E_BAR'] = f['nuae_data']['rms_energy'][:, 1]
        simtab['RMS_NU_X'] = f['nux_data']['rms_energy'][:, 1]

        simtab['ALPHA_NU_E'] = (2.0 * simtab['E_NU_E'] ** 2 - simtab['RMS_NU_E'] ** 2) / (simtab['RMS_NU_E'] ** 2 - simtab['E_NU_E'] ** 2)
        simtab['ALPHA_NU_E_BAR'] = (2.0 * simtab['E_NU_E_BAR'] ** 2 - simtab['RMS_NU_E_BAR'] ** 2) / (simtab['RMS_NU_E_BAR'] ** 2 - simtab['E_NU_E_BAR'] ** 2)
        simtab['ALPHA_NU_X'] = (2.0 * simtab['E_NU_X'] ** 2 - simtab['RMS_NU_X'] ** 2) / (simtab['RMS_NU_X'] ** 2 - simtab['E_NU_X'] ** 2)

        # Set model metadata.
        self.filename = os.path.basename(filename)
        self.EOS = eos
        self.progenitor_mass = float(filename.split('_')[-1].strip('m%.h5')) * u.Msun
        self.turbmixing_param = float(filename.split('_')[-2].strip('a%'))

        # Get grid of model times.
        self.time = simtab['TIME'] * u.s

        # Set up dictionary of luminosity, mean energy and shape parameter
        # alpha, keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        for flavor in Flavor:
            # Note: file only contains NU_E, NU_E_BAR, and NU_X, so double up
            # the use of NU_X for NU_X_BAR.
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor

            self.luminosity[flavor] = simtab['L_{}'.format(_flav.name)] * u.erg/u.s
            self.meanE[flavor] = simtab['E_{}'.format(_flav.name)] * u.MeV
            self.pinch[flavor] = simtab['ALPHA_{}'.format(_flav.name)]

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time

    def get_initialspectra(self, t, E):
        """Get neutrino spectra/luminosity curves before oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Warren_2020 Model: {}\n'.format(self.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Turb. mix param : {}'.format(self.turbmixing_param),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Warren_2020 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
             '|Turb. mixing param. | {}|'.format(self.turbmixing_param),
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)


class Kuroda_2020(SupernovaModel):
    """Set up a model based on simulations from Kuroda et al. (2020)."""

    def __init__(self, filename, eos='LS220'):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        eos : string
            Equation of state used in simulation
        """
        # Load up model metadata.
        self.filename = filename
        self.EOS = eos

        # Read ASCII data.
        simtab = Table.read(filename, format='ascii')

        # Get grid of model times.
        self.time = (simtab['Tpb[ms]'] * u.ms).to('s')

        # Set up dictionary of luminosity, mean energy and shape parameter
        # alpha, keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}

        for flavor in Flavor:
            # Note: file only contains NU_E, NU_E_BAR, and NU_X, so double up
            # the use of NU_X for NU_X_BAR.
            _flav = Flavor.NU_X if flavor == Flavor.NU_X_BAR else flavor
            if _flav.is_neutrino:
                _fkey = _flav.name.lower()
            else:
                _fkey = _flav.name.strip('%_BAR').lower().replace('_', '_a')

            self.luminosity[flavor] = simtab['<L{}>'.format(_fkey)] * 1e51 * u.erg/u.s
            self.meanE[flavor] = simtab['<E{}>'.format(_fkey)] * u.MeV

            # There is no pinch parameter so use alpha=2.0.
            self.pinch[flavor] = np.full_like(self.meanE[flavor].value, 2.)

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.time

    def get_initialspectra(self, t, E):
        """Get neutrino spectra/luminosity curves before oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid division by zero in energy PDF below.
        E[E==0] = np.finfo(float).eps * E.unit

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to('erg').value

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in Flavor:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Kuroda_2020 Model: {}\n'.format(self.filename)
        s = ['Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Kuroda_2020 Model**: {}\n\n'.format(self.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|EOS | {}|'.format(self.EOS)
             ]
        return mod + '\n'.join(s)


class Fornax_2021_2D(SupernovaModel):
    """Model based on axisymmetric simulations from A. Burrows and D.  Vartanyan, Nature 589:29, 2021. Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.2d/.
    """

    def __init__(self, filename):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """
        # Set up model metadata.
        self.progenitor_mass = float(filename.split('_')[-3][:-1]) * u.Msun

        # Conversion of flavor to key name in the model HDF5 file.
        self._flavorkeys = { Flavor.NU_E : 'nu0',
                             Flavor.NU_E_BAR : 'nu1',
                             Flavor.NU_X : 'nu2',
                             Flavor.NU_X_BAR : 'nu2' }

        # Open HDF5 data file.
        self._h5file = h5py.File(filename, 'r')

        # Get grid of model times.
        self.time = self._h5file['nu0'].attrs['time'] * u.s

        # Compute luminosity by integrating over model energy bins.
        self.luminosity = {}
        for flavor in Flavor:
            key = self._flavorkeys[flavor]
            dE = np.asarray(self._h5file[key]['degroup'])
            n = len(dE[0])
            dLdE = np.zeros((len(self.time), n), dtype=float)
            for i in range(n):
                dLdE[:,i] = self._h5file[key]["g{}".format(i)]

            # Note factor of 0.5 in nu_x and nu_x_bar.
            factor = 1e50*u.erg/u.s if flavor.is_electron else 0.5*1e50*u.erg/u.s
            self.luminosity[flavor] = np.sum(dLdE*dE, axis=1) * factor

    def get_time(self):
        return self.time

    def get_initialspectra(self, t, E):
        """Get neutrino spectra/luminosity curves after oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid "division by zero" in retrieval of the spectrum.
        E[E == 0] = np.finfo(float).eps * E.unit
        logE = np.log10(E.to('MeV').value)

        # Make sure the input time uses the same units as the model time grid.
        # Convert input time to a time index.
        t = t.to(self.time.unit)
        j = (np.abs(t - self.time)).argmin()

        for flavor in Flavor:
            key = self._flavorkeys[flavor]

            # Energy in units of MeV.
            _E = self._h5file[key]['egroup'][j]

            # Pad log(E) array with values where flux is fixed to zero.
            _logE = np.log10(_E)
            _dlogE = np.diff(_logE)
            _logEbins = np.insert(_logE, 0, np.log10(np.finfo(float).eps))
            _logEbins = np.append(_logEbins, _logE[-1] + _dlogE[-1])

            # Spectrum in units of 1e50 erg/s/MeV.
            # Pad with values where flux is fixed to zero.
            _dLdE = np.asarray([0.] + [self._h5file[key]['g{}'.format(i)][j] for i in range(12)] + [0.])

            # Linear interpolation in flux.
            factor = 1e50*u.MeV
            if not flavor.is_electron:
                # Model flavors (internally) are nu_e, nu_e_bar, and nu_x, which stands
                # for nu_mu(_bar) and nu_tau(_bar), making the flux 4x higher than nu_e and nu_e_bar.
                # Since we separate NU_X and NU_X_BAR here, multiply by 2x, not 4x.
                factor *= 0.5

            initialspectra[flavor] = np.interp(logE, _logEbins, _dLdE) * factor

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Fornax 2D Model: {}\n'.format(self._h5file.filename)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass)
            ]
        return mod + '\n'.join(s)

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Fornax 2D Model**: {}\n\n'.format(self._h5file.filename)
        s = ['|Parameter|Value|',
             '|:---------|:-----:|',
             '|Progenitor mass | ${0.value:g}$ {0.unit:latex}|'.format(self.progenitor_mass),
            ]
        return mod + '\n'.join(s)


class SNOwGLoBES:
    """A model that does not inherit from SupernovaModel (yet) and imports a group of SNOwGLoBES files."""

    def __init__(self, tarfilename):
        """Initialize model from a tar archive.

        Parameters
        ----------
        tarfilename: str
            Absolute or relative path to tar archive with SNOwGLoBES files.
        """
        self.tfname = tarfilename
        tf = tarfile.open(self.tfname)

        # For now just pull out the "NoOsc" files.
        datafiles = sorted([f.name for f in tf if '.dat' in f.name])
        noosc = [df for df in datafiles if 'NoOsc' in df]
        noosc.sort(key=len)

        # Loop through the noosc files and pull out the number fluxes.
        self.time = []
        self.energy = None
        self.flux = {}
        self.fmin = 1e99
        self.fmax = -1e99

        for nooscfile in noosc:
            with tf.extractfile(nooscfile) as f:
                logging.debug('Reading {}'.format(nooscfile))
                meta = f.readline()
                metatext = meta.decode('utf-8')
                t = float(metatext.split('TBinMid=')[-1].split('sec')[0])
                dt = float(metatext.split('tBinWidth=')[-1].split('s')[0])
                dE = float(metatext.split('eBinWidth=')[-1].split('MeV')[0])

                data = Table.read(f, format='ascii.commented_header', header_start=-1)
                data.meta['t'] = t
                data.meta['dt'] = dt
                data.meta['dE'] = dE

                self.time.append(t)
                if self.energy is None:
                    self.energy = (data['E(GeV)'].data*1000).tolist()

            for flavor in ['NuE', 'NuMu', 'NuTau', 'aNuE', 'aNuMu', 'aNuTau']:
                if flavor in self.flux:
                    self.flux[flavor].append(data[flavor].data.tolist())
                else:
                    self.flux[flavor] = [data[flavor].data.tolist()]

        # We now have a table with rows=times and columns=energies. Transpose
        # so that rows=energy and cols=time.
        for k, v in self.flux.items():
            self.flux[k] = np.transpose(self.flux[k])
            self.fmin = np.minimum(self.fmin, np.min(self.flux[k]))
            self.fmax = np.maximum(self.fmax, np.max(self.flux[k]))

    def get_fluence(self, t):
        """Return the fluence at a given time t.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        fluence : dict
            A dictionary giving fluence at time t, keyed by flavor.
        """
        idx = get_closest(self.time, t)

        fluence = {}
        for k, fl in self.flux.items():
            fluence[k] = fl[:,idx]

        return fluence
