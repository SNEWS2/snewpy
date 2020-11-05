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

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
from scipy.special import loggamma

import os
import logging
import re

import tarfile
import h5py

from .neutrino import Flavor
from .flavor_transformation import *


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

    def get_oscillatedspectra(self, t, E ):
        """Get neutrino spectra after applying oscillation.

        Parameters
        ----------
        t : float
            Time to evaluate initial and oscillated spectra.
        E : float or ndarray
            Energies to evaluate the initial and oscillated spectra.

        Returns
        -------
        oscillatedspectra : dict
            Dictionary of oscillated spectra, keyed by neutrino flavor.
        """
        initialspectra = self.get_initialspectra(t, E)
        oscillatedspectra = {}

        oscillatedspectra[Flavor.NU_E] = \
            self.FT.prob_ee(t, E) * initialspectra[Flavor.NU_E] + \
            self.FT.prob_ex(t, E) * initialspectra[Flavor.NU_X]

        oscillatedspectra[Flavor.NU_X] = \
            self.FT.prob_xe(t, E) * initialspectra[Flavor.NU_E] + \
            self.FT.prob_xx(t, E) * initialspectra[Flavor.NU_X] 

        oscillatedspectra[Flavor.NU_E_BAR] = \
            self.FT.prob_eebar(t, E) * initialspectra[Flavor.NU_E_BAR] + \
            self.FT.prob_exbar(t, E) * initialspectra[Flavor.NU_X_BAR]

        oscillatedspectra[Flavor.NU_X_BAR] = \
            self.FT.prob_xebar(t, E) * initialspectra[Flavor.NU_E_BAR] + \
            self.FT.prob_xxbar(t, E) * initialspectra[Flavor.NU_X_BAR] 

        return oscillatedspectra    


class Nakazato_2013(SupernovaModel):
    """Set up a model based on simulations from Nakazato et al., ApJ S 205:2,
    2013 and ApJ 804:75, 2015. See also http://asphwww.ph.noda.tus.ac.jp/snn/.
    """

    def __init__(self, filename, flavor_xform):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        flavor_xform : FlavorTransformation
            Flavor transformation object with survival probabilities.
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
        self.fitsfile = os.path.basename(filename)

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
        self.FT = flavor_xform

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
            L  = np.interp(t, self.time, self.luminosity[flavor].to('erg/s'))
            Ea = np.interp(t, self.time, self.meanE[flavor].to('erg'))
            a  = np.interp(t, self.time, self.pinch[flavor])

            # For numerical stability, evaluate log PDF and then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a) 
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra

    def __repr__(self):
        """Default representation of the model.
        """
        mod = 'Nakazato_2013 Model: {}\n'.format(self.fitsfile)
        s = ['Progenitor mass : {}'.format(self.progenitor_mass),
             'Metallicity     : {}'.format(self.metallicity),
             'Revival time    : {}'.format(self.revival_time),
             'Eq. of state    : {}'.format(self.EOS)
             ]
        return mod + '\n'.join(s)
        
    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = '**Nakazato_2013 Model**: {}\n\n'.format(self.fitsfile)
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

    def __init__(self, filename, flavor_xform):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        flavor_xform : FlavorTransformation
            Flavor transformation object with survival probabilities.
        """
        self.file = Table.read(filename)
        self.filename = filename
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d(self.get_time(), self.get_luminosity(flavor))
            self.meanE[flavor] = interp1d(self.get_time(), self.get_mean_energy(flavor))
            self.pinch[flavor] = interp1d(self.get_time(), self.get_pinch_param(flavor))
        self.FT = flavor_xform
            
    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.file['TIME']
    
    def get_luminosity(self, flavor):
        """Get model luminosity L_nu.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        luminosity : ndarray
            Grid of luminosity values (erg/s) for this flavor.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['L_{}'.format(flavor.name.upper())]
        
    def get_mean_energy(self, flavor):
        """Get model mean energy <E_nu>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['E_{}'.format(flavor.name.upper())]
    
    def get_pinch_param(self, flavor):
        """Get spectral pinch parameter alpha(t).

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        alpha : ndarray
            Grid of alpha versus time.
        """
        if (flavor == Flavor.NU_X_BAR):
            flavor = Flavor.NU_X
        return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
    def get_EOS(self):
        """Model equation of state.

        Returns
        -------
        eos : str
            Model equation of state.
        """
        return self.filename.split('-')[1]
    
    def get_progenitor_mass(self):
        """Progenitor mass.

        Returns
        -------
        mass : float
            Progenitor mass, in units of solar mass.
        """
        return float(self.split('-')[-1].split('.')[0].strip('s'))

    def get_initialspectra(self,t,E):
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
        for flavor in Flavor:
            L = self.luminosity[flavor](t)
            Ea = self.meanE[flavor](t)          # <E_nu(t)>
            Ea = Ea*1e6 * 1.60218e-12
            a = self.pinch[flavor](t)           # alpha_nu(t)
            E[E==0] = np.finfo(float).eps       # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a) 
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra


class Bollig_2016(SupernovaModel):
    """Set up a model based on simulations from Bollig et al. (2016). Models were taken, with permission from the Garching Supernova Archive
    """
    def __init__(self, filename, FlavorTransformation, eos="LS220"):
        """Initialize model.

        Parameters
        ----------
        filename : str
        Absolute or relative path to file prefix, we add nue/nuebar/nux
        flavor_xform : FlavorTransformation
        Flavor transformation object with survival probabilities.
        eos : string
        Equation of state used in simulation
        """
        nue = Table.read(filename+"_"+eos+"_nue",names=["TIME","L_NU_E","E_NU_E","MS_NU_E"],format='ascii')
        nue['ALPHA_NU_E'] = (2.0*nue['E_NU_E']**2 - nue['MS_NU_E'])/(nue['MS_NU_E'] - nue['E_NU_E']**2)
        nue["L_NU_E"] *= 1e51
        nuebar = Table.read(filename+"_"+eos+"_nuebar",names=["TIME","L_NU_E_BAR","E_NU_E_BAR","MS_NU_E_BAR"],format='ascii')
        nuebar['ALPHA_NU_E_BAR'] = (2.0*nuebar['E_NU_E_BAR']**2 - nuebar['MS_NU_E_BAR'])/(nuebar['MS_NU_E_BAR'] - nuebar['E_NU_E_BAR']**2)
        nuebar["L_NU_E_BAR"] *= 1e51
        nux = Table.read(filename+"_"+eos+"_nux",names=["TIME","L_NU_X","E_NU_X","MS_NU_X"],format='ascii')
        nux['ALPHA_NU_X'] = (2.0*nux['E_NU_X']**2 - nux['MS_NU_X'])/(nux['MS_NU_X'] - nux['E_NU_X']**2)
        nux["L_NU_X"] *= 1e51
        tmptable = join(nuebar,nux,keys="TIME",join_type="left")
        self.file = join(nue,tmptable,keys="TIME",join_type="left")
        self.filename = filename
        self.shortname = "Bollig2016_"+filename.split('/')[-1]
        self.eos = eos
        self.luminosity={}
        self.meanE={}
        self.pinch={}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d( self.get_time() , self.get_luminosity(flavor) )
            self.meanE[flavor] = interp1d( self.get_time() , self.get_mean_energy(flavor) )
            self.pinch[flavor] = interp1d( self.get_time() , self.get_pinch_param(flavor) )
        self.FT=FlavorTransformation

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.file['TIME']

    def get_luminosity(self, flavor):
        """Get model luminosity L_nu.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        luminosity : ndarray
            Grid of luminosity values (erg/s) for this flavor.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['L_{}'.format(flavor.name.upper())]

    def get_mean_energy(self, flavor):
        """Get model mean energy <E_nu>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['E_{}'.format(flavor.name.upper())]

    def get_meansq_energy(self, flavor):
        """Get model mean squared energy <E_nu^2>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean squared energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['MS_{}'.format(flavor.name.upper())]

    def get_pinch_param(self, flavor):
        """Get spectral pinch parameter alpha(t).

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        alpha : ndarray
            Grid of alpha versus time.
        """
        if (flavor == Flavor.NU_X_BAR):
            flavor = Flavor.NU_X
        return self.file['ALPHA_{}'.format(flavor.name.upper())]

    def get_EOS(self):
        """Model equation of state.

        Returns
        -------
        eos : str
            Model equation of state.
        """
        return self.eos

    def get_progenitor_mass(self):
        """Progenitor mass.

        Returns
        -------
        mass : float
            Progenitor mass, in units of solar mass. This strips other model information, use filename for full model
        """
        return re.sub("[^0123456789\.]","",filename)

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
        for flavor in Flavor:
            L = self.luminosity[flavor](t)
           
            Ea = self.meanE[flavor](t)          # <E_nu(t)>
            Ea = Ea*1e6 * 1.60218e-12
            a = self.pinch[flavor](t)           # alpha_nu(t)
            E[E==0] = np.finfo(float).eps       # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra


class OConnor_2015(SupernovaModel):
    """Set up a model based on the black hole formation simulation in O'Connor (2015). 
    """
    def __init__(self, filename, FlavorTransformation, eos="LS220"):
        """Initialize model.

        Parameters
        ----------
        filename : str
        Absolute or relative path to file prefix, we add nue/nuebar/nux
        flavor_xform : FlavorTransformation
        Flavor transformation object with survival probabilities.
        eos : string
        Equation of state used in simulation
        """
        table = Table.read(filename,names=["TIME","L_NU_E","L_NU_E_BAR","L_NU_X","E_NU_E","E_NU_E_BAR","E_NU_X","RMS_NU_E","RMS_NU_E_BAR","RMS_NU_X"],format='ascii')
        table['ALPHA_NU_E'] = (2.0*table['E_NU_E']**2 - table['RMS_NU_E']**2)/(table['RMS_NU_E']**2 - table['E_NU_E']**2)
        table['ALPHA_NU_E_BAR'] = (2.0*table['E_NU_E_BAR']**2 - table['RMS_NU_E_BAR']**2)/(table['RMS_NU_E_BAR']**2 - table['E_NU_E_BAR']**2)
        table['ALPHA_NU_X'] = (2.0*table['E_NU_X']**2 - table['RMS_NU_X']**2)/(table['RMS_NU_X']**2 - table['E_NU_X']**2)
        table['L_NU_X'] /= 4.0
        self.file = table
        self.filename = filename
        self.shortname = "OConnor2015_s40WH07_LS220"
        self.eos = eos
        self.luminosity={}
        self.meanE={}
        self.pinch={}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d( self.get_time() , self.get_luminosity(flavor) )
            self.meanE[flavor] = interp1d( self.get_time() , self.get_mean_energy(flavor) )
            self.pinch[flavor] = interp1d( self.get_time() , self.get_pinch_param(flavor) )
        self.FT=FlavorTransformation

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.file['TIME']

    def get_luminosity(self, flavor):
        """Get model luminosity L_nu.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        luminosity : ndarray
            Grid of luminosity values (erg/s) for this flavor.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['L_{}'.format(flavor.name.upper())]

    def get_mean_energy(self, flavor):
        """Get model mean energy <E_nu>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['E_{}'.format(flavor.name.upper())]

    def get_rootmeansq_energy(self, flavor):
        """Get model mean squared energy <E_nu^2>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean squared energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['RMS_{}'.format(flavor.name.upper())]

    def get_pinch_param(self, flavor):
        """Get spectral pinch parameter alpha(t).

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        alpha : ndarray
            Grid of alpha versus time.
        """
        if (flavor == Flavor.NU_X_BAR):
            flavor = Flavor.NU_X
        return self.file['ALPHA_{}'.format(flavor.name.upper())]

    def get_EOS(self):
        """Model equation of state.

        Returns
        -------
        eos : str
            Model equation of state.
        """
        return self.eos

    def get_progenitor_mass(self):
        """Progenitor mass.

        Returns
        -------
        mass : float
            Progenitor mass, in units of solar mass. This strips other model information, use filename for full model
        """
        return 40.0

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
        for flavor in Flavor:
            L = self.luminosity[flavor](t)
            Ea = self.meanE[flavor](t)          # <E_nu(t)>
            Ea = Ea*1e6 * 1.60218e-12
            a = self.pinch[flavor](t)           # alpha_nu(t)
            E[E==0] = np.finfo(float).eps       # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra


class Warren_2020(SupernovaModel):
    "Set up a model based on simulations from Warren et al. (2020)"

    def __init__(self, filename, flavor_xform, eos = 'LS220'):
        """Initialize model.

            Parameters
            ----------
            filename : str
                Absolute or relative path to file prefix, we add nue/nuebar/nux
            flavor_xform : FlavorTransformation
                Flavor transformation object with survival probabilities.
            eos : string
                Equation of state used in simulation
        """
        f = h5py.File(filename, 'r')
        table = Table()

        for i in range(len(f['nue_data']['lum'])):
            if f['sim_data']['shock_radius'][i][1] > 0.00001:
                bounce = f['sim_data']['shock_radius'][i][0]
                break

        table['TIME'] = f['nue_data']['lum'][:, 0] - bounce
        table['L_NU_E'] = f['nue_data']['lum'][:, 1] * 1e51
        table['L_NU_E_BAR'] = f['nuae_data']['lum'][:, 1] * 1e51
        table['L_NU_X'] = f['nux_data']['lum'][:, 1] * 1e51 / 4.0
        table['E_NU_E'] = f['nue_data']['avg_energy'][:, 1]
        table['E_NU_E_BAR'] = f['nuae_data']['avg_energy'][:, 1]
        table['E_NU_X'] = f['nux_data']['avg_energy'][:, 1]
        table['RMS_NU_E'] = f['nue_data']['rms_energy'][:, 1]
        table['RMS_NU_E_BAR'] = f['nuae_data']['rms_energy'][:, 1]
        table['RMS_NU_X'] = f['nux_data']['rms_energy'][:, 1]
        table['ALPHA_NU_E'] = (2.0 * table['E_NU_E'] ** 2 - table['RMS_NU_E'] ** 2) / (
                    table['RMS_NU_E'] ** 2 - table['E_NU_E'] ** 2)
        table['ALPHA_NU_E_BAR'] = (2.0 * table['E_NU_E_BAR'] ** 2 - table['RMS_NU_E_BAR'] ** 2) / (
                    table['RMS_NU_E_BAR'] ** 2 - table['E_NU_E_BAR'] ** 2)
        table['ALPHA_NU_X'] = (2.0 * table['E_NU_X'] ** 2 - table['RMS_NU_X'] ** 2) / (
                    table['RMS_NU_X'] ** 2 - table['E_NU_X'] ** 2)
        self.file = table
        self.filename = filename
        self.eos = eos
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d(self.get_time(), self.get_luminosity(flavor))
            self.meanE[flavor] = interp1d(self.get_time(), self.get_mean_energy(flavor))
            self.pinch[flavor] = interp1d(self.get_time(), self.get_pinch_param(flavor))
        self.FT = flavor_xform

    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.file['TIME']

    def get_luminosity(self, flavor):
        """Get model luminosity L_nu.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        luminosity : ndarray
            Grid of luminosity values (erg/s) for this flavor.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['L_{}'.format(flavor.name.upper())]

    def get_mean_energy(self, flavor):
        """Get model mean energy <E_nu>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['E_{}'.format(flavor.name.upper())]

    def get_rootmeansq_energy(self, flavor):
        """Get model mean squared energy <E_nu^2>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean squared energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['RMS_{}'.format(flavor.name.upper())]

    def get_pinch_param(self, flavor):
        """Get spectral pinch parameter alpha(t).

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        alpha : ndarray
            Grid of alpha versus time.
        """
        if (flavor == Flavor.NU_X_BAR):
            flavor = Flavor.NU_X
        return self.file['ALPHA_{}'.format(flavor.name.upper())]

    def get_EOS(self):
        """Model equation of state.

        Returns
        -------
        eos : str
            Model equation of state.
        """
        return self.eos

    def get_progenitor_mass(self):
        """Progenitor mass.

        Returns
        -------
        mass : float
            Progenitor mass, in units of solar mass. This strips other model information, use filename for full model
        """
        return float(self.filename.split('_')[3].strip('m').strip('.h5'))

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
        for flavor in Flavor:
            L = self.luminosity[flavor](t)
          
            Ea = self.meanE[flavor](t)  # <E_nu(t)>
            Ea = Ea*1e6 * 1.60218e-12
            a = self.pinch[flavor](t)  # alpha_nu(t)
            E[E == 0] = np.finfo(float).eps  # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2 + a) * np.log(Ea) + (1 + a) * np.log(1 + a)
                       - loggamma(1 + a) + a * np.log(E) - (1 + a) * (E / Ea))

        return initialspectra

class Janka(SupernovaModel):
    """Set up a model based on simulations from Janka, I'll have to update this descriptioin later because I dont know where this is from
    """

    def __init__(self, filename, flavor_xform):
        """Initialize model.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        flavor_xform : FlavorTransformation
            Flavor transformation object with survival probabilities.
        """
        self.file = Table.read(filename)
        self.filename = filename
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d(self.get_time(), self.get_luminosity(flavor))
            self.meanE[flavor] = interp1d(self.get_time(), self.get_mean_energy(flavor))
            self.pinch[flavor] = interp1d(self.get_time(), self.get_pinch_param(flavor))
        self.FT = flavor_xform
            
    def get_time(self):
        """Get grid of model times.

        Returns
        -------
        time : ndarray
            Grid of times used in the model.
        """
        return self.file['TIME']
    
    def get_luminosity(self, flavor):
        """Get model luminosity L_nu.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        luminosity : ndarray
            Grid of luminosity values (erg/s) for this flavor.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['L_{}'.format(flavor.name.upper())]
        
    def get_mean_energy(self, flavor):
        """Get model mean energy <E_nu>.

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        energy : ndarray
            Grid of mean energy versus time.
        """
        if flavor == Flavor.NU_X_BAR:
            flavor = Flavor.NU_X
        return self.file['E_{}'.format(flavor.name.upper())]
    
    def get_pinch_param(self, flavor):
        """Get spectral pinch parameter alpha(t).

        Parameters
        ----------
        flavor : Flavor
            Neutrino flavor type.

        Returns
        -------
        alpha : ndarray
            Grid of alpha versus time.
        """
        if (flavor == Flavor.NU_X_BAR):
            flavor = Flavor.NU_X
        return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
    def get_EOS(self):
        """Model equation of state.

        Returns
        -------
        eos : str
            Model equation of state.
        """
        return self.filename.split('-')[1]
    
    def get_progenitor_mass(self):
        """Progenitor mass.

        Returns
        -------
        mass : float
            Progenitor mass, in units of solar mass.
        """
        return float(self.split('-')[-1].split('.')[0].strip('s'))

    def get_initialspectra(self,t,E):
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
            Dictionary of model spectra, keyed by neutrino flavor."""
    
        initialspectra = {}
        for flavor in Flavor:
            L = self.luminosity[flavor](t)
            Ea = self.meanE[flavor](t)          # <E_nu(t)>
            Ea = Ea*1e6 * 1.60218e-12
            a = self.pinch[flavor](t)           # alpha_nu(t)
            E[E==0] = np.finfo(float).eps       # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a) 
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra

'''
# class Fornax2019(SupernovaModel):
    
#     def __init__(self, filename):
#         self.file = Table.read(filename)
#         self.filename = filename
        
#     def get_time(self):
#         return self.file['TIME']
    
#     def get_luminosity(self, flavor):
#         if flavor == Flavor.NU_X_BAR:
#             flavor = Flavor.NU_X
#         return self.file['L_{}'.format(flavor.name.upper())]
        
#     def get_mean_energy(self, flavor):
#         if flavor == Flavor.NU_X_BAR:
#             flavor = Flavor.NU_X
#         return self.file['E_{}'.format(flavor.name.upper())]
    
#     def get_pinch_param(self, flavor):
#         if (flavor == Flavor.NU_X_BAR):
#             flavor = Flavor.NU_X
#         return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
#     def get_EOS(self):
#         return self.filename.split('-')[1]
    
#     def get_progenitor_mass(self):
#         return float(self.split('-')[-1].split('.')[0].strip('s'))

class OConnor_2013(SupernovaModel):
       
    eos = 'LS220'
    mass = 30
    def __init__(self, filename, FlavorTransformation):
        spectra_tuple = get_spectra(mass, eos)
        self.eos = eos
        self.luminosity = get_luminosity
        self.meanE = spectra_list(2)
        self.luminosity = get_luminosity(mass, eos)
        self.
    
    def get_luminosity(mass, eos):
        # Open luminosity file.
        tf = tarfile.open('{}_timeseries.tar.gz'.format(eos))
    
        # Extract luminosity data.
        dataname = 's{:d}_{}_timeseries.dat'.format(mass, eos)
        datafile = tf.extractfile(dataname)
        lumdata = ascii.read(datafile, names=['t', 'Le', 'Lae', 'Lx',
                                              'Ee_avg', 'Eae_avg', 'Ex_avg',
                                              'Ee_rms', 'Eae_rms', 'Ex_rms'])
        return lumdata
    
    
    def get_spectra(mass, eos):
        # Open spectra file.
        tf = tarfile.open('{}_timeseries_spectra.tar.gz'.format(eos))
    
        # Extract luminosity data.
        dataname = 's{:d}_{}_timeseries_spectra.dat'.format(mass, eos)
        datafile = tf.extractfile(dataname)
    
        t = []
        E = []
        spectra = []
    
        for line in datafile:
            tokens = [float(x) for x in line.strip().split()]
            if len(tokens) == 1:
                if t:
                    E = _E
                    spectra.append(Table([_Fe, _Fae, _Fx],
                                         names=['Fe', 'Fae', 'Fx'],
                                         meta={'t' : t[-1]}))
                t.append(tokens[0])
                _E, _Fe, _Fae, _Fx = [[] for _ in range(4)]
            elif len(tokens) == 4:
                _E.append(tokens[0])
                _Fe.append(tokens[1])
                _Fae.append(tokens[2])
                _Fx.append(tokens[3])
    
        spectra.append(Table([_Fe, _Fae, _Fx],
                             names=['Fe', 'Fae', 'Fx'],
                             meta={'t' : t[-1]}))
    
        return t, E, spectra
'''

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
