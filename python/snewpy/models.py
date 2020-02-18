# -*- coding: utf-8 -*-
"""
A submodule with classes used for supernova model files stored on disk. Based
on the models_class.py script started by Navya Uberoi and updated by James
Kneller.

The changes, made by James Kneller, are:

#. The Enum class inherited by the Flavor class was changed to IntEnum and the
assigned values altered so that the flavor enums could be used as list indicii.
The member functions should still work as before.
#. Extra member data was added to the SupernovaModel class to store the
interpolations of the luminosity, mean energy and pinch parameter as functions
of time.
#. The SupernovaModel class now requires instantation with a
FlavorTransofrmation class.  The FlavorTransofrmation class - defined in
FlavorTransformation.py - has two members which return the electron neutrino
survival probability p and electron antineutrino survival probability pbar.
#. Extra methods were added to the SupernovaModel class as its descendents to
return the initial and oscillated spectra at a given time and neutrino energy.
"""

from abc import abstractmethod, ABC
from enum import IntEnum

import astropy
from astropy.table import Table, join

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
from scipy.special import loggamma

import logging
import re

import tarfile

from snewpy.FlavorTransformation import *

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

class Flavor(IntEnum):
    """Enumeration of CCSN Neutrino flavors.
    """
    nu_e = 2
    nu_e_bar = 1
    nu_x = 3
    nu_x_bar = 0
    
    def to_tex(self):
        """LaTeX-compatible string representations of flavor.
        """
        
        if '_bar' in self.name:
            return r'$\overline{{\nu}}_{0}$'.format(self.name[3])
        return r'$\{0}$'.format(self.name)

    @property
    def is_electron(self):
        return self.value in (Flavor.nu_e.value, Flavor.nu_e_bar.value)

    @property
    def is_neutrino(self):
        return self.value in (Flavor.nu_e.value, Flavor.nu_x.value)

    @property
    def is_antineutrino(self):
        return self.value in (Flavor.nu_e_bar.value, Flavor.nu_x_bar.value)


class SupernovaModel(ABC):
    """Base class defining an interface to a supernova model."""
    
    def __init__(self):
        pass
    
    @abstractmethod
    def get_time(self):
        pass
    
    @abstractmethod
    def get_mean_energy(self, flavor):
        pass
    
    @abstractmethod
    def get_luminosity(self, flavor):
        pass
    
    @abstractmethod
    def get_pinch_param(self, flavor):
        pass
    
    @abstractmethod
    def get_EOS(self):
        pass
    
    @abstractmethod
    def get_progenitor_mass(self):
        pass

    @abstractmethod
    def get_initialspectra(self, t, E):
        """Get neutrino spectra/luminosity curves at the source.

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

    def get_oscillatedspectra(self, t, E):
        """Get neutrino spectra/luminosity curves after oscillation.

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
        oscillatedspectra[Flavor.nu_e] = \
            self.FT.p() * initialspectra[Flavor.nu_e] + \
            (1. - self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_x] = \
            (1. - self.FT.p()) * initialspectra[Flavor.nu_e] + \
            (1. + self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_e_bar] = \
            self.FT.pbar() * initialspectra[Flavor.nu_e_bar] + \
            (1. - self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        oscillatedspectra[Flavor.nu_x_bar] = \
            (1. - self.FT.p()) * initialspectra[Flavor.nu_e_bar] + \
            (1. + self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        return oscillatedspectra    


class Nakazato2013(SupernovaModel):
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
        self.file = Table.read(filename)
        self.filename = filename
        self.luminosity={}
        self.meanE={}
        self.pinch={}
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if (flavor == Flavor.nu_x_bar):
            flavor = Flavor.nu_x
        return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
    def get_EOS(self):
        """Model equation of state.

        Returns
        -------
        eos : str
            Model equation of state.
        """
        return self.filename.split('-')[1].upper()
    
    def get_progenitor_mass(self):
        """Progenitor mass.

        Returns
        -------
        mass : float
            Progenitor mass, in units of solar mass.
        """
        return float(self.filename.split('-')[-1].strip('s%.0.fits'))
    
    def get_revival_time(self):
        """Revival time of model explosion; specific to the Nakazato models.

        Returns
        -------
        time : float
            Revival time, in ms.
        """
        return float(self.filename.split('-')[-2].strip('t_rev%ms'))

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


class Sukhbold2015(SupernovaModel):
    """Set up a model based on simulations from Sukhbold et al., ApJ 821:38,
    2016. Models were shared privately by email.
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if (flavor == Flavor.nu_x_bar):
            flavor = Flavor.nu_x
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
            Ea = Ea*1e6
            a = self.pinch[flavor](t)           # alpha_nu(t)
            E[E==0] = np.finfo(float).eps       # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a) 
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra

class Bollig2016(SupernovaModel):
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if (flavor == Flavor.nu_x_bar):
            flavor = Flavor.nu_x
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
            L=L*6.242e11
            Ea = self.meanE[flavor](t)          # <E_nu(t)>
            Ea = Ea*1e6
            a = self.pinch[flavor](t)           # alpha_nu(t)
            E[E==0] = np.finfo(float).eps       # Avoid division by zero.

            # For numerical stability, evaluate log PDF then exponentiate.
            initialspectra[flavor] = \
                np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                       - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea))

        return initialspectra

class OConnor2015(SupernovaModel):
    """Set up a model based on the black hole formation simulation in
    O'Connor (2015). 
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
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
        if (flavor == Flavor.nu_x_bar):
            flavor = Flavor.nu_x
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
            L=L*6.242e11
            Ea = self.meanE[flavor](t)          # <E_nu(t)>
            Ea = Ea*1e6
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
#         if flavor == Flavor.nu_x_bar:
#             flavor = Flavor.nu_x
#         return self.file['L_{}'.format(flavor.name.upper())]
        
#     def get_mean_energy(self, flavor):
#         if flavor == Flavor.nu_x_bar:
#             flavor = Flavor.nu_x
#         return self.file['E_{}'.format(flavor.name.upper())]
    
#     def get_pinch_param(self, flavor):
#         if (flavor == Flavor.nu_x_bar):
#             flavor = Flavor.nu_x
#         return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
#     def get_EOS(self):
#         return self.filename.split('-')[1]
    
#     def get_progenitor_mass(self):
#         return float(self.split('-')[-1].split('.')[0].strip('s'))
'''

class SNOwGLoBES:
    """A model that does not inherit from SupernovaModel (yet) and imports a
    group of SNOwGLoBES files."""

    def __init__(self, tarfilename):
        """Initialize model from a tar archive.

        Parameters
        ----------
        tarfilename: str
            Absolute or relative path to tar archive with SNOwGLoBES files.
        """
        self.tfname = tarfilename
        tf = tarfile.open(self.tfname)
        logging.debug('Opening {}'.format(tarfilename))

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
