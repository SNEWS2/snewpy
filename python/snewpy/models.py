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
from astropy.table import Table

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
import math

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
        pass

    def get_oscillatedspectra(self,t,E):
        initialspectra = self.get_initialspectra(t, E)
        oscillatedspectra = {}
        oscillatedspectra[Flavor.nu_e] = self.FT.p() * initialspectra[Flavor.nu_e] + (1. - self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_x] = (1. - self.FT.p()) * initialspectra[Flavor.nu_e] + (1. + self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_e_bar] = self.FT.pbar() * initialspectra[Flavor.nu_e_bar] + (1. - self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        oscillatedspectra[Flavor.nu_x_bar] = (1. - self.FT.p()) * initialspectra[Flavor.nu_e_bar] + (1. + self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        return oscillatedspectra    


class Nakazato2013(SupernovaModel):
    
    def __init__(self, filename, FlavorTransformation):
        self.file = Table.read(filename)
        self.filename = filename
        self.luminosity={}
        self.meanE={}
        self.pinch={}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d( self.get_time() , self.get_luminosity(flavor) )
            self.meanE[flavor] = interp1d( self.get_time() , self.get_mean_energy(flavor) )
            self.pinch[flavor] = interp1d( self.get_time() , self.get_pinch_param(flavor) )
        self.FT=FlavorTransformation
        
    def get_time(self):
        return self.file['TIME']
    
    def get_luminosity(self, flavor):
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
        return self.file['L_{}'.format(flavor.name.upper())]
        
    def get_mean_energy(self, flavor):
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
        return self.file['E_{}'.format(flavor.name.upper())]
    
    def get_pinch_param(self, flavor):
        if (flavor == Flavor.nu_x_bar):
            flavor = Flavor.nu_x
        return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
    def get_EOS(self):
        return self.filename.split('-')[1].upper()
    
    def get_progenitor_mass(self):
        return float(self.filename.split('-')[-1].strip('s%.0.fits'))
    
    def get_revival_time(self):
        return float(self.filename.split('-')[-2].strip('t_rev%ms'))

    def get_initialspectra(self,t,E):
        initialspectra={}
        for flavor in Flavor:
            L=self.luminosity[flavor](t)
            ME = self.meanE[flavor](t)
            ME = ME*1e6 * 1.60218e-12            
            alpha=self.pinch[flavor](t)       
            initialspectra[flavor] = L / ME * np.float_power(alpha+1.,alpha+1.)/ME/math.gamma(alpha+1.)*np.float_power(E/ME,alpha)*np.exp(-(alpha+1.)*E/ME)
        return initialspectra


class Sukhbold2015(SupernovaModel):
    
    def __init__(self, filename, FlavorTransformation):
        self.file = Table.read(filename)
        self.filename = filename
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}
        for flavor in Flavor:
            self.luminosity[flavor] = interp1d(self.get_time(), self.get_luminosity(flavor))
            self.meanE[flavor] = interp1d(self.get_time(), self.get_mean_energy(flavor))
            self.pinch[flavor] = interp1d(self.get_time(), self.get_pinch_param(flavor))
        self.FT = FlavorTransformation
            
    def get_time(self):
        return self.file['TIME']
    
    def get_luminosity(self, flavor):
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
        return self.file['L_{}'.format(flavor.name.upper())]
        
    def get_mean_energy(self, flavor):
        if flavor == Flavor.nu_x_bar:
            flavor = Flavor.nu_x
        return self.file['E_{}'.format(flavor.name.upper())]
    
    def get_pinch_param(self, flavor):
        if (flavor == Flavor.nu_x_bar):
            flavor = Flavor.nu_x
        return self.file['ALPHA_{}'.format(flavor.name.upper())]
    
    def get_EOS(self):
        return self.filename.split('-')[1]
    
    def get_progenitor_mass(self):
        return float(self.split('-')[-1].split('.')[0].strip('s'))

    def get_initialspectra(self,t,E):
        initialspectra={}
        for flavor in Flavor:
            L = self.luminosity[flavor](t)
            ME = self.meanE[flavor](t)
            ME = ME*1e6 * astropy.units.eV            
            alpha = self.pinch[flavor](t)                
            initialspectra[flavor] = L / ME * np.float_power(alpha+1.,alpha+1.)/ME/math.gamma(alpha+1.)*np.float_power(E/ME,alpha)*np.exp(-(alpha+1.)*E/ME)
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
