
'''
  This file is part of SNEWPY.

  This file is a modified version of the code uploaded by nuberoi.  
  The changes, made by James Kneller, are:
      1) The Enum class inherited by the Flavor class was changed to IntEnum 
      and the assigned values altered so that the flavor enums could be 
      used as list indicii. The member functions should still work as before.
      2) Extra member data was added to the SupernovaModel class to store the interpolations
      of the luminosity, mean energy and pinch parameter as functions of time.
      3) The SupernovaModel class now requires instantation with a FlavorTransofrmation class.
      The FlavorTransofrmation class - defined in FlavorTransformation.py - has two members 
      which return the electron neutrino survival probability p and electron antineutrino survival 
      probability pbar.
      4) Extra methods were added to the SupernovaModel class as its descendents 
      to return the initial and oscillated spectra at a given time and neutrino energy.
'''

from abc import abstractmethod, ABC
from enum import IntEnum

import astropy
from astropy.table import Table

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
import math

from snewpy.FlavorTransformation import *

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
    def get_initialspectra(self,t,E):
        pass

    @abstractmethod
    def get_oscillatedspectra(self,t,E):
        pass


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
            ME=ME*1e6 * 1.60218e-12            
            alpha=self.pinch[flavor](t)       
            initialspectra[flavor] = L / ME * np.float_power(alpha+1.,alpha+1.)/ME/math.gamma(alpha+1.)*np.float_power(E/ME,alpha)*np.exp(-(alpha+1.)*E/ME)
        return initialspectra

    def get_oscillatedspectra(self,t,E):
        initialspectra=self.get_initialspectra(t,E)
        oscillatedspectra={}
        oscillatedspectra[Flavor.nu_e]=self.FT.p() * initialspectra[Flavor.nu_e] + (1.-self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_x]=(1.-self.FT.p()) * initialspectra[Flavor.nu_e] + (1.+self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_e_bar]=self.FT.pbar() * initialspectra[Flavor.nu_e_bar] + (1.-self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        oscillatedspectra[Flavor.nu_x_bar]=(1.-self.FT.p()) * initialspectra[Flavor.nu_e_bar] + (1.+self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        return oscillatedspectra


class Sukhbold2015(SupernovaModel):
    
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
        return self.filename.split('-')[1]
    
    def get_progenitor_mass(self):
        return float(self.split('-')[-1].split('.')[0].strip('s'))

    def get_initialspectra(self,t,E):
        initialspectra={}
        for flavor in Flavor:
            L=self.luminosity[flavor](t)
            ME = self.meanE[flavor](t)
            ME=ME*1e6 * astropy.units.eV            
            alpha=self.pinch[flavor](t)                
            initialspectra[flavor] = L / ME * np.float_power(alpha+1.,alpha+1.)/ME/math.gamma(alpha+1.)*np.float_power(E/ME,alpha)*np.exp(-(alpha+1.)*E/ME)
        return initialspectra
    
    def get_oscillatedspectra(self,t,E):
        initialspectra=self.get_initialspectra(t,E)
        oscillatedspectra={}
        oscillatedspectra[Flavor.nu_e]=self.FT.p() * initialspectra[Flavor.nu_e] + (1.-self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_x]=(1.-self.FT.p()) * initialspectra[Flavor.nu_e] + (1.+self.FT.p()) * initialspectra[Flavor.nu_x]
        oscillatedspectra[Flavor.nu_e_bar]=self.FT.pbar() * initialspectra[Flavor.nu_e_bar] + (1.-self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        oscillatedspectra[Flavor.nu_x_bar]=(1.-self.FT.p()) * initialspectra[Flavor.nu_e_bar] + (1.+self.FT.pbar()) * initialspectra[Flavor.nu_x_bar]
        return oscillatedspectra    

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

