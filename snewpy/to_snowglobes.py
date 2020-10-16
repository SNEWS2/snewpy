# -*- coding: utf-8 -*-
"""
snewpy.scripts.to_snowglobes
============================

Convert an arbitrary model to SNOwGLoBES format. Based on SNEWPY.py script by
E. O'Connor and J. P. Kneller. Improved by Segev Benzvi. 
Arkin Worlikar added the neutrino decay classes

This version will subsample the times in a supernova model, produce energy
tables expected by SNOwGLoBES, and compress the output into a tarfile.
"""

import numpy as np

import os
import io
import tarfile

import logging

from models import *
from flavor_transformation import *

def generate(model_path, model_file, model_type, transformation_type, transformation_parameters, d, output_filename, ntbins, deltat):

    # Chooses model format. model_format_dict associates the model format name with it's class
    model_class_dict = {'Nakazato_2013':Nakazato_2013, 'Sukhbold_2015':Sukhbold_2015, 'Bollig_2016':Bollig_2016, 'OConnor_2015':OConnor_2015, 'Janka':Janka, 'Warren_2020':Warren_2020}
    model_class = model_class_dict[model_type]
    
    # chooses flavor transformation, works in the same way as model_format
    flavor_transformation_dict = {'NoTransformation':NoTransformation, 'AdiabaticMSW_NMO':AdiabaticMSW_NMO, 'AdiabaticMSW_IMO':AdiabaticMSW_IMO, 'NonAdiabaticMSWH_NMO':NonAdiabaticMSWH_NMO, 'NonAdiabaticMSWH_IMO':NonAdiabaticMSWH_IMO, 'TwoFlavorDecoherence':TwoFlavorDecoherence, 'ThreeFlavorDecoherence':ThreeFlavorDecoherence, 'NeutrinoDecay_NMO':NeutrinoDecay_NMO, 'NeutrinoDecay_IMO':NeutrinoDecay_IMO}
    flavor_transformation = flavor_transformation_dict[transformation_type]

    snmodel = model_class(model_path+model_file, flavor_transformation(transformation_parameters))

    # Subsample the model time. Default to 30 time slices.
    tmin = snmodel.get_time()[0]
    tmax = snmodel.get_time()[-1]
    if deltat is not None:
        dt = deltat
    elif ntbins is not None:
        dt = (tmax - tmin) / (ntbins+1)
    else:
        ntbins=30
        dt = (tmax - tmin) / (ntbins+1)

    tedges = np.arange(tmin, tmax, dt)
    times = 0.5*(tedges[1:] + tedges[:-1])
    
    # Generate output.
    if output_filename is not None:
        tfname = output_filename
    elif '.fits' in model_file:
        tfname = model_file.replace('.fits', '.' + transformation_type + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(tmin,tmax,ntbins,d) + 'kpc.tar.bz2')
    elif '.dat' in model_file:
        tfname = model_file.replace('.dat', '.' + transformation_type  + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(tmin,tmax,ntbins,d) + 'kpc.tar.bz2')
    elif '.h5' in model_file:
        tfname = model_file.replace('.h5', '.' + transformation_type  + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(tmin,tmax,ntbins,d) + 'kpc.tar.bz2')
    else:
        tfname = model_file + '.' + transformation_type  + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(tmin,tmax,ntbins,d) + 'kpc.tar.bz2'

    with tarfile.open(model_path + tfname, 'w:bz2') as tf:
        #creates file in tar archive that gives information on parameters
        output = '\n'.join(map(str,transformation_parameters)).encode('ascii')
        tf.addfile(tarfile.TarInfo(name='parameterinfo'), io.BytesIO(output))        
        
        keV = 1e3 * 1.60218e-12        # eV to erg
        MeV = 1e6 * 1.60218e-12
        GeV = 1e9 * 1.60218e-12

        energy = np.linspace(0, 100, 501) * MeV
        
        # Loop over sampled times.
        for i, t in enumerate(times):
            osc_spectra = snmodel.get_oscillatedspectra(t, energy)
            
            osc_fluence = {}
            table = []

            table.append('# TBinMid={:g}sec TBinWidth={:g}s EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2'.format(t, dt))
            table.append('# E(GeV)	NuE	NuMu	NuTau	aNuE	aNuMu	aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 200.*keV  / (4.*np.pi*(d*1000*3.086e+18)**2)

                s = '{:17.8E}'.format(E/GeV)
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_e])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_e_bar])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x_bar])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x_bar])
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            extension = ".dat"
            filename = model_file.replace('.fits', '.tbin{:01d}.'.format(i+1) + transformation_type + '.{:.3f},{:.3f},{:01d}-{:.1f}kpc{}'.format(tmin,tmax,ntbins,d,extension) )
            info = tarfile.TarInfo(name=filename)
            info.size = len(output)
            tf.addfile(info, io.BytesIO(output))

    return tfname






