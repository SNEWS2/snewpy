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

from snewpy.models import *
from .neutrino import MassHierarchy
from snewpy.flavor_transformation import *
from astropy import units as u


def generate_output_name(output_filename, model_file,
                         transformation_type, tmin,
                         tmax, ntbins, d):
    if output_filename is not None:
        tfname = output_filename + 'kpc.tar.bz2'
    elif '.fits' in model_file:
        tfname = model_file.replace('.fits', f'.{transformation_type}{tmin:.3f},{tmax:.3f},{ntbins:d}-{d:.1f}kpc.tar.bz2')
    elif '.dat' in model_file:
        tfname = model_file.replace('.dat', f'.{transformation_type}.{tmin:.3f},{tmax:.3f},{ntbins:d}-{d:.1f}kpc.tar.bz2')
    elif '.h5' in model_file:
        tfname = model_file.replace('.h5', f'.{transformation_type}.{tmin:.3f},{tmax:.3f},{ntbins:d}-{d:.1f}kpc.tar.bz2')
    else:
        tfname = model_file + f'.{transformation_type}.{tmin:.3f},{tmax:.3f},{ntbins:d}-{d:.1f}kpc.tar.bz2'
    return tfname 

def get_model(model_type):
    """ Chooses model format. model_format_dict associates
         the model format name with it's class
    """
    model_class_dict = {'Nakazato_2013':Nakazato_2013, 
                        'Sukhbold_2015':Sukhbold_2015, 
                        'Bollig_2016':Bollig_2016, 
                        'OConnor_2015':OConnor_2015, 
                        'Fornax_2021':Fornax_2021_2D, 
                        'Warren_2020':Warren_2020, 
                        'Analytic3Species':Analytic3Species, 
                        'Zha_2021':Zha_2021, 
                        'Tamborra_2014':Tamborra_2014, 
                        'Walk_2018':Walk_2018, 
                        'Walk_2019':Walk_2019}
    model_class = model_class_dict[model_type]
    return model_class

def get_flavor(transformation_type):
    """ Chooses flavor transformation, 
        works in the same way as model_format
    """
    flavor_transformation_dict = {'NoTransformation':NoTransformation(), 
                                  'AdiabaticMSW_NMO':AdiabaticMSW(mh=MassHierarchy.NORMAL), 
                                  'AdiabaticMSW_IMO':AdiabaticMSW(mh=MassHierarchy.INVERTED), 
                                  'NonAdiabaticMSWH_NMO':NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 
                                  'NonAdiabaticMSWH_IMO':NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 
                                  'TwoFlavorDecoherence':TwoFlavorDecoherence(), 
                                  'ThreeFlavorDecoherence':ThreeFlavorDecoherence(), 
                                  'NeutrinoDecay_NMO':NeutrinoDecay(mh=MassHierarchy.NORMAL), 
                                  'NeutrinoDecay_IMO':NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    flavor_transformation = flavor_transformation_dict[transformation_type]
    return flavor_transformation

def generate_time_series(model_path, 
                         model_file, 
                         model_type, 
                         transformation_type, 
                         transformation_parameters, 
                         d, 
                         output_filename, 
                         ntbins = 30, 
                         deltat = None):
    """
    model_path : str
        Location of the SN models
    model_file : str
        Name of the model file at model_path
    model_type : str
        The model to be used from snewpy.models, the options are;
        Nakazato_2013, Sukhbold_2015, Bollig_2016, OConnor_2015, Fornax_2021, 
        Warren_2020, Analytic3Species, Zha_2021, Tamborra_2014, Walk_2018, Walk_2019
    transformation_type : str
        Transformation to be used select from
        NoTransformation, AdiabaticMSW_NMO, AdiabaticMSW_IMO, 
        NonAdiabaticMSWH_NMO, NonAdiabaticMSWH_IMO, TwoFlavorDecoherence, 
        ThreeFlavorDecoherence, NeutrinoDecay_NMO, NeutrinoDecay_IMO
    transformation_parameters : 
    d : float
        distance to supernova in kpc
    output_filename : str
        output file name
    ntbins: int, optional
        desired number of time bins
    deltat : float, optional
        the time step to divide timeseries. If given number of bins is recalculated
        and ntbins is overwritten.
    """
    model_class = get_model(model_type)
    flavor_transformation = get_flavor(transformation_type)

    # sets the SN model obj
    snmodel = model_class(f"{model_path}/{model_file}")

    # Subsample the model time. Default to 30 time slices.
    # snmodel.get_time() brings the timeseries from the model
    tmin = snmodel.get_time()[0]
    tmax = snmodel.get_time()[-1]
    if deltat is not None:                     # if delta-t is given, 
        dt = deltat                      
        ntbins = len(np.arange(tmin,tmax,dt))  # calculate nt-bins and overwrite
    else:                                       # if delta-t is not given, use default nt-bins
        dt = (tmax - tmin) / (ntbins+1)

    tedges = np.arange(tmin, tmax, dt)
    times = 0.5*(tedges[1:] + tedges[:-1])
    
    # Generate output.
    tfname = generate_output_name(output_filename, model_file,
                                  transformation_type, tmin,
                                  tmax, ntbins, d)

    with tarfile.open(model_path + tfname, 'w:bz2') as tf:
        # creates file in tar archive that gives information on parameters
        output = '\n'.join(map(str,transformation_type)).encode('ascii')
        tf.addfile(tarfile.TarInfo(name='parameterinfo'), io.BytesIO(output))        

        MeV = 1.60218e-6 * u.erg
        energy = np.linspace(0, 100, 501) * MeV #1MeV
        
        # Loop over sampled times.
        for i, t in enumerate(times):
            osc_spectra = snmodel.get_oscillatedspectra(t, energy,flavor_transformation)
            
            osc_fluence = {}
            table = []

            table.append(f'# TBinMid={t:g}sec TBinWidth={dt:g}s EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2')
            table.append('# E(GeV)	NuE	NuMu	NuTau	aNuE	aNuMu	aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 0.2 * MeV / (4.*np.pi*(d*1000*3.086e+18)**2)


                s = f'{E/(1e3 * MeV):17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_E]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_E_BAR]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X_BAR]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X_BAR]:17.8E}'
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            extension = ".dat"
            filename = model_file.replace('.fits', f'.tbin{i+1:01d}.{transformation_type}.{tmin:.3f},{tmax:.3f},{ntbins:01d}-{d:.1f}kpc{extension}')
            info = tarfile.TarInfo(name=filename)
            info.size = len(output)
            tf.addfile(info, io.BytesIO(output))
    return tfname

def generate_fluence(model_path, model_file, model_type, transformation_type, d, output_filename, tstart=None, tend=None):
    model_class = get_model(model_type)
    flavor_transformation = get_flavor(transformation_type)
    snmodel = model_class(model_path+"/"+model_file)

    #set the timings up
    #default if inputs are None: full time window of the model
    if tstart is None:
        tstart = snmodel.get_time()[0]
        tend = snmodel.get_time()[-1]

    try:
        if len(tstart/u.s) > 0:
            t0 = tstart[0]
            t1 = tend[-1]
            nbin = len(tstart/u.s)
    except:
        t0 = tstart
        t1 = tend
        nbin = 1

    times = 0.5*(tstart + tend)

    model_times = snmodel.get_time()
    model_tstart = model_times*1.0
    model_tend = model_times*1.0

    model_tstart[0] = model_times[0]
    for i in range(1,len(model_times),1):
        model_tstart[i] = 0.5*(model_times[i]+model_times[i-1])
        model_tend[i-1] = model_tstart[i]
    model_tend[len(model_times)-1] = model_times[-1]
    
    if nbin>1:
        starting_index = np.zeros(len(times),dtype=np.int64)
        ending_index = np.zeros(len(times),dtype=np.int64)
        for i in range(len(tstart)):
            starting_index[i] = next(j for j, t in enumerate(model_tend) if t>tstart[i])
            ending_index[i] = next(j for j, t in enumerate(model_tend) if t>=tend[i])
    else:
        starting_index = [next(j for j, t in enumerate(model_tend) if t>tstart)]
        ending_index = [next(j for j, t in enumerate(model_tend) if t>=tend)]

    
    # Generate output.
    tfname = generate_output_name(output_filename, model_file,
                                  transformation_type, tmin,
                                  tmax, ntbins, d)

    with tarfile.open(model_path + tfname, 'w:bz2') as tf:
        #creates file in tar archive that gives information on parameters
        output = '\n'.join(map(str,transformation_type)).encode('ascii')
        tf.addfile(tarfile.TarInfo(name='parameterinfo'), io.BytesIO(output))

        MeV = 1.60218e-6 * u.erg
        energy = np.linspace(0, 100, 501) * MeV

        # Loop over sampled times.
        for i in range(nbin):

            if nbin > 1:
                ta = tstart[i]
                tb = tend[i]
                t = times[i]
                dt = tb-ta
            else:
                ta = tstart
                tb = tend
                t = times
                dt = tb-ta

            #first time bin of model in requested interval
            osc_spectra = snmodel.get_oscillatedspectra(model_times[starting_index[i]], energy,flavor_transformation)

            if dt < model_tend[starting_index[i]]-ta:
                dt=dt
            else:
                for flavor in Flavor:
                    osc_spectra[flavor] *= (model_tend[starting_index[i]]-ta)

                #intermediate time bins of model in requested interval
                for j in range(starting_index[i]+1,ending_index[i],1):
                    temp_spectra = snmodel.get_oscillatedspectra(model_times[j], energy,flavor_transformation)
                    for flavor in Flavor:
                        osc_spectra[flavor] += temp_spectra[flavor]*(model_tend[j]-model_tstart[j])

                #last time bin of model in requested interval
                temp_spectra = snmodel.get_oscillatedspectra(model_times[ending_index[i]], energy,flavor_transformation)
                for flavor in Flavor:
                    osc_spectra[flavor] += temp_spectra[flavor]*(tb-model_tstart[ending_index[i]])

                for flavor in Flavor:
                    osc_spectra[flavor] /= (tb-ta)

            osc_fluence = {}
            table = []

            table.append(f'# TBinMid={t:g}sec TBinWidth={dt:g}s EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2')
            table.append('# E(GeV)	NuE	NuMu	NuTau	aNuE	aNuMu	aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 0.2 * MeV  / (4.*np.pi*(d*1000*3.086e+18)**2)

                s = '{E/(1e3 * MeV):17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_E]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_E_BAR]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X_BAR]:17.8E}'
                s = f'{s}{osc_fluence[Flavor.NU_X_BAR]:17.8E}'
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            extension = ".dat"
            if output_filename is not None:
                if nbin>1:
                    filename = output_filename+"_"+str(i)+extension
                else:
                    filename = output_filename+extension
            elif '.fits' in model_file:
                filename = model_file.replace('.fits', f'.tbin{i+1:01d}.{transformation_type}.{t0:.3f},{t1:.3f},{nbin:01d}-{d:.1f}kpc{extension}')
            elif '.dat' in model_file:
                filename = model_file.replace('.dat', f'.tbin{i+1:01d}.{transformation_type}.{t0:.3f},{t1:.3f},{nbin:01d}-{d:.1f}kpc{extension}')
            elif '.h5' in model_file:
                filename = model_file.replace('.h5', f'.tbin{i+1:01d}.{transformation_type}.{t0:.3f},{t1:.3f},{nbin:01d}-{d:.1f}kpc{extension}')
            else:
                filename = model_file + f'.tbin{i+1:01d}.{transformation_type}.{t0:.3f},{t1:.3f},{nbin:01d}-{d:.1f}kpc{extension}'

            info = tarfile.TarInfo(name=filename)
            info.size = len(output)
            tf.addfile(info, io.BytesIO(output))

    return tfname





