# -*- coding: utf-8 -*-
"""The ``snewpy.snowglobes`` module contains functions for interacting with SNOwGLoBES.

`SNOwGLoBES <https://github.com/SNOwGLoBES/snowglobes>`_ can estimate detected
event rates from a given input supernova neutrino flux. It supports many
different neutrino detectors, detector materials and interaction channels.
There are three basic steps to using SNOwGLoBES from SNEWPY:

* **Generating input files for SNOwGLoBES:**
    There are two ways to do this, either generate a time series or a fluence file. This is done taking as input the supernova simulation model.
    The first will evaluate the neutrino flux at each time step, the latter will compute the integrated neutrino flux (fluence) in the time bin.
    The result is a compressed .tar file containing all individual input files.
* **Running SNOwGLoBES:**
    This step convolves the fluence generated in the previous step with the cross-sections for the interaction channels happening in various detectors supported by SNOwGLoBES.
    It takes into account the effective mass of the detector as well as a smearing matrix describing the energy-dependent detection efficiency.
    The output gives the number of events detected as a function of energy for each interaction channel, integrated in a given time window (or time bin), or in a snapshot in time.
* **Collating SNOwGLoBES outputs:**
    This step puts together all the interaction channels and time bins evaluated by SNOwGLoBES in a single file (for each detector and for each time bin).
    The output tables allow to build the detected neutrino energy spectrum and neutrino time distribution, for each reaction channel or the sum of them.
"""

from __future__ import unicode_literals

import fnmatch
import io
import logging
import os
import re
import tarfile
import zipfile
from pathlib import Path
import subprocess
import traceback
import itertools
import pdb

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

import snewpy.models
from snewpy.flavor_transformation import *
from snewpy.neutrino import Flavor, MassHierarchy

mpl.use('Agg')

logger = logging.getLogger(__name__)

def generate_time_series(model_path, model_type, transformation_type, d, output_filename=None, ntbins=30, deltat=None):
    """Generate time series files in SNOwGLoBES format.

    This version will subsample the times in a supernova model, produce energy
    tables expected by SNOwGLoBES, and compress the output into a tarfile.

    Parameters
    ----------
    model_path : str
        Input file containing neutrino flux information from supernova model.
    model_type : str
        Format of input file. Matches the name of the corresponding class in :py:mod:`snewpy.models`.
    transformation_type : str
        Name of flavor transformation. See snewpy.flavor_transformation documentation for possible values.
    d : int or float
        Distance to supernova in kpc.
    output_filename : str or None
        Name of output file. If ``None``, will be based on input file name.
    ntbins : int
        Number of time slices. Will be ignored if ``deltat`` is also given.
    deltat : astropy.Quantity or None
        Length of time slices.

    Returns
    -------
    str
        Path of compressed .tar file with neutrino flux data.
    """
    model_class = getattr(snewpy.models, model_type)

    # Choose flavor transformation. Use dict to associate the transformation name with its class.
    flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    flavor_transformation = flavor_transformation_dict[transformation_type]

    model_dir, model_file = os.path.split(os.path.abspath(model_path))
    snmodel = model_class(model_path)

    # Subsample the model time. Default to 30 time slices.
    tmin = snmodel.get_time()[0]
    tmax = snmodel.get_time()[-1]
    if deltat is not None:
        dt = deltat
        ntbins = int((tmax-tmin)/dt)
    else:
        dt = (tmax - tmin) / (ntbins+1)

    tedges = np.arange(tmin/u.s, tmax/u.s, dt/u.s)*u.s
    times = 0.5*(tedges[1:] + tedges[:-1])

    # Generate output.
    if output_filename is not None:
        tfname = output_filename + 'kpc.tar.bz2'
    else:
        model_file_root, _ = os.path.splitext(model_file)  # strip extension (if present)
        tfname = model_file_root + '.' + transformation_type + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(tmin, tmax, ntbins, d) + 'kpc.tar.bz2'

    with tarfile.open(os.path.join(model_dir, tfname), 'w:bz2') as tf:
        #creates file in tar archive that gives information on parameters
        output = '\n'.join(map(str, transformation_type)).encode('ascii')
        tf.addfile(tarfile.TarInfo(name='parameterinfo'), io.BytesIO(output))

        MeV = 1.60218e-6 * u.erg
        energy = np.linspace(0, 100, 501) * MeV  # 1MeV

        # Loop over sampled times.
        for i, t in enumerate(times):
            osc_spectra = snmodel.get_transformed_spectra(t, energy, flavor_transformation)

            osc_fluence = {}
            table = []

            table.append('# TBinMid={:g}sec TBinWidth={:g}s EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2'.format(t, dt))
            table.append('# E(GeV)	NuE	NuMu	NuTau	aNuE	aNuMu	aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 0.2 * MeV / (4.*np.pi*(d*1000*3.086e+18)**2)

                s = '{:17.8E}'.format(E/(1e3 * MeV))
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_E])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_E_BAR])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X_BAR])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X_BAR])
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            extension = ".dat"
            model_file_root, _ = os.path.splitext(model_file)
            filename = model_file_root + '.tbin{:01d}.'.format(i+1) + transformation_type + \
                '.{:.3f},{:.3f},{:01d}-{:.1f}kpc{}'.format(tmin/u.s, tmax/u.s, ntbins, d, extension)

            info = tarfile.TarInfo(name=filename)
            info.size = len(output)
            tf.addfile(info, io.BytesIO(output))

    return os.path.join(model_dir, tfname)


def generate_fluence(model_path, model_type, transformation_type, d, output_filename=None, tstart=None, tend=None):
    """Generate fluence files in SNOwGLoBES format.

    This version will subsample the times in a supernova model, produce energy
    tables expected by SNOwGLoBES, and compress the output into a tarfile.

    Parameters
    ----------
    model_path : str
        Input file containing neutrino flux information from supernova model.
    model_type : str
        Format of input file. Matches the name of the corresponding class in :py:mod:`snewpy.models`.
    transformation_type : str
        Name of flavor transformation. See snewpy.flavor_transformation documentation for possible values.
    d : int or float
        Distance to supernova in kpc.
    output_filename : str or None
        Name of output file. If ``None``, will be based on input file name.
    tstart : astropy.Quantity or None
        Start of time interval to integrate over, or list of start times of the time series bins.
    tend : astropy.Quantity or None
        End of time interval to integrate over, or list of end times of the time series bins.

    Returns
    -------
    str
        Path of compressed .tar file with neutrino flux data.
    """
    model_class = getattr(snewpy.models, model_type)

    # Choose flavor transformation. Use dict to associate the transformation name with its class.
    flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    flavor_transformation = flavor_transformation_dict[transformation_type]

    model_dir, model_file = os.path.split(os.path.abspath(model_path))
    snmodel = model_class(model_path)

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
    for i in range(1, len(model_times), 1):
        model_tstart[i] = 0.5*(model_times[i]+model_times[i-1])
        model_tend[i-1] = model_tstart[i]
    model_tend[len(model_times)-1] = model_times[-1]

    if nbin > 1:
        starting_index = np.zeros(len(times), dtype=np.int64)
        ending_index = np.zeros(len(times), dtype=np.int64)
        for i in range(len(tstart)):
            starting_index[i] = next(j for j, t in enumerate(model_tend) if t > tstart[i])
            ending_index[i] = next(j for j, t in enumerate(model_tend) if t >= tend[i])
    else:
        starting_index = [next(j for j, t in enumerate(model_tend) if t > tstart)]
        ending_index = [next(j for j, t in enumerate(model_tend) if t >= tend)]

    # Generate output.
    if output_filename is not None:
        tfname = output_filename+'.tar.bz2'
    else:
        model_file_root, _ = os.path.splitext(model_file)  # strip extension (if present)
        tfname = model_file_root + '.' + transformation_type + '.{:.3f},{:.3f},{:d}-{:.1f}'.format(t0, t1, nbin, d) + 'kpc.tar.bz2'

    with tarfile.open(os.path.join(model_dir, tfname), 'w:bz2') as tf:
        #creates file in tar archive that gives information on parameters
        output = '\n'.join(map(str, transformation_type)).encode('ascii')
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
            osc_spectra = snmodel.get_transformed_spectra(model_times[starting_index[i]], energy, flavor_transformation)

            if dt < model_tend[starting_index[i]]-ta:
                dt = dt
            else:
                for flavor in Flavor:
                    osc_spectra[flavor] *= (model_tend[starting_index[i]]-ta)

                #intermediate time bins of model in requested interval
                for j in range(starting_index[i]+1, ending_index[i], 1):
                    temp_spectra = snmodel.get_transformed_spectra(model_times[j], energy, flavor_transformation)
                    for flavor in Flavor:
                        osc_spectra[flavor] += temp_spectra[flavor]*(model_tend[j]-model_tstart[j])

                #last time bin of model in requested interval
                temp_spectra = snmodel.get_transformed_spectra(
                    model_times[ending_index[i]], energy, flavor_transformation)
                for flavor in Flavor:
                    osc_spectra[flavor] += temp_spectra[flavor]*(tb-model_tstart[ending_index[i]])

                for flavor in Flavor:
                    osc_spectra[flavor] /= (tb-ta)

            osc_fluence = {}
            table = []

            table.append('# TBinMid={:g}sec TBinWidth={:g}s EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2'.format(t, dt))
            table.append('# E(GeV)	NuE	NuMu	NuTau	aNuE	aNuMu	aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 0.2 * MeV / (4.*np.pi*(d*1000*3.086e+18)**2)

                s = '{:17.8E}'.format(E/(1e3 * MeV))
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_E])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_E_BAR])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X_BAR])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.NU_X_BAR])
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            extension = ".dat"
            if output_filename is not None:
                if nbin > 1:
                    filename = output_filename+"_"+str(i)+extension
                else:
                    filename = output_filename+extension
            else:
                model_file_root, _ = os.path.splitext(model_file)  # strip extension (if present)
                filename = model_file_root + '.tbin{:01d}.'.format(i+1) + transformation_type + \
                    '.{:.3f},{:.3f},{:01d}-{:.1f}kpc{}'.format(t0, t1, nbin, d, extension)

            info = tarfile.TarInfo(name=filename)
            info.size = len(output)
            tf.addfile(info, io.BytesIO(output))

    return os.path.join(model_dir, tfname)

det_materials = { 'icecube':'water',
                  'wc100kt30prct':'water',
                  'wc100kt15prct':'water',
                  'hyperk30prct':'water',
                  'km3net':'water',
                  'ar40kt':'argon',
                  'novaND':'nova_soup',
                  'novaFD':'nova_soup',
                  'scint20kt':'scint',
                  'halo1':'lead',
                  'halo2':'lead'
                  }

categories_map = {'water': ["nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16"],
                  'lead' : ["nc_", "*_e_", "nue_Pb208_1n", "nue_Pb208_2n"],
                  'argon': ["nc_", "*_e_", "nue_Ar40", "nuebar_Ar40"],
                  'scint': ["nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12", "nue_C13"],
                  'nova_soup':["nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12"]
                 }

def load_datafile(f):
    """ Read data from  snowglobes output file"""
    return np.loadtxt(f,comments=['---','Total','#'], dtype=[('energy','f8'),('events','f8')])

def load_channels(fname):
    t = np.loadtxt(fname, dtype=[('name','U100'),('number','u4'),('parity','U1'),('flavor','U1'),('weight','f8')])
    return t

def load_target_masses(fname):
    t = np.loadtxt(fname, dtype=[('name','U100'),('mass','f8'),('factor','f8')])
    tgt_mass = t['mass']*t['factor']
    return dict(zip(t['name'],tgt_mass))


def simulate(SNOwGLoBESdir, tarball_path, detector_input="all", verbose=False):
    """Takes as input the neutrino flux files and configures and runs the supernova script inside SNOwGLoBES, which outputs calculated event rates expected for a given (set of) detector(s). These event rates are given as a function of the neutrino energy and time, for each interaction channel.

    Parameters
    ----------
    SNOwGLoBESdir : str
        Path to directory where SNOwGLoBES is installed.
    tarball_path : str
        Path of compressed .tar file produced e.g. by ``generate_time_series()`` or ``generate_fluence()``.
    detector_input : str
        Name of detector. If ``"all"``, will use all detectors supported by SNOwGLoBES.
    verbose : bool
        Whether to generate verbose output, e.g. for debugging.
    """
    
    sng = Path(SNOwGLoBESdir)
    #prepare output folder
    output_dir = sng/'out'
    output_dir.mkdir(exist_ok=True)
    flux_dir = sng/'fluxes'
    #extract flux files
    with tarfile.open(tarball_path,'r') as f:
        flux_files = f.getnames() 
        f.extractall(path=flux_dir)
     
    flux_files = [Path(f) for f in flux_files if f.endswith('.dat')]
    logger.debug(f'Flux files = {flux_files}')
    #regexps to use
    re_flux = re.compile('flux_file\s?=.*\n')
    re_tgt  = re.compile('target_mass\s?=.*\n')

    tgt_masses = load_target_masses(sng/'detector_configurations.dat')

    def prepare_globes_file(output, flux_fname, channels, detector_name):
        """Function to write the snowglobes configuration to file """
        logger.debug(f'Prepare globes run for [flux={flux_fname}, det={detector_name}]')

        def cat(fname):
            with open(fname,'r') as f:
                return f.read()
        if flux_fname.is_file()==False:
            raise FileNotFoundError(flux_fname)
        #write preamble
        output.write(cat(sng/'glb/preamble.glb'))
        #write flux, substituted with input argument
        s = cat(sng/'glb/flux.glb')
        s = re_flux.sub(f'flux_file=   "{flux_fname.relative_to(sng)}"\n',s)
        output.write(s)

        # smearing for each channel
        for c in channels:
            output.write(f'include \"{sng}/smear/smear_{c["name"]}_{detector_name}.dat"\n')

        # Writes modified DETECTOR to GLOBESFILE; replace the mass with the calculated TargetMass
        s = cat(sng/'glb/detector.glb')
        mass = tgt_masses[detector_name] 
        s = re_tgt.sub(f'target_mass=  {mass:13.6f}\n',s)
        output.write(s)

        # Now the cross-sections.  Note that some of these are repeated even though
        # it is not necessary (xscns for several flavors can be in the same file).
        # This is just to make a consistent loop over channels.
        output.write("\n\n /******** Cross-sections *********/\n \n")
        for c in channels:
            output.write(f'cross(#{c["name"]})< @cross_file= "{sng}/xscns/xs_{c["name"]}.dat" >\n')

        output.write("\n /******** Channels *********/\n \n")
        for c in channels:
            output.write(f'channel(#{c["name"]}_signal)<\n')
            output.write('      @channel= #supernova_flux:  {0}:    {1}:     {1}:    #{2}:    #{2}_smear\n'.format(c["parity"],c["flavor"],c["name"]))
            # Get the post-smearing efficiencies by channel
            try: 
                with open(sng/f'effic/effic_{c["name"]}_{detector_name}.dat') as f:
                    for l in f:
                        output.write('       @post_smearing_efficiencies = '+l)
            except IOError as e:
                logger.warning("Failed to open efficiency: "+traceback.format_exc())
            finally:
                output.write("\n>\n\n")
        #write postamble
        output.write(cat(sng/'glb/postamble.glb'))

    def apply_weights(detector_name, flux_fname, channels, do_smear):
        for c in channels:
            input_fname = sng/'out/{fl}_{ch}_{det}_events{sm}_unweighted.dat'.format(
                                        fl=flux_fname.stem,ch=c["name"], det=detector_name,
                                        sm="_smeared" if do_smear else "")
            output_fname = str(input_fname).replace('unweighted','weighted')
            data = load_datafile(input_fname)
            if len(data):
                data['events']*=c["weight"]
                np.savetxt(output_fname, data, fmt='%11g')

    def run_simulation(flux_fname, detector_name, do_apply_weight=False):
        """this function runs supernova"""
        material = det_materials[detector_name]
        channel_fname = sng/f'channels/channels_{material}.dat'
        channels = load_channels(channel_fname)
        globes_fname = sng/'supernova.glb'
        with open(globes_fname,'w') as output:
            prepare_globes_file(output,flux_fname,channels, detector_name)
        # Runs supernova
        os.chdir(sng)
        subprocess.run(['bin/supernova', flux_fname.stem,channel_fname, detector_name])
        # Calls apply_weights to format the output of supernova
        if do_apply_weight:
            logger.debug("Applying channel weighting factors to output")
            for do_smear in (True,False):
                apply_weights(detector_name,flux_fname,channels,do_smear)

    #THE DRIVER
    #This is the place where you can comment out any detectors you don't want to run
    #This runs the entire module, for each detector configuration
  
    if detector_input == "all":
        detector_input=list(det_materials.keys())
    if isinstance(detector_input,str):
        detector_input=[detector_input]
    logger.info(f'detector_input={detector_input}')
    for flux_fname in flux_files:
        for det in detector_input:
            try:
                run_simulation(flux_dir/flux_fname, det, "weight")
            except Exception as e:
                logger.exception(f'Failed processing "{det}" with "{flux_fname}": '+traceback.format_exc())

def collate(SNOwGLoBESdir, tarball_path, detector_input="all", skip_plots=False, verbose=False, remove_generated_files=True):
    """Collates SNOwGLoBES output files and generates plots or returns a data table.

    Parameters
    ----------
    SNOwGLoBESdir : str
        Path to directory where SNOwGLoBES is installed.
    tarball_path : str
        Path of compressed .tar file produced e.g. by ``generate_time_series()`` or ``generate_fluence()``.
    detector_input : str
        Name of detector. If ``"all"``, will use all detectors supported by SNOwGLoBES.
    skip_plots: bool
        If False, it gives as output the plot of the energy distribution for each time bin and for each interaction channel.
    verbose : bool
        Whether to generate verbose output, e.g. for debugging.
    remove_generated_files: bool
        Remove the output files from SNOwGLoBES, collated files, and .png's made for this snewpy run. 

    Returns
    -------
    dict
        Dictionary of data tables: One table per time bin; each table contains in the first column the energy bins, in the remaining columns the number of events for each interaction channel in the detector.
    """
    sng = Path(SNOwGLoBESdir)
    homebase = sng/'out'
    #Determine the files list
    with tarfile.open(tarball_path,'r') as f:
        flux_files = f.getnames() 
    flux_files = [f for f in flux_files if f.endswith('.dat')]


    FilesToCleanup = []

    def add_funct(flux, detector, smear, *arg):
        """Add_funct sums up relevant files from output generated above"""
        #Defining different variables for the combinatoric iterations
        #And reformatting some variables for different naming uses
        if smear.endswith('_w'):
            weight_value = "weighted"
        else:
            weight_value = "unweighted"

        if smear == "_smeared_w" or smear == "_smeared_u":
            smear_value = "smeared"
            smear_title = "detected"
            x_label = "Detected Energy (GeV)"
            y_label = "Events"
        else:
            smear_value = "unsmeared"
            smear_title = "interaction"
            x_label = "Neutrino Energy (GeV)"
            y_label = "Interaction Events"

        #a map to redefine detector labels
        det_labels_map={
            "ar40kt_eve":"ar40kt",
            "wc100kt30prct_eve":"wc100kt30prct"
            }
        det_label = det_labels_map.get(detector,detector)
            
        colors = ["k", "r", "g", "y", "b", "m", "c"]  # colors of graphed values
        compile_dict = {}

        events_all = []
        for input_val in arg:
            pattern = "{0}_{1}*{2}*{3}*.dat".format(flux,input_val, detector, smear)
            #Loop over files, corresponding to given pattern
            matching_files = list(homebase.glob(pattern))
            if matching_files:
                #Read all the files into a 2D array (Npoins * Nfiles)
                data=np.stack([load_datafile(f) for f in matching_files])
                #check that all energies are equal along the axis
                Es = data['energy']
                assert np.all(Es==Es[0])
                #sum up values and store them in list
                events_all.append(data['events'].sum(axis=0))
        if not events_all:
            return
        events_all = [Es[0].T]+events_all
        events_all = np.stack(events_all, axis=-1)
        #Creates the condensed data file & applies formatting
        # this part making new files with only useful info
        condensed_file = sng/'out/Collated_{0}_{1}_events_{2}_{3}.dat'.format(
                             flux, det_label, smear_value, weight_value)
        FilesToCleanup.append(condensed_file)
        header = 'Energy(GeV) '+' '.join(f'{a:20s}' for a in arg)+'\n'+'-'*100
        np.savetxt(condensed_file,events_all, header=header, fmt='%23.15g')

        if (skip_plots is False):
            r = 0

            if detector_input == "all":  # aka all detectors are being run, not a specific one
                relevant_list = arg
            else:
                relevant_list = sum_categories

            while r < len(relevant_list):  # Graph labels going back up here isn't working
                new_dict = {}
                for key in compile_dict:
                    new_dict[key] = compile_dict[key][r]

                if arg[r] == "ibd":  # Color labels on graph determined by name, given from value passed through flux details
                    name = "Inverse Beta Decay"
                elif arg[r] == "*_e_":
                    name = r'${\nu}_x+e^-$'
                elif arg[r] == "nc_":
                    name = "Neutral Current"
                elif arg[r] == "nue_Pb208_1n":
                    name = r'${\nu}_e {}^{208}Pb$' + " 1n"
                elif arg[r] == "nue_Pb208_2n":
                    name = r'${\nu}_e {}^{208}Pb$' + " 2n"
                elif arg[r] == "nue_O16":
                    name = r'${\nu}_e$' + " " + r'${}^{16}O$'
                elif arg[r] == "nuebar_O16":
                    name = r'$\bar{{\nu}_e}$' + " " + r'$^{16}O$'
                elif arg[r] == "nue_C12":
                    name = r'${\nu}_e$' + " " + r'${}^{12}C$'
                elif arg[r] == "nue_C13":
                    name = r'${\nu}_e$' + " " + r'${}^{13}C$'
                elif arg[r] == "nuebar_C12":
                    name = r'$\bar{{\nu}_e}$' + " " + r'$^{12}C$'
                elif arg[r] == "nue_Ar40":
                    name = r'${\nu}_e$' + " " + r'${}^{40}Ar$'
                elif arg[r] == "nuebar_Ar40":
                    name = r'$\bar{{\nu}_e}$' + " " + r'${}^{40}Ar$'
                else:
                    print("Invalid Input")

                if sum(new_dict.values()) != 0:  # plotting each color of bars individually

                    #Creates the plots
                    all_values = list(new_dict.values())
                    max_val = max(all_values)  # ensures only values greater than 0.1 show up on the plot and legend
                    if max_val > 0.1:
                        plt.plot(list(new_dict.keys()), list(new_dict.values()),
                                 linewidth=1, drawstyle='steps', color=colors[r], label=name)
                r += 1

            axes = plt.gca()
            axes.set_xlim([None, 0.10])
            axes.set_ylim([0.10, None])
            axes.set_yscale('log')
            plt.legend(bbox_to_anchor=(0.5, 0.5, 0.5, 0.5), loc='best', borderaxespad=0)  # formats complete graph
            plt.ylabel(y_label)
            plt.xlabel(x_label)
            plt.title(str(flux).capitalize() + " " + str(detector_label).capitalize() + " " +
                      str(weight_value).capitalize() + " " + str(smear_title).capitalize() + " Events")
            plt.savefig("{0}/out/{1}_{2}_{3}_{4}_log_plot.png".format(SNOwGLoBESdir, flux,
                        detector, smear_value, weight_value), dpi=300, bbox_inches='tight')
            FilesToCleanup.append("{0}/out/{1}_{2}_{3}_{4}_log_plot.png".format(SNOwGLoBESdir,
                                  flux, detector, smear_value, weight_value))
            plt.clf()

    #flux details                 # detector, smeared/unsmeared, *arg reaction placeholders
    smearvals = ["_smeared_w", "_smeared_u", "nts_w", "nts_u"]
    #weightvals = ["_weighted", "unweighted"]

    #The driver, runs through the addition function for every detector
    #The values at the end correspond to the different categories that each file is being summed into
    if detector_input == "all":
        detector_input=list(det_materials.keys())
    if isinstance(detector_input,str):
        detector_input=[detector_input]
    logger.info(f'detector_input={detector_input}')
    for flux_file in flux_files:
        for det in detector_input:
            for smearval in smearvals:
                add_funct(Path(flux_file).stem, det, smearval,*categories_map[det_materials[det]])

    #Now create tarball output
    #Makes a tarfile with the condensed data files and plots
    output_fname = tarball_path.rsplit('.tar',1)[0]+'_SNOprocessed'
    outname = Path(output_fname).name #name in archive
    with tarfile.open(output_fname+'.tar.gz', "w:gz") as tar:
        for file in homebase.glob('*Collated*'):
            tar.add(file, arcname=outname+'/'+str(file))
        for file in homebase.glob('*.png'):
            tar.add(file, arcname=outname+'/'+str(file))

    #collect the returned tables
    returned_tables = {}
    for file in homebase.glob('*Collated*'):
            data = np.loadtxt(file, skiprows=2, unpack=True)
            with open(file, 'r') as f:
                header= f.readline()
                header=header.lstrip('#')
            returned_tables[file.name] = {'data':data,'header':header}

    #removes all snowglobes output files, collated files, and .png's made for this snewpy run
    if remove_generated_files:
        for file in FilesToCleanup:
            os.remove(file)

    #Removes all the fluxfiles unzipped from the tarfile
    for file in flux_files:
        (sng/'fluxes'/file).unlink()
    try:
        os.remove(SNOwGLoBESdir + "/fluxes/parameterinfo")
    except OSError:
        pass

    return returned_tables
