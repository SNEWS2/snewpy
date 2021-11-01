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
from subprocess import call

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

import snewpy.models
from snewpy.flavor_transformation import *
from snewpy.neutrino import Flavor, MassHierarchy

mpl.use('Agg')


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

    #Extracts data from tarfile and sets up lists of paths and fluxfilenames for later use

    if tarfile.is_tarfile(tarball_path):  # extracts tarfile
        if (verbose):
            print("Valid Tar file")
        tar = tarfile.open(tarball_path)
        TarballFileNames = tar.getnames()
        tar.extractall(path=SNOwGLoBESdir + "/fluxes")
        tar.close()
        if (verbose):
            print("Tar file data extracted")
    elif zipfile.is_zipfile(tarball_path):  # extracts zipfile
        if (verbose):
            print("Valid Zip file")
        zip = zipfile.ZipFile(tarball_path)
        TarballFileNames = zip.namelist()
        for fileName in TarballFileNames:
            if str(fileName)[0] == "." or str(fileName)[0] == "_":
                continue
            else:
                zip.extract(fileName, path=SNOwGLoBESdir + "/fluxes")
                if (verbose):
                    print("zipfile extracted")
        zip.close()
        if (verbose):
            print("Zip file data extracted")
    else:
        print("Invalid Tarfile")  # doesn't accept anything other than tar and zip files

    #Creates a new out folder if one does not exist

    out_folder = "{0}/out".format(SNOwGLoBESdir)
    try:
        os.mkdir(out_folder)
    except OSError as e:
        if (verbose):
            print("out folder already exists")
    else:
        if (verbose):
            print("out folder made")

    #Creates directories for each flux file within the out folder

    FluxFileNameStems = []
    #extension = '.dat'
    for IndividualFile in TarballFileNames:
        extension_beginning = str(IndividualFile).rfind(".")
        if extension_beginning > -1:
            extension = str(IndividualFile[extension_beginning:])
            ExtensionFound = str(IndividualFile).rfind(extension)
            # gets rid of extension at the end of flux files
            FluxFileNameStems.append(str(IndividualFile)[0:ExtensionFound])
        else:
            continue

    ##################################################################################################
    #The main function that reformats data from a bunch of different files in order to run supernova.c

    # this function runs supernova
    def format_globes_for_supernova(FluxNameStem, ChannelName, ExpConfigName, NoWeight=1, extension=extension):
        ExeName = "bin/supernova"  # Supernova
        ChannelFileName = "{0}/channels/channels_{1}.dat".format(SNOwGLoBESdir, ChannelName)

        GlobesFileName = "{0}/supernova.glb".format(SNOwGLoBESdir)
        GLOBESFILE = open(GlobesFileName, 'w')
        GLOBESFILE.truncate(0)

        # Write PREAMBLE to GLOBESFILE

        PREAMBLE = open("{0}/glb/preamble.glb".format(SNOwGLoBESdir))

        for line in PREAMBLE:

            GLOBESFILE.write(line)

        PREAMBLE.close()

        # Checks FluxFileNamePath existence

        FluxFileNamePath = "fluxes/" + FluxNameStem + extension
        if not os.path.isfile(SNOwGLoBESdir+"/"+FluxFileNamePath):
            print(("Flux file name {0} not found".format(FluxFileNamePath)))
            return "Incomplete"

        # Write modified FLUX to GLOBESFILE; replace the flux file name with the input argument

        FLUX = open("{0}/glb/flux.glb".format(SNOwGLoBESdir))

        for line in FLUX:
            if "flux_file" in line:
                line = "        @flux_file=  \"{0}\"\n".format(FluxFileNamePath)
            GLOBESFILE.write(line)

        FLUX.close()

        # Now go through the channels and put in the relevant lines

        # First, smearing for each channel

        # Checks channel existence
        if not os.path.isfile(ChannelFileName):
            print(("Channel file name {0} not found".format(ChannelFileName)))
            return "Incomplete"

        # Grabs channel name from each line in CHANFILE and write it (modified) to GLOBESFILE
        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')
            ChanName = TempArray[0]

            GLOBESFILE.write("include \"{0}/smear/smear_{1}_{2}.dat\"\n".format(SNOwGLoBESdir, ChanName, ExpConfigName))

        CHANFILE.close()

        # Checks DetectorFileNamePath existence

        DetectorFileName = "detector_configurations.dat"
        DetectorFileNamePath = "{0}/{1}".format(SNOwGLoBESdir, DetectorFileName)

        if not os.path.isfile(DetectorFileNamePath):
            print(("Detector file name {0} not found".format(DetectorFileName)))
            return "Incomplete"

        # Defines Masses and NormFact (hashes, A.K.A. dictionaries) for later use

        DETFILENAME = open(DetectorFileNamePath)
        Masses = {}
        NormFact = {}
        for line in DETFILENAME:

            if line[-1] == "\n":
                line = line[:-1]

            if "#" in line:
                continue

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')
            DetName = TempArray[0]
            Masses[DetName] = TempArray[1]
            NormFact[DetName] = TempArray[2]

            if (DetName == "" or Masses == "" or NormFact == ""):
                continue
        DETFILENAME.close()

        # Mass and target normalization by species

        # Checks Masses[ExpConfigName] existence

        if not Masses.get(ExpConfigName):
            print("Error: please enter a valid experiment configuration")
            return "Incomplete"

        TargetMass = format(float(Masses[ExpConfigName]) * float(NormFact[ExpConfigName]), "13.6f")

        # Writes modified DETECTOR to GLOBESFILE; replace the mass with the calculated TargetMass

        DETECTOR = open("{0}/glb/detector.glb".format(SNOwGLoBESdir))

        for line in DETECTOR:

            if "mass" in line:
                line = "$target_mass=  {0}".format(TargetMass)

            GLOBESFILE.write(line)

        DETECTOR.close()
        # Now the cross-sections.  Note that some of these are repeated even though
        # it is not necessary (xscns for several flavors can be in the same file).
        # This is just to make a consistent loop over channels.

        GLOBESFILE.write("\n\n /******** Cross-sections *********/\n \n")

        # Writes a modified CHANFILE to GLOBESFILE

        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')

            ChanName = TempArray[0]

            GLOBESFILE.write("cross(#{0})<\n".format(ChanName))
            GLOBESFILE.write("      @cross_file= \"{0}/xscns/xs_{1}.dat\"\n".format(SNOwGLoBESdir, ChanName))
            GLOBESFILE.write(">\n")

        CHANFILE.close()
        # Now, the channel definitions

        GLOBESFILE.write("\n /******** Channels *********/\n \n")

        # Checks ChannelFileName existence

        if not os.path.isfile(ChannelFileName):
            print(("Channel file name {0} not found".format(ChannelFileName)))
            return "Incomplete"

        # Writes a modified CHANFILE and EFF_FILE to GLOBESFILE
        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')

            ChanName = TempArray[0]
            CpState = TempArray[2]
            InFlav = TempArray[3]

            GLOBESFILE.write("channel(#{0}_signal)<\n".format(ChanName))
            GLOBESFILE.write("      @channel= #supernova_flux:  {0}:    {1}:     {2}:    #{3}:    #{4}_smear\n".format(
                CpState, InFlav, InFlav,                                                                         ChanName, ChanName))

            # Get the post-smearing efficiencies by channel

            EffFile = "{0}/effic/effic_{1}_{2}.dat".format(SNOwGLoBESdir, ChanName, ExpConfigName)
            try:
                EFF_FILE = open(EffFile)
            except IOError:
                GLOBESFILE.write("\n>\n\n")
                continue
            else:
                for line in EFF_FILE:
                    GLOBESFILE.write("       @post_smearing_efficiencies = {0}".format(line))

                GLOBESFILE.write("\n>\n\n")

                EFF_FILE.close()
        CHANFILE.close()

        # End matter

        # Writes POSTAMBLE to GLOBESFILE

        POSTAMBLE = open("{0}/glb/postamble.glb".format(SNOwGLoBESdir))

        for line in POSTAMBLE:

            GLOBESFILE.write(line)

        POSTAMBLE.close()
        GLOBESFILE.close()

        # Runs supernova
        os.chdir(SNOwGLoBESdir)
        ComString = "{0} {1} {2} {3}".format(ExeName, FluxNameStem, ChannelFileName, ExpConfigName)
        call(ComString, shell=True)

        # Calls apply_weights to format the output of supernova

        if NoWeight != 1:
            if (verbose):
                print("Applying channel weighting factors to output")

            apply_weights("unsmeared", FluxNameStem, ChannelFileName, ExpConfigName)
            apply_weights("smeared", FluxNameStem, ChannelFileName, ExpConfigName)

        return "Complete"

    #######################################################################
    #This function applies the weights to the weightedfilename via tempfile

    def apply_weights(smearing, FluxNameStem, ChannelFileName, ExpConfigName):
        # Opens CHANFILE in order for the code to be able to find all of the output files (from supernova)

        #formats for naming of the files
        if smearing == "unsmeared":
            smear_input = ""  # snowglobes1.2 does not insert "_unsmeared" into its output filenames
        elif smearing == "smeared":
            smear_input = "_smeared"
        else:
            print("Invalid smearing input")

        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')

            ChanName = TempArray[0]
            NumTargetFactordraft = TempArray[4]
            if NumTargetFactordraft[-1] == ".":
                NumTargetFactor = float(NumTargetFactordraft[:-1])  # the weighting factor
            else:
                NumTargetFactor = float(NumTargetFactordraft)

            # Write weighted file
            # This applies the weighting factor to unweighted data and reformats a few lines for consistency
            UnweightFileName = "{0}/out/{1}_{2}_{3}_events{4}_unweighted.dat".format(
                SNOwGLoBESdir, FluxNameStem, ChanName, ExpConfigName, smear_input)
            WeightedFileName = "{0}/out/{1}_{2}_{3}_events{4}_weighted.dat".format(
                SNOwGLoBESdir, FluxNameStem, ChanName, ExpConfigName, smear_input)

            with open(UnweightFileName, "r") as f:
                WEIGHT = open(WeightedFileName, "w")
                for line in f:

                    clean_line = re.sub('\s*', '', line)
                    clean_line = re.sub('\s+', ' ', line)

                    if "---" in line:
                        WEIGHT.write(line + "\n")
                    elif "Total:" in line:
                        continue
                    else:
                        TempArray2 = clean_line.split(' ')
                        Enbin = TempArray2[0]
                        Evrate = TempArray2[1]
                        if Enbin != "":  # writing the new weighted lines to a file
                            newline = "{0}            {1} ".format(Enbin, float(Evrate) * float(NumTargetFactor))
                            WEIGHT.write(newline + "\n")
                        elif Enbin == "" and Evrate == "0.05" or Evrate == "0.051":
                            newline = "{0}            {1} ".format(
                                Evrate, float(TempArray2[2]) * float(NumTargetFactor))
                            WEIGHT.write(newline + "\n")
                        else:
                            continue
                WEIGHT.close()
        CHANFILE.close()

    fluxes_list = []  # for a list of IndividualFluxFileName files

    #THE DRIVER
    #This is the place where you can comment out any detectors you don't want to run
    #This runs the entire module, for each detector configuration

    for IndividualFluxFile in FluxFileNameStems:
        if (verbose):
            print(IndividualFluxFile)
        position = FluxFileNameStems.index(IndividualFluxFile) + 1
        total = len(FluxFileNameStems)
        percent_calc = (position/total)*100
        percentage_done = str(round(percent_calc, 2))

        if detector_input == "all":
            if (verbose):
                print("Running all detectors")
            if format_globes_for_supernova(IndividualFluxFile, "water", "icecube", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||icecube".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "water", "wc100kt30prct", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||wc100kt30prct".format(IndividualFluxFile))
            #if format_globes_for_supernova(IndividualFluxFile, "water", "wc100kt30prct_he", "weight") == "Complete":
                # if verbose: print("Finished {0}||wc100kt30prct_he".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "water", "wc100kt15prct", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||wc100kt15prct".format(IndividualFluxFile))
            #if format_globes_for_supernova(IndividualFluxFile, "water", "hyperk30prct", "weight") == "Complete": #not working
                # if verbose: print("Finished {0}||hyperk30prct".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "argon", "ar40kt", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||ar40kt".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "lead", "halo1", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||halo1".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "lead", "halo2", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||halo2".format(IndividualFluxFile))
            #if format_globes_for_supernova(IndividualFluxFile, "argon", "ar40kt_he", "weight") == "Complete":
                # if verbose: print("Finished {0}||ar40kt_he".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "scint", "scint20kt", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||scint20kt".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "nova_soup", "novaND", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||novaND".format(IndividualFluxFile))
            if format_globes_for_supernova(IndividualFluxFile, "nova_soup", "novaFD", "weight") == "Complete":
                if verbose:
                    print("Finished {0}||novaFD".format(IndividualFluxFile))
            fluxes_list.append(IndividualFluxFile)
            if verbose:
                print("\n\n\nCalculations are " + percentage_done + "% completed.\n\n\n")
        else:  # This is called if you choose to input a single detector in SNEWPY.py, and just run that one
            if verbose:
                print("Running selected detector:", detector_input)
            if detector_input in ("icecube", "wc100kt30prct", "wc100kt30prct_he", "wc100kt15prct", "hyperk30prct", "km3net"):
                detector_material = "water"
            elif detector_input in ("ar40kt_he", "ar40kt"):
                detector_material = "argon"
            elif detector_input in ("novaND", "novaFD"):
                detector_material = "nova_soup"
            elif detector_input in ("scint20kt",):
                detector_material = "scint"
            elif detector_input in ("halo1", "halo2"):
                detector_material = "lead"
            else:
                print("Unanticipated value for detector input")

            detector_output = detector_input
            if verbose:
                print("Detector material:", detector_material)
            if format_globes_for_supernova(IndividualFluxFile, detector_material, detector_input, "weight") == "Complete":
                if verbose:
                    print("Finished {0}||{1}".format(IndividualFluxFile, detector_output))
            fluxes_list.append(IndividualFluxFile)
            if (verbose):
                print('\n'*3)
            if (int(percent_calc) % 10 == 0):
                print("Calculations are " + percentage_done + "% completed. (", position, " of ", total, ")")
            if (verbose):
                print('\n'*3)


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
    model_dir, tarball = os.path.split(os.path.abspath(tarball_path))

    #Determines type of input file

    if ".tar.bz2" in str(tarball):
        outputnamestem = tarball[0:str(tarball).rfind(".tar.bz2")]
        tar = tarfile.open(tarball_path)
        TarballFileNames = tar.getnames()
        tar.close()
    elif ".tar.gz" in str(tarball):
        outputnamestem = tarball[0:str(tarball).rfind(".tar.gz")]
        tar = tarfile.open(tarball_path)
        TarballFileNames = tar.getnames()
        tar.close()
    elif ".zip" in str(tarball):
        outputnamestem = tarball[0:str(tarball).rfind(".zip")]
        zip = zipfile.ZipFile(tarball_path)
        TarballFileNames = zip.namelist()
        zip.close()
    else:
        print("Invalid Tar file")

    #Determining the configuration of the files and directories within the tarfile and defining variables for uniformity

    FluxFileNameStems = []
    extension = '.dat'
    for IndividualFile in TarballFileNames:
        extension_beginning = str(IndividualFile).rfind(".")
        if extension_beginning > -1:
            extension = str(IndividualFile[extension_beginning:])
            ExtensionFound = str(IndividualFile).rfind(extension)
            # gets rid of extension at the end of flux files
            FluxFileNameStems.append(str(IndividualFile)[0:ExtensionFound])
        else:
            continue

    smearvals = ["_smeared_w", "_smeared_u", "nts_w", "nts_u"]
    #weightvals = ["_weighted", "unweighted"]

    FilesToCleanup = []

    #Add_funct sums up relevant files from output generated above
    def add_funct(flux, detector, smear, *arg):
        homebase = SNOwGLoBESdir + "/out/"

        #Defining different variables for the combinatoric iterations
        #And reformatting some variables for different naming uses
        if smear == "_smeared_w" or smear == "nts_w":
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

        if detector == "ar40kt_eve":
            detector_label = "ar40kt"
            detector_label2 = "ar40kt"
        elif detector == "ar40kt_he":
            detector_label = "ar40kt he"
            detector_label2 = "ar40kt_he"
        elif detector == "wc100kt30prct_eve":
            detector_label = "wc100kt30prct"
            detector_label2 = "wc100kt30prct"
        elif detector == "wc100kt30prct_he":
            detector_label = "wc100kt30prct he"
            detector_label2 = "wc100kt30prct_he"
        else:
            detector_label = detector
            detector_label2 = detector

        colors = ["k", "r", "g", "y", "b", "m", "c"]  # colors of graphed values
        compile_dict = {}

        #Splits values from each file into Energy and Events, then sums appropriate ones
        for input_val in arg:
            flux_filenames = os.listdir(homebase)
            final_dict = {}
            temp_dict = {}
            for afile in flux_filenames:  # Does the actual summing of essential values in each file
                name_iteration = "{0}_{1}*{2}*{3}*".format(flux, input_val, detector, smear)
                if fnmatch.fnmatch(afile, name_iteration):  # and "Collated" not in str(afile):
                    FilesToCleanup.append(homebase+"/"+afile)
                    with open("{0}/{1}".format(homebase, afile)) as fileinsides:
                        #Splits each line & converts Energy and Events into floats
                        for line in fileinsides.readlines():
                            content = line.split()
                            if not content:
                                continue
                            elif content[0].startswith('---'):
                                break
                            energy = float(content[0])
                            events = float(content[1])

                            temp_dict[energy] = events
                            if len(final_dict) < len(temp_dict):
                                final_dict[energy] = 0

                    for energy in temp_dict:
                        final_dict[energy] = final_dict[energy] + temp_dict[energy]  # sums appropriate values #here

                if len(compile_dict) < len(final_dict):
                    for energy in final_dict:
                        compile_dict[energy] = []

            for k, v in list(final_dict.items()):
                # This is the dictionary with energy bins and lists of events corresponding to interaction type
                compile_dict[k].append(v)

        #Creates the condensed data file & applies formatting
        # this part making new files with only useful info
        condensed_file = "{0}/out/Collated_{1}_{2}_events_{3}_{4}.dat".format(
            SNOwGLoBESdir, flux, detector_label2, smear_value, weight_value)
        FilesToCleanup.append(condensed_file)
        new_f = open(condensed_file, "w")
        if len(arg) == 4:
            new_f.write("Energy(GeV)            {0}                    {1}                    {2}                {3}".format(
                arg[0], arg[1], arg[2], arg[3]))
        elif len(arg) == 5:
            new_f.write("Energy(GeV)            {0}                    {1}                    {2}                {3}                    {4}".format(
                arg[0], arg[1], arg[2], arg[3], arg[4]))
        else:
            new_f.write("Energy(GeV)            {0}                    {1}                    {2}                {3}                    {4}                    {5}".format(
                arg[0], arg[1], arg[2], arg[3], arg[4], arg[5]))
        new_f.write("\n")
        new_f.write("-"*100)
        new_f.write("\n")

        #converts the dictionary into readable info
        for key in compile_dict:
            first_space_num = 23 - len(str(key))
            first_space = " " * first_space_num
            value = list(compile_dict[key])
            if len(value) != 1:
                i = 0
                while i < len(value):  # Formats the lines so that there are sufficient spaces between each value for them to be readable
                    if len(str(value[i])) > 16:
                        separation_space_num = 23-len(str(round(value[i], 16)))
                        separation_space = " " * separation_space_num
                        value[i] = str(round(value[i], 16)) + separation_space
                        i += 1
                    else:
                        separation_space_num = 23-len(str(round(value[i], 15)))
                        separation_space = " " * separation_space_num
                        value[i] = str(round(value[i], 15)) + separation_space
                        i += 1
                if len(arg) == 4:
                    new_f.write("{0}{1}{2}{3}{4}{5}".format(key, first_space, value[0], value[1], value[2], value[3]))
                elif len(arg) == 5:
                    # for normal, next is for spreadsheet counting
                    new_f.write("{0}{1}{2}{3}{4}{5}{6}".format(key, first_space,
                                value[0], value[1], value[2], value[3], value[4]))
                else:
                    # for normal, next is for spreadsheet counting
                    new_f.write("{0}{1}{2}{3}{4}{5}{6}{7}".format(key, first_space,
                                value[0], value[1], value[2], value[3], value[4], value[5]))
            new_f.write('\n')
        new_f.close()

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

    #The driver, runs through the addition function for every detector
    #The values at the end correspond to the different categories that each file is being summed into
    if detector_input == "all":
        for single_flux in FluxFileNameStems:
            for smearval in smearvals:
                add_funct(str(single_flux), "wc100kt15prct", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")  # everything after smearval corresponds to a *arg value
                add_funct(str(single_flux), "wc100kt30prct", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct(str(single_flux), "wc100kt30prct_he", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct(str(single_flux), "halo1", smearval, "nc_", "*_e_", "nue_Pb208_1n", "nue_Pb208_2n")
                add_funct(str(single_flux), "halo2", smearval, "nc_", "*_e_", "nue_Pb208_1n", "nue_Pb208_2n")
                add_funct(str(single_flux), "ar40kt_eve", smearval, "nc_", "*_e_", "nue_Ar40", "nuebar_Ar40")
                add_funct(str(single_flux), "ar40kt_he", smearval, "nc_", "*_e_", "nue_Ar40", "nuebar_Ar40")
                add_funct(str(single_flux), "icecube", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct(str(single_flux), "hyperk30prct", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct(str(single_flux), "scint20kt", smearval, "nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12", "nue_C13")
                add_funct(str(single_flux), "novaND", smearval, "nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12")
                add_funct(str(single_flux), "novaFD", smearval, "nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12")

            flux_position = FluxFileNameStems.index(single_flux) + 1
            total_fluxes = len(FluxFileNameStems)
            plot_percent_calc = (flux_position/total_fluxes)*100
            plot_percentage_done = str(round(plot_percent_calc, 2))
            if (verbose):
                print('\n'*3)
            if (verbose):
                print("Plots are " + plot_percentage_done + "% completed.")
            if (verbose):
                print('\n'*3)
    else:  # if just one detector is inputted in SNEWPY.py
        for single_flux in FluxFileNameStems:
            for smearval in smearvals:
                if (verbose):
                    print("Running inputted detector")
                if detector_input in ("wc100kt15prct", "wc100kt30prct", "wc100kt30prct_he", "icecube", "hyperk30prct", "km3net"):
                    sum_categories = ["nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16"]
                elif detector_input in ("halo1", "halo2"):
                    sum_categories = ["nc_", "*_e_", "nue_Pb208_1n", "nue_Pb208_2n"]
                elif detector_input in ("ar40kt", "ar40kt_he"):
                    sum_categories = ["nc_", "*_e_", "nue_Ar40", "nuebar_Ar40"]
                elif detector_input in ("scint20kt"):
                    sum_categories = ["nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12", "nue_C13"]
                elif detector_input in ("novaND", "novaFD"):
                    sum_categories = ["nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12"]
                else:
                    print("Unanticipated value for detector input")
                argList1 = [str(single_flux), detector_input, smearval]
                argList1.extend(sum_categories)
                add_funct(*argList1)

                flux_position = FluxFileNameStems.index(single_flux) + 1
                total_fluxes = len(FluxFileNameStems)
                plot_percent_calc = (flux_position/total_fluxes)*100
                plot_percentage_done = str(round(plot_percent_calc, 2))
            if (verbose):
                print('\n'*3)
            if (verbose):
                print("Plots are " + plot_percentage_done + "% completed.")
            if (verbose):
                print('\n'*3)

    #Now create tarball output
    #Makes a tarfile with the condensed data files and plots
    tar = tarfile.open(model_dir + "/" + outputnamestem + "_SNOprocessed.tar.gz", "w:gz")
    for file in os.listdir(SNOwGLoBESdir + "/out"):
        if "Collated" in str(file):
            tar.add(SNOwGLoBESdir + "/out/" + file, arcname=outputnamestem+'_SNOprocessed/'+file)
        elif ".png" in str(file):
            tar.add(SNOwGLoBESdir + "/out/" + file, arcname=outputnamestem+'_SNOprocessed/'+file)
        else:
            continue
    tar.close()

    returned_tables = {}
    for file in os.listdir(SNOwGLoBESdir + "/out"):
        if "Collated" in str(file):
            returned_tables[file] = {}
            with open(SNOwGLoBESdir + "/out/"+file, 'r') as fstream:
                returned_tables[file]['header'] = fstream.readline()
            returned_tables[file]['data'] = np.loadtxt(SNOwGLoBESdir + "/out/" + file, skiprows=2, unpack=True)

    #removes all snowglobes output files, collated files, and .png's made for this snewpy run
    if remove_generated_files:
        for file in FilesToCleanup:
            os.remove(file)

    #Removes all the fluxfiles unzipped from the tarfile
    for file in FluxFileNameStems:
        os.remove(SNOwGLoBESdir + "/fluxes/" + file + extension)
    try:
        os.remove(SNOwGLoBESdir + "/fluxes/parameterinfo")
    except OSError:
        pass

    return returned_tables
