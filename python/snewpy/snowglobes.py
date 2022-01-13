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

import io
import logging
import os
import re
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from tqdm.auto import tqdm
from warnings import warn

import snewpy.models
from snewpy.flavor_transformation import *
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.snowglobes_interface import SNOwGLoBES, SimpleRate

logger = logging.getLogger(__name__)

def generate_time_series(model_path, model_type, transformation_type, d, output_filename=None, ntbins=30, deltat=None, snmodel_dict={}):
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
    snmodel_dict : dict
        Additional arguments for setting up the supernova model. See documentation of relevant ``SupernovaModel`` subclass for available options. (Optional)

    Returns
    -------
    str
        Path of compressed .tar file with neutrino flux data.
    """
    model_class = getattr(snewpy.models.ccsn, model_type)

    # Choose flavor transformation. Use dict to associate the transformation name with its class.
    flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    flavor_transformation = flavor_transformation_dict[transformation_type]

    model_dir, model_file = os.path.split(os.path.abspath(model_path))
    snmodel = model_class(model_path, **snmodel_dict)

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


def generate_fluence(model_path, model_type, transformation_type, d, output_filename=None, tstart=None, tend=None, snmodel_dict={}):
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
    snmodel_dict : dict
        Additional arguments for setting up the supernova model. See documentation of relevant ``SupernovaModel`` subclass for available options. (Optional)

    Returns
    -------
    str
        Path of compressed .tar file with neutrino flux data.
    """
    model_class = getattr(snewpy.models.ccsn, model_type)

    # Choose flavor transformation. Use dict to associate the transformation name with its class.
    flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    flavor_transformation = flavor_transformation_dict[transformation_type]

    model_dir, model_file = os.path.split(os.path.abspath(model_path))
    snmodel = model_class(model_path, **snmodel_dict)

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

def simulate(SNOwGLoBESdir, tarball_path, detector_input="all", verbose=False, detector_effects=True):
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
        [DEPRECATED, DO NOT USE.]
    detector_effects : bool
         Whether to account for detector smearing and efficiency.
    """
    if verbose:
        warn(f"The 'verbose' parameter to 'snewpy.snowglobes.simulate()' is deprecated and should not be used.", FutureWarning)
    
    sng = SNOwGLoBES(SNOwGLoBESdir) if detector_effects else SimpleRate(SNOwGLoBESdir)

    if detector_input == 'all':
        detector_input = list(sng.detectors)
        detector_input.remove('d2O')
    elif isinstance(detector_input,str):
        detector_input = [detector_input]
    
    result = {}
    #Extracts data from tarfile and sets up lists of paths and fluxfilenames for later use
    with TemporaryDirectory(prefix='snowglobes') as tempdir:
        with tarfile.open(tarball_path) as tar:
            tar.extractall(tempdir)

        flux_files = list(Path(tempdir).glob('*.dat'))
        if len(detector_input)>0:
            detector_input = tqdm(detector_input, desc='Detectors', leave=False)
        for det in detector_input:
            res=sng.run(flux_files, det)
            result[det]=dict(zip((f.stem for f in flux_files),res))

    # save result to file for re-use in collate()
    cache_file = tarball_path[:tarball_path.rfind('.tar')] + '.npy'
    logging.info(f'Saving simulation results to {cache_file}')
    np.save(cache_file, result)
    return result 

re_chan_label = re.compile(r'nu(e|mu|tau)(bar|)_([A-Z][a-z]*)(\d*)_?(.*)')
def get_channel_label(c):
    mapp = {'nc':'NeutralCurrent',
            'ibd':'Inverse Beta Decay',
            'e':r'${\nu}_x+e^-$'}
    def gen_label(m):
        flv,bar,Nuc,num,res = m.groups()
        if flv!='e':
            flv='\\'+flv
        if bar:
            bar='\\'+bar
        s = f'${bar}{{\\nu}}_{flv}$ '+f'${{}}^{{{num}}}{Nuc}$ '+res
        return s

    if c in mapp:
        return mapp[c]
    else: 
        return re_chan_label.sub(gen_label, c) 

def collate(SNOwGLoBESdir, tarball_path, detector_input="all", skip_plots=False, verbose=False, remove_generated_files=True, smearing = True):
    """Collates SNOwGLoBES output files and generates plots or returns a data table.

    Parameters
    ----------
    SNOwGLoBESdir : str
        [DEPRECATED, DO NOT USE.]
    tarball_path : str
        Path of compressed .tar file produced e.g. by ``generate_time_series()`` or ``generate_fluence()``.
    detector_input : str
        [DEPRECATED, DO NOT USE. SNEWPY will use all detectors included in the tarball.]
    skip_plots: bool
        If False, it gives as output the plot of the energy distribution for each time bin and for each interaction channel.
    verbose : bool
        [DEPRECATED, DO NOT USE.]
    remove_generated_files: bool
        [DEPRECATED, DO NOT USE.]
    smearing: bool
        Also consider results with smearing effects.

    Returns
    -------
    dict
        Dictionary of data tables: One table per time bin; each table contains in the first column the energy bins, in the remaining columns the number of events for each interaction channel in the detector.
    """
    if verbose:
        warn(f"The 'verbose' parameter to 'snewpy.snowglobes.collate()' is deprecated and should not be used.", FutureWarning)
    if detector_input != "all":
        warn(f"The 'detector_input' parameter to 'snewpy.snowglobes.simulate()' is deprecated and should not be used.", FutureWarning)
    if not remove_generated_files:
        warn(f"The 'remove_generated_files' parameter to 'snewpy.snowglobes.simulate()' is deprecated and should not be used.", FutureWarning)

    def aggregate_channels(table, **patterns):
        #rearrange the table to have only channel column
        levels = list(table.columns.names)
        levels.remove('channel')
        t = table.stack(levels)
        for name,pattern in patterns.items():
            #get channels which contain `like`
            t_sel = t.filter(like=pattern)
            #sum over them and save to a separate column
            t_agg = t_sel.sum(axis='columns')
            #drop processed channels
            t.drop(t_sel.columns, axis='columns',inplace=True)
            t[name]=t_agg #fill the column
        #return table with the original levels order
        t = t.unstack(levels)
        t = t.reorder_levels(table.columns.names, axis=1)
        return t
        
    def do_plot(table, params):
        #plotting the events from given table
        flux,det,weighted,smeared = params
        for c in table.columns:
            if table[c].max() > 0.1:
                plt.plot(table[c],drawstyle='steps',label=get_channel_label(c), lw=1)
        plt.xlim(right=0.10)
        plt.ylim(bottom=0.10)
        plt.yscale('log')
        plt.legend(bbox_to_anchor=(0.5, 0.5, 0.5, 0.5), loc='best', borderaxespad=0)  # formats complete graph
        smear_title = 'Interaction' if smeared=='unsmeared' else 'Detected'
        plt.title(f'{flux} {det.capitalize()} {weighted.capitalize()} {smear_title} Events')
        if smeared=='smeared':
            plt.xlabel('Detected Energy (GeV)')
            plt.ylabel('Events')  
        else:
            plt.xlabel('Neutrino Energy (GeV)')
            plt.ylabel('Interaction Events')  

    #read the results from storage
    cache_file = tarball_path[:tarball_path.rfind('.tar')] + '.npy'
    logging.info(f'Reading tables from {cache_file}')
    tables = np.load(cache_file, allow_pickle=True).tolist()
    #This output is similar to what produced by:
    #tables = simulate(SNOwGLoBESdir, tarball_path,detector_input)

    #dict for old-style results, for backward compatibiity
    results = {}
    smearing_options = ['smeared','unsmeared'] if smearing else ['unsmeared']
    #save collated files:
    with TemporaryDirectory(prefix='snowglobes') as tempdir:
        tempdir = Path(tempdir)
        for det in tables:
            results[det] = {}
            for flux,t in tables[det].items():
                t = aggregate_channels(t,nc='nc_',e='_e')
                for w in ['weighted','unweighted']:
                    for s in smearing_options:
                        table = t[w][s]
                        filename_base = f'{flux}_{det}_events_{s}_{w}'
                        filename = tempdir/f'Collated_{filename_base}.dat'
                        #save results to text files
                        with open(filename,'w') as f:
                            f.write(table.to_string(float_format='%23.15g'))
                        #format the results for the output
                        header = 'Energy '+' '.join(list(table.columns))
                        data = table.to_numpy().T
                        index = table.index.to_numpy()
                        data = np.concatenate([[index],data])
                        results[filename.name] = {'header':header,'data':data}
                        #optionally plot the results
                        if skip_plots is False:
                            plt.figure(dpi=300)
                            do_plot(table,(flux,det,w,s))
                            filename = tempdir/f'{filename_base}_log_plot.png'
                            plt.savefig(filename.with_suffix('.png'), dpi=300, bbox_inches='tight')
        #Make a tarfile with the condensed data files and plots
        output_name = Path(tarball_path).stem
        output_name = output_name[:output_name.rfind('.tar')]+'_SNOprocessed'
        output_path = Path(tarball_path).parent/(output_name+'.tar.gz')
        with tarfile.open(output_path, "w:gz") as tar:
            for file in tempdir.iterdir():
                tar.add(file,arcname=output_name+'/'+file.name)
        logging.info(f'Created archive: {output_path}')
    return results 
