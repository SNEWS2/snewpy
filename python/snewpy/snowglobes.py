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

import logging
import re
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from tqdm.auto import tqdm

import snewpy.models
from snewpy.flavor_transformation import *
from snewpy.neutrino import Flavor, MassHierarchy
from contextlib import contextmanager
from snewpy.flux import Flux
from snewpy.snowglobes_interface import SNOwGLoBES

logger = logging.getLogger(__name__)

def _load_model(model_path, model_type):
    model_class = getattr(snewpy.models.ccsn, model_type)
    return model_class(model_path)

def _load_transformation(transformation_type: str):
    # Choose flavor transformation. Use dict to associate the transformation name with its class.
    flavor_transformation_dict = {'NoTransformation': NoTransformation(), 'AdiabaticMSW_NMO': AdiabaticMSW(mh=MassHierarchy.NORMAL), 'AdiabaticMSW_IMO': AdiabaticMSW(mh=MassHierarchy.INVERTED), 'NonAdiabaticMSWH_NMO': NonAdiabaticMSWH(mh=MassHierarchy.NORMAL), 'NonAdiabaticMSWH_IMO': NonAdiabaticMSWH(mh=MassHierarchy.INVERTED), 'TwoFlavorDecoherence': TwoFlavorDecoherence(), 'ThreeFlavorDecoherence': ThreeFlavorDecoherence(), 'NeutrinoDecay_NMO': NeutrinoDecay(mh=MassHierarchy.NORMAL), 'NeutrinoDecay_IMO': NeutrinoDecay(mh=MassHierarchy.INVERTED)}
    return flavor_transformation_dict[transformation_type]


@contextmanager
def _archive_dir(tar_filename):
    """ Context manager to create a temporary dir, return it to the user.
    An exit of the context: pack all files in the directory to the tar archive and delete the dir
    """
    with TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        #return the path to user
        yield tempdir
        #pack and close
        with tarfile.open(tar_filename, 'w:bz2') as tar:
            for f in tempdir.iterdir():
                tar.add(f, arcname=f.name)


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

    model_path = Path(model_path)
    snmodel = _load_model(model_path, model_type)
    flavor_transformation = _load_transformation(transformation_type)

    # Subsample the model time. Default to 30 time slices.
    tmin = snmodel.get_time()[0].to_value('s')
    tmax = snmodel.get_time()[-1].to_value('s')
    if deltat is not None:
        tedges = np.arange(tmin, tmax, deltat) << u.s
    else:
        tedges = np.linspace(tmin, tmax, ntbins+2) << u.s
    times = 0.5*(tedges[1:] + tedges[:-1])
    dt = np.diff(tedges)

    energy = np.linspace(0, 100, 501) * u.MeV
    flux = snmodel.get_transformed_flux(times, energy, flavor_transformation, d*u.kpc)
    #multiply flux by the bin sizes
    flux.array *= 0.2*dt[0].to_value('s')
    # Save all to tar file
    if output_filename is None:
        output_filename = f'{model_path.stem}.{transformation_type}.{tmin:.3f} s,{tmax:.3f} s,{ntbins:d}-{d:.1f}'
    output_filename = model_path.parent/f'{output_filename}kpc.tar.bz2'

    with _archive_dir(output_filename) as tempdir:
        #creates file in tar archive that gives information on parameters
        with open(tempdir/'parameterinfo', 'w') as f:
            f.write(transformation_type)
        filename = f'{model_path.stem}.tbin{{n_time:01d}}.{transformation_type}'+\
                   f'.{tmin:.3f},{tmax:.3f},{ntbins:01d}-{d:.1f}kpc.dat'
        header=f'TBinMid={{time}}sec TBinWidth={dt[0]:g}s '+\
                ' EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2\n'
        flux.to_snowglobes(tempdir/filename,header)
    return output_filename


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

    model_path = Path(model_path)
    snmodel = _load_model(model_path,model_type)
    flavor_transformation = _load_transformation(transformation_type)

    #set the timings up
    #default if inputs are None: full time window of the model
    if tstart is None:
        tstart = snmodel.get_time()[0]
    if tend is None:
        tend = snmodel.get_time()[-1]
    if np.isscalar(tstart.value): tstart = [tstart]
    if np.isscalar(tend.value):   tend = [tend]
    nbin = len(tstart)
    tmin = min(tstart)
    tmax = max(tend)

    energy = np.linspace(0, 100, 501) * u.MeV
    times = snmodel.get_time()
    flux = snmodel.get_transformed_flux(times,energy,flavor_transformation, d*u.kpc)

    #multiply flux by the bin sizes
    flux.array*=0.2
    fluence = flux.integral('time',list(zip(tstart,tend)))
    if nbin==1: fluence=[fluence]
    # Generate output
    if output_filename is None:
        output_filename =f'{model_path.stem}.{transformation_type}.{tmin:.3f},{tmax:.3f},{nbin:d}-{d:.1f}kpc'
    output_filename = model_path.parent/f'{output_filename}.tar.bz2'

    with _archive_dir(output_filename) as tempdir:
        #creates file in tar archive that gives information on parameters
        with open(tempdir/'parameterinfo','w') as f:
            f.write(transformation_type)

        for i,(ta,tb) in enumerate(zip(tstart,tend)):
            t = 0.5*(tb+ta)
            dt = tb-ta
            filename = f'{model_path.stem}.tbin{i:01d}.{transformation_type}'+\
                       f'.{tmin:.3f},{tmax:.3f},{nbin:01d}-{d:.1f}kpc.dat'
            header=f'TBinMid={t:g}sec TBinWidth={dt:g}s '+\
                    'EBinWidth=0.2MeV Fluence at Earth for this timebin in neutrinos per cm^2\n'
            fluence[i].to_snowglobes(tempdir/filename,header)
 
    return output_filename


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
    
    sng = SNOwGLoBES(SNOwGLoBESdir)
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

re_chan_label = re.compile('nu(e|mu|tau)(bar|)_([A-Z][a-z]*)(\d*)_?(.*)')
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
    #save collated files:
    with TemporaryDirectory(prefix='snowglobes') as tempdir:
        tempdir = Path(tempdir)
        for det in tables:
            results[det] = {}
            for flux,t in tables[det].items():
                t = aggregate_channels(t,nc='nc_',e='_e')
                for w in ['weighted','unweighted']:
                    for s in ['smeared','unsmeared']:
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

       
