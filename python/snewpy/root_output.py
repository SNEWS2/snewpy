# -*- coding: utf-8 -*-
""" The root_output module allows to collate all the SNOwGLoBES outputs in a single ROOT file with 2D histograms
    * **Collating SNOwGLoBES outputs:**
    This step puts together all the interaction channels and time bins evaluated by SNOwGLoBES in a single ROOT file (for each detector).
    The output 2D histograms allow to build the detected neutrino energy spectrum and neutrino time distribution, for each reaction channel or the sum of them.
"""

import io
import logging
import os
import re
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from tqdm.auto import tqdm

import snewpy.models
from snewpy.flavor_transformation import *
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.snowglobes_interface import SNOwGLoBES

logger = logging.getLogger(__name__)

import ROOT

def collate(tarball_path, detector_input="all", verbose=False):
    """Collates SNOwGLoBES output files and returns a ROOT file with 2D histograms.

    Parameters
    ----------
    tarball_path : str
        Path of compressed .tar file produced e.g. by ``generate_time_series()`` or ``generate_fluence()``.
    detector_input : str
        Name of detector. If ``"all"``, will use all detectors supported by SNOwGLoBES.
    verbose : bool
        Whether to generate verbose output, e.g. for debugging.

    Returns
    -------
    string
        Output file name
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
        

    #read the results from storage
    cache_file = tarball_path[:tarball_path.rfind('.tar')] + '.npy'
    logging.info(f'Reading tables from {cache_file}')
    tables = np.load(cache_file, allow_pickle=True).tolist()
    #This output is similar to what produced by:

    #dict for old-style results, for backward compatibiity
    results = {}
    outputnames = {}
    hist_dict = {}
    tmin,tmax,ntbins = 0,0,0
    #save collated files:
    with TemporaryDirectory(prefix='snowglobes') as tempdir:
        tempdir = Path(tempdir)
        for det in tables:
            results[det] = {}
            first = True
            for flux,t in tables[det].items():
                # Get time binning
                if first:
                    # Time boundaries are formatted to be %0.3f in generate_time_series
                    print(flux)
                    tmin,tmax,ntbins = flux.split(',')
                    tmin = float(".".join(tmin.split('.')[-2:]))
                    tmax = float(tmax)
                    ntbins = int(ntbins.split('-')[0])
                # Find general name for flux file as well as value of time bin
                flux_start, flux_end = flux.split('tbin')
                flux_split = flux_end.split('.')
                time_bin = int(flux_split[0])
                time_val = (time_bin - 0.5) * (tmax - tmin)/ntbins + tmin
                flux_end = ".".join(flux_split[1:])
                flux_root = flux_start + flux_end
                # Aggregate channels
                t = aggregate_channels(t,nc='nc_',e='_e')
                for w in ['weighted','unweighted']:
                    for s in ['smeared','unsmeared']:
                        # Get rates
                        table = t[w][s]
                        data = table.to_numpy()
                        energies = table.index.to_numpy()
                        # Initialize histogram if necessary
                        file_root = f"{flux_root}_{det}_events_{s}_{w}"
                        file_name = f"{Path(tarball_path).parent}/Collated_{file_root}.root"
                        energy_step = 0.5 * (energies[1] - energies[0])
                        if file_name not in outputnames.keys():
                            rootfile = ROOT.TFile.Open(file_name, "RECREATE")
                            outputnames[file_name] = rootfile
                            for channel in table.columns:
                                histname = f"{file_root}_{channel}" 
                                hist_dict[histname] = ROOT.TH2D(channel, f"Rates for {channel} interaction", 
                                                                 len(energies), energies[0] - energy_step, 
                                                                energies[-1] + energy_step, 
                                                                ntbins, tmin, tmax)
                                hist_dict[histname].SetDirectory(rootfile)
                        # Fill histogram with table
                        for energ,line in zip(energies,data):
                            for l,c in zip(line,table.columns):
                                histname = f"{file_root}_{c}"
                                hist_dict[histname].Fill(energ,time_val,l)
        # Save data and close ROOT files
        for fname in outputnames:
            outputnames[fname].Write()
            outputnames[fname].Close()
    return hist_dict
