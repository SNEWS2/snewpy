#!/usr/bin/env python
""" This module allows to gather the SNOwGLoBES outputs into a ROOT files with 2D histograms.
    These 2D histograms show the rate as a function of time and energy for any given interaction channel.
    There is one ROOT file per detector, per smearing option, and per weighting option.
    Generate a tarball with the input files by running
            snowglobes.generate_time_series(...)
            snowglobes.simulate(...)
            snowglobes.collate(...)

    Then, to run the code, use:
        python Convert_to_ROOT.py <path to the tar file with collated outputs>
"""

import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse

import numpy as np

import ROOT

# Get the path to input files (command line argument). Default is current folder
parser = argparse.ArgumentParser()
parser.add_argument('path', help='Path to tarball containing collated SNEWPY outputs')
args = parser.parse_args()

tarball_path = args.path

# Iterate over collated SNEWPY files and gather all time bins
outputnames = {}
hist_dict = {}
tmin,tmax,ntbins = 0,0,0
column_names = []
with TemporaryDirectory(prefix='snowglobes') as tempdir:
    #Extracts data from tarfile and sets up lists of paths and fluxfilenames for later use
    with tarfile.open(tarball_path) as tar:
        tar.extractall(tempdir)

    flux_files = list(Path(tempdir).glob('*/Collated*.dat'))
    first = True
    for flux_file in flux_files:
        flux_root = flux_file.stem
        if first:
            # Time boundaries are formatted to be %0.3f in generate_time_series
            tmin,tmax,ntbins = flux_root.split(',')
            tmin = float(".".join(tmin.split('.')[-2:]))
            tmax = float(tmax)
            ntbins = int(ntbins.split('-')[0])
            print(f"Time binning is: tmin={tmin}, tmax={tmax}, ntbins={ntbins}")
            first = False
        # Get column names
        with open(flux_file) as fin:
            column_names = fin.readline().split()[1:]
        # Find general name for flux file as well as value of time bin
        flux_start, flux_end = flux_root.split('tbin')
        flux_split = flux_end.split('.')
        time_bin = int(flux_split[0])
        time_val = (time_bin - 0.5) * (tmax - tmin)/ntbins + tmin
        flux_end = ".".join(flux_split[1:])
        flux_root = flux_start + flux_end
        # Get data from file
        data = np.loadtxt(flux_file, skiprows=2)
        energies = data[:, 0]
        data = data[:, 1:]
        # Initialize histograms and open files
        file_name = f"{Path(tarball_path).parent}/{flux_root}.root"
        energy_step = 0.5 * (energies[1] - energies[0])
        if file_name not in outputnames.keys():
            print(f"Creating file {file_name}")
            rootfile = ROOT.TFile.Open(file_name, "RECREATE")
            outputnames[file_name] = rootfile
            for channel in column_names:
                histname = f"{flux_root}_{channel}" 
                hist_dict[histname] = ROOT.TH2D(channel, f"Rates for {channel} interaction", 
                                                    len(energies), energies[0] - energy_step, 
                                                energies[-1] + energy_step, 
                                                ntbins, tmin, tmax)
                hist_dict[histname].SetDirectory(rootfile)
        # Fill histogram with table
        for energ,line in zip(energies,data):
            for l,c in zip(line,column_names):
                histname = f"{flux_root}_{c}"
                hist_dict[histname].Fill(energ,time_val,l)
# Save data and close ROOT files
print("Closing files...")
for fname in outputnames:
    outputnames[fname].Write()
    outputnames[fname].Close()
