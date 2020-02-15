# -*- coding: utf-8 -*-
"""
snewpy.scripts.to_snowglobes
============================

Convert an arbitrary model to SNOwGLoBES format. Based on SNEWPY.py script by
E. O'Connor and J. P. Kneller.

This version will subsample the times in a supernova model, produce energy
tables expected by SNOwGLoBES, and compress the output into a tarfile.
"""

import numpy as np
from argparse import ArgumentParser

import os
import io
import tarfile

import logging

from snewpy.models import *
from snewpy.FlavorTransformation import *

def main(options=None):
    # Parse command-line arguments.
    p = ArgumentParser(description='Convert to SNOwGLoBES format.')
    p.add_argument('infile', nargs=1,
                   help='Supernova model input file (Nakazato only).')
    p.add_argument('-o', '--output', default=None,
                   help='Output tarfile name (if customization desired)')

    tbingroup = p.add_mutually_exclusive_group()
    tbingroup.add_argument('-n', '--ntbins', type=int,
                           help='Number of bins used to sample model')
    tbingroup.add_argument('-t', '--deltat', type=float,
                           help='Time binning used to sample model [sec]')

    p.add_argument('-v', '--verbose', action='store_true', default=False,
                   help='Activate verbose log for debugging')

    if options is None:
        args = p.parse_args()
    else:
        args = p.parse_args(options)

    # Set verbosity of the log.
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # Load up the model. To do: support more than Nakazato format.
    infile = args.infile[0]
    snmodel = Nakazato2013(infile, NoTransformation())

    # Subsample the model time. Default to 30 time slices.
    tmin = snmodel.get_time()[0]
    tmax = snmodel.get_time()[-1]
    if args.deltat is not None:
        dt = args.deltat
    elif args.ntbins is not None:
        dt = (tmax - tmin) / (args.ntbins+1)
    else:
        dt = (tmax - tmin) / 31

    tedges = np.arange(tmin, tmax, dt)
    times = 0.5*(tedges[1:] + tedges[:-1])

    # Generate output.
    if args.output is not None:
        tfname = args.output
    else:
        tfname = infile.replace('.fits', '.SNOformat.tar.bz2')

    with tarfile.open(tfname, 'w:bz2') as tf:
        d = 10. *1000.*3.086e+18       # luminosity to fluence
        keV = 1e3 * 1.60218e-12        # eV to erg
        MeV = 1e6 * 1.60218e-12
        GeV = 1e9 * 1.60218e-12

        energy = np.linspace(0, 100, 501) * MeV

        # Loop over sampled times.
        for i, t in enumerate(times):
            osc_spectra = snmodel.get_oscillatedspectra(t, energy)
            osc_fluence = {}
            table = []

            table.append('# GCD@TBinMid={:g}sec@(tBinWidth={:g}s)(eBinWidth=0.2MeV) Flux in Number Neutrinos per cm^2'.format(t, dt))
            table.append('# E(GeV)            NuE               NuMu NuTau             aNuE              aNuMu             aNuTau')

            # Generate energy + number flux table.
            for j, E in enumerate(energy):
                for flavor in Flavor:
                    osc_fluence[flavor] = osc_spectra[flavor][j] * dt * 200.*keV  / (4.*np.pi*d**2)
                
                s = '{:17.8E}'.format(E/GeV)
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_e])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x]/2)
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x]/2)
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_e_bar])
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x_bar]/2)
                s = '{}{:17.8E}'.format(s, osc_fluence[Flavor.nu_x_bar]/2)
                table.append(s)
                logging.debug(s)

            # Encode energy/flux table and output to file in tar archive.
            output = '\n'.join(table).encode('ascii')

            infoname = '{:02d}Tbins/{}-tbin{:02d}.NoOsc.dat'.format(
                           len(times),
                           os.path.basename(infile).replace('.fits', ''),
                           i + 1)
            info = tarfile.TarInfo(name=infoname)
            info.size = len(output)

            logging.info('Time {:g} s; writing {} to {}'.format(t, infoname, tfname))
            tf.addfile(info, io.BytesIO(output))

