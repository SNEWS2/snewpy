# -*- coding: utf-8 -*-
"""
snewpy.scripts.to_snowglobes
============================

Convert an arbitrary model to SNOwGLoBES format. Based on SNEWPY.py script by
E. O'Connor and J. P. Kneller.
"""

import numpy as np
from argparse import ArgumentParser

from snewpy.models_class import *
from snewpy.FlavorTransformation import *

def OutputInSnowGlobesFormat(SN):
    d=10. *1000.*3.086e+18
    keV=1e3 * 1.60218e-12
    MeV=1e6 * 1.60218e-12
    GeV=1e9 * 1.60218e-12
    
    times=SN.get_time()
    timebin=0.01
      
    for j in range(0,1000):
        t=timebin*j
        
        if(t<=times[0] or t>=times[-1]):
            break
            
        print(t)
        
        filename  =  SN.filename.split(".fits")[0] + "-t=" + str(t) + "-dt=" + str(timebin) + ".SNOformat.dat" 
        file = open(filename,"w")

        E=np.linspace(0,100,501)*MeV
        OscillatedSpectra = SN.get_oscillatedspectra(t,E) 

        OscillatedFluence={}
        for i in range(0,501):
            for flavor in Flavor:
                OscillatedFluence[flavor]=OscillatedSpectra[flavor][i] * timebin * 200.*keV  / (4.*np.pi*d**2)

            file.write("{0:.4f}".format(E[i]/GeV))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_e]))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x]/2.))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x]/2.))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_e_bar]))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x_bar]/2.))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x_bar]/2.)+"\n")
        file.close()


def main(options=None):
    p = ArgumentParser(description='Convert to SNOwGLoBES format.')
    p.add_argument('infile', nargs=1,
                   help='Supernova model input file (Nakazato only).')

    if options is None:
        args = p.parse_args()
    else:
        args = p.parse_args(options)

    SNModel = Nakazato2013(args.infile[0], AdiabaticMSW_NMO())

    OutputInSnowGlobesFormat(SNModel)
