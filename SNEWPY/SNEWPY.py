'''
  Copyright (c) 2020 James Kneller

  This file is part of SNEWPY.

  SNEWPY is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  SNEWPY is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
'''

from abc import abstractmethod, ABC
from enum import Enum

from astropy.table import Table

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import scipy
import math

from models_class import *
from FlavorTransformation import *

SNModel = Nakazato2013('nakazato-LS220-BH-z0.004-s30.0.fits',AdiabaticMSW_NMO())

def OutputInSnowGlobesFormat(SN):
    d=10. *1000.*3.086e+18
    keV=1e3 * 1.60218e-12
    GeV=1e9 * 1.60218e-12
    times=SN.get_time()
        
    for j in range(0,1000):
        timebin=0.01
        t=timebin*j
        
        if(t<=times[0] or t>=times[-1]):
            break
            
        print(t)
        
        filename  =  SN.filename.split(".fits")[0] + "-t=" + str(t) + "-dt=" + str(timebin) + ".SNOformat.dat" 
        file = open(filename,"w")

        OscillatedFluence={}
        for i in range(0,501):
            E=i*200.*keV
            OscillatedSpectra = SN.get_oscillatedspectra(t,E)
            
            for flavor in Flavor:
                OscillatedFluence[flavor]=OscillatedSpectra[flavor] * 200.*keV  / (4.*np.pi*d**2)

            file.write("{0:.4f}".format(E/GeV))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_e]))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x]/2.))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x]/2.))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_e_bar]))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x_bar]/2.))
            file.write("\t"+str(OscillatedFluence[Flavor.nu_x_bar]/2.)+"\n")
        file.close()

OutputInSnowGlobesFormat(SNModel)


