#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 22:37:52 2022
Script used for plotting individual time slice data
@author: phyics
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.models import Nakazato_2013,Walk_2018
from snewpy.flavor_transformation import NoTransformation # just use NoTransformation for now to keep things simple
from snowglobes import generate_time_series
import os
import data_handlers as handlers

import sys
sys.path.insert(0,'./SURF2020fork')

#simulation details
modelFilePathBase = "../..//models/Nakazato_2013/"
modelFilePath = modelFilePathBase + "nakazato-shen-z0.004-t_rev100ms-s20.0.fits"
model = Nakazato_2013(modelFilePath)
model_type="Nakazato_2013"
step_size = 1 # 0.04 for better results
deltat=step_size*u.s
detector = 'wc100kt30prct'
d = 10 # in pc, distance to SN
snowglobes_out_name="snowglobes-output"
snowglobes_dir = os.environ['SNOWGLOBES']
print(os.environ['SNOWGLOBES'])
smearing = 'smeared'
model
transform = "NoTransformation"
print(f'Timeframe of model is from {model.time[0]} to {model.time[len(model.time)-1]}')

profiles = handlers.build_detector_profiles()

generate_time_series(modelFilePath,model_type,transform,d,log_bins=False)



# l_plot_data, l_raw_data, l_l_data = t.create_detector_event_scatter(
#     modelFilePath,model_type,
#     detector,
#     model,
#     deltat=deltat,
#     transformation=transform,
#     data_calc=profiles[detector]['handler'],
#     use_cache=True,
#     log_bins=True
#     )
    
    # t.create_regular_plot(l_raw_data, ['ibd','nue+es','nc'],
    #                       '{model} {detector} {transform} Logged Bins'.format(model=model_type,detector=detector,transform=transform),
    #                       ylab="Event Counts",xlab="Logged Time Bin No",use_x_log=False,save=True,show=True)
    # t.create_default_detector_plot(l_plot_data, profiles[detector]['axes'](),
    #                                f'Nakazato {detector} {transform} dt={str(deltat)} Logged Bins',show=True)
    # # also create regular plot for comparison
    # plot_data, raw_data, l_data = t.create_detector_event_scatter(
    #     modelFilePath,model_type,
    #     detector,
    #     model,
    #     deltat=deltat,
    #     transformation=transform,
    #     data_calc=profiles[detector]['handler'],
    #     use_cache=True
    #     )
    # t.create_default_detector_plot(plot_data, profiles[detector]['axes'](), f'Nakazato {detector} {transform} dt={str(deltat)}',save=True)