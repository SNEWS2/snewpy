#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 23:16:27 2022
This contains data handlers for handling individual model data
@author: phyics
"""

import numpy as np
import os

supported_detectors = [
    'scint20kt',
    'ar40kt',
    'wc100kt30prct'
    ]

def build_detector_profiles():
    '''
    Builds a profile for a detector that includes a dictionary of data handlers
    and the corresponding data labels (used for labeling axes)

    Returns
    -------
    detector_defs : dictionary
        The dictionary of profiles

    '''
    detector_defs = {}
    module_scope = __import__(__name__)
    for detector in supported_detectors:
        detector_defs[detector] = {
            'handler':getattr(module_scope, f'h_{detector}'),
            'axes':getattr(module_scope,f'axes_{detector}')
            }
    return detector_defs

'''
What follows are the calculations that are made for each time bin when collated
data comes through. Each time bin has event rates per energy bin. Order should
be nc, nu_e+es, then inverse beta decay events.

Treat all the NC as nux,  all the nue-nucleus and ES as nue, and IBD+nuebar-nucleus as nuebar
'''

def h_scint20kt(data):
    # must return a list of a, b, c
    nue_plus_es=np.sum(data['nue_C12'])+np.sum(data['nue_C13']+data['e']) # nue
    ibd = np.sum(data['ibd'])
    nc = np.sum(data['nc'])
    return [nc,nue_plus_es,ibd]

def h_ar40kt(data):
    nue_plus_es = np.sum(data['nue_Ar40'])+np.sum(data['e'])
    nc = np.sum(data['nc'])
    nue_bar = np.sum(data['nuebar_Ar40'])
    return [nc,nue_plus_es,nue_bar]

def h_wc100kt30prct(data):
    ibd = np.sum(data['ibd'])
    nue_plus_es=np.sum(data['nue_O16'])+np.sum(data['e'])
    nc = np.sum(data['nc'])
    return [nc,nue_plus_es,ibd]

'''
Data handler axes labels definitions
'''

def same_axes():
    return [r'$\nu_x$',r'$\nu_e$',r'$\bar{\nu_e}$']

def axes_scint20kt():
    return same_axes()

def axes_ar40kt():
    return same_axes()

def axes_wc100kt30prct():
    return same_axes()