#!/usr/bin/env python

from argparse import ArgumentParser
import to_snowglobes
import run_snowglobes
import from_snowglobes
import tarfile

SNOwGLoBES_path = "/location/of/snowglobes/" #where snowglobes is located
models_dir = "/location/of/models/Nakazato_2013/" #where models (aka input for to_snowglobes) is located
file_name = 'nakazato-LS220-BH-z0.004-s30.0.fits'

#input necessary for to_snowglobes
modeltype = 'Nakazato_2013'
transformation = 'AdiabaticMSW_NMO'
output = None
ntbins =  None
deltat = None

theta12 = 33.
theta13 = 9.
theta23 = 45.

m=None
tau=None

d=10

mixing_parameters = [theta12,theta13,theta23]    
decay_parameters = [m,tau,d]
parameters = mixing_parameters

#Running the modules

outfile = to_snowglobes.generate(models_dir, file_name, modeltype, transformation, parameters, output, ntbins, deltat, d) #runs to_snowglobes

#outfile = "ezyzip.zip"

run_snowglobes.go(SNOwGLoBES_path, models_dir, outfile)

from_snowglobes.collate(SNOwGLoBES_path, models_dir, outfile)

#an optional fourth argument in both run_ and from_ is the detector name, if one wants to only run 1 detector, rather than all of them
#fourth argument is the same except in the case of wc100kt30prct and ar40kt
#For wc100kt30prct: run_snowglobes takes in "wc100kt30prct", while from_snowglobes takes in "wc100kt30prct_eve"
#For ar40kt: run_snowglobes takes in "ar40kt", while from_snowglobes takes in "ar40kt_eve"
#All other times they are identical

