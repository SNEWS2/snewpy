
from snewpy import snowglobes
import numpy as np
from astropy import units as u

SNOwGLoBES_path = "/location/of/snowglobes/" #where snowglobes is located
SNEWPY_models_base = "/location/of/models/" #where models (aka input for to_snowglobes) is located
output_path = "/path/to/output/" #where the output files will be located

#set distance in kpc
distance=10

#set SNOwGLoBES detector to use
detector = "icecube"

#set SNEWPY model type and filename
modeltype = 'Tamborra_2014'
modeldir = SNEWPY_models_base+"/"+modeltype+"/"
model = 's20.0c_3D_dir1'
filename = model

#set desired flavor transformation prescription
transformation = 'AdiabaticMSW_NMO'

#snewpy.snowglobes creates a tarred file of snowglobes fluences
#this is stored in the model types directory in snewpy/models
outfile = modeltype+"_"+model+"_"+transformation

#Specify sequence of time intervals, one fluence
#file is made for each elements with dt=tstart[i]-tend[i]
window_tstart = 0.001
window_tend = 0.331
window_bins = 330
tstart = np.linspace(window_tstart,window_tend,window_bins,endpoint=False)*u.s
tend = tstart + (window_tend-window_tstart)/window_bins*u.s
tmid = (tstart+tend)*0.5

#Generate fluence file for SNOwGLoBES (there are two options, here we use generate_fluence)
print("Preparing fluences...")
tarredfile = snowglobes.generate_fluence(modeldir+filename, modeltype, transformation, distance, outfile,tstart,tend)
print("Done fluences...")

print("Running snowglobes...")
#now run SNOwGLoBES, this will loop over all the fluence files in `tarredfile`
snowglobes.simulate(SNOwGLoBES_path, tarredfile, detector_input=detector)
print("Done snowglobes...")

#now collate results of output of SNOwGLoBES
print("Collating...")
tables = snowglobes.collate(SNOwGLoBES_path, tarredfile, skip_plots=True)

#read results from SNOwGLoBES and put lightcurve in output file for snewpdag
fout = open(output_path+"snewpy_output_"+detector+"_"+modeltype+"_"+filename+"_1msbin.txt", "a")
nevents = np.zeros(len(tmid))
for i in range(len(tmid)):
    key = "Collated_"+outfile+"_"+str(i)+"_"+detector+"_events_smeared_weighted.dat"
    for j in range(1,len(tables[key]['header'].split())):
        nevents[i] += sum(tables[key]['data'][j])
    print(i, "\t" , nevents[i], file=fout)
