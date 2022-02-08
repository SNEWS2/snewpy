#!/usr/bin/env python
from snewpy import snowglobes

SNOwGLoBES_path = "/path/to/snowglobes/"  # directory where SNOwGLoBES is located
SNEWPY_model_dir = "/path/to/snewpy/models/"  # directory containing model input files

distance = 10  # Supernova distance in kpc
detector = "wc100kt30prct" #SNOwGLoBES detector for water Cerenkov
modeltype = 'Bollig_2016' # Model type from snewpy.models
model = 's11.2c' # Name of model
transformation = 'AdiabaticMSW_NMO' # Desired flavor transformation

# Construct file system path of model file and name of output file
model_path = SNEWPY_model_dir + "/" + modeltype + "/" + model
outfile = modeltype + "_" + model + "_" + transformation

# Now, do the main work:
print("Generating fluence files ...")
tarredfile = snowglobes.generate_fluence(model_path, modeltype, transformation, distance, outfile)

print("Simulating detector effects with SNOwGLoBES ...")
snowglobes.simulate(SNOwGLoBES_path, tarredfile, detector_input=detector)

print("Collating results ...")
tables = snowglobes.collate(SNOwGLoBES_path, tarredfile, skip_plots=True)


# Use results to print the number of events in different interaction channels
key = f"Collated_{outfile}_{detector}_events_smeared_weighted.dat"
total_events = 0
for i, channel in enumerate(tables[key]['header'].split()):
    if i == 0:
        continue
    n_events = sum(tables[key]['data'][i])
    total_events += n_events
    print(f"{channel:10}: {n_events:.3f} events")

#Super-K has 32kT inner volume
print("Total events in Super-K-like detector:",0.32*total_events)

