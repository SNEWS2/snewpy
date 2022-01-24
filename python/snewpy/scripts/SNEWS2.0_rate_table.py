#!/usr/bin/env python

import numpy as np
import os
from astropy.io import ascii
from snewpy import snowglobes

home_directory = os.getcwd()
SNOwGLoBES_path = "/path/to/snowglobes/"  # directory where SNOwGLoBES is located
SNEWPY_models_base = "/path/to/snewpy/models/"  # directory containing model input files

d = 10  # distance of supernova in kpc

dets = ["wc100kt30prct","ar40kt","halo1","halo2","scint20kt","novaFD","icecube"]
ref_mass = {"wc100kt30prct":100,"ar40kt":40,"halo1":0.079,"halo2":1,"scint20kt":20,"novaFD":14,"icecube":51600,"km3net":69366*3}

models = {}
models['s11.2'] = {'type':'Bollig_2016','file_name':'s11.2c'}
models['s27.0'] = {'type':'Bollig_2016','file_name':'s27.0c'}
models['s40'] = {'type':'OConnor_2015','file_name':'M1_neutrinos.dat'}

transformations = ['AdiabaticMSW_NMO','AdiabaticMSW_IMO']

total_events = {}

have_data_saved = False
if (have_data_saved is False):
    #Running the modules
    for model in models:
        total_events[model] = {}
        for transformation in transformations:
            total_events[model][transformation] = {}
            file_name = models[model]['file_name']
            modeltype = models[model]['type']
            outfile = modeltype+"_"+model+"_summed_"+transformation
            model_dir = SNEWPY_models_base+"/"+modeltype+"/"

            tarredfile = snowglobes.generate_fluence(model_dir + file_name, modeltype, transformation, d, outfile)
            for det in dets:
                snowglobes.simulate(SNOwGLoBES_path, tarredfile, detector_input=det)
                tables = snowglobes.collate(SNOwGLoBES_path, tarredfile, skip_plots=True)

                #for our table, interesting number is the smeared total number of events
                key = "Collated_"+outfile+"_"+det+"_events_smeared_weighted.dat"
                total_events[model][transformation][det+"smeared"]=0
                for j in range(1,len(tables[key]['header'].split())):
                    total_events[model][transformation][det+"smeared"] += sum(tables[key]['data'][j])

                key = "Collated_"+outfile+"_"+det+"_events_unsmeared_weighted.dat"
                total_events[model][transformation][det+"unsmeared"]=0
                for j in range(1,len(tables[key]['header'].split())):
                    total_events[model][transformation][det+"unsmeared"] += sum(tables[key]['data'][j])

    os.chdir(home_directory)
    np.save("SNEWS2.0_whitepaper_table_data.npy",total_events)
else:
    total_events = np.load("SNEWS2.0_whitepaper_table_data.npy",allow_pickle=True).tolist()
 

    
#Now lets make the table:
def round_to_2(x):
    if x==0: 
        return 0
    else:
        return round(x, -int(np.floor(np.log10(np.abs(x))))+1)

det_maps = {"Super-K":"wc100kt30prct","Hyper-K":"wc100kt30prct","IceCube":"icecube","LVD":"scint20kt",
            "KamLAND":"scint20kt","Borexino":"scint20kt","JUNO":"scint20kt","SNO+":"scint20kt","NO{$\\nu$}A":"novaFD",
            "HALO":"halo1","HALO-1kT":"halo2","DUNE":"ar40kt","MicroBooNe":"ar40kt","SBND":"ar40kt","Baksan":"scint20kt"}
    
data = {}
data['Experiment'] = ['Super-K','Hyper-K','IceCube','LVD','KamLAND','Borexino',
                       'JUNO','SNO+','NO{$\\nu$}A','Baksan','HALO','HALO-1kT','DUNE','MicroBooNe','SBND']

data['Type'] = ["H$_2$O/$\\bar{\\nu}_e$","H$_2$O/$\\bar{\\nu}_e$","String/$\\bar{\\nu}_e$",
                'C$_n$H$_{2n}$/$\\bar{\\nu}_e$','C$_n$H$_{2n}$/$\\bar{\\nu}_e$',
                'C$_n$H$_{2n}$/$\\bar{\\nu}_e$','C$_n$H$_{2n}$/$\\bar{\\nu}_e$',
                'C$_n$H$_{2n}$/$\\bar{\\nu}_e$','C$_n$H$_{2n}$/$\\bar{\\nu}_e$','C$_n$H$_{2n}$/$\\bar{\\nu}_e$',
                "Lead/$\\nu_e$","Lead/$\\nu_e$","Ar/$\\nu_e$","Ar/$\\nu_e$","Ar/$\\nu_e$"]
data['Mass [kt]'] = [32,220,51600,1,1,0.278,20,0.78,14,0.240,0.079,1,40,0.09,0.12]

data['Location'] = ["Japan","Japan","South Pole","Italy","Japan","Italy","China",
                    "Canada","USA","Russia","Canada","Italy","USA","USA","USA"]

data['\\SI{11.2}{\solarmass}'] = []
data['\\SI{27.0}{\solarmass}'] = []
data['\\SI{40.0}{\solarmass}'] = []

for experiment in range(len(data['Experiment'])):
    mass = data['Mass [kt]'][experiment]
    dettype = det_maps[data['Experiment'][experiment]]
    base_mass = ref_mass[dettype]
    
    counts_LCN = int(total_events['s11.2']['AdiabaticMSW_NMO'][dettype+"smeared"]*mass/base_mass)
    counts_LCI = int(total_events['s11.2']['AdiabaticMSW_IMO'][dettype+"smeared"]*mass/base_mass)
    counts_MCN = int(total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"smeared"]*mass/base_mass)
    counts_MCI = int(total_events['s27.0']['AdiabaticMSW_IMO'][dettype+"smeared"]*mass/base_mass)
    counts_HCN = int(total_events['s40']['AdiabaticMSW_NMO'][dettype+"smeared"]*mass/base_mass)
    counts_HCI = int(total_events['s40']['AdiabaticMSW_IMO'][dettype+"smeared"]*mass/base_mass)

    post = ['','','','','','']
    if counts_LCN>10000: 
        counts_LCN = int(counts_LCN/1000.0+0.5)
        post[0] = 'K'
    if counts_MCN>10000: 
        counts_MCN = int(counts_MCN/1000+0.5)
        post[2] = 'K'
    if counts_HCN>10000: 
        counts_HCN = int(counts_HCN/1000+0.5)
        post[4] = 'K'
    if counts_LCI>10000: 
        counts_LCI = int(counts_LCI/1000+0.5)
        post[1] = 'K'
    if counts_MCI>10000: 
        counts_MCI = int(counts_MCI/1000+0.5)
        post[3] = 'K'
    if counts_HCI>10000: 
        counts_HCI = int(counts_HCI/1000+0.5)
        post[5] = 'K'
        
    data['\\SI{11.2}{\solarmass}'].append(str(round_to_2(counts_LCN))+post[0]+"/"+str(round_to_2(counts_LCI))+post[1])
    data['\\SI{27.0}{\solarmass}'].append(str(round_to_2(counts_MCN))+post[2]+"/"+str(round_to_2(counts_MCI))+post[3])
    data['\\SI{40.0}{\solarmass}'].append(str(round_to_2(counts_HCN))+post[4]+"/"+str(round_to_2(counts_HCI))+post[5])


    
#Hacking in different numbers

#For IceCube, the effective mass in SNOwGLoBES is artificially high.  This is because the
#non-standard energy dependence is handled through the efficiencies.  To get an effective
#mass we take the ratio of the total weighted events to the unweighted events and multiply
#the unweighted mass (the entry in SNOwGLoBES), see below for details.  Here we take the
#effective mass of the s27 normal scenario and discuss the range in the table caption.

dettype='icecube'
mass=51600
data['Mass [kt]'][2] = "~"+str(100*int(0.01*(mass*total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"smeared"]/
      total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"unsmeared"]+50)))+"*"

#DarkSide-20k & Ar/any $\nu$ & 0.02 & Italy & \\
#For DarkSide-20k, the text on the SNEWS2.0 paper says the s27 will produce 250 events at 10kpc with a
#detector mass of 47 tonnes
s27_epT = 5319.1489361703

mass = 0.047
data['Experiment'].append('DarkSide-20k')
data['Type'].append('Ar/any $\\nu$')
data['Mass [kt]'].append(mass)
data['Location'].append('Italy')
data['\\SI{11.2}{\solarmass}'].append("-")
data['\\SI{27.0}{\solarmass}'].append(int(mass*s27_epT))
data['\\SI{40.0}{\solarmass}'].append("-")


#For XENONnT, LZ, and PandaX-4T, we take the values from Lang et al. [Phys. Rev. D 94 (2016) no.10, 103009].
#We take the S2-only, 60PE threshold numbers for the 11.2 msun and 27 msun models as they
#are the ones used here. We take 9.4 events/tonne for the 11.2 model and 17.6 events/tonne for the 27.0 model.
s11_epT = 9400.
s27_epT = 17600.

mass = 0.008
data['Experiment'].append('XENONnT')
data['Type'].append('Xe/any $\\nu$')
data['Mass [kt]'].append(mass)
data['Location'].append('Italy')
data['\\SI{11.2}{\solarmass}'].append(int(mass*s11_epT))
data['\\SI{27.0}{\solarmass}'].append(int(mass*s27_epT))
data['\\SI{40.0}{\solarmass}'].append("-")

mass = 0.007
data['Experiment'].append('LZ')
data['Type'].append('Xe/any $\\nu$')
data['Mass [kt]'].append(mass)
data['Location'].append('USA')
data['\\SI{11.2}{\solarmass}'].append(int(mass*s11_epT))
data['\\SI{27.0}{\solarmass}'].append(int(mass*s27_epT))
data['\\SI{40.0}{\solarmass}'].append("-")

mass = 0.004
data['Experiment'].append('PandaX-4T')
data['Type'].append('Xe/any $\\nu$')
data['Mass [kt]'].append(mass)
data['Location'].append('China')
data['\\SI{11.2}{\solarmass}'].append(int(mass*s11_epT))
data['\\SI{27.0}{\solarmass}'].append(int(mass*s27_epT))
data['\\SI{40.0}{\solarmass}'].append("-")


ascii.write(data,Writer=ascii.Latex)
