import numpy as np
import os

from snewpy.rate_calculator import RateCalculator
from snewpy.models.ccsn import Bollig_2016, OConnor_2015
from snewpy.neutrino import MassHierarchy, MixingParameters, ThreeFlavorMixingParameters
from snewpy.flavor_transformation import AdiabaticMSW

from astropy import units as u

#Select output format, Mathjax or LaTeX
while True:
    outputFormat = input('Enter 1 for Mathjax format or 0 for LaTeX format: ')
    if outputFormat == '1':
        useMathjax = True
        break
    elif outputFormat == '0':
        useMathjax = False
        break
    else:
        print('Please enter 1 or 0.')

models = { 's11.2' : Bollig_2016(progenitor_mass=11.2<<u.Msun), 
           's27.0' : Bollig_2016(progenitor_mass=27<<u.Msun),  
           's40' : OConnor_2015(progenitor_mass=40<<u.Msun) }

transformations = { 'AdiabaticMSW_NMO' : AdiabaticMSW(MixingParameters('NORMAL')), 
                    'AdiabaticMSW_IMO' : AdiabaticMSW(MixingParameters('INVERTED')) }

detectors = ["wc100kt30prct", "ar40kt", "halo1", "halo2", "scint20kt", "novaFD", 
             "icecube", "km3net", "ds20", "xent", "lz", "pandax"]
        
detector_masses = {"wc100kt30prct": 100, "ar40kt": 40, "halo1": 0.079, "halo2": 1, "scint20kt": 20, "novaFD": 14,
                   "icecube": 51600, "km3net": 69366 * 3, "ds20": 0.0386, "xent": 0.006, "lz": 0.007, "pandax": 0.004}

detector_effects = {'smeared' : True, 'unsmeared' : False}

rc=RateCalculator()

total_events = {}

energies = np.linspace(0,100,501)<<u.MeV
distance = 10*u.kpc
        
# Running the modules
for effects in detector_effects:
    total_events[effects] = {}
    for model in models:
        times = models[model].get_time()
        total_events[effects][model] = {}
        for transformation in transformations:
            total_events[effects][model][transformation] = {}
            #get the flux from the model
            flux = model.get_flux(t=times, E=energies, distance=distance, flavor_xform=transformation)
            fluence = flux.integrate('time')
            for detector in detectors:
                events = rc.run(fluence, detectors[detector], detector_effects=detector_effects[effects])
                total_events[effects][model][transformation][detector] = sum([chan.integrate_or_sum('energy').array.squeeze().value for chan in events.values()])     

home_directory = os.getcwd()
os.chdir(home_directory)
np.savez("SNEWS2.0_whitepaper_table_data.npz", total_events)

# Now lets make the table:
def round_to_2(x):
    if x == 0:
        return 0
    else:
        return round(x, -int(np.floor(np.log10(np.abs(x)))) + 1)


detector_maps = {"Super-K": "wc100kt30prct", "Hyper-K": "wc100kt30prct", "IceCube": "icecube", "KM3NeT":"km3net",
                 "LVD": "scint20kt", "KamLAND": "scint20kt", "Borexino": "scint20kt", "JUNO": "scint20kt",
                 "SNO+": "scint20kt", "NO${\\nu}$A": "novaFD", "HALO": "halo1", "HALO-1kT": "halo2", "DUNE": "ar40kt",
                 "MicroBooNe": "ar40kt", "SBND": "ar40kt", "Baksan": "scint20kt", "DarkSide-20k": "ds20", "XENONnT": "xent",
                 "LZ": "lz", "PandaX-4T": "pandax"}

data = {}
data['Experiment'] = ['Super-K', 'Hyper-K', 'IceCube', 'KM3NeT', 'LVD', 'KamLAND', 'Borexino',
                      'JUNO', 'SNO+', 'NO${\\nu}$A', 'Baksan', 'HALO', 'HALO-1kT', 'DUNE', 'MicroBooNe', 'SBND',
                      'DarkSide-20k', 'XENONnT', 'LZ', 'PandaX-4T']

data['Type'] = ["\\text{H$_2$O$/\\bar{\\nu}_e$}", "\\text{H$_2$O$/\\bar{\\nu}_e$}", "\\text{String}/\\bar{\\nu}_e",
                "\\text{String}/\\bar{\\nu}_e", '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}', '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}',
                '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}', '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}',
                '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}', '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}', '\\text{C$_n$H$_{2n}/\\bar{\\nu}_e$}',
                "\\text{Lead/$\\nu_e$}", "\\text{Lead/$\\nu_e$}", "\\text{Ar/$\\nu_e$}", "\\text{Ar/$\\nu_e$}", "\\text{Ar/$\\nu_e$}", "\\text{Ar/any $\\nu$}", "\\text{Xe/any $\\nu$}",
                "\\text{Xe/any $\\nu$}", "\\text{Xe/any $\\nu$}"]
data['Mass [kt]'] = [32, 220, 51600, 69366*3, 1, 1, 0.278, 20, 0.78, 14, 0.240, 0.079, 1, 40, 0.09, 0.12, 0.0386, 0.006,
                     0.007, 0.004]

data['Location'] = ["Japan", "Japan", "South Pole", "Italy", "Italy", "Japan", "Italy", "China",
                    "Canada", "USA", "Russia", "Canada", "Italy", "USA", "USA", "USA", "Italy", "Italy", "USA", "China"]

data['11.2 M$_\\odot$'] = []
data['27.0 M$_\\odot$'] = []
data['40.0 M$_\\odot$'] = []


for experiment in range(len(data['Experiment'])):
    mass = data['Mass [kt]'][experiment]
    detector_type = detector_maps[data['Experiment'][experiment]]
    base_mass = detector_masses[detector_type]

    counts_LCN = int(total_events['smeared']['s11.2']['AdiabaticMSW_NMO'][detector_type] * mass / base_mass)
    counts_LCI = int(total_events['smeared']['s11.2']['AdiabaticMSW_IMO'][detector_type] * mass / base_mass)
    counts_MCN = int(total_events['smeared']['s27.0']['AdiabaticMSW_NMO'][detector_type] * mass / base_mass)
    counts_MCI = int(total_events['smeared']['s27.0']['AdiabaticMSW_IMO'][detector_type] * mass / base_mass)
    counts_HCN = int(total_events['smeared']['s40']['AdiabaticMSW_NMO'][detector_type] * mass / base_mass)
    counts_HCI = int(total_events['smeared']['s40']['AdiabaticMSW_IMO'][detector_type] * mass / base_mass)

    post = ['', '', '', '', '', '']
    if counts_LCN > 10000:
        counts_LCN = int(counts_LCN / 1000.0 + 0.5)
        post[0] = '\\text{K}'
    if counts_MCN > 10000:
        counts_MCN = int(counts_MCN / 1000 + 0.5)
        post[2] = '\\text{K}'
    if counts_HCN > 10000:
        counts_HCN = int(counts_HCN / 1000 + 0.5)
        post[4] = '\\text{K}'
    if counts_LCI > 10000:
        counts_LCI = int(counts_LCI / 1000 + 0.5)
        post[1] = '\\text{K}'
    if counts_MCI > 10000:
        counts_MCI = int(counts_MCI / 1000 + 0.5)
        post[3] = '\\text{K}'
    if counts_HCI > 10000:
        counts_HCI = int(counts_HCI / 1000 + 0.5)
        post[5] = '\\text{K}'

    data['11.2 M$_\\odot$'].append(
        str(round_to_2(counts_LCN)) + post[0] + "/" + str(round_to_2(counts_LCI)) + post[1])
    data['27.0 M$_\\odot$'].append(
        str(round_to_2(counts_MCN)) + post[2] + "/" + str(round_to_2(counts_MCI)) + post[3])
    data['40.0 M$_\\odot$'].append(
        str(round_to_2(counts_HCN)) + post[4] + "/" + str(round_to_2(counts_HCI)) + post[5])

# For IceCube & KM3NeT, the effective masses in SNOwGLoBES are artificially high.  This is because the
# non-standard energy dependence is handled through the efficiencies.  To get an effective
# mass we take the ratio of the total weighted events to the unweighted events and multiply
# the unweighted mass (the entry in SNOwGLoBES), see below for details.  Here we take the
# effective mass of the s27 normal scenario and discuss the range in the table caption.

detector_type = 'icecube'
mass = 51600
data['Mass [kt]'][2] = "~"+str(int(round(mass*total_events['smeared']['s27.0']['AdiabaticMSW_NMO'][detector_type]/
      total_events['unsmeared']['s27.0']['AdiabaticMSW_NMO'][detector_type], -2)))+"*"

detector_type = 'km3net'
mass = 69366 * 3
data['Mass [kt]'][3] = "~"+str(int(round(mass*total_events['smeared']['s27.0']['AdiabaticMSW_NMO'][detector_type]/
      total_events['unsmeared']['s27.0']['AdiabaticMSW_NMO'][detector_type], -1)))+"*"

# Formatting the dictionary to be compatible with LaTeX & MathJax (useful for html)
def dictArray(dictionary):
    table = r'\begin{array} {|r|c|r|c|c|c|c|}\hline '

    # Writes & formats the header row
    is_first_column = True
    for columnHeader in dictionary:
        if not is_first_column:
            table = table + ' & '
        table = table + '\\text{' + columnHeader + '}'
        is_first_column = False
    table = table+r'\\ \hline '

    # Writes the center rows
    tempVal = list(dictionary)[0]   # Calls an arbitrary key's stored data (list)
    numDetectors = len(dictionary[tempVal])  # Uses the list's length to calculate the number of rows needed (num detectors)

    for det in range(numDetectors):    # Iterates through all data in a row/ for each detector
        is_first_column = True
        columnPos = 0   #tracks which column of the row the loop is in
        for col in dictionary:  # Iterates through all column titles/ dictionary keys
            if not is_first_column:
                table = table + ' & '
            if columnPos == 0 or columnPos == 3:
                table = table + '\\text{' + str(dictionary[col][det]) + '}'
            if columnPos != 0 and columnPos != 3:
                table = table + str(dictionary[col][det])
            is_first_column = False
            columnPos += 1
        table = table + r'\\ \hline '

    table = table + r' \end{array}'

    return table

# Prints out code for html that puts the data into a table (MathJax array)
if useMathjax:
    print('<body>')
    print('  <script id="MathJax-script" async')
    print('          src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">\n  </script>\n<p>')
    print(dictArray(data))
    print('\n</p>\n</body>')
# Just copy & paste the output into any html page to insert the table!

# Prints out data in a LaTeX table (array) script
else:
    print("\\usepackage{amsmath}\n\\begin{document}")   #loads forever without amsmath package
    print(f"$$\n{dictArray(data)}\n$$")
    print("\\end{document}")
#Just copy & paste the output into LaTeX!
