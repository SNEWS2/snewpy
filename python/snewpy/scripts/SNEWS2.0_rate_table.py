import numpy as np
import os
from snewpy import snowglobes

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

home_directory = os.getcwd()
SNOwGLoBES_path = None  # change to SNOwGLoBES directory if using a custom detector configuration
SNEWPY_models_base = "Your/SNEWPY/Models/Path/Here"  # local directory containing model input files ("SNEWPY_models")

d = 10  # distance of supernova in kpc

dets = ["wc100kt30prct", "ar40kt", "halo1", "halo2", "scint20kt", "novaFD", "icecube", "km3net", "ds20", "xent", "lz",
        "pandax"]
ref_mass = {"wc100kt30prct": 100, "ar40kt": 40, "halo1": 0.079, "halo2": 1, "scint20kt": 20, "novaFD": 14,
            "icecube": 51600, "km3net": 69366 * 3, "ds20": 0.0386, "xent": 0.006, "lz": 0.007, "pandax": 0.004}

models = {}
models['s11.2'] = {'type': 'Bollig_2016', 'file_name': 's11.2c'}
models['s27.0'] = {'type': 'Bollig_2016', 'file_name': 's27.0c'}
models['s40'] = {'type': 'OConnor_2015', 'file_name': 'M1_neutrinos.dat'}

transformations = ['AdiabaticMSW_NMO', 'AdiabaticMSW_IMO']

total_events = {}

have_data_saved = False
if (have_data_saved is False):
    # Running the modules
    for model in models:
        total_events[model] = {}
        for transformation in transformations:
            total_events[model][transformation] = {}
            file_name = models[model]['file_name']
            modeltype = models[model]['type']
            outfile = modeltype + "_" + model + "_summed_" + transformation
            model_dir = SNEWPY_models_base + "/" + modeltype + "/"

            tarredfile = snowglobes.generate_fluence(model_dir + file_name, modeltype, transformation, d, outfile)
            for det in dets:
                snowglobes.simulate(SNOwGLoBES_path, tarredfile, detector_input = det)
                tables = snowglobes.collate(SNOwGLoBES_path, tarredfile, skip_plots = True)

                # for our table, interesting number is the smeared total number of events
                key = "Collated_" + outfile + "_" + det + "_events_smeared_weighted.dat"
                total_events[model][transformation][det + "smeared"] = 0
                for j in range(1, len(tables[key]['header'].split())):
                    total_events[model][transformation][det + "smeared"] += sum(tables[key]['data'][j])

                key = "Collated_" + outfile + "_" + det + "_events_unsmeared_weighted.dat"
                total_events[model][transformation][det + "unsmeared"] = 0
                for j in range(1, len(tables[key]['header'].split())):
                    total_events[model][transformation][det + "unsmeared"] += sum(tables[key]['data'][j])

    os.chdir(home_directory)
    np.save("SNEWS2.0_whitepaper_table_data.npy", total_events)
else:
    total_events = np.load("SNEWS2.0_whitepaper_table_data.npy", allow_pickle = True).tolist()


# Now lets make the table:
def round_to_2(x):
    if x == 0:
        return 0
    else:
        return round(x, -int(np.floor(np.log10(np.abs(x)))) + 1)


det_maps = {"Super-K": "wc100kt30prct", "Hyper-K": "wc100kt30prct", "IceCube": "icecube", "KM3NeT":"km3net",
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
    dettype = det_maps[data['Experiment'][experiment]]
    base_mass = ref_mass[dettype]

    counts_LCN = int(total_events['s11.2']['AdiabaticMSW_NMO'][dettype + "smeared"] * mass / base_mass)
    counts_LCI = int(total_events['s11.2']['AdiabaticMSW_IMO'][dettype + "smeared"] * mass / base_mass)
    counts_MCN = int(total_events['s27.0']['AdiabaticMSW_NMO'][dettype + "smeared"] * mass / base_mass)
    counts_MCI = int(total_events['s27.0']['AdiabaticMSW_IMO'][dettype + "smeared"] * mass / base_mass)
    counts_HCN = int(total_events['s40']['AdiabaticMSW_NMO'][dettype + "smeared"] * mass / base_mass)
    counts_HCI = int(total_events['s40']['AdiabaticMSW_IMO'][dettype + "smeared"] * mass / base_mass)

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

dettype = 'icecube'
mass = 51600
data['Mass [kt]'][2] = "~"+str(int(round(mass*total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"smeared"]/
      total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"unsmeared"], -2)))+"*"

dettype = 'km3net'
mass = 69366 * 3
data['Mass [kt]'][3] = "~"+str(int(round(mass*total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"smeared"]/
      total_events['s27.0']['AdiabaticMSW_NMO'][dettype+"unsmeared"], -1)))+"*"

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