#!/usr/bin#!/usr/bin/env python
from __future__ import unicode_literals
import os
import tarfile
import zipfile
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import fnmatch
import logging
import numpy as np

logging.getLogger().setLevel(logging.CRITICAL)

# Anne Graf
# Written & modulized Summer 2020
# Takes in input flux files and configures and runs supernova (which outputs calculated rates)
            
def collate(Branch, Model_Path, Tarball, detector_input = all, skip_plots=False, return_tables=False, verbose=False, remove_generated_files=True):
    #Determines type of input file

    if ".tar.bz2" in str(Tarball):
        outputnamestem = Tarball[0:str(Tarball).rfind(".tar.bz2")]
        tar = tarfile.open(Model_Path+"/"+Tarball)
        TarballFileNames = tar.getnames()
        tar.close()
    elif ".tar.gz" in str(Tarball):
        outputnamestem = Tarball[0:str(Tarball).rfind(".tar.gz")]
        tar = tarfile.open(Model_Path+"/"+Tarball)
        TarballFileNames = tar.getnames()
        tar.close()
    elif ".zip" in str(Tarball):
        outputnamestem = Tarball[0:str(Tarball).rfind(".zip")]
        zip = zipfile.ZipFile(Model_Path+"/"+Tarball)
        TarballFileNames = zip.namelist()
        zip.close()
    else:
        print("Invalid Tar file")


    #Determining the configuration of the files and directories within the tarfile and defining variables for uniformity

    FluxFileNameStems = []
    extension = '.dat'
    for IndividualFile in TarballFileNames:
        extension_beginning = str(IndividualFile).rfind(".")
        if extension_beginning > -1:
            extension = str(IndividualFile[extension_beginning:])
            ExtensionFound = str(IndividualFile).rfind(extension)           
            FluxFileNameStems.append( str(IndividualFile)[0:ExtensionFound] ) # gets rid of extension at the end of flux files            
        else:
            continue

    smearvals = ["_smeared_w", "_smeared_u", "nts_w", "nts_u"]
    #weightvals = ["_weighted", "unweighted"]

    FilesToCleanup = []

    #Add_funct sums up relevant files from output generated above
    def add_funct(flux, detector, smear, *arg):
        homebase = Branch + "/out/"
        
        #Defining different variables for the combinatoric iterations
        #And reformatting some variables for different naming uses
        if smear == "_smeared_w" or smear == "nts_w":
            weight_value = "weighted"
        else:
            weight_value = "unweighted"
        
        if smear == "_smeared_w" or smear == "_smeared_u":
            smear_value = "smeared"
            smear_title = "detected"
            x_label = "Detected Energy (GeV)"
            y_label = "Events"
        else:
            smear_value = "unsmeared"
            smear_title = "interaction"
            x_label = "Neutrino Energy (GeV)"
            y_label = "Interaction Events"

        if detector == "ar40kt_eve":
            detector_label = "ar40kt"
            detector_label2 = "ar40kt"
        elif detector == "ar40kt_he":
            detector_label = "ar40kt he"
            detector_label2 = "ar40kt_he"
        elif detector == "wc100kt30prct_eve":
            detector_label = "wc100kt30prct"
            detector_label2 = "wc100kt30prct"
        elif detector == "wc100kt30prct_he":
            detector_label = "wc100kt30prct he"
            detector_label2 = "wc100kt30prct_he"
        else:
            detector_label = detector
            detector_label2 = detector

        colors = [ "k", "r", "g", "y", "b", "m", "c"]        #colors of graphed values
        compile_dict = {}

        #Splits values from each file into Energy and Events, then sums appropriate ones
        for input_val in arg:
            flux_filenames = os.listdir(homebase)
            final_dict= {}    
            temp_dict = {}
            for afile in flux_filenames:                #Does the actual summing of essential values in each file
                name_iteration = "{0}_{1}*{2}*{3}*".format(flux, input_val, detector, smear)
                if fnmatch.fnmatch(afile, name_iteration): # and "Collated" not in str(afile):
                    FilesToCleanup.append(homebase+"/"+afile)
                    fileinsides = open("{0}/{1}".format(homebase, afile)) #fine
                    
                    lines = fileinsides.readlines()
                    lines.pop(-1) #removes empty lines
                    lines.pop(-1)
                    lines.pop(-1)
                    if "unweight" in str(afile):
                        lines.pop(0)
                        
                    #Splits each line & converts Energy and Events into floats
                    for line in lines:
                        if "." in line:
                        
                            #Determining which characters in the line compose the energy bin and corresponding event rate
                            first_char = line.find("0") #finds first character in a line
                            end_energy = line.find(" ", first_char, -1) #finds the last character of the energy value
                            if end_energy == -1:
                                end_energy = line.find("\t", first_char, -1) #alternative last character of the energy value
                            event_found = line.find(".", end_energy, -1) #determines approximately where event value is located in the line
                            begin_event = line.rfind(" ", end_energy, event_found) #first character of event value
                            end_event = line.find(" ", event_found, -1) #last character of event value
                            if end_event == -1:
                                end_event = len(line) #alternative last character of event value
                                begin_event = end_energy + 1 #alternative first character of event value
                            
                            Energy = float((line[first_char:end_energy])) #Complete energy characters
                            Events = float((line[begin_event + 1: (end_event)])) #Complete event characters
                            temp_dict[Energy] = Events
                            if len(final_dict) < len(temp_dict):
                                final_dict[Energy] = 0
                            else:
                                continue
                        else:
                            print("empty line")
                            
                    for value in temp_dict:
                        final_dict[value] = final_dict[value] + temp_dict[value] #sums appropriate values #here
                    fileinsides.close()
                if len(compile_dict) < len(final_dict):
                    for k, v in list(final_dict.items()):
                        compile_dict[k] = []
                else:
                    continue     
            for k, v in list(final_dict.items()):
                compile_dict[k].append(v)  #This is the dictionary with energy bins and lists of events corresponding to interaction type

        #Creates the condensed data file & applies formatting
        condensed_file ="{0}/out/Collated_{1}_{2}_events_{3}_{4}.dat".format(Branch, flux, detector_label2, smear_value, weight_value) #this part making new files with only useful info
        FilesToCleanup.append(condensed_file)
        new_f = open(condensed_file,"w")
        if len(arg) == 4:  
            new_f.write("Energy(GeV)            {0}                    {1}                    {2}                {3}".format(arg[0], arg[1], arg[2], arg[3]))
        elif len(arg) == 5:
            new_f.write("Energy(GeV)            {0}                    {1}                    {2}                {3}                    {4}".format(arg[0], arg[1], arg[2], arg[3], arg[4]))
        else:
            new_f.write("Energy(GeV)            {0}                    {1}                    {2}                {3}                    {4}                    {5}".format(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5]))
        new_f.write("\n")
        new_f.write("-"*100)
        new_f.write("\n") 
                
        #converts the dictionary into readable info
        for key in compile_dict:
            first_space_num = 23 - len(str(key))
            first_space = " " * first_space_num
            value = list(compile_dict[key])
            if len(value) != 1:
                i = 0
                while i < len(value): #Formats the lines so that there are sufficient spaces between each value for them to be readable
                    if len(str(value[i])) > 16:
                        separation_space_num = 23-len(str(round(value[i], 16)))
                        separation_space = " " * separation_space_num
                        value[i] = str(round(value[i], 16)) + separation_space
                        i +=1
                    else:
                        separation_space_num = 23-len(str(round(value[i], 15)))
                        separation_space = " " * separation_space_num
                        value[i] = str(round(value[i], 15)) + separation_space
                        i +=1
                if len(arg) == 4:
                    new_f.write("{0}{1}{2}{3}{4}{5}".format(key, first_space, value[0], value[1], value[2], value[3]))
                elif len(arg) == 5:
                    new_f.write("{0}{1}{2}{3}{4}{5}{6}".format(key, first_space, value[0], value[1], value[2], value[3], value[4])) #for normal, next is for spreadsheet counting
                else:
                    new_f.write("{0}{1}{2}{3}{4}{5}{6}{7}".format(key, first_space, value[0], value[1], value[2], value[3], value[4], value[5])) #for normal, next is for spreadsheet counting
            new_f.write('\n')
        new_f.close()

        if (skip_plots is False):
            r = 0

            if detector_input == all: #aka all detectors are being run, not a specific one
                relevant_list = arg
            else:
                relevant_list = sum_categories

            while r < len(relevant_list):              #Graph labels going back up here isn't working 
                new_dict = {}
                for key in compile_dict:
                    new_dict[key] = compile_dict[key][r]
            
                if arg[r] == "ibd":                #Color labels on graph determined by name, given from value passed through flux details
                    name = "Inverse Beta Decay"
                elif arg[r] == "*_e_":
                    name = r'${\nu}_x+e^-$'
                elif arg[r] == "nc_":
                    name = "Neutral Current"
                elif arg[r] == "Pb208_1n":
                    name = r'${}^{208}Pb$' + " 1n"
                elif arg[r] == "Pb208_2n":
                    name = r'${}^{208}Pb$' + " 2n"
                elif arg[r] == "nue_O16":
                    name = r'${\nu}_e$' + " " + r'${}^{16}O$'
                elif arg[r] == "nuebar_O16":
                    name = r'$\bar{{\nu}_e}$' + " " + r'$^{16}O$'
                elif arg[r] == "nue_C12":
                    name = r'${\nu}_e$' + " " + r'${}^{12}C$'
                elif arg[r] == "nue_C13":
                    name = r'${\nu}_e$' + " " + r'${}^{13}C$'
                elif arg[r] == "nuebar_C12":
                    name = r'$\bar{{\nu}_e}$' + " " + r'$^{12}C$'
                elif arg[r] == "nue_Ar40":
                    name = r'${\nu}_e$' + " " + r'${}^{40}Ar$'
                elif arg[r] == "nuebar_Ar40":
                    name = r'$\bar{{\nu}_e}$' + " " + r'${}^{40}Ar$'
                else:
                    print("Invalid Input")

                if sum(new_dict.values()) != 0:         #plotting each color of bars individually
                
                    #Creates the plots
                    all_values = list(new_dict.values())
                    max_val = max(all_values) #ensures only values greater than 0.1 show up on the plot and legend
                    if max_val > 0.1:
                        plt.plot(list(new_dict.keys()), list(new_dict.values()), linewidth = 1, drawstyle = 'steps', color = colors[r], label = name)
                r += 1

            axes = plt.gca()
            axes.set_xlim([None,0.10])
            axes.set_ylim([0.10, None]) 
            axes.set_yscale('log')
            plt.legend(bbox_to_anchor = (0.5, 0.5, 0.5, 0.5), loc='best', borderaxespad=0)        #formats complete graph
            plt.ylabel (y_label)
            plt.xlabel (x_label)
            plt.title (str(flux).capitalize() + " " + str(detector_label).capitalize() + " " + str(weight_value).capitalize() + " " + str(smear_title).capitalize() + " Events")
            plt.savefig("{0}/out/{1}_{2}_{3}_{4}_log_plot.png".format(Branch,flux, detector, smear_value, weight_value), dpi = 300, bbox_inches = 'tight')
            FilesToCleanup.append("{0}/out/{1}_{2}_{3}_{4}_log_plot.png".format(Branch,flux, detector, smear_value, weight_value))
            plt.clf()
        
    #flux details                 # detector, smeared/unsmeared, *arg reaction placeholders
    
    #The driver, runs through the addition function for every detector
    #The values at the end correspond to the different categories that each file is being summed into
    if detector_input == all:
        for single_flux in FluxFileNameStems:
            for smearval in smearvals:
                add_funct (str(single_flux), "wc100kt15prct", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16") #everything after smearval corresponds to a *arg value
                add_funct (str(single_flux), "wc100kt30prct", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct (str(single_flux), "wc100kt30prct_he", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct (str(single_flux), "halo1", smearval, "nc_", "*_e_", "Pb208_1n", "Pb208_2n") 
                add_funct (str(single_flux), "halo2", smearval, "nc_", "*_e_", "Pb208_1n", "Pb208_2n")
                add_funct (str(single_flux), "ar40kt_eve", smearval,"nc_", "*_e_", "nue_Ar40", "nuebar_Ar40")
                add_funct (str(single_flux), "ar40kt_he", smearval, "nc_", "*_e_", "nue_Ar40", "nuebar_Ar40")
                add_funct (str(single_flux), "icecube", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct (str(single_flux), "hyperk30prct", smearval, "nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16")
                add_funct (str(single_flux), "scint20kt", smearval, "nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12","nue_C13")
                add_funct (str(single_flux), "novaND", smearval, "nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12")
                add_funct (str(single_flux), "novaFD", smearval, "nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12")

            flux_position = FluxFileNameStems.index(single_flux) + 1
            total_fluxes = len(FluxFileNameStems)
            plot_percent_calc = (flux_position/total_fluxes)*100
            plot_percentage_done = str(round(plot_percent_calc,2))
            if (verbose): print('\n'*3)
            if (verbose): print("Plots are " + plot_percentage_done + "% completed.")
            if (verbose): print('\n'*3)
    else: #if just one detector is inputted in SNEWPY.py
        for single_flux in FluxFileNameStems:
            for smearval in smearvals:
                if (verbose):  print("Running inputted detector")
                if detector_input in ("wc100kt15prct", "wc100kt30prct", "wc100kt30prct_he", "icecube", "hyperk30prct", "km3net"):
                    sum_categories = ["nc_", "*_e_", "ibd", "nue_O16", "nuebar_O16"]
                elif detector_input in ("halo1", "halo2"):
                    sum_categories = ["nc_", "*_e_", "Pb208_1n", "Pb208_2n"]
                elif detector_input in ("ar40kt", "ar40kt_he"):
                    sum_categories = ["nc_", "*_e_", "nue_Ar40", "nuebar_Ar40"]
                elif detector_input in ("scint20kt"):
                    sum_categories = ["nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12", "nue_C13"]
                elif detector_input in ("novaND", "novaFD"):
                    sum_categories = ["nc_", "*_e_", "ibd", "nue_C12", "nuebar_C12"]
                else:
                    print("Unanticipated value for detector input")
                argList1 = [str(single_flux), detector_input, smearval]
                argList1.extend(sum_categories)
                add_funct(*argList1)

                
                flux_position = FluxFileNameStems.index(single_flux) + 1
                total_fluxes = len(FluxFileNameStems)
                plot_percent_calc = (flux_position/total_fluxes)*100
                plot_percentage_done = str(round(plot_percent_calc,2))
            if (verbose): print('\n'*3)
            if (verbose): print("Plots are " + plot_percentage_done + "% completed.")
            if (verbose): print('\n'*3)


    #Now create tarball output
    #Makes a tarfile with the condensed data files and plots
    tar=tarfile.open(Model_Path + "/" + outputnamestem + "_SNOprocessed.tar.gz", "w:gz")
    for file in os.listdir(Branch + "/out"):
        if "Collated" in str(file):
            tar.add(Branch + "/out/" + file,arcname=outputnamestem+'_SNOprocessed/'+file)
        elif ".png" in str(file):
            tar.add(Branch + "/out/" + file,arcname=outputnamestem+'_SNOprocessed/'+file)
        else:
            continue
    tar.close()

    if (return_tables is True):
        returned_tables = {}
        for file in os.listdir(Branch + "/out"):
            if "Collated" in str(file):
                returned_tables[file] = {}
                fstream = open(Branch + "/out/"+file, 'r')
                returned_tables[file]['header'] = fstream.readline()
                fstream.close()
                returned_tables[file]['data'] = np.loadtxt(Branch + "/out/" + file,skiprows=2,unpack=True)

    #removes all snowglobes output files, collated files, and .png's made for this snewpy run 
    if (remove_generated_files==True):
        for file in FilesToCleanup:
            os.remove(file)
                
    #Removes all the fluxfiles unzipped from the tarfile
    for file in FluxFileNameStems:
        os.remove(Branch + "/fluxes/" + file + extension)
    try:
        os.remove(Branch + "/fluxes/parameterinfo")
    except OSError:
        print("")

    if (return_tables is True):
        return returned_tables
