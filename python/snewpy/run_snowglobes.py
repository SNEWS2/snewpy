#!/usr/bin/env python
from __future__ import unicode_literals
import os
import tarfile
import zipfile
import re
import matplotlib as mpl
mpl.use('Agg')
import shutil 
from subprocess import call

# Tomer K. Goldhagen
# Translated from perl to python June 2019
# Anne Graf
# Made it Better Summer 2020
# Takes in input flux files and configures and runs supernova (which outputs calculated rates)

#for tar file in folder, call go & replace go in master
           
def go(SNOwGLoBESdir, Models_Path, Tarball, detector_input = all, verbose=False):  #function for entire file  
    #Extracts data from tarfile and sets up lists of paths and fluxfilenames for later use

    if tarfile.is_tarfile(Models_Path+"/"+Tarball): #extracts tarfile
        if (verbose): print("Valid Tar file")
        tar = tarfile.open(Models_Path+"/"+Tarball)
        TarballFileNames = tar.getnames()
        tar.extractall(path = SNOwGLoBESdir + "/fluxes")
        tar.close()
        if (verbose): print("Tar file data extracted")
    elif zipfile.is_zipfile(Models_Path+"/"+Tarball): #extracts zipfile
        if (verbose): print("Valid Zip file")
        zip = zipfile.ZipFile(Models_Path+"/"+Tarball)
        TarballFileNames = zip.namelist()
        for fileName in TarballFileNames:
            if str(fileName)[0] == "." or str(fileName)[0] == "_":
                continue
            else:        
                zip.extract(fileName, path = SNOwGLoBESdir + "/fluxes")
                if (verbose): print("zipfile extracted")
        zip.close()
        if (verbose): print("Zip file data extracted")
    else:
        print("Invalid Tarfile") #doesn't accept anything other than tar and zip files
       
    #Creates a new out folder if one does not exist
        
    out_folder = "{0}/out".format(SNOwGLoBESdir)     
    try:
        os.mkdir(out_folder)
    except OSError as e:
        if (verbose): print("out folder already exists")
    else:
        if (verbose): print("out folder made")
        
    #Creates directories for each flux file within the out folder

    FluxFileNameStems = []
    #extension = '.dat'
    for IndividualFile in TarballFileNames:
        extension_beginning = str(IndividualFile).rfind(".")
        if extension_beginning > -1:
            extension = str(IndividualFile[extension_beginning:])
            ExtensionFound = str(IndividualFile).rfind(extension)           
            FluxFileNameStems.append( str(IndividualFile)[0:ExtensionFound] ) # gets rid of extension at the end of flux files            
        else:
            continue

    ##################################################################################################
    #The main function that reformats data from a bunch of different files in order to run supernova.c

    def format_globes_for_supernova(FluxNameStem, ChannelName, ExpConfigName, NoWeight = 1, extension = extension): #this function runs supernova
        ExeName = "bin/supernova" #Supernova
        ChannelFileName = "{0}/channels/channels_{1}.dat".format(SNOwGLoBESdir, ChannelName)

        GlobesFileName = "{0}/supernova.glb".format(SNOwGLoBESdir)
        GLOBESFILE = open(GlobesFileName, 'w')
        GLOBESFILE.truncate(0)

        # Write PREAMBLE to GLOBESFILE

        PREAMBLE = open("{0}/glb/preamble.glb".format(SNOwGLoBESdir))

        for line in PREAMBLE:

            GLOBESFILE.write(line)

        PREAMBLE.close()

        # Checks FluxFileNamePath existence

        FluxFileNamePath = "fluxes/" + FluxNameStem + extension
        if not os.path.isfile(SNOwGLoBESdir+"/"+FluxFileNamePath): 
            print(("Flux file name {0} not found".format(FluxFileNamePath)))
            return "Incomplete"

        # Write modified FLUX to GLOBESFILE; replace the flux file name with the input argument

        FLUX = open("{0}/glb/flux.glb".format(SNOwGLoBESdir))

        for line in FLUX:
            if "flux_file" in line: line = "        @flux_file=  \"{0}\"\n".format(FluxFileNamePath)
            GLOBESFILE.write(line)

        FLUX.close()

        # Now go through the channels and put in the relevant lines

        # First, smearing for each channel

        # Checks channel existence
        if not os.path.isfile(ChannelFileName):
            print(("Channel file name {0} not found".format(ChannelFileName)))
            return "Incomplete"

        # Grabs channel name from each line in CHANFILE and write it (modified) to GLOBESFILE
        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')
            ChanName = TempArray[0];

            GLOBESFILE.write("include \"{0}/smear/smear_{1}_{2}.dat\"\n".format(SNOwGLoBESdir, ChanName, ExpConfigName))

        CHANFILE.close()
        
        # Checks DetectorFileNamePath existence

        DetectorFileName = "detector_configurations.dat"
        DetectorFileNamePath = "{0}/{1}".format(SNOwGLoBESdir, DetectorFileName)

        if not os.path.isfile(DetectorFileNamePath):
            print(("Detector file name {0} not found".format(DetectorFileName)))
            return "Incomplete"

        # Defines Masses and NormFact (hashes, A.K.A. dictionaries) for later use

        DETFILENAME = open(DetectorFileNamePath)
        Masses = {}
        NormFact = {}    
        for line in DETFILENAME:

            if line[-1] == "\n": line = line[:-1]

            if "#" in line: continue

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')
            DetName = TempArray[0]
            Masses[DetName] = TempArray[1]
            NormFact[DetName] = TempArray[2]

            if (DetName == "" or Masses == "" or NormFact == ""): continue
        DETFILENAME.close()

        # Mass and target normalization by species

        # Checks Masses[ExpConfigName] existence

        if not Masses.get(ExpConfigName):
            print("Error: please enter a valid experiment configuration")
            return "Incomplete"

        TargetMass = format(float(Masses[ExpConfigName]) * float(NormFact[ExpConfigName]), "13.6f") 

        # Writes modified DETECTOR to GLOBESFILE; replace the mass with the calculated TargetMass 

        DETECTOR = open("{0}/glb/detector.glb".format(SNOwGLoBESdir))

        for line in DETECTOR:

            if "mass" in line: line = "$target_mass=  {0}".format(TargetMass);
    
            GLOBESFILE.write(line)

        DETECTOR.close()
        # Now the cross-sections.  Note that some of these are repeated even though 
        # it is not necessary (xscns for several flavors can be in the same file).
        # This is just to make a consistent loop over channels.

        GLOBESFILE.write("\n\n /******** Cross-sections *********/\n \n")

        # Writes a modified CHANFILE to GLOBESFILE

        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')

            ChanName = TempArray[0]

            GLOBESFILE.write("cross(#{0})<\n".format(ChanName))
            GLOBESFILE.write("      @cross_file= \"{0}/xscns/xs_{1}.dat\"\n".format(SNOwGLoBESdir, ChanName))
            GLOBESFILE.write(">\n")

        CHANFILE.close()
        # Now, the channel definitions

        GLOBESFILE.write("\n /******** Channels *********/\n \n")

        # Checks ChannelFileName existence

        if not os.path.isfile(ChannelFileName):
            print(("Channel file name {0} not found".format(ChannelFileName)))
            return "Incomplete"

        # Writes a modified CHANFILE and EFF_FILE to GLOBESFILE
        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:
    
            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ') 

            ChanName = TempArray[0]
            CpState = TempArray[2]
            InFlav = TempArray[3]
    
            GLOBESFILE.write("channel(#{0}_signal)<\n".format(ChanName))
            GLOBESFILE.write("      @channel= #supernova_flux:  {0}:    {1}:     {2}:    #{3}:    #{4}_smear\n".format(CpState, InFlav, InFlav,                                                                         ChanName, ChanName))
    
            # Get the post-smearing efficiencies by channel
    
            EffFile = "{0}/effic/effic_{1}_{2}.dat".format(SNOwGLoBESdir, ChanName, ExpConfigName)
            try:
                EFF_FILE = open(EffFile)
            except IOError:
                GLOBESFILE.write("\n>\n\n")
                continue
            else:
                for line in EFF_FILE:
                    GLOBESFILE.write("       @post_smearing_efficiencies = {0}".format(line))

                GLOBESFILE.write("\n>\n\n")

                EFF_FILE.close()
        CHANFILE.close()

        # End matter

        # Writes POSTAMBLE to GLOBESFILE

        POSTAMBLE = open("{0}/glb/postamble.glb".format(SNOwGLoBESdir))

        for line in POSTAMBLE:
    
            GLOBESFILE.write(line)

        POSTAMBLE.close()
        GLOBESFILE.close()

        # Runs supernova
        os.chdir(SNOwGLoBESdir)
        ComString = "{0} {1} {2} {3}".format(ExeName, FluxNameStem, ChannelFileName, ExpConfigName)
        call(ComString, shell = True)

        # Calls apply_weights to format the output of supernova

        if NoWeight != 1:
            if (verbose): print("Applying channel weighting factors to output")

            apply_weights("unsmeared", FluxNameStem, ChannelFileName, ExpConfigName)
            apply_weights("smeared", FluxNameStem, ChannelFileName, ExpConfigName)
    
        return "Complete"
        
    #######################################################################
    #This function applies the weights to the weightedfilename via tempfile
    
    def apply_weights(smearing, FluxNameStem, ChannelFileName, ExpConfigName):   
        # Opens CHANFILE in order for the code to be able to find all of the output files (from supernova)

        #formats for naming of the files
        if smearing == "unsmeared":
            smear_input = "" # snowglobes1.2 does not insert "_unsmeared" into its output filenames
        elif smearing == "smeared":
            smear_input = "_smeared"
        else:
            print("Invalid smearing input")
        
        CHANFILE = open(ChannelFileName)

        for line in CHANFILE:

            clean_line = re.sub('\s*', '', line)
            clean_line = re.sub('\s+', ' ', line)

            TempArray = clean_line.split(' ')

            ChanName = TempArray[0]
            NumTargetFactordraft = TempArray[4]
            if NumTargetFactordraft[-1] == ".":
                NumTargetFactor = float(NumTargetFactordraft[:-1]) #the weighting factor
            else:
                NumTargetFactor = float(NumTargetFactordraft)


            # Writes modified UNWEIGHT to WEIGHTED

            UnweightFileName = "{0}/out/{1}_{2}_{3}_events{4}_unweighted.dat".format(SNOwGLoBESdir, FluxNameStem, ChanName, ExpConfigName, smear_input)
            WeightedFileName = "{0}/out/{1}_{2}_{3}_events{4}_weighted.dat".format(SNOwGLoBESdir, FluxNameStem, ChanName, ExpConfigName, smear_input)

            UNWEIGHT = open(UnweightFileName)
            
            #This applies the weighting factor, and reformates a few lines for consistency
            with open(UnweightFileName,"r") as f:
                WEIGHT = open(WeightedFileName, "w")
                for line in f:
            
                    clean_line = re.sub('\s*', '', line)
                    clean_line = re.sub('\s+', ' ', line)
                
                    if "---" in line:
                        WEIGHT.write(line + "\n")
                    elif "Total:" in line:
                        continue
                    else:
                        TempArray2 = clean_line.split(' ')
                        Enbin = TempArray2[0]
                        Evrate = TempArray2[1] 
                        if Enbin != "": #writing the new weighted lines to a file
                            newline = "{0}            {1} ".format(Enbin, float(Evrate) * float(NumTargetFactor))
                            WEIGHT.write(newline + "\n")
                        elif Enbin == "" and Evrate == "0.05" or Evrate == "0.051":
                            newline = "{0}            {1} ".format(Evrate, float(TempArray2[2]) * float(NumTargetFactor))
                            WEIGHT.write(newline + "\n")
                        else:
                            continue
                WEIGHT.close()      
        CHANFILE.close()

    fluxes_list = [] #for a list of IndividualFluxFileName files

    #THE DRIVER
    #This is the place where you can comment out any detectors you don't want to run
    #This runs the entire module, for each detector configuration

    for IndividualFluxFile in FluxFileNameStems:
        if (verbose): print(IndividualFluxFile)
        position = FluxFileNameStems.index(IndividualFluxFile) + 1
        total = len(FluxFileNameStems)
        percent_calc = (position/total)*100
        percentage_done = str(round(percent_calc,2))
        
        if detector_input == all:
            if (verbose): print("Running all detectors")
            if format_globes_for_supernova(IndividualFluxFile, "water", "icecube", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||icecube"'.format(IndividualFluxFile))
                if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "water", "wc100kt30prct", "weight") == "Complete": # each of these calls the subroutine with varying
                if (verbose): os.system('echo "Finished {0}||wc100kt30prct"'.format(IndividualFluxFile))          # settings and echoes a completion confirmation
                if (verbose): print('')
            #if format_globes_for_supernova(IndividualFluxFile, "water", "wc100kt30prct_he", "weight") == "Complete": 
                #if (verbose): os.system('echo "Finished {0}||wc100kt30prct_he"'.format(IndividualFluxFile)) 
                #if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "water", "wc100kt15prct", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||wc100kt15prct"'.format(IndividualFluxFile))
                if (verbose): print('')
            #if format_globes_for_supernova(IndividualFluxFile, "water", "hyperk30prct", "weight") == "Complete": #not working
                #if (verbose): os.system('echo "Finished {0}||hyperk30prct"'.format(IndividualFluxFile))
                #if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "argon", "ar40kt", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||ar40kt"'.format(IndividualFluxFile))
                if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "lead", "halo1", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||halo1"'.format(IndividualFluxFile))
                if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "lead", "halo2", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||halo2"'.format(IndividualFluxFile))
                if (verbose): print('')
            #if format_globes_for_supernova(IndividualFluxFile, "argon", "ar40kt_he", "weight") == "Complete":
                #if (verbose): os.system('echo "Finished {0}||ar40kt_he"'.format(IndividualFluxFile))
                #if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "scint", "scint20kt", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||scint20kt"'.format(IndividualFluxFile))
                if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "nova_soup", "novaND", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||novaND"'.format(IndividualFluxFile))
                if (verbose): print('')
            if format_globes_for_supernova(IndividualFluxFile, "nova_soup", "novaFD", "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||novaFD"'.format(IndividualFluxFile))
                if (verbose): print('')
            fluxes_list.append(IndividualFluxFile)
            if (verbose): print('\n'*3)
            if (verbose): print("Calculations are " + percentage_done + "% completed.")
            if (verbose): print('\n'*3)
        else: #This is called if you choose to input a single detector in SNEWPY.py, and just run that one
            if (verbose): print("Running inputted detector")
            if (verbose): print(detector_input)
            if detector_input in ("icecube", "wc100kt30prct", "wc100kt30prct_he", "wc100kt15prct", "hyperk30prct", "km3net"):
                detector_material = "water"
            elif detector_input in  ("ar40kt_he", "ar40kt"):
                detector_material = "argon"
            elif detector_input in ("novaND", "novaFD"):
                detector_material = "nova_soup"
            elif detector_input in ("scint20kt"):
                detector_material = "scint"
            elif detector_input in ("halo1", "halo2"):
                detector_material = "lead"
            else:
                print("Unanticipated value for detector input")
                
            detector_output = detector_input
            if (verbose): print(detector_material)
            if format_globes_for_supernova(IndividualFluxFile, detector_material, detector_input, "weight") == "Complete":
                if (verbose): os.system('echo "Finished {0}||{1}"'.format(IndividualFluxFile, detector_output))
                if (verbose): print('')
            fluxes_list.append(IndividualFluxFile)
            if (verbose): print('\n'*3)
            if (int(percent_calc)%10==0): print("Calculations are " + percentage_done + "% completed. (",position," of ",total,")")
            if (verbose): print('\n'*3)

