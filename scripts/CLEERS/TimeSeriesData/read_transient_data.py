'''
    \file read_transient_data.py
    \brief Script to read in all CLEERS data of a particular folder
    \details Python script to read in CLEERS transient data for
                for a particular folder or folders. This script works
                in tandem with the transient_data.py script which
                contains the necessary objects for storing and operating
                on a set of time series data. What this script does is
                create new objects and methods for dealing with folders
                filled with similar transient data and performing a series
                of like actions on all that data and creating output files
                in sub-folders with the newly processed data.
    \author Austin Ladshaw
    \date 02/24/2020
    \copyright This software was designed and built at the Oak Ridge National
                    Laboratory (ORNL) National Transportation Research Center
                    (NTRC) by Austin Ladshaw for research in the catalytic
                    reduction of NOx. Copyright (c) 2020, all rights reserved.
'''

#Generally speaking, we do not know whether or not there is bypass data in the
#   folders that may need to be paired, so we import all objects from transient_data
from transient_data import TransientData, PairedTransientData
import os, sys, getopt

## ------- Start: TransientDataFolder Object -------- ##
class TransientDataFolder(object):
    def __init__(self,folder):
        if os.path.isdir(folder) == False:
            print("Error! Given argument is not a folder!")
            return
        self.folder_name = folder
        self.file_names = []                    #List of file names for non-bypass files
        self.bypass_names = []                  #List of file names for bypass files
        self.file_pairs = {}                    #Map of paired files: key ==> file_base_name --> (bypass_file, result_file)
        self.unpaired_data = {}                 #Map of TransientData objects by file_name
        self.paired_data = {}                   #Map of PairedTransientData objects by file_name
        self.has_unpaired = False
        self.has_paired = False
        self.unread = []                        #List to correlate with all file_names to determine if a file has been read or not
        #Iterate through the files in the folder and store specific file names
        for file in os.listdir(self.folder_name):
            #Check to make sure the files have the correct extension
            if len(file.split(".")) < 2:
                print("Error! No file extension was given!")
                print("\t"+file+" has been skipped...")
            else:
                if file.split(".")[1] != "dat":
                    print("Error! Unexpected file extension or object is not a file!")
                    print("\t"+file+" has been skipped...")
                else:
                    if "-bp." in file:
                        self.bypass_names.append(file)
                        self.file_pairs[file.split("-bp.")[0]] = []
                    else:
                        self.file_names.append(file)
                        self.unread.append(True)

        #Determine what to do when no bypass files are provided
        if len(self.bypass_names) == 0:
            print("Warning! No bypass data provided...")
            print("\tPutting all data into unpaired_data structure as TransientData objects...")
            self.has_unpaired = True
            for file in self.file_names:
                print("\nReading the following...")
                print("\t"+file)
                self.unpaired_data[file] = TransientData(self.folder_name+"/"+file)
                self.unpaird_data[file].compressColumns()
        #Check for file_names that should correspond to bypass_names
        else:
            self.has_paired = True
            for base in self.file_pairs:
                i=0
                j=0
                for file in self.file_names:
                    if base in file:
                        self.file_pairs[base].append( (base+"-bp.dat",file) )
                        print("\nReading and pairing the following...")
                        print("\t"+self.file_pairs[base][i][0])
                        print("\t"+self.file_pairs[base][i][1])
                        self.paired_data[file] = PairedTransientData(self.folder_name+"/"+self.file_pairs[base][i][0],self.folder_name+"/"+self.file_pairs[base][i][1])
                        self.paired_data[file].compressColumns()
                        self.paired_data[file].alignData()
                        self.unread[j] = False
                        i+=1
                    j+=1

        #Check the unread list for any files that were unread/unpaired and read them as unpaired files
        j=0
        for check in self.unread:
            if check == True:
                self.unpaired_data[self.file_names[j]] = TransientData(self.folder_name+"/"+self.file_names[j])
            j+=1

    def __str__(self):
        message =  "\n ---- Folder: " + str(self.folder_name) + " ---- "
        message += "\n    Number of result files: " + str(len(self.file_names))
        message += "\n    Number of bypass files: " + str(len(self.bypass_names))
        if self.has_paired == True:
            message += "\n ---- Paired Data ---- "
            for base in self.file_pairs:
                message += "\n\t" + self.file_pairs[base][0][0] + " pairs with ..."
                for pair in self.file_pairs[base]:
                    message += "\n\t\t" + str(pair[1])
                message += "\n"
        if self.has_unpaired == True:
            message += "\n ---- Unpaired Data ---- "
            for obj in self.unpaired_data:
                message += "\n\t" + obj
        return message

## ------ Testing ------ ##
test01 = TransientDataFolder("BASFCuSSZ13-700C4h-NH3storage")
print(test01)
