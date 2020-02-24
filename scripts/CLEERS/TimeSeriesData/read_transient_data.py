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
'''
    This object will iteratively read all data files in a given folder and store that
    information into TransientData or PairedTransientData objects depending on the
    typical file name flags used to distinguish between bypass runs and actual data
    runs. The data in each file of the folder can have a bypass pairing or can be
    solitary. HOWEVER, all data files should have all the same common column names.

    NOTE: ALL data files in the folder should have the same column name formatting!
            This is so that we can process entire folders of data en masse by performing
            the same set of actions across all data sets to a particular column.
            IF DATA HAS DIFFERENT COLUMN NAMES, IT SHOULD BE IN DIFFERENT FOLDER!

    For instance: If we want to perform a data reduction analysis and extract only
                    a specific sub-set of columns from all the files, then those
                    column names need to be all the same. That way we can use a single
                    function call/operation to perform the analysis on all data.
'''
class TransientDataFolder(object):
    #Initialize the object by reading in all readable data
    def __init__(self,folder,addNoise = True):
        if os.path.isdir(folder) == False:
            print("Error! Given argument is not a folder!")
            return
        self.folder_name = folder
        self.file_names = []                    #List of file names for non-bypass files
        self.bypass_names = []                  #List of file names for bypass files
        self.file_pairs = {}                    #Map of paired files: key ==> file_base_name --> (bypass_file, result_file)
        self.unpaired_data = {}                 #Map of TransientData objects by file_name
        self.paired_data = {}                   #Map of PairedTransientData objects by file_name
        self.has_unpaired = False               #Flag to denote whether or not folder contains unpaired data
        self.has_paired = False                 #Flag to denote whether or not folder contains paired data (note: CAN have both)
        self.unread = []                        #List to correlate with all file_names to determine if a file has been read or not
        self.was_compressed = False             #Flag to denote whether or not user has requested row compression of data sets
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
                        self.paired_data[file].alignData(addNoise)
                        self.unread[j] = False
                        i+=1
                    j+=1

        #Check the unread list for any files that were unread/unpaired and read them as unpaired files
        j=0
        for check in self.unread:
            if check == True:
                self.unpaired_data[self.file_names[j]] = TransientData(self.folder_name+"/"+self.file_names[j])
            j+=1

    #Print object attributes to the console
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

    #Display names of all columns for all data sets
    #       NOTE: ASSUMES ALL DATA SETS HAVE SAME COMMON COLUMN NAMES
    def displayColumnNames(self):
        if self.has_paired == True:
            print("\n ---- Paired Data Column Names ---- \n")
            self.paired_data[list(self.paired_data.keys())[-1]].displayColumnNames()
        if self.has_unpaired == True:
            print("\n ---- Unpaired Data Column Names ---- \n")
            self.unpaired_data[list(self.unpaired_data.keys())[-1]].displayColumnNames()

    #Delete specific columns from all data sets that we do not need
    def deleteColumns(self, column_list):
        #For unpaired data, just call the corresponding sub-class method
        for item in self.unpaired_data:
            self.unpaired_data[item].deleteColumns(column_list)
        if self.has_paired == True  and self.was_compressed == True:
            print("Warning! Deleting columns after running compressRows() may cause errors in paired data!")
            print("\tYou should run deleteColumns() before running compressRows()...")
        #For the paired data, we need to expand the list to include the [bypass] suffixed info in results_trans_obj
        extended_column_list = []
        if type(column_list) is list:
            for name in column_list:
                extended_column_list.append(name)
                if "[bypass]" not in name:
                    extended_column_list.append(name+"[bypass]")
        else:
            extended_column_list.append(column_list)
            if "[bypass]" not in column_list:
                extended_column_list.append(column_list+"[bypass]")
        for item in self.paired_data:
            self.paired_data[item].deleteColumns(extended_column_list)

    #Retain on specific columns from all data sets
    def retainOnlyColumns(self, column_list):
        #For unpaired data, just call the corresponding sub-class method
        for item in self.unpaired_data:
            self.unpaired_data[item].retainOnlyColumns(column_list)
        if self.has_paired == True  and self.was_compressed == True:
            print("Warning! Deleting columns after running compressRows() may cause errors in paired data!")
            print("\tYou should run deleteColumns() before running compressRows()...")
        #For the paired data, we need to expand the list to include the [bypass] suffixed info in results_trans_obj
        extended_column_list = []
        if type(column_list) is list:
            for name in column_list:
                extended_column_list.append(name)
                if "[bypass]" not in name:
                    if name != self.paired_data[list(self.paired_data.keys())[-1]].result_trans_obj.time_key:
                        extended_column_list.append(name+"[bypass]")
        else:
            extended_column_list.append(column_list)
            if "[bypass]" not in column_list:
                if column_list != self.paired_data[list(self.paired_data.keys())[-1]].result_trans_obj.time_key:
                    extended_column_list.append(column_list+"[bypass]")
        for item in self.paired_data:
            #NOTE: After data has been aligned, we only need to operate on result_trans_obj
            self.paired_data[item].result_trans_obj.retainOnlyColumns(extended_column_list)
            self.paired_data[item].bypass_trans_obj.retainOnlyColumns(column_list) #May be unnecessary

    #Function to return an instance of the TransientData or PairedTransientData object given a file name
    #   NOTE: This function returns different data types depending on the argument it gets
    def grabDataObj(self,file):
        if file in self.paired_data.keys():
            return self.paired_data[file]
        elif file in self.unpaired_data.keys():
            return self.unpaired_data[file]
        else:
            print("Error! No such file was read by this object!")
            return

    #Function to compress all rows of data for each data object according to size of rows
    #       NOTE: This function should be called before printing to a file, but after everything else
    #
    #       This function accepts 2 optional arguments: (i) num_rows_target and (ii) max_factor
    #               num_rows_target is used to represent the number of rows of data we
    #                   want to compress the data to. For instance, num_rows_target = 1000
    #                   means that after compression, we want the data to fit within about
    #                   1000 rows of data.
    #               max_factor puts a cap on the maximum level of compression we will allow,
    #                   regardless of the num_rows_target. For intance, a max_factor of 10
    #                   means that the number of data rows will be reduced by a factor of 10
    #                   at most, and no more than that.
    def compressAllRows(self, num_rows_target = 1000, max_factor = 10):
        #Check the number of rows for each set of data and determine an appropriate compression factor
        for file in self.unpaired_data:
            num_rows = len(self.unpaired_data[file].data_map[self.unpaired_data[file].time_key])
            factor = int(num_rows/num_rows_target)
            if factor > max_factor:
                factor = max_factor
            self.unpaired_data[file].compressRows(factor)
        for file in self.paired_data:
            num_rows = len(self.paired_data[file].result_trans_obj.data_map[self.paired_data[file].result_trans_obj.time_key])
            factor = int(num_rows/num_rows_target)
            if factor > max_factor:
                factor = max_factor
            self.paired_data[file].compressRows(factor)
        self.was_compressed = True

    #Function to perform a math operation to all like columns of data in all data objects
    def mathOperations(self, column_name, operator, value_or_column, append_new = False, append_name = ""):
        return

    #Function to compute a mass retention integral for all columns of the given name
    #   NOTE: This function only works with PAIRED data
    def calculatePairedRetentionIntegrals(self, column_name, normalized = False, conv_factor = 1):
        if self.was_compressed == True:
            print("Warning! Computing integrals after data compression may cause inaccuracies in results!")

    #Function to compute a mass retention integral for all columns of the given name
    #   NOTE: This function only works with UNPAIRED data
    def calculateUnpairedRetentionIntegrals(self, column_name, normalized = False, conv_factor = 1):
        if self.was_compressed == True:
            print("Warning! Computing integrals after data compression may cause inaccuracies in results!")

    #Function to print all results (paired and unpaired) to a series of output files in a sub-directory
    def printAlltoFile(self, subdir = ""):
        return

## ------ Testing ------ ##
test01 = TransientDataFolder("BASFCuSSZ13-700C4h-NH3storage")
#test01.deleteColumns(["Time",'CO2% (20)'])
test01.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)', 'TC bot sample in (C)', 'TC bot sample mid 1 (C)', 'TC bot sample mid 2 (C)', 'TC bot sample out 1 (C)', 'TC bot sample out 2 (C)', 'P bottom in (bar)', 'P bottom out (bar)'])
#test01.displayColumnNames()
#print(test01)
#print(test01.grabDataObj("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-200C.dat"))

test01.compressAllRows()
