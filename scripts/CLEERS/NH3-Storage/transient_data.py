'''
    \file transient_data.py
    \brief Read TransientData from CLEERS team
    \details Python script to read in CLEERS transient data for
                NH3 storage on Cu-SSZ-13. This script will store the
                orginal data as is and provide other functions to
                redistribute, print, or parse that data as needed.
    \author Austin Ladshaw
    \date 02/07/2020
    \copyright This software was designed and built at the Oak Ridge National
                    Laboratory (ORNL) National Transportation Research Center
                    (NTRC) by Austin Ladshaw for research in the catalytic
                    reduction of NOx. Copyright (c) 2020, all rights reserved.
'''

    #NOTE: The CLEERS data files are very, very large, so I am saving them as *.dat files.
    #       The reasoning behind this is so that I can direct 'git' to ignore files that
    #       end with a *.dat file extension. This prevents the repository from becoming
    #       bloated. The *.dat files behave exactly like regular text files.

import math
import os
from statistics import mean, stdev
import random

## ---------------- Begin: Definition of TransientData object ------------
class TransientData(object):
    # Initialize data object by passing the current file to it
    # Each key in the data_map represents a column label
    #       Each key maps to a data list
    # The input file should have a specific convention for naming a file
    #       e.g., 20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-150C.dat
    #
    #       Each important piece of information is split by a "-" character. We then use this to parse
    #   the file name to obtain particular information. However, that file name convention is not necessarily
    #   consistent for all files. So we are limited in what can be interpreted from the names. The most important
    #   information is as follows...
    #
    #       item[2] = Name of the catalyst material
    #       item[3] = Aging condition
    #       item[5] = flow rate information (in kilo-volumes per hour)
    #       item[-1] = either a temperature or "bp"
    #                   temperature is irrelevant, but bp indicates that this is inlet information
    #                   (Also note, item[-1] will carry the file extension with it)
    def __init__(self, file):
        #Parse the file name to gain specific information
        self.input_file_name = file
        file_name_info = file.split("-")
        self.material_name = file_name_info[2]
        if file_name_info[3] == "700C4h":
            self.aging_condition = "De-greened"
        else:
            self.aging_condition = file_name_info[3]
        try:
            self.flow_rate = float(file_name_info[5].split("k")[0])*1000
            self.have_flow_rate = True
        except:
            self.flow_rate = file_name_info[5]
            self.have_flow_rate = False
        if file_name_info[-1].split(".")[0] == "bp":
            self.inlet_data = True
            self.isothermal_temp = 25
        else:
            self.inlet_data = False
            try:
                self.isothermal_temp = float(file_name_info[-1].split(".")[0].split("C")[0])
            except:
                self.isothermal_temp = file_name_info[-1].split(".")[0]

        self.data_file = open(file,"r") #Contains data file we are digitizing
        statinfo = os.stat(file)
        self.exp_header = ''            #Contains the first line of the data file
        self.data_map = {}              #Contains a map of all the data by column
        self.num_rows = 0               #Contains the number of rows of data
        self.ordered_key_list = []      #Contains an ordered list of column names
        self.change_time = []           #Contains an ordered list of the times when experimental inputs changed
                                        #   NOTE: this is just used to initialize data in the map, it does not
                                        #           change or update if the map changes
        self.input_change = {}          #Contains a map of the inputs values that correspond to change_time
                                        #   Keys of this map are modifications to the keys of data_map that
                                        #   correspond to output values corresponding to the input values given

        self.time_key = ""              #Contains the name of the time key for the data_map
                                        #   Key must contain "Elapsed Time (unit)" in the data file
                                        #   unit can be any time units, but the first part of the string
                                        #   must alwas exist as "Elapsed Time (..."
        if statinfo.st_size >= 10000000:
            print("Reading large file. Please wait...")
        self.readFile()
        if statinfo.st_size >= 10000000:
            print("Finished!")

    def __str__(self):
        message = "\nFile Name:\t" + self.data_file.name
        message += "\nMaterial:\t" + self.material_name
        message += "\nAging Condition:\t" + self.aging_condition
        message += "\nFlow Rate (hr^-1):\t" + str(self.flow_rate)
        message += "\nInlet Conditions:\t"
        if self.inlet_data == True:
            message += "True"
        else:
            message += "False"
            message += "\nIsothermal Temp (C):\t" + str(self.isothermal_temp)
        message += "\nFile Header:\t" + self.exp_header
        message += "\nNumber of Columns:\t" + str(len(self.data_map))
        message += "\nNumber of Rows:\t" + str(self.num_rows)
        return message + "\n"

    def displayColumnNames(self):
        print(self.data_map.keys())

    def readFile(self):
        i = 0
        first_data_line = False
        for line in self.data_file:
            line_list = line.split('\t')
            #Ignore the first line of the data file and stores as header
            if (i == 0):
                self.exp_header = line

            # This is the first real header of the data
            if (i == 1):
                for item in line_list:
                    #Ignore the ending character
                    if item != '\n':
                        #Find the time_key by parsing and checking item
                        if item.strip().split("(")[0] == "Elapsed Time ":
                            self.time_key = item.strip()
                        #NOTE: The registered keys will be the column names stripped of leading and trailing whitespaces
                        self.data_map[item.strip()] = []
                        self.ordered_key_list.append(item.strip())
                    else:
                        # Force a new column for input conditions
                        self.change_time.append(0)

            if i == 2:
                first_data_line = True
            else:
                first_data_line = False
            if (i > 1): #or greater
                n = 0
                for item in line_list:
                    #Check item to make sure that we do not get another column key
                    # If we do, then do not append data to map, instead record
                    #   the data_map[time_key] value from prior as the point
                    #   when input conditions changed
                    #NOTE: Check the stripped items for keys
                    if item.strip() in self.data_map.keys():
                        if (n == 0):
                            self.change_time.append(float(self.data_map[self.time_key][-1]))
                    # If we don't, then record the data into the map
                    else:
                        #Ignore ending character
                        if item != '\n':
                            #Attempt to store data as a number
                            try:
                                self.data_map[self.ordered_key_list[n]].append(float(item))
                            #Store data as string if it is not a number
                            except:
                                if first_data_line == True:
                                    self.data_map[self.ordered_key_list[n]].append(item)
                                else:
                                    if len(self.data_map[self.ordered_key_list[n]]) > 0:
                                        if type(self.data_map[self.ordered_key_list[n]][-1]) is not int and type(self.data_map[self.ordered_key_list[n]][-1]) is not float:
                                            self.data_map[self.ordered_key_list[n]].append(item)
                                        #else:
                                            #print("Warnging")
                                    #else:
                                        #print("Warning")
                                #self.data_map[self.ordered_key_list[n]].append(item)
                                #if i == 2:
                                    #self.data_map[self.ordered_key_list[n]].append(item)
                                #else:
                                    #if type(self.data_map[self.ordered_key_list[n]][i-3]) is not int and type(self.data_map[self.ordered_key_list[n]][i-3]) is not float:
                                        #self.data_map[self.ordered_key_list[n]].append(item)
                                    #else:
                                        #print("Warning! Some data lost...")
                    n+=1
            i+=1
        #END of line loop
        self.num_rows = len(self.data_map[self.time_key])
        self.closeFile()

    def closeFile(self):
        self.data_file.close()

    #This function will add a column to the data map given the column name and associated data
    #   NOTE: Appending the column does NOT copy the column into the map. It merely directs
    #       the map to point to the given data_set. If you change the data_set that you pass
    #       to this function, then the data in this object's map will also change.
    def appendColumn(self, column_name, data_set):
        #First, check to make sure the map has been prepared
        if (len(self.data_map) == 0):
            print("Error! File has not been read and stored!")
            return

        #Make sure the data_set is actually a list of data
        if type(data_set) is not list:
            print("Error! The data_set must be a list of data!")
            return

        #Next, check to make sure the column_name is not already in the map
        if column_name in self.data_map.keys():
            if len(self.data_map[column_name]) != len(data_set):
                print("Error! That data column is already in the structure and sizes are mis-matched!")
                return

        #Lastly, check to make sure the length of the data_set matches the length of the time series data
        if len(self.data_map[self.time_key]) != len(data_set):
            print("Error! The data set size does not match the existing data set size!")
            return

        self.data_map[column_name] = data_set

    #This function will extract column sets from the data_map and return a new, reduced map
    #   NOTE: Column list must be a list of valid keys in the data_map
    def extractColumns(self, column_list):
        new_map = {}
        if type(column_list) is list:
            for item in column_list:
                #Check to make sure the item is a key in data_map
                if item in self.data_map.keys():
                    new_map[item] = self.data_map[item]
                else:
                    print("Error! Invalid Key!")
            #End item loop
        else:
            if column_list in self.data_map.keys():
                new_map[column_list] = self.data_map[column_list]
            else:
                print("Error! Invalid Key!")
        return new_map

    #This function will extract a row of data (or set of rows) based on the value of Elapsed time provided
    def extractRows(self, min_time, max_time):
        new_map = {}
        for item in self.data_map:
            new_map[item] = []
        #loop through all time data
        n = 0
        for time in self.data_map[self.time_key]:
            if time > max_time:
                break
            if time >= min_time:
                for item in self.data_map:
                    if len(self.data_map[item]) > n:
                        new_map[item].append(self.data_map[item][n])
            n+=1

        return new_map

    # This function will get a particular data point based on the given elapsed time and column name
    def getDataPoint(self, time_value, column_name):
        point = 0
        if column_name not in self.data_map.keys():
            print("Error! Invalid Key!")
            return 0
        n = 0
        start_time = 0
        start_point = 0
        end_time = 0
        end_point = 0
        if time_value >= self.data_map[self.time_key][-1]:
            return self.data_map[column_name][-1]
        #What should we do if we reach beyond the total time?
        for time in self.data_map[self.time_key]:
            if time > time_value:
                if n == 0:
                    start_time = 0
                    start_point = self.data_map[column_name][n]
                else:
                    start_time = self.data_map[self.time_key][n-1]
                    start_point = self.data_map[column_name][n-1]
                end_time = time
                end_point = self.data_map[column_name][n]
                break
            n+=1
        #Perform linear interpolation between start_point and end_point
        if type(end_point) is not int and type(end_point) is not float:
            point = end_point
        else:
            try:
                point = (end_point - start_point)/(end_time - start_time)*(time_value - start_time) + start_point
            except:
                point = end_point
        return point

    #This function will iterate through all columns to find data that can be compressed or eliminated
    #       For instance,   if a column contains no data, then delete it
    #                       if there are multiple columns that carrier similar info, then combine them
    def compressColumns(self):
        #Iterate through the data map to find information to delete immediately
        list_to_del = []
        for item in self.data_map:
            if len(self.data_map[item]) < len(self.data_map[self.time_key]):
                list_to_del.append(item)
        for item in list_to_del:
            del self.data_map[item]

        #Iterate through the map and find compressible columns
        frac_keys = {}  #Map of a list of columns that can compress
        for item in self.data_map:
            first = item.split("(")
            if len(first) > 1:
                last = first[1].split(")")
            else:
                last = ""
            if first[0] in frac_keys.keys():
                #Value check is used to determine whether or not the duplicated key
                #       is significant or can be neglected. If it is significant, then
                #       the value will be a number. Otherwise, the columns can't be merged.
                try:
                    val = float(last[0])
                    frac_keys[first[0]].append(item)
                except:
                    val = 0
            else:
                frac_keys[first[0]] = []
                frac_keys[first[0]].append(item)
        #Any frac_keys[key] ==> list that has a length > 1 can be compressed
        #       Need to add new key to data_map to hold compressed data, then
        #       delete the old keys to reduce the map size.
        for item in frac_keys:
            if len(frac_keys[item]) > 1:
                val_list = []
                new_name = item +"("
                j=0
                for sub_key in frac_keys[item]:
                    first = sub_key.split("(")
                    last = first[1].split(")")
                    val_list.append(float(last[0]))
                    if j == 0:
                        new_name += str(int(val_list[j]))
                    else:
                        new_name += "," + str(int(val_list[j]))
                    j+=1
                new_name += ")"
                self.data_map[new_name] = []
                #Now loop through all rows and insert proper data into new column
                n = 0
                for value in self.data_map[frac_keys[item][0]]:
                    #value we are on here corresponds to the val_list[0] limit
                    #Loop through all other columns that can compress and
                    #find value to register based on limits in val_list
                    dist = []
                    dist.append(val_list[0] - value)
                    for i in range(1,len(frac_keys[item])):
                        dist.append(val_list[i] - self.data_map[frac_keys[item][i]][n])

                    #The index of the smallest postive dist is the index i for the data to append
                    reg_index = 0
                    old_d = dist[0]
                    for i in range(1,len(frac_keys[item])):
                        if dist[i] > 0 and old_d > 0:
                            if dist[i] < old_d:
                                reg_index = i
                                old_d = dist[i]
                        elif old_d < 0:
                            old_d = dist[i]
                            reg_index = i
                    self.data_map[new_name].append(self.data_map[frac_keys[item][reg_index]][n])
                    n+=1
                #End value loop

                #Now, delete the original columns
                for sub_key in frac_keys[item]:
                    del self.data_map[sub_key]

            #End if


    #This function will delete the given columns from the map
    def deleteColumns(self, column_list):
        #NOTE: column_list is either a list of columns or a single column_name
        if type(column_list) is list:
            #iterate through list and delete any columns that can be deleted
            for item in column_list:
                if item in self.data_map.keys():
                    del self.data_map[item]
                else:
                    print("Error! No such column exists!")
        else:
            #delete only the given column
            if column_list in self.data_map.keys():
                del self.data_map[column_list]
            else:
                print("Error! No such column exists!")

    #This function will delete all columns in the data_map except for the ones specified to retain
    def retainOnlyColumns(self, column_list):
        #First, check to make sure that the columns named in the list are valid
        keep = {}
        if type(column_list) is list:
            for item in column_list:
                if item not in self.data_map.keys():
                    print("Error! Invalid Column Name! No further action taken to delete columns...")
                    return
                else:
                    keep[item] = 0
        else:
            if column_list not in self.data_map.keys():
                print("Error! Invalid Column Name! No further action taken to delete columns...")
                return
            else:
                keep[column_list] = 0

        #Create list of columns to delete
        list_to_del = []
        for item in self.data_map:
            if item not in keep.keys():
                list_to_del.append(item)

        #Delete items in the list
        for item in list_to_del:
            del self.data_map[item]

    #This function is used to print processed data to an output file
    def printAlltoFile(self, file_name = ""):
        if file_name == "":
            file_name = self.input_file_name.split(".")[0]+"-output.dat"
        file = open(file_name,'w')
        file.write(str(self))
        file.write("\n")
        i=0
        first = ""
        for item in self.data_map:
            if i == 0:
                file.write(str(item))
                first = str(item)
            else:
                file.write("\t"+str(item))
            i+=1
        file.write("\n")
        j=0
        for value in self.data_map[first]:
            i = 0
            for item in self.data_map:
                if i == 0:
                    file.write(str(self.data_map[item][j]))
                else:
                    file.write("\t"+str(self.data_map[item][j]))
                i+=1
            file.write("\n")
            j+=1
        file.close()

    #This function is used to print select columns of data to a file
    def printColumnstoFile(self, column_list, file_name = ""):
        if type(column_list) is not list:
            print("Error! You must provide a list of columns to print to a file!")
            return
        for name in column_list:
            if name not in self.data_map.keys():
                print("Error! Invalid column names given...")
                return
        if file_name == "":
            file_name = self.input_file_name.split(".")[0]+"-SelectedOutput.dat"
        file = open(file_name,'w')
        file.write(str(self))
        file.write("\n")
        i=0
        for name in column_list:
            if i == 0:
                file.write(str(name))
            else:
                file.write("\t"+str(name))
            i+=1
        file.write("\n")
        j=0
        for value in self.data_map[column_list[0]]:
            i = 0
            for name in column_list:
                if i == 0:
                    file.write(str(self.data_map[name][j]))
                else:
                    file.write("\t"+str(self.data_map[name][j]))
                i+=1
            file.write("\n")
            j+=1
        file.close()

    #This function compresses the rows of the data map
    #   New columns are created to store compressed data and old columns are deleted to reserve space
    #   The 'factor' argument is optional and is used to determine how much compression to use
    #       Default is 2x compression: ==>  Cuts number of rows in half
    def compressRows(self, factor = 2):
        new_key_list = {}   #Map that links old key values (item) to the new keys
        factor = int(factor)
        for item in self.data_map:
            new_key = str(item) + "-" + str(int(factor)) + "x Compression"
            new_key_list[item] = new_key

        #time_key is forced to change here
        self.time_key += "-" + str(int(factor)) + "x Compression"

        #Compression will work by reading the first 'factor' # of values from a list and averaging them
        #   That average value becomes the new value to place into the new data_map with a new key
        #   This cycle repeats until no data is left to compress
        for item in new_key_list:
            self.data_map[new_key_list[item]] = []
            i = 0
            reset = 1
            avg = 0
            for value in self.data_map[item]:
                if type(value) is not int and type(value) is not float:
                    if reset == factor:
                        reset = 1
                        self.data_map[new_key_list[item]].append(value)
                    else:
                        reset+=1
                else:
                    if reset == factor:
                        avg += value
                        avg = avg/float(factor)
                        reset = 1
                        self.data_map[new_key_list[item]].append(avg)
                        avg = 0
                    else:
                        avg+=value
                        reset+=1
                i+=1
            del self.data_map[item]
        self.num_rows = len(self.data_map[self.time_key])


    #This function will take in a data_map key name and a list of changed values
    #   to add to the map of input_change for each input that corresponds to an
    #   output in data_map for the size of the change_time list
    #
    #       NOTE: Make sure you give the same units for the data in the input_list
    #               as the units provided in the corresponding output in data_map
    def registerChangedInput(self, data_key, input_list):
        #First, check to make sure that the data_key is valid
        if data_key not in self.data_map.keys():
            print("Error! No corresponding output value exists in data_map!")
            return

        #Next, check to make sure that the input_list size is the same as the change_time size
        if len(input_list) != len(self.change_time):
            print("Error! List of given changed inputs does not match length of change_time")
            return

        self.input_change[data_key+'[input]'] = input_list

    #This function automates the above function by utiliizing the corresponding
    #   output information of the given data_key to automatically approximate the
    #   input data for each change_time. That input data is estimated by averaging
    #   the last few output data points within the corresponding time range.
    #   By default, the last few data points are taken as the last 10 data points,
    #   however, you can override this by simply calling this function with a
    #   different value for avg_points. You may also specify whether or not the
    #   inlet conditions for this data set should be non-negative (e.g., for things
    #   such as inlet concentrations or molefractions)
    def autoregChangedInput(self, data_key, avg_points = 10, non_neg = True):
        #First, check to make sure that the data_key is valid
        if data_key not in self.data_map.keys():
            print("Error! No corresponding output value exists in data_map!")
            return

        #Loop through the data_map[time_key] data in REVERSE to get last data first
        points = 0
        change_loc = len(self.change_time)-1
        value_sum = 0
        time_index = len(self.data_map[self.time_key])-1
        has_calc = False
        #Initialize the list of inlet conditions (because we fill it in backwards)
        avg_list = [0.0]*len(self.change_time)
        for time in reversed(self.data_map[self.time_key]):
            if points >= avg_points:
                #Run a calculation here
                if has_calc == False:
                    value_sum = value_sum/avg_points
                    has_calc = True
                    if value_sum < 0 and non_neg == True:
                        value_sum = 0
                    avg_list[change_loc] = value_sum
                value_sum = 0
                #Update only if we are in next time bin
                if self.change_time[change_loc] >= time:
                    change_loc = change_loc - 1
                    #this is always the first point to grab, unless it is the first iteration
                    value_sum += self.data_map[data_key][time_index]
                    has_calc = False
                    points = 1
            else:
                #grab more data
                value_sum += self.data_map[data_key][time_index]
                points+=1
            time_index = time_index - 1
        #End time loop
        self.registerChangedInput(data_key, avg_list)

    #This function will create a new column in the data_map by creating a step-wise set of
    #   input data based on the change_time and corresponding input_change information.
    #   When calling this function, it is unnecessary to call registerChangedInput as
    #   this function will automatically perform the associated actions of that function.
    #       NOTE: data_key can be a single column or a list of columns
    def createStepChangeInputData(self, data_key, avg_points = 10, non_neg = True):
        if type(data_key) is not list:
            self.autoregChangedInput(data_key,avg_points,non_neg)
        else:
            for key in data_key:
                self.autoregChangedInput(key,avg_points,non_neg)
        #Iterate through the input_change map
        for new_key in self.input_change:
            i=0
            self.data_map[new_key] = []
            for time in self.data_map[self.time_key]:
                try:
                    time_limit = self.change_time[i+1]
                except:
                    time_limit = self.change_time[-1]+time
                if time > time_limit:
                    i+=1
                self.data_map[new_key].append(self.input_change[new_key][i])


    #This function will perform a trapezoid rule integration between two given curves
    #   in the data_map versus the time_key set of data. The first value of the integrated
    #   curve is always assumed to be zero. The units for the given columns do not matter
    #   as this method will specifically produce a normalized integrated curve. Generally,
    #   this function is used to create a data column for Mass Retained in the catalyst.
    #
    #   The following relationship is assumed...
    #
    #   d(MR)/dt = Q*(Min - Mout)
    #
    #   Min = Mass in (given data column to represent inlet mass)
    #   Mout = Mass out (given data column to represent outlet mass)
    #   Q = flow rate (usually as space velocity [ hr^-1 ])
    #   MR = Mass retained in the catalyst (Representative of adsorbed mass)
    def calculateRetentionNormalizedIntegral(self, inlet_column, outlet_column):
        if inlet_column not in self.data_map.keys():
            print("Error! No corresponding inlet column exists in data_map!")
            return
        if outlet_column not in self.data_map.keys():
            print("Error! No corresponding outlet column exists in data_map!")
            return
        if self.have_flow_rate == False:
            print("Error! Cannot calculate retention integral without a flow rate. Flow rate is expected in the file name in kilo-volumes per hour...")
            return

        ret_key = inlet_column.split()[0]+"-Retained (normalized)"
        self.data_map[ret_key] = []
        self.data_map[ret_key].append(0)
        MR_old = 0
        max_value = abs(MR_old)
        time_old = self.data_map[self.time_key][0]
        Min_old = self.data_map[inlet_column][0]
        Mout_old = self.data_map[outlet_column][0]
        i=1
        while i<len(self.data_map[self.time_key]):
            time_new = self.data_map[self.time_key][i]
            Min_new = self.data_map[inlet_column][i]
            Mout_new = self.data_map[outlet_column][i]
            MR_new = MR_old + (time_new-time_old)*self.flow_rate*( (Min_old+Min_new)/2 - (Mout_old+Mout_new)/2 )
            if abs(MR_new) > max_value:
                max_value = abs(MR_new)
            self.data_map[ret_key].append(MR_new)
            time_old = time_new
            Min_old = Min_new
            Mout_old = Mout_new
            MR_old = MR_new
            i+=1

        #Loop one last time to normalize the integrated curve
        i=0
        for value in self.data_map[ret_key]:
            self.data_map[ret_key][i] = self.data_map[ret_key][i]/max_value
            i+=1




## ---------------- End: Definition of TransientData object ------------



## ---------------- Begin: Definition of PairedTransientData object ------------
''' This object is used when you want to 'pair' inlet and outlet data sets together,
    as well as perform some post-processing such as integrals over data sets. To initialize
    the data set, you must pass a data file that contains the inlet data. For the
    CLEERS data sets, an inlet data file is denoted by a "-bp" at the end of the file name
    as opposed to an isothermal temperature.
'''


## ---------------- End: Definition of PairedTransientData object ------------
class PairedTransientData(object):
    #When creating an instance of this object, you must pass to files to the constructor
    #   (i)     A file for the input data (or by-pass run data)
    #   (ii)    A file for the output data (or non-by-pass run data)
    #The constructor will check the file names to make sure that the given information aligns
    #   so that the given by-pass and non-by-pass data sets should actually be paired.
    #   Because of this check, maintaining the same file name conventions as before is necessary.
    def __init__(self, bypass_file, result_file):

        # The constructor for the TransientData will automatically read the files
        self.bypass_trans_obj = TransientData(bypass_file)
        self.result_trans_obj = TransientData(result_file)
        self.aligned = False        #Flag used to determine whether or not the data sets are aligned in time

        # Check some specific file information to make sure there are no errors
        self.file_errors = False
        if self.bypass_trans_obj.inlet_data == False:
            print("Error! Given by-pass data file is not recognized as by-pass data from the given file name...")
            self.file_errors = True
        if self.result_trans_obj.inlet_data == True:
            print("Error! Given results data file is not recognized as results data from the given file name...")
            self.file_errors = True
        if self.bypass_trans_obj.material_name != self.result_trans_obj.material_name:
            print("Error! Files given indicate they are for 2 different materials. They cannot be paired...")
            self.file_errors = True
        if self.bypass_trans_obj.aging_condition != self.result_trans_obj.aging_condition:
            print("Error! Files given indicate they are for 2 different aging conditions. They cannot be paired...")
            self.file_errors = True
        if self.bypass_trans_obj.flow_rate != self.result_trans_obj.flow_rate:
            print("Error! Files given indicate they are for 2 different flow rates. They cannot be paired...")
            self.file_errors = True

        #Check to make sure that the column names are the same for the by-pass and results data
        for item in self.bypass_trans_obj.data_map:
            if item not in self.result_trans_obj.data_map:
                print("Error! By-pass file and Results file must have all the same columns and column names...")
                self.file_errors = True
                break

    #Function to print file information message to console
    def __str__(self):
        message = "\nBy-pass Information\n"
        message += "-------------------\n"
        message += str(self.bypass_trans_obj)
        message += "\nResults Information\n"
        message += "-------------------\n"
        message += str(self.result_trans_obj)
        if self.aligned == False:
            message += "\n\tWARNING: Data sets are unaligned currently!\n\t\tRun alignData() before processing information...\n"
        return message

    #Function to display the column names to console
    def displayColumnNames(self):
        self.bypass_trans_obj.displayColumnNames()

    #Function to compress columns in both data sets
    def compressColumns(self):
        self.bypass_trans_obj.compressColumns()
        self.result_trans_obj.compressColumns()

    #Function to append a column to by-pass data
    def appendBypassColumn(self, column_name, data_set):
        self.bypass_trans_obj.appendColumn(column_name, data_set)

    #Function to append a column to result data
    def appendResultColumn(self, column_name, data_set):
        self.result_trans_obj.appendColumn(column_name, data_set)

    #Function to extract a map of columns from the by-pass data
    def extractBypassColumns(self, column_list):
        return self.bypass_trans_obj.extractColumns(column_list)

    #Function to extract a map of columns from the result data
    def extractResultColumns(self, column_list):
        return self.result_trans_obj.extractColumns(column_list)

    #Function to extract a set of rows from the by-pass data
    def extractBypassRows(self, min_time, max_time):
        return self.bypass_trans_obj.extractRows(min_time,max_time)

    #Function to extract a set of rows from the result data
    def extractResultRows(self, min_time, max_time):
        return self.result_trans_obj.extractRows(min_time,max_time)

    #Function to get a data point from the by-pass data
    def getBypassDataPoint(self, time_value, column_name):
        return self.bypass_trans_obj.getDataPoint(time_value, column_name)

    #Function to get a data point from the result data
    def getResultDataPoint(self, time_value, column_name):
        return self.result_trans_obj.getDataPoint(time_value, column_name)

    #Function to delete a set of columns from both data sets
    def deleteColumns(self, column_list):
        self.bypass_trans_obj.deleteColumns(column_list)
        self.result_trans_obj.deleteColumns(column_list)

    #Function to delete a set all columns from both data sets, except for the specified set
    def retainOnlyColumns(self, column_list):
        self.bypass_trans_obj.retainOnlyColumns(column_list)
        self.result_trans_obj.retainOnlyColumns(column_list)

    #Function to align the bypass and results data so each has the same number of rows at the appropriate time values
    #NOTE: This function is riddled with comment lines and print statements that are used for debugging
    #       DO NOT REMOVE ANY COMMENTS UNLESS THE FUNCTION IS FULLY TESTED AND APPROVED
    #           This function is exceedingly complicated, do not modify unless you know what you are doing
    def alignData(self, addNoise = True):
        if self.aligned == True:
            print("Data already aligned. Cannot re-align...")
            return
        # Identify which data set will be considered the 'master_set'
        #   Master set will be the data set that has the most data available
        #   The sub-set will be forced to conform to the master set
        # Our goal is to expand the data in the sub-set to match the master set through
        #   linear interpolations or other means. Thus, the master set of data won't be
        #   changed, only the sub-set of data will change.
        master_set = {}             #Point to master data
        sub_set = {}                #Point to non-master data
        new_set = {}                #New map to make to override non-master data when finished
        master_change_time = []     #Point to list of master change_time
        sub_change_time = []        #Point to list of non-master change_time
        isResultMaster = False      #True if result_trans_obj is master, False if bypass_trans_obj is master
        time_column = self.bypass_trans_obj.time_key
        # Find master and point master's data_map to master_set
        #       This is so we only have to write code to iterate on master set
        #       Do not copy, only point. Master set will not change, so a pointer works best
        # Do same with sub_set, point sub_set to the non-master data set
        #       However, we also need a new data set to create. That new data set
        #       will eventually override the non-master original data
        if len(self.bypass_trans_obj.data_map[time_column]) > len(self.result_trans_obj.data_map[time_column]):
            isResultMaster = False
            master_set = self.bypass_trans_obj.data_map
            for i in range(1,len(self.bypass_trans_obj.change_time)):
                master_change_time.append((self.bypass_trans_obj.change_time[i-1],self.bypass_trans_obj.change_time[i]))
            master_change_time.append((self.bypass_trans_obj.change_time[-1],self.bypass_trans_obj.data_map[time_column][-1]))

            sub_set = self.result_trans_obj.data_map
            for i in range(1,len(self.result_trans_obj.change_time)):
                sub_change_time.append((self.result_trans_obj.change_time[i-1],self.result_trans_obj.change_time[i]))
            sub_change_time.append((self.result_trans_obj.change_time[-1],self.result_trans_obj.data_map[time_column][-1]))
        else:
            isResultMaster = True
            master_set = self.result_trans_obj.data_map
            for i in range(1,len(self.result_trans_obj.change_time)):
                master_change_time.append((self.result_trans_obj.change_time[i-1],self.result_trans_obj.change_time[i]))
            master_change_time.append((self.result_trans_obj.change_time[-1],self.result_trans_obj.data_map[time_column][-1]))

            sub_set = self.bypass_trans_obj.data_map
            for i in range(1,len(self.bypass_trans_obj.change_time)):
                sub_change_time.append((self.bypass_trans_obj.change_time[i-1],self.bypass_trans_obj.change_time[i]))
            sub_change_time.append((self.bypass_trans_obj.change_time[-1],self.bypass_trans_obj.data_map[time_column][-1]))

        #Check to make sure that the change times are already aligned
        if len(master_change_time) != len(sub_change_time):
            print("Error! Cannot align data if the number of 'change_time' instances are not the same!")
            return

        #Initialize the new_set with same number of columns and column names,
        #   but with the number of data points matching master_set
        for item in sub_set:
            #Check the data type of the master list at the given item and initialize new_set accordingly
            if type(master_set[item][0]) is not int and type(master_set[item][0]) is not float:
                new_set[item] = [""]*len(master_set[item])
            else:
                new_set[item] = [0.0]*len(master_set[item])

        #First data point in the new_set should auto-align to the first master_set time,
        #   regardless of discrepencies in the first time point. This is because
        #   we have no other information yet to perform any necessary interpolation.
        for item in new_set:
            new_set[item][0] = sub_set[item][0]
        new_set[time_column][0] = master_set[time_column][0]

        #Data from sub_set should be either compressed or elongated depending on the change_times
        '''
        i=0
        for tuple in sub_change_time:
            sub_dist = tuple[1]-tuple[0]
            master_dist = master_change_time[i][1] - master_change_time[i][0]
            if master_dist >= sub_dist:
                print("Elongate Sub (or keep same)")
            else:
                print("Compress Sub")
            i+=1
        '''

        #Want to iterate through all rows in the new data set
        print("\nData Alignment has started... Please wait...")
        tup_i = 0       #Index for tuple frame
        i = 1           #Index for master row
        j = 1           #Index for sub row
        #The time shifted values below represent how far in time we have traversed for the given tuple frame
        master_time_shifted = 0     #Start with no time shift for first reference frame
        sub_time_shifted = 0        #Same reason as above
        points_in_frame = 0         #Counter to keep track of the total number of points in the sub_set
                                    #   that exist in the current master_frame
        new_frame = False
        avg = {}
        std = {}
        hasCalcMean = False
        while i < len(new_set[time_column]):
            #Enfore matching time values
            new_set[time_column][i] = master_set[time_column][i]

            #Find when to update the tuple frame
            if master_set[time_column][i] > master_change_time[tup_i][1]:
                # ------------ Update tuple frame -------------
                new_frame = True
                tup_i +=1
                #When the tuple frame updates, we are now forced to update j
                #   j gets updated to the first indexed time value within the tuple frame
                for n in range(j,len(sub_set[time_column])):
                    if sub_change_time[tup_i][0] <= sub_set[time_column][n]:
                        j = n+1
                        break

                #At this point, i and j are aligned to their own tuple frames
                #       However, each may be misaligned with each other
                #Each time we move into a tuple frame, we should use a time_shited
                #       value for the actual time to represent distance traveled in
                #       the given frame of reference.
                #This resets the time shifting at each tuple frame update
                master_time_shifted = master_set[time_column][i-1] - master_change_time[tup_i][0]
                sub_time_shifted = sub_set[time_column][j-1] - sub_change_time[tup_i][0]
            else:
                new_frame = False
            #END Update for tuple frame

            #Calculate distance in time of our position in the master set of data from the previous point
            master_dt = master_set[time_column][i] - master_set[time_column][i-1]
            #Calculate the master time within the frame (not the actual time value)
            master_time_shifted += master_dt

            #Calculate target time for sub_set based on time travel and current tuple frame
            target_time = sub_time_shifted + master_time_shifted
            #This line was needed for additional alignment because our first frame starts at 0, but first data starts  > 0
            if tup_i == 0:
                target_time += sub_set[time_column][j-1]

            #Check to see if target time for sub_set is outside of the correct tuple frame
            #   If it is outside of the frame, then extrapolate values
            if target_time > (sub_change_time[tup_i][1]-sub_change_time[tup_i][0]):
                # --------------- Target Out of Frame ------------------
                # When the target is out of the tuple frame, we need to supplement the data to
                #   fill in any gaps. Simplest way to fill the gap is just to repeat the last
                #   relevant data point from the previous frame. However, This is likely
                #   unrealistic and may produce error.
                #
                #   Other options to consider:
                #       (i) Grab the last point and extend it to the next frame
                #       (ii) Average the last few points of the current frame to use as an extension
                #       (iii) Use a running average with artifical "noise" added for realism
                for item in new_set:
                    if item is not time_column:
                        if type(new_set[item][i]) is not int and type(new_set[item][i]) is not float:
                            #If the data in the column is non-numeric, then just grab the prior time value
                            new_set[item][i] = new_set[item][i-1]
                        else:
                            #If the data in the column is numeric, then just grab the prior time value
                            #new_set[item][i] = new_set[item][i-1]

                            #If the data in the column is numeric, then just grab average the last 1/4 of elements
                            if addNoise == False:
                                new_set[item][i] = new_set[item][i-1]
                                if hasCalcMean == False:
                                    avg[item] = mean(new_set[item][i-int(points_in_frame/4):i])
                                new_set[item][i] = avg[item]
                            #If the data in the column is numeric, then average the last few points and
                            #   and in some random noise based on the standard deviation of the last few points
                            else:
                                new_set[item][i] = new_set[item][i-1]
                                if hasCalcMean == False:
                                    avg[item] = mean(new_set[item][i-int(points_in_frame/4):i])
                                    std[item] = stdev(new_set[item][i-int(points_in_frame/4):i])
                                new_set[item][i] = random.gauss(avg[item],std[item])
                    else:
                        new_set[item][i] = master_set[item][i]
                hasCalcMean = True
            #   If it is inside of the frame, then interpolate or extract values
            else:
                #----------------- Data Mapping ---------------------
                #master_set[time_column][i]  ---> target_time+sub_change_time[tup_i][0]
                # Current time in master set maps to the target time in the sub_set
                #       plus the off-set dictated by the tuple frame of the sub_set
                hasCalcMean = False
                if new_frame == False:
                    points_in_frame += 1
                else:
                    points_in_frame = 1
                for item in new_set:
                    # Perform a linear interpolation for the correct data_map and place into the new_set
                    if item is not time_column:
                        if isResultMaster == True:
                            new_set[item][i] = self.bypass_trans_obj.getDataPoint(target_time+sub_change_time[tup_i][0],item)
                        else:
                            new_set[item][i] = self.result_trans_obj.getDataPoint(target_time+sub_change_time[tup_i][0],item)
                    else:
                        new_set[item][i] = master_set[item][i]
            i+=1
        #END while loop through all times

        #At this point, the new_set should contain all of the data for the sub_set that needs alignment
        #   Delete all existing data in the sub_set and replace with the new_set
        for item in new_set:
            del sub_set[item]
        for item in new_set:
            sub_set[item] = []
            for value in new_set[item]:
                sub_set[item].append(value)

        #If we made it this far without error, then the data should be aligned (fingers crossed, no jynx)
        print("Complete!")
        self.aligned = True
        if isResultMaster == True:
            self.bypass_trans_obj.num_rows = len(sub_set[time_column])
        else:
            self.result_trans_obj.num_rows = len(sub_set[time_column])

        #Lastly, for book-keeping purposes, point all results columns to the corresponding bypass
        #       columns by appending the columns with the appropriate column names
        #After this method has been called, result and bypass structures no longer
        #       have identical numbers of columns and column names.
        for item in self.bypass_trans_obj.data_map:
            if item != self.bypass_trans_obj.time_key:
                appendName = item+"[bypass]"
                self.result_trans_obj.appendColumn(appendName, self.bypass_trans_obj.data_map[item])
    #End alignData()


    #NOTE: Before you can compress the rows of data, each data set must be aligned
    #Function will compress the rows of data based on the given compression factor
    def compressRows(self, factor = 2):
        if self.aligned == False:
            print("Error! Cannot compress the rows of data until data is aligned! Run alignData() first...")
            return

        self.bypass_trans_obj.compressRows(factor)
        self.result_trans_obj.compressRows(factor)

    #Function will compute a normalized mass retention integral for the given column
    #   The align data function must have already been called. That function
    #   ensures the data sets for bypass and results are aligned and appends the
    #   aligned data from bypass to the results with "[bypass]" added to the
    #   end of the file name.
    def calculateRetentionNormalizedIntegral(self, column_name):
        #Check for alignment
        if self.aligned == False:
            print("Error! Data must be aligned first! Call alignData()...")
            return
        appendName = column_name+"[bypass]"
        if appendName not in self.result_trans_obj.data_map.keys():
            print("Error! The appended column_name does not match any by-pass/result data keys!")
            return
        if column_name not in self.result_trans_obj.data_map.keys():
            print("Error! The column_name does not match any result data keys!")
            return
        self.result_trans_obj.calculateRetentionNormalizedIntegral(appendName, column_name)

    #Function will print out the results to an sinle output file
    #   Data must already be aligned. This allows us to only print out information
    #   in the results object, since after alignment the results object is linked
    #   to the columns in the bypass object.
    def printAlltoFile(self, file_name = ""):
        if self.aligned == False:
            print("Error! Data must be aligned first! Call alignData()...")
            print("\t\t(NOTE: if you want to print out the bypass and result data separately,")
            print("\t\t       then call each object's s respective printAlltoFile() functions)")
            return
        if file_name == "":
            file_name = self.result_trans_obj.input_file_name.split(".")[0]+"-AllPairedOutput.dat"
        self.result_trans_obj.printAlltoFile(file_name)

    #Function to print only specific columns to an output file
    def printColumnstoFile(self, column_list, file_name = ""):
        if self.aligned == False:
            print("Error! Data must be aligned first! Call alignData()...")
            print("\t\t(NOTE: if you want to print out the bypass and result data separately,")
            print("\t\t       then call each object's s respective printAlltoFile() functions)")
            return
        if file_name == "":
            file_name = self.result_trans_obj.input_file_name.split(".")[0]+"-SelectedPairedOutput.dat"
        self.result_trans_obj.printAlltoFile(column_list, file_name)



## ------ Testing ------

# Testing with the NH3 TPDs at constant H2O concentration
'''
test01 = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-150C.dat")

test01.compressColumns()
#NOTE: After compressing the columns, some column names with change to reflect the combination of data
#test01.displayColumnNames()
test01.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)', 'TC bot sample in (C)', 'TC bot sample mid 1 (C)', 'TC bot sample mid 2 (C)', 'TC bot sample out 1 (C)', 'TC bot sample out 2 (C)', 'P bottom in (bar)', 'P bottom out (bar)'])
test01.createStepChangeInputData('NH3 (300,3000)')
test01.calculateRetentionNormalizedIntegral('NH3 (300,3000)[input]','NH3 (300,3000)')
#NOTE: Consider using current number of rows to determine how much compression to use
#print(test01.num_rows)
test01.compressRows(10)
test01.printAlltoFile()
'''

# Testing with the Competition between H2O and NH3 TPDs
'''
test02 = TransientData("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-150C.dat")
test02.compressColumns()
#test02.displayColumnNames()
test02.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)', 'TC bot sample in (C)', 'TC bot sample mid 1 (C)', 'TC bot sample mid 2 (C)', 'TC bot sample out 1 (C)', 'TC bot sample out 2 (C)', 'P bottom in (bar)', 'P bottom out (bar)'])
test02.createStepChangeInputData(['NH3 (300,3000)','H2O% (20)'],40)
test02.calculateRetentionNormalizedIntegral('NH3 (300,3000)[input]','NH3 (300,3000)')
test02.calculateRetentionNormalizedIntegral('H2O% (20)[input]','H2O% (20)')
#print(test01.num_rows)
#print("H2O-150")
#print(test02.change_time)
test02.compressRows(2)
test02.printAlltoFile()
'''


#Testing with by-pass data
'''
test03 = TransientData("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-bp.dat")
test03.compressColumns()
test03.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)'])
#print(test03.change_time)
test03.compressRows(1)
test03.printAlltoFile()

test04 = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-bp.dat")
test04.compressColumns()
test04.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)'])
#print(test04.change_time)
test04.compressRows(1)
test04.printAlltoFile()
'''


#Testing Paired data
'''
h2o_comp = False
if h2o_comp == True:
    test05 = PairedTransientData("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-bp.dat","20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-150C.dat")
else:
    test05 = PairedTransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-bp.dat","20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-150C.dat")
test05.compressColumns()
#test05.displayColumnNames()

#NOTE: If you are going to only retain specific columns, you should run that function
#       prior to running alignData(). This saves significant computational effort.
test05.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)'])
test05.alignData()

#NOTE: Always perform data processing BEFORE compressing the rows!
test05.calculateRetentionNormalizedIntegral('NH3 (300,3000)')
test05.calculateRetentionNormalizedIntegral('H2O% (20)')

if h2o_comp == True:
    test05.compressRows(2)
else:
    test05.compressRows(10)

test05.printAlltoFile()
'''

## ----- End Testing -----
