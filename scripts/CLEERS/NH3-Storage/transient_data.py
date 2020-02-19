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
    #       item[5] = flow rate information (in volumes per hour)
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
                                self.data_map[self.ordered_key_list[n]].append(item)
                    n+=1
            i+=1
        #END of line loop
        self.num_rows = len(self.data_map[self.time_key])
        self.closeFile()

    def closeFile(self):
        self.data_file.close()

    #This function will add a column to the data map given the column name and associated data
    def appendColumn(self, column_name, data_set):
        #First, check to make sure the map has been prepared
        if (len(self.data_map) == 0):
            print("Error! File has not been read and stored!")
            return

        #Next, check to make sure the column_name is not already in the map
        if column_name in self.data_map.keys():
            print("Error! That data column is already in the structure!")
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
            point = (end_point - start_point)/(end_time - start_time)*(time_value - start_time) + start_point
        return point

    #This function will iterate through all columns to find data that can be compressed or eliminated
    #       For instance,   if a column contains no data, then delete it
    #                       if there are multiple columns that carrier similar info, then combine them
    def compressColumns(self):
        #Iterate through the data map to find information to delete immediately
        list_to_del = []
        for item in self.data_map:
            if len(self.data_map[item]) == 0:
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

        #Check for max values that are very small and make some kind of correction
        if max_value == 0:
            print("here")
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


## ------ Testing ------

# Testing with the NH3 TPDs at constant H2O concentration
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

# Testing with the Competition between H2O and NH3 TPDs
test02 = TransientData("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-150C.dat")
test02.compressColumns()
#test02.displayColumnNames()
test02.retainOnlyColumns(['Elapsed Time (min)','NH3 (300,3000)', 'H2O% (20)', 'TC bot sample in (C)', 'TC bot sample mid 1 (C)', 'TC bot sample mid 2 (C)', 'TC bot sample out 1 (C)', 'TC bot sample out 2 (C)', 'P bottom in (bar)', 'P bottom out (bar)'])
test02.createStepChangeInputData(['NH3 (300,3000)','H2O% (20)'])
test02.calculateRetentionNormalizedIntegral('NH3 (300,3000)[input]','NH3 (300,3000)')
test02.calculateRetentionNormalizedIntegral('H2O% (20)[input]','H2O% (20)')
#print(test01.num_rows)
test02.compressRows(2)
test02.printAlltoFile()
print(test02)

## ----- End Testing -----
