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
        self.exp_header = ''            #Contains the first line of the data file
        self.data_map = {}              #Contains a map of all the data by column
        self.ordered_key_list = []      #Contains an ordered list of column names
        self.change_time = []           #Contains an ordered list of the times when experimental inputs changed
                                        #   NOTE: this is just used to initialize data in the map, it does not
                                        #           change or update if the map changes
        self.input_change = {}          #Contains a map of the inputs values that correspond to change_time
                                        #   Keys of this map are modifications to the keys of data_map that
                                        #   correspond to output values corresponding to the input values given
        self.readFile()

    def __str__(self):
        message = "\nFile Name: " + self.data_file.name
        message += "\nMaterial: " + self.material_name
        message += "\nAging Condition: " + self.aging_condition
        message += "\nFlow Rate (hr^-1): " + str(self.flow_rate)
        message += "\nInlet Conditions: "
        if self.inlet_data == True:
            message += "True"
        else:
            message += "False"
            message += "\nIsothermal Temp (C): " + str(self.isothermal_temp)
        message += "\nFile Header: " + self.exp_header
        message += "\nNumber of Columns: " + str(len(self.data_map))
        message += "\nNumber of Rows: " + str(len(self.data_map['Elapsed Time (min)']))
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
                        self.data_map[item] = []
                        self.ordered_key_list.append(item)
                    else:
                        # Force a new column for input conditions
                        self.change_time.append(0)

            if (i > 1): #or greater
                n = 0
                for item in line_list:
                    #Check item to make sure that we do not get another column key
                    # If we do, then do not append data to map, instead record
                    #   the 'Elapsed Time (min)' value from prior as the point
                    #   when input conditions changed
                    if item in self.data_map.keys():
                        if (n == 0):
                            self.change_time.append(float(self.data_map['Elapsed Time (min)'][-1]))
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
        if len(self.data_map['Elapsed Time (min)']) != len(data_set):
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
        for time in self.data_map['Elapsed Time (min)']:
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
        for time in self.data_map['Elapsed Time (min)']:
            if time > time_value:
                if n == 0:
                    start_time = 0
                    start_point = self.data_map[column_name][n]
                else:
                    start_time = self.data_map['Elapsed Time (min)'][n-1]
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
                for sub_key in frac_keys[item]:
                    first = sub_key.split("(")
                    last = first[1].split(")")
                    val_list.append(float(last[0]))
                self.data_map[item] = []
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
                    #print(dist)
                    #print(reg_index)
                    self.data_map[item].append(self.data_map[frac_keys[item][reg_index]][n])
                    n+=1
                #End value loop
                print(self.data_map[item])
            #End if

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

        #Loop through the 'Elapsed Time (min)' data in REVERSE to get last data first
        points = 0
        change_loc = len(self.change_time)-1
        value_sum = 0
        time_index = len(self.data_map['Elapsed Time (min)'])-1
        has_calc = False
        #Initialize the list of inlet conditions (because we fill it in backwards)
        avg_list = [0.0]*len(self.change_time)
        for time in reversed(self.data_map['Elapsed Time (min)']):
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

## ---------------- End: Definition of TransientData object ------------


## ------ Testing ------
#test01 = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-150C.dat")
test01 = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-bp.dat")
print(test01)
#test01.displayColumnNames()
#print(test01.extractColumns("NH3 (3000) "))
#print(test01.extractRows(1.95,2))
print(test01.getDataPoint(1.27,"NH3 (3000) "))

test01.compressColumns()
#test01.displayColumnNames()

## ----- End Testing -----
