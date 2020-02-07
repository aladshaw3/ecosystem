''' Python script to read in CLEERS transient data for
    NH3 storage on Cu-SSZ-13. This script will store the
    orginal data as is and provide other functions to
    redistribute, print, or parse that data as needed. '''

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
    def __init__(self, file):
        self.data_file = open(file,"r") #Contains data file we are digitizing
        self.exp_header = ''            #Contains the first line of the data file
        self.data_map = {}              #Contains a map of all the data by column
        self.ordered_key_list = []      #Contains an ordered list of column names
        self.change_time = []           #Contains an ordered list of the times when experimental inputs changed
        self.input_change = {}          #Contains a map of the inputs values that correspond to change_time
                                        #   Keys of this map are modifications to the keys of data_map that
                                        #   correspond to output values corresponding to the input values given
        self.readFile()
        self.closeFile()

    def __str__(self):
        message = "\nFile Name: " + self.data_file.name
        message += "\nFile Header: " + self.exp_header
        message += "Number of Columns: " + str(len(self.data_map))
        message += "\nNumber of Rows: " + str(len(self.data_map['Elapsed Time (min)']))
        return message + "\n"

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
        for item in column_list:
            #Check to make sure the item is a key in data_map
            if item in self.data_map.keys():
                new_map[item] = self.data_map[item]
            else:
                print("Error! Invalid Key!")
        #End item loop
        return new_map

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
test01 = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-150C.dat")
print(test01)
#print(test01.change_time)
#test01.registerChangedInput('NH3 (3000) ',[0,1000,800,600,400,200,100,50,25,12.5,0])
#print(test01.input_change['NH3 (3000) [input]'])
test01.autoregChangedInput('NH3 (3000) ')
print(test01.change_time)
print(test01.input_change)

#test02 = TransientData("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-150C.dat")
#print(test02)

'''
set = []
i = 0
for item in test02.data_map['Elapsed Time (min)']:
    set.append(i)
    i+=1
test02.appendColumn('new',set)
print(test02.data_map['new'])
'''
#print(test02.data_map.keys())

#map = test02.extractColumns( ['Elapsed Time (min)','NH3 (3000) ','H2O% (20) '] )

#print(map)

#Direct access to class data is allowed
#print(test02.data_map['Elapsed Time (min)'][-1])
## ----- End Testing -----
