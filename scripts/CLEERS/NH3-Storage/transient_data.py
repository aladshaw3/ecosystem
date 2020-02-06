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

        return new_map


## ---------------- End: Definition of TransientData object ------------


## ------ Testing ------
test01 = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-150C.dat")
print(test01)

test02 = TransientData("20160209-CLRK-BASFCuSSZ13-700C4h-NH3H2Ocomp-30k-0_2pctO2-11-3pctH2O-400ppmNH3-150C.dat")
print(test02)

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

map = test02.extractColumns( ['Elapsed Time (min)','NH3 (3000) ','H2O% (20) '] )

#print(map)

#Direct access to class data is allowed
#print(test02.data_map['Elapsed Time (min)'][-1])
## ----- End Testing -----
