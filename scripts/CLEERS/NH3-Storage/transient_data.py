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
        message = "File Name: " + self.data_file.name
        message += "\nFile Header: " + self.exp_header
        message += "Number of Columns: " + str(len(self.data_map))
        message += "\nNumber of Rows: " + str(len(self.data_map['Elapsed Time (min)']))
        return message

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

## ---------------- End: Definition of TransientData object ------------


## ------ Testing ------
test = TransientData("20160205-CLRK-BASFCuSSZ13-700C4h-NH3DesIsoTPD-30k-0_2pctO2-5pctH2O-150C.dat")
print(test)

## ----- End Testing -----
