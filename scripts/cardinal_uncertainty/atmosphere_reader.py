## Python script to read atmosphere files and apply changes ##
## Run python scripts using Python 3.5 or newer ##

''' Atmosphere script:
    ----------------
    Object-Oriented approach to reading the atmosphere data files used by
    Cardinal. This object will read in the file and create an editable map
    of the file. That map can then be used to change the atmospheric conditions
    automatically for the uncertainty/sensitivity analyses.
    
    Author:     Austin Ladshaw
    Date:       05/14/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.
    '''
from enum import Enum

class AtmParam(Enum):
    Temperature  = 0
    Pressure = 1
    Humidity = 2

class Atmosphere(object):
    def __init__(self):
        self.map = {}

    def __str__(self):
        string = ""
        i = 0
        for alt in self.map:
            if i != 0:
                string += "\n"
            string += str(alt) + "\t"
            for item in self.map[alt]:
                string += str(item) + "\t"
            i += 1
        return string

    def readFile(self, file):
        for line in open(file,"r"):
            i = 0
            for item in line.split():
                if i == 0:
                    try:
                        alt = float(item)
                        self.map[alt] = []
                    except:
                        alt = item
                        self.map[alt] = []
                else:
                    try:
                        self.map[alt].append(float(item))
                    except:
                        self.map[alt].append(item)
                i += 1
        #End Loop over all lines

    def print2file(self, filename):
        file = open(filename, 'w')
        info = str(self)
        file.write(info)
        file.close()

    def editValue_PercentChange(self, e_num, per):
        i = 0
        for key in self.map:
            if i != 0:
                if e_num == AtmParam.Temperature:
                    self.map[key][0] = self.map[key][0]*(1.0 + (per/100.0))
                if e_num == AtmParam.Pressure:
                    self.map[key][1] = self.map[key][1]*(1.0 + (per/100.0))
                if e_num == AtmParam.Humidity:
                    self.map[key][2] = self.map[key][2]*(1.0 + (per/100.0))
            i += 1

    def editValue_LinearChange(self, e_num, value):
        i = 0
        for key in self.map:
            if i != 0:
                if e_num == AtmParam.Temperature:
                    self.map[key][0] = self.map[key][0] + value
                if e_num == AtmParam.Pressure:
                    self.map[key][1] = self.map[key][1] + value
                if e_num == AtmParam.Humidity:
                    self.map[key][2] = self.map[key][2] + value
            i += 1

    def clearMap(self):
        self.map.clear()

## Testing ##
'''
test = Atmosphere()
test.readFile("DefaultAtmosphere.txt")

test.editValue_PercentChange(AtmParam.Humidity, 10.0)
test.editValue_LinearChange(AtmParam.Temperature, -10.0)

test.print2file("DefaultAtmosphere_Mod.txt")

test2 = Atmosphere()
test2.readFile("DefaultAtmosphere_Mod.txt")
print(test2)
'''

