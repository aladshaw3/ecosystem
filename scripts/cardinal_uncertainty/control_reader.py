## Python script to read cardinal control files and make changes ##
## Run python scripts using Python 3.5 or newer ##

''' Control file script:
    ----------------
    Object-Oriented approach to reading the Cardinal control files
    This object will use the previously established python yaml reader
    to read the control files, then make changes to the python yaml map
    and record those changes in the new control file.
    
    Author:     Austin Ladshaw
    Date:       05/14/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.
    '''
import sys, os
sys.path.insert(0, '../yaml-reader')
import yaml_reader as yr
from enum import Enum

class ControlParam(Enum):
    Yield  = 0
    Height = 1
    Level = 2

class ControlFile(yr.YAML):
    def __init__(self):
        yr.YAML.__init__(self)

    def printMap(self):
        print(self.map)


## Testing ##
control = ControlFile()
control.readFile("1979-Test-Case.txt")
#print(control)
control.printMap()

