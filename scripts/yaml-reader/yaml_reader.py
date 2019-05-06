## Python script to read yaml files using ecosystem library functions ##
## Run python scripts using Python 3.5 or newer ##

''' YAML script:
    ----------------
    Object-Oriented approach to interfacing with the C/C++ functions and objects
    built into the ecosystem library (libeco.so). This script will provide a
    Python interface to already previously developed C++ set of objects for
    reading and storing a digital record of a yaml formatted file. It was specifically
    designed for the purpose of maintaining the same yaml styling as what is used
    by the ecosystem library.
    
    Author:     Austin Ladshaw
    Date:       05/06/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.'''

from ctypes import *
lib = cdll.LoadLibrary('../../libeco.so')

''' Python Class for interfacing with the C++ yaml wrapper
    ------------------------------------------------------
    This class provides the necessary interface for reading
    yaml formatted files. It uses the yaml wrapper class
    previously developed in the ecosystem project. This is
    done to ensure consistency between usage and application
    of yaml based input files within the entire project.
    
    NOTE: You need to direct python what the return types
        and argument types are for each C-style function.
        This is especially important when dealing with
        C++ objects, as you will need to use void * as
        arguments and return types. 
    
        Pass const char* as arg.encode()
        Pass instances of C++ structures as c_void_p'''
class YAML(object):
    def __init__(self):
        initialize = lib.New_YAML
        initialize.restype = c_void_p
        self.obj = initialize()

    def readFile(self, file):
        readFunc = lib.YAML_executeYamlRead
        readFunc.restype = c_int
        readFunc.argtypes = [c_void_p, c_char_p]
        return readFunc(self.obj,file.encode())

    def displayContents(self):
        dispFunc = lib.YAML_DisplayContents
        dispFunc.argtypes = [c_void_p]
        return dispFunc(self.obj)

    def docKeys(self):
        docKeys_func = lib.YAML_DocumentKeys
        docKeys_func.restype = c_char_p
        docKeys_func.argtypes = [c_void_p]
        _res = docKeys_func(self.obj)
        result = str(_res)
        return result

yaml = YAML()
yaml.readFile("1979-Test-Case.txt")
yaml.displayContents()
keys = yaml.docKeys()
print(keys)
start = False
s_keys = ""
for char in keys:
    if char == "'" and start == False:
        start = True
    if start == True and char != "'":
        s_keys += char
print(s_keys)
