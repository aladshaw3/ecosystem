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
        self.map = {}

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
        docKeys_size = lib.YAML_DocumentKeys_Size
        docKeys_func = lib.YAML_DocumentKeys
        docKeys_size.restype = c_int
        docKeys_size.argtypes = [c_void_p]
        docKeys_func.argtypes = [c_void_p, c_char_p]
        size = docKeys_size(self.obj)
        keys = create_string_buffer(size)
        docKeys_func(self.obj, keys)
        long_keys = str(keys.value,'utf-8')
        list_keys = long_keys.split(":")
        for name in list_keys:
            self.map[name] = {}

    def headKeys(self):
        for doc in self.map:
            headKeys_size = lib.YAML_HeaderKeys_Size
            headKeys_func = lib.YAML_HeaderKeys
            headKeys_size.restype = c_int
            headKeys_size.argtypes = [c_void_p, c_char_p]
            headKeys_func.argtypes = [c_void_p, c_char_p, c_char_p]
            size = headKeys_size(self.obj,doc.encode())
            keys = create_string_buffer(size)
            headKeys_func(self.obj,doc.encode(),keys)
            long_keys = str(keys.value,'utf-8')
            list_keys = long_keys.split(":")
            for name in list_keys:
                # Don't add a map if there is no key
                if name != '':
                    self.map[doc][name] = {}

    def subKeys(self):
        for doc in self.map:
            for head in self.map[doc]:
                subKeys_size = lib.YAML_SubHeaderKeys_Size
                subKeys_func = lib.YAML_SubHeaderKeys
                subKeys_size.restype = c_int
                subKeys_size.argtypes = [c_void_p, c_char_p, c_char_p]
                subKeys_func.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p]
                size = subKeys_size(self.obj, doc.encode(), head.encode())
                keys = create_string_buffer(size)
                subKeys_func(self.obj,doc.encode(),head.encode(),keys)
                long_keys = str(keys.value,'utf-8')
                list_keys = long_keys.split(":")
                for name in list_keys:
                    # Don't add a map if there is no key
                    if name != '':
                        self.map[doc][head][name] = {}

    def docData(self):
        for doc in self.map:
            docData_size = lib.YAML_DocumentData_Size
            docData_func = lib.YAML_DocumentData
            docData_size.restype = c_int
            docData_size.argtypes = [c_void_p, c_char_p]
            docData_func.argtypes = [c_void_p, c_char_p, c_char_p]
            size = docData_size(self.obj, doc.encode())
            key_values = create_string_buffer(size)
            docData_func(self.obj, doc.encode(), key_values)
            long_keys = str(key_values.value,'utf-8')
            list_keys = long_keys.split("*")
            for items in list_keys:
                key = items.split(":")
                if key[0] != '':
                    self.map[doc][key[0]] = key[1]

    def headData(self):
        for doc in self.map:
            for head in self.map[doc]:
                headData_size = lib.YAML_HeaderData_Size
                headData_func = lib.YAML_HeaderData
                headData_size.restype = c_int
                headData_size.argtypes = [c_void_p, c_char_p, c_char_p]
                headData_func.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p]
                size = headData_size(self.obj, doc.encode(), head.encode())
                key_values = create_string_buffer(size)
                headData_func(self.obj, doc.encode(), head.encode(), key_values)
                long_keys = str(key_values.value,'utf-8')
                list_keys = long_keys.split("*")
                for items in list_keys:
                    key = items.split(":")
                    if key[0] != '':
                        self.map[doc][head][key[0]] = key[1]

    def subData(self):
        for doc in self.map:
            for head in self.map[doc]:
                for sub in self.map[doc][head]:
                    subData_size = lib.YAML_SubHeaderData_Size
                    subData_func = lib.YAML_SubHeaderData
                    subData_size.restype = c_int
                    subData_size.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p]
                    subData_func.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p, c_char_p]
                    size = subData_size(self.obj, doc.encode(), head.encode(), sub.encode())
                    key_values = create_string_buffer(size)
                    subData_func(self.obj, doc.encode(), head.encode(), sub.encode(), key_values)
                    long_keys = str(key_values.value,'utf-8')
                    list_keys = long_keys.split("*")
                    for items in list_keys:
                        key = items.split(":")
                        if key[0] != '':
                            self.map[doc][head][sub][key[0]] = key[1]

## END YAML Class ##


yaml = YAML()
#yaml.readFile("1979-Test-Case.txt")
yaml.readFile("test_input.yml")
#yaml.displayContents()
yaml.docKeys()
#print(yaml.map)
#print(yaml.map["ODE_Options"]) #this works
yaml.headKeys()
#print(yaml.map)

yaml.subKeys()
#print(yaml.map)

#yaml.docData()

#print(yaml.map)
#yaml.map["TestDoc3"]["apple"] = 5.0   # works
#print(yaml.map["TestDoc3"]["apple"]) # works

yaml.subData()
yaml.headData()
yaml.docData()  #Must call in reverse order
print(yaml.map)

