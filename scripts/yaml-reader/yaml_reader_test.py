## Python script to read yaml files using ecosystem library functions ##
## Run python scripts using Python 3.5 or newer ##

''' YAML test script:
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

import yaml_reader

yaml = yaml_reader.YAML()
#yaml.readFile("1979-Test-Case.txt")
yaml.readFile("test_input.yml")
#yaml.displayContents()
#yaml.docKeys()
#print(yaml.map)
#print(yaml.map["ODE_Options"]) #this works
#yaml.headKeys()
#print(yaml.map)

#yaml.subKeys()
#print(yaml.map)

#yaml.docData()

#print(yaml.map)
#yaml.map["TestDoc3"]["apple"] = 5.0   # works
#print(yaml.map["TestDoc3"]["apple"]) # works

#yaml.subData()
#yaml.headData()
#yaml.docData()  #Must call in reverse order
#print(yaml.map)
#yaml.formMap()
#print(yaml)
yaml.print2file("test_out.txt")

yaml2 = yaml_reader.YAML()
yaml2.readFile("test_out.txt")
#yaml2.formMap()
#print(yaml2)

yaml3 = yaml_reader.YAML()
yaml3.map["Doc1"] = {}
yaml3.map["Doc2"] = {}
yaml3.map["Doc1"]["key1"] = "val1"
yaml3.map["Doc1"]["key2"] = False
yaml3.map["Doc2"]["head1"] = {}
yaml3.map["Doc2"]["head1"]["head_key"] = 10
yaml3.map["Doc2"]["head1"]["sub"] = {}
yaml3.map["Doc2"]["head1"]["sub"]["key"] = 5.5

yaml3.print2file("test2.txt")

# Convert to int
print( int(yaml.map["TestDoc1"]["scenario"]["numvar"]) + 10)

# Convert to float/double
print( float(yaml.map["TestDoc1"]["scenario"]["pH"])*10 )

# Convert to bool#
print( (yaml.map["TestDoc1"]["scenario"]["steadystate"]) )

print( type(yaml.map["TestDoc1"]["scenario"]["numvar"]) )  ## Yes

print( (yaml.map["TestDoc1"]["scenario"]["t_out"]) + 5.5 )

# Now automatically converts to types: int, double/float, bool, string

