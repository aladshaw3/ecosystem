## @example Example of the usage of the yaml_reader
#
# Using the reader...
#
# import yaml_reader \n
# yaml = yaml_reader.YAML() \n
# yaml.readFile("test_input.yml") \n
# yaml.print2file("test_out.txt") \n
#
#
# Manually creating a yaml map and printing to a file...
#
# yaml3 = yaml_reader.YAML() \n
# yaml3.map["Doc1"] = {} \n
# yaml3.map["Doc2"] = {} \n
# yaml3.map["Doc1"]["key1"] = "val1" \n
# yaml3.map["Doc1"]["key2"] = False \n
# yaml3.map["Doc2"]["head1"] = {} \n
# yaml3.map["Doc2"]["head1"]["head_key"] = 10 \n
# yaml3.map["Doc2"]["head1"]["sub"] = {} \n
# yaml3.map["Doc2"]["head1"]["sub"]["key"] = 5.5 \n

yaml3.print2file("test2.txt")

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
