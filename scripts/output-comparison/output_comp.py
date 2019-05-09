## Python script to read output files and compare them ##
## Run python scripts using Python 3.5 or newer ##

''' YAML script:
    ----------------
    Script for reading in output produced from an executable and comparing it
    to another output file. Useful for creating small unit tests for code validations
    or making evaluations on how changes in parameters change predicted outcomes. Thus,
    this script will be utilized for the purpose of performing sensitivity analyses.
    
    Author:     Austin Ladshaw
    Date:       05/09/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.'''

import difflib

### Class Object for Comparing Files ###

class FileCompare(object):
    # Initialization constructor must take in strings for the gold file and test file
    def __init__(self, gold, test):
        self.num_diff = 0.0                 # Values closest to zero represent smallest differences
        self.str_ratio = 1.0                # Values closest to one represent smallest differences
        self.gold_file = open(gold,"r")
        self.test_file = open(test,"r")

    def closeFiles(self):
        self.gold_file.close()
        self.test_file.close()

### End Class Object ###

''' Notes 
    -----
    read()   gives all text in file
    read(n)  gives first n characters in file
    readline()   gives first line in file      (**Each instance reads the next line)
    readline(n)  gives nth line in the file
    readlines()  gives list of all lines in file
    for line in file:   loops over all lines in the file
    string.split(ch)    produces list of sub-strings parsed by the ch character
                        leaving ch blank defaults to splitting on blank spaces
    use == to compare two strings
        e.g.,  str1 == str2   -->   True if same, False if different
'''

### Testing ###
comp = FileCompare("gold_out.txt","test_out.txt")
line1 = comp.gold_file.readline().split()
linea = comp.test_file.readline().split()
print(line1)
print(len(line1))
i = 0
sum = 0
for str in line1:
    if str != linea[i]:
        sum += 1
    i += 1

print(sum)

d = difflib.Differ()
diff = d.compare(line1, linea)
print(list(diff))

hd = difflib.HtmlDiff()
hdiff = hd.make_file(line1, linea)
print(hdiff)

### Perform a comparison of string lines and report a numeric difference ratio ###
sd = difflib.SequenceMatcher()
sd.set_seqs(line1,linea)
### if ratio == 1 (SAME lines) if ratio == 0 (Completely Different)
# ratio = 2.0*M/T  where M = number of matching elements and T = total number of elements (in both).
print(line1)
print(linea)
print(sd.ratio())

''' Plan
    ----
    (1) Iterate through each line in each file
            How can we tell the total number of lines? (iterate through smallest)
            
    (2) Split each line into individual strings
            Parse on spaces
    
    (3) Iterate (nested) over each individual string in a lint
            (loop over the smallest string list)
    
    (4) Check each string for types
            a)  compare numbers/floats (try) and developed numeric different (sq. error)
            b)  compare strings/bools (except) using SequenceMatcher and produce ratio
            c)  Summate all errors keeping track of error types
    
    (5) Continue through files and report levels of differences/similarities 
            Store results in the object 
            Use to determine how much simulation results have changed between runs
'''


comp.closeFiles()


