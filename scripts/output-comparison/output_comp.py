## @package output_comp
#
#   @brief Python script to read output files and compare them
#
#   @details Script for reading in output produced from an executable and comparing it
#            to another output file. Useful for creating small unit tests for code validations
#            or making evaluations on how changes in parameters change predicted outcomes. Thus,
#            this script will be utilized for the purpose of performing sensitivity analyses.
#
#    @author     Austin Ladshaw
#
#    @date       05/09/2019
#
#    @copyright  This software was designed and built at the Georgia Institute
#                of Technology by Austin Ladshaw for research in the area of
#                radioactive particle decay and transport. Copyright (c) 2019,
#                all rights reserved.

import difflib
import math

## Class Object for Comparing Files
#
# This object will read in a pair of files line by line, parse the line into
#   a set of words or numbers, and compare all items in the files to find differences.
#   This is useful when you are doing unit testing (i.e., want to see if a simulation
#   produced the same output after updating some code) or can be used for a large
#   scale uncertainty or sensitivity analysis between a 'gold' standard simulation
#   and a 'test' case.
#
#   Method:
#
#    (1) Iterate through each line in each file \n
#            (iterate through smallest)
#
#    (2) Split each line into individual strings \n
#            Parse on spaces
#
#    (3) Iterate (nested) over each individual string in a line \n
#            (loop over the smallest string list)
#
#    (4) Check each string for types \n
#            a)  compare numbers/floats (try) and developed numeric different (sq. error) \n
#            b)  compare strings/bools (except) using SequenceMatcher and produce ratio \n
#            c)  Summate all errors keeping track of error types
#
#    (5) Continue through files and report levels of differences/similarities \n
#            Store results in the object \n
#            Use to determine how much simulation results have changed between runs
#
#   Notes:
#
#    read()   gives all text in file
#
#    read(n)  gives first n characters in file
#
#    readline()   gives first line in file      (**Each instance reads the next line)
#
#    readline(n)  gives nth line in the file
#
#    readlines()  gives list of all lines in file
#
#    for line in file:   loops over all lines in the file
#
#    string.split(ch)    produces list of sub-strings parsed by the ch character \n
#                        leaving ch blank defaults to splitting on blank spaces
#
#    use == to compare two strings
#
#        e.g.,  str1 == str2   -->   True if same, False if different
#
class FileCompare(object):
    ## Initialization constructor must take in strings for the gold file and test file
    #
    #   @param gold name of the file that you are comparing against
    #   @param test name of the file being tested for changes against the gold file
    def __init__(self, gold, test):
        ## Variable to keep track of numerical differences in the files
        #Values closest to zero represent smallest differences
        self.num_diff = 0.0
        ## Variable to keep track of string differences in the files
        #Values closest to zero represent smallest differences
        self.str_diff = 0.0
        ## Number of total numbers in both files
        self.total_num = 0
        ## Number of total words in both files
        self.total_word = 0
        self.hasBeenRead = False
        ## Number of lines in gold file
        self.lnum_gold = len(open(gold,"r").readlines())
        ## Number of lines in test file
        self.lnum_test = len(open(test,"r").readlines())
        self.gold_file = open(gold,"r")
        self.test_file = open(test,"r")
        self.computeErrors()

    ##Function to display results of the comparison to the console
    def __str__(self):
        if self.hasBeenRead == False:
            return "Must call computeErrors() first!"
        else:
            message = "\nComparison between " + self.gold_file.name + " and " + self.test_file.name
            message += "\n\nResults\n"
            message += "-------\n"
            message += "Avg Word Diff (%) =\t" + str(self.str_diff*100.0) + "\n"
            message += "Total Words       =\t" + str(self.total_word) + "\n"
            message += "Avg Num Diff (%)  =\t" + str(self.num_diff*100.0) + "\n"
            message += "Total Numbers     =\t" + str(self.total_num) + "\n"
            return message

    ## Function to read through each file line by line and compare each element
    #   and entry in each file with each other. If the two files are the same,
    #   then there will be no differences recorded between the two files.
    def computeErrors(self):
        if self.lnum_test > self.lnum_gold:
            short_file = self.gold_file
            long_file = self.test_file
        else:
            short_file = self.test_file
            long_file = self.gold_file

        #Iterate through the short file
        for s_line in short_file:
            l_line = long_file.readline()
            s_line_list = s_line.split()
            l_line_list = l_line.split()

            if len(s_line_list) > len(l_line_list):
                short_list = l_line_list
                long_list = s_line_list
            else:
                short_list = s_line_list
                long_list = l_line_list

            i = 0
            #Iteration through the short list
            for s in short_list:
                l = long_list[i]
                try:
                    if abs(float(s)) > abs(float(l)):
                        self.num_diff += abs(float(s) - float(l))/abs(float(s))
                    elif abs(float(s)) < abs(float(l)):
                        self.num_diff += abs(float(s) - float(l))/abs(float(l))
                    else:
                        self.num_diff += 0.0

                    self.total_num += 1
                except:
                    s = difflib.SequenceMatcher(None, s, l)
                    self.str_diff += s.ratio()
                    self.total_word += 1
                i += 1

            #List iteration complete

            #Finish the long list
            for lf in long_list[i:]:
                try:
                    self.num_diff += 1.0
                    self.total_num += 1
                except:
                    s = difflib.SequenceMatcher(None, "", lf)
                    self.str_diff += s.ratio()
                    self.total_word += 1
            #End Finishing
        #File Iteration complete

        #Finish reading the long file
        for remaining in long_file:
            for item in remaining.split():
                try:
                    self.num_diff += 1.0
                    self.total_num += 1
                except:
                    s = difflib.SequenceMatcher(None, "", item)
                    self.str_diff += s.ratio()
                    self.total_word += 1
        #End finishing

        #Final Editing Steps
        self.num_diff = (self.num_diff/self.total_num)
        self.str_diff = 1.0 - (self.str_diff/self.total_word)
        self.hasBeenRead = True
        self.closeFiles()

    ##Function to close the open files (called automatically by computeErrors())
    def closeFiles(self):
        self.gold_file.close()
        self.test_file.close()

### End Class Object ###

### Testing ###
#comp = FileCompare("gold_out.txt","test_out.txt")
#print(comp)
