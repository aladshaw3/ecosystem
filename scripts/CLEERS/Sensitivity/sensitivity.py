'''
    \file sensitivity.py
    \brief Simple Sensitivity Analysis
    \details Python script to perform a sensitivity analysis on a simple model
                by taking finite partial derivatives of that model given a set of
                parameters.
    \author Austin Ladshaw
    \date 02/12/2020
    \copyright This software was designed and built at the Oak Ridge National
                    Laboratory (ORNL) National Transportation Research Center
                    (NTRC) by Austin Ladshaw for research in the catalytic
                    reduction of NOx. Copyright (c) 2020, all rights reserved.
'''

import math
import sys

# ----------- Start: Definition of Class object for Sensitivity --------------------
class Sensitivity(object):
    #Constructor for the object
    def __init__(self, func, func_params, func_conds):
        self.errors = False
        # Check to make sure the arguments are dictionaries
        if type(func_params) is not dict:
            print("Error! func_params must be a dictionary object!")
            self.errors = True
            return
        if type(func_conds) is not dict:
            print("Error! func_conds must be a dictionary object!")
            self.errors = True
            return
        # Check to make sure func is a Function
        if callable(func) == False:
            print("Error! func must be a callable function object!")
            self.errors = True
            return

        #NOTE: func_params and func_conds should be map objects whose keys are the names of parameters or conditions of interest
        self.func = func                    #A function that produces a single output given a set of parameters and conditions
        self.func_params = func_params      #A set of parameters that the sensitivity analysis will be performed on
        self.func_conds = func_conds        #A set of other conditions or information the model needs to use
        self.partials = {}                  #Computed set of partial derivatives or percent changes
        self.relative_sensitivity = False   #When set to True, the partials are computed based on a percent change to the variable
        self.percent_change = 0.0           #Percent change to apply to parameters when relative_sensitivity = True
        self.partials_computed = False
        self.lowest_sensitivity = ""
        self.highest_sensitivity = ""
        self.sorted_param_sensitivity = {}  #Stores a sorted map of most to least sensitve parameters

        #Initialize the list of partial derivatives and look for errors
        for item in self.func_params:
            if item in self.partials.keys():
                print("Error! Duplicate parameter name!")
                self.errors = True
                return
            #Check the func_params at the item to make sure it is a number
            if type(self.func_params[item]) is not int and type(self.func_params[item]) is not float:
                print("Error! Parameters must be numbers!")
                self.errors = True
                return
            self.partials[item] = 0.0

    #Function to print the results of the analysis to the console
    def __str__(self):
        if self.partials_computed == False:
            return "Analysis not yet performed! Call function 'compute_partials' first...\n"
        else:
            if self.relative_sensitivity == False:
                message = " ------------ Sensitivity Analysis Performed using Standard Finite Differences ------------ \n"
            else:
                message = " ------------ Sensitivity Analysis Performed using " + str(self.percent_change)+ "% Change in Parameters ------------ \n"
            message += "Function parameters:\n"
            for item in self.func_params:
                message += "\t" + str(item) + "\t" + str(self.func_params[item]) + "\n"
            message += "Function conditions:\n"
            for item in self.func_conds:
                message += "\t" + str(item) + "\t" + str(self.func_conds[item]) + "\n"
            message += "Partial derivatives:\n"
            for item in self.partials:
                message += "\t" + str(item) + "\t" + str(self.partials[item]) + "\n"
            message += "Sorted Sensitivity:\n"
            for item in self.sorted_param_sensitivity:
                message += "\t" + str(item) + "\t" + str(self.sorted_param_sensitivity[item]) + "\n"
            message += "Most sensitive parameter: " + self.highest_sensitivity + "\n"
            message += "Least sensitive parameter: " + self.lowest_sensitivity + "\n"
            return message

    #Function to call the users function with their parameters and conditions
    def eval_func(self):
        if self.errors == False:
            return self.func(self.func_params,self.func_conds)
        else:
            print("Errors exist in the object!")
            return 0.0

    #Function to compute all partials
    def compute_partials(self, relative = False, per = 1):
        if self.errors == True:
            print("Errors exist in the object!")
            return
        old_func = self.eval_func()
        self.relative_sensitivity = relative
        self.percent_change = per
        i = 0
        smallest = 0
        largest = 0
        for param_name in self.partials:
            old_value = self.func_params[param_name]
            if relative == False:
                new_value = old_value + math.sqrt(sys.float_info.epsilon)
            else:
                new_value = (1 + (per/100))*old_value
            self.func_params[param_name] = new_value
            new_func = self.eval_func()
            self.func_params[param_name] = old_value
            try:
                if relative == False:
                    self.partials[param_name] = (new_func-old_func)/(new_value-old_value)
                else:
                    try:
                        self.partials[param_name] = (new_func-old_func)/old_func*100
                    except:
                        try:
                            self.partials[param_name] = (new_func-old_func)/new_func*100
                        except:
                            self.partials[param_name] = 0
            except:
                new_value = 1.01*old_value
                self.func_params[param_name] = new_value
                new_func = self.eval_func()
                self.func_params[param_name] = old_value
                if relaive == False:
                    self.partials[param_name] = (new_func-old_func)/(new_value-old_value)
                else:
                    try:
                        self.partials[param_name] = (new_func-old_func)/old_func*100
                    except:
                        try:
                            self.partials[param_name] = (new_func-old_func)/new_func*100
                        except:
                            self.partials[param_name] = 0
            if i == 0:
                smallest = abs(self.partials[param_name])
                self.lowest_sensitivity = param_name
                largest = abs(self.partials[param_name])
                self.highest_sensitivity = param_name
            else:
                if abs(self.partials[param_name]) > largest:
                    largest = abs(self.partials[param_name])
                    self.highest_sensitivity = param_name
                if abs(self.partials[param_name]) < smallest:
                    smallest = abs(self.partials[param_name])
                    self.lowest_sensitivity = param_name
            i += 1
        self.partials_computed = True
        #Sort the parameters from most to least sensitive
        self.sorted_param_sensitivity = {k: v for k, v in sorted(self.partials.items(), key=lambda item: abs(item[1]), reverse=True)}

# ----------- End: Definition of Class object for Sensitivity --------------------

#Helper function to iterate through all permutations of conditions
#       Function will return True after the last permutation has been made
# cond_value = current list of values of conditions that needs updating
# cond_limit_lower = list of the upper limits of the conditions
# cond_limit_upper = list of the upper limits of the conditions
def update_cond(cond_value, cond_limit_lower, cond_limit_upper):
    complete = False
    i = 0
    for value in cond_value:
        dist = cond_limit_upper[i] - cond_limit_lower[i]
        if abs(dist) <= 10*math.sqrt(sys.float_info.epsilon):
            #print("Skipped this analysis... Spacing between bounds is too small\n")
            cond_value[i] = cond_limit_lower[i]
            if i == len(cond_limit_upper)-1:
                complete = True
        else:
            #base is calculated based on the number of conditions to vary
            #       min_base = 2 (can't do less than this) -->  max_base = 10 (don't need more than this)
            base = int(math.pow(10,math.log10(100000)/len(cond_value)))
            if base > 10:
                base = 10
            if base < 2:
                base = 2
            update = dist/base
            if value+update <= cond_limit_upper[i]:
                cond_value[i] = value+update
                break
            else:
                cond_value[i] = cond_limit_lower[i]
                if i == len(cond_limit_upper)-1:
                    complete = True
        i+=1
    return complete


# ----------- Start: Definition of Class object for SensitivitySweep --------------------
''' The SensitivitySweep object is an object that uses the Sensitivity object to
    calculation partials or changes in a model with changes in parameters, but
    also repeats this process for a ranged of conditions to produce sensitivity
    matrices that are output to a file. This is necessary for complex models as
    it is possible that a model will not be sensitive to a certain parameter under
    certain conditions, but becomes more sensitive as the conditions change. '''
class SensitivitySweep(object):
    #Constructor
    def __init__(self, func, func_params, func_conds_tuples):
        '''
        NOTE: func_conds_tuples must be a dictionary whose keys are the simulation/model
                condtions and whose values are tuples representing the lower and upper
                bounds of the conditions, respectively.

                e.g.,  func_conds_tuples["Temp"] = (273, 373)
                            a condition for temperature that spans 100 degrees
        '''
        self.errors = False
        start_conditions = {}
        if type(func_conds_tuples) is not dict:
            print("Error! func_conds_tuples must be a dictionary object!")
            self.errors = True
            return
        for item in func_conds_tuples:
            if type(func_conds_tuples[item]) is not tuple:
                print("Error! Items in func_conds_tuples must be tuples!")
                self.errors = True
                return
            start_conditions[item] = func_conds_tuples[item][0]
        self.cond_tuples = func_conds_tuples
        self.sens_obj = Sensitivity(func, func_params, start_conditions)

    #Run the Sensitivity Sweep Analysis and print results to a file
    #       User may also specify whether or not to use relative parameter changes and the percent to change
    def run_sweep(self, sensitivity_file_name = "SensitivitySweepAnalysis.txt", relative = False, per = 1):
        file = open(sensitivity_file_name,'w')
        #Write out header information
        file.write("-------------- Sensitivity Sweep Results ---------------\n")
        if relative == False:
            file.write("\tAnalysis performed using Finite Differences\n")
        else:
            file.write("\tAnalysis performed using " + str(per) + "% Change in baseline parameter values\n")
        file.write("\nSystem Parameters:\n")
        for item in self.sens_obj.func_params:
            file.write(str(item) + "\t" + str(self.sens_obj.func_params[item]) + "\n")
        file.write("\nSystem Baseline Conditions:\n")
        for item in self.sens_obj.func_conds:
            file.write(str(item) + "\t" + str(self.sens_obj.func_conds[item]) + "\n")

        #Begin looping through conditions, running sensitivity analyses, and printing results in a matrix
        i = 0
        for cond in self.sens_obj.func_conds:
            file.write("\n")
            file.write("Sensitivity Matrix for " + str(cond) + "\n-----------------------------\n")
            file.write(str(cond))
            for param in self.sens_obj.func_params:
                if relative == False:
                    file.write("\tdf/d("+str(param)+")")
                else:
                    file.write("\tdf(%)@d("+str(param)+")")
            file.write("\n")

            og_cond = self.sens_obj.func_conds[cond]
            #Check to see what type of data the condition is (should be number or boolean)
            if type(self.sens_obj.func_conds[cond]) is int or type(self.sens_obj.func_conds[cond]) is float:
                dist = self.cond_tuples[cond][-1] - self.cond_tuples[cond][0]
                if abs(dist) <= 10*math.sqrt(sys.float_info.epsilon):
                    file.write("Skipped this analysis... Spacing between bounds is too small\n")
                else:
                    dx = dist/10
                    for n in range(0,11):
                        update = self.cond_tuples[cond][0] + n*dx
                        self.sens_obj.func_conds[cond] = update
                        self.sens_obj.compute_partials(relative, per)
                        file.write(str(self.sens_obj.func_conds[cond]))
                        for param in self.sens_obj.func_params:
                            file.write("\t"+str(self.sens_obj.partials[param]))
                        file.write("\n")
            elif type(self.sens_obj.func_conds[cond]) is boolean:
                #Run simulation as is
                self.sens_obj.compute_partials(relative, per)
                file.write(str(self.sens_obj.func_conds[cond]))
                for param in self.sens_obj.func_params:
                    file.write("\t"+str(self.sens_obj.partials[param]))
                file.write("\n")
                #Run simulation changing the boolean
                if self.sens_obj.func_conds[cond] == True:
                    self.sens_obj.func_conds[cond] = False
                else:
                    self.sens_obj.func_conds[cond] = True
                self.sens_obj.compute_partials(relative, per)
                file.write(str(self.sens_obj.func_conds[cond]))
                for param in self.sens_obj.func_params:
                    file.write("\t"+str(self.sens_obj.partials[param]))
                file.write("\n")
            else:
                print("Error! Invalid condition type in self.sens_obj.func_conds")
                file.write("Skipped due to type error\n")

            #Recover original condition before proceeding to next
            self.sens_obj.func_conds[cond] = og_cond
            i += 1
        #End func_conds loop
        file.close()

    #Function to perform an Exhaustive Sensitivity Analysis
    def run_exhaustive_sweep(self, sensitivity_file_name = "ExhaustiveSensitivitySweepAnalysis.dat", relative = False, per = 1):
        #Enfore a .dat file extension. This is used to flag the file as very large so that
        #   our git tracking ignores that file. May be unnecessary, but better safe than sorry.
        set = sensitivity_file_name.split('.')
        sensitivity_file_name = set[0] + ".dat"
        file = open(sensitivity_file_name,'w')
        #Write out header information
        file.write("-------------- Exhaustive Sensitivity Sweep Results ---------------\n")
        if relative == False:
            file.write("\tAnalysis performed using Finite Differences\n")
        else:
            file.write("\tAnalysis performed using " + str(per) + "% Change in baseline parameter values\n")
        file.write("\nSystem Parameters:\n")
        for item in self.sens_obj.func_params:
            file.write(str(item) + "\t" + str(self.sens_obj.func_params[item]) + "\n")

        cond_list = []
        cond_value = []
        cond_limit_upper = []
        cond_limit_lower = []
        for cond in self.sens_obj.func_conds:
            if type(self.sens_obj.func_conds[cond]) is not int and type(self.sens_obj.func_conds[cond]) is not float:
                print("Non-numeric conditions are not currently supported in Exhaustive Search method...")
                print("Variations of the " + str(cond) + " condition will be skipped...")
            else:
                cond_list.append(cond)
                cond_value.append(self.cond_tuples[cond][0])
                cond_limit_lower.append(self.cond_tuples[cond][0])
                cond_limit_upper.append(self.cond_tuples[cond][-1])

        #NOTE: This limitation keeps the total permutation count around 100,000.
        #       Beyond this point, the analysis is too complex for this script
        if len(cond_list) > 15:
            print("Too many conditions to vary! Number of permutations exceeds limits!")
            print("Try using run_sweep() instead or reduce the number of conditions to below 15...")
            self.errors = True
            file.write("\nToo many conditions to vary! Number of permutations exceeds limits!\n")
            file.write("Try using run_sweep() instead or reduce the number of conditions to below 15...\n")
            file.close()
            return
        #Print out a header line for conditions and parameters
        file.write("\n")
        n = 0
        for cond in self.sens_obj.func_conds:
            file.write("\tcond_"+str(n))
            n+=1
        n = 0
        for param in self.sens_obj.func_params:
            file.write("\tpartial_"+str(n))
            n+=1
        file.write("\nPermutation")
        for cond in self.sens_obj.func_conds:
            file.write("\t"+str(cond))
        for param in self.sens_obj.func_params:
            if relative == False:
                file.write("\tdf/d("+str(param)+")")
            else:
                file.write("\tdf(%)@d("+str(param)+")")
        file.write("\n")

        #Changing conditions by going through all permuntations
        complete = False
        perm = 0
        while complete == False:
            #Run current case
            i = 0
            for cond in cond_list:
                self.sens_obj.func_conds[cond] = cond_value[i]
                i += 1
            self.sens_obj.compute_partials(relative, per)
            file.write(str(perm))
            for cond in self.sens_obj.func_conds:
                file.write("\t"+str(self.sens_obj.func_conds[cond]))
            for param in self.sens_obj.func_params:
                file.write("\t"+str(self.sens_obj.partials[param]))
            file.write("\n")
            #Update Values
            complete = update_cond(cond_value, cond_limit_lower, cond_limit_upper)
            perm+=1
        file.close()


# ----------- End: Definition of Class object for SensitivitySweep --------------------

# ---------- Testing -------------


def test_func(params, conds):
    return params["B"]*params["B"]*params["A"] + conds["X"]*params["A"] + conds["Y"]+conds["Z"]*params["B"]


test_params = {}
test_params["A"] = 20
test_params["B"] = 3

test_conds = {}
test_conds["X"] = 0
test_conds["Y"] = 0
test_conds["Z"] = 0

test_tuples = {}

for item in test_conds:
    test_tuples[item] = (test_conds[item], test_conds[item]+2)
#test_tuples["Y"] = (0,0)
#print(test_tuples)

test_obj = SensitivitySweep(test_func,test_params, test_tuples)
test_obj.run_exhaustive_sweep(sensitivity_file_name = "ExhaustiveSensitivitySweepAnalysis.dat", relative = True, per = 10)
