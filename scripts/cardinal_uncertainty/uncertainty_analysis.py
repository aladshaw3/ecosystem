## Python script to perform uncertainty analysis for cardinal ##
## Run python scripts using Python 3.5 or newer ##

''' Uncertainty script:
    ----------------
    Object-Oriented approach to running an uncertainty analysis on the
    simulation results from the CARDINAL module of ecosystem. Uses various
    control file structures to change simulation behavior based on modifying
    input parameters. Output results are compared to assess how much the 
    simulation changes as a result of the input changes. Uses the knapsack
    optimizer to find sub-set of parameters that result in highest levels
    of uncertainty (or lowest levels)
    
    Author:     Austin Ladshaw
    Date:       05/15/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.
    '''
import sys, os
sys.path.insert(0, '../output-comparison')
sys.path.insert(0, '../knapsack-optimizer')
import atmosphere_reader as ar
import control_reader as cr
import output_comp as oc
import knapsack as ns
from ctypes import *
lib = cdll.LoadLibrary('../../libeco.so')

class CardinalParam(object):
    tag = 0
    # Inputs:   Name of the parameter, value to change parameter by, whether or not value is a percentage
    def __init__(self, name, change, percentage):
        #Names are specific and must mimic Enums from Atmosphere and Control
        # Addition name includes FPY for Fission Product Yields
        # Note: FPY has its own values that it changes by (-1, 0, 1)
        self.param_name = name
        self.param_change = change
        self.usePercent = percentage
        self.id = CardinalParam.tag
        CardinalParam.tag += 1

    def __str__(self):
        return self.param_name + " : " + str(self.param_change) + " : " + str(self.usePercent)

# Constraint Function
def subset_size(list, size, args):
    # args is a Tuple and the first value is the size limit
    if size > args[0]:
        return False
    else:
        return True
#End Constraint Function

# Objective Function
## args is a Tuple whose values are as follows:
#           args[0] = cloud_rise file, args[1] = cloud_growth file, and args[2] = nuclides file
#           args[3] = True (for maximizing) and False (for minimizing)
#           args[4] = Control_Input file, args[5] = Atmosphere_Input file, and args[6] = database path
def output_change(list, size, args):
    numeric_difference = 0.0
    #Compute numeric difference by running a simulation and comparing output
    #   Setup function
    card_func = lib.cardinal_simulation
    card_func.restype = c_int
    #                    sim_file,  atm_file, nuc_path, unc_opt, rise_out, growth_out, nuc_out
    card_func.argtypes = [c_char_p, c_char_p, c_char_p, c_int, c_char_p, c_char_p, c_char_p]
    
    # Read input files to prepare edits
    control = cr.ControlFile()
    control.readFile(args[4])
    
    atm = ar.Atmosphere()
    atm.readFile(args[5])
    
    # Change the input files and prepare output files
    c_file = args[4].split(".")[0]
    a_file = args[5].split(".")[0]
    fpy_opt = 0
    rise_out = args[0].split(".")[0]
    growth_out = args[1].split(".")[0]
    nuc_out = args[2].split(".")[0]
    
    for param in list:
        c_file += "_" + str(param.id)
        a_file += "_" + str(param.id)
        rise_out += "_" + str(param.id)
        growth_out += "_" + str(param.id)
        nuc_out += "_" + str(param.id)
        if param.param_name == "Temperature" and param.usePercent == True:
            atm.editValue_PercentChange(ar.AtmParam.Temperature, param.param_change)
        if param.param_name == "Pressure" and param.usePercent == True:
            atm.editValue_PercentChange(ar.AtmParam.Pressure, param.param_change)
        if param.param_name == "Humidity" and param.usePercent == True:
            atm.editValue_PercentChange(ar.AtmParam.Humidity, param.param_change)
        if param.param_name == "Yield" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.Yield, param.param_change)
        if param.param_name == "Height" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.Height, param.param_change)
        if param.param_name == "Level" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.Level, param.param_change)
        if param.param_name == "Wind" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.Wind, param.param_change)
        if param.param_name == "Mass" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.Mass, param.param_change)
        if param.param_name == "FissExtent" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.FissExtent, param.param_change)
        if param.param_name == "FissYield" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.FissYield, param.param_change)
        if param.param_name == "PartDist" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.PartDist, param.param_change)
        if param.param_name == "Casing" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.Casing, param.param_change)
        if param.param_name == "FracMod" and param.usePercent == True:
            control.editValue_PercentChange(cr.ControlParam.FracMod, param.param_change)
        if param.param_name == "FPY":
            fpy_opt = param.param_change
        if param.param_name == "Temperature" and param.usePercent == False:
            atm.editValue_LinearChange(ar.AtmParam.Temperature, param.param_change)
        if param.param_name == "Pressure" and param.usePercent == False:
            atm.editValue_LinearChange(ar.AtmParam.Pressure, param.param_change)
        if param.param_name == "Humidity" and param.usePercent == False:
            atm.editValue_LinearChange(ar.AtmParam.Humidity, param.param_change)
        if param.param_name == "Yield" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.Yield, param.param_change)
        if param.param_name == "Height" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.Height, param.param_change)
        if param.param_name == "Level" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.Level, param.param_change)
        if param.param_name == "Wind" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.Wind, param.param_change)
        if param.param_name == "Mass" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.Mass, param.param_change)
        if param.param_name == "FissExtent" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.FissExtent, param.param_change)
        if param.param_name == "FissYield" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.FissYield, param.param_change)
        if param.param_name == "PartDist" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.PartDist, param.param_change)
        if param.param_name == "Casing" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.Casing, param.param_change)
        if param.param_name == "FracMod" and param.usePercent == False:
            control.editValue_LinearChange(cr.ControlParam.FracMod, param.param_change)

    c_file += ".txt"
    a_file += ".txt"
    rise_out += ".txt"
    growth_out += ".txt"
    nuc_out += ".txt"

    # Only print to file if file doesn't exist
    if c_file != args[4]:
        control.print2file(c_file)
    if a_file != args[5]:
        atm.print2file(a_file)

    # Call the solver
    success = success = card_func(c_file.encode(), a_file.encode(), args[6].encode(), fpy_opt, rise_out.encode(), growth_out.encode(), nuc_out.encode())
    if success != 0:
        print("Error has occurred...")
        return 0.0
    
    # Make comparisons and record results
    comp_rise = oc.FileCompare("output/"+args[0],"output/"+rise_out)
    comp_growth = oc.FileCompare("output/"+args[1],"output/"+growth_out)
    comp_nuc = oc.FileCompare("output/"+args[2],"output/"+nuc_out)
    numeric_difference += comp_rise.num_diff + comp_growth.num_diff + comp_nuc.num_diff
    
    if args[3] == True:
        return numeric_difference
    else:
        return -numeric_difference
#End Objective Function

## ----------- Run Script -------------- ##

# Choose what list of params to use and how to change them
param_list = []
param_list.append(CardinalParam("Temperature",10.0,True))
param_list.append(CardinalParam("Pressure",10.0,True))
param_list.append(CardinalParam("Humidity",10.0,True))
param_list.append(CardinalParam("Yield",10.0,True))
param_list.append(CardinalParam("Height",10.0,True))
param_list.append(CardinalParam("Level",10.0,True))
param_list.append(CardinalParam("Wind",10.0,True))
param_list.append(CardinalParam("Mass",10.0,True))
param_list.append(CardinalParam("FissExtent",10.0,True))
param_list.append(CardinalParam("FissYield",10.0,True))
param_list.append(CardinalParam("PartDist",10.0,True))
param_list.append(CardinalParam("Casing",10.0,True))
param_list.append(CardinalParam("FracMod",10.0,True))
param_list.append(CardinalParam("FPY",1.0,True))

test = []
test.append(CardinalParam("Height",10.0,False))
test.append(CardinalParam("Humidity",10.0,True))

obj_args = ("py_cloud_rise_base.txt","py_cloud_growth_base.txt","py_nuclide_base.txt",True,"1979-Test-Case.txt","DefaultAtmosphere.txt","../../database/")
size_limit = 1
#value = output_change(test, len(test), obj_args)

optimizer = ns.ZeroOneKnapsack()
optimizer.register_objective_func(output_change, *obj_args)
optimizer.register_constraints(subset_size, size_limit)
optimizer.exhaustive_search(False)
#optimizer.eval_objective_func(test)
#optimizer.eval_constraints(test)
(val, new_list, status) = optimizer.Optimize(test) #Appears to be working, needs more tests

print((val, new_list, status))
#print(value)
for obj in new_list:
    print(obj)


