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
        string = "Name =\t" + self.param_name + "\t"
        if self.usePercent == True:
            string += "Percent Change =\t"
        else:
            string += "Value Change =\t"
        string += str(self.param_change) + "\tID Tag =\t" + str(self.id)
        return string

# Constraint Function
def subset_size(list, size, args):
    # args is a Tuple and the first value is the size limit (args[0] = upper limit && args[1] = lower limit)
    if size > args[0] or size < args[1]:
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

    # Skip solver if files are same
    if c_file == args[4] and a_file == args[5]:
        return 0.0

    # Modify the path to the input files to keep directory clean
    new_dir = args[4].split(".")[0]
    if os.path.exists(new_dir) == False:
        os.mkdir(new_dir)
    c_file = new_dir + "/" + c_file
    a_file = new_dir + "/" + a_file
    out_dir = "output/"+new_dir
    if os.path.exists(out_dir) == False:
        os.mkdir(out_dir)
    rise_out = new_dir + "/" + rise_out
    growth_out = new_dir + "/" + growth_out
    nuc_out = new_dir + "/" + nuc_out

    # Only print to file if file doesn't exist
    if c_file != args[4]:
        control.print2file(c_file)
    if a_file != args[5]:
        atm.print2file(a_file)
    
    # Call the solver
    success = success = card_func(c_file.encode(), a_file.encode(), args[6].encode(), int(fpy_opt), rise_out.encode(), growth_out.encode(), nuc_out.encode())
    if success != 0:
        print("Error has occurred...")
        return 0.0
    
    # Make comparisons and record results
    comp_rise = oc.FileCompare("output/"+args[0],"output/"+rise_out)
    comp_growth = oc.FileCompare("output/"+args[1],"output/"+growth_out)
    comp_nuc = oc.FileCompare("output/"+args[2],"output/"+nuc_out)
    numeric_difference += comp_rise.num_diff + comp_growth.num_diff + comp_nuc.num_diff
    
    if args[3] == True:
        return numeric_difference/3.0
    else:
        return -numeric_difference/3.0
#End Objective Function



## ----------- Run Script -------------- ##

### ----------------------- User May Interface Below -----------------------###
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
param_list.append(CardinalParam("FPY",1,True))

# Define the Control input file, Atmosphere file, and database path
control_file = "1979-Test-Case.txt"
atm_file = "DefaultAtmosphere.txt"
data_path = "../../database/"

# Pick whether or not to maximize or minimize uncertainty and the size of the parameter subset you will seek
Maximize = True
size_toplimit = 1
size_bottomlimit = 0

### ----------------------- User May Interface Above -----------------------###



### ----------------------- DO NOT CHANGE BELOW -----------------------###

# Automatically defining output files (Names come from control file name)
rise_out = control_file.split(".")[0] + "_cloud_rise_out.txt"
growth_out = control_file.split(".")[0] + "_cloud_growth_out.txt"
nuc_out = control_file.split(".")[0] + "_nuclides_out.txt"

# Call Solver routine for the base input files
card_func = lib.cardinal_simulation
card_func.restype = c_int
#                    sim_file,  atm_file, nuc_path, unc_opt, rise_out, growth_out, nuc_out
card_func.argtypes = [c_char_p, c_char_p, c_char_p, c_int, c_char_p, c_char_p, c_char_p]
success = card_func(control_file.encode(), atm_file.encode(), data_path.encode(), 0, rise_out.encode(), growth_out.encode(), nuc_out.encode())
if success != 0:
    print("Error in simulation setup! Cannot proceed...\n")
    exit()

#Args List  (cloud_rise_out base file, cloud_growth_out base file, nuclide_out base file, True=Maximize, Control input file, Atmosphere input file, path/to/database)
obj_args = (rise_out, growth_out, nuc_out, Maximize, control_file, atm_file, data_path)

# Setup the uncertainty analyzer and call the optimization routine
optimizer = ns.ZeroOneKnapsack()
optimizer.register_objective_func(output_change, *obj_args)
optimizer.register_constraints(subset_size, size_toplimit, size_bottomlimit)
optimizer.exhaustive_search(False)

# Open a results file to print some info
result_file = open(control_file.split(".")[0] + "_OptimizationResults.txt",'w')
result_file.write("Optimization Goal:\t")
if Maximize == True:
    result_file.write("Find params that Maximize uncertainty\n")
else:
    result_file.write("Find params that Minimize uncertainty\n")
result_file.write("Parameter set size restrictions:\t" + str(size_bottomlimit) + " <= size <= " + str(size_toplimit) + "\n")
result_file.write("\nSet of parameters, and their change values, being considered...\n")
result_file.write("---------------------------------------------------------------\n")
result_file.write("Name\tID\tChangeType\tChangeValue\n")
for param in param_list:
    result_file.write(param.param_name + "\t" + str(param.id) + "\t")
    if param.usePercent == True:
        result_file.write("Percent\t")
    else:
        result_file.write("Linear\t")
    result_file.write(str(param.param_change)+"\n")

# Call the optimizer
(val, new_list, status) = optimizer.Optimize(param_list)

# Record the results
result_file.write("\nOptimized set of parameters...\n")
result_file.write("--------------------------------\n")
result_file.write("Name\tID\tChangeType\tChangeValue\n")
suffix = ""
for param in new_list:
    result_file.write(param.param_name + "\t" + str(param.id) + "\t")
    suffix += "_" + str(param.id)
    if param.usePercent == True:
        result_file.write("Percent\t")
    else:
        result_file.write("Linear\t")
    result_file.write(str(param.param_change)+"\n")

result_file.write("\n\nControl Input file:\t" + control_file.split(".")[0] + "/" + control_file.split(".")[0] + suffix + ".txt")
result_file.write("\nAtmosphere Input file:\t" + control_file.split(".")[0] + "/" + atm_file.split(".")[0] + suffix + ".txt")
result_file.write("\n\nCloud Rise output file:\t" + "output/" + control_file.split(".")[0] + "/" + rise_out.split(".")[0] + suffix + ".txt")
result_file.write("\nCloud Growth output file:\t" + "output/" + control_file.split(".")[0] + "/" + growth_out.split(".")[0] + suffix + ".txt")
result_file.write("\nNuclides output file:\t" + "output/" + control_file.split(".")[0] + "/" + nuc_out.split(".")[0] + suffix + ".txt")
result_file.close()

print("\nAveraged Percent Error in Response =\t" + str(val))
print("Was converged? =\t" + str(status))
print("\nChoosen list of parameters\n---------------")
for obj in new_list:
    print(obj)
