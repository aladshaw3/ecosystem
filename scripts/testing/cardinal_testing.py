## Python script to test running eco functions ##
## Run python scripts using Python 3.5 or newer ##

# Add path to optimiser to python and import methods
import sys, os
sys.path.insert(0, '../knapsack-optimizer')
sys.path.insert(0, '../output-comparison')
import knapsack as ns
import output_comp as oc
from ctypes import *
lib = cdll.LoadLibrary('../../libeco.so')

#Setup function
card_func = lib.cardinal_simulation
card_func.restype = c_int
#                    sim_file,  atm_file, nuc_path, unc_opt, rise_out, growth_out, nuc_out
card_func.argtypes = [c_char_p, c_char_p, c_char_p, c_int, c_char_p, c_char_p, c_char_p]

#Setup function args
sim_file = "1979-Test-Case.txt"
atm_file = "DefaultAtmosphere.txt"
nuc_path = "../../database/"
unc_opt = 0
rise_out = "py_cloud_rise_base.txt"
growth_out = "py_cloud_growth_base.txt"
nuc_out = "py_nuclide_base.txt"

#Call routine
success = card_func(sim_file.encode(), atm_file.encode(), nuc_path.encode(), unc_opt, rise_out.encode(), growth_out.encode(), nuc_out.encode())

#Setup new function args (changed the burst height by 10 m)
sim_file = "1979-Test-Case.txt"
atm_file = "DefaultAtmosphere_15.txt"
nuc_path = "../../database/"
unc_opt = 0
rise_out = "py_cloud_rise_mod.txt"
growth_out = "py_cloud_growth_mod.txt"
nuc_out = "py_nuclide_mod.txt"

#Call routine again
success = card_func(sim_file.encode(), atm_file.encode(), nuc_path.encode(), unc_opt, rise_out.encode(), growth_out.encode(), nuc_out.encode())

#All output is in a sub-directory (/output/) from this directory
comp_rise = oc.FileCompare("output/py_cloud_rise_base.txt","output/py_cloud_rise_mod.txt")
print(comp_rise)

comp_growth = oc.FileCompare("output/py_cloud_growth_base.txt","output/py_cloud_growth_mod.txt")
print(comp_growth)

comp_nuc = oc.FileCompare("output/py_nuclide_base.txt","output/py_nuclide_mod.txt")
print(comp_nuc)
