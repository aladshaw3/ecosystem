## Python script to test running eco functions ##
## Run python scripts using Python 3.5 or newer ##
import numpy as np

# Add path to optimiser to python
import sys, os
sys.path.insert(0, '../knapsack-optimizer')

import knapsack as ns

from ctypes import *
lib = cdll.LoadLibrary('../../libeco.so')
print(lib)
#success = lib.blah()

func = lib.obj_func
func.restype = c_double
func.argtypes = [POINTER(c_double), c_int]

data = POINTER(c_double)(c_double())
data[0] = 2.0
data[1] = 3.0
len = 2
args = POINTER(c_double)(c_double())
args[0] = 2.0
args[1] = 1.0
set = (len, args)

print(func(data,*set))

data2 = []

# Check to see if variable is a python list
if isinstance(data,(list,)):
    print('This is a python list')
else:
    print('This is NOT a python list')
    #Convert c_double pointer to a python list
    data2 = np.fromiter(data, dtype=np.float, count=len).tolist()

if isinstance(data2,(list,)):
    print('This is NOW a python list')
    print(data2)
else:
    print('You screwed up')
    print(data2)

# NOTE: Must use lists in our python script for the knapsack problem
# HOWEVER: Cannot use lists as arguments for calling C-style functions
# SOLUTION: Create a python function that accepts the list and converts it
#           then sends the converted information to the C-style function

def py_func(list, size, args):
    # Input is a python list, size, and a tuple of other args
    
    #Use syntax below to create a double array of specific size for C
    data_set = (c_double * size)()
    arg_set = (c_double * size)()
    py_set = (size, arg_set)
    
    # NOTE: You are now able to use any python types within both list
    #       and args, so long as you convert into the expected types
    #       for the C-style function call
    i = 0
    for val in list:
        data_set[i] = val
        i += 1
    i = 0
    if size > 0:
        for val in args[0]:
            if i < size:
                arg_set[i] = val
            i += 1
    return func(data_set,*py_set)


print(py_func([2,3], 2, ([2,1],)))   #Works!!!
print(py_func([3], 1, ([1],)))   #Works!!!
print(py_func([3,10,4,5], 4, ([1,2,3,4],)))   #Works!!!

def pick_one(list, size, args):
    if size <= 1:
        return True
    else:
        return False

# Test knapsack with C function
list = [1,2,5,3,10,12,4,6]
args = [0,1,2,3,4,5,6,7]
prob = ns.ZeroOneKnapsack()
prob.register_objective_func(py_func, args)
prob.register_constraints(pick_one)
prob.exhaustive_search(False)
(val, new_list, status) = prob.Optimize(list)

print('\nAfter Optimization: Test 01 Results')
print(val)
print(status)
print('\nTaken')
for obj in new_list:
    print(obj)

# Run command line programs in current directory
os.system("ls -a")
#os.system("cp test.py test_copy.py") # WORKS
# MUST make an output directory (CARDNINAL EXPECTS IT TO EXIST!!!)
os.system("mkdir output")
os.system("../../eco -e cardinal -i ../../input_files/CARDINAL/1979-Test-Case.txt ../../input_files/CARDINAL/DefaultAtmosphere.txt ../../database/")
