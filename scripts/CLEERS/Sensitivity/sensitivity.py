''' Python script to perform a sensitivity analysis on a simple model
    by taking finite partial derivatives of that model given a set of
    parameters. '''

import math
import sys

# ----------- Start: Definition of Class object for Sensitivity --------------------
class Sensitivity(object):
    #Constructor for the object
    def __init__(self, func, func_params, func_consts):
        self.errors = False
        # Check to make sure the arguments are dictionaries
        if type(func_params) is not dict:
            print("Error! func_params must be a dictionary object!")
            self.errors = True
            return
        if type(func_consts) is not dict:
            print("Error! func_consts must be a dictionary object!")
            self.errors = True
            return
        # Check to make sure func is a Function
        if callable(func) == False:
            print("Error! func must be a callable function object!")
            self.errors = True
            return

        #NOTE: func_params and func_consts should be map objects whose keys are the names of parameters or constants of interest
        self.func = func                    #A function that produces a single output given a set of parameters and constants
        self.func_params = func_params      #A set of parameters that the sensitivity analysis will be performed on
        self.func_consts = func_consts      #A set of other constants or information the model needs to use
        self.partials = {}                  #Computed set of partial derivatives
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
            message = "Function parameters: " + str(self.func_params) + "\n"
            message += "Partial derivatives: " + str(self.partials) + "\n"
            message += "Sorted Sensitivity: " + str(self.sorted_param_sensitivity) + "\n"
            message += "Most sensitive parameter: " + self.highest_sensitivity + "\n"
            message += "Least sensitive parameter: " + self.lowest_sensitivity + "\n"
            return message

    #Function to call the users function with their parameters and constants
    def eval_func(self):
        if self.errors == False:
            return self.func(self.func_params,self.func_consts)
        else:
            print("Errors exist in the object!")
            return 0.0

    #Function to compute a partial derivative for the given parameter
    def calc_partial(self, param_name):
        if self.errors == True:
            print("Errors exist in the object!")
            return
        if param_name not in self.partials.keys():
            print("Error! Parameter name not registered!")
            self.errors = True
            return
        old_value = self.func_params[param_name]
        new_value = old_value + math.sqrt(sys.float_info.epsilon)
        old_func = self.eval_func()
        self.func_params[param_name] = new_value
        new_func = self.eval_func()
        self.func_params[param_name] = old_value
        self.partials[param_name] = (new_func-old_func)/(new_value-old_value)

    #Function to compute all partials
    def compute_partials(self):
        if self.errors == True:
            print("Errors exist in the object!")
            return
        old_func = self.eval_func()
        i = 0
        smallest = 0
        largest = 0
        for param_name in self.partials:
            #NOTE: For efficiency, we will neglect the constant parameter checking
            #if param_name not in self.partials.keys():
            #    print("Error! Parameter name not registered!")
            #    self.errors = True
            #    return
            old_value = self.func_params[param_name]
            new_value = old_value + math.sqrt(sys.float_info.epsilon)
            self.func_params[param_name] = new_value
            new_func = self.eval_func()
            self.func_params[param_name] = old_value
            self.partials[param_name] = (new_func-old_func)/(new_value-old_value)
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
        self.sorted_param_sensitivity = {k: v for k, v in sorted(self.partials.items(), key=lambda item: item[1], reverse=True)}

# ----------- End: Definition of Class object for Sensitivity --------------------


# ---------- Testing -------------
def test_func(params, consts):
    return params["C"]*params["C"]*params["A"] + consts["B"]


test_params = {}
test_params["A"] = 2
test_params["C"] = 3

test_consts = {}
test_consts["B"] = 5

test_obj = Sensitivity(test_func,test_params,test_consts)
test_obj.compute_partials()

print(test_obj)
