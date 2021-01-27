## This framework file creates an object used to perform multiple parameter
# searches for the reaction_model_framework by randomly sampling starting
# parameter values within a global neighborhood of values, then performing
# a set of gradient searches (using 'ipopt') within local neighborhoods.

# Import the pyomo environment
from pyomo.environ import *

from reaction_model_framework import ReactionModel
import random
import os

# Import the IDAES environment (includes a custom 'ipopt' library)
# NOTE: This script piggy-backs off of the 'ipopt' solver factory available
#       through the IDAES framework. You must have installed IDAES through
#       the 'conda' environment/package manager.
#
#       For IDAES installation instructions:
#           https://idaes-pse.readthedocs.io/en/stable/advanced_user_guide/advanced_install/index.html#advanced-user-installation

# IDAES is not installed, then the script will search for any other available 'ipopt' library
if os.environ['CONDA_DEFAULT_ENV'] == "idaes":
    from idaes.core import *


# Object to perform Monte Carlo analysis
class MCSearch_ReactionModel(object):
    # Default Constructor
    def __init__(self):
        self.RxnModel = ReactionModel()   # Reaction model object
        self.ready = False
        self.best_result_id = 0           # index value containing the best results
        self.best_obj = -1                # Value of the best objective function
        self.best_params = {}             # Dictionary for the best parameters
        self.best_params_start_window = {} # Dictionary for the starting window
        self.best_params_start_frac = {}    # Dictionary for the starting fraction

    # Read the control file
    def read_yaml_simfile(self, file):
        self.RxnModel.read_yaml_simfile(file)

        if self.RxnModel.simulate_only == True:
            print("\nError! Key 'Simulate_Only' must be False to run parameter analysis!\n")
            return
        else:
            self.ready = True

        for rxn in self.RxnModel.instance.rxns:
            self.best_params[rxn] = {}
            self.best_params_start_window[rxn] = {}
            self.best_params_start_frac[rxn] = {}
            self.best_params[rxn]["A"] = value(self.RxnModel.instance.A[rxn])
            self.best_params[rxn]["E"] = value(self.RxnModel.instance.E[rxn])
            self.best_params[rxn]["B"] = value(self.RxnModel.instance.B[rxn])

            if self.RxnModel.instance.A[rxn].lb == None:
                self.RxnModel.instance.A[rxn].setlb(value(result.RxnModel.instance.A[rxn])*0.5)
            if self.RxnModel.instance.A[rxn].ub == None:
                self.RxnModel.instance.A[rxn].setub(value(result.RxnModel.instance.A[rxn])*2.0)

            if self.RxnModel.instance.E[rxn].lb == None:
                self.RxnModel.instance.E[rxn].setlb(value(result.RxnModel.instance.E[rxn])*0.5)
            if self.RxnModel.instance.E[rxn].ub == None:
                self.RxnModel.instance.E[rxn].setub(value(result.RxnModel.instance.E[rxn])*2.0)

            if self.RxnModel.instance.B[rxn].lb == None:
                self.RxnModel.instance.B[rxn].setlb(value(result.RxnModel.instance.B[rxn])*0.5)
            if self.RxnModel.instance.B[rxn].ub == None:
                self.RxnModel.instance.B[rxn].setub(value(result.RxnModel.instance.B[rxn])*2.0)

            self.best_params_start_window[rxn]["A"] = (self.RxnModel.instance.A[rxn].lb, self.RxnModel.instance.A[rxn].ub)
            self.best_params_start_window[rxn]["E"] = (self.RxnModel.instance.E[rxn].lb, self.RxnModel.instance.E[rxn].ub)
            self.best_params_start_window[rxn]["B"] = (self.RxnModel.instance.B[rxn].lb, self.RxnModel.instance.B[rxn].ub)

            try:
                self.best_params_start_frac[rxn]["A"] = (self.RxnModel.instance.A[rxn].ub/value(self.RxnModel.instance.A[rxn]))-1
            except:
                self.best_params_start_frac[rxn]["A"] = 0
            try:
                self.best_params_start_frac[rxn]["E"] = (self.RxnModel.instance.E[rxn].ub/value(self.RxnModel.instance.E[rxn]))-1
            except:
                self.best_params_start_frac[rxn]["E"] = 0
            try:
                self.best_params_start_frac[rxn]["B"] = (self.RxnModel.instance.B[rxn].ub/value(self.RxnModel.instance.B[rxn]))-1
            except:
                self.best_params_start_frac[rxn]["B"] = 0

    # Run a number of simulations
    def run_monte_carlo_analysis(self, iter):
        if self.ready == False:
            print("\nError! Model improperly setup! Call read_yaml_simfile() first...\n")
            return

        if not isinstance(iter, int):
            print("\nError! Model improperly setup! Iterations must be an int argument!\n")
            return

        # Run an initial simulation first to establish the current objective value
        self.RxnModel.run_model()
        self.best_obj = self.RxnModel.obj_value

        for rxn in self.RxnModel.instance.rxns:
            self.best_params[rxn]["A"] = value(self.RxnModel.instance.A[rxn])
            self.best_params[rxn]["E"] = value(self.RxnModel.instance.E[rxn])
            self.best_params[rxn]["B"] = value(self.RxnModel.instance.B[rxn])

        for i in range(iter):
            # Select random parameters in window
            for rxn in self.RxnModel.instance.rxns:
                Aval = random.uniform(self.best_params_start_window[rxn]["A"][0], self.best_params_start_window[rxn]["A"][1])
                Eval = random.uniform(self.best_params_start_window[rxn]["E"][0], self.best_params_start_window[rxn]["E"][1])
                Bval = random.uniform(self.best_params_start_window[rxn]["B"][0], self.best_params_start_window[rxn]["B"][1])

                #self.RxnModel.set_A(rxn, Aval, False, 0.5*self.best_params_start_frac[rxn]["A"], None)
                #self.RxnModel.set_E(rxn, Eval, False, 0.5*self.best_params_start_frac[rxn]["E"], None)
                #self.RxnModel.set_B(rxn, Bval, False, 0.5*self.best_params_start_frac[rxn]["B"], None)
                self.RxnModel.instance.A[rxn].set_value(Aval)
                self.RxnModel.instance.E[rxn].set_value(Eval)
                self.RxnModel.instance.B[rxn].set_value(Bval)

            # Run new simulation
            self.RxnModel.run_model()

            # If results have improved, then save them
            if self.RxnModel.obj_value < self.best_obj:
                self.best_obj = self.RxnModel.obj_value
                for rxn in self.RxnModel.instance.rxns:
                    self.best_params[rxn]["A"] = value(self.RxnModel.instance.A[rxn])
                    self.best_params[rxn]["E"] = value(self.RxnModel.instance.E[rxn])
                    self.best_params[rxn]["B"] = value(self.RxnModel.instance.B[rxn])
        #End range loop

        # After the loop is complete, perform one more optimization to finalize
        self.RxnModel.simulate_only = True
        for rxn in self.RxnModel.instance.rxns:
            self.RxnModel.set_A(rxn, self.best_params[rxn]["A"], False, 0.5*self.best_params_start_frac[rxn]["A"], None)
            self.RxnModel.set_E(rxn, self.best_params[rxn]["E"], False, 0.5*self.best_params_start_frac[rxn]["E"], None)
            self.RxnModel.set_B(rxn, self.best_params[rxn]["B"], False, 0.5*self.best_params_start_frac[rxn]["B"], None)

        self.RxnModel.run_model()
        self.RxnModel.print_to_file("MCSearch")
