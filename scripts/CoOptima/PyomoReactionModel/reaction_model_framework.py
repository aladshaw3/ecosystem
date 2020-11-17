## This is the framework set of objects and helper functions for the inhibition kinetic model
# contactor using Pyomo model objects and the 'ipopt' solver library
#
# NOTE: This script piggy-backs off of the 'ipopt' solver factory available
#       through the IDAES framework. You must have installed IDAES through
#       the 'conda' environment/package manager.
#
#       For IDAES installation instructions:
#           https://idaes-pse.readthedocs.io/en/stable/advanced_user_guide/advanced_install/index.html#advanced-user-installation

# Import the pyomo environment
from pyomo.environ import *

# Import the IDAES environment (includes SolverFactory() and 'ipopt')
from idaes.core import *

# Other import statements
import yaml
import os.path
from os import path

## Return string value for a valid weight method
def interpret_weight_method(name):
    option = name.lower()
    if option == "default":
        return "default"
    elif option == "equalizer":
        return "equalizer"
    elif option == "relative":
        return "relative"
    else:
        return "default"

##Calculation of rate constant
# E = activation energy in J/mol
# A = pre-exponential factor (units depend on reaction)
# B = power of temperature (usually = 0)
# T = temperature of system in K
def rate_const(A, B, E, T):
    return A*T**B*exp(-E/8.3145/T)

##Class object for the reaction model
class ReactionModel(object):
    #Default constructor
    # Pass a list of species names (Specs) for the mass balances
    def __init__(self):
        self.model = AbstractModel()
        self.model.tau = Param(domain=NonNegativeReals, initialize=1, mutable=True)
        self.model.eps = Param(domain=NonNegativeReals, initialize=0.5, mutable=True)

        self.built = False
        self.mb_set = False
        self.rxns_set = False
        self.temp_set = False
        self.simulate_only = True
        self.solved = False
        self.run_seriel = True   #Run simulations at each temperature in series
        self.data_set = False
        self.sub_dir = ""
        self.weight_method = "default"
        self.obj_value = 0
        self.total_var = 0

        self.temperature_set = {}  #Used to simulate at various temperatures
        self.Solution = {}  #Used to store solutions
        self.Data = {}  #Used to store data for optimization
        self.Weight_Factors = {} #Dictionary of weight factors for data fitting
        self.fix_dict = {}  #Dictionary of kinetic parameters to fix/lock

    # Read a data file
    #
    # Data format for data file
    #   Temp     \t Species[0]        \t Species[1]     \t .....
    #   T_val[0] \t Species_val[0][0] \t Species[1][0]  \t....
    #   .....
    def read_data_file(self, file):
        i = 0
        ordered_list = []
        for line in open(file, "r"):
            # Read in header that contains species names
            if i==0:
                j=0
                for item in line.split():
                    if j>0:
                        self.Data[item] = {}
                        self.Weight_Factors[item] = {}
                        ordered_list.append(item)
                    j+=1
            # Read in data for each species
            else:
                j=0
                for item in line.split():
                    if j==0:
                        self.temperature_set["T"+str(i-1)] = float(item)+273.15
                    else:
                        self.Data[ordered_list[j-1]]["T"+str(i-1)] = float(item)
                        self.Weight_Factors[ordered_list[j-1]]["T"+str(i-1)] = 1.0
                    j+=1
            i+=1

    # Read a yaml control file to setup the problem
    def read_yaml_simfile(self, file):
        dir_list = file.split("/")
        for item in dir_list:
            if "." not in item:
                self.sub_dir += item+"/"

        stream = open(file, "r")
        file_contents = yaml.safe_load_all(stream)
        for doc in file_contents:
            if isinstance(doc["Run_Seriel"], bool):
                self.run_seriel = doc["Run_Seriel"]
            else:
                print("\nError! Key 'Run_Seriel' must be a boolean!\n")
                return

            if isinstance(doc["Simulate_Only"], bool):
                self.simulate_only = doc["Simulate_Only"]
            else:
                print("\nError! Key 'Simulate_Only' must be a boolean!\n")
                return

            try:
                self.weight_method = interpret_weight_method( doc["Weight_Method"] )
            except:
                self.weight_method = "default"

            i=0
            temp_list = []
            data_spec = []
            if self.simulate_only == True:
                if type(doc["Temperature_Set"]) is not list:
                    print("\nError! Key 'Temperature_Set' must be a list!\n")
                    return

                for temp in doc["Temperature_Set"]:
                    self.temperature_set["T"+str(i)] = temp
                    temp_list.append("T"+str(i))
                    self.Solution["T"+str(i)] = {}
                    i+=1
            else:
                self.run_seriel = False
                if path.exists( self.sub_dir+doc["Data_File"] ):
                    self.read_data_file( self.sub_dir+doc["Data_File"] )
                    for temp in self.temperature_set:
                        temp_list.append(temp)
                        self.Solution[temp] = {}
                    for spec in self.Data:
                        data_spec.append(spec)
                else:
                    print("Error! Given file in Key 'Data_File' does not exist or the directory is wrong!")
                    return

            if self.run_seriel == True:
                self.add_temperatures( [temp_list[0]] )
            else:
                self.add_temperatures( temp_list )

            if type(doc["Chemical_Species"]) is not list:
                print("\nError! Key 'Chemical_Species' must be a list!\n")
                return
            self.add_species( doc["Chemical_Species"] )

            self.add_dataset( data_spec, doc["Chemical_Species"] )

            for species in doc["Chemical_Species"]:
                for temp in temp_list:
                    self.Solution[temp][species] = 0

            if type(doc["Rxn_Keys"]) is not list:
                print("\nError! Key 'Rxn_Keys' must be a list!\n")
                return
            self.add_reactions( doc["Rxn_Keys"] )

            self.build_instance()

            self.set_weight_factors()

            for species in self.Data:
                for temp in self.Data[species]:
                    self.set_data(species,temp,self.Data[species][temp])

            for temp in temp_list:
                self.set_temperature(temp,self.temperature_set[temp])
                if self.run_seriel == True:
                    break
            if isinstance(doc["Space_Velocity"], (int,float)):
                self.set_tau( doc["Space_Velocity"] )
            else:
                print("\nError! Key 'Space_Velocity' must be a int or float!\n")
                return

            if isinstance(doc["Void_Fraction"], (int,float)):
                self.set_eps( doc["Void_Fraction"] )
            else:
                print("\nError! Key 'Void_Fraction' must be a int or float!\n")
                return

            for species in doc["Inlet_Conc"]:
                self.set_initial(species,doc["Inlet_Conc"][species])
                self.set_inlet(species,doc["Inlet_Conc"][species])
            for rxn in doc["Reaction_Params"]:
                unbound_A = False
                unbound_E = False
                unbound_B = False
                delta_A = None
                delta_E = None
                delta_B = None
                frac_A = 0.2
                frac_E = 0.2
                frac_B = 0.2
                lock_A = False
                lock_B = False
                lock_E = False

                try:
                    if isinstance(doc["Reaction_Params"][rxn]["lock_A"], bool):
                        self.fix_dict[rxn]["A"] = doc["Reaction_Params"][rxn]["lock_A"]
                except:
                    self.fix_dict[rxn]["A"] = False
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["lock_B"], bool):
                        self.fix_dict[rxn]["B"] = doc["Reaction_Params"][rxn]["lock_B"]
                except:
                    self.fix_dict[rxn]["B"] = False
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["lock_E"], bool):
                        self.fix_dict[rxn]["E"] = doc["Reaction_Params"][rxn]["lock_E"]
                except:
                    self.fix_dict[rxn]["E"] = False
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["unbound_A"], bool):
                        unbound_A = doc["Reaction_Params"][rxn]["unbound_A"]
                except:
                    unbound_A = False
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["unbound_E"], bool):
                        unbound_E = doc["Reaction_Params"][rxn]["unbound_E"]
                except:
                    unbound_E = False
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["unbound_B"], bool):
                        unbound_B = doc["Reaction_Params"][rxn]["unbound_B"]
                except:
                    unbound_B = False

                try:
                    if isinstance(doc["Reaction_Params"][rxn]["delta_A"], (int,float)):
                        delta_A = doc["Reaction_Params"][rxn]["delta_A"]
                except:
                    delta_A = None
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["delta_E"], (int,float)):
                        delta_E = doc["Reaction_Params"][rxn]["delta_E"]
                except:
                    delta_E = None
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["delta_B"], (int,float)):
                        delta_B = doc["Reaction_Params"][rxn]["delta_B"]
                except:
                    delta_B = None

                try:
                    if isinstance(doc["Reaction_Params"][rxn]["frac_A"], (int,float)):
                        frac_A = doc["Reaction_Params"][rxn]["frac_A"]
                except:
                    frac_A = 0.2
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["frac_E"], (int,float)):
                        frac_E = doc["Reaction_Params"][rxn]["frac_E"]
                except:
                    frac_E = 0.2
                try:
                    if isinstance(doc["Reaction_Params"][rxn]["frac_B"], (int,float)):
                        frac_B = doc["Reaction_Params"][rxn]["frac_B"]
                except:
                    frac_B = 0.2

                try:
                    if isinstance(doc["All_Rxn_Unbound"], bool):
                        unbound_A = doc["All_Rxn_Unbound"]
                        unbound_E = doc["All_Rxn_Unbound"]
                        unbound_B = doc["All_Rxn_Unbound"]
                except:
                    pass

                try:
                    if isinstance(doc["All_Rxn_Search_Fraction"], (int, float)):
                        frac_A = doc["All_Rxn_Search_Fraction"]
                        frac_E = doc["All_Rxn_Search_Fraction"]
                        frac_B = doc["All_Rxn_Search_Fraction"]
                        delta_A = None
                        delta_E = None
                        delta_B = None
                except:
                    pass

                self.set_A(rxn,float(doc["Reaction_Params"][rxn]["A"]),unbound_A,frac_A,delta_A)
                self.set_E(rxn,float(doc["Reaction_Params"][rxn]["E"]),unbound_E,frac_E,delta_E)
                self.set_B(rxn,float(doc["Reaction_Params"][rxn]["B"]),unbound_B,frac_B,delta_B)
            for rxn in doc["Reactant_Powers"]:
                for species in doc["Reactant_Powers"][rxn]:
                    self.set_species_rxn_power(species,rxn,float(doc["Reactant_Powers"][rxn][species]))
            for species in doc["Molar_Contributions_to_MBs"]:
                for rxn in doc["Molar_Contributions_to_MBs"][species]:
                    self.set_scale(species,rxn,float(doc["Molar_Contributions_to_MBs"][species][rxn]))

            break


    # Display pprint() of the model instance
    def Display(self):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.pprint()

    # Print results to a file
    def print_to_file(self, name=""):
        if self.solved == False:
            print("Model is not solved!")
            print("\tMust call run_model() first...")
            return
        if name == "" or name == None:
            name = "ReactionModel"
        if self.simulate_only == True:
            name+="-SimulationResults.txt"
        else:
            name+="-OptimizationResults.txt"
        file = open(self.sub_dir+name,"w")
        file.write(" ---- Reactions ---- \n")

        for rxn in self.instance.rxns:
            file.write("\t"+str(rxn)+" =\tk_"+str(rxn))
            for species in self.instance.MBs:
                if abs(value(self.instance.powers[species,rxn])) > 0:
                    if abs(value(self.instance.powers[species,rxn])) > 1:
                        file.write(" * "+str(species)+"^"+str(value(self.instance.powers[species,rxn])))
                    else:
                        file.write(" * "+str(species))
            file.write("\n")

        file.write("\n")
        file.write("\n ---- Rxn Contributions to MBs ----\n")
        for species in self.instance.MBs:
            file.write("\tM_"+str(species)+":")
            i=0
            for rxn in self.instance.rxns:
                if abs(value(self.instance.scale[species,rxn])) > 0:
                    if abs(value(self.instance.scale[species,rxn])) > 1:
                        if i==0:
                            file.write("\t"+str(value(self.instance.scale[species,rxn]))+"*"+str(rxn))
                            i+=1
                        else:
                            file.write(" + "+str(value(self.instance.scale[species,rxn]))+"*"+str(rxn))
                            i+=1
                    else:
                        if i==0:
                            file.write("\t"+str(rxn))
                            i+=1
                        else:
                            file.write(" + "+str(rxn))
                            i+=1
            file.write("\n")

        file.write("\n ---- Mass Balances ---- \n")
        for species in self.instance.MBs:
            file.write("\ttau*("+str(species)+" - "+str(species)+"_in) = -(1 - eps)*M_"+str(species))
            file.write("\n")

        file.write("\n ---- Parameters ---- \n")
        file.write("\ttau =\t"+str(value(self.instance.tau))+"\n")
        file.write("\teps =\t"+str(value(self.instance.eps))+"\n")
        file.write("\tInlet Conc:\n")
        for species in self.instance.MBs:
            file.write("\t\t"+str(species)+"_in =\t"+str(value(self.instance.Cin[species]))+"\n")
        file.write("\tArrhenius Rxn Param:\n")
        for rxn in self.instance.rxns:
            file.write("\t\tA_"+str(rxn)+" =\t"+str(value(self.instance.A[rxn]))+"\n")
            file.write("\t\tB_"+str(rxn)+" =\t"+str(value(self.instance.B[rxn]))+"\n")
            file.write("\t\tE_"+str(rxn)+" =\t"+str(value(self.instance.E[rxn]))+"\n")
        file.write("\n")

        file.write("T[C]")
        for species in self.instance.MBs:
            file.write("\t"+species)
        file.write("\n")
        for temp in self.temperature_set:
            file.write(str(value(self.temperature_set[temp])-273.15))
            for species in self.instance.MBs:
                file.write("\t"+str(value(self.Solution[temp][species])))
            file.write("\n")
        file.close()


    # Create model instance
    def build_instance(self):
        self.model.obj = Objective(rule=self.obj_func)
        self.model.cons = Constraint(self.model.T_set, self.model.MBs, rule=self.mb_constraint)
        self.instance = self.model.create_instance()
        self.built = True

    # Define an initial guess to the concentrations
    def initial_guess(self):
        # Regardless of whether or not we are just simulating, we must fix the kinetics for initial guess
        self.fix_kinetics()
        self.instance.del_component(self.instance.obj)
        self.instance.obj = Objective(expr=0)
        i=0
        for temp_solve in self.instance.T_set:
            if i!=0:
                # Establish initial guesses for next steps
                for species in self.instance.MBs:
                    self.instance.C[species,temp_solve].set_value(value(self.instance.C[species,"T"+str(i-1)]))

            for temp_hold in self.instance.T_set:
                if temp_solve != temp_hold:
                    for species in self.instance.MBs:
                        self.instance.C[species,temp_hold].fix()
                        self.instance.cons[temp_hold,species].deactivate()
                    #End species in hold loop
            #End hold loop
            solver = SolverFactory('ipopt')
            solver.options['print_user_options'] = 'yes'
            solver.options['print_level'] = 5 # Will print out warnings and final analysis only
            solver.options['tol'] = 1e-6
            solver.options['acceptable_tol'] = 1e-6
            solver.options['compl_inf_tol'] = 1e-6
            solver.options['max_iter'] = 20*self.total_var
            solver.options['obj_scaling_factor'] = 1 #Set scaling factor to value similar to tol?
            results = solver.solve(self.instance, tee=True, load_solutions=False)
            try:
                self.instance.solutions.load_from(results)
            except:
                print("\n ---------------- ERROR! -----------------------\n")

            for temp_hold in self.instance.T_set:
                if temp_solve != temp_hold:
                    for species in self.instance.MBs:
                        self.instance.C[species,temp_hold].unfix()
                        self.instance.cons[temp_hold,species].activate()
            i+=1
        #End solve loop
        self.instance.del_component(self.instance.obj)
        self.instance.obj = Objective(rule=self.obj_func)
        if self.simulate_only == False:
            self.unfix_kinetics()

    # Run the model
    def run_model(self):
        # Fix kinetic variables if we are just simulating
        if self.simulate_only == True:
            self.fix_kinetics()

        if self.run_seriel == False:
            # May need to run in seriel mode first to establish a guess
            self.initial_guess()
            solver = SolverFactory('ipopt')
            solver.options['print_user_options'] = 'yes'
            solver.options['print_level'] = 5 # Will print out warnings and final analysis only
            solver.options['tol'] = 1e-6
            solver.options['acceptable_tol'] = 1e-6
            solver.options['compl_inf_tol'] = 1e-6
            solver.options['max_iter'] = 30*self.total_var
            solver.options['obj_scaling_factor'] = 1 #Set scaling factor to value similar to tol?
            if self.simulate_only == False:
                solver.options['diverging_iterates_tol'] = 1e100
            results = solver.solve(self.instance, tee=True, load_solutions=False)
            try:
                self.instance.solutions.load_from(results)
            except:
                for temp in self.instance.T_set:
                    for species in self.instance.MBs:
                        self.Solution[temp][species] = value(self.instance.C[species,temp])
            self.solved = True
            self.obj_value = value(self.instance.obj)

            for temp in self.instance.T_set:
                for species in self.instance.MBs:
                    self.Solution[temp][species] = value(self.instance.C[species,temp])
        else:
            i=0
            for temp in self.temperature_set:
                self.set_temperature("T0",self.temperature_set[temp])
                for species in self.instance.MBs:
                    self.set_initial(species, value(self.instance.C[species,"T0"]))

                solver = SolverFactory('ipopt')
                solver.options['print_user_options'] = 'yes'
                solver.options['print_level'] = 5 # Will print out warnings and final analysis only
                solver.options['tol'] = 1e-6
                solver.options['acceptable_tol'] = 1e-6
                solver.options['compl_inf_tol'] = 1e-6
                solver.options['max_iter'] = 20*self.total_var
                solver.options['obj_scaling_factor'] = 1 #Set scaling factor to value similar to tol?
                results = solver.solve(self.instance, tee=True, load_solutions=False)
                try:
                    self.instance.solutions.load_from(results)
                except:
                    for species in self.instance.MBs:
                        self.Solution["T"+str(i-1)][species] = value(self.instance.C[species,"T0"])
                self.solved = True

                for species in self.instance.MBs:
                    self.Solution[temp][species] = value(self.instance.C[species,"T0"])
                i+=1


    # Add temperatures to simulate
    def add_temperatures(self, Temps):
        self.model.T_set = Set(initialize=Temps)
        self.model.T = Param(self.model.T_set, domain=NonNegativeReals, initialize=298, mutable=True)
        self.temp_set = True

    # Add species list to the model
    def add_species(self, Specs):
        if self.temp_set == False:
            print("Number of temperatures has not been set!")
            print("\tMust call add_temperatures() first...")
            return
        self.model.MBs = Set(initialize=Specs)
        self.total_var += len(Specs)
        self.model.C = Var(self.model.MBs, self.model.T_set, domain=NonNegativeReals, bounds=(1e-20,1e20), initialize=0)
        #self.model.C = Var(self.model.MBs, self.model.T_set, domain=Reals, initialize=0)
        self.model.Cin = Param(self.model.MBs, domain=NonNegativeReals, initialize=0, mutable=True)
        self.mb_set = True

    # Add data to the model
    def add_dataset(self, Specs, RealSpecs):
        if self.mb_set == False:
            print("Number of mass balances has not been set!")
            print("\tMust call add_species() first...")
            return
        for item in Specs:
            if item not in RealSpecs:
                print("Error! Given data value not applicable!")
                return
        self.model.MB_data = Set(initialize=Specs)
        self.model.C_data = Param(self.model.MB_data, self.model.T_set, domain=Reals, initialize=0, mutable=True)
        self.model.wf = Param(self.model.MB_data, self.model.T_set, domain=NonNegativeReals, initialize=1, mutable=True)
        self.data_set = True

    # Add reactions to the system (each MB gets a set of reactions)
    def add_reactions(self, rxns):
        if self.data_set == False:
            print("Number of mass balances and data_set has not been set!")
            print("\tMust call add_dataset() first...")
            return
        self.model.rxns = Set(initialize=rxns)
        self.total_var += 3*len(rxns)
        self.model.scale = Param(self.model.MBs, self.model.rxns, domain=Reals, initialize=0, mutable=True)
        self.model.A = Var(self.model.rxns, domain=NonNegativeReals, initialize=0)
        self.model.B = Var(self.model.rxns, domain=Reals, initialize=0)
        self.model.E = Var(self.model.rxns, domain=Reals, initialize=0)
        self.model.powers = Param(self.model.MBs, self.model.rxns, domain=Any, initialize=0, mutable=True)
        self.rxns_set = True

        for r in rxns:
            self.fix_dict[r] = {}
            self.fix_dict[r]["A"] = False
            self.fix_dict[r]["B"] = False
            self.fix_dict[r]["E"] = False

    # Set initial conditions for a species
    def set_initial(self, species, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        if value <= 0:
            value = 1e-6
        for temp in self.instance.T_set:
            self.instance.C[species,temp].set_value(value)

    # Set inlet condition for a species
    def set_inlet(self, species, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        if value <= 0:
            value = 1e-6
        self.instance.Cin[species].set_value(value)

    # Set data values by temp and species
    def set_data(self, species, temp, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.C_data[species,temp].set_value(value)

    # Set the scale factors for each species, rxn pair
    def set_scale(self, species, rxn, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.scale[species,rxn].set_value(value)

    # Set the pre-exponential A for each reaction
    def set_A(self, rxn, value, unbounded=False, factor=0.2, delta=None):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.A[rxn].set_value(value)
        if factor < 0:
            factor = 0.1
        if factor > 1:
            factor = 0.9
        if unbounded == True:
            return
        if delta == None:
            self.instance.A[rxn].setlb(value*(1.0-factor))
            self.instance.A[rxn].setub(value*(1.+factor))
        else:
            if (value-delta) <= 0:
                self.instance.A[rxn].setlb(1e-16)
            else:
                self.instance.A[rxn].setlb(value-delta)
            self.instance.A[rxn].setub(value+delta)

    # Set the power coefficient B for each reaction
    def set_B(self, rxn, value, unbounded=False, factor=0.2, delta=None):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.B[rxn].set_value(value)
        if factor < 0:
            factor = 0.1
        if factor > 1:
            factor = 0.9
        if unbounded == True:
            return
        if delta == None:
            self.instance.B[rxn].setlb(value*(1.0-factor))
            self.instance.B[rxn].setub(value*(1.+factor))
        else:
            self.instance.B[rxn].setlb(value-delta)
            self.instance.B[rxn].setub(value+delta)

    # Set the activation energy for each reaction
    def set_E(self, rxn, value, unbounded=False, factor=0.2, delta=None):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.E[rxn].set_value(value)
        if factor < 0:
            factor = 0.1
        if factor > 1:
            factor = 0.9
        if unbounded == True:
            return
        if delta == None:
            self.instance.E[rxn].setlb(value*(1.0-factor))
            self.instance.E[rxn].setub(value*(1.+factor))
        else:
            self.instance.E[rxn].setlb(value-delta)
            self.instance.E[rxn].setub(value+delta)

    # Function to fix the kinetic parameters
    def fix_kinetics(self):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        for rxn in self.instance.rxns:
            self.instance.A[rxn].fix()
            self.instance.B[rxn].fix()
            self.instance.E[rxn].fix()
            if value(self.instance.A[rxn]) > 0:
                if self.fix_dict[rxn]["A"] == False:
                    self.total_var = self.total_var - 1
            if abs(value(self.instance.B[rxn])) > 0:
                if self.fix_dict[rxn]["B"] == False:
                    self.total_var = self.total_var - 1
            if abs(value(self.instance.E[rxn])) > 0:
                if self.fix_dict[rxn]["E"] == False:
                    self.total_var = self.total_var - 1

    # Function to unfix kinetics
    def unfix_kinetics(self):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        for rxn in self.instance.rxns:
            if value(self.instance.A[rxn]) > 0:
                if self.fix_dict[rxn]["A"] == False:
                    self.instance.A[rxn].unfix()
                    self.total_var = self.total_var + 1
            if abs(value(self.instance.B[rxn])) > 0:
                if self.fix_dict[rxn]["B"] == False:
                    self.instance.B[rxn].unfix()
                    self.total_var = self.total_var + 1
            if abs(value(self.instance.E[rxn])) > 0:
                if self.fix_dict[rxn]["E"] == False:
                    self.instance.E[rxn].unfix()
                    self.total_var = self.total_var + 1

    # Set the tau value
    def set_tau(self, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.tau.set_value(value)

    # Set the eps value
    def set_eps(self, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        if value > 1 or value < 0:
            print("Porosity value must be between (0,1)...")
            return
        self.instance.eps.set_value(value)

    # Set the temperature value
    def set_temperature(self, temp, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.T[temp].set_value(value)

    # Set reaction species list
    def set_species_rxn_power(self, species, rxn, value):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        self.instance.powers[species,rxn].set_value(value)

    # Function to determine weight factors
    def set_weight_factors(self):
        if self.built == False:
            print("Model is not constructed!")
            print("\tMust call build_instance() first...")
            return
        if self.weight_method == "default":
            return
        elif self.weight_method == "equalizer":
            self.weight_equalizer()
        elif self.weight_method == "relative":
            self.weight_relative()
        else:
            return

        for temp in self.instance.T_set:
            for spec in self.instance.MB_data:
                self.instance.wf[spec, temp].set_value(self.Weight_Factors[spec][temp])

    # Weight Equalization Method
    def weight_equalizer(self):
        for temp in self.instance.T_set:
            sum=0
            for spec in self.instance.MB_data:
                sum+=self.Data[spec][temp]
            for spec in self.instance.MB_data:
                self.Weight_Factors[spec][temp] = (sum/(self.Data[spec][temp]+1e-4))**2

    # Weight Equalization Method
    def weight_relative(self):
        for temp in self.instance.T_set:
            for spec in self.instance.MB_data:
                self.Weight_Factors[spec][temp] = (1.0/(self.Data[spec][temp]+1e-4))**2

    # Return reaction rate
    def comp_rate(self, rxn, model, temp):
        k = rate_const(model.A[rxn], model.B[rxn], model.E[rxn], model.T[temp])
        for name in model.MBs:
            k=k*model.C[name,temp]**model.powers[name,rxn]
        return k

    # Return rhs of the mb
    def rhs_mb(self, species, model, temp):
        return model.tau*(model.C[species, temp] - model.Cin[species])

    # Return the lhs of the mb
    def lhs_mb(self, species, model, temp):
        sum = 0
        for rxn in model.rxns:
            sum += model.scale[species,rxn]*self.comp_rate(rxn, model, temp)
        return sum*-(1 - model.eps)

    # Iterable Constraint Function
    def mb_constraint(self, model, temp, species):
        return self.lhs_mb(species, model, temp) == self.rhs_mb(species, model, temp)

    # Objective function for all mass balances
    def obj_func(self, model):
        sum = 0
        for spec in model.MB_data:
            for temp in model.T_set:
                sum+= model.wf[spec, temp]*(model.C[spec, temp] - model.C_data[spec, temp])**2
        return sum
