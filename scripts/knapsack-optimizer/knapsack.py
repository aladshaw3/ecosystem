''' knapsack script:
    ----------------
    Object-Oriented approach to solving the quantitized optimization problem
    with generic constraints and generic objective functions. User specifies
    constraints and objective functions and passes those instructions to the
    ZeroOneKnapsack object. Then, by running the optimizer, the object will
    find the sub-set of items, given an available list of items, that will
    maximize the given objective function, subject to the user's constraints. 
    
    Uses: This script will be used as part of the uncertainty analyses needed
            to assess data needs and/or model sensitivity to a discrete 
            number of control parameters. 
    
    Author:     Austin Ladshaw
    Date:       05/06/2019
    Copyright:  This software was designed and built at the Georgia Institute
                of Technology by Austin Ladshaw for research in the area of
                radioactive particle decay and transport. Copyright (c) 2019,
                all rights reserved.'''

## Create Class Objects ##

class Item(object):
    #Item object: each item as a name, value, and unique id
    tag = 0
    def __init__(self):
        self.name = " "             #Name of the object
        self.value = 0              #Value associated with that object
        self.id = Item.tag          #Item id tag (unique for each instance of the object)
        Item.tag += 1
    def __str__(self):
        return self.name + ' (' + str(self.id) + '): <' + str(self.value) + '>'
    def set_name(self, name):
        self.name = name
    def set_value(self, value):
        self.value = value
    def get_name(self):
        return self.name
    def get_value(self):
        return self.value
    def get_id(self):
        return self.id

class ZeroOneKnapsack(object):
    #ZeroOneKnapsack object: solves the 0/1 knapsack for a
    #                   generic set of objectives and constraints
    def __init__(self):
        self.obj_func = []          #Storage for the objective function (passed by user)
        self.obj_args = []          #Arguments tuple for the objective function (passed by user)
        self.constraint_func = []   #Storage for the constraint function (passed by user)
        self.constraint_args = []   #Arguments tuple for the constraint function (passed by user)
        self.whole_tree = True      #Boolean to force exhaustive search of tree
    
    # Function to set the exhaustive search boolean
    def exhaustive_search(self, opt):
        self.whole_tree = opt

    # Function to register the objective function and it's arguments it needs
    def register_objective_func(self, obj_func, *obj_args):
        if (len(self.obj_func) == 0):
            self.obj_func.append(obj_func)
            self.obj_args.append(obj_args)
        else:
            self.obj_func[0] = obj_func
            self.obj_args[0] = obj_args

    # Function to register the constraints function and it's arguments it needs
    def register_constraints(self, con_fun, *con_args):
        if (len(self.constraint_func) == 0):
            self.constraint_func.append(con_fun)
            self.constraint_args.append(con_args)
        else:
            self.constraint_func[0] = con_fun
            self.constraint_args[0] = con_args

    ''' NOTE: Both the constraint function and objective function supplied by the user must accept 
                two arguments: (i) list of items and (ii) tuple of other any other arguments.
        
              Return types for each function are to be expected as follows:
                    Constraint Function --> returns boolean (True = constraints met)
                    Objective Function --> returns scalar 
                                                Postive Scalar == Maximization
                                                Negative Scalar == Minimization 
        '''
    
    def eval_constraints(self, list):
        return self.constraint_func[0](list, len(list), *self.constraint_args)

    def eval_objective_func(self, list):
        return self.obj_func[0](list, len(list), *self.obj_args)

    # Function to solve the optimization problem given the original list of items to consider
    ##      NOTE: User should not pass any arguments for 'taken' as this is used only internally
    def Optimize(self, list, taken=None):
        if taken is None:
            taken = []
        #Check constraints, compute value, view list
        status = self.eval_constraints(taken)
        #value = self.eval_objective_func(taken)
        value = 0.0
        result = (value, taken, status)
        
        # Check for empty list
        if list == []:
            value = self.eval_objective_func(taken)
            result = (value, taken, status)
        else:
            # If constraints are violated or search is exhaustive, force a move to left and right
            if status == False or self.whole_tree == True:
                violation = True
                i = 0
                #Continue to loop while constraints are still violated
                while violation:
                    try:
                        nextItem = list[i]
                    except:
                        break
                    i += 1
                    #Create a copy of current taken list before adding new item
                    not_taken = taken.copy()
                    taken.append(nextItem)

                    #Move down left branch
                    (left_val, left_take, left_status) = self.Optimize(list[i:],taken)

                    #Move down right branch
                    (right_val, right_take, right_status) = self.Optimize(list[i:],not_taken)
            
                    #Check the results
                    if (left_status == True) and (right_status == True):
                        violation = False
                        # Take largest value
                        if right_val > left_val:
                            result = (right_val, right_take, right_status)
                        else:
                            result = (left_val, left_take, left_status)
                    elif (left_status == True) and (right_status == False):
                        violation = False
                        result = (left_val, left_take, left_status)
                    elif (left_status == False) and (right_status == True):
                        violation = False
                        result = (right_val, right_take, right_status)
                    else:
                        # If constraints are still violated, then force another move
                        violation = True
                        if right_val > left_val:
                            result = (right_val, right_take, right_status)
                        else:
                            result = (left_val, left_take, left_status)
                #END While Loop

            # If constraints are not violated and search is shallow, check left and check right
            else:
                #Create copies to pass to arguments
                not_taken = taken.copy()
                taken.append(list[0])
            
                #Peak left and right to look for a move violation
                left_status = self.eval_constraints(taken)
                right_status = self.eval_constraints(not_taken)
                
                # Right move valid, but left isn't, then move right
                if left_status == False and right_status == True:
                    (right_val, right_take, right_status) = self.Optimize(list[1:],not_taken)
                    result = (right_val, right_take, right_status)
                
                # Left move valid, but right isn't, then move left
                elif left_status == True and right_status == False:
                    (left_val, left_take, left_status) = self.Optimize(list[1:],taken)
                    result = (left_val, left_take, left_status)

                # Both moves valid, then explore both branches
                elif left_status == True and right_status == True:
                    (left_val, left_take, left_status) = self.Optimize(list[1:],taken)
                    (right_val, right_take, right_status) = self.Optimize(list[1:],not_taken)
                
                    # Take largest value
                    if right_val > left_val:
                        result = (right_val, right_take, right_status)
                    else:
                        result = (left_val, left_take, left_status)

                # No moves valid, return current result
                else:
                    value = self.eval_objective_func(taken)
                    result = (value, taken, status)

        return result
    #End Optimization Routine
