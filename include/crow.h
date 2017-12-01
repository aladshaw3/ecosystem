/*!
 *  \file crow.h
 *	\brief Coupled Reaction Object Workspace
 *	\details This file creates objects and subroutines for setting up and solving systems
 *			of reaction driven equations using the DOVE (see dove.h) solver. It combines
 *			a generalized description of chemical reaction mathematics with a comprehensive
 *			input file framework to allow systems of equations to be developed on the fly
 *			and solved with reasonable accuracy and efficiency. 
 *
 *			Mathematical description of ConstReaction (single reaction):
 *			------------------------------------------------------------
 *			du_i/dt = v_i*( k1*Product(j,u_j^v_j) - k2*Product(l,u_l^v_l) )
 *
 *			where i,j,l are indices of variables, k1,k2 are reaction constants, and Product(i,arg)
 *			is the product of all args in i.
 *			NOTE: use sign of stoichiometric coefficients to distinguish reactants from products
 *					(+) = product  and   (-) = reactant
 *
 *			Mathematical description of MultiConstReaction (many reactions):
 *			----------------------------------------------------------------
 *			For a rate expression derived from multiple reaction mechanisms
 *			du_i/dt = SUM( du_i/dt for each reaction in the mechanism )
 *
 *			NOTE: See above for description of single reaction rate functions
 *
 *			\warning This kernel is still under active development. Use with caution!
 *
 *
 *  \author Austin Ladshaw
 *	\date 11/27/2017
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2017, all
 *             rights reserved.
 */

#include "dove.h"

#ifndef crow_h
#define crow_h

/// Enumeration for the list of valid CROW function types
/** This enumeration will define all the function types that have so far been created in CROW. So far,
	the list of valid options is as follows...
 
	\param CONSTREACTION ConstReaction objects for basic chemical reaction mechanisms
	\param MULTICONSTREACTION MultiConstReaction objects for advanced mechanisms with const coeffs
	\param INVALID Default used to denote when a type was not correctly defined*/
typedef enum {CONSTREACTION,
				MULTICONSTREACTION, INVALID} func_type;

/// ConstReaction class is an object for information and functions associated with the Generic Reaction
/** This is a C++ style object designed to store and operate on the generic representation of a 
	reaction mechanism. In this object, the reaction parameters are treated as constants and do
	not change with temperature. This object can be inherited from to add functionality that will
	compute reaction parameters as a function of system temperature. In addition, this object will
	contribute user functions (rate functions and jacobians) to DOVE, which will be
	interfaced through CROW_DATA (see below the class definition). */
class ConstReaction
{
public:
	ConstReaction();							///< Default constructor
	~ConstReaction();							///< Default destructor
	
	void InitializeSolver(Dove &Solver);		///< Function to initialize the ConstReaction object from the Dove object
	void SetIndex(int index);					///< Function to set the index of the function/variable of interest
	void SetForwardRate(double rate);			///< Function to se the forward rate constant of the reaction
	void SetReverseRate(double rate);			///< Function to se the reverse rate constant of the reaction
	void InsertStoichiometry(int i, int v);		///< Insert a Stoichiometric value to the existing map
	
	double getForwardRate();						///< Function to return the forward reaction rate constant
	double getReverseRate();						///< Function to return the reverse reaction rate constant
	int getIndex();									///< Function to return the primary variable index
	std::map<int, int> & getStoichiometryMap();	///< Function to return reference to the Stoichiometry map object
	
protected:
	double forward_rate;						///< Reaction rate constant associated with reactants
	double reverse_rate;						///< Reaction rate constant associated with products
	int main_index;								///< Variable index for the variable of interest
	std::map<int, int> stoic;					///< Map of Stoichiometric coefficients for the reaction (access by var index)
	Dove *SolverInfo;							///< Pointer to the Dove Object
private:
	
};

/// Rate function for the ConstReaction Object
/** This function defines the reaction rate function that will be used in the Dove Object to solve
 the system of equations. Arguments passed to this function are standard and are required in order
 to have this function registered in the Dove object itself. Parameters of this function are as
 follows...
 
 \param i index of the non-linear variable for which this function applies
 \param u matrix of all non-linear variables in the Dove system
 \param t value of time in the current simulation (user must define units)
 \param data pointer to a data structure used to delineate parameters of the function
 \param dove reference to the Dove object itself for access to specific functions*/
double rate_func_ConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove);

/// Jacobi function for the ConstReaction Object
/** This function defines the Jacobian functions that will be used in the Dove Object to precondition
 the system of equations. Arguments passed to this function are standard and are required in order
 to have this function registered in the Dove object itself. Parameters of this function are as
 follows...
 
 \param i index of the non-linear variable for which this function applies
 \param j index of the non-linear variable for which we are taking the derivative with respect to
 \param u matrix of all non-linear variables in the Dove system
 \param t value of time in the current simulation (user must define units)
 \param data pointer to a data structure used to delineate parameters of the function
 \param dove reference to the Dove object itself for access to specific functions*/
double jacobi_func_ConstReaction(int i, int j, const Matrix<double> &u, double t, const void *data, const Dove &dove);

/// MultiConstReaction class is an object associated with rate functions deriving from multiple reactions
/** This is a C++ style object designed to store and operate on the generic representation of multiple
	reaction mechanisms. In this object, the reaction parameters are treated as constants and do
	not change with temperature. This object can be inherited from to add functionality that will
	compute reaction parameters as a function of system temperature. In addition, this object will
	contribute user functions (rate functions and jacobians) to DOVE, which will be
	interfaced through CROW_DATA (see below the class definition).*/
class MultiConstReaction
{
public:
	MultiConstReaction();							///< Default constructor
	~MultiConstReaction();							///< Default destructor
	
	void SetNumberReactions(unsigned int i);		///< Set the number of reactions for the rate function
	void InitializeSolver(Dove &Solver);			///< Function to initialize each ConstReaction object from the Dove object
	void SetIndex(int index);						///< Function to set the index of the primary variable species (same for all reactions)
	
	void SetForwardRate(int react, double rate);	///< Function to set the forward reaction rate for the indicated reaction
	void SetReverseRate(int react, double rate);	///< Function to set the reverse reaction rate for the indicated reaction
	void InsertStoichiometry(int react, int i, int v);	///< Function to insert stoichiometry for the given reaction
	
	int getNumReactions();							///< Return the number of reactions in the object
	ConstReaction &getReaction(int i);				///< Return reference to the ith ConstReaction object
	
protected:
	int num_reactions;									///< Number of reaction objects
	std::vector<ConstReaction> reactions;				///< List of reaction objects associated with MultiConstReaction
	
private:
	
};

/// Rate function for the MultiConstReaction Object
/** This function defines the reaction rate function that will be used in the Dove Object to solve
 the system of equations. Arguments passed to this function are standard and are required in order
 to have this function registered in the Dove object itself. Parameters of this function are as
 follows...
 
 \param i index of the non-linear variable for which this function applies
 \param u matrix of all non-linear variables in the Dove system
 \param t value of time in the current simulation (user must define units)
 \param data pointer to a data structure used to delineate parameters of the function
 \param dove reference to the Dove object itself for access to specific functions*/
double rate_func_MultiConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove);

/// Jacobi function for the MultiConstReaction Object
/** This function defines the Jacobian functions that will be used in the Dove Object to precondition
 the system of equations. Arguments passed to this function are standard and are required in order
 to have this function registered in the Dove object itself. Parameters of this function are as
 follows...
 
 \param i index of the non-linear variable for which this function applies
 \param j index of the non-linear variable for which we are taking the derivative with respect to
 \param u matrix of all non-linear variables in the Dove system
 \param t value of time in the current simulation (user must define units)
 \param data pointer to a data structure used to delineate parameters of the function
 \param dove reference to the Dove object itself for access to specific functions*/
double jacobi_func_MultiConstReaction(int i, int j, const Matrix<double> &u, double t, const void *data, const Dove &dove);

/// Primary data structure for the CROW routines
/** This is a c-style data structure used to house all CROW data necessary to perform simulations on
	systems of non-linear reaction equations. It is the primary structure that will interface with DOVE
	and be passed to the DOVE object as the user_data structure. Nested within this structure will be
	all the parameter information necessary to delineate and evaluate the functions registered in DOVE. */
typedef struct CROW_DATA
{
	std::unordered_map<int, ConstReaction> const_reacts;	///< Map of constant reaction objects used in CROW (access by variable index)
	std::unordered_map<int, MultiConstReaction> multi_const_reacts;	///< Map of constant reaction objects used in CROW (access by variable index)
	Dove SolverInfo;										///< Dove object that holds all information associated with the solver
	FILE *OutputFile;										///< Pointer to the output file for CROW
	bool FileOutput = true;									///< Boolean to determine whether or not to print results to a file
	yaml_cpp_class yaml_object;								///< yaml object to read and access digitized yaml documents (see yaml_wrapper.h)
} CROW_DATA;

///Function to print header information about CROW to output file
void print2file_crow_header(CROW_DATA *dat);

///Function to validate solver choice
//** Returns true for Linear and false for Nonlinear */
bool solver_choice(std::string &choice);

///Function to validate linesearch choice
linesearch_type linesearch_choice(std::string &choice);

///Function to validate linear solver choice
krylov_method linearsolver_choice(std::string &choice);

///Function to determine whether or not to precondition
bool use_preconditioning(std::string &choice);

///Function to validate preconditioning choice
precond_type preconditioner_choice(std::string &choice);

///Function to validate timestepper choice
timestep_type timestepper_choice(std::string &choice);

///Function to validate integration method choice
integrate_subtype integration_choice(std::string &choice);

///Function to validate Function type
func_type function_choice(std::string &choice);

///Function to read a yaml input file to setup a CROW simulation
int read_crow_input(CROW_DATA *dat);

///Function to intialize system information
int read_crow_system(CROW_DATA *dat);

///Function to read the header files for each variable
int read_crow_functions(CROW_DATA *dat);

///Function to initialize ConstReaction information
int read_crow_ConstReaction(int index, std::unordered_map<int, ConstReaction> &map, Document &info, Dove &solver);

///Function to initialize MultiConstReaction information
int read_crow_MultiConstReaction(int index, std::unordered_map<int, MultiConstReaction> &map, Document &info, Dove &solver);

///Run CROW scenario
int CROW_SCENARIO(const char *yaml_input);

///Run the CROW test
int CROW_TESTS();

#endif /* crow_h */
