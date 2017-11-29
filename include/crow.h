/*!
 *  \file crow.h
 *	\brief Coupled Reaction Object Workspace
 *	\details This file creates objects and subroutines for setting up and solving systems
 *			of reaction driven equations using the DOVE (see dove.h) solver. It combines
 *			a generalized description of chemical reaction mathematics with a comprehensive
 *			input file framework to allow systems of equations to be developed on the fly
 *			and solved with reasonable accuracy and efficiency. 
 *
 *			Mathematical description of generic reaction:
 *			---------------------------------------------
 *			a_i*du_i/dt = k1*Product(j,u_j^v_j) - k2*Product(l,u_l^v_l)
 *
 *			where i,j,l are indices of variables, k1,k2 are reaction constants, a is a time
 *			coefficient, and Product(i,arg) is the product of all args in i.
 *			NOTE: a_i = 1/v_i for a single reaction mechanism equation
 *			NOTE: use sign of stoichiometric coefficients to distinguish reactants from products
 *					(+) = product  and   (-) = reactant
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

/// ConstReaction class is an object for information and functions associated with the Generic Reaction
/** This is a C++ style object designed to store and operate on the generic representation of a 
	reaction mechanism. In this object, the reaction parameters are treated as constants and do
	not change with temperature. This object can be inherited from to add functionality that will
	compute reaction parameters as a function of system temperature. In addition, this object will
	contribute user functions (time coeffs, rate functions, and jacobians) to DOVE, which will be
	interfaced through CROW_DATA. */
class ConstReaction
{
public:
	ConstReaction();							///< Default constructor
	~ConstReaction();							///< Default destructor
	
	void InitializeSolver(Dove &Solver);		///< Function to initialize the ConstReaction object from the Dove object
	void SetIndex(int index);					///< Function to set the index of the function/variable of interest
	void SetForwardRate(double rate);			///< Function to se the forward rate constant of the reaction
	void SetReverseRate(double rate);			///< Function to se the reverse rate constant of the reaction
	void InsertStoichiometry(int i, int v);	///< Insert a Stoichiometric value to the existing map
	
	void ComputeTimeCoeff();					///< Function to compute and store the time coefficient
	
	double getForwardRate();						///< Function to return the forward reaction rate constant
	double getReverseRate();						///< Function to return the reverse reaction rate constant
	double getTimeCoeff();							///< Function to return the time coefficient constant
	int getIndex();									///< Function to return the primary variable index
	std::map<int, int> & getStoichiometryMap();	///< Function to return reference to the Stoichiometry map object
	
protected:
	double forward_rate;						///< Reaction rate constant associated with reactants
	double reverse_rate;						///< Reaction rate constant associated with products
	double time_coeff;							///< Coefficient for the time derivative of the variable
	int main_index;								///< Variable index for the variable of interest
	std::map<int, int> stoic;				///< Map of Stoichiometric coefficients for the reaction (access by var index)
	Dove *SolverInfo;							///< Pointer to the Dove Object
private:
	
};

/// Coefficient function for the ConstReaction Object
/** This function defines the time coefficient function that will be used in the Dove Object to solve
 the system of equations. Arguments passed to this function are standard and are required in order
 to have this function registered in the Dove object itself. Parameters of this function are as
 follows...
 
 \param i index of the non-linear variable for which this function applies
 \param u matrix of all non-linear variables in the Dove system
 \param t value of time in the current simulation (user must define units)
 \param data pointer to a data structure used to delineate parameters of the function
 \param dove reference to the Dove object itself for access to specific functions*/
double coeff_func_ConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove);

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

/// Primary data structure for the CROW routines
/** This is a c-style data structure used to house all CROW data necessary to perform simulations on
	systems of non-linear reaction equations. It is the primary structure that will interface with DOVE
	and be passed to the DOVE object as the user_data structure. Nested within this structure will be
	all the parameter information necessary to delineate and evaluate the functions registered in DOVE. */
typedef struct CROW_DATA
{
	std::vector<ConstReaction> const_reacts;	///< List of constant reaction objects used in CROW
	Dove SolverInfo;							///< Dove object that holds all information associated with the solver
	FILE *OutputFile;							///< Pointer to the output file for CROW
	bool FileOutput = true;						///< Boolean to determine whether or not to print results to a file
} CROW_DATA;

///Function to print header information about CROW to output file
void print2file_crow_header(CROW_DATA *dat);

///Run CROW scenario
int CROW_SCENARIO(const char *yaml_input);

///Run the CROW test
int CROW_TESTS();

#endif /* crow_h */
