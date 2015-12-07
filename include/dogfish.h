/*!
 *  \file dogfish.h dogfish.cpp
 *	\brief Diffusion Object Governing Fiber Interior Sorption History
	\details This set of objects and functions is used to numerically solve linear or non-linear
			diffusion physics of aqueous ions into cylindrical adsorbent fibers. Boundary conditions
			for this problem could be a film mass transfer, reaction, or dirichlet condition depending
			on the type of problem being solve.
 
 *	\warning Functions and methods in this file are still under construction.
 *  \author Austin Ladshaw
 *	\date 04/09/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#ifndef DOGFISH_HPP_
#define DOGFISH_HPP_

#include "finch.h"
#include "mola.h"

/// Data structure for species-specific parameters
/** C-style object to hold information on all adsorbing species. Parameters are 
	given descriptive names to indicate what each is for. */
typedef struct
{
	//Default parameters to go with default functions
	double intraparticle_diffusion;						///< Units: um^2/hr
	double film_transfer_coeff;							///< Units: um/hr
	double surface_concentration;						///< Units: mol/kg
	double initial_sorption;							///< Units: mol/kg
	
	//Additional info
	double sorbed_molefraction;							///< Molefraction of sorbed species
	Molecule species;									///< Adsorbed species Molecule Object
	
}DOGFISH_PARAM;

/// Primary data structure for running the DOGFISH application
/** C-style object to hold information for the adsorption simulations. Contains function
	pointers and other data structures. This information is passed around to other functions
	used to simulate the fiber diffusion physics. */
typedef struct
{
	unsigned long int total_steps = 0;		///< Total number of solver steps taken
	double time_old = 0.0;					///< Old value of time (hrs)
	double time = 0.0;						///< Current value of time (hrs)
	bool Print2File = true;					///< True = results to .txt; False = no printing
	bool Print2Console = true;				///< True = results to console; False = no printing
	bool DirichletBC = false;				///< False = uses film mass transfer for BC, True = Dirichlet BC
	bool NonLinear = false;					///< False = Solve directly, True = Solve iteratively
	double t_counter = 0.0;					///< Counter for the time output
	double t_print;							///< Print output at every t_print time (hrs)
	
	int NumComp;							///< Number of species to track
	double end_time;						///< Units: hours
	double total_sorption_old;				///< Per mass or volume of single fiber
	double total_sorption;					///< Per mass or volume of single fiber
	double fiber_length;					///< Units: um
	double fiber_diameter;					///< Units: um
	double fiber_specific_area;				///< Units: m^2/kg
	
	FILE *OutputFile;						///< Output file pointer to the output file for postprocesses and results
	
	//Function pointers for the parameters to be evaluated in FINCH for DOGFISH

	double (*eval_R) (int i, int l, const void *data);	///< Function pointer to evaluate retardation coefficient
	double (*eval_DI) (int i, int l, const void *data);	///< Function pointer to evaluate intraparticle diffusivity
	double (*eval_kf) (int i, const void *data);		///< Function pointer to evaluate film mass transfer coefficient
	double (*eval_qs) (int i, const void *data);		///< Function pointer to evaluate fiber surface concentration
	const void *user_data;					///< Data structure for users info to calculate all parameters
	std::vector<FINCH_DATA> finch_dat;		///< Data structure for FINCH_DATA objects
	std::vector<DOGFISH_PARAM> param_dat;	///< Data structure for DOGFISH_PARAM objects
	
}DOGFISH_DATA;

/// Function to print a species based header for the output file
void print2file_species_header(FILE *Output, DOGFISH_DATA *dog_dat, int i);

/// Function to print a time and space header for the output file
void print2file_DOGFISH_header(DOGFISH_DATA *dog_dat);

/// Function to print out the old time results for the output file
void print2file_DOGFISH_result_old(DOGFISH_DATA *dog_dat);

/// Function to print out the new time results for the output file
void print2file_DOGFISH_result_new(DOGFISH_DATA *dog_dat);

//----- Default functions for DOGFISH parameters -----

/// Default function for the retardation coefficient
/** The default retardation coefficient for this problem is 1.0 for all time and space. Therefore,
	this function will only ever return a 1.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param data pointer to the DOGFISH_DATA structure*/
double default_Retardation(int i, int l, const void *data);

/// Default function for the intraparticle diffusion coefficient
/** The default intraparticle diffusivity is to assume that each species i has a constant diffusivity.
	Therefore, this function returns the value of the parameter intraparticle_diffusion from the 
	DOGFISH_PARAM structure for each adsorbing species i. Each species may have a different diffusivity.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param data pointer to the DOGFISH_DATA structure*/
double default_IntraDiffusion(int i, int l, const void *data);

/// Default function for the film mass transfer coefficient
/** The default film mass transfer coefficient will be to assume that this value is a constant for
	each species i. Therefore, this function returns the parameter value of film_transfer_coeff from
	the DOGFISH_PARAM structure for each adsorbing species i.
 
	\param i index for the ith adsorbing species
	\param data pointer to the DOGFISH_DATA structure */
double default_FilmMTCoeff(int i, const void *data);

/// Default function for the fiber surface concentration
/** The default fiber surface concentration will be to assume that this value is a constant for
	each species i. Therefore, this function returns the parameter value of surface_concentration from
	the DOGFISH_PARAM structure for each adsorbing species i.
 
	\param i index for the ith adsorbing species
	\param data pointer to the DOGFISH_DATA structure*/
double default_SurfaceConcentration(int i, const void *data);
//----- END Default function definitions -------------

/// \cond
double LangmuirSurfaceConcentration(int i, const void *data);
/// \endcond

/// Function will set up the memory and pointers for use in the DOGFISH simulations
/** The pointers to the output file, parameter functions, and data structures are passed into
	this function to setup the problem in memory. This function must always be called prior to
	calling any other DOGFISH routine and after the DOGFISH_DATA structure has been initialized.
 
	\param file pointer to the output file to print out results
	\param eval_R function pointer for the retardation coefficient function
	\param eval_DI function pointer for the intraparticle diffusion function
	\param eval_kf function pointer for the film mass transfer function
	\param eval_qs function pointer for the surface concentration function
	\param user_data pointer for the user's own data structure (only if using custom functions)
	\param dog_dat pointer for the DOGFISH_DATA structure*/
int setup_DOGFISH_DATA(FILE *file, double (*eval_R) (int i, int l, const void *user_data),
					   double (*eval_DI) (int i, int l, const void *user_data),
					   double (*eval_kf) (int i, const void *user_data),
					   double (*eval_qs) (int i, const void *user_data),
					   const void *user_data, DOGFISH_DATA *dog_dat);

/// Function to serially call all other functions need to solve the system at one time step
/** This function will call the DOGFISH_preprocesses function, followed by the FINCH solver
	functions for each species i, then call the DOGFISH_postprocesses function. After completion,
	this would have solved the diffusion physics for a single time step. */
int DOGFISH_Executioner(DOGFISH_DATA *dog_dat);

/// Function called to evaluate the initial conditions for the time dependent problem
/** This function will use information in DOGFISH_DATA to setup the initial conditions,
	initial parameter values, and initial sorption averages for each species. This function
	always assumes a constant initial condition for the sorption of each species. */
int set_DOGFISH_ICs(DOGFISH_DATA *dog_dat);

/// Function sets the time step size for the next step forward in the simulation
/** This function will set the next time step size based on the spatial discretization 
	of the fiber. Maximum time step size is locked at 0.5 hours. */
int set_DOGFISH_timestep(DOGFISH_DATA *dog_dat);

/// Function to perform preprocess actions to be used before calling any solver
/** This function will call all of the parameter functions in order to establish boundary
	condition parameter values prior to calling the FINCH solvers. */
int DOGFISH_preprocesses(DOGFISH_DATA *dog_dat);

/// Function to calculate the values of all parameters for all species at all nodes
/** This function is passed to the FINCH_DATA data structure and set as the setparams function
	pointer. FINCH calls this function during it's solver routine to setup the non-linear form
	of the problem and solve the non-linear system. 
 
	\param user_data this is actually the DOGFISH_DATA structure, but is passed anonymously to FINCH */
int set_DOGFISH_params(const void *user_data);

/// Function to perform post-solve actions such as printing out results
/** This function increments the total_steps counter in DOGFISH_DATA to keep a running total
	of all solver steps taken. Additionally, it prints out the results of the current time
	simulation to the output file. */
int DOGFISH_postprocesses(DOGFISH_DATA *dog_dat);

/// Function to reset the matrices and vectors and prepare for next time step
/** This function will reset the matrix and vector information of DOGFISH_DATA and FINCH_DATA
	to prepare for the next simulation step in time. */
int DOGFISH_reset(DOGFISH_DATA *dog_dat);

/// Function performs all necessary steps to step the diffusion simulation through time
/** This function calls the initial conditions, set time step, executioner, and reset functions 
	to step the simulation through time. It will only exit when the simulation time is reached
	or if an error occurs. */
int DOGFISH(DOGFISH_DATA *dog_dat);


/// Running DOGFISH tests
/** This function is called from the UI to run a test simulation of DOGFISH. Ouput is stored
	in a DOGFISH_TestOutput.txt file in a sub-directory "output" from the directory in which
	the executable was called. */
int DOGFISH_TESTS();

#endif
