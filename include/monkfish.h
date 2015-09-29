/*!
 *  \file monkfish.h monkfish.cpp
 *	\brief Multi-fiber wOven Nest Kernel For Interparticle Sorption History
 *	\details This file contains structures and functions associated with modeling
 *			the sorption characteristics of woven fiber bundles used to recover 
 *			uranium from seawater. It is coupled with the DOGFISH kernel that
 *			determines the sorption of individual fibers. This kernel will resolve
 *			the interparticle diffusion between bundles of individual fibers in a
 *			woven ball-like domain.
 *
 *	\warning Functions and methods in this file are still under construction.
 *  \author Austin Ladshaw
 *	\date 04/14/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#ifndef MONKFISH_HPP_
#define MONKFISH_HPP_

#include "dogfish.h"

/// Data structure for species specific information and parameters
/** C-style object to hold information associated with the different species
	present in the interparticle diffusion problem. Each species may have different
	diffusivities, mass transfer coefficients, etc. Average adsorption for each species
	will be held in matrix objects. */
typedef struct
{
	double interparticle_diffusion;			///< Units: cm^2/hr
	double exterior_concentration;			///< Units: mol/L
	double exterior_transfer_coeff;			///< Units: cm/hr
	double sorbed_molefraction;				///< Units: -
	double initial_sorption;				///< Units: mg/g
	double sorption_bc;						///< Units: mg/g
	
	double intraparticle_diffusion;			///< Units: um^2/hr
	double film_transfer_coeff;				///< Units: um/hr
	
	Matrix<double> avg_sorption;			///< Units: mg/g
	Matrix<double> avg_sorption_old;		///< Units: mg/g
	
	Molecule species;						///< Species in the liquid phase
	
}MONKFISH_PARAM;

/// Primary data structure for running MONKFISH
/** C-style object holding simulation information for MONKFISH as well as common system
	parameters like fiber density, fiber diameter, fiber length, etc. This object also 
	contains function pointers to different parameter evaluation functions that can be
	changed to suit a particular problem. Default functions will be given, so not every
	user needs to override these functions. This structure also contains vectors of other
	objects including FINCH and DOGFISH objects to resolve the diffusion physics at both
	the macro- and micro-scale. */
typedef struct
{
	unsigned long int total_steps = 0;		///< Total number of steps taken by the algorithm (iterations and time steps)
	double time_old = 0.0;					///< Old value of time in the simulation (hrs)
	double time = 0.0;						///< Current value of time in the simulation (hrs)
	bool Print2File = true;					///< True = results to .txt; False = no printing
	bool Print2Console = true;				///< True = results to console; False = no printing
	bool DirichletBC = true;				///< False = uses film mass transfer for BC, True = Dirichlet BC
	bool NonLinear = false;					///< False = Solve directly, True = Solve iteratively
	bool haveMinMax = false;				///< True = know min and max fiber density, False = only know avg density (Used in ICs)
	bool MultiScale = true;					///< True = solve single fiber model at nodes, False = solve equilibrium at nodes
	int level = 2;							///< Level of coupling between multiple scales (default = 2)
	double t_counter = 0.0;					///< Counter for the time output
	double t_print;							///< Print output at every t_print time (hrs)
	
	int NumComp;							///< Number of species to track
	double end_time;						///< Units: hours
	double total_sorption_old;				///< Old total adsorption per mass of woven nest (mg/g)
	double total_sorption;					///< Current total adsorption per mass woven nest (mg/g)
	double single_fiber_density;			///< Units: g/L
	double avg_fiber_density;				///< Units: g/L (Used in ICs)
	double max_fiber_density;				///< Units: g/L (Used in ICs)
	double min_fiber_density;				///< Units: g/L (Used in ICs)
	double max_porosity;					///< Units: -
	double min_porosity;					///< Units: -
	double domain_diameter;					///< Nominal diameter of the woven fiber ball - Units: cm
	
	FILE *Output;							///< Output file pointer for printing to text file
	
	//Function pointers for the parameters to be evaluated in FINCH for MONKFISH
	
	/// Function pointer to evaluate the porosity of the woven bundle of fibers
	double (*eval_eps) (int i, int l, const void *user_data);
	
	/// Function pointer to evaluate the fiber density in the domain
	double (*eval_rho) (int i, int l, const void *user_data);
	
	/// Function pointer to evaluate the interparticle diffusivity
	double (*eval_Dex) (int i, int l, const void *user_data);
	
	/// Function pointer to evaluate the adsorption strength for the macro-scale
	double (*eval_ads) (int i, int l, const void *user_data);
	
	/// Function pointer to evaluate the retardation coefficient for the macro-scale
	double (*eval_Ret) (int i, int l, const void *user_data);
	
	/// Function pointer to evaluate the exterior concentration for the domain
	double (*eval_Cex) (int i, const void *user_data);
	
	/// Function pointer to evalutate the film mass transfer coefficient for the macro-scale
	double (*eval_kf) (int i, const void *user_data);
	
	const void *user_data;					///< User supplied data function to evaluate the function pointers (Default = MONKFISH_DATA)
	std::vector<FINCH_DATA> finch_dat;		///< FINCH data structures to solve each species interparticle diffusion equation
	std::vector<MONKFISH_PARAM> param_dat;	///< MONKFISH parameter data structure for each species adsorbing
	std::vector<DOGFISH_DATA> dog_dat;		///< DOGFISH data structures for each node in the macro-scale problem
	
}MONKFISH_DATA;

//--------Default Parameter Functions---------------

/// Default porosity function for MONKFISH
/** This function assumes a linear relationship between the maximum porosity at the center of
	the woven fibers and the minimum porosity at the edge of the woven fiber bundle. 
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_porosity(int i, int l, const void *user_data);

/// Default density function for MONKFISH
/** This function calls the porosity function and uses the single fiber density to provide an
	estimate of the bulk fiber density locally in the woven fiber bundle.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_density(int i, int l, const void *user_data);

/// Default interparticle diffusion function
/** This function assumes that the interparticle diffusivity is a contant and returns that 
	diffusivity multiplied by the domain porosity to form the effective diffusion coefficient
	in the domain.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_interparticle_diffusion(int i, int l, const void *user_data);

/// Default adsorption strength function
/** This function will either use the default equilibrium function or the DOGFISH simulation
	result to produce the approximate adsorption strength using perturbation theory.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_monk_adsorption(int i, int l, const void *user_data);

/// Default equilibirium adsorption function in mg/g
/** This function uses the exterior species' concentration (mol/L), the species' molecular 
	weight (g/mol), and the bulk fiber density (g/L) to calculate the adsorption equilibrium
	in mg/g. It assumes that the exterior concentration represents the moles of species per 
	liter of solution that is being sorbed.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_monk_equilibrium(int i, int l, const void *user_data);

/// Default retardation coefficient function
/** This function calls the porosity, density, and adsorption functions to evaluate the 
	retardation coefficient of the diffusing material.
 
	\param i index for the ith adsorbing species
	\param l index for the lth node in the domain
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_monkfish_retardation(int i, int l, const void *user_data);

/// Default exterior concentratio function
/** This function assumes that the exterior concentration for sorption is just equal to 
	the value of exterior_concentration given in MONKFISH_PARAM.
 
	\param i index for the ith adsorbing species
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_exterior_concentration(int i, const void *user_data);

/// Default film mass transfer function
/** This function assumes that the film mass transfer coefficient is just equal to the
	value of the film_transfer_coeff in MONKFISH_PARAM.
 
	\param i index for the ith adsorbing species
	\param user_data pointer to the MONKFISH_DATA structure*/
double default_film_transfer(int i, const void *user_data);
//--------End Default Parameter Functions-----------

/// Setup function to allocate memory and setup function pointers for the MONKFISH simulation
/** This function will allocate memory and setup the MONKFISH problem. To specify use of the
	default functions in MONKFISH, pass NULL args for all function pointers and the user_data
	data structure. Otherwise, pass in your own custom arguments. The MONKFISH_DATA pointer
	must always be passed to this function. 
 
	\param file pointer to the output file to print out results
	\param eval_porosity function pointer for the bulk domain porosity function
	\param eval_density function pointer for the bulk domain density function
	\param eval_ext_diff function pointer for the interparticle diffusion function
	\param eval_adsorb function pointer for the adsorption strength function
	\param eval_retard function pointer for the retardation coefficient function
	\param eval_ext_conc function pointer for the external concentration function
	\param eval_ext_film function pointer for the external film mass transfer function
	\param dog_diffusion function pointer for the DOGFISH diffusion function (see dogfish.h)
	\param dog_ext_film function pointer for the DOGFISH film mass transfer (see dogfish.h)
	\param dog_surf_conc function pointer for the DOGFISH surface concentration (see dogfish.h)
	\param user_data pointer for the user's own data structure (only if using custom functions)
	\param monk_dat pointer for the MONKFISH_DATA structure*/
int setup_MONKFISH_DATA(FILE *file,
						double (*eval_porosity) (int i, int l, const void *user_data),
						double (*eval_density) (int i, int l, const void *user_data),
						double (*eval_ext_diff) (int i, int l, const void *user_data),
						double (*eval_adsorb) (int i, int l, const void *user_data),
						double (*eval_retard) (int i, int l, const void *user_data),
						double (*eval_ext_conc) (int i, const void *user_data),
						double (*eval_ext_film) (int i, const void *user_data),
						double (*dog_diffusion) (int i, int l, const void *user_data),
						double (*dog_ext_film) (int i, const void *user_data),
						double (*dog_surf_conc) (int i, const void *user_data),
						const void *user_data, MONKFISH_DATA *monk_dat);


/// Function to run tests on the MONKFISH algorithms
/** This function currently does nothing and is not callable from the UI. */
int MONKFISH_TESTS();

#endif
