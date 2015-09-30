/*!
 *  \file scopsowl_opt.h scopsowl_opt.cpp
 *	\brief Optimization Routine for Surface Diffusivities in SCOPSOWL
 *	\details This file contains structures and functions associated with performing non-linear
 *			least-squares optimization of the SCOPSOWL simulation results against actual kinetic
 *			adsorption data. The optimization routine here allows you to run data comparisons and
 *			optimizations in three forms: (i) Rough optimizations - cheaper operations, but less
 *			accurate, (ii) Exact optmizations - much more expensive, but greater accuracy, and (iii)
 *			data/model comparisons - no optimization, just using system parameters to compare
 *			simulation results agains a set of data. 
 *
 *			Depending on the level of optimization desired, this routine could take several minutes
 *			or several hours. The optimization/comparisons are printed out in two files: (i) a parameter
 *			file, which contains the simulation partial pressures and temperatures and the optimized 
 *			diffusivities with the euclidean norm of the fitting and (ii) a comparison file that 
 *			shows the model value and data value at each time step for each kinetic curve. 
 *
 *			The optimized diffusion parameters are given for each individual kinetic data curve. Each
 *			data curve will have a different pairing of partial pressure and temperature. Because of 
 *			this, you will get a list of different diffusivities for each data curve. To get the 
 *			optimum kinetic parameters from this list of diffusivities, you must fit the diffusion
 *			parameter values to the following diffusion function model... \n
 *
 *			D_opt = D_ref * exp(-E / (R*T) ) * pow(p , (T_ref/T) - B ) \n
 *			
 *			where D_ref is the Reference Diffusivity (um^2/hr), E is the activation energy for adsorption
 *			(J/mol), R is the gas law constant (J/K/mol), T is the system temperature (K), p is the 
 *			partial pressure of the adsorbing species (kPa), T_ref is the Reference Temperature (K), and
 *			B is the Affinity constant. This algorithm does not automatically produce these parameters
 *			for you, but gives you everything you need to produce them yourself. \n
 *
 *			This routine allows you to optimize multiple kinetic curves at one time. However, all
 *			data must be for the same adsorbent-adsorbate system. In other words, the adsorbent and
 *			adsorbate pair must be the same for each kinetic curve analyzed. Also, each experiment
 *			must have been done in a thin bed or continuous flow system where the adsorbents were
 *			exposed to a nearly constant outside partial pressure for all time steps and the gas
 *			velocity of that system is assumed constant for all experiments. This experimental setup
 *			is very typical for studying adsorption kinetics for gas-solid systems.
 *
 *  \author Austin Ladshaw
 *	\date 05/14/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "scopsowl.h"

#ifndef SCOPSOWL_OPT_HPP_
#define SCOPSOWL_OPT_HPP_

/// Data structure for the SCOPSOWL optmization routine
/** C-style object holding information about the optimization routine as well as the standard
	SCOPSOwl_DATA structure for SCOPSOWL simulations.*/
typedef struct
{
	int num_curves;					///< Number of adsorption curves to analyze
	int evaluation;					///< Number of times the eval function has been called for a single curve
	unsigned long int total_eval;	///< Total number of evaluations needed for completion
	int current_points;				///< Number of points in the current curve
	int num_params = 1;				///< Number of adjustable parameters for the current curve (currently only supports 1)
	int diffusion_type;				///< Flag to identify type of diffusion function to use
	int adsorb_index;				///< Component index for adsorbable species
	int max_guess_iter = 20;		///< Maximum allowed guess iterations (default = 20)
	
	bool Optimize;					///< True = run optimization, False = run a comparison
	bool Rough;						///< True = use only a rough estimate, False = run full optimization
	
	double current_temp;			///< Temperature for current curve
	double current_press;			///< Partial pressure for current curve
	double current_equil;			///< Equilibrium data point for the current curve
	double simulation_equil;		///< Equilibrium simulation point for the current curve
	
	double max_bias;				///< Positive maximum bias plausible for fitting
	double min_bias;				///< Negative minimum bias plausible for fitting
	double e_norm;					///< Euclidean norm of current fit
	double f_bias;					///< Function bias of current fit
	double e_norm_old;				///< Euclidean norm of the previous fit
	double f_bias_old;				///< Function bias of the previous fit
	double param_guess;				///< Parameter guess for the surface/crystal diffusivity
	double param_guess_old;			///< Parameter guess for the previous curve
	double rel_tol_norm = 0.01;		///< Tolerance for convergence of the guess norm
	double abs_tol_bias = 1.0;		///< Tolerance for convergence of the guess bias
	
	std::vector<double> y_base;		///< Gas phase mole fractions in absense of adsorbing species
	
	std::vector<double> q_data;		///< Amount adsorbed at a particular point in current curve
	std::vector<double> q_sim;		///< Amount adsorbed based on the simulation
	std::vector<double> t;			///< Time points in the current curve
	
	FILE *ParamFile;				///< Output file for parameter results
	FILE *CompareFile;				///< Output file for comparison of results
	
	SCOPSOWL_DATA owl_dat;			///< Data structure for the SCOPSOWL simulation
	
}SCOPSOWL_OPT_DATA;

/// Function to set the rest of the gas phase mole fractions based on current mole fraction of adsorbing gas
/** This function takes the current mole fraction of the adsorbing gas and calculates the gas mole fractions of
	the other gases in the sytem based on the standard inlet gas composition given in the scenario file. */
int SCOPSOWL_OPT_set_y(SCOPSOWL_OPT_DATA *owl_opt);

/// Function to set up an initial guess for the surface diffusivity parameter in SCOPSOWL
/** This function performs the Rough optimization on the surface diffusivity based on the idea of reducing or
	eliminating function bias between data and simulation. A positive function bias means that the simulation
	curve is "higher" than the data curve and a negative function bias means that the simulation curve is 
	"lower" than the data curve. We use this information to incrementally adjust the rate of surface diffusion
	until this bias is near zero. When bias is near zero, the simulation is nearly optimized, but further
	refinement may be necessary to find the true minimum solution. */
int initial_guess_SCOPSOWL(SCOPSOWL_OPT_DATA *owl_opt);

/// Function that works in conjunction with the lmfit routine to minimize the euclidean norm between function and data
/** This function will run the SCOPSOWL simulation at a given value of surface diffusivity and produce residuals that
	feed into the Levenberg-Marquardt's algorithm for non-linear least-squares regression. The form of this function is
	specific to the format required by the lmfit routine. 
 
	\param par array of parameters that are to be optimized
	\param m_dat number of data points or functions to evaluate
	\param data user supplied data structure holding information necessary to form the residuals
	\param fvec array of residuals computed at the current parameter values
	\param info integer pointer denoting whether or not the user requests to end a particular simulation*/
void eval_SCOPSOWL_Uptake(const double *par, int m_dat, const void *data, double *fvec, int *info);

/// Function called to perform the optimization routine given a specific set of information and data
/** This is the function that is callable by the UI. The user must provide 5 input files to the routine in order to
	establish simulation conditions, adsborbent properties, component properties, adsorbate equilibrium parameters,
	and the set of data that we are comparing the simulations to. Each input file has a very specific structure and
	order to the information that it contains. The structure here is DIFFERENT than the structure for just running
	standard SCOPSOWL simulations (see scopsowl.h). 
 
	\param scene Sceneario Input File
	\param sorbent Adsorbent Input File
	\param comp Component Input File
	\param sorbate Adsorbate Input File
	\param data Kinetic Adsorption Data File*/

/** \note Much of the structure of these input files are "similar" to that of the input files used in SCOPSOWL_SCENARIOS
	(see scopsowl.h), but with some notable differences. Below gives the format for each input file with an example. Make
	sure your input files follow this format before calling this routine from the UI. \n
 
	Scenario Input File
	-------------------
	Optimization? (0 = false, 1 = true) [tab] Rough Optimization? (0 = false, 1 = true) \n
	Surf. Diff. (0 = constant, 1 = simple Darken, 2 = theoretical Darken) [tab] BC Type (0 = Neumann, 1 = Dirichlet) \n
	Total Pressure (kPa) [tab] Gas Velocity (cm/s) \n
	Number of Gaseous Species \n
	Initial Adsorption Total (mol/kg) \n
	Name [tab] Adsorbable? (0 = false, 1 = true) [tab] Inlet Gas Mole Fraction [tab] Initial Adsorbed Mole Fraction \n
	(NOTE: The above line is repeated for all species in gas phase. Also, this algorithm only allows you to consider 
	one adsorbable gas component. Inlet gas mole fractions must be non-zero for all non-adsorbing gases and must sum
	to 1.) \n
 
	Example Scenario Input
	----------------------
	1	0	\n
	0	0	\n
	101.35	0.36	\n
	5	\n
	0.0	\n
	N2	0	0.7825	0.0	\n
	O2	0	0.2081	0.0	\n
	Ar	0	0.009	0.0	\n
	CO2	0	0.0004	0.0	\n
	H2O	1	0.0		0.0	\n
 
	Above example is for running optimizations on data collected with a gas stream at 0.36 cm/s with 5 gas species in
	the mixture, only H2O of which is adsorbing. The "base line" or "inlet gas" without H2O has a composition of N2 at
	0.7825, O2 at 0.2081, Ar at 0.009, and CO2 at 0.0004. \n
 
	Adsorbent Input File
	--------------------
	Heterogeneous Pellet? (0 = false, 1 = true) \n
	Macro Coord. (2 = spherical, 1 = cylindrical) { [tab] Char. Length (cm) (i.e., cylinder length) } \n
	(NOTE: Char. Length is only needed if problem is not spherical) \n
	Pellet Radius (cm) [tab] Pellet Density (kg/L) [tab] Porosity (vol. void / vol. binder) [tab] Pore Radius (cm) \n
	(Below is only needed if pellet is Heterogeneous) \n
	Micro Coord. (2 = spherical, 1 = cylindrical) { [tab] Char. Length (cm) (i.e., cylinder length) } \n
	Crystal Radius (um) [tab] Binder Fraction (vol. binder / vol. pellet) \n
 
	Example Adsorbent Input
	-----------------------
	1	\n
	2	\n
	0.118	1.69	0.272	3.5E-6	\n
	2	\n
	2.0	0.175 \n
 
	Above example is nearly identical to the file given in the SCOPSOWL_SCENARIO example (see scopsowl.h). However,
	here we do not give an integer flag denoting whether or not we are considering surface diffusion as a mechanism.
	This is because we automatically assume that surface diffusion is a mechanism in the system, since that is the
	unknown parameter that we are performing the optimizations for. \n
 
	Component Input File
	--------------------
	Molar Weight of ith species (g/mol) [tab] Specific Heat of ith species (J/g/K) \n
	Sutherland Viscosity (g/cm/s) [tab] Sutherland Temperature (K) [tab] Sutherland Constant (K) of ith species \n
	(repeat above for all species in same order they appeared in the Scenario Input File) \n
 
	Example Component Input
	-----------------------
	28.016	1.04	\n
	0.0001781	300.55	111.0	\n
	32.0	0.919	\n
	0.0002018	292.25	127.0	\n
	39.948	0.522	\n
	0.0002125	273.11	144.4	\n
	44.009	0.846	\n
	0.000148	293.15	240.0	\n
	18.0	1.97	\n
	0.0001043	298.16	784.72	\n
 
	Above example is exactly the same as in the SCOPSOWL_SCENARIO example (see scopsowl.h). There is no difference
	in the input file formats for this input. Keep in mind that the order is VERY important! All species information
	must be in the same order that the species appeared in the Scenario input file. \n
	
	Adsorbate Input File
	--------------------
	Reference Diffusivity (um^2/hr) [tab] Activation Energy (J/mol) of ith adsorbable species\n
	Reference Temperature (K) [tab] Affinity Constant (-) of ith adsorbable species \n
	van der Waals Volume (cm^3/mol) of ith species \n
	GSTA adsorption capacity (mol/kg) of ith species \n
	Number of GSTA parameters of ith species \n
	Enthalpy (J/mol) of nth site       [tab]       Entropy of nth site (J/K/mol)       of ith species \n
	(repeat enthalpy and entropy for all n sites in species i) \n
	(repeat above for all species i) \n
 
	Example Adsorbate Input
	-----------------------
	0		0	\n
	0		0	\n
	13.91	\n
	11.67	\n
	4	\n
	-46597.5	-53.6994	\n
	-125024		-221.073	\n
	-193619		-356.728	\n
	-272228		-567.459	\n
 
	Above example gives the equilibrium parameters associated with the H2O-MS3A single component adsorption
	system. Note that the kinetic parameters (Ref. Diff., Act. Energy, Ref. Temp., and Affinity) were all
	given a value of zero. These values are irrelavent if we are running an optimization because they will
	be replaced with a single estimate for the diffusivity that is being optimization for. However, if we
	wanted to run this routine with comparisons and not do any optmization, then you would need to provide
	non-zero values for these parameters (at least for Ref. Diff.). \n
 
	Data Input File
	---------------
	Number of Kinetic Data Curves \n
	Number of data points in the ith curve \n
	Temperature (K) [tab]  Partial Pressure (kPa) [tab] Equilibrium Adsorption (mol/kg) all of ith curve \n
	Time point 1 (hrs) [tab] Adsorption 1 (mol/kg) of ith curve \n
	Time point 1 (hrs) [tab] Adsorption 2 (mol/kg) of ith curve \n
	... (Repeat for all time-adsorption data points) \n
	(Repeat above for all curves i) \n
 
	Example Data Input
	------------------
	40 \n
	2990 \n
	298.15	0.000310922	2.9 \n
	0	0 \n
	0.166666667	0.001834419 \n
	0.333611111	0.004880247 \n
	0.5	0.008306803 \n
	... \n
	2789 \n
	298.15	0.00055189	5 \n
	0	0 \n
	0.166944444	0.003350185 \n
	0.333611111	0.007418267 \n
	0.5	0.009930906 \n
	0.666666667	0.014597236 \n
	0.833611111	0.021377373 \n
	.... \n
 
	Above is a partial example for a data set of 40 kinetic curves. The first curve contains 2990 data points 
	and has temperature of 298.15 K, partial pressure of 0.000310922 kPa, and an equilibrium adsorption of 
	2.9. Each first time point should start from 0 hours and each initial adsorption should correspond to the
	value of initial adsorption indicated in the Scenario input file. Then, this structure is repeated for all
	adsorptio curves. \n
 */
int SCOPSOWL_OPTIMIZE(const char *scene, const char *sorbent, const char *comp, const char *sorbate, const char *data);

#endif
