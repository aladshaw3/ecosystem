//----------------------------------------
//  Created by Austin Ladshaw on 6/10/14
//  Copyright (c) 2014
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*!
 *  \file finch.h finch.cpp
 *	\brief Flux-limiting Implicit Non-oscillatory Conservative High-resolution scheme
 *	\details This is a conservative finite differences scheme based on the Kurganov and Tadmoor (2000)
 *		MUSCL scheme for non-linear conservation laws. It can solve 1-D conservation law problems in
 *		three different coordinate systems: (i) Cartesian - axial, (ii) Cylindrical - radial, and (iii)
 *		Spherical - radial. It is the backbone algorithm behind all 1-D PDE problems in the ecosystem
 *		software.
 *
 *		The form of the general conservation law problem that FINCH solves is...
 *
 *		z^d*R*du/dt = d/dz(z^d*D*du/dz) - d/dz(z^d*v*u) - z^d*k*u + z^d*S
 *
 *		where R, D, v, k, and S are the parameters of the problem and d, z, and u are the coordinates,
 *		spatial dimension, and conserved quantities, respectively. The parameter R is a retardation
 *		coefficient, D is a diffusion coefficient, v is a velocity, k is a reaction coefficient, and S
 *		is a forcing function or source/sink term. 
 *
 *		FINCH supports the use of both Dirichlet and Neuman boundary conditions as the input/inlet
 *		condition and uses the No Flux (or Natural) boundary condition for the output/outlet of the
 *		domain. For radial problems, the outlet is always taken to the the center of the cylindrical
 *		or spherical particle. This enforces the symmetry of the problem. For axial problems, the outlet
 *		is determined by the sign of the velocity term and is therefore choosen by the routine based on
 *		the actual flow direction in the domain. 
 *
 *		Parameters of the problem can be coupled to the variable u and also be functions of space and
 *		time. The coupling of the parameters with the variable forces the problem to become non-linear,
 *		which requires iteration to solve. The default iterative method is a built-in Picard's method.
 *		This method is equivalent to an inexact Newton method, because we use the Linear Solve of this
 *		system as a weak approximation to the non-linear solve. Generally, this method is sufficient and
 *		is the most efficient. However, if a problem is particularly difficult to solve, then we can call
 *		some of the non-linear solvers developed in LARK. If PJFNK is used, then the Linear Solve for
 *		the FINCH problem is used as the Preconditioner for the Linear Solve in PJFNK.
 *
 *		This algorithm comes packaged with three different slope limiter functions to stabilize the 
 *		velocity term for highly advectively dominate problems. The available slope limiters are: (i)
 *		minmod, (ii) van Albada, and (iii) ospre. By default, the FINCH setup function will set the
 *		slope limiter to ospre, because this method provides a reasonable compromise between accuracy
 *		and efficiency. 
 *		
 *		Slope Limiter Stats:
 *		--------------------
 *		minmod -> Highest Accuracy, Lowest Efficiency \n
 *		van Albada -> Lowest Accuracy, Highest Efficiency \n
 *		ospre -> Average Accuracy, Average Efficiency \n
 *
 *  \author Austin Ladshaw
 *	\date 01/29/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#ifndef FINCH_HPP_
#define FINCH_HPP_

#include "macaw.h"
#include "lark.h"

/// List of enum options to define the solver type in FINCH
typedef enum { FINCH_Picard, LARK_Picard, LARK_PJFNK } finch_solve_type;

/// List of enum options to define the coordinate system in FINCH
typedef enum { Cartesian, Cylindrical, Spherical } finch_coord_type;

/// Data structure for the FINCH object
/** C-style object that holds data, functions, and other structures necessary to discretize
	and solve a FINCH problem. All of this information must be overriden or initialized prior
	to running a FINCH simulation. Many, many default functions are provided to make it easier
	to incorporate FINCH into other problems. The main function to override will be the setparams
	function. This will be a function that the user provides to tell the FINCH simulation how
	the parameters of the problem vary in time and space and whether or not they are coupled
	the the variable u. All functions are overridable and several can be skipped entirely, or
	called directly at different times in the execution of a particular routine. This make FINCH
	extremely flexible to the user. */

/** \note All parameters and dimensions do not carry any units with them. The user is required to
	keep track of all their own units in their particular problem and ensure that units will cancel
	and be consistent in their own physical model. */
typedef struct
{
	//All parameters are given default values prior to build and run
  	int d = 0;				///< Dimension of the problem: 0 = cartesian, 1 = cylindrical, 2 = spherical
	double dt = 0.0125;		///< Time step
	double dt_old = 0.0125;	///< Previous time step
	double T = 1.0;			///< Total time
	double dz = 0.1;		///< Space step
	double L = 1.0;			///< Total space
  	double s = 1.0;			///< Char quantity (spherical = 1, cylindrical = length, cartesian = area)
	double t = 0.0;			///< Current Time
	double t_old = 0.0;		///< Previous Time
	
  	double uT = 0.0;		///< Total amount of conserved quantity in domain
  	double uT_old = 0.0;	///< Old Total amount of conserved quantity
  	double uAvg = 0.0;		///< Average amount of conserved quantity in domain
  	double uAvg_old = 0.0;	///< Old Average amount of conserved quantity
	double uIC = 0.0;		///< Initial condition of Conserved Quantity
	double vIC = 1.0;		///< Initial condition of Velocity
	double DIC = 1.0;		///< Initial condition of Dispersion
  	double kIC = 1.0;		///< Initial condition of Reaction
  	double RIC = 1.0;		///< Initial condition of the Time Coefficient
	double uo = 1.0;		///< Boundary Value of Conserved Quantity
	double vo = 1.0;		///< Boundary Value of Velocity
	double Do = 1.0;		///< Boundary Value of Dispersion
  	double ko = 1.0;		///< Boundary Value of Reaction
  	double Ro = 1.0;		///< Boundary Value of Time Coefficient
	double kfn = 1.0;		///< Film mass transfer coefficient Old
  	double kfnp1 = 1.0;		///< Film mass transfer coefficient New
  	double lambda_I;        ///< Boundary Coefficient for Implicit Neumann (Calculated at Runtime)
  	double lambda_E;        ///< Boundary Coefficient for Explicit Neumann (Calculated at Runtime)
	
	int LN = 10;				///< Number of nodes
	bool CN = true;				///< True if Crank-Nicholson, false if Implicit, never use explicit
	bool Update = false;		///< Flag to check if the system needs updating
	bool Dirichlet = false;		///< Flag to indicate use of Dirichlet or Neumann starting boundary
  	bool CheckMass = false;		///< Flag to indicate whether or not mass is to be checked
  	bool ExplicitFlux = false;	///< Flag to indicate whether or not to use fully explicit flux limiters
  	bool Iterative = true;		///< Flag to indicate whether to solve directly, or iteratively
  	bool SteadyState = false;	///< Flag to determine whether or not to solve the steady-state problem
  	bool NormTrack = true;		///< Flag to determine whether or not to track the norms during simulation
	
	double beta = 0.5;			///< Scheme type indicator: 0.5=CN & 1.0=Implicit; all else NULL
  	double tol_rel = 1e-6;		///< Relative Tolerance for Convergence
	double tol_abs = 1e-6;		///< Absolute Tolerance for Convergence
  	int max_iter = 20;			///< Maximum number of iterations allowed
  	int total_iter = 0;			///< Total number of iterations made
	int nl_method = FINCH_Picard;///< Non-linear solution method - default = FINCH_Picard
	
	//Coefficients For Semi-Discrete Form
	std::vector<double> CL_I;	///< Left side, implicit coefficients (Calculated at Runtime)
	std::vector<double> CL_E;	///< Left side, explicit coefficients (Calculated at Runtime)
	std::vector<double> CC_I;	///< Centered, implicit coefficients (Calculated at Runtime)
	std::vector<double> CC_E;	///< Centered, explicit coefficients (Calculated at Runtime)
	std::vector<double> CR_I;	///< Right side, implicit coefficients (Calculated at Runtime)
	std::vector<double> CR_E;	///< Right side, explicit coefficients (Calculated at Runtime)
	//Fluxes for Semi-Discrete Form
	std::vector<double> fL_I;	///< Left side, implicit fluxes (Calculated at Runtime)
	std::vector<double> fL_E;	///< Left side, explicit fluxes (Calculated at Runtime)
	std::vector<double> fC_I;	///< Centered, implicit fluxes (Calculated at Runtime)
	std::vector<double> fC_E;	///< Centered, explicit fluxes (Calculated at Runtime)
	std::vector<double> fR_I;	///< Right side, implicit fluxes (Calculated at Runtime)
	std::vector<double> fR_E;	///< Right side, explicit fluxes (Calculated at Runtime)
	
	//For CN Fully Discrete Form
	std::vector<double> OI;		///< Implicit upper diagonal matrix elements (Calculated at Runtime)
	std::vector<double> OE;		///< Explicit upper diagonal matrix elements (Calculated at Runtime)
	std::vector<double> NI;		///< Implicit diagonal matrix elements (Calculated at Runtime)
	std::vector<double> NE;		///< Explicit diagonal matrix elements (Calculated at Runtime)
	std::vector<double> MI;		///< Implicit lower diagonal matrix elements (Calculated at Runtime)
	std::vector<double> ME;		///< Explicit lower diagonal matrix elements (Calculated at Runtime)
  
  	//Small Vectors used in Slope Reconstructions (MAX Size = 3)
  	std::vector<double> uz_l_I, uz_lm1_I, uz_lp1_I;	///< Implicit local slopes (Calculated at Runtime)
  	std::vector<double> uz_l_E, uz_lm1_E, uz_lp1_E;	///< Explicit local slopes (Calculated at Runtime)
	
	//comp
	Matrix<double> unm1;			///< Conserved Quantity Older
	Matrix<double> un;				///< Conserved Quantity Old
	Matrix<double> unp1;			///< Conserved Quantity New
	Matrix<double> u_star;			///< Conserved Quantity Projected New
	Matrix<double> ubest;			///< Best found solution if solving iteratively
	//sys
	Matrix<double> vn;				///< Velocity Old
	Matrix<double> vnp1;			///< Velocity New
	Matrix<double> Dn;				///< Dispersion Old
	Matrix<double> Dnp1;			///< Dispersion New
  	Matrix<double> kn;				///< Reaction Old
  	Matrix<double> knp1;			///< Reaction New
	Matrix<double> Sn;				///< Forcing Function Old
	Matrix<double> Snp1;			///< Forcing Function New
  	Matrix<double> Rn;				///< Time Coeff Old
  	Matrix<double> Rnp1;			///< Time Coeff New
	//comp
	Matrix<double> Fn;				///< Flux Limiter Old
	Matrix<double> Fnp1;			///< Flux Limiter New
	Matrix<double> gI;				///< Implicit Side Boundary Conditions
	Matrix<double> gE;				///< Explicit Side Boundary Conditions
  	//Iteration Info
  	Matrix<double> res;				///< Current residual
  	Matrix<double> pres;			///< Current search direction
	
	//Function pointers for user defined functions (defaults are provided)
	//NOTE: if defaults are used, user_data must be a FINCH_DATA object
	//Otherwise, user can use any data structure that contains a FINCH_DATA object
	int (*callroutine) (const void *user_data); ///< Function pointer to executioner (DEFAULT = default_execution)
	int (*setic) (const void *user_data);		///< Function pointer to initial conditions (DEFAULT = default_ic)
	int (*settime) (const void *user_data);		///< Function pointer to set time step (DEFAULT = default_timestep)
	int (*setpreprocess) (const void *user_data);///< Function pointer to preprocesses (DEFAULT = default_preprocess)
	int (*solve) (const void *user_data);		///< Function pointer to the solver (DEFAULT = default_solve)
	int (*setparams) (const void *user_data);	///< Function pointer to set parameters (DEFAULT = default_params)
	int (*discretize) (const void *user_data);	///< Function pointer to discretization (DEFAULT = ospre_discretization)
	int (*setbcs) (const void *user_data);		///< Function pointer to set boundary conditions (DEFAULT = default_bcs)
	/// Function pointer to the residual function (DEFAULT = default_res)
	int (*evalres) (const Matrix<double>& x, Matrix<double>& res, const void *user_data);
	/// Function pointer to the preconditioning function (DEFAULT = default_precon)
	int (*evalprecon) (const Matrix<double>& b, Matrix<double>& p, const void *user_data);
	int (*setpostprocess) (const void *user_data);	///< Function pointer to the postprocesses (DEFAULT = default_postprocess)
	int (*resettime) (const void *user_data);		///< Function pointer to reset time (DEFAULT = default_reset)
	
	//LARK data structures
	PICARD_DATA picard_dat;			///< Data structure for PICARD method (no need to use this)
	PJFNK_DATA pjfnk_dat;			///< Data structure for PJFNK method (more rigours method)
	const void *param_data;		///< User's data structure used to evaluate the parameter function (Must override if setparams is overriden)
	
}FINCH_DATA;

/// Function returns the maximum in a list of values
double max(std::vector<double> &values);

/// Function returns the minimum in a list of values
double min(std::vector<double> &values);

/// Function returns the result of the minmod function acting on a list of values
double minmod(std::vector<double> &values);

/// Function integrates the conserved quantity to return it's total in the domain
int uTotal(FINCH_DATA *dat);

/// Function integrates the conserved quantity to reture it's average in the domain
int uAverage(FINCH_DATA *dat);

/// Function checks the unp1 vector for negative values and will adjust if needed
/** This function can be turned off or on in the FINCH_DATA structure. Typically, you
	will want to leave this on so that the routine does not return negative values for u.
	However, if you want to get negative values of u, then turn this option off. */
int check_Mass(FINCH_DATA *dat);

/// Function solves the discretized FINCH problem directly by assuming it is linear
int l_direct(FINCH_DATA *dat);

/// Function to perform the necessary LARK Picard iterative method (not typically used)
int lark_picard_step(const Matrix<double> &x, Matrix<double> &G, const void *data);

/// Function to solve the discretized FINCH problem iteratively by assuming it is non-linear
/** \note If the problem is actually linear, then this will solve it in one iteration. So it
	may be best to always assume the problem is non-linear. */
int nl_picard(FINCH_DATA *dat);

/// Function to setup memory and set user defined functions into the FINCH object
/** This function MUST be called prior to running any FINCH based simulation. However, you are
	only every required to provide this function with the FINCH_DATA pointer. It is recommended,
	however, that you do provide the user_setparams and param_data pointers, as these will likely
	vary significantly from problem to problem. 
 
	After the problem is setup in memory, you do not technically have to have FINCH call all of
	it's own functions. You can write your own executioner, initial conditions, and other functions
	and decided how and when everything is called. Then just call the solve function in FINCH_DATA
	when you want to use the FINCH solver. This is how FINCH is used in SKUA, SCOPSOWL, DOGFISH, and
	MONKFISH. */
int setup_FINCH_DATA( int (*user_callroutine) (const void *user_data),
					  int (*user_setic) (const void *user_data),
					  int (*user_timestep) (const void *user_data),
					  int (*user_preprocess) (const void *user_data),
					  int (*user_solve) (const void *user_data),
					  int (*user_setparams) (const void *user_data),
					  int (*user_discretize) (const void *user_data),
					  int (*user_bcs) (const void *user_data),
					  int (*user_res) (const Matrix<double>&x, Matrix<double>& res, const void *user_data),
					  int (*user_precon) (const Matrix<double>&b, Matrix<double>& p, const void *user_data),
					  int (*user_postprocess) (const void *user_data),
					  int (*user_reset) (const void *user_data),
					  FINCH_DATA *dat, const void *param_data);

/// Function will print out a dimension header for FINCH output
void print2file_dim_header(FILE *Output, FINCH_DATA *dat);

/// Function will print out a time header for FINCH output
void print2file_time_header(FILE *Output, FINCH_DATA *dat);

/// Function will print out the old results to the variable u
void print2file_result_old(FILE *Output, FINCH_DATA *dat);

/// Function will print out the new results to the variable u
void print2file_result_new(FILE *Output, FINCH_DATA *dat);

/// Function will force print out a blank line
void print2file_newline(FILE *Output, FINCH_DATA *dat);

/// Function will force print out a tab
void print2file_tab(FILE *Output, FINCH_DATA *dat);

//Default Functions in FINCH--------------------------------------------

/// Default executioner function for FINCH
int default_execution(const void *user_data);

int default_ic(const void *user_data);

int default_timestep(const void *user_data);

int default_preprocess(const void *user_data);

int default_solve(const void *user_data);

int default_params(const void *user_data);

int minmod_discretization(const void *user_data);

int vanAlbada_discretization(const void *user_data);

int ospre_discretization(const void *user_data);

int default_bcs(const void *user_data);

int default_res(const Matrix<double> &x, Matrix<double> &res, const void *user_data);

int default_precon(const Matrix<double> &b, Matrix<double> &p, const void *user_data);

int default_postprocess(const void *user_data);

int default_reset(const void *user_data);
//END Default Functions------------------------------------------------

//Specific Test Functions----------------------------------------------
/// \cond
int buckley_leverett_ic(const void *user_data);

int buckley_leverett_params(const void *user_data);

int burgers_ic(const void *user_data);

int burgers_params(const void *user_data);

int burgers_bcs(const void *user_data);
/// \endcond
//END Specific Functions-----------------------------------------------

int FINCH_TESTS();

#endif
