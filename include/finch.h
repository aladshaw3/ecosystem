//----------------------------------------
//  Created by Austin Ladshaw on 6/10/14
//  Copyright (c) 2014
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef FINCH_HPP_
#define FINCH_HPP_

#include "macaw.h"
#include "lark.h"

//List of macros to define the solver type in FINCH
#define FINCH_Picard 0
#define LARK_Picard 1
#define LARK_PJFNK 2
#define Cartesian 0
#define Cylindrical 1
#define Spherical 2

typedef struct
{
	//All parameters are given default values prior to build and run
  	int d = 0;				//Dimension of the problem: 0 = cartesian, 1 = cylindrical, 2 = spherical
	double dt = 0.0125;		//Time step
	double dt_old = 0.0125;	//Previous time step
	double T = 1.0;			//Total time
	double dz = 0.1;		//Space step
	double L = 1.0;			//Total space
  	double s = 1.0;			//Char quantity (spherical = 1, cylindrical = length, cartesian = area)
	double t = 0.0;			//Current Time
	double t_old = 0.0;		//Previous Time
	
  	double uT = 0.0;		//Total amount of conserved quantity in domain
  	double uT_old = 0.0;	//Old Total amount of conserved quantity
  	double uAvg = 0.0;		//Average amount of conserved quantity in domain
  	double uAvg_old = 0.0;	//Old Average amount of conserved quantity
	double uIC = 0.0;		//Initial condition of Conserved Quantity
	double vIC = 1.0;		//Initial condition of Velocity
	double DIC = 1.0;		//Initial condition of Dispersion
  	double kIC = 1.0;		//Initial condition of Reaction
  	double RIC = 1.0;		//Initial condition of the Time Coefficient
	double uo = 1.0;		//Boundary Value of Conserved Quantity
	double vo = 1.0;		//Boundary Value of Velocity
	double Do = 1.0;		//Boundary Value of Dispersion
  	double ko = 1.0;		//Boundary Value of Reaction
  	double Ro = 1.0;		//Boundary Value of Time Coefficient
	double kfn = 1.0;		//Film mass transfer coefficient Old
  	double kfnp1 = 1.0;		//Film mass transfer coefficient New
  	double lambda_I;        //Boundary Coefficient for Implicit Neumann (Calculated at Runtime)
  	double lambda_E;        //Boundary Coefficient for Explicit Neumann (Calculated at Runtime)
	
	int LN = 10;				//Number of nodes
	bool CN = true;				//True if Crank-Nicholson, false if Implicit, never use explicit
	bool Update = false;		//Flag to check if the system needs updating
	bool Dirichlet = false;		//Flag to indicate use of Dirichlet or Neumann starting boundary
  	bool CheckMass = false;		//Flag to indicate whether or not mass is to be checked
  	bool ExplicitFlux = false;	//Flag to indicate whether or not to use fully explicit flux limiters
  	bool Iterative = true;		//Flag to indicate whether to solve directly, or iteratively
  	bool SteadyState = false;	//Flag to determine whether or not to solve the steady-state problem
  	bool NormTrack = true;		//Flag to determine whether or not to track the norms during simulation
	
	double beta = 0.5;			//Scheme type indicator: 0.5=CN & 1.0=Implicit; all else NULL
  	double tol_rel = 1e-6;		//Relative Tolerance for Convergence
	double tol_abs = 1e-6;		//Absolute Tolerance for Convergence
  	int max_iter = 20;			//Maximum number of iterations allowed
  	int total_iter = 0;			//Total number of iterations made
	int nl_method = FINCH_Picard;//Non-linear solution method - default = FINCH_Picard
	
	//Coefficients For Semi-Discrete Form (fill out in Parameter function)
	std::vector<double> CL_I;
	std::vector<double> CL_E;
	std::vector<double> CC_I;
	std::vector<double> CC_E;
	std::vector<double> CR_I;
	std::vector<double> CR_E;
	//Fluxes for Semi-Discrete Form (fill out in Parameter function)
	std::vector<double> fL_I;
	std::vector<double> fL_E;
	std::vector<double> fC_I;
	std::vector<double> fC_E;
	std::vector<double> fR_I;
	std::vector<double> fR_E;
	
	//For CN Fully Discrete Form
	std::vector<double> OI;		//Implicit l+1 coefficients
	std::vector<double> OE;		//Explicit l+1 coefficients
	std::vector<double> NI;		//Implicit l coefficients
	std::vector<double> NE;		//Explicit l coefficients
	std::vector<double> MI;		//Implicit l-1 coefficients
	std::vector<double> ME;		//Explicit l-1 coefficients
  
  	//Small Vectors used in Slope Reconstructions (MAX Size = 3)
  	std::vector<double> uz_l_I, uz_lm1_I, uz_lp1_I;
  	std::vector<double> uz_l_E, uz_lm1_E, uz_lp1_E;
	
	//comp
	Matrix unm1;			//Conserved Quantity Older
	Matrix un;				//Conserved Quantity Old
	Matrix unp1;			//Conserved Quantity New
	Matrix u_star;			//Conserved Quantity Projected New
	Matrix ubest;			//Best found solution if solving iteratively
	//sys
	Matrix vn;				//Velocity Old
	Matrix vnp1;			//Velocity New
	Matrix Dn;				//Dispersion Old
	Matrix Dnp1;			//Dispersion New
  	Matrix kn;				//Reaction Old
  	Matrix knp1;			//Reaction New
	Matrix Sn;				//Forcing Function Old
	Matrix Snp1;			//Forcing Function New
  	Matrix Rn;				//Time Coeff Old
  	Matrix Rnp1;			//Time Coeff New
	//comp
	Matrix Fn;				//Flux Limiter Old
	Matrix Fnp1;			//Flux Limiter New
	Matrix gI;				//Implicit Side BC
	Matrix gE;				//Explicit Side BC
  	//Iteration Info
  	Matrix res;				//Current residual
  	Matrix pres;			//Current search direction
	
	//Function pointers for user defined functions (defaults are provided)
	//NOTE: if defaults are used, user_data must be a FINCH_DATA object
	//Otherwise, user can use any data structure that contains a FINCH_DATA object
	int (*callroutine) (const void *user_data); 
	int (*setic) (const void *user_data);
	int (*settime) (const void *user_data);
	int (*setpreprocess) (const void *user_data);
	int (*solve) (const void *user_data);
	int (*setparams) (const void *user_data);
	int (*discretize) (const void *user_data);
	int (*setbcs) (const void *user_data);
	int (*evalres) (const Matrix& x, Matrix& res, const void *user_data);
	int (*evalprecon) (const Matrix& b, Matrix& p, const void *user_data);
	int (*setpostprocess) (const void *user_data);
	int (*resettime) (const void *user_data);
	
	//LARK data structures
	PICARD_DATA picard_dat;			//Only used in Testing (better version available in nl_picard)
	PJFNK_DATA pjfnk_dat;			//Used for Newton-Krylov iterative methods
	const void *param_data;			//This is the user data structure used to evaluate the parameter function
	
}FINCH_DATA;

double max(std::vector<double> &values);

double min(std::vector<double> &values);

double minmod(std::vector<double> &values);

int uTotal(FINCH_DATA *dat);

int uAverage(FINCH_DATA *dat);

int check_Mass(FINCH_DATA *dat);

int l_direct(FINCH_DATA *dat);

int lark_picard_step(const Matrix &x, Matrix &G, const void *data);

int nl_picard(FINCH_DATA *dat);

//Function to setup memory and set user defined functions into the FINCH object
int setup_FINCH_DATA( int (*user_callroutine) (const void *user_data),
					  int (*user_setic) (const void *user_data),
					  int (*user_timestep) (const void *user_data),
					  int (*user_preprocess) (const void *user_data),
					  int (*user_solve) (const void *user_data),
					  int (*user_setparams) (const void *user_data),
					  int (*user_discretize) (const void *user_data),
					  int (*user_bcs) (const void *user_data),
					  int (*user_res) (const Matrix&x, Matrix& res, const void *user_data),
					  int (*user_precon) (const Matrix&b, Matrix& p, const void *user_data),
					  int (*user_postprocess) (const void *user_data),
					  int (*user_reset) (const void *user_data),
					  FINCH_DATA *dat, const void *param_data);

void print2file_dim_header(FILE *Output, FINCH_DATA *dat);

void print2file_time_header(FILE *Output, FINCH_DATA *dat);

void print2file_result_old(FILE *Output, FINCH_DATA *dat);

void print2file_result_new(FILE *Output, FINCH_DATA *dat);

void print2file_newline(FILE *Output, FINCH_DATA *dat);

void print2file_tab(FILE *Output, FINCH_DATA *dat);

//Default Functions in FINCH--------------------------------------------
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

int default_res(const Matrix &x, Matrix &res, const void *user_data);

int default_precon(const Matrix &b, Matrix &p, const void *user_data);

int default_postprocess(const void *user_data);

int default_reset(const void *user_data);
//END Default Functions------------------------------------------------

//Specific Test Functions----------------------------------------------
int buckley_leverett_ic(const void *user_data);

int buckley_leverett_params(const void *user_data);

int burgers_ic(const void *user_data);

int burgers_params(const void *user_data);

int burgers_bcs(const void *user_data);
//END Specific Functions-----------------------------------------------

int FINCH_TESTS();

#endif
