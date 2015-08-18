//----------------------------------------
//  Created by Austin Ladshaw on 1/26/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "finch.h"				//FINCH handles the physics solver and discretization of the problem
#include "magpie.h"				//MAGPIE handles the adsorption equilibria equations
#include "egret.h"				//EGRET handles the parameter estimation for gas phase properties
	
#ifndef SKUA_HPP_
#define SKUA_HPP_

#ifndef D_inf					//Empirical Estimate of infinite dilution diffusivity (um^2/hr)
#define D_inf(Dref,Tref,B,p,T) ( Dref * pow(p+sqrt(DBL_EPSILON),(Tref/T)-B) )
#endif

#ifndef D_o						//Arrhenius Diffusivity (um^2/hr)
#define D_o(Diff,E,T) ( Diff * exp(-E/(Rstd*T)) )
#endif

#ifndef D_c
#define D_c(Diff,phi) ( Diff * (1.0/((1.0+1.1E-6)-phi) ) )
#endif

typedef struct
{
	double activation_energy;				//Activation energy for surface diffusion (J/mol)
	double ref_diffusion;					//Reference state diffusivity (um^2/hr)
	double ref_temperature;					//Reference temperature for empirical adjustments (K)
	double affinity;						//Affinity parameter used in empirical adjustments (-)
	double ref_pressure;					//Reference pressure for empirical adjustments (kPa)
	double film_transfer;					//Film mass transfer coeff (um/hr)
	
	double xIC;								//Inside initial mole fractions of each component (-)
	double y_eff;							//Effective interior gas mole fraction based on adsorption (-)
	
	double Qstn;							//Old heat of adsorption (J/mol)
	double Qstnp1;							//New heat of adsorption (J/mol)
	double xn;								//Old adsorbed mole fraction (-)
	double xnp1;							//New adsorbed mole fraction (-)
	
	bool Adsorbable;						//Boolean to identify with components are adsorbable
	
	std::string speciesName;				//String to hold the name of each species
}SKUA_PARAM;

typedef struct
{
	unsigned long int total_steps;			//Running total of all calculation steps
	int coord;								//Used to determine the coordinates of the problem
	double sim_time;						//Stopping time for the simulation (hrs)
	double t_old;							//Old time of the simulations (hrs)
	double t;								//Current time of the simulations (hrs)
	double t_counter = 0.0;					//Counts for print times for output (hrs)
	double t_print;							//Prints out every t_print time (hrs)
	double qTn;								//Old total amounts adsorbed (mol/kg)
	double qTnp1;							//New total amounts adsorbed (mol/kg)
	bool Print2File = true;					//True = results to .txt; False = no printing
	bool Print2Console = true;				//True = results to console; False = no printing
	
	double gas_velocity;					//Superficial Gas Velocity arount pellet (cm/s)
	double pellet_radius;					//Nominal radius of the pellet/crystal (um)
	double char_measure;					//Length or Area if in Cylindrical or Cartesian coordinates (um or um^2)
	bool DirichletBC = true;				//True = Dirichlet BC; False = Neumann BC
	bool NonLinear = true;					//True = Non-linear solver; False = Linear solver
	std::vector<double> y;					//Outside mole fractions of each component (-)
	
	FILE *OutputFile;						//Output file pointer to the output file for postprocesses
	double (*eval_diff) (int i, int l, const void *user_data);	//Function pointer for evaluating Dc
	double (*eval_kf) (int i, const void *user_data);			//Function pointer for evaluating kf
	const void *user_data;					//Data structure for users info to calculate Dc and kf
	MAGPIE_DATA magpie_dat;					//Data structure for adsorption equilibria
	MIXED_GAS *gas_dat;						//Pointer to the MIXED_GAS data structure (may or may not be used)
	std::vector<FINCH_DATA> finch_dat;		//Data structure for adsorption kinetics
	std::vector<SKUA_PARAM> param_dat;		//Data structure for SKUA specific parameters 
}SKUA_DATA;

void print2file_species_header(FILE *Output, SKUA_DATA *skua_dat, int i);

void print2file_SKUA_time_header(FILE *Output, SKUA_DATA *skua_dat, int i);

void print2file_SKUA_header(SKUA_DATA *skua_dat);

void print2file_SKUA_results_old(SKUA_DATA *skua_dat);

void print2file_SKUA_results_new(SKUA_DATA *skua_dat);

//--------- Default Parameter Functions -------------
double default_Dc(int i, int l, const void *data);

double default_kf(int i, const void *data);
//---------------------------------------------------

//-------- Test Parameter Functions -----------------
double const_Dc(int i, int l, const void *data);

double simple_darken_Dc(int i, int l, const void *data);

//double average_darken_Dc(int i, int l, const void *data);

double theoretical_darken_Dc(int i, int l, const void *data);

//double average_theoretical_darken_Dc(int i, int l, const void *data);

//double const_darken_Dc(int i, int l, const void *data);

double empirical_kf(int i, const void *data);

double const_kf(int i, const void *data);
//---------------------------------------------------

int molefractionCheck(SKUA_DATA *skua_dat);

int setup_SKUA_DATA(FILE *file, double (*eval_Dc) (int i, int l, const void *user_data),
					double (*eval_Kf) (int i, const void *user_data), const void *user_data,
					MIXED_GAS *gas_data, SKUA_DATA *skua_dat);

int SKUA_Executioner(SKUA_DATA *skua_dat);

int set_SKUA_ICs(SKUA_DATA *skua_dat);

int set_SKUA_timestep(SKUA_DATA *skua_dat);

int SKUA_preprocesses(SKUA_DATA *skua_dat);

int set_SKUA_params(const void *user_data);

int SKUA_postprocesses(SKUA_DATA *skua_dat);

int SKUA_reset(SKUA_DATA *skua_dat);

int SKUA(SKUA_DATA *skua_dat);

//-------- Running Specific Tests ------------
int SKUA_CYCLE_TEST01(SKUA_DATA *skua_dat);

int SKUA_CYCLE_TEST02(SKUA_DATA *skua_dat);

int SKUA_LOW_TEST03(SKUA_DATA *skua_dat);

int SKUA_MID_TEST04(SKUA_DATA *skua_dat);
//--------------------------------------------

int SKUA_SCENARIOS(const char *scene, const char *sorbent, const char *comp, const char *sorbate);

int SKUA_TESTS();

#endif
