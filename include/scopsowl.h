//----------------------------------------
//  Created by Austin Ladshaw on 1/29/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "egret.h"							//EGRET handles the estimation of gas-phase properties
#include "skua.h"							//SKUA couples MAGPIE and Surface kinetics together

#ifndef SCOPSOWL_HPP_
#define SCOPSOWL_HPP_

#ifndef Dp
#define Dp(Dm,ep) (ep*ep*Dm)					//Estimate of Pore Diffusivity (cm^2/s)
#endif

#ifndef Dk
#define Dk(rp,T,MW) (9700.0*rp*pow((T/MW),0.5))	//Estimate of Knudsen Diffusivity (cm^2/s)
#endif

#ifndef avgDp
#define avgDp(Dp,Dk) (pow(((1/Dp)+(1/Dk)),-1.0))//Estimate of Average Pore Diffusion (cm^2/s)
#endif

typedef struct
{
	Matrix qAvg;								//Average adsorbed amount for a species at each node (mol/kg)
	Matrix qAvg_old;							//Old Average adsorbed amount for a species at each node (mol/kg)
	
	Matrix Qst;									//Heat of adsorption for all nodes (J/mol)
	Matrix Qst_old;								//Old Heat of adsorption for all nodes (J/mol)
	
	Matrix dq_dc;							//Storage vector for current adsorption slope (dq/dc) (L/kg)
	
	double xIC;									//Initial conditions for adsorbed molefractions
	
	double qIntegralAvg;						//Integral average of adsorption over the entire pellet (mol/kg)
	double qIntegralAvg_old;					//Old Integral average of adsorption over the entire pellet (mol/kg)
	
	double QstAvg;								//Integral average heat of adsorption (J/mol)
	double QstAvg_old;							//Old integral average heat of adsorption (J/mol)
	
	double qo;									//Boundary value of adsorption if using Dirichlet BCs (mol/kg)
	double Qsto;								//Boundary value of adsorption heat if using Dirichlet BCs (J/mol)
	double dq_dco;								//Boundary value of adsorption slope for Dirichelt BCs (L/kg)
	
	double pore_diffusion;					//Value for constant pore diffusion (cm^2/hr)
	double film_transfer;					//Value for constant film mass transfer (cm/hr)
	
	double activation_energy;				//Activation energy for surface diffusion (J/mol)
	double ref_diffusion;					//Reference state diffusivity (um^2/hr)
	double ref_temperature;					//Reference temperature for empirical adjustments (K)
	double affinity;						//Affinity parameter used in empirical adjustments (-)
	double ref_pressure;					//Reference pressure used in empirical adjustments (kPa)
	
	bool Adsorbable;							//True = Adsorbable; False = Non-adsorbable
	
	std::string speciesName;					//String to hold the name of each species
	
}SCOPSOWL_PARAM_DATA;

typedef struct
{
	unsigned long int total_steps;			//Running total of all calculation steps
	int coord_macro;						//Coordinate system for large pellet
	int coord_micro;						//Coordinate system for small crystal (if any)
	int level = 2;							//Level of coupling between the different scales (default = 2)
	double sim_time;						//Stopping time for the simulation (hrs)
	double t_old;							//Old time of the simulations (hrs)
	double t;								//Current time of the simulations (hrs)
	double t_counter = 0.0;					//Counter for the time output
	double t_print;							//Print output at every t_print time (hrs)
	
	bool Print2File = true;					//True = results to .txt; False = no printing
	bool Print2Console = true;				//True = results to console; False = no printing
	bool SurfDiff = true;					//True = includes SKUA; False = No SKUA
	bool Heterogeneous = true;				//True = pellet is made of binder and crystals, False = all one phase
	
	double gas_velocity;					//Superficial Gas Velocity arount pellet (cm/s)
	double total_pressure;					//Gas phase total pressure (kPa)
	double gas_temperature;					//Gas phase temperature (K)
	double pellet_radius;					//Nominal radius of the pellet (cm)
	double crystal_radius;					//Nominal radius of the crystal (um)
	double char_macro;						//Characteristic size for macro scale (cm or cm^2)
	double char_micro;						//Characteristic size for micro scale (um or um^2)
	double binder_fraction;					//Volume of binder per total volume of pellet (-)
	double binder_porosity;					//Volume of pores per volume of binder (-)
	double binder_poresize;					//Nominal radius of the binder pores (cm)
	double pellet_density;					//Mass of the pellet per volume of pellet (kg/L)
	
	bool DirichletBC = false;				//True = Dirichlet BC; False = Neumann BC
	bool NonLinear = true;					//True = Non-linear solver; False = Linear solver
	std::vector<double> y;					//Outside mole fractions of each component (-)
	std::vector<double> tempy;				//Temporary place holder for gas mole fractions in other locations (-)
	
	FILE *OutputFile;											//Output file pointer to the output file for postprocesses
	double (*eval_ads) (int i, int l, const void *user_data);	//Function pointer for evaluating adsorption (mol/kg)
	double (*eval_retard) (int i, int l, const void *user_data);//Function pointer for evaluating retardation (-)
	double (*eval_diff) (int i, int l, const void *user_data);	//Function pointer for evaluating pore diffusion (cm^2/hr)
	double (*eval_surfDiff) (int i, int l, const void *user_data);//Function pointer for evaluating surface diffusion (um^2/hr)
	double (*eval_kf) (int i, const void *user_data);			//Function pointer for evaluating film mass transfer (cm/hr)
	
	const void *user_data;						//Data structure for users info to calculate parameters
	MIXED_GAS *gas_dat;							//Pointer to the MIXED_GAS data structure (may or may not be used)
	MAGPIE_DATA magpie_dat;						//Data structure for a magpie problem (to be used if not using skua)
	std::vector<FINCH_DATA> finch_dat;			//Data structure for pore adsorption kinetics for all species (u in mol/L)
	std::vector<SCOPSOWL_PARAM_DATA> param_dat;	//Data structure for parameter info for all species
	
	std::vector<SKUA_DATA> skua_dat;			//Data structure holding a skua object for all nodes (each skua has an object for each species)
	
}SCOPSOWL_DATA;

void print2file_species_header(FILE *Output, SCOPSOWL_DATA *owl_dat, int i);

void print2file_SCOPSOWL_time_header(FILE *Output, SCOPSOWL_DATA *owl_dat, int i);

void print2file_SCOPSOWL_header(SCOPSOWL_DATA *owl_dat);

void print2file_SCOPSOWL_result_old(SCOPSOWL_DATA *owl_dat);

void print2file_SCOPSOWL_result_new(SCOPSOWL_DATA *owl_dat);

double default_adsorption(int i, int l, const void *user_data);

double default_retardation(int i, int l, const void *user_data);

double default_pore_diffusion(int i, int l, const void *user_data);

double default_surf_diffusion(int i, int l, const void *user_data);

double default_effective_diffusion(int i, int l, const void *user_data);

double const_pore_diffusion(int i, int l, const void *user_data);

double default_filmMassTransfer(int i, const void *user_data);

double const_filmMassTransfer(int i, const void *user_data);

int setup_SCOPSOWL_DATA(FILE *file,
						double (*eval_sorption) (int i, int l, const void *user_data),
						double (*eval_retardation) (int i, int l, const void *user_data),
						double (*eval_pore_diff) (int i, int l, const void *user_data),
						double (*eval_filmMT) (int i, const void *user_data),
						double (*eval_surface_diff) (int i, int l, const void *user_data),
						const void *user_data,MIXED_GAS *gas_data,SCOPSOWL_DATA *owl_data);

int SCOPSOWL_Executioner(SCOPSOWL_DATA *owl_dat);

int set_SCOPSOWL_ICs(SCOPSOWL_DATA *owl_dat);

int set_SCOPSOWL_timestep(SCOPSOWL_DATA *owl_dat);

int SCOPSOWL_preprocesses(SCOPSOWL_DATA *owl_dat);

int set_SCOPSOWL_params(const void *user_data);

int SCOPSOWL_postprocesses(SCOPSOWL_DATA *owl_dat);

int SCOPSOWL_reset(SCOPSOWL_DATA *owl_dat);

int SCOPSOWL(SCOPSOWL_DATA *owl_dat);

int LARGE_CYCLE_TEST01(SCOPSOWL_DATA *owl_dat);

int SMALL_CYCLE_TEST02(SCOPSOWL_DATA *owl_dat);

int CURVE_TEST03(SCOPSOWL_DATA *owl_dat);

int CURVE_TEST04(SCOPSOWL_DATA *owl_dat);

int CURVE_TEST05(SCOPSOWL_DATA *owl_dat);

int SCOPSOWL_SCENARIOS(const char *scene, const char *sorbent, const char *comp, const char *sorbate);

int SCOPSOWL_TESTS();

#endif
