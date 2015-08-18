//----------------------------------------
//  Created by Austin Ladshaw on 5/11/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "skua.h"

#ifndef SKUA_OPT_HPP_
#define SKUA_OPT_HPP_

typedef struct
{
	int num_curves;					//Number of adsorption curves to analyze
	int evaluation;					//Number of times the eval function has been called for a single curve
	unsigned long int total_eval;	//Total number of evaluations needed for completion
	int current_points;				//Number of points in the current curve
	int num_params = 1;					//Number of adjustable parameters for the current curve
	int diffusion_type;				//Flag to identify type of diffusion function to use
	int adsorb_index;				//Component index for adsorbable species
	int max_guess_iter = 20;		//Maximum allowed guess iterations (default = 20)
	
	bool Optimize;					//True = run optimization, False = run a comparison
	bool Rough;						//True = use only a rough estimate, False = run full optimization
	
	double current_temp;			//Temperature for current curve
	double current_press;			//Partial pressure for current curve
	double current_equil;			//Equilibrium data point for the current curve
	double simulation_equil;		//Equilibrium simulation point for the current curve
	
	double max_bias;				//Positive maximum bias plausible for fitting
	double min_bias;				//Negative minimum bias plausible for fitting
	double e_norm;					//Euclidean norm of current fit
	double f_bias;					//Function bias of current fit
	double e_norm_old;				//Euclidean norm of the previous fit
	double f_bias_old;				//Function bias of the previous fit
	double param_guess;				//Parameter guess for the surface/crystal diffusivity
	double param_guess_old;			//Parameter guess for the previous curve
	double rel_tol_norm = 0.1;		//Tolerance for convergence of the guess norm
	double abs_tol_bias = 0.1;		//Tolerance for convergence of the guess bias
	
	std::vector<double> y_base;		//Gas phase mole fractions in absense of adsorbing species
	
	std::vector<double> q_data;		//Amount adsorbed at a particular point in current curve
	std::vector<double> q_sim;		//Amount adsorbed based on the simulation
	std::vector<double> t;			//Time points in the current curve
	
	FILE *ParamFile;				//Output file for parameter results
	FILE *CompareFile;				//Output file for comparison of results
	
	SKUA_DATA skua_dat;				//Data structure for the SKUA simulation
	
}SKUA_OPT_DATA;

int SKUA_OPT_set_y(SKUA_OPT_DATA *skua_opt);

int initial_guess_SKUA(SKUA_OPT_DATA *skua_opt);

void eval_SKUA_Uptake(const double *par, int m_dat, const void *data, double *fvec, int *info);

int SKUA_OPTIMIZE(const char *scene, const char *sorbent, const char *comp, const char *sorbate, const char *data);

#endif
