//----------------------------------------
//  Created by Austin Ladshaw on 4/14/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef MONKFISH_HPP_
#define MONKFISH_HPP_

#include "dogfish.h"

//Data structure for species specific info
typedef struct
{
	double interparticle_diffusion;			//Units: cm^2/hr
	double exterior_concentration;			//Units: mol/L
	double exterior_transfer_coeff;			//Units: cm/hr
	double sorbed_molefraction;				//Units: -
	double initial_sorption;				//Units: mg/g
	double sorption_bc;						//Units: mg/g
	
	double intraparticle_diffusion;			//Units: um^2/hr
	double film_transfer_coeff;				//Units: um/hr
	
	Matrix avg_sorption;					//Units: mg/g
	Matrix avg_sorption_old;				//Units: mg/g
	
	Molecule species;						//Species in the liquid phase
	
}MONKFISH_PARAM;

//Primary data structure for running MONKFISH
typedef struct
{
	unsigned long int total_steps = 0;
	double time_old = 0.0;
	double time = 0.0;
	bool Print2File = true;					//True = results to .txt; False = no printing
	bool Print2Console = true;				//True = results to console; False = no printing
	bool DirichletBC = true;				//False = uses film mass transfer for BC, True = Dirichlet BC
	bool NonLinear = false;					//False = Solve directly, True = Solve iteratively
	bool haveMinMax = false;				//True = know min and max fiber density, False = only know avg density (Used in ICs)
	bool MultiScale = true;					//True = solve single fiber model at nodes, False = solve equilibrium at nodes
	int level = 2;							//Level of coupling between multiple scales (default = 2)
	double t_counter = 0.0;					//Counter for the time output
	double t_print;							//Print output at every t_print time (hrs)
	
	int NumComp;							//Number of species to track
	double end_time;						//Units: hours
	double total_sorption_old;				//Per mass or volume of woven nest
	double total_sorption;					//Per mass or volume of woven nest
	double single_fiber_density;			//Units: g/L
	double avg_fiber_density;				//Units: g/L (Used in ICs)
	double max_fiber_density;				//Units: g/L (Used in ICs)
	double min_fiber_density;				//Units: g/L (Used in ICs)
	double max_porosity;					//Units: -
	double min_porosity;					//Units: -
	double domain_diameter;					//Units: cm
	
	FILE *Output;							//Output file pointer
	//Function pointers for the parameters to be evaluated in FINCH for MONKFISH
	double (*eval_eps) (int i, int l, const void *user_data);
	double (*eval_rho) (int i, int l, const void *user_data);
	double (*eval_Dex) (int i, int l, const void *user_data);
	double (*eval_ads) (int i, int l, const void *user_data);
	double (*eval_Ret) (int i, int l, const void *user_data);
	double (*eval_Cex) (int i, const void *user_data);
	double (*eval_kf) (int i, const void *user_data);
	
	const void *user_data;					//User supplied data function
	std::vector<FINCH_DATA> finch_dat;		//FINCH data structure (species)
	std::vector<MONKFISH_PARAM> param_dat;	//MONKFISH parameter data (species)
	std::vector<DOGFISH_DATA> dog_dat;		//DOGFISH data structure (nodal)
	
}MONKFISH_DATA;

//--------Default Parameter Functions---------------
double default_porosity(int i, int l, const void *user_data);

double default_density(int i, int l, const void *user_data);

double default_interparticle_diffusion(int i, int l, const void *user_data);

double default_monk_adsorption(int i, int l, const void *user_data);

double default_monk_equilibrium(int i, int l, const void *user_data);

double default_monkfish_retardation(int i, int l, const void *user_data);

double default_exterior_concentration(int i, const void *user_data);

double default_film_transfer(int i, const void *user_data);
//--------End Default Parameter Functions-----------

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



int MONKFISH_TESTS();

#endif
