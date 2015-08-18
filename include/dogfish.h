//----------------------------------------
//  Created by Austin Ladshaw on 4/9/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef DOGFISH_HPP_
#define DOGFISH_HPP_

#include "finch.h"
#include "mola.h"

//Data structure for species-specific parameters
typedef struct
{
	//Default parameters to go with default functions
	double intraparticle_diffusion;						//Units: um^2/hr
	double film_transfer_coeff;							//Units: um/hr
	double surface_concentration;						//Units: mg/g
	double initial_sorption;							//Units: mg/g
	
	//Additional info
	double sorbed_molefraction;
	Molecule species;									//Adsorbed species
	
}DOGFISH_PARAM;

//Primary data structure for running the DOGFISH application
typedef struct
{
	unsigned long int total_steps = 0;
	double time_old = 0.0;
	double time = 0.0;
	bool Print2File = true;					//True = results to .txt; False = no printing
	bool Print2Console = true;				//True = results to console; False = no printing
	bool DirichletBC = false;				//False = uses film mass transfer for BC, True = Dirichlet BC
	bool NonLinear = false;					//False = Solve directly, True = Solve iteratively
	double t_counter = 0.0;					//Counter for the time output
	double t_print;							//Print output at every t_print time (hrs)
	
	int NumComp;							//Number of species to track
	double end_time;						//Units: hours
	double total_sorption_old;				//Per mass or volume of single fiber
	double total_sorption;					//Per mass or volume of single fiber
	double fiber_length;					//Units: um
	double fiber_diameter;					//Units: um
	
	FILE *OutputFile;						//Output file pointer to the output file for postprocesses
	//Function pointers for the parameters to be evaluated in FINCH for DOGFISH
	double (*eval_R) (int i, int l, const void *data);
	double (*eval_DI) (int i, int l, const void *data);
	double (*eval_kf) (int i, const void *data);
	double (*eval_qs) (int i, const void *data);
	const void *user_data;					//Data structure for users info to calculate Dc and kf
	std::vector<FINCH_DATA> finch_dat;		//Data structure for adsorption kinetics
	std::vector<DOGFISH_PARAM> param_dat;	//Data structure for SKUA specific parameters
	
}DOGFISH_DATA;

void print2file_species_header(FILE *Output, DOGFISH_DATA *dog_dat, int i);

void print2file_DOGFISH_header(DOGFISH_DATA *dog_dat);

void print2file_DOGFISH_result_old(DOGFISH_DATA *dog_dat);

void print2file_DOGFISH_result_new(DOGFISH_DATA *dog_dat);

//----- Default functions for DOGFISH parameters -----
double default_Retardation(int i, int l, const void *data);

double default_IntraDiffusion(int i, int l, const void *data);

double default_FilmMTCoeff(int i, const void *data);

double default_SurfaceConcentration(int i, const void *data);
//----- END Default function definitions -------------

int setup_DOGFISH_DATA(FILE *file, double (*eval_R) (int i, int l, const void *user_data),
					   double (*eval_DI) (int i, int l, const void *user_data),
					   double (*eval_kf) (int i, const void *user_data),
					   double (*eval_qs) (int i, const void *user_data),
					   const void *user_data, DOGFISH_DATA *dog_dat);

int DOGFISH_Executioner(DOGFISH_DATA *dog_dat);

int set_DOGFISH_ICs(DOGFISH_DATA *dog_dat);

int set_DOGFISH_timestep(DOGFISH_DATA *dog_dat);

int DOGFISH_preprocesses(DOGFISH_DATA *dog_dat);

int set_DOGFISH_params(const void *user_data);

int DOGFISH_postprocesses(DOGFISH_DATA *dog_dat);

int DOGFISH_reset(DOGFISH_DATA *dog_dat);

int DOGFISH(DOGFISH_DATA *dog_dat);


//Running DOGFISH tests
int DOGFISH_TESTS();

#endif
