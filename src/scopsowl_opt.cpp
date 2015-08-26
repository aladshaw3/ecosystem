//----------------------------------------
//  Created by Austin Ladshaw on 5/14/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------


#include "scopsowl_opt.h"

//Function to setup the gas phase mole fractions
int SCOPSOWL_OPT_set_y(SCOPSOWL_OPT_DATA *owl_opt)
{
	int success = 0;
	
	double yad_sum = 0.0;
	double rat_frac = 0.0;
	
	//Fill out IC for non-adsorbing species
	int first = 0;
	bool changed = false;
	for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
	{
		if (owl_opt->owl_dat.param_dat[i].Adsorbable == true)
		{
			yad_sum = yad_sum + owl_opt->owl_dat.y[i];
		}
	}
	for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
	{
		if (owl_opt->owl_dat.param_dat[i].Adsorbable == false)
		{
			if (changed == false)
			{
				first = i; changed = true;
				rat_frac = rat_frac + (owl_opt->y_base[i] / owl_opt->y_base[first]);
			}
			else
			{
				rat_frac = rat_frac + (owl_opt->y_base[i] / owl_opt->y_base[first]);
			}
		}
	}
	
	owl_opt->owl_dat.y[first] = (1 - yad_sum) / rat_frac;
	yad_sum = 0.0;
	for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
	{
		if (owl_opt->owl_dat.param_dat[i].Adsorbable == false)
		{
			owl_opt->owl_dat.y[i] = owl_opt->owl_dat.y[first] * (owl_opt->y_base[i] / owl_opt->y_base[first]);
			owl_opt->owl_dat.magpie_dat.gpast_dat[i].y = 0.0;
		}
		else
		{
			owl_opt->owl_dat.magpie_dat.gpast_dat[i].y = owl_opt->owl_dat.y[i];
		}
		yad_sum = yad_sum + owl_opt->owl_dat.y[i];
	}
	if (yad_sum > (1.0 + 1e-6) || yad_sum < (1.0 - 1e-6)) {mError(invalid_gas_sum); return -1;}
	
	
	return success;
}

//Initial guessing algorithm to roughly approximate the diffusion parameter
int initial_guess_SCOPSOWL(SCOPSOWL_OPT_DATA *owl_opt)
{
	int success = 0;
	if (owl_opt->Optimize == false) return 0;
	owl_opt->evaluation = 0;
	
	//Loop over max iter
	int k = 0;
	double base = 1.0;
	double best_par = owl_opt->param_guess;
	double par_old = owl_opt->param_guess;
	double par_new = owl_opt->param_guess;
	double best_norm = 1.0;
	double par_min = DBL_MIN, par_max = DBL_MAX;
	bool Increasing = false;
	do
	{
		//Set the parameters based on optimization call
		owl_opt->owl_dat.param_dat[owl_opt->adsorb_index].ref_diffusion = owl_opt->param_guess;
		std::cout << "\nD_ref Guess (um^2/hr) = " << owl_opt->param_guess;
		
		//Establish the Initial conditions
		owl_opt->owl_dat.total_steps = 0;
		success = set_SCOPSOWL_ICs(&owl_opt->owl_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Set initial points and time steps
		owl_opt->q_sim[0] = owl_opt->owl_dat.param_dat[owl_opt->adsorb_index].qIntegralAvg;
		for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
		{
			owl_opt->owl_dat.finch_dat[i].t_old = owl_opt->t[0];
			owl_opt->owl_dat.finch_dat[i].dt = owl_opt->t[1] - owl_opt->t[0];
			owl_opt->owl_dat.finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].dt + owl_opt->owl_dat.finch_dat[i].t_old;
			
			if (owl_opt->owl_dat.SurfDiff == true && owl_opt->owl_dat.Heterogeneous == true)
			{
				for (int l=0; l<owl_opt->owl_dat.finch_dat[i].LN; l++)
				{
					owl_opt->owl_dat.skua_dat[l].finch_dat[i].dt = owl_opt->owl_dat.finch_dat[i].dt;
					owl_opt->owl_dat.skua_dat[l].finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].t;
					owl_opt->owl_dat.skua_dat[l].t_old = owl_opt->owl_dat.finch_dat[i].t_old;
					owl_opt->owl_dat.skua_dat[l].t = owl_opt->owl_dat.finch_dat[i].t;
				}
			}
		}
		owl_opt->owl_dat.t = owl_opt->t[1];
		owl_opt->owl_dat.t_old = owl_opt->t[0];
		
		//Messaging to user
		std::cout << "\n--------------------\nInitial Guess Operation " << (owl_opt->evaluation+1) << std::endl;
		int ppd = owl_opt->current_points / 10;
		int dot = 0;
		
		//Begin looping to complete simulation
		owl_opt->f_bias = 0.0;
		owl_opt->e_norm = 0.0;
		for (int n=1; n<owl_opt->current_points; n++)
		{
			//Printout messages as the job builds to completion
			if ((dot*ppd) < n && (dot*10) < 100)
			{
				std::cout << "[" << (dot*10) << "%]\n";
				dot++;
			}
			
			//Conditional update statements
			if (owl_opt->owl_dat.finch_dat[owl_opt->adsorb_index].Update == true)
			{
				success = SCOPSOWL_reset(&owl_opt->owl_dat);
				if (success != 0) {mError(simulation_fail); return -1;}
			}
			
			//Set the timestep
			for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
			{
				owl_opt->owl_dat.finch_dat[i].dt = owl_opt->t[n] - owl_opt->t[n-1];
				owl_opt->owl_dat.finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].dt + owl_opt->owl_dat.finch_dat[i].t_old;
				
				if (owl_opt->owl_dat.SurfDiff == true && owl_opt->owl_dat.Heterogeneous == true)
				{
					for (int l=0; l<owl_opt->owl_dat.finch_dat[i].LN; l++)
					{
						owl_opt->owl_dat.skua_dat[l].finch_dat[i].dt = owl_opt->owl_dat.finch_dat[i].dt;
						owl_opt->owl_dat.skua_dat[l].finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].t;
						owl_opt->owl_dat.skua_dat[l].t_old = owl_opt->owl_dat.finch_dat[i].t_old;
						owl_opt->owl_dat.skua_dat[l].t = owl_opt->owl_dat.finch_dat[i].t;
					}
				}
			}
			owl_opt->owl_dat.t = owl_opt->owl_dat.finch_dat[0].t;
			owl_opt->owl_dat.t_old = owl_opt->owl_dat.finch_dat[0].t_old;
			
			//Call skua simulation
			success = SCOPSOWL_Executioner(&owl_opt->owl_dat);
			owl_opt->owl_dat.total_steps++;
			if (success == 0)
			{
				for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
					owl_opt->owl_dat.finch_dat[i].Update = true;
			}
			else {mError(simulation_fail); owl_opt->owl_dat.finch_dat[owl_opt->adsorb_index].Update = false; return -1;}
			
			//Store solution and form residuals
			owl_opt->q_sim[n] = owl_opt->owl_dat.param_dat[owl_opt->adsorb_index].qIntegralAvg;
			owl_opt->f_bias = owl_opt->f_bias + ((owl_opt->q_sim[n] / owl_opt->simulation_equil) - (owl_opt->q_data[n] / owl_opt->current_equil));
			owl_opt->e_norm = owl_opt->e_norm + pow(((owl_opt->q_sim[n] / owl_opt->simulation_equil) - (owl_opt->q_data[n] / owl_opt->current_equil)), 2.0);
		}
		
		owl_opt->e_norm = sqrt(owl_opt->e_norm);
		if (k == 0)
		{
			base = owl_opt->e_norm;
		}
		owl_opt->evaluation++;
		std::cout << "[100%] Complete!\n--------------------" << std::endl;
		std::cout << "Max Bias\tCurnt Bias\tMin Bias\tE.Norm\n";
		std::cout << owl_opt->max_bias <<" >=\t"<< owl_opt->f_bias <<" >=\t"<< owl_opt->min_bias <<"\t"<< owl_opt->e_norm << std::endl;
		owl_opt->total_eval = owl_opt->total_eval + owl_opt->owl_dat.total_steps;
		
		//Check convergence
		if ((owl_opt->e_norm/base) <= owl_opt->rel_tol_norm)
		{
			std::cout << "Euclidean Norm is within tolerence!\n";
			best_norm = owl_opt->e_norm;
			best_par = owl_opt->param_guess;
			break;
		}
		if (fabs(owl_opt->f_bias) <= owl_opt->abs_tol_bias)
		{
			std::cout << "Function Bias is within tolerence!\n";
			best_norm = owl_opt->e_norm;
			best_par = owl_opt->param_guess;
			break;
		}
		
		//Alter the parameter guess based on biases
		if (k == 0)
		{
			if (owl_opt->f_bias > 0.0)
			{
				par_max = par_new;
				par_new = par_new*(1.0 - (owl_opt->f_bias/owl_opt->max_bias));
				Increasing = false;
			}
			else
			{
				par_min = par_new;
				par_new = par_new/(1.0 - (owl_opt->f_bias/owl_opt->min_bias));
				Increasing = true;
			}
			best_par = owl_opt->param_guess;
			best_norm = owl_opt->e_norm;
		}
		else
		{
			//Store the best results
			if (owl_opt->e_norm < best_norm)
			{
				best_norm = owl_opt->e_norm;
				best_par = owl_opt->param_guess;
			}
			
			//Positive bias
			if (owl_opt->f_bias > 0.0 && owl_opt->f_bias_old > 0.0)
			{
				par_max = par_new;
				//Showing improvement
				if (fabs(owl_opt->f_bias) < fabs(owl_opt->f_bias_old))
				{
					if (Increasing == false)
					{
						if ((owl_opt->f_bias/owl_opt->f_bias_old) <= 0.5)
							par_new = par_new*(1.0 - (owl_opt->f_bias/owl_opt->max_bias));
						else
							par_new = par_new*((1.0 - (owl_opt->f_bias/owl_opt->max_bias))/(2.0*(owl_opt->f_bias/owl_opt->f_bias_old)));
						if (par_new <= par_min)
						{
							par_new = (0.5 * par_min) + (0.5 * par_old);
						}
					}
					else
					{
						if ((owl_opt->f_bias/owl_opt->f_bias_old) <= 0.5)
							par_new = par_new/(1.0 - (owl_opt->f_bias/owl_opt->max_bias));
						else
							par_new = par_new/((1.0 - (owl_opt->f_bias/owl_opt->max_bias))/(2.0*(owl_opt->f_bias/owl_opt->f_bias_old)));
						if (par_new >= par_max)
						{
							par_new = (0.5 * par_max) + (0.5 * par_old);
						}
					}
				}
				//Getting worse
				else
				{
					if (Increasing == false)
					{
						par_new = par_old/((1.0 - (owl_opt->f_bias/owl_opt->max_bias)));
						if (par_new >= par_max)
						{
							par_new = (0.5 * par_max) + (0.5 * par_old);
						}
						Increasing = true;
					}
					else
					{
						par_new = par_old*((1.0 - (owl_opt->f_bias/owl_opt->max_bias)));
						if (par_new <= par_min)
						{
							par_new = (0.5 * par_min) + (0.5 * par_old);
						}
						Increasing = false;
					}
				}
			}
			//Negative bias
			else if (owl_opt->f_bias < 0.0 && owl_opt->f_bias_old < 0.0)
			{
				par_min = par_new;
				//Showing improvement
				if (fabs(owl_opt->f_bias) < fabs(owl_opt->f_bias_old))
				{
					if (Increasing == true)
					{
						if ((owl_opt->f_bias/owl_opt->f_bias_old) <= 0.5)
							par_new = par_new/(1.0 - (owl_opt->f_bias/owl_opt->min_bias));
						else
							par_new = par_new/((1.0 - (owl_opt->f_bias/owl_opt->min_bias))/(2.0*(owl_opt->f_bias/owl_opt->f_bias_old)));
						if (par_new >= par_max)
						{
							par_new = (0.5 * par_max) + (0.5 * par_old);
						}
					}
					else
					{
						if ((owl_opt->f_bias/owl_opt->f_bias_old) <= 0.5)
							par_new = par_new*(1.0 - (owl_opt->f_bias/owl_opt->min_bias));
						else
							par_new = par_new*((1.0 - (owl_opt->f_bias/owl_opt->min_bias))/(2.0*(owl_opt->f_bias/owl_opt->f_bias_old)));
						if (par_new <= par_min)
						{
							par_new = (0.5 * par_min) + (0.5 * par_old);
						}
					}
				}
				//Getting worse
				else
				{
					if (Increasing == false)
					{
						par_new = par_old/((1.0 - (owl_opt->f_bias/owl_opt->min_bias)));
						if (par_new >= par_max)
						{
							par_new = (0.5 * par_max) + (0.5 * par_old);
						}
						Increasing = true;
					}
					else
					{
						par_new = par_old*((1.0 - (owl_opt->f_bias/owl_opt->min_bias)));
						if (par_new <= par_min)
						{
							par_new = (0.5 * par_min) + (0.5 * par_old);
						}
						Increasing = false;
					}
				}
			}
			//Bias swap
			else
			{
				if (owl_opt->f_bias > 0.0)
				{
					par_max = par_new;
				}
				else
				{
					par_min	= par_new;
				}
				double total_bias = fabs(owl_opt->f_bias) + fabs(owl_opt->f_bias_old);
				if (par_old > par_new)
					Increasing = true;
				else
					Increasing = false;
				par_new = ((fabs(owl_opt->f_bias)/total_bias)*par_old) + ((fabs(owl_opt->f_bias_old)/total_bias)*par_new);
			}
		}
		
		//Update for next iteration
		par_old = owl_opt->param_guess;
		owl_opt->param_guess = fabs(par_new);
		owl_opt->e_norm_old = owl_opt->e_norm;
		owl_opt->f_bias_old = owl_opt->f_bias;
		
		k++;
	} while(k < owl_opt->max_guess_iter);
	
	owl_opt->param_guess = best_par;
	owl_opt->e_norm = best_norm;
	std::cout << "\nBest Guess for LMFIT Iterations: D_ref = " << best_par << " um^2/hr" << std::endl;
	owl_opt->evaluation = 0;
	return success;
}

//LMIFT function to evaluate the SCOPSOWL simulations
void eval_SCOPSOWL_Uptake(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	int success = 0;
	SCOPSOWL_OPT_DATA *owl_opt = (SCOPSOWL_OPT_DATA *) data;
	
	//Set the parameters based on optimization call
	if (owl_opt->Optimize == true)
	{
		owl_opt->owl_dat.param_dat[owl_opt->adsorb_index].ref_diffusion = fabs(par[0]);
	}
		
	//Establish the Initial conditions
	owl_opt->owl_dat.total_steps = 0;
	success = set_SCOPSOWL_ICs(&owl_opt->owl_dat);
	if (success != 0) {mError(simulation_fail); *info = -1; return;}
		
	//Set initial points and time steps
	owl_opt->q_sim[0] = owl_opt->owl_dat.param_dat[owl_opt->adsorb_index].qIntegralAvg;
	for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
	{
		owl_opt->owl_dat.finch_dat[i].dt = owl_opt->t[1] - owl_opt->t[0];
		owl_opt->owl_dat.finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].dt + owl_opt->owl_dat.finch_dat[i].t_old;
			
		if (owl_opt->owl_dat.SurfDiff == true && owl_opt->owl_dat.Heterogeneous == true)
		{
			for (int l=0; l<owl_opt->owl_dat.finch_dat[i].LN; l++)
			{
				owl_opt->owl_dat.skua_dat[l].finch_dat[i].dt = owl_opt->owl_dat.finch_dat[i].dt;
				owl_opt->owl_dat.skua_dat[l].finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].t;
				owl_opt->owl_dat.skua_dat[l].t_old = owl_opt->owl_dat.finch_dat[i].t_old;
				owl_opt->owl_dat.skua_dat[l].t = owl_opt->owl_dat.finch_dat[i].t;
			}
		}
	}
	owl_opt->owl_dat.t = owl_opt->owl_dat.finch_dat[0].t;
	owl_opt->owl_dat.t_old = owl_opt->owl_dat.finch_dat[0].t_old;
		
	//Messaging to user
	std::cout << "\n--------------------\nLMFIT Operation " << (owl_opt->evaluation+1) << std::endl;
	int ppd = owl_opt->current_points / 10;
	int dot = 0;
		
	//Begin looping to complete simulation
	fvec[0] = 0.0;
	owl_opt->f_bias = 0.0;
	for (int n=1; n<owl_opt->current_points; n++)
	{
		//Printout messages as the job builds to completion
		if ((dot*ppd) < n && (dot*10) < 100)
		{
			std::cout << "[" << (dot*10) << "%]\n";
			dot++;
		}
			
		//Conditional update statements
		if (owl_opt->owl_dat.finch_dat[owl_opt->adsorb_index].Update == true)
		{
			success = SCOPSOWL_reset(&owl_opt->owl_dat);
			if (success != 0) {mError(simulation_fail); *info = -1; return;}
		}
			
		//Set the timestep
		for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
		{
			owl_opt->owl_dat.finch_dat[i].dt = owl_opt->t[n] - owl_opt->t[n-1];
			owl_opt->owl_dat.finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].dt + owl_opt->owl_dat.finch_dat[i].t_old;
				
			if (owl_opt->owl_dat.SurfDiff == true && owl_opt->owl_dat.Heterogeneous == true)
			{
				for (int l=0; l<owl_opt->owl_dat.finch_dat[i].LN; l++)
				{
					owl_opt->owl_dat.skua_dat[l].finch_dat[i].dt = owl_opt->owl_dat.finch_dat[i].dt;
					owl_opt->owl_dat.skua_dat[l].finch_dat[i].t = owl_opt->owl_dat.finch_dat[i].t;
					owl_opt->owl_dat.skua_dat[l].t_old = owl_opt->owl_dat.finch_dat[i].t_old;
					owl_opt->owl_dat.skua_dat[l].t = owl_opt->owl_dat.finch_dat[i].t;
				}
			}
		}
		owl_opt->owl_dat.t = owl_opt->owl_dat.finch_dat[0].t;
		owl_opt->owl_dat.t_old = owl_opt->owl_dat.finch_dat[0].t_old;
			
		//Call skua simulation
		success = SCOPSOWL_Executioner(&owl_opt->owl_dat);
		owl_opt->owl_dat.total_steps++;
		if (success == 0)
		{
			for (int i=0; i<owl_opt->owl_dat.magpie_dat.sys_dat.N; i++)
				owl_opt->owl_dat.finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); owl_opt->owl_dat.finch_dat[owl_opt->adsorb_index].Update = false; *info = -1; return;}
			
		//Store solution and form residuals
		owl_opt->q_sim[n] = owl_opt->owl_dat.param_dat[owl_opt->adsorb_index].qIntegralAvg;
		fvec[n] = (owl_opt->q_sim[n] / owl_opt->simulation_equil) - (owl_opt->q_data[n] / owl_opt->current_equil);
		owl_opt->f_bias = owl_opt->f_bias + ((owl_opt->q_sim[n] / owl_opt->simulation_equil) - (owl_opt->q_data[n] / owl_opt->current_equil));
		
	}
		
	owl_opt->evaluation++;
	std::cout << "[100%] Complete!\n--------------------" << std::endl;
	owl_opt->total_eval = owl_opt->total_eval + owl_opt->owl_dat.total_steps;
	if (owl_opt->Optimize == false) *info = -1;
	if (owl_opt->Rough == true) *info = -1;
		
}

//Optimization routine for SCOPSOWL
int SCOPSOWL_OPTIMIZE(const char *scene, const char *sorbent, const char *comp, const char *sorbate, const char *data)
{
	int success = 0;
	std::string sceneName, sorbentName, compName, sorbateName, dataName;
	
	//Check the input files
	if (scene == NULL || sorbent == NULL || comp == NULL || sorbate == NULL || data == NULL)
	{
		std::cout << "Enter name of Scenario File: ";
		std::cin >> sceneName;
		std::cout << "Enter name of Adsorbent Property File: ";
		std:: cin >> sorbentName;
		std::cout << "Enter name of Component Property File: ";
		std::cin >> compName;
		std::cout << "Enter name of Adsorbate Property File: ";
		std::cin >> sorbateName;
		std::cout << "Enter name of Adsorption Data File: ";
		std::cin >> dataName;
		std::cout << "\n";
		
		scene = sceneName.c_str();
		sorbent = sorbentName.c_str();
		comp = compName.c_str();
		sorbate = sorbateName.c_str();
		data = dataName.c_str();
	}
	
	std::ifstream sceneFile ( scene );
	std::ifstream sorbentFile ( sorbent );
	std::ifstream compFile ( comp );
	std::ifstream sorbateFile ( sorbate );
	std::ifstream dataFile ( data );
	
	if (sceneFile.good()==false || sorbentFile.good()==false || compFile.good()==false || sorbateFile.good()==false || dataFile.good()==false)
	{
		mError(file_dne);
		return -1;
	}
	
	//Declarations
	SCOPSOWL_OPT_DATA owl_opt;
	MIXED_GAS mixture;
	lm_status_struct status;
	lm_control_struct control;
	int i_read;
	double d_read;
	std::string s_read;
	double time;
	FILE *ParamResults, *Comparison;
	
	//Initializations
	time = clock();
	ParamResults = fopen("output/SCOPSOWL_OPT_ParamFile.txt", "w+");
	Comparison = fopen("output/SCOPSOWL_OPT_CompareFile.txt", "w+");
	if (ParamResults == nullptr)
	{
		system("mkdir output");
		ParamResults = fopen("output/SCOPSOWL_OPT_ParamFile.txt", "w+");
		Comparison = fopen("output/SCOPSOWL_OPT_CompareFile.txt", "w+");
	}
	
	owl_opt.owl_dat.total_steps = 0;
	owl_opt.total_eval = 0;
	owl_opt.evaluation = 0;
	owl_opt.owl_dat.Print2Console = false;
	owl_opt.owl_dat.Print2File = false;
	owl_opt.owl_dat.NonLinear = true;
	owl_opt.owl_dat.SurfDiff = true;
	owl_opt.owl_dat.level = 1;
	owl_opt.num_params = 1;
	control = lm_control_double;
	control.printflags = 3;
	status.nfev = 0;
	
	// Read (1) sceneFile
	sceneFile >> i_read;
	if (i_read == 0) owl_opt.Optimize = false;
	else if (i_read == 1) owl_opt.Optimize = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read;
	if (i_read == 0) owl_opt.Rough = false;
	else if (i_read == 1) owl_opt.Rough = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read; owl_opt.diffusion_type = i_read;
	if (i_read == 0) owl_opt.owl_dat.eval_surfDiff = (*default_Dc);//GOOD
	else if (i_read == 1) owl_opt.owl_dat.eval_surfDiff = (*simple_darken_Dc);//GOOD
	else if (i_read == 2) owl_opt.owl_dat.eval_surfDiff = (*theoretical_darken_Dc);//GOOD
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read;
	if (i_read == 0) owl_opt.owl_dat.DirichletBC = false;
	else if (i_read == 1) owl_opt.owl_dat.DirichletBC = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> d_read; owl_opt.owl_dat.magpie_dat.sys_dat.PT = d_read;
	owl_opt.owl_dat.total_pressure = d_read;
	sceneFile >> d_read; owl_opt.owl_dat.gas_velocity = d_read;
	sceneFile >> i_read; owl_opt.owl_dat.magpie_dat.sys_dat.N = i_read;
	owl_opt.owl_dat.magpie_dat.gsta_dat.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	owl_opt.owl_dat.magpie_dat.gpast_dat.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	owl_opt.owl_dat.magpie_dat.mspd_dat.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	owl_opt.owl_dat.param_dat.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	owl_opt.owl_dat.finch_dat.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	owl_opt.owl_dat.y.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	owl_opt.y_base.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	sceneFile >> d_read; owl_opt.owl_dat.magpie_dat.sys_dat.qT = d_read;
	int number_adsorbable = 0;
	for (int i=0; i<owl_opt.owl_dat.magpie_dat.sys_dat.N; i++)
	{
		sceneFile >> s_read; owl_opt.owl_dat.param_dat[i].speciesName = s_read;
		sceneFile >> i_read;
		if (i_read == 0) owl_opt.owl_dat.param_dat[i].Adsorbable = false;
		else if (i_read == 1)
		{
			owl_opt.owl_dat.param_dat[i].Adsorbable = true;
			owl_opt.adsorb_index = i;
			number_adsorbable++;
		}
		else {mError(invalid_boolean); return -1;}
		sceneFile >> d_read; owl_opt.y_base[i] = d_read;
		owl_opt.owl_dat.y[i] = 0.0;
		sceneFile >> d_read; owl_opt.owl_dat.magpie_dat.gpast_dat[i].x = d_read;
		owl_opt.owl_dat.param_dat[i].xIC = d_read;
		
		//Additional magpie initializations
		owl_opt.owl_dat.magpie_dat.mspd_dat[i].eta.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
		owl_opt.owl_dat.magpie_dat.gpast_dat[i].gama_inf.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
		owl_opt.owl_dat.magpie_dat.gpast_dat[i].po.resize(owl_opt.owl_dat.magpie_dat.sys_dat.N);
	}
	if (number_adsorbable > 1 || number_adsorbable == 0) {mError(invalid_components); return -1;}
	sceneFile.close();
	
	//Initialize gas mixture data
	success = initialize_data(owl_opt.owl_dat.magpie_dat.sys_dat.N, &mixture);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	// Read (2) sorbentFile
	sorbentFile >> i_read;
	if (i_read == 0) owl_opt.owl_dat.Heterogeneous = false;
	else if (i_read == 1) owl_opt.owl_dat.Heterogeneous = true;
	else {mError(invalid_boolean); return -1;}
	
	//IMPORTANT!!! If Heterogeneous is false, you must change the function pointer
	if (owl_opt.owl_dat.Heterogeneous == false)
		owl_opt.owl_dat.eval_surfDiff = (*default_surf_diffusion);
	
	sorbentFile >> i_read; owl_opt.owl_dat.coord_macro = i_read;
	if (i_read > 2 || i_read < 0) {mError(invalid_boolean); return -1;}
	if (i_read == 0 || i_read == 1)
	{
		sorbentFile >> d_read; owl_opt.owl_dat.char_macro = d_read;
	}
	else
		owl_opt.owl_dat.char_macro = 1.0;
	sorbentFile	>> d_read; owl_opt.owl_dat.pellet_radius = d_read;
	sorbentFile >> d_read; owl_opt.owl_dat.pellet_density = d_read;
	sorbentFile >> d_read; owl_opt.owl_dat.binder_porosity = d_read;
	sorbentFile >> d_read; owl_opt.owl_dat.binder_poresize = d_read;
	if (owl_opt.owl_dat.Heterogeneous == true)
	{
		sorbentFile >> i_read; owl_opt.owl_dat.coord_micro = i_read;
		if (i_read > 2 || i_read < 0) {mError(invalid_boolean); return -1;}
		if (i_read == 0 || i_read == 1)
		{
			sorbentFile >> d_read; owl_opt.owl_dat.char_micro = d_read;
		}
		else
			owl_opt.owl_dat.char_micro = 1.0;
	}
	sorbentFile >> d_read; owl_opt.owl_dat.crystal_radius = d_read;
	sorbentFile >> d_read; owl_opt.owl_dat.binder_fraction = d_read;
	sorbentFile.close();
	
	// Read (3) compFile
	for (int i=0; i<owl_opt.owl_dat.magpie_dat.sys_dat.N; i++)
	{
		compFile >> d_read; mixture.species_dat[i].molecular_weight = d_read;
		compFile >> d_read; mixture.species_dat[i].specific_heat = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Viscosity = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Temp = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Const = d_read;
	}
	compFile.close();
	
	// Read (4) sorbateFile
	for (int i=0; i<owl_opt.owl_dat.magpie_dat.sys_dat.N; i++)
	{
		//Only read in data for adsorbable components in order of appearence
		if (owl_opt.owl_dat.param_dat[i].Adsorbable == true)
		{
			sorbateFile >> d_read; owl_opt.owl_dat.param_dat[i].ref_diffusion = d_read;			//um^2/hr
			sorbateFile >> d_read; owl_opt.owl_dat.param_dat[i].activation_energy = d_read;		//J/mol
			sorbateFile >> d_read; owl_opt.owl_dat.param_dat[i].ref_temperature = d_read;		//K
			sorbateFile >> d_read; owl_opt.owl_dat.param_dat[i].affinity = d_read;				//-
			if (owl_opt.Optimize == true)
			{
				owl_opt.owl_dat.param_dat[i].ref_diffusion = 1.0;
				owl_opt.owl_dat.param_dat[i].activation_energy = 0.0;
				owl_opt.owl_dat.param_dat[i].ref_temperature = 0.0;
				owl_opt.owl_dat.param_dat[i].affinity = 0.0;
			}
			sorbateFile >> d_read; owl_opt.owl_dat.magpie_dat.mspd_dat[i].v = d_read;				//cm^3/mol
			sorbateFile >> d_read; owl_opt.owl_dat.magpie_dat.gsta_dat[i].qmax = d_read;			//mol/kg
			sorbateFile >> i_read; owl_opt.owl_dat.magpie_dat.gsta_dat[i].m = i_read;				//-
			owl_opt.owl_dat.magpie_dat.gsta_dat[i].dHo.resize(owl_opt.owl_dat.magpie_dat.gsta_dat[i].m);
			owl_opt.owl_dat.magpie_dat.gsta_dat[i].dSo.resize(owl_opt.owl_dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<owl_opt.owl_dat.magpie_dat.gsta_dat[i].m; n++)
			{
				sorbateFile >> d_read; owl_opt.owl_dat.magpie_dat.gsta_dat[i].dHo[n] = d_read;	//J/mol
				sorbateFile >> d_read; owl_opt.owl_dat.magpie_dat.gsta_dat[i].dSo[n] = d_read;	//J/K/mol
			}
		}
		//Otherwise, set values to zeros
		else
		{
			owl_opt.owl_dat.param_dat[i].ref_diffusion = 0.0;				//um^2/hr
			owl_opt.owl_dat.param_dat[i].activation_energy = 0.0;			//J/mol
			owl_opt.owl_dat.param_dat[i].ref_temperature = 0.0;				//K
			owl_opt.owl_dat.param_dat[i].affinity = 0.0;					//-
			owl_opt.owl_dat.magpie_dat.mspd_dat[i].v = 0.0;				//cm^3/mol
			owl_opt.owl_dat.magpie_dat.gsta_dat[i].qmax = 0.0;			//mol/kg
			owl_opt.owl_dat.magpie_dat.gsta_dat[i].m = 1;					//-
			owl_opt.owl_dat.magpie_dat.gsta_dat[i].dHo.resize(owl_opt.owl_dat.magpie_dat.gsta_dat[i].m);
			owl_opt.owl_dat.magpie_dat.gsta_dat[i].dSo.resize(owl_opt.owl_dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<owl_opt.owl_dat.magpie_dat.gsta_dat[i].m; n++)
			{
				owl_opt.owl_dat.magpie_dat.gsta_dat[i].dHo[n] = 0.0;	//J/mol
				owl_opt.owl_dat.magpie_dat.gsta_dat[i].dSo[n] = 0.0;	//J/K/mol
			}
		}
	}
	sorbateFile.close();
	
	//Call the setup function
	success = setup_SCOPSOWL_DATA(NULL, default_adsorption, default_retardation, default_pore_diffusion, default_filmMassTransfer, owl_opt.owl_dat.eval_surfDiff, (void *)&owl_opt.owl_dat, &mixture, &owl_opt.owl_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	// Read (5) top of data file
	dataFile >> i_read; owl_opt.num_curves = i_read;
	
	//Make a header for the parameter file
	fprintf(ParamResults, "Curve#\tT[K]\tp[kPa]\ttheta\tD_ref[um^2/hr]\tE.Norm\n");

	//Start Optimization Loop
	int curve = 0;
	double par[owl_opt.num_params];
	owl_opt.param_guess_old = 0.001;
	do
	{
		// Continue reading data file
		dataFile >> i_read; owl_opt.current_points = i_read;
		owl_opt.q_data.resize(owl_opt.current_points);
		owl_opt.q_sim.resize(owl_opt.current_points);
		owl_opt.t.resize(owl_opt.current_points);
		dataFile >> d_read; owl_opt.owl_dat.magpie_dat.sys_dat.T = d_read;
		owl_opt.owl_dat.gas_temperature = d_read;
		owl_opt.current_temp = d_read;
		dataFile >> d_read; owl_opt.owl_dat.y[owl_opt.adsorb_index] = d_read / owl_opt.owl_dat.magpie_dat.sys_dat.PT;
		owl_opt.current_press = d_read;
		dataFile >> d_read; owl_opt.current_equil = d_read;
		
		//Set up the gas molefractions
		success = SCOPSOWL_OPT_set_y(&owl_opt);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		std::cout << "\nPerforming operations for curve " << (curve+1) << " . Stand by:\n";
		std::cout << "Gas Composition\n-----------------\n";
		for (int i=0; i<owl_opt.owl_dat.magpie_dat.sys_dat.N; i++)
		{
			std::cout << "y[" << i << "] = " << owl_opt.owl_dat.y[i] << std::endl;
		}
		
		//Establish the final equilibrium point
		owl_opt.owl_dat.magpie_dat.sys_dat.Recover = false;
		if (owl_opt.owl_dat.magpie_dat.sys_dat.N > 1)
			owl_opt.owl_dat.magpie_dat.sys_dat.Carrier = true;
		else
			owl_opt.owl_dat.magpie_dat.sys_dat.Carrier = false;
		owl_opt.owl_dat.magpie_dat.sys_dat.Output = false;
		success = MAGPIE(&owl_opt.owl_dat.magpie_dat);
		if (success < 0 || success > 5) {mError(simulation_fail); return -1;}
		else success = 0;
		owl_opt.total_eval = owl_opt.total_eval + owl_opt.owl_dat.magpie_dat.sys_dat.total_eval;
		owl_opt.simulation_equil = owl_opt.owl_dat.magpie_dat.sys_dat.qT;
		
		std::cout << "Ads. Equil. = " << owl_opt.simulation_equil << std::endl;
		
		//Read in the curve data points
		double maxsum = 0.0, minsum = 0.0;
		for (int n=0; n<owl_opt.current_points; n++)
		{
			dataFile >> d_read; owl_opt.t[n] = d_read;
			dataFile >> d_read; owl_opt.q_data[n] = d_read;
			if (n > 0)
			{
				maxsum = maxsum + (1.0 - (owl_opt.q_data[n]/owl_opt.current_equil));
				minsum = minsum + (0.0 - (owl_opt.q_data[n]/owl_opt.current_equil));
			}
		}
		owl_opt.max_bias = maxsum;
		owl_opt.min_bias = minsum;
		
		//Perform initial guess routine
		owl_opt.param_guess = 0.001;
		success = initial_guess_SCOPSOWL(&owl_opt);
		if (success != 0) {mError(simulation_fail); return -1;}
		par[0] = owl_opt.param_guess;
		status.fnorm = owl_opt.e_norm;
		
		//Call the LMFIT routine
		owl_opt.evaluation = 0;
		lmmin(owl_opt.num_params, par, owl_opt.current_points, (void *)&owl_opt, eval_SCOPSOWL_Uptake, &control, &status, lm_printout_std);
		owl_opt.total_eval = owl_opt.total_eval + status.nfev;
		owl_opt.param_guess_old = fabs(par[0]);
		printf("\nStatus after %d function evaluations:\n  %s\n\n",status.nfev,lm_infmsg[status.info]);
		std::cout << "F.Bias:\t" << owl_opt.f_bias << std::endl;
		std::cout << "E.Norm:\t" << status.fnorm << std::endl << std::endl;
		
		//Header file for curve results
		if (owl_opt.Optimize == true)
			fprintf(ParamResults, "%i\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n",(curve+1),owl_opt.current_temp,owl_opt.current_press,owl_opt.simulation_equil/owl_opt.owl_dat.magpie_dat.gsta_dat[owl_opt.adsorb_index].qmax,fabs(par[0]),status.fnorm);
		else
			fprintf(ParamResults, "%i\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n",(curve+1),owl_opt.current_temp,owl_opt.current_press,owl_opt.simulation_equil/owl_opt.owl_dat.magpie_dat.gsta_dat[owl_opt.adsorb_index].qmax,owl_opt.owl_dat.param_dat[owl_opt.adsorb_index].ref_diffusion,status.fnorm);
		fprintf(Comparison,"Curve\t%i\np[kPa]\t%.6g\tT[K]\t%.6g\n",(curve+1),owl_opt.current_press,owl_opt.current_temp);
		fprintf(Comparison,"Time[hr]\tData(Raw)\tData(Normal)\tModel(Raw)\tModel(Normal)\n");
		
		//Print out current results
		for (int n=0; n<owl_opt.current_points; n++)
		{
			fprintf(Comparison,"%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t",owl_opt.t[n],owl_opt.q_data[n],(owl_opt.q_data[n] / owl_opt.current_equil),owl_opt.q_sim[n],(owl_opt.q_sim[n]/owl_opt.simulation_equil));
			fprintf(Comparison,"\n");
		}
		fprintf(Comparison,"------------------------------------------------------\n\n");
		
		owl_opt.evaluation = 0;
		curve++;
	} while (curve < owl_opt.num_curves);
	//End Loop
	dataFile.close();
	
	//END of optimization
	fclose(ParamResults);
	fclose(Comparison);
	time = clock() - time;
	std::cout << "Simulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Evaluations: " << owl_opt.total_eval << "\n";
	std::cout << "Evaluations/sec: " << owl_opt.total_eval/(time / CLOCKS_PER_SEC) << "\n";
	
	return success;
}