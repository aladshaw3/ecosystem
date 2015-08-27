//----------------------------------------
//  Created by Austin Ladshaw on 5/11/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "skua_opt.h"

//Function to set the gas phase mole fractions correctly
int SKUA_OPT_set_y(SKUA_OPT_DATA *skua_opt)
{
	int success = 0;
	
	double yad_sum = 0.0;
	double rat_frac = 0.0;
	
	//Fill out IC for non-adsorbing species
	int first = 0;
	bool changed = false;
	for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
	{
		if (skua_opt->skua_dat.param_dat[i].Adsorbable == true)
		{
			yad_sum = yad_sum + skua_opt->skua_dat.y[i];
		}
	}
	for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
	{
		if (skua_opt->skua_dat.param_dat[i].Adsorbable == false)
		{
			if (changed == false)
			{
				first = i; changed = true;
				rat_frac = rat_frac + (skua_opt->y_base[i] / skua_opt->y_base[first]);
			}
			else
			{
				rat_frac = rat_frac + (skua_opt->y_base[i] / skua_opt->y_base[first]);
			}
		}
	}
	
	skua_opt->skua_dat.y[first] = (1 - yad_sum) / rat_frac;
	yad_sum = 0.0;
	for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
	{
		if (skua_opt->skua_dat.param_dat[i].Adsorbable == false)
		{
			skua_opt->skua_dat.y[i] = skua_opt->skua_dat.y[first] * (skua_opt->y_base[i] / skua_opt->y_base[first]);
			skua_opt->skua_dat.magpie_dat.gpast_dat[i].y = 0.0;
		}
		else
		{
			skua_opt->skua_dat.magpie_dat.gpast_dat[i].y = skua_opt->skua_dat.y[i];
		}
		yad_sum = yad_sum + skua_opt->skua_dat.y[i];
	}
	if (yad_sum > (1.0 + 1e-6) || yad_sum < (1.0 - 1e-6)) {mError(invalid_gas_sum); return -1;}

	
	return success;
}

//Form a good initial guess for the optimization routine
int initial_guess_SKUA(SKUA_OPT_DATA *skua_opt)
{
	int success = 0;
	if (skua_opt->Optimize == false) return 0;
	skua_opt->evaluation = 0;
	
	//Loop over max iter
	int k = 0;
	double base = 1.0;
	double best_par = skua_opt->param_guess;
	double par_old = skua_opt->param_guess;
	double par_new = skua_opt->param_guess;
	double best_norm = 1.0;
	do
	{
		//Establish the Initial conditions
		skua_opt->skua_dat.total_steps = 0;
		success = set_SKUA_ICs(&skua_opt->skua_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	
		//Set initial points and time steps
		skua_opt->q_sim[0] = skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].uAvg;
		for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
		{
			skua_opt->skua_dat.finch_dat[i].dt = skua_opt->t[1] - skua_opt->t[0];
			skua_opt->skua_dat.finch_dat[i].t = skua_opt->skua_dat.finch_dat[i].dt + skua_opt->skua_dat.finch_dat[i].t_old;
		}
		skua_opt->skua_dat.t = skua_opt->skua_dat.finch_dat[0].t;
		skua_opt->skua_dat.t_old = skua_opt->skua_dat.finch_dat[0].t_old;
	
		//Messaging to user
		std::cout << "\n--------------------\nInitial Guess Operation " << (skua_opt->evaluation+1) << std::endl;
		int ppd = skua_opt->current_points / 10;
		int dot = 0;
	
		//Begin looping to complete simulation
		skua_opt->f_bias = 0.0;
		skua_opt->e_norm = 0.0;
		for (int n=1; n<skua_opt->current_points; n++)
		{
			//Printout messages as the job builds to completion
			if ((dot*ppd) < n && (dot*10) < 100)
			{
				std::cout << "[" << (dot*10) << "%]\n";
				dot++;
			}
		
			//Conditional update statements
			if (skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].Update == true)
			{
				success = SKUA_reset(&skua_opt->skua_dat);
				if (success != 0) {mError(simulation_fail); return -1;}
			}
		
			//Set the timestep
			for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
			{
				skua_opt->skua_dat.finch_dat[i].dt = skua_opt->t[n] - skua_opt->t[n-1];
				skua_opt->skua_dat.finch_dat[i].t = skua_opt->skua_dat.finch_dat[i].dt + skua_opt->skua_dat.finch_dat[i].t_old;
			}
			skua_opt->skua_dat.t = skua_opt->skua_dat.finch_dat[0].t;
			skua_opt->skua_dat.t_old = skua_opt->skua_dat.finch_dat[0].t_old;
		
			//Set the parameters based on optimization call
			if (skua_opt->Optimize == true)
			{
				skua_opt->skua_dat.param_dat[skua_opt->adsorb_index].ref_diffusion = skua_opt->param_guess;
			}
		
			//Call skua simulation
			success = SKUA_Executioner(&skua_opt->skua_dat);
			skua_opt->skua_dat.total_steps++;
			if (success == 0)
			{
				for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
					skua_opt->skua_dat.finch_dat[i].Update = true;
			}
			else {mError(simulation_fail); skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].Update = false; return -1;}
		
			//Store solution and form residuals
			skua_opt->q_sim[n] = skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].uAvg;
			skua_opt->f_bias = skua_opt->f_bias + ((skua_opt->q_sim[n] / skua_opt->simulation_equil) - (skua_opt->q_data[n] / skua_opt->current_equil));
			skua_opt->e_norm = skua_opt->e_norm + pow(((skua_opt->q_sim[n] / skua_opt->simulation_equil) - (skua_opt->q_data[n] / skua_opt->current_equil)), 2.0);
		}
		skua_opt->e_norm = sqrt(skua_opt->e_norm);
		if (k == 0)
		{
			base = skua_opt->e_norm;
		}
		skua_opt->evaluation++;
		std::cout << "[100%] Complete!\n--------------------" << std::endl;
		std::cout << "Max Bias\tCurnt Bias\tMin Bias\tE.Norm\n";
		std::cout << skua_opt->max_bias <<" >=\t"<< skua_opt->f_bias <<" >=\t"<< skua_opt->min_bias <<"\t"<< skua_opt->e_norm << std::endl;
		skua_opt->total_eval = skua_opt->total_eval + skua_opt->skua_dat.total_steps;
		
		//Check convergence
		if ((skua_opt->e_norm/base) <= skua_opt->rel_tol_norm)
		{
			std::cout << "Euclidean Norm is within tolerence!\n";
			best_norm = skua_opt->e_norm;
			best_par = skua_opt->param_guess;
			break;
		}
		if (fabs(skua_opt->f_bias) <= skua_opt->abs_tol_bias)
		{
			std::cout << "Function Bias is within tolerence!\n";
			best_norm = skua_opt->e_norm;
			best_par = skua_opt->param_guess;
			break;
		}
		
		//Alter the parameter guess based on biases
		if (k == 0)
		{
			if (skua_opt->f_bias > 0.0)
			{
				par_new = par_new*(1.0 - (skua_opt->f_bias/skua_opt->max_bias));
			}
			else
			{
				par_new = par_new/(1.0 - (skua_opt->f_bias/skua_opt->min_bias));
			}
			best_par = skua_opt->param_guess;
			best_norm = skua_opt->e_norm;
		}
		else
		{
			//Store the best results
			if (skua_opt->e_norm < best_norm)
			{
				best_norm = skua_opt->e_norm;
				best_par = skua_opt->param_guess;
			}
			
			//Showing improvement
			if (skua_opt->e_norm < (skua_opt->e_norm_old+skua_opt->rel_tol_norm))
			{
				//Positive bias
				if (skua_opt->f_bias > 0.0 && skua_opt->f_bias_old > 0.0)
				{
					//Showing improvement
					if (fabs(skua_opt->f_bias) < fabs(skua_opt->f_bias_old))
					{
						if ((skua_opt->f_bias/skua_opt->f_bias_old) <= 0.5)
							par_new = par_new*(1.0 - (skua_opt->f_bias/skua_opt->max_bias));
						else
							par_new = par_new*((1.0 - (skua_opt->f_bias/skua_opt->max_bias))/(2.0*(skua_opt->f_bias/skua_opt->f_bias_old)));
					}
					//Getting worse
					else
					{
						par_new = par_new*((1.0 - (skua_opt->f_bias/skua_opt->max_bias))/1.5);
					}
				}
				//Negative bias
				else if (skua_opt->f_bias < 0.0 && skua_opt->f_bias_old < 0.0)
				{
					//Showing improvement
					if (fabs(skua_opt->f_bias) < fabs(skua_opt->f_bias_old))
					{
						if ((skua_opt->f_bias/skua_opt->f_bias_old) <= 0.5)
							par_new = par_new/(1.0 - (skua_opt->f_bias/skua_opt->min_bias));
						else
							par_new = par_new/((1.0 - (skua_opt->f_bias/skua_opt->min_bias))/(2.0*(skua_opt->f_bias/skua_opt->f_bias_old)));
					}
					//Getting worse
					else
					{
						par_new = par_new/((1.0 - (skua_opt->f_bias/skua_opt->min_bias))/1.5);
					}
				}
				//Bias swap
				else
				{
					double total_bias = fabs(skua_opt->f_bias) + fabs(skua_opt->f_bias_old);
					par_new = ((fabs(skua_opt->f_bias)/total_bias)*par_old) + ((fabs(skua_opt->f_bias_old)/total_bias)*par_new);
				}
			}
			//Bias swap occured
			else if ( (skua_opt->f_bias < 0.0 && skua_opt->f_bias_old > 0.0) || (skua_opt->f_bias > 0.0 && skua_opt->f_bias_old < 0.0) )
			{
				double total_bias = fabs(skua_opt->f_bias) + fabs(skua_opt->f_bias_old);
				par_new = ((fabs(skua_opt->f_bias)/total_bias)*par_old) + ((fabs(skua_opt->f_bias_old)/total_bias)*par_new);
			}
			//Getting worse
			else
			{
				skua_opt->e_norm = skua_opt->e_norm_old;
				skua_opt->f_bias = skua_opt->f_bias_old;
				std::cout << "Euclidean Norm stopped improving!" << std::endl;
				break;
			}
		}
	
		//Update for next iteration
		par_old = skua_opt->param_guess;
		skua_opt->param_guess = par_new;
		skua_opt->e_norm_old = skua_opt->e_norm;
		skua_opt->f_bias_old = skua_opt->f_bias;
	
		k++;
	} while(k < skua_opt->max_guess_iter);
		
	skua_opt->param_guess = best_par;
	skua_opt->e_norm = best_norm;
	std::cout << "\nBest Guess for LMFIT Iterations: D_ref = " << best_par << " um^2/hr" << std::endl;
	skua_opt->evaluation = 0;
	return success;
}

//Evaluation function used by LMFIT to perform optimization
void eval_SKUA_Uptake(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	int success = 0;
	SKUA_OPT_DATA *skua_opt = (SKUA_OPT_DATA *) data;
	
	//Establish the Initial conditions
	skua_opt->skua_dat.total_steps = 0;
	success = set_SKUA_ICs(&skua_opt->skua_dat);
	if (success != 0) {mError(simulation_fail); *info = -1; return;}
	
	//Set initial points and time steps
	fvec[0] = 0.0;
	skua_opt->q_sim[0] = skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].uAvg;
	for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
	{
		skua_opt->skua_dat.finch_dat[i].dt = skua_opt->t[1] - skua_opt->t[0];
		skua_opt->skua_dat.finch_dat[i].t = skua_opt->skua_dat.finch_dat[i].dt + skua_opt->skua_dat.finch_dat[i].t_old;
	}
	skua_opt->skua_dat.t = skua_opt->skua_dat.finch_dat[0].t;
	skua_opt->skua_dat.t_old = skua_opt->skua_dat.finch_dat[0].t_old;
	
	//Messaging to user
	std::cout << "\n--------------------\nLMFIT Operation " << (skua_opt->evaluation+1) << std::endl;
	int ppd = m_dat / 10;
	int dot = 0;
	
	//Begin looping to complete simulation
	for (int n=1; n<m_dat; n++)
	{
		//Printout messages as the job builds to completion
		if ((dot*ppd) < n && (dot*10) < 100)
		{
			std::cout << "[" << (dot*10) << "%]\n";
			dot++;
		}
		
		//Conditional update statements
		if (skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].Update == true)
		{
			success = SKUA_reset(&skua_opt->skua_dat);
			if (success != 0) {mError(simulation_fail); *info = -1; return;}
		}
		
		//Set the timestep
		for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
		{
			skua_opt->skua_dat.finch_dat[i].dt = skua_opt->t[n] - skua_opt->t[n-1];
			skua_opt->skua_dat.finch_dat[i].t = skua_opt->skua_dat.finch_dat[i].dt + skua_opt->skua_dat.finch_dat[i].t_old;
		}
		skua_opt->skua_dat.t = skua_opt->skua_dat.finch_dat[0].t;
		skua_opt->skua_dat.t_old = skua_opt->skua_dat.finch_dat[0].t_old;
		
		//Set the parameters based on optimization call
		if (skua_opt->Optimize == true)
			skua_opt->skua_dat.param_dat[skua_opt->adsorb_index].ref_diffusion = fabs(par[0]);
		
		//Call skua simulation
		success = SKUA_Executioner(&skua_opt->skua_dat);
		skua_opt->skua_dat.total_steps++;
		if (success == 0)
		{
			for (int i=0; i<skua_opt->skua_dat.magpie_dat.sys_dat.N; i++)
				skua_opt->skua_dat.finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].Update = false; *info = -1; return;}
		
		//Store solution and form residuals
		skua_opt->q_sim[n] = skua_opt->skua_dat.finch_dat[skua_opt->adsorb_index].uAvg;
		fvec[n] = (skua_opt->q_sim[n] / skua_opt->simulation_equil) - (skua_opt->q_data[n] / skua_opt->current_equil);
	}
	skua_opt->evaluation++;
	std::cout << "[100%] Complete!\n--------------------" << std::endl;
	skua_opt->total_eval = skua_opt->total_eval + skua_opt->skua_dat.total_steps;
	if (skua_opt->Optimize == false) *info = -1;
	if (skua_opt->Rough == true) *info = -1;
}

//Main optimization function
int SKUA_OPTIMIZE(const char *scene, const char *sorbent, const char *comp, const char *sorbate, const char *data)
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
	SKUA_OPT_DATA skua_opt;
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
	ParamResults = fopen("output/SKUA_OPT_ParamFile.txt", "w+");
	Comparison = fopen("output/SKUA_OPT_CompareFile.txt", "w+");
	if (ParamResults == nullptr)
	{
		system("mkdir output");
		ParamResults = fopen("output/SKUA_OPT_ParamFile.txt", "w+");
		Comparison = fopen("output/SKUA_OPT_CompareFile.txt", "w+");
	}
	skua_opt.skua_dat.total_steps = 0;
	skua_opt.total_eval = 0;
	skua_opt.evaluation = 0;
	skua_opt.skua_dat.Print2Console = false;
	skua_opt.skua_dat.Print2File = false;
	skua_opt.skua_dat.NonLinear = true;
	skua_opt.num_params = 1;
	control = lm_control_double;
	control.printflags = 3;
	status.nfev = 0;
	
	// Read (1) sceneFile
	sceneFile >> i_read;
	if (i_read == 0) skua_opt.Optimize = false;
	else if (i_read == 1) skua_opt.Optimize = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read;
	if (i_read == 0) skua_opt.Rough = false;
	else if (i_read == 1) skua_opt.Rough = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read; skua_opt.diffusion_type = i_read;
	if (i_read == 0) skua_opt.skua_dat.eval_diff = (*default_Dc);
	else if (i_read == 1) skua_opt.skua_dat.eval_diff = (*simple_darken_Dc);
	else if (i_read == 2) skua_opt.skua_dat.eval_diff = (*theoretical_darken_Dc);
	else {mError(invalid_boolean); return -1;}
	sceneFile >> i_read;
	if (i_read == 0) skua_opt.skua_dat.DirichletBC = false;
	else if (i_read == 1) skua_opt.skua_dat.DirichletBC = true;
	else {mError(invalid_boolean); return -1;}
	sceneFile >> d_read; skua_opt.skua_dat.magpie_dat.sys_dat.PT = d_read;
	sceneFile >> d_read; skua_opt.skua_dat.gas_velocity = d_read;
	sceneFile >> i_read; skua_opt.skua_dat.magpie_dat.sys_dat.N = i_read;
	skua_opt.skua_dat.magpie_dat.gsta_dat.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	skua_opt.skua_dat.magpie_dat.gpast_dat.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	skua_opt.skua_dat.magpie_dat.mspd_dat.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	skua_opt.skua_dat.param_dat.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	skua_opt.skua_dat.finch_dat.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	skua_opt.skua_dat.y.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	skua_opt.y_base.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	sceneFile >> d_read; skua_opt.skua_dat.magpie_dat.sys_dat.qT = d_read;
	int number_adsorbable = 0;
	for (int i=0; i<skua_opt.skua_dat.magpie_dat.sys_dat.N; i++)
	{
		sceneFile >> s_read; skua_opt.skua_dat.param_dat[i].speciesName = s_read;
		sceneFile >> i_read;
		if (i_read == 0) skua_opt.skua_dat.param_dat[i].Adsorbable = false;
		else if (i_read == 1)
		{
			skua_opt.skua_dat.param_dat[i].Adsorbable = true;
			skua_opt.adsorb_index = i;
			number_adsorbable++;
		}
		else {mError(invalid_boolean); return -1;}
		sceneFile >> d_read; skua_opt.y_base[i] = d_read;
		skua_opt.skua_dat.y[i] = 0.0;
		sceneFile >> d_read; skua_opt.skua_dat.param_dat[i].xIC = d_read;
		
		//Additional magpie initializations
		skua_opt.skua_dat.magpie_dat.mspd_dat[i].eta.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
		skua_opt.skua_dat.magpie_dat.gpast_dat[i].gama_inf.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
		skua_opt.skua_dat.magpie_dat.gpast_dat[i].po.resize(skua_opt.skua_dat.magpie_dat.sys_dat.N);
	}
	if (number_adsorbable > 1 || number_adsorbable == 0) {mError(invalid_components); return -1;}
	sceneFile.close();
	
	//Initialize gas mixture data
	success = initialize_data(skua_opt.skua_dat.magpie_dat.sys_dat.N, &mixture);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	// Read (2) sorbentFile
	sorbentFile >> i_read; skua_opt.skua_dat.coord = i_read;
	if (i_read > 2 || i_read < 0) {mError(invalid_boolean); return -1;}
	if (i_read == 0 || i_read == 1)
	{
		sorbentFile >> d_read; skua_opt.skua_dat.char_measure = d_read;
	}
	else
		skua_opt.skua_dat.char_measure = 1.0;
	sorbentFile	>> d_read; skua_opt.skua_dat.pellet_radius = d_read;
	sorbentFile.close();
	
	// Read (3) compFile
	for (int i=0; i<skua_opt.skua_dat.magpie_dat.sys_dat.N; i++)
	{
		compFile >> d_read; mixture.species_dat[i].molecular_weight = d_read;
		compFile >> d_read; mixture.species_dat[i].specific_heat = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Viscosity = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Temp = d_read;
		compFile >> d_read; mixture.species_dat[i].Sutherland_Const = d_read;
	}
	compFile.close();
	
	// Read (4) sorbateFile
	for (int i=0; i<skua_opt.skua_dat.magpie_dat.sys_dat.N; i++)
	{
		//Only read in data for adsorbable components in order of appearence
		if (skua_opt.skua_dat.param_dat[i].Adsorbable == true)
		{
			sorbateFile >> d_read; skua_opt.skua_dat.param_dat[i].ref_diffusion = d_read;			//um^2/hr
			sorbateFile >> d_read; skua_opt.skua_dat.param_dat[i].activation_energy = d_read;		//J/mol
			sorbateFile >> d_read; skua_opt.skua_dat.param_dat[i].ref_temperature = d_read;			//K
			sorbateFile >> d_read; skua_opt.skua_dat.param_dat[i].affinity = d_read;				//-
			if (skua_opt.Optimize == true)
			{
				skua_opt.skua_dat.param_dat[i].activation_energy = 0.0;
				skua_opt.skua_dat.param_dat[i].ref_diffusion = 1.0;
				skua_opt.skua_dat.param_dat[i].ref_temperature = 0.0;
				skua_opt.skua_dat.param_dat[i].affinity = 0.0;
			}
			sorbateFile >> d_read; skua_opt.skua_dat.magpie_dat.mspd_dat[i].v = d_read;				//cm^3/mol
			sorbateFile >> d_read; skua_opt.skua_dat.magpie_dat.gsta_dat[i].qmax = d_read;			//mol/kg
			sorbateFile >> i_read; skua_opt.skua_dat.magpie_dat.gsta_dat[i].m = i_read;				//-
			skua_opt.skua_dat.magpie_dat.gsta_dat[i].dHo.resize(skua_opt.skua_dat.magpie_dat.gsta_dat[i].m);
			skua_opt.skua_dat.magpie_dat.gsta_dat[i].dSo.resize(skua_opt.skua_dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<skua_opt.skua_dat.magpie_dat.gsta_dat[i].m; n++)
			{
				sorbateFile >> d_read; skua_opt.skua_dat.magpie_dat.gsta_dat[i].dHo[n] = d_read;	//J/mol
				sorbateFile >> d_read; skua_opt.skua_dat.magpie_dat.gsta_dat[i].dSo[n] = d_read;	//J/K/mol
			}
		}
		//Otherwise, set values to zeros
		else
		{
			skua_opt.skua_dat.param_dat[i].ref_diffusion = 0.0;				//um^2/hr
			skua_opt.skua_dat.param_dat[i].activation_energy = 0.0;			//J/mol
			skua_opt.skua_dat.magpie_dat.mspd_dat[i].v = 0.0;				//cm^3/mol
			skua_opt.skua_dat.magpie_dat.gsta_dat[i].qmax = 0.0;			//mol/kg
			skua_opt.skua_dat.magpie_dat.gsta_dat[i].m = 1;					//-
			skua_opt.skua_dat.magpie_dat.gsta_dat[i].dHo.resize(skua_opt.skua_dat.magpie_dat.gsta_dat[i].m);
			skua_opt.skua_dat.magpie_dat.gsta_dat[i].dSo.resize(skua_opt.skua_dat.magpie_dat.gsta_dat[i].m);
			for (int n=0; n<skua_opt.skua_dat.magpie_dat.gsta_dat[i].m; n++)
			{
				skua_opt.skua_dat.magpie_dat.gsta_dat[i].dHo[n] = 0.0;	//J/mol
				skua_opt.skua_dat.magpie_dat.gsta_dat[i].dSo[n] = 0.0;	//J/K/mol
			}
		}
	}
	sorbateFile.close();
	
	//Call the setup function
	success = setup_SKUA_DATA(NULL, skua_opt.skua_dat.eval_diff, empirical_kf, (void *)&skua_opt.skua_dat, &mixture, &skua_opt.skua_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	// Read (5) top of data file
	dataFile >> i_read; skua_opt.num_curves = i_read;
	
	//Make a header for the parameter file
	if (skua_opt.Optimize == true)
		fprintf(ParamResults, "Curve#\tT[K]\tp[kPa]\tD_ref[um^2/hr]\tE.Norm\n");
	else
		fprintf(ParamResults, "Curve#\tT[K]\tp[kPa]\tD_ref[um^2/hr]\tE[J/mol]\tE.Norm\n");
	
	// Begin optimization loop
	int curve = 0;
	double par[skua_opt.num_params];
	skua_opt.param_guess_old = 1.0;
	do
	{
		// Continue reading data file
		dataFile >> i_read; skua_opt.current_points = i_read;
		skua_opt.q_data.resize(skua_opt.current_points);
		skua_opt.q_sim.resize(skua_opt.current_points);
		skua_opt.t.resize(skua_opt.current_points);
		dataFile >> d_read; skua_opt.skua_dat.magpie_dat.sys_dat.T = d_read;
		skua_opt.current_temp = d_read;
		dataFile >> d_read; skua_opt.skua_dat.y[skua_opt.adsorb_index] = d_read / skua_opt.skua_dat.magpie_dat.sys_dat.PT;
		skua_opt.current_press = d_read;
		dataFile >> d_read; skua_opt.current_equil = d_read;
		
		//Set up the gas molefractions
		success = SKUA_OPT_set_y(&skua_opt);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Establish the final equilibrium point
		skua_opt.skua_dat.magpie_dat.sys_dat.Recover = false;
		if (skua_opt.skua_dat.magpie_dat.sys_dat.N > 1)
			skua_opt.skua_dat.magpie_dat.sys_dat.Carrier = true;
		else
			skua_opt.skua_dat.magpie_dat.sys_dat.Carrier = false;
		success = MAGPIE(&skua_opt.skua_dat.magpie_dat);
		if (success < 0 || success > 3) {mError(simulation_fail); return -1;}
		else success = 0;
		skua_opt.total_eval = skua_opt.total_eval + skua_opt.skua_dat.magpie_dat.sys_dat.total_eval;
		skua_opt.simulation_equil = skua_opt.skua_dat.magpie_dat.sys_dat.qT;
		
		//Read in the curve data points
		double maxsum = 0.0, minsum = 0.0;
		for (int n=0; n<skua_opt.current_points; n++)
		{
			dataFile >> d_read; skua_opt.t[n] = d_read;
			dataFile >> d_read; skua_opt.q_data[n] = d_read;
			if (n > 0)
			{
				maxsum = maxsum + (1.0 - (skua_opt.q_data[n]/skua_opt.current_equil));
				minsum = minsum + (0.0 - (skua_opt.q_data[n]/skua_opt.current_equil));
			}
		}
		skua_opt.max_bias = maxsum;
		skua_opt.min_bias = minsum;
		
		//Perform initial guess routine
		std::cout << "\nPerforming operations for curve " << (curve+1) << " . Stand by:\n";
		skua_opt.param_guess = skua_opt.param_guess_old;
		success = initial_guess_SKUA(&skua_opt);
		if (success != 0) {mError(simulation_fail); return -1;}
		par[0] = skua_opt.param_guess;
		status.fnorm = skua_opt.e_norm;
		
		//Call the LMFIT routine
		lmmin(skua_opt.num_params, par, skua_opt.current_points, (void *)&skua_opt, eval_SKUA_Uptake, &control, &status, lm_printout_std);
		skua_opt.total_eval = skua_opt.total_eval + status.nfev;
		skua_opt.param_guess_old = fabs(par[0]);
		printf("\nStatus after %d function evaluations:\n  %s\n\n",status.nfev,lm_infmsg[status.info]);
		std::cout << "E.Norm: "<< status.fnorm << std::endl << std::endl;
		if (status.info > 3 && status.info != 11) {mError(simulation_fail); return -1;}
		
		//Header file for curve results
		if (skua_opt.Optimize == true)
			fprintf(ParamResults, "%i\t%.6g\t%.6g\t%.6g\t%.6g\n",(curve+1),skua_opt.current_temp,skua_opt.current_press,fabs(par[0]),status.fnorm);
		else
			fprintf(ParamResults, "%i\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n",(curve+1),skua_opt.current_temp,skua_opt.current_press,skua_opt.skua_dat.param_dat[skua_opt.adsorb_index].ref_diffusion,skua_opt.skua_dat.param_dat[skua_opt.adsorb_index].activation_energy,status.fnorm);
		fprintf(Comparison,"Curve\t%i\np[kPa]\t%.6g\tT[K]\t%.6g\n",(curve+1),skua_opt.current_press,skua_opt.current_temp);
		fprintf(Comparison,"Time[hr]\tData(Raw)\tData(Normal)\tModel(Raw)\tModel(Normal)\n");
		
		//Print out current results
		for (int n=0; n<skua_opt.current_points; n++)
		{
			fprintf(Comparison,"%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t",skua_opt.t[n],skua_opt.q_data[n],(skua_opt.q_data[n] / skua_opt.current_equil),skua_opt.q_sim[n],(skua_opt.q_sim[n]/skua_opt.simulation_equil));
			fprintf(Comparison,"\n");
		}
		fprintf(Comparison,"------------------------------------------------------\n\n");
		
		skua_opt.evaluation = 0;
		curve++;
	} while (curve < skua_opt.num_curves /*1*/);
	dataFile.close();
	
	//END of optimization
	fclose(ParamResults);
	fclose(Comparison);
	time = clock() - time;
	std::cout << "Optimization Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Evaluations: " << skua_opt.total_eval << "\n";
	std::cout << "Evaluations/sec: " << skua_opt.total_eval/(time / CLOCKS_PER_SEC) << "\n";
	
	return success;
}
