//----------------------------------------
//  Created by Austin Ladshaw on 4/9/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 *			DOGFISH = Diffusion Object Governing Fiber Interior Sorption History
 *
 */

#include "dogfish.h"

//Header file for output
void print2file_species_header(FILE *Output, DOGFISH_DATA *dog_dat, int i)
{
	const char *name = dog_dat->param_dat[i].species.MolecularFormula().c_str();
	fprintf(Output, "%s\t", name);
	if (dog_dat->finch_dat[i].Dirichlet == true)
		fprintf(Output,"-\t");
	for (int l=0; l<dog_dat->finch_dat[i].LN+2; l++)
	{
		fprintf(Output,"-\t");
	}
}

//Print out the standard DOGFISH header to the output file
void print2file_DOGFISH_header(DOGFISH_DATA *dog_dat)
{
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		print2file_species_header(dog_dat->OutputFile, dog_dat, i);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
	}
	print2file_newline(dog_dat->OutputFile, NULL);
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		print2file_dim_header(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
	}
	print2file_newline(dog_dat->OutputFile, NULL);
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		print2file_time_header(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
	}
	print2file_newline(dog_dat->OutputFile, NULL);
}

//Printout old time result to file
void print2file_DOGFISH_result_old(DOGFISH_DATA *dog_dat)
{
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		print2file_result_old(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
	}
	print2file_newline(dog_dat->OutputFile, NULL);
}

//Printout new time result to file
void print2file_DOGFISH_result_new(DOGFISH_DATA *dog_dat)
{
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		print2file_result_new(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
		print2file_tab(dog_dat->OutputFile, &dog_dat->finch_dat[i]);
	}
	print2file_newline(dog_dat->OutputFile, NULL);
}

//Default Retardation function
double default_Retardation(int i, int l, const void *data)
{
	return 1.0;
}

//Default Intraparticle Diffusion function
double default_IntraDiffusion(int i, int l, const void *data)
{
	DOGFISH_DATA *dat = (DOGFISH_DATA *) data;
	return dat->param_dat[i].intraparticle_diffusion;
}

//Default Film Mass Transfer function
double default_FilmMTCoeff(int i, const void *data)
{
	DOGFISH_DATA *dat = (DOGFISH_DATA *) data;
	return dat->param_dat[i].film_transfer_coeff;
}

//Default Surface Concentration function
double default_SurfaceConcentration(int i, const void *data)
{
	DOGFISH_DATA *dat = (DOGFISH_DATA *) data;
	return dat->param_dat[i].surface_concentration;
}

//Setup function
int setup_DOGFISH_DATA(FILE *file, double (*eval_R) (int i, int l, const void *user_data),
					   double (*eval_DI) (int i, int l, const void *user_data),
					   double (*eval_kf) (int i, const void *user_data),
					   double (*eval_qs) (int i, const void *user_data),
					   const void *user_data, DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	if (file == NULL)
		dog_dat->Print2File = false;
	else
	{
		dog_dat->Print2File = true;
		dog_dat->OutputFile = file;
	}
	if (eval_R == NULL)
		dog_dat->eval_R = (*default_Retardation);
	else
		dog_dat->eval_R = (*eval_R);
	if (eval_DI == NULL)
		dog_dat->eval_DI = (*default_IntraDiffusion);
	else
		dog_dat->eval_DI = (*eval_DI);
	if (eval_kf == NULL)
		dog_dat->eval_kf = (*default_FilmMTCoeff);
	else
		dog_dat->eval_kf = (*eval_kf);
	if (eval_qs == NULL)
		dog_dat->eval_qs = (*default_SurfaceConcentration);
	else
		dog_dat->eval_qs = (*eval_qs);
	dog_dat->user_data = user_data;
	
	//Setup FINCH data and remaining info
	dog_dat->time = 0.0;
	dog_dat->time_old = 0.0;
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		dog_dat->finch_dat[i].t = dog_dat->time;
		dog_dat->finch_dat[i].t_old = dog_dat->time_old;
		dog_dat->finch_dat[i].L = dog_dat->fiber_diameter / 2.0;		//um
		dog_dat->finch_dat[i].T = dog_dat->end_time;					//hrs
		dog_dat->finch_dat[i].Dirichlet = dog_dat->DirichletBC;
		dog_dat->finch_dat[i].LN = 10;
		dog_dat->finch_dat[i].SteadyState = false;
		dog_dat->finch_dat[i].CheckMass = false;
		dog_dat->finch_dat[i].Iterative = dog_dat->NonLinear;
		dog_dat->finch_dat[i].NormTrack = dog_dat->Print2Console;
		dog_dat->finch_dat[i].nl_method = 0;
		dog_dat->finch_dat[i].d = 1;
		dog_dat->finch_dat[i].s = dog_dat->fiber_length;
		dog_dat->finch_dat[i].t = 0.0;
		dog_dat->finch_dat[i].t_old = 0.0;
		dog_dat->finch_dat[i].Ro = 1.0;
		dog_dat->finch_dat[i].RIC = 1.0;
		dog_dat->finch_dat[i].vo = 0.0;
		dog_dat->finch_dat[i].vIC = 0.0;
		dog_dat->finch_dat[i].kIC = 0.0;
		dog_dat->finch_dat[i].ko = 0.0;
		
		success = setup_FINCH_DATA(NULL,NULL,NULL,NULL,NULL,set_DOGFISH_params,NULL,NULL,NULL,NULL,NULL,NULL,&dog_dat->finch_dat[i],user_data);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	
	return success;
}

int DOGFISH_Executioner(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	//Perform Preprocess Actions
	success = DOGFISH_preprocesses(dog_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Loop for all components
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		//Solve the system
		success = (*dog_dat->finch_dat[i].solve) ((void *)&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
			
		//Check for negative concentrations
		success = check_Mass(&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
			
		//Form Totals
		success = uTotal(&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
			
		//Form Averages
		success = uAverage(&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	
	//Perform Postprocess Actions
	success = DOGFISH_postprocesses(dog_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	return success;
}

int set_DOGFISH_ICs(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	//This is where we could call any speciation functions to establish starting point
	
	//Set ICs for FINCH sub-problems
	dog_dat->total_sorption = dog_dat->total_sorption_old;
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		//Non-function based ICs
		dog_dat->finch_dat[i].Sn.ConstantICFill(0.0);
		dog_dat->finch_dat[i].Snp1.ConstantICFill(0.0);
		dog_dat->finch_dat[i].vn.ConstantICFill(0.0);
		dog_dat->finch_dat[i].vnp1.ConstantICFill(0.0);
		dog_dat->finch_dat[i].kn.ConstantICFill(0.0);
		dog_dat->finch_dat[i].knp1.ConstantICFill(0.0);
		dog_dat->finch_dat[i].uo = dog_dat->param_dat[i].initial_sorption;
		dog_dat->param_dat[i].sorbed_molefraction = dog_dat->param_dat[i].initial_sorption / dog_dat->total_sorption;
		dog_dat->finch_dat[i].un.ConstantICFill(dog_dat->finch_dat[i].uo);
		dog_dat->finch_dat[i].unm1 = dog_dat->finch_dat[i].un;
		dog_dat->finch_dat[i].unp1 = dog_dat->finch_dat[i].un;
	}
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		//Function based ICs
		dog_dat->finch_dat[i].Ro = (*dog_dat->eval_R) (i, -1, dog_dat->user_data);
		dog_dat->finch_dat[i].Do = (*dog_dat->eval_DI) (i, -1, dog_dat->user_data);
		dog_dat->finch_dat[i].kfn = (*dog_dat->eval_kf) (i, dog_dat->user_data);
		for (int l=0; l<dog_dat->finch_dat[i].LN; l++)
		{
			dog_dat->finch_dat[i].Rn(l,0) = (*dog_dat->eval_R) (i, l, dog_dat->user_data);
			dog_dat->finch_dat[i].Dn(l,0) = (*dog_dat->eval_DI) (i, l, dog_dat->user_data);
		}
		dog_dat->finch_dat[i].Rnp1 = dog_dat->finch_dat[i].Rn;
		dog_dat->finch_dat[i].Dnp1 = dog_dat->finch_dat[i].Dn;
		dog_dat->finch_dat[i].kfnp1 = dog_dat->finch_dat[i].kfn;
		
		//Form Totals
		success = uTotal(&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Form Averages
		success = uAverage(&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Update averages and totals
		dog_dat->finch_dat[i].uT_old = dog_dat->finch_dat[i].uT;
		dog_dat->finch_dat[i].uAvg_old = dog_dat->finch_dat[i].uAvg;
	}
	
	return success;
}

int set_DOGFISH_timestep(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		dog_dat->finch_dat[i].dt = dog_dat->finch_dat[i].dz / 4.0;
		if (dog_dat->finch_dat[i].dt >= 1.0)
			dog_dat->finch_dat[i].dt = 0.5;
		dog_dat->finch_dat[i].t = dog_dat->finch_dat[i].dt + dog_dat->finch_dat[i].t_old;
	}
	dog_dat->time_old = dog_dat->finch_dat[0].t_old;
	dog_dat->time = dog_dat->finch_dat[0].t;
	
	return success;
}

int DOGFISH_preprocesses(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	//Want to use this function to establish BCs
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		dog_dat->finch_dat[i].uo = (*dog_dat->eval_qs) (i, dog_dat->user_data);
		dog_dat->finch_dat[i].kfnp1 = (*dog_dat->eval_kf) (i, dog_dat->user_data);
		dog_dat->finch_dat[i].Ro = (*dog_dat->eval_R) (i, -1, dog_dat->user_data);
		dog_dat->finch_dat[i].Do = (*dog_dat->eval_DI) (i, -1, dog_dat->user_data);
	}
	
	return success;
}

int set_DOGFISH_params(const void *user_data)
{
	int success = 0;
	DOGFISH_DATA *dog_dat = (DOGFISH_DATA *) user_data;
	
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		for (int l=0; l<dog_dat->finch_dat[i].LN; l++)
		{
			dog_dat->finch_dat[i].Rnp1(l,0) = (*dog_dat->eval_R) (i, l, user_data);
			dog_dat->finch_dat[i].Dnp1(l,0) = (*dog_dat->eval_DI) (i, l, user_data);
		}
	}
	
	return success;
}

int DOGFISH_postprocesses(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	dog_dat->total_sorption = 0.0;
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		dog_dat->total_sorption+=dog_dat->finch_dat[i].uAvg;
		dog_dat->total_steps+=dog_dat->finch_dat[i].total_iter;
		dog_dat->finch_dat[i].total_iter = 0;
	}
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		if (dog_dat->total_sorption == 0.0)
			dog_dat->param_dat[i].sorbed_molefraction = 0.0;
		else
			dog_dat->param_dat[i].sorbed_molefraction = dog_dat->finch_dat[i].uAvg / dog_dat->total_sorption;
	}
	
	//Print results to output file
	if (dog_dat->Print2File == true)
	{
		dog_dat->t_counter = dog_dat->t_counter + dog_dat->finch_dat[0].dt;
		if (dog_dat->t_counter >= dog_dat->t_print)
		{
			print2file_DOGFISH_result_new(dog_dat);
			dog_dat->t_counter = 0.0;
		}
	}
	
	return success;
}

int DOGFISH_reset(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	for (int i=0; i<dog_dat->NumComp; i++)
	{
		success = (*dog_dat->finch_dat[i].resettime) ((void *)&dog_dat->finch_dat[i]);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	dog_dat->total_sorption_old = dog_dat->total_sorption;
	if (dog_dat->time_old == 0.0)
		dog_dat->time_old = 0.0;
	else
		dog_dat->time_old = dog_dat->time;
	
	return success;
}

int DOGFISH(DOGFISH_DATA *dog_dat)
{
	int success = 0;
	
	//Print to file
	if (dog_dat->Print2File == true)
	{
		print2file_DOGFISH_header(dog_dat);
	}
	
	//Set Initial Conditions
	success = set_DOGFISH_ICs(dog_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Print ICs
	if (dog_dat->Print2File == true)
	{
		print2file_DOGFISH_result_old(dog_dat);
	}

	
	//Solve a series of time steps implicitly
	do
	{
		if (dog_dat->finch_dat[0].Update == true)
		{
			success = DOGFISH_reset(dog_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}
		success = set_DOGFISH_timestep(dog_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
		std::cout << "Evaluating time: " << dog_dat->time << " hours..." << std::endl;
		success = DOGFISH_Executioner(dog_dat);
		if (success == 0)
		{
			std::cout << "Simulation Successful!\n" << std::endl;
			for (int i=0; i<dog_dat->NumComp; i++)
				dog_dat->finch_dat[i].Update = true;
		}
		else {mError(simulation_fail); dog_dat->finch_dat[0].Update = false; return -1;}
		dog_dat->total_steps++;
	} while (dog_dat->time < dog_dat->end_time);
	
	return success;
}


//Testing Suite for DOGFISH
int DOGFISH_TESTS()
{
	int success = 0;
	double time;
	DOGFISH_DATA dog_dat;
	FILE *TestOutput;
	std::cout << "Start of DOGFISH Tests\n\n";
	
	//Standard initializations
	time = clock();
	TestOutput = fopen("output/DOGFISH_TestOutput.txt","w+");
	if (TestOutput == nullptr)
	{
		system("mkdir output");
		TestOutput = fopen("output/DOGFISH_TestOutput.txt","w+");
	}
	dog_dat.total_steps = 0;
	dog_dat.time_old = 0.0;
	dog_dat.time = 0.0;
	
	//This is where some input file would be read
	dog_dat.NumComp = 1;
	dog_dat.DirichletBC = true;
	dog_dat.NonLinear = true;
	dog_dat.param_dat.resize(dog_dat.NumComp);
	dog_dat.finch_dat.resize(dog_dat.NumComp);
	dog_dat.end_time = 1500.0;												//hours
	dog_dat.t_print = dog_dat.end_time / 1000.0;
	dog_dat.total_sorption_old = 0.0;										//mg/g (Total IC)
	dog_dat.fiber_length = 1000.0;											//um
	dog_dat.fiber_diameter = 76.5 * 2.0;									//um
	double check = 0.0;
	for (int i=0; i<dog_dat.NumComp; i++)
	{
		dog_dat.param_dat[i].sorbed_molefraction = 1.0/dog_dat.NumComp;
		dog_dat.param_dat[i].intraparticle_diffusion = 2.148;				//um^2/hr
		dog_dat.param_dat[i].film_transfer_coeff = 1.0;						//um/hr
		dog_dat.param_dat[i].surface_concentration = 1.0;					//mg/g
		check+=dog_dat.param_dat[i].sorbed_molefraction;					//-
		dog_dat.param_dat[i].initial_sorption = 0.0;						//mg/g (individual IC)
	}
	if (check > 1.0+1e-6 || check < 1.0-1e-6)
	{
		std::cout << "check sould be 1, but equals " << check << std::endl;
		mError(invalid_solid_sum);
		return -1;
	}
	
	//Testing Purposes (formula's to be determined by input file or speciation)
	dog_dat.param_dat[0].species.Register("UO2 2+ (aq)");
	
	//Call the setup function
	success = setup_DOGFISH_DATA(TestOutput, default_Retardation, default_IntraDiffusion, default_FilmMTCoeff, default_SurfaceConcentration, (void *)&dog_dat, &dog_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Call the routine
	success = DOGFISH(&dog_dat);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	
	fclose(TestOutput);
	time = clock() - time;
	std::cout << "Simulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Evaluations: " << dog_dat.total_steps << "\n";
	std::cout << "Evaluations/sec: " << dog_dat.total_steps/(time / CLOCKS_PER_SEC) << "\n";
	std::cout << "\nEnd of DOGFISH Tests\n\n";
	return success;
}