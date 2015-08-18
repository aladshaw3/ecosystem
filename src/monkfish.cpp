//----------------------------------------
//  Created by Austin Ladshaw on 4/14/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 *			MONKFISH = Multi-fiber wOven Nest Kernel For Interparticle Sorption History
 *
 */

#include "monkfish.h"

//Default porosity function
double default_porosity(int i, int l, const void *user_data)
{
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	if (l < 0)
		return monk_dat->max_porosity;
	else
		return monk_dat->min_porosity + ( ((double)l*monk_dat->finch_dat[i].dz) * ( (monk_dat->max_porosity-monk_dat->min_porosity) / (monk_dat->domain_diameter/2.0) ) );
}

//Default density function
double default_density(int i, int l, const void *user_data)
{
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	return monk_dat->single_fiber_density * (1.0 - (*monk_dat->eval_eps) (i, l, user_data));
}

//Default interparticle diffusion function
double default_interparticle_diffusion(int i, int l, const void *user_data)
{
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	double eps = (*monk_dat->eval_eps) (i, l, user_data);
	return eps*eps*eps*monk_dat->param_dat[i].interparticle_diffusion;
}

//Deafault adsorption function {NOTE: All of this will probably have to change!!!}
double default_monk_adsorption(int i, int l, const void *user_data)
{
	int success = 0;
	double dq_dc = 0.0;
	double eps = sqrt(DBL_EPSILON);
	double qEPS = 0.0;
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	
	//Boundary value
	if (l < 0)
	{
		//Store real exterior concentration locally, then perturb with eps
		double uo_temp = monk_dat->finch_dat[i].uo;
		monk_dat->finch_dat[i].uo = uo_temp + eps;
		
		//Call equilibrium/kinetics at perturbation and store result
		qEPS = default_monk_equilibrium(i, -1, user_data);
		
		//Restore boundary condition and recall equilibrium
		monk_dat->finch_dat[i].uo = uo_temp;
		monk_dat->param_dat[i].sorption_bc = default_monk_equilibrium(i, -1, user_data);
		
		//Approximate adsorption strength
		dq_dc = (qEPS-monk_dat->param_dat[i].sorption_bc)/eps;
	}
	//Nodal value
	else
	{
		//Based only on equilibrium
		if (monk_dat->MultiScale == false)
		{
			//Store real concentration locally, then perturb with eps
			double u_temp = monk_dat->finch_dat[i].unp1(l,0);
			monk_dat->finch_dat[i].unp1.edit(l,0,u_temp + eps);
			
			//Call equilibrium/kinetics at perturbation and store result
			qEPS = default_monk_equilibrium(i, l, user_data);
			
			//Restore boundary condition and recall equilibrium
			monk_dat->finch_dat[i].unp1.edit(l, 0, u_temp);
			monk_dat->param_dat[i].avg_sorption.edit(l, 0, default_monk_equilibrium(i, l, user_data));
			
			//Approximate adsorption strength
			dq_dc = (qEPS-monk_dat->param_dat[i].avg_sorption(l,0))/eps;
		}
		//Based on DOGFISH
		else
		{
			//Store real concentration locally, then perturb with eps
			double u_temp = monk_dat->finch_dat[i].unp1(l,0);
			monk_dat->finch_dat[i].unp1.edit(l,0,u_temp + eps);
			
			//Establish information necessary to run a DOGFISH simulation
			monk_dat->dog_dat[l].param_dat[i].intraparticle_diffusion = monk_dat->param_dat[i].intraparticle_diffusion;
			monk_dat->dog_dat[l].param_dat[i].film_transfer_coeff = monk_dat->param_dat[i].film_transfer_coeff;
			monk_dat->dog_dat[l].param_dat[i].surface_concentration = default_monk_equilibrium(i, l, user_data);
			
			//Call DOGFISH and store result
			monk_dat->dog_dat[l].total_steps = 0;
			success = DOGFISH_Executioner(&monk_dat->dog_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
			monk_dat->total_steps = monk_dat->total_steps + monk_dat->dog_dat[l].total_steps;
			qEPS = monk_dat->dog_dat[l].finch_dat[i].uAvg;
			
			//Restore concentration and setup DOGFISH
			monk_dat->finch_dat[i].unp1.edit(l, 0, u_temp);
			monk_dat->dog_dat[l].param_dat[i].intraparticle_diffusion = monk_dat->param_dat[i].intraparticle_diffusion;
			monk_dat->dog_dat[l].param_dat[i].film_transfer_coeff = monk_dat->param_dat[i].film_transfer_coeff;
			monk_dat->dog_dat[l].param_dat[i].surface_concentration = default_monk_equilibrium(i, l, user_data);
			
			//Call DOGFISH and store result
			monk_dat->dog_dat[l].total_steps = 0;
			success = DOGFISH_Executioner(&monk_dat->dog_dat[l]);
			if (success != 0) {mError(simulation_fail); return -1;}
			monk_dat->total_steps = monk_dat->total_steps + monk_dat->dog_dat[l].total_steps;
			monk_dat->param_dat[i].avg_sorption.edit(l, 0, monk_dat->dog_dat[l].finch_dat[i].uAvg);
			
			//Approximate adsorption strength
			dq_dc = (qEPS-monk_dat->param_dat[i].avg_sorption(l,0))/eps;
		}
	}
	
	return dq_dc;
}

//Default adsorption equilibrium function (returns the equilibrium adsorption in mg/g) {Will have to be rewritten}
double default_monk_equilibrium(int i, int l, const void *user_data)
{
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	double equ = 0.0;
	
	if (l < 0)
		equ = monk_dat->finch_dat[i].uo / ((*monk_dat->eval_rho) (i, -1, user_data)) * monk_dat->param_dat[i].species.MolarWeight() * 1000.0;
	else
		equ = monk_dat->finch_dat[i].uo / ((*monk_dat->eval_rho) (i, l, user_data)) * monk_dat->param_dat[i].species.MolarWeight() * 1000.0;
	
	return equ;
}

//Default retardation fuction
double default_monkfish_retardation(int i, int l, const void *user_data)
{
	double Ret = 0.0;
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	if (l < 0)
		Ret = (*monk_dat->eval_eps) (i, l, user_data);
	else
		Ret = (*monk_dat->eval_eps) (i, l, user_data) + ( (*monk_dat->eval_rho) (i, l, user_data) * (*monk_dat->eval_ads) (i, l, user_data));
	
	if (Ret < 0.0)
	{
		mError(unstable_matrix);
		Ret = fabs(Ret);
	}
	
	return Ret;
}

//Default exterior concentration function
double default_exterior_concentration(int i, const void *user_data)
{
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	return monk_dat->param_dat[i].exterior_concentration;
}

//Default film mass transfer function
double default_film_transfer(int i, const void *user_data)
{
	MONKFISH_DATA *monk_dat = (MONKFISH_DATA *) user_data;
	return monk_dat->param_dat[i].exterior_transfer_coeff;
}

//Testing of MONKFISH
int MONKFISH_TESTS()
{
	int success = 0;
	
	return success;
}
