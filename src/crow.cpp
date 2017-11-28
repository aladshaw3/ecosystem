/*!
 *  \file crow.h
 *	\brief Coupled Reaction Object Workspace
 *  \author Austin Ladshaw
 *	\date 11/27/2017
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2017, all
 *             rights reserved.
 */

#include "crow.h"

/*
 *								Start: ConstReaction
 *	-------------------------------------------------------------------------------------
 */

//Constructor
ConstReaction::ConstReaction()
{
	forward_rate = 1.0;
	reverse_rate = 1.0;
	time_coeff = 1.0;
	main_index = 0;
}

//Destructor
ConstReaction::~ConstReaction()
{
	stoic.clear();
}

//Function to initialize solver
void ConstReaction::InitializeSolver(Dove &Solver)
{
	this->SolverInfo = &Solver;
}

//Set function/variable index
void ConstReaction::SetIndex(int index)
{
	this->main_index = index;
}

//Set forward rate
void ConstReaction::SetForwardRate(double rate)
{
	this->forward_rate = rate;
}

//Set reverse rate
void ConstReaction::SetReverseRate(double rate)
{
	this->reverse_rate = rate;
}

//Insert stoichiometry
void ConstReaction::InsertStoichiometry(int i, int v)
{
	this->stoic[i] = v;
}

//Compute time coeff
void ConstReaction::ComputeTimeCoeff()
{
	try
	{
		this->time_coeff = 1.0/((double)this->stoic.at(this->main_index));
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		this->time_coeff = 1.0;
	}
}

//Return forward rate
double ConstReaction::getForwardRate()
{
	return this->forward_rate;
}

//Return reverse rate
double ConstReaction::getReverseRate()
{
	return this->reverse_rate;
}

//Return time coefficient
double ConstReaction::getTimeCoeff()
{
	return this->time_coeff;
}

//Return variable index
int ConstReaction::getIndex()
{
	return this->main_index;
}

//Return reference to map
std::map<int,int> & ConstReaction::getStoichiometryMap()
{
	return this->stoic;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: ConstReaction
 */

//Time coeff function
double coeff_func_ConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	CROW_DATA *dat = (CROW_DATA *) data;
	return dat->const_reacts[i].getTimeCoeff();
}

//Rate function
double rate_func_ConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double rate = 0.0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, int>::iterator it;
	
	//Iterate forward through the map
	int count = 0;
	double forward = 0.0, reverse = 0.0;
	bool first_prod = true, first_reac = true;
	for (it = dat->const_reacts[i].getStoichiometryMap().begin(); it != dat->const_reacts[i].getStoichiometryMap().end(); it++)
	{
		if (it->second > 0)
		{
			if (first_prod == true)
			{
				reverse = pow(u(it->first,0), (double)it->second);
				first_prod = false;
			}
			else
				reverse = reverse * pow(u(it->first,0), (double)it->second);
		}
		else if (it->second < 0)
		{
			if (first_reac == true)
			{
				forward = pow(u(it->first,0), (double)abs(it->second));
				first_reac = false;
			}
			else
				forward = forward * pow(u(it->first,0), (double)abs(it->second));
		}
		else
		{
			//Do Nothing
		}
		count++;
	}
	
	//Check the sign of the variable of interest's stoichiometry
	if (dat->const_reacts[i].getStoichiometryMap()[dat->const_reacts[i].getIndex()] > 0)
	{
		rate = dat->const_reacts[i].getForwardRate()*forward - dat->const_reacts[i].getReverseRate()*reverse;
	}
	else if (dat->const_reacts[i].getStoichiometryMap()[dat->const_reacts[i].getIndex()] < 0)
	{
		rate = dat->const_reacts[i].getReverseRate()*reverse - dat->const_reacts[i].getForwardRate()*forward;
	}
	else
	{
		rate = 0.0;
	}
	
	return rate;
}

//Jacobi function
double jacobi_func_ConstReaction(int i, int j, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double jac = 0.0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, int>::iterator it;
	
	//Assume j is a variable related to the ith function
	try
	{
		//Check stoichiometry of j
		bool prod = false;
		if (dat->const_reacts[i].getStoichiometryMap().at(j) > 0)
			prod = true;
		
		//Iterate forward through the map
		int count = 0;
		double forward = 0.0, reverse = 0.0;
		bool first_prod = true, first_reac = true;
		
		if (prod == true)
		{
			for (it = dat->const_reacts[i].getStoichiometryMap().begin(); it != dat->const_reacts[i].getStoichiometryMap().end(); it++)
			{
				if (it->second > 0)
				{
					if (first_prod == true)
					{
						if (it->first == j)
							reverse = pow(u(it->first,0), (double)(it->second-1))*(double)(it->second);
						else
							reverse = pow(u(it->first,0), (double)it->second);
						first_prod = false;
					}
					else
						if (it->first == j)
							reverse = reverse * pow(u(it->first,0), (double)(it->second-1))*(double)(it->second);
						else
							reverse = reverse * pow(u(it->first,0), (double)it->second);
				}
				count++;
			}
		}
		else
		{
			for (it = dat->const_reacts[i].getStoichiometryMap().begin(); it != dat->const_reacts[i].getStoichiometryMap().end(); it++)
			{
				if (it->second < 0)
				{
					if (first_reac == true)
					{
						if (it->first == j)
							forward = pow(u(it->first,0), (double)(abs(it->second)-1))*(double)abs(it->second);
						else
							forward = pow(u(it->first,0), (double)abs(it->second));
						first_reac = false;
					}
					else
						if (it->first == j)
							forward = forward * pow(u(it->first,0), (double)(abs(it->second)-1))*(double)abs(it->second);
						else
							forward = forward * pow(u(it->first,0), (double)abs(it->second));
				}
				count++;
			}
		}
		
		//Compute Jac from relevant info
		if (dat->const_reacts[i].getStoichiometryMap()[dat->const_reacts[i].getIndex()] > 0)
		{
			jac = dat->const_reacts[i].getForwardRate()*forward - dat->const_reacts[i].getReverseRate()*reverse;
		}
		else if (dat->const_reacts[i].getStoichiometryMap()[dat->const_reacts[i].getIndex()] < 0)
		{
			jac = dat->const_reacts[i].getReverseRate()*reverse - dat->const_reacts[i].getForwardRate()*forward;
		}
		else
			jac = 0.0;

	}
	//If it is not related, then return 0
	catch (std::out_of_range)
	{
		jac = 0.0;
	}
	
	return jac;
}

//Execute CROW
int CROW_SCENARIO(const char *yaml_input)
{
	int success = 0;
	
	std::cout << "This is an executable\n";
	
	return success;
}

//Crow test
int CROW_TESTS()
{
	int success = 0;
	
	FILE *file;
	file = fopen("output/CROW_Tests.txt", "w+");
	if (file == nullptr)
	{
		system("mkdir output");
		file = fopen("output/CROW_Tests.txt", "w+");
	}
	
	CROW_DATA test01;
	test01.SolverInfo.set_numfunc(1);
	test01.SolverInfo.set_outputfile(file);
	
	test01.const_reacts.resize(1);
	test01.const_reacts[0].InitializeSolver(test01.SolverInfo);
	
	test01.SolverInfo.set_output(true);
	test01.SolverInfo.set_NonlinearOutput(true);
	test01.SolverInfo.set_endtime(1);
	test01.SolverInfo.set_userdata((void*)&test01);
	test01.SolverInfo.set_tolerance(1e-6);
	test01.SolverInfo.set_fileoutput(true);
	test01.SolverInfo.set_timestepmax(0.5);
	test01.SolverInfo.set_timestepper(CONSTANT);
	test01.SolverInfo.set_variableName(0, "u");
	
	test01.const_reacts[0].SetIndex(0);
	test01.const_reacts[0].SetForwardRate(0);
	test01.const_reacts[0].SetReverseRate(1);
	
	test01.const_reacts[0].InsertStoichiometry(0, 1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[0].ComputeTimeCoeff(); //MUST CALL THIS BEFORE SOLVE!!!
	
	test01.SolverInfo.registerCoeff(0, coeff_func_ConstReaction);
	test01.SolverInfo.registerFunction(0, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(0, 0, jacobi_func_ConstReaction);
	
	test01.SolverInfo.set_timestep(0.05);
	test01.SolverInfo.set_integrationtype(BE);
	test01.SolverInfo.set_initialcondition("u", 1);
	
	test01.SolverInfo.solve_all();
	
	std::cout << "This is a test\n";
	
	return success;
}
