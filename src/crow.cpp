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
		this->time_coeff = 1.0/fabs((double)this->stoic.at(this->main_index));
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

//Print header to file
void print2file_crow_header(CROW_DATA *dat)
{
	if (dat->FileOutput != true)
		return;
	
	fprintf(dat->OutputFile, "------------------- System Information ------------------\n\n");
	fprintf(dat->OutputFile, "Number of Variables/Functions = \t%i\n",dat->SolverInfo.getNumFunc());
	fprintf(dat->OutputFile, "List of Variable Names:");
	for (int i=0; i<dat->SolverInfo.getNumFunc(); i++)
		fprintf(dat->OutputFile, "\n\t%s",dat->SolverInfo.getVariableName(i).c_str());
	fprintf(dat->OutputFile, "\n");
	fprintf(dat->OutputFile,"\nDOVE Solver Options:\n");
	fprintf(dat->OutputFile, "\tINTEGRATION TYPE = \t");
	switch (dat->SolverInfo.getIntegrationType())
	{
		case EXPLICIT:
			fprintf(dat->OutputFile, "EXPLICIT\n");
			break;
			
		case IMPLICIT:
			fprintf(dat->OutputFile, "IMPLICIT\n");
			break;
	}
	fprintf(dat->OutputFile, "\tINTEGRATION METHOD = \t");
	switch (dat->SolverInfo.getIntegrationMethod())
	{
		case BE:
			fprintf(dat->OutputFile, "Backwards-Euler\n");
			break;
			
		case FE:
			fprintf(dat->OutputFile, "Forwards-Euler\n");
			break;
			
		case CN:
			fprintf(dat->OutputFile, "Crank-Nicholson\n");
			break;
			
		case BDF2:
			fprintf(dat->OutputFile, "Backwards Differentiation Formula: 2nd Order\n");
			break;
			
		case RK4:
			fprintf(dat->OutputFile, "Runge-Kutta: 4th Order\n");
			break;
			
		case RKF:
			fprintf(dat->OutputFile, "Runge-Kutta-Fehlberg\n");
			break;
	}
	fprintf(dat->OutputFile, "\tTIMESTEPPER METHOD = \t");
	switch (dat->SolverInfo.getTimeStepper())
	{
		case CONSTANT:
			fprintf(dat->OutputFile, "Constant\n");
			break;
			
		case ADAPTIVE:
			fprintf(dat->OutputFile, "Adaptive\n");
			break;
			
		case FEHLBERG:
			fprintf(dat->OutputFile, "Fehlberg\n");
			break;
			
		case RATEBASED:
			fprintf(dat->OutputFile, "Ratebased\n");
			break;
	}
	fprintf(dat->OutputFile, "\tLINESEARCH METHOD = \t");
	switch (dat->SolverInfo.getLinesearchMethod())
	{
		case BT:
			fprintf(dat->OutputFile, "Backtracking\n");
			break;
			
		case ABT:
			fprintf(dat->OutputFile, "Adaptive Backtracking\n");
			break;
			
		case NO_LS:
			fprintf(dat->OutputFile, "None\n");
			break;
	}
	fprintf(dat->OutputFile, "\tPRECONDITIONING METHOD = \t");
	if (dat->SolverInfo.isPreconditioned() == false)
		fprintf(dat->OutputFile, "None\n");
	else
	{
		switch (dat->SolverInfo.getPreconditioner())
		{
			case JACOBI:
				fprintf(dat->OutputFile, "Jacobi\n");
				break;
			
			case TRIDIAG:
				fprintf(dat->OutputFile, "Tridiagonal\n");
				break;
			
			case UGS:
				fprintf(dat->OutputFile, "Upper Gauss-Seidel\n");
				break;
				
			case LGS:
				fprintf(dat->OutputFile, "Lower Gauss-Seidel\n");
				break;
				
			case SGS:
				fprintf(dat->OutputFile, "Symmetric Gauss-Seidel\n");
				break;
		}
	}
	fprintf(dat->OutputFile, "\tLINEAR SOLVER = \t");
	switch (dat->SolverInfo.getLinearMethod())
	{
		case GMRESLP:
			fprintf(dat->OutputFile, "Generalized Minimum Residuals (GMRESLP)\n");
			break;
			
		case PCG:
			fprintf(dat->OutputFile, "Conjugate Gradients (PCG)\n");
			break;
			
		case BiCGSTAB:
			fprintf(dat->OutputFile, "Biconjugate Gradients Stabilized (BiCGSTAB)\n");
			break;
			
		case CGS:
			fprintf(dat->OutputFile, "Conjugate Gradients Squared (CGS)\n");
			break;
			
		case FOM:
			fprintf(dat->OutputFile, "Full Orthogonalization (FOM)\n");
			break;
			
		case GMRESRP:
			fprintf(dat->OutputFile, "Generalized Minimum Residuals (GMRESRP)\n");
			break;
			
		case GCR:
			fprintf(dat->OutputFile, "Generalized Conjugate Residuals (GCR)\n");
			break;
			
		case GMRESR:
			fprintf(dat->OutputFile, "Recursive Generalized Minimum Residuals (GMRESR)\n");
			break;
			
		case KMS:
			fprintf(dat->OutputFile, "Kylov Multi-space (KMS)\n");
			break;
			
		case QR:
			fprintf(dat->OutputFile, "QR Factorization (QR)\n");
			break;
	}
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
	double time;
	
	// ---------------------------- Initializations ---------------------------------
	time = clock();
	CROW_DATA test01;
	FILE *file;
	file = fopen("output/CROW_Tests.txt", "w+");
	if (file == nullptr)
	{
		system("mkdir output");
		file = fopen("output/CROW_Tests.txt", "w+");
	}
	test01.SolverInfo.set_outputfile(file);
	test01.OutputFile = file;
	
	// --------------------------- Input Information (Done on Read) -------------------------------
	
	//---System Info---
	int N = 3;
	test01.SolverInfo.set_numfunc(N);
	test01.SolverInfo.set_variableName(0, "I2");
	test01.SolverInfo.set_variableName(1, "Ag0");
	test01.SolverInfo.set_variableName(2, "AgI");
	test01.SolverInfo.set_initialcondition("I2", 1);
	test01.SolverInfo.set_initialcondition("Ag0", 1);
	test01.SolverInfo.set_initialcondition("AgI", 0);
	
	//---Solver Info (may consolidate with system info)---
	test01.FileOutput = true;
	test01.SolverInfo.set_output(true);
	test01.SolverInfo.set_NonlinearOutput(true);
	test01.SolverInfo.set_LinearOutput(false);
	test01.SolverInfo.set_endtime(20);
	test01.SolverInfo.set_userdata((void*)&test01);
	test01.SolverInfo.set_tolerance(1e-6);
	test01.SolverInfo.set_LinearRelTol(1e-6);
	test01.SolverInfo.set_fileoutput(true);
	test01.SolverInfo.set_timestepmax(0.5);
	test01.SolverInfo.set_timestepper(RATEBASED);
	test01.SolverInfo.set_timestep(0.05);
	test01.SolverInfo.set_integrationtype(BE);
	test01.SolverInfo.set_LinearMethod(GMRESRP);
	test01.SolverInfo.set_preconditioner(SGS);
	test01.SolverInfo.set_Preconditioning(true);
	
	//---Function Info (Specific to ConstReaction)---
	
	//Use emplace_back for each instance of a variable (C++11)
	test01.const_reacts.emplace_back();
	test01.const_reacts[0].InitializeSolver(test01.SolverInfo);
	test01.const_reacts[0].SetIndex(0); //Index should match on this call
	test01.const_reacts[0].SetForwardRate(1);
	test01.const_reacts[0].SetReverseRate(0);
	test01.const_reacts[0].InsertStoichiometry(0, -1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[0].InsertStoichiometry(1, -2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[0].InsertStoichiometry(2, 2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[0].ComputeTimeCoeff(); //MUST CALL THIS BEFORE SOLVE!!!
	test01.SolverInfo.registerCoeff(0, coeff_func_ConstReaction);
	test01.SolverInfo.registerFunction(0, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(0, 0, jacobi_func_ConstReaction); //automate
	test01.SolverInfo.registerJacobi(0, 1, jacobi_func_ConstReaction); //automate
	
	test01.const_reacts.emplace_back();
	test01.const_reacts[1].InitializeSolver(test01.SolverInfo);
	test01.const_reacts[1].SetIndex(1); //Index should match on this call
	test01.const_reacts[1].SetForwardRate(1);
	test01.const_reacts[1].SetReverseRate(0);
	test01.const_reacts[1].InsertStoichiometry(0, -1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[1].InsertStoichiometry(1, -2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[1].InsertStoichiometry(2, 2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[1].ComputeTimeCoeff(); //MUST CALL THIS BEFORE SOLVE!!!
	test01.SolverInfo.registerCoeff(1, coeff_func_ConstReaction);
	test01.SolverInfo.registerFunction(1, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(1, 0, jacobi_func_ConstReaction);
	test01.SolverInfo.registerJacobi(1, 1, jacobi_func_ConstReaction);
	
	test01.const_reacts.emplace_back();
	test01.const_reacts[2].InitializeSolver(test01.SolverInfo);
	test01.const_reacts[2].SetIndex(2); //Index should match on this call
	test01.const_reacts[2].SetForwardRate(1);
	test01.const_reacts[2].SetReverseRate(0);
	test01.const_reacts[2].InsertStoichiometry(0, -1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[2].InsertStoichiometry(1, -2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[2].InsertStoichiometry(2, 2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[2].ComputeTimeCoeff(); //MUST CALL THIS BEFORE SOLVE!!!
	test01.SolverInfo.registerCoeff(2, coeff_func_ConstReaction);
	test01.SolverInfo.registerFunction(2, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(2, 0, jacobi_func_ConstReaction);
	test01.SolverInfo.registerJacobi(2, 1, jacobi_func_ConstReaction);
	
	
	//---Call solver---
	print2file_crow_header(&test01);
	test01.SolverInfo.solve_all();
	
	//---Exit Messages and cleanup---
	time = clock() - time;
	std::cout << "\nCROW Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
