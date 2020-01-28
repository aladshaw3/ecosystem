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
	if (rate < 0)
		this->forward_rate = fabs(rate);
	else
		this->forward_rate = rate;
}

//Set reverse rate
void ConstReaction::SetReverseRate(double rate)
{
	if (rate < 0)
		this->reverse_rate = fabs(rate);
	else
		this->reverse_rate = rate;
}

//Insert stoichiometry
void ConstReaction::InsertStoichiometry(int i, double v)
{
	this->stoic[i] = v;
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

//Return variable index
int ConstReaction::getIndex()
{
	return this->main_index;
}

//Return reference to map
std::map<int,double> & ConstReaction::getStoichiometryMap()
{
	return this->stoic;
}

//Rate function
double rate_func_ConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double rate = 0.0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, double>::iterator it;

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
				reverse = pow(u(it->first,0), it->second);
				first_prod = false;
			}
			else
				reverse = reverse * pow(u(it->first,0), it->second);
		}
		else if (it->second < 0)
		{
			if (first_reac == true)
			{
				forward = pow(u(it->first,0), abs(it->second));
				first_reac = false;
			}
			else
				forward = forward * pow(u(it->first,0), abs(it->second));
		}
		else
		{
			//Do Nothing
		}
		count++;
	}

	//Check the sign of the variable of interest's stoichiometry
	int index = dat->const_reacts[i].getIndex();
	double stoic = dat->const_reacts[i].getStoichiometryMap()[index];
	if (stoic > 0)
	{
		rate = dat->const_reacts[i].getForwardRate()*forward - dat->const_reacts[i].getReverseRate()*reverse;
	}
	else if (stoic < 0)
	{
		rate = dat->const_reacts[i].getReverseRate()*reverse - dat->const_reacts[i].getForwardRate()*forward;
	}
	else
	{
		rate = 0.0;
	}

	return rate*abs(stoic);
}

//Jacobi function
double jacobi_func_ConstReaction(int i, int j, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double jac = 0.0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, double>::iterator it;
	double stoic;

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
							reverse = pow(u(it->first,0), (it->second-1))*(it->second);
						else
							reverse = pow(u(it->first,0), it->second);
						first_prod = false;
					}
					else
						if (it->first == j)
							reverse = reverse * pow(u(it->first,0), (it->second-1))*(it->second);
						else
							reverse = reverse * pow(u(it->first,0), it->second);
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
							forward = pow(u(it->first,0), (abs(it->second)-1))*abs(it->second);
						else
							forward = pow(u(it->first,0), abs(it->second));
						first_reac = false;
					}
					else
						if (it->first == j)
							forward = forward * pow(u(it->first,0), (abs(it->second)-1))*abs(it->second);
						else
							forward = forward * pow(u(it->first,0), abs(it->second));
				}
				count++;
			}
		}

		//Compute Jac from relevant info
		int index = dat->const_reacts[i].getIndex();
		stoic = dat->const_reacts[i].getStoichiometryMap()[index];
		if (stoic > 0)
		{
			jac = dat->const_reacts[i].getForwardRate()*forward - dat->const_reacts[i].getReverseRate()*reverse;
		}
		else if (stoic < 0)
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

	return jac*(double)abs(stoic);
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: ConstReaction
 */

/*
 *								Start: MultiConstReaction
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
MultiConstReaction::MultiConstReaction()
:
reactions(1)
{
	num_reactions = 0;
}

//Default destructor
MultiConstReaction::~MultiConstReaction()
{
	reactions.clear();
}

//Set number of reactions
void MultiConstReaction::SetNumberReactions(unsigned int i)
{
	this->num_reactions = i;
	reactions.resize(i);
}

//Initialize solver
void MultiConstReaction::InitializeSolver(Dove &Solver)
{
	for (int i=0; i<this->num_reactions; i++)
		this->reactions[i].InitializeSolver(Solver);
}

//Set index
void MultiConstReaction::SetIndex(int index)
{
	for (int i=0; i<this->num_reactions; i++)
		this->reactions[i].SetIndex(index);
}

//Set forward rate
void MultiConstReaction::SetForwardRate(int react, double rate)
{
	this->getReaction(react).SetForwardRate(rate);
}

//Set reverse rate
void MultiConstReaction::SetReverseRate(int react, double rate)
{
	this->getReaction(react).SetReverseRate(rate);
}

//Set stoichiometry
void MultiConstReaction::InsertStoichiometry(int react, int i, double v)
{
	this->getReaction(react).InsertStoichiometry(i, v);
}

//Get number of reactions
int MultiConstReaction::getNumReactions()
{
	return this->num_reactions;
}

//Get ConstReaction object
ConstReaction & MultiConstReaction::getReaction(int i)
{
	try
	{
		return this->reactions[i];
	}
	catch (std::out_of_range)
	{
		mError(out_of_bounds);
		return this->reactions[0];
	}
}

/// Rate function for the MultiConstReaction Object
double rate_func_MultiConstReaction(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double rate = 0.0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, double>::iterator it;

	//Loop through all reactions
	for (int r=0; r<dat->multi_const_reacts[i].getNumReactions(); r++)
	{
		//Iterate forward through the map
		int count = 0;
		double forward = 0.0, reverse = 0.0;
		bool first_prod = true, first_reac = true;
		for (it = dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().end(); it++)
		{
			if (it->second > 0)
			{
				if (first_prod == true)
				{
					reverse = pow(u(it->first,0), it->second);
					first_prod = false;
				}
				else
					reverse = reverse * pow(u(it->first,0), it->second);
			}
			else if (it->second < 0)
			{
				if (first_reac == true)
				{
					forward = pow(u(it->first,0), abs(it->second));
					first_reac = false;
				}
				else
					forward = forward * pow(u(it->first,0), abs(it->second));
			}
			else
			{
				//Do Nothing
			}
			count++;
		}

		//Check the sign of the variable of interest's stoichiometry
		int index = dat->multi_const_reacts[i].getReaction(r).getIndex();
		double stoic = dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap()[index];
		if (stoic > 0)
		{
			rate = rate + abs(stoic)*(dat->multi_const_reacts[i].getReaction(r).getForwardRate()*forward - dat->multi_const_reacts[i].getReaction(r).getReverseRate()*reverse);
		}
		else if (stoic < 0)
		{
			rate = rate + abs(stoic)*(dat->multi_const_reacts[i].getReaction(r).getReverseRate()*reverse - dat->multi_const_reacts[i].getReaction(r).getForwardRate()*forward);
		}
		else
		{
			rate = rate + 0.0;
		}
	}

	return rate;
}

/// Jacobi function for the MultiConstReaction Object
double jacobi_func_MultiConstReaction(int i, int j, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double jac = 0.0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, double>::iterator it;
	double stoic;

	//Loop through all reactions
	for (int r=0; r<dat->multi_const_reacts[i].getNumReactions(); r++)
	{
		//Assume j is a variable related to the ith function
		try
		{
			//Check stoichiometry of j
			bool prod = false;
			if (dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().at(j) > 0)
				prod = true;

			//Iterate forward through the map
			int count = 0;
			double forward = 0.0, reverse = 0.0;
			bool first_prod = true, first_reac = true;

			if (prod == true)
			{
				for (it = dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().end(); it++)
				{
					if (it->second > 0)
					{
						if (first_prod == true)
						{
							if (it->first == j)
								reverse = pow(u(it->first,0), (it->second-1))*(it->second);
							else
								reverse = pow(u(it->first,0), it->second);
							first_prod = false;
						}
						else
							if (it->first == j)
								reverse = reverse * pow(u(it->first,0), (it->second-1))*(it->second);
							else
								reverse = reverse * pow(u(it->first,0), it->second);
					}
					count++;
				}
			}
			else
			{
				for (it = dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap().end(); it++)
				{
					if (it->second < 0)
					{
						if (first_reac == true)
						{
							if (it->first == j)
								forward = pow(u(it->first,0), (abs(it->second)-1))*abs(it->second);
							else
								forward = pow(u(it->first,0), abs(it->second));
							first_reac = false;
						}
						else
							if (it->first == j)
								forward = forward * pow(u(it->first,0), (abs(it->second)-1))*abs(it->second);
							else
								forward = forward * pow(u(it->first,0), abs(it->second));
					}
					count++;
				}
			}

			//Compute Jac from relevant info
			int index = dat->multi_const_reacts[i].getReaction(r).getIndex();
			stoic = dat->multi_const_reacts[i].getReaction(r).getStoichiometryMap()[index];
			if (stoic > 0)
			{
				jac += abs(stoic)*(dat->multi_const_reacts[i].getReaction(r).getForwardRate()*forward - dat->multi_const_reacts[i].getReaction(r).getReverseRate()*reverse);
			}
			else if (stoic < 0)
			{
				jac += abs(stoic)*(dat->multi_const_reacts[i].getReaction(r).getReverseRate()*reverse - dat->multi_const_reacts[i].getReaction(r).getForwardRate()*forward);
			}
			else
				jac += 0.0;

		}
		//If it is not related, then return 0
		catch (std::out_of_range)
		{
			jac += 0.0;
		}
	}

	return jac;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: MultiConstReaction
 */

/*
 *								Start: InfiniteBath
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
InfiniteBath::InfiniteBath()
{
	Value = 0;
	main_index = 0;
}

//Default destructor
InfiniteBath::~InfiniteBath()
{
	weights.clear();
}

//Initialize solver
void InfiniteBath::InitializeSolver(Dove &Solver)
{
	this->SolverInfo = &Solver;
}

//Set index
void InfiniteBath::SetIndex(int index)
{
	this->main_index = index;
}

//Set value
void InfiniteBath::SetValue(double val)
{
	this->Value = val;
}

//Insert weight
void InfiniteBath::InsertWeight(int i, double w)
{
	if (w < 0)
		w = 0;
	this->weights[i] = w;
}

//Return value
double InfiniteBath::getValue()
{
	return this->Value;
}

//Return index
int InfiniteBath::getIndex()
{
	return this->main_index;
}

//Return reference to map
std::map<int, double> & InfiniteBath::getWeightMap()
{
	return this->weights;
}

/// Rate function for the InfiniteBath Object
double rate_func_InfiniteBath(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, double>::iterator it;

	//Loop over all mapped values
	for (it = dat->inf_bath[i].getWeightMap().begin(); it != dat->inf_bath[i].getWeightMap().end(); it++)
	{
		res += u(it->first,0)*it->second;
	}
	res = dat->inf_bath[i].getValue() - res;

	return res;
}

/// Jacobi function for the InfiniteBath Object
double jacobi_func_InfiniteBath(int i, int j, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double jac = 0;
	CROW_DATA *dat = (CROW_DATA *) data;
	std::map<int, double>::iterator it;

	//Assume jacobian is only registered for species involved
	try
	{
		//Loop over all mapped values
		for (it = dat->inf_bath[i].getWeightMap().begin(); it != dat->inf_bath[i].getWeightMap().end(); it++)
		{
			if (it->first == j)
				jac += it->second;
			else
				jac += 0;
		}
	}
	catch (std::out_of_range)
	{
		jac = 0.0;
	}

	return jac;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: InfiniteBath
 */

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
	fprintf(dat->OutputFile, "\tSOLVER TYPE = \t");
	if (dat->SolverInfo.isLinear() == true)
		fprintf(dat->OutputFile, "Linear\n");
	else
	{
		if (dat->SolverInfo.isPreconditioned() == false || dat->SolverInfo.getLinearMethod() == FOM)
			fprintf(dat->OutputFile, "Jacobian-Free Newton-Krylov (JFNK)\n");
		else if (dat->SolverInfo.getLinearMethod() == QR)
			fprintf(dat->OutputFile, "Direct Newton Method\n");
		else
			fprintf(dat->OutputFile, "Preconditioned Jacobian-Free Newton-Krylov (PJFNK)\n");
	}
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
	fprintf(dat->OutputFile, "\tPRECONDITION METHOD = \t");
	if (dat->SolverInfo.isPreconditioned() == false || dat->SolverInfo.getLinearMethod() == FOM || dat->SolverInfo.getLinearMethod() == QR)
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
	fprintf(dat->OutputFile, "\tLINEAR SOLVER METHOD = \t");
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
	fprintf(dat->OutputFile, "\tRESTART ITERATION LEVEL =\t%i\n", dat->SolverInfo.getRestartLevel());
	fprintf(dat->OutputFile, "\tRECURSION ITERATION LEVEL =\t%i\n", dat->SolverInfo.getRecursionLevel());
	fprintf(dat->OutputFile, "\tMAXIMUM LINEAR ITERATIONS =\t%i\n", dat->SolverInfo.getMaxLinearIterations());
	fprintf(dat->OutputFile, "\tLINEAR TOLERANCE (Absolute) =\t%.6g\n", dat->SolverInfo.getLinearToleranceABS());
	fprintf(dat->OutputFile, "\tLINEAR TOLERANCE (Relative) =\t%.6g\n", dat->SolverInfo.getLinearToleranceREL());
	fprintf(dat->OutputFile, "\tMAXIMUM NONLINEAR ITERATIONS =\t%i\n", dat->SolverInfo.getMaxNonlinearIterations());
	fprintf(dat->OutputFile, "\tNONLINEAR TOLERANCE (Absolute) =\t%.6g\n", dat->SolverInfo.getNonlinearToleranceABS());
	fprintf(dat->OutputFile, "\tNONLINEAR TOLERANCE (Relative) =\t%.6g\n", dat->SolverInfo.getNonlinearToleranceREL());
	fprintf(dat->OutputFile,"\n-----------------------------------------------------------\n\n");

	//Loop over all ConstReaction Objects
	int i = 0;
	for (auto &x: dat->const_reacts)
	{
		if (i == 0)
			fprintf(dat->OutputFile, "---------------- Constant Reaction Objects ----------------\n\n");

		fprintf(dat->OutputFile, "Variable:\t%s\nEquation:\td(%s)/dt = ", dat->SolverInfo.getVariableName(x.first).c_str(), dat->SolverInfo.getVariableName(x.first).c_str());

		//Reaction Equation Information
		std::map<int, double>::iterator it;
		//Forward rate / positive reactants
		if (dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()] > 0)
		{
			//Forward part
			if (dat->const_reacts[x.first].getForwardRate() != 0)
			{
				if (abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]) > 1)
					fprintf(dat->OutputFile, "%.6g*", abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]));
				fprintf(dat->OutputFile, "kf");
				for (it = dat->const_reacts[x.first].getStoichiometryMap().begin(); it != dat->const_reacts[x.first].getStoichiometryMap().end(); it++)
				{
					if (it->second < 0)
					{
						if (abs(it->second) == 1)
							fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
						else
							fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
					}
					else
					{
						//Do Nothing
					}
				}
			}

			//Reverse part
			if (dat->const_reacts[x.first].getReverseRate() != 0)
			{
				if (abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]) > 1)
					fprintf(dat->OutputFile, " - %.6g*kr", abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]));
				else
					fprintf(dat->OutputFile, " - kr");
				for (it = dat->const_reacts[x.first].getStoichiometryMap().begin(); it != dat->const_reacts[x.first].getStoichiometryMap().end(); it++)
				{
					if (it->second > 0)
					{
						if (abs(it->second) == 1)
							fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
						else
							fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
					}
					else
					{
						//Do Nothing
					}
				}
			}
		}
		//Reverse rate / negative reactants
		else if (dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()] < 0)
		{
			//Reverse part
			if (dat->const_reacts[x.first].getReverseRate() != 0)
			{
				if (abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]) > 1)
					fprintf(dat->OutputFile, "%.6g*", abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]));
				fprintf(dat->OutputFile, "kr");
				for (it = dat->const_reacts[x.first].getStoichiometryMap().begin(); it != dat->const_reacts[x.first].getStoichiometryMap().end(); it++)
				{
					if (it->second > 0)
					{
						if (abs(it->second) == 1)
							fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
						else
							fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
					}
					else
					{
						//Do Nothing
					}
				}
			}

			//Forward part
			if (dat->const_reacts[x.first].getForwardRate() != 0)
			{
				if (abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]) > 1)
					fprintf(dat->OutputFile, " - %.6g*kf", abs(dat->const_reacts[x.first].getStoichiometryMap()[dat->const_reacts[x.first].getIndex()]));
				else
					fprintf(dat->OutputFile, " - kf");
				for (it = dat->const_reacts[x.first].getStoichiometryMap().begin(); it != dat->const_reacts[x.first].getStoichiometryMap().end(); it++)
				{
					if (it->second < 0)
					{
						if (abs(it->second) == 1)
							fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
						else
							fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
					}
					else
					{
						//Do Nothing
					}
				}
			}
		}
		else
		{
			//Do Nothing
		}
		fprintf(dat->OutputFile, "\n");

		//Parameter Information
		if (dat->const_reacts[x.first].getForwardRate() != 0)
			fprintf(dat->OutputFile, "Parameter:\tkf =\t%.6g\n", dat->const_reacts[x.first].getForwardRate());
		if (dat->const_reacts[x.first].getReverseRate() != 0)
			fprintf(dat->OutputFile, "Parameter:\tkr =\t%.6g\n", dat->const_reacts[x.first].getReverseRate());

		fprintf(dat->OutputFile, "\n");
		if (i == dat->const_reacts.size()-1)
			fprintf(dat->OutputFile,"-----------------------------------------------------------\n\n");
		i++;
	}//END ConstReaction Loop

	//Loop over all MultiConstReaction Objects
	i = 0;
	for (auto &x: dat->multi_const_reacts)
	{
		if (i == 0)
			fprintf(dat->OutputFile, "------------ Multiple Constant Reaction Objects ------------\n\n");

		fprintf(dat->OutputFile, "Variable:\t%s\nEquation:\td(%s)/dt = ", dat->SolverInfo.getVariableName(x.first).c_str(), dat->SolverInfo.getVariableName(x.first).c_str());

		//Loop over all reactions in object
		for (int r=0; r<dat->multi_const_reacts[x.first].getNumReactions(); r++)
		{
			//Reaction Equation Information
			std::map<int, double>::iterator it;
			//Forward rate / positive reactants
			if (dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()] > 0)
			{
				//Forward part
				if (dat->multi_const_reacts[x.first].getReaction(r).getForwardRate() != 0)
				{
					if (r == 0)
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, "%.6g*", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]));
						fprintf(dat->OutputFile, "kf_%i",r);
					}
					else
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, "+ %.6g*kf_%i", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]),r);
						else
							fprintf(dat->OutputFile, "+ kf_%i",r);
					}
					for (it = dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().end(); it++)
					{
						if (it->second < 0)
						{
							if (abs(it->second) == 1)
								fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
							else
								fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
						}
						else
						{
							//Do Nothing
						}
					}
				}

				//Reverse part
				if (dat->multi_const_reacts[x.first].getReaction(r).getReverseRate() != 0)
				{
					if (r == 0)
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, " - %.6g*kr_%i", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]),r);
						else
							fprintf(dat->OutputFile, " - kr_%i",r);
					}
					else
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, "- %.6g*kr_%i", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]),r);
						else
							fprintf(dat->OutputFile, "- kr_%i",r);
					}
					for (it = dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().end(); it++)
					{
						if (it->second > 0)
						{
							if (abs(it->second) == 1)
								fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
							else
								fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
						}
						else
						{
							//Do Nothing
						}
					}
				}
			}
			//Reverse rate / negative reactants
			else if (dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()] < 0)
			{
				//Reverse part
				if (dat->multi_const_reacts[x.first].getReaction(r).getReverseRate() != 0)
				{
					if (r == 0)
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, "%.6g*", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]));
						fprintf(dat->OutputFile, "kr_%i",r);
					}
					else
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, "+ %.6g*kr_%i", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]),r);
						else
							fprintf(dat->OutputFile, "+ kr_%i",r);
					}
					for (it = dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().end(); it++)
					{
						if (it->second > 0)
						{
							if (abs(it->second) == 1)
								fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
							else
								fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
						}
						else
						{
							//Do Nothing
						}
					}
				}

				//Forward part
				if (dat->multi_const_reacts[x.first].getReaction(r).getForwardRate() != 0)
				{
					if (r == 0)
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, " - %.6g*kf_%i", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]),r);
						else
							fprintf(dat->OutputFile, " - kf_%i",r);
					}
					else
					{
						if (abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]) > 1)
							fprintf(dat->OutputFile, "- %.6g*kf_%i", abs(dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap()[dat->multi_const_reacts[x.first].getReaction(r).getIndex()]),r);
						else
							fprintf(dat->OutputFile, "- kf_%i",r);
					}
					for (it = dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().begin(); it != dat->multi_const_reacts[x.first].getReaction(r).getStoichiometryMap().end(); it++)
					{
						if (it->second < 0)
						{
							if (abs(it->second) == 1)
								fprintf(dat->OutputFile, "*(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
							else
								fprintf(dat->OutputFile, "*(%s)^%.6g", dat->SolverInfo.getVariableName(it->first).c_str(), abs(it->second));
						}
						else
						{
							//Do Nothing
						}
					}
				}
			}
			else
			{
				//Do Nothing
			}
			fprintf(dat->OutputFile, " ");

		}//END Reaction Loop
		fprintf(dat->OutputFile, "\n");

		//Loop over all reactions in object
		for (int r=0; r<dat->multi_const_reacts[x.first].getNumReactions(); r++)
		{
			//Parameter Information
			if (dat->multi_const_reacts[x.first].getReaction(r).getForwardRate() != 0)
				fprintf(dat->OutputFile, "Parameter:\tkf_%i =\t%.6g\n",r,dat->multi_const_reacts[x.first].getReaction(r).getForwardRate());
			if (dat->multi_const_reacts[x.first].getReaction(r).getReverseRate() != 0)
				fprintf(dat->OutputFile, "Parameter:\tkr_%i =\t%.6g\n",r,dat->multi_const_reacts[x.first].getReaction(r).getReverseRate());
		}

		fprintf(dat->OutputFile, "\n");
		if (i == dat->multi_const_reacts.size()-1)
			fprintf(dat->OutputFile,"-----------------------------------------------------------\n\n");
		i++;
	}//END MultiConstReaction Loop

	//Loop over all InfiniteBath Objects
	i = 0;
	for (auto &x: dat->inf_bath)
	{
		if (i == 0)
			fprintf(dat->OutputFile, "------------------ Infinite Bath Objects ------------------\n\n");

		fprintf(dat->OutputFile, "Variable:\t%s\nEquation:\t%.6g = ", dat->SolverInfo.getVariableName(x.first).c_str(),x.second.getValue());

		//Infinite Bath Equation Information
		std::map<int, double>::iterator it;
		for (it = dat->inf_bath[x.first].getWeightMap().begin(); it != dat->inf_bath[x.first].getWeightMap().end(); it++)
		{
			if (it == dat->inf_bath[x.first].getWeightMap().begin())
			{
				if (it->second > 1)
					fprintf(dat->OutputFile, "%.6g*(%s)", it->second, dat->SolverInfo.getVariableName(it->first).c_str());
				else
					fprintf(dat->OutputFile, "(%s)", dat->SolverInfo.getVariableName(it->first).c_str());
			}
			else
			{
				if (it->second > 1)
					fprintf(dat->OutputFile, " + %.6g*(%s)", it->second, dat->SolverInfo.getVariableName(it->first).c_str());
				else
					fprintf(dat->OutputFile, " + (%s)", dat->SolverInfo.getVariableName(it->first).c_str());
			}
		}

		fprintf(dat->OutputFile, "\n\n");
		if (i == dat->inf_bath.size()-1)
			fprintf(dat->OutputFile,"-----------------------------------------------------------\n\n");
		i++;
	}//END InfiniteBath Loop

}

///Function to validate Function type
func_type function_choice(std::string &choice)
{
	func_type type = INVALID;

	std::string copy = choice;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);

	if (copy == "constreaction")
		type = CONSTREACTION;
	else if (copy == "multiconstreaction")
		type = MULTICONSTREACTION;
	else if (copy == "infinitebath")
		type = INFINITEBATH;
	else
		type = INVALID;

	return type;
}

//Read input file
int read_crow_input(CROW_DATA *dat)
{
	int success = 0;

	success = read_crow_system(dat);
	if (success != 0) {mError(read_error); return -1;}

	success = read_crow_functions(dat);
	if (success != 0) {mError(read_error); return -1;}

	return success;
}

//Function to intialize system information
int read_crow_system(CROW_DATA *dat)
{
	int success = 0;

	//Find all required info from System Document
	try
	{
		dat->SolverInfo.set_numfunc(dat->yaml_object.getYamlWrapper()("System")["num_var"].getInt());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}
	int size = 0;
	try
	{
		size = dat->yaml_object.getYamlWrapper()("System")("var_names").getDataMap().size();
		if (size != dat->SolverInfo.getNumFunc())
		{
			mError(missing_information);
			return -1;
		}
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}
	try
	{
		for (auto &x: dat->yaml_object.getYamlWrapper()("System")("var_names").getDataMap().getMap())
		{
			int index = atoi(x.first.c_str());
			dat->SolverInfo.set_variableName(index, x.second.getString());
		}
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	//Solver Options
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("solve_opts")["solver_type"].getString();
		dat->SolverInfo.set_LinearStatus(solver_choice(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_LinearStatus(false);
	}
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("solve_opts")["line_search"].getString();
		dat->SolverInfo.set_LineSearchMethod(linesearch_choice(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_LineSearchMethod(NO_LS);
	}
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("solve_opts")["linear_solver"].getString();
		dat->SolverInfo.set_LinearMethod(linearsolver_choice(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_LinearMethod(QR);
	}
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("solve_opts")["preconditioning"].getString();
		dat->SolverInfo.set_Preconditioning(use_preconditioning(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_Preconditioning(false);
	}
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("solve_opts")["preconditioning"].getString();
		dat->SolverInfo.set_preconditioner(preconditioner_choice(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_preconditioner(JACOBI);
	}
	try
	{
		dat->SolverInfo.set_RecursionLevel(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["recursion"].getInt());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_RecursionLevel(1);
	}
	try
	{
		dat->SolverInfo.set_RestartLimit(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["restart"].getInt());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_RestartLimit(20);
	}
	try
	{
		dat->SolverInfo.set_MaxLinearIterations(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["max_lin_it"].getInt());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_MaxLinearIterations(100);
	}
	try
	{
		dat->SolverInfo.set_MaxNonLinearIterations(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["max_nl_it"].getInt());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_MaxNonLinearIterations(10);
	}
	try
	{
		dat->SolverInfo.set_LinearRelTol(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["lin_rel_tol"].getDouble());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_LinearRelTol(0.01);
	}
	try
	{
		dat->SolverInfo.set_LinearAbsTol(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["lin_abs_tol"].getDouble());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_LinearAbsTol(0.0001);
	}
	try
	{
		dat->SolverInfo.set_NonlinearRelTol(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["nl_rel_tol"].getDouble());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_NonlinearRelTol(1e-6);
	}
	try
	{
		dat->SolverInfo.set_NonlinearAbsTol(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["nl_abs_tol"].getDouble());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_NonlinearAbsTol(0.0001);
	}
	try
	{
		dat->SolverInfo.set_fileoutput(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["file_output"].getBool());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_fileoutput(true);
	}
	try
	{
		dat->SolverInfo.set_output(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["console_output"].getBool());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_output(true);
	}
	try
	{
		dat->SolverInfo.set_NonlinearOutput(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["nl_output"].getBool());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_NonlinearOutput(true);
	}
	try
	{
		dat->SolverInfo.set_LinearOutput(dat->yaml_object.getYamlWrapper()("System")("solve_opts")["lin_output"].getBool());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_LinearOutput(false);
	}

	//Runtime Options
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("run_time")["timestepper"].getString();
		dat->SolverInfo.set_timestepper(timestepper_choice(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_timestepper(CONSTANT);
	}
	try
	{
		std::string choice = dat->yaml_object.getYamlWrapper()("System")("run_time")["integration"].getString();
		dat->SolverInfo.set_integrationtype(integration_choice(choice));
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_integrationtype(BE);
	}
	try
	{
		dat->SolverInfo.set_timestep(dat->yaml_object.getYamlWrapper()("System")("run_time")["dt"].getDouble());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}
	try
	{
		dat->SolverInfo.set_endtime(dat->yaml_object.getYamlWrapper()("System")("run_time")["end_time"].getDouble());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}
	try
	{
		dat->SolverInfo.set_timestepmax(dat->yaml_object.getYamlWrapper()("System")("run_time")["dtmax"].getDouble());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_timestepmax(100.0);
	}
	try
	{
		dat->SolverInfo.set_timestepmin(dat->yaml_object.getYamlWrapper()("System")("run_time")["dtmin"].getDouble());
	}
	catch (std::out_of_range)
	{
		dat->SolverInfo.set_timestepmin(0.0);
	}

	return success;
}

///Function to read the header files for each variable
int read_crow_functions(CROW_DATA *dat)
{
	int success = 0;

	//Check the number of headers
	int headers = (int)dat->yaml_object.getYamlWrapper().getDocMap().size();
	if (headers != dat->SolverInfo.getNumFunc()+1)
	{
		mError(missing_information);
		return -1;
	}

	//Create space for all functions first
	int valid_names = 0;
	try
	{
		for (auto &x: dat->yaml_object.getYamlWrapper().getDocMap())
		{
			if (dat->SolverInfo.isValidName(x.first))
			{
				std::string choice = dat->yaml_object.getYamlWrapper()(x.first)["func_type"].getString();
				valid_names++;
			}
		}

	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}
	if (valid_names != dat->SolverInfo.getNumFunc())
	{
		mError(missing_information);
		return -1;
	}

	//Loop through all document headers and call the appropriate read function
	int index = 0;
	try
	{
		for (auto &x: dat->yaml_object.getYamlWrapper().getDocMap())
		{
			if (dat->SolverInfo.isValidName(x.first))
			{
				index = dat->SolverInfo.getVariableIndex(x.first);
				std::string choice = dat->yaml_object.getYamlWrapper()(x.first)["func_type"].getString();
				func_type type = function_choice(choice);

				switch (type)
				{
					case CONSTREACTION:
						success = read_crow_ConstReaction(index, dat->const_reacts, dat->yaml_object.getYamlWrapper().getDocument(x.first), dat->SolverInfo);
						if (success == -1) {mError(read_error);return -1;}
						break;

					case MULTICONSTREACTION:
						success = read_crow_MultiConstReaction(index, dat->multi_const_reacts, dat->yaml_object.getYamlWrapper().getDocument(x.first), dat->SolverInfo);
						if (success == -1) {mError(read_error);return -1;}
						break;

					case INFINITEBATH:
						success = read_crow_InfiniteBath(index, dat->inf_bath, dat->yaml_object.getYamlWrapper().getDocument(x.first), dat->SolverInfo);
						if (success == -1) {mError(read_error);return -1;}
						break;

					default:
						mError(invalid_type);
						return -1;
						break;
				}
			}
		}

	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	return success;
}

//Function to initialize ConstReaction information
int read_crow_ConstReaction(int index, std::unordered_map<int, ConstReaction> &map, Document &info, Dove &solver)
{
	int success = 0;

	//try to set initial condition
	try
	{
		solver.set_initialcondition(index, info["initial_cond"].getDouble());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	map[index].InitializeSolver(solver);
	map[index].SetIndex(index);

	//Register appropriate functions
	solver.registerCoeff(index, default_coeff);
	solver.registerFunction(index, rate_func_ConstReaction);

	int num_rate = 0;
	//Read forward and reverse rate info (at least one is required)
	try
	{
		map[index].SetForwardRate(info["forward_rate"].getDouble());
		num_rate++;
	}
	catch (std::out_of_range)
	{
		map[index].SetForwardRate(0);
	}
	try
	{
		map[index].SetReverseRate(info["reverse_rate"].getDouble());
		num_rate++;
	}
	catch (std::out_of_range)
	{
		map[index].SetReverseRate(0);
	}
	if (num_rate == 0)
	{
		mError(missing_information);
		return -1;
	}

	//Check the size of the map
	try
	{
		if (info("stoichiometry").getDataMap().size() < 1)
		{
			mError(missing_information);
			return -1;
		}
	}
	catch (std::out_of_range)
	{
		mError(read_error);
		return -1;
	}

	//Loop through all stoichiometry
	for (auto &x: info("stoichiometry").getDataMap())
	{
		if (solver.isValidName(x.first) == false)
		{
			mError(read_error);
			return -1;
		}

		//Grab stoichiometry
		int species = solver.getVariableIndex(x.first);
		try
		{
			map[index].InsertStoichiometry(species, x.second.getDouble());
		}
		catch (std::out_of_range)
		{
			mError(read_error);
			return -1;
		}

		//Register appropriate jacobian functions
		solver.registerJacobi(index, species, jacobi_func_ConstReaction);
	}

	return success;
}

///Function to initialize MultiConstReaction information
int read_crow_MultiConstReaction(int index, std::unordered_map<int, MultiConstReaction> &map, Document &info, Dove &solver)
{
	int success = 0;

	//Must determine number of reactions for object FIRST!
	try
	{
		map[index].SetNumberReactions(info["num_reactions"].getInt());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	//try to set initial condition
	try
	{
		solver.set_initialcondition(index, info["initial_cond"].getDouble());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	map[index].InitializeSolver(solver);
	map[index].SetIndex(index);

	//Register appropriate functions
	solver.registerCoeff(index, default_coeff);
	solver.registerFunction(index, rate_func_MultiConstReaction);

	//Check number of headers to ensure it is correct
	int headers = (int)info.getHeadMap().size();
	if (headers != map[index].getNumReactions())
	{
		mError(missing_information);
		return -1;
	}

	//Loop through all headers to setup reactions (NOTE: header names don't matter)
	int react = 0;
	for (auto &y: info.getHeadMap())
	{
		int num_rate = 0;
		//Read forward and reverse rate info (at least one is required)
		try
		{
			map[index].SetForwardRate(react, info(y.first)["forward_rate"].getDouble());
			num_rate++;
		}
		catch (std::out_of_range)
		{
			map[index].SetForwardRate(react,0);
		}
		try
		{
			map[index].SetReverseRate(react, info(y.first)["reverse_rate"].getDouble());
			num_rate++;
		}
		catch (std::out_of_range)
		{
			map[index].SetReverseRate(react,0);
		}
		if (num_rate == 0)
		{
			mError(missing_information);
			return -1;
		}

		//Check the size of the map
		try
		{
			if (info(y.first)("stoichiometry").getMap().size() < 1)
			{
				mError(missing_information);
				return -1;
			}
		}
		catch (std::out_of_range)
		{
			mError(read_error);
			return -1;
		}

		//Loop through all stoichiometry
		for (auto &x: info(y.first)("stoichiometry").getMap())
		{
			if (solver.isValidName(x.first) == false)
			{
				mError(read_error);
				return -1;
			}

			//Grab stoichiometry
			int species = solver.getVariableIndex(x.first);
			try
			{
				map[index].InsertStoichiometry(react, species, x.second.getDouble());
			}
			catch (std::out_of_range)
			{
				mError(read_error);
				return -1;
			}

			//Register appropriate jacobian functions
			solver.registerJacobi(index, species, jacobi_func_ConstReaction);
		}
		react++;
	}

	return success;
}

///Function to initialize InfiniteBath information
int read_crow_InfiniteBath(int index, std::unordered_map<int, InfiniteBath> &map, Document &info, Dove &solver)
{
	int success = 0;

	map[index].InitializeSolver(solver);
	map[index].SetIndex(index);

	//Register appropriate functions
	solver.registerCoeff(index, default_coeff);
	solver.registerFunction(index, rate_func_InfiniteBath);

	//Try to set constant value
	try
	{
		map[index].SetValue(info["const_value"].getDouble());
		solver.set_initialcondition(index, info["const_value"].getDouble());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	//Try to set initial condition
	try
	{
		solver.set_initialcondition(index, info["initial_cond"].getDouble());
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	//Check the size of the map
	try
	{
		if (info("weights").getDataMap().size() < 1)
		{
			mError(missing_information);
			return -1;
		}
	}
	catch (std::out_of_range)
	{
		mError(read_error);
		return -1;
	}

	//Loop through all weights
	bool found_index = false;
	for (auto &x: info("weights").getDataMap())
	{
		if (solver.isValidName(x.first) == false)
		{
			mError(read_error);
			return -1;
		}

		int species = solver.getVariableIndex(x.first);
		if (index == species)
			found_index = true;

		map[index].InsertWeight(species, x.second.getDouble());
		solver.registerJacobi(index, species, jacobi_func_InfiniteBath);
	}

	//Check to ensure non-singularity
	if (found_index == false)
	{
		mError(missing_information);
		mError(singular_matrix);
		return -1;
	}

	solver.set_variableSteadyState(index);

	return success;
}

//Execute CROW
int CROW_SCENARIO(const char *yaml_input)
{
	int success = 0;
	double time;

	// ---------------------------- Initializations ---------------------------------
	time = clock();
	CROW_DATA crow;
	FILE *file;
	file = fopen("output/CROW_Results.txt", "w+");
	if (file == nullptr)
	{
		success = system("mkdir output");
		file = fopen("output/CROW_Results.txt", "w+");
	}
	crow.SolverInfo.set_outputfile(file);
	crow.OutputFile = file;
	crow.SolverInfo.set_userdata((void*)&crow);

	//Read input file
	success = crow.yaml_object.executeYamlRead(yaml_input);
	if (success != 0) {mError(file_dne); return -1;}
	success = read_crow_input(&crow);
	if (success != 0) {mError(read_error); return -1;}

	//Execute solver functions
	print2file_crow_header(&crow);
	crow.SolverInfo.solve_all();

	//Exit Messages
	time = clock() - time;
	std::cout << "\nCROW Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";

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
		success = system("mkdir output");
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
	test01.SolverInfo.set_LinearMethod(QR);
	test01.SolverInfo.set_preconditioner(SGS);
	test01.SolverInfo.set_Preconditioning(true);

	//---Function Info (Specific to ConstReaction)---

	//Use emplace_back for each instance of a variable (C++11)
	/*
	//test01.const_reacts.emplace();
	test01.const_reacts[0].InitializeSolver(test01.SolverInfo);
	test01.const_reacts[0].SetIndex(0); //Index are the same!!!
	test01.const_reacts[0].SetForwardRate(1);
	test01.const_reacts[0].SetReverseRate(1);
	test01.const_reacts[0].InsertStoichiometry(0, -1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[0].InsertStoichiometry(1, -2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[0].InsertStoichiometry(2, 2); //Species 0, stoic 1 (+) for products
	test01.SolverInfo.registerCoeff(0, default_coeff); //Variable index
	test01.SolverInfo.registerFunction(0, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(0, 0, jacobi_func_ConstReaction); //automate
	test01.SolverInfo.registerJacobi(0, 1, jacobi_func_ConstReaction); //automate
	 */

	test01.multi_const_reacts[0].SetNumberReactions(1);//MUST BE CALLED FIRST!!!
	test01.multi_const_reacts[0].InitializeSolver(test01.SolverInfo);
	test01.multi_const_reacts[0].SetIndex(0);

	test01.multi_const_reacts[0].SetForwardRate(0, 1);
	test01.multi_const_reacts[0].SetReverseRate(0, 1);
	test01.multi_const_reacts[0].InsertStoichiometry(0, 0, -1); //Species 0, stoic 1 (+) for products
	test01.multi_const_reacts[0].InsertStoichiometry(0, 1, -2); //Species 0, stoic 1 (+) for products
	test01.multi_const_reacts[0].InsertStoichiometry(0, 2, 2); //Species 0, stoic 1 (+) for products

	test01.SolverInfo.registerCoeff(0, default_coeff);
	test01.SolverInfo.registerFunction(0, rate_func_MultiConstReaction);
	test01.SolverInfo.registerJacobi(0, 0, jacobi_func_MultiConstReaction);
	test01.SolverInfo.registerJacobi(0, 1, jacobi_func_MultiConstReaction);

	//test01.const_reacts.emplace();
	test01.const_reacts[1].InitializeSolver(test01.SolverInfo);
	test01.const_reacts[1].SetIndex(1); //Index are the same!!!
	test01.const_reacts[1].SetForwardRate(1);
	test01.const_reacts[1].SetReverseRate(1);
	test01.const_reacts[1].InsertStoichiometry(0, -1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[1].InsertStoichiometry(1, -2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[1].InsertStoichiometry(2, 2); //Species 0, stoic 1 (+) for products
	test01.SolverInfo.registerCoeff(1, default_coeff);
	test01.SolverInfo.registerFunction(1, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(1, 0, jacobi_func_ConstReaction);
	test01.SolverInfo.registerJacobi(1, 1, jacobi_func_ConstReaction);

	/*
	test01.const_reacts.emplace();
	test01.const_reacts[2].InitializeSolver(test01.SolverInfo);
	test01.const_reacts[2].SetIndex(2); //Index are the same!!!
	test01.const_reacts[2].SetForwardRate(1);
	test01.const_reacts[2].SetReverseRate(1);
	test01.const_reacts[2].InsertStoichiometry(0, -1); //Species 0, stoic 1 (+) for products
	test01.const_reacts[2].InsertStoichiometry(1, -2); //Species 0, stoic 1 (+) for products
	test01.const_reacts[2].InsertStoichiometry(2, 2); //Species 0, stoic 1 (+) for products
	test01.SolverInfo.registerCoeff(2, default_coeff);
	test01.SolverInfo.registerFunction(2, rate_func_ConstReaction);
	test01.SolverInfo.registerJacobi(2, 0, jacobi_func_ConstReaction);
	test01.SolverInfo.registerJacobi(2, 1, jacobi_func_ConstReaction);
	*/


	//test01.multi_const_reacts.emplace();
	test01.multi_const_reacts[2].SetNumberReactions(2);//MUST BE CALLED FIRST!!!
	test01.multi_const_reacts[2].InitializeSolver(test01.SolverInfo);
	test01.multi_const_reacts[2].SetIndex(2);

	test01.multi_const_reacts[2].SetForwardRate(0, 0);
	test01.multi_const_reacts[2].SetReverseRate(0, 1);
	test01.multi_const_reacts[2].InsertStoichiometry(0, 0, -1); //Species 0, stoic 1 (+) for products
	test01.multi_const_reacts[2].InsertStoichiometry(0, 1, -2); //Species 0, stoic 1 (+) for products
	test01.multi_const_reacts[2].InsertStoichiometry(0, 2, 2); //Species 0, stoic 1 (+) for products

	test01.multi_const_reacts[2].SetForwardRate(1, 1);
	test01.multi_const_reacts[2].SetReverseRate(1, 0);
	test01.multi_const_reacts[2].InsertStoichiometry(1, 0, -1); //Species 0, stoic 1 (+) for products
	test01.multi_const_reacts[2].InsertStoichiometry(1, 1, -2); //Species 0, stoic 1 (+) for products
	test01.multi_const_reacts[2].InsertStoichiometry(1, 2, 2); //Species 0, stoic 1 (+) for products


	test01.SolverInfo.registerCoeff(2, default_coeff);
	test01.SolverInfo.registerFunction(2, rate_func_MultiConstReaction);
	test01.SolverInfo.registerJacobi(2, 0, jacobi_func_MultiConstReaction);
	test01.SolverInfo.registerJacobi(2, 1, jacobi_func_MultiConstReaction);

	//test01.SolverInfo.set_variableSteadyStateAll();


	//---Call solver---
	print2file_crow_header(&test01);
	test01.SolverInfo.solve_all();

	//---Exit Messages and cleanup---
	time = clock() - time;
	std::cout << "\nCROW Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";

	return success;
}
