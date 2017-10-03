/*!
 *  \file dove.h
 *	\brief Dynamic Ode solver with Various Established methods
 *  \author Austin Ladshaw
 *	\date 09/25/2017
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2017, all
 *             rights reserved.
 */

#include "dove.h"

/*
 *								Start: Dove Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
Dove::Dove()
{
	dt = 0.0;
	dt_old = 0.0;
	time_end = 0.0;
	time = 0.0;
	time_old = 0.0;
	time_older = 0.0;
	dtmin = sqrt(DBL_EPSILON);
	dtmax = 100.0;
	tolerance = 1e-6;
	int_type = IMPLICIT;
	int_sub = BE;
	timestepper = CONSTANT;
	Output = nullptr;
	num_func = 1;
	Converged = true;
	user_data = NULL;
	residual = residual_BE;
	precon = NULL;
	newton_dat.LineSearch = true;
	newton_dat.NL_Output = false;
	DoveOutput = false;
}

//Default destructor
Dove::~Dove()
{
	this->user_jacobi.clear();
	if (this->Output != nullptr)
		fclose(Output);
}

//Set number of functions
void Dove::set_numfunc(int i)
{
	if (i < 1)
	{
		mError(invalid_size);
		this->un.set_size(1, 1);
		this->unp1.set_size(1, 1);
		this->unm1.set_size(1, 1);
		this->user_func.set_size(1, 1);
		this->user_coeff.set_size(1, 1);
		this->user_jacobi.resize(1);
		this->num_func = 1;
	}
	else
	{
		this->un.set_size(i, 1);
		this->unp1.set_size(i, 1);
		this->unm1.set_size(i, 1);
		this->user_func.set_size(i, 1);
		this->user_coeff.set_size(i, 1);
		this->user_jacobi.resize(i);
		this->num_func = i;
	}
	this->set_defaultCoeffs();
	this->set_defaultJacobis();
}

//Set time step to value
void Dove::set_timestep(double d)
{
	if (d <= sqrt(DBL_EPSILON))
		this->dt = sqrt(DBL_EPSILON);
	else
		this->dt = d;
}

//Set min dt
void Dove::set_timestepmin(double dmin)
{
	if (dmin <= sqrt(DBL_EPSILON))
		this->dtmin = sqrt(DBL_EPSILON);
	else
		this->dtmin = dmin;
}

//Set max dt
void Dove::set_timestepmax(double dmax)
{
	if (dmax <= sqrt(DBL_EPSILON))
		this->dtmax = sqrt(DBL_EPSILON);
	else
		this->dtmax = dmax;
}

//Set the end time
void Dove::set_endtime(double e)
{
	if (e <= sqrt(DBL_EPSILON))
		this->time_end = sqrt(DBL_EPSILON);
	else
		this->time_end = e;
}

//Set the type of integration scheme to use
void Dove::set_integrationtype(integrate_subtype type)
{
	this->int_sub = type;
	switch (type)
	{
		case BE:
			this->int_type = IMPLICIT;
			this->residual = residual_BE;
			break;
			
		case FE:
			this->int_type = EXPLICIT;
			this->residual = nullptr;
			break;
			
		case CN:
			this->int_type = IMPLICIT;
			this->residual = residual_CN;
			break;
			
		case BDF2:
			this->int_type = IMPLICIT;
			this->residual = residual_BDF2;
			break;
			
		case RK4:
			this->int_type = EXPLICIT;
			this->residual = nullptr;
			break;
			
		case RKF:
			this->int_type = EXPLICIT;
			this->residual = nullptr;
			break;
			
		default:
			this->int_type = IMPLICIT;
			this->residual = residual_BE;
			break;
	}
}

//Set the type of time stepper
void Dove::set_timestepper(timestep_type type)
{
	this->timestepper = type;
}

//Set the output file
void Dove::set_outputfile(FILE *file)
{
	this->Output = file;
}

//Set user data
void Dove::set_userdata(const void *data)
{
	this->user_data = data;
}

//Set the initial condition
void Dove::set_initialcondition(int i, double ic)
{
	this->un.edit(i, 0, ic);
	this->unp1.edit(i, 0, ic);
	this->unm1.edit(i, 0, ic);
}

//Set output conditions for Dove
void Dove::set_output(bool choice)
{
	this->DoveOutput = choice;
}

//Set the tolerance
void Dove::set_tolerance(double tol)
{
	if (tol < MIN_TOL)
		tol = MIN_TOL;
	this->tolerance = tol;
	this->newton_dat.nl_tol_abs = tol;
}

//Set nl_abs_tol
void Dove::set_NonlinearAbsTol(double tol)
{
	this->newton_dat.nl_tol_abs = tol;
}

//Set nl_rel_tol
void Dove::set_NonlinearRelTol(double tol)
{
	this->newton_dat.nl_tol_rel	= tol;
}

//Set l_abs_tol
void Dove::set_LinearAbsTol(double tol)
{
	this->newton_dat.lin_tol_abs = tol;
}

//Set l_rel_tol
void Dove::set_LinearRelTol(double tol)
{
	this->newton_dat.lin_tol_rel = tol;
}

//Set the default coeffs
void Dove::set_defaultCoeffs()
{
	for (int i=0; i<this->num_func; i++)
		this->registerCoeff(i, default_coeff);
}

//Set the default jacobi
void Dove::set_defaultJacobis()
{
	for (int i=0; i<this->num_func; i++)
		this->registerJacobi(i, i, default_jacobi);
}

//Set NL output
void Dove::set_NonlinearOutput(bool choice)
{
	this->newton_dat.NL_Output = choice;
}

//Set L output
void Dove::set_LinearOutput(bool choice)
{
	this->newton_dat.L_Output = choice;
}

//Set linear method
void Dove::set_LinearMethod(krylov_method choice)
{
	this->newton_dat.linear_solver = choice;
}

//Set the line search method
void Dove::set_LineSearchMethod(linesearch_type choice)
{
	switch (choice)
	{
		case BT:
			this->newton_dat.LineSearch = true;
			this->newton_dat.Bounce = false;
			break;
			
		case ABT:
			this->newton_dat.LineSearch = true;
			this->newton_dat.Bounce = true;
			break;
			
		case NO_LS:
			this->newton_dat.LineSearch = false;
			this->newton_dat.Bounce = false;
			break;
			
		default:
			this->newton_dat.LineSearch = false;
			this->newton_dat.Bounce = false;
			break;
	}
}

//Register user function
void Dove::registerFunction(int i, double (*func) (int i, const Matrix<double> &u, double t, const void *data) )
{
	if ((*func) == NULL)
	{
		mError(nullptr_func);
		this->user_func.edit(i, 0, default_func);
	}
	else
		this->user_func.edit(i, 0, func);
}

//Register time coeff functions
void Dove::registerCoeff(int i, double (*coeff) (int i, const Matrix<double> &u, double t, const void *data) )
{
	if ((*coeff) == NULL)
		this->user_coeff.edit(i, 0, default_coeff);
	else
		this->user_coeff.edit(i, 0, coeff);
}

//Register jacobians
void Dove::registerJacobi(int i, int j, double (*jac) (int i, int j, const Matrix<double> &u, double t, const void *data) )
{
	if ((*jac) == NULL)
	{
		this->user_jacobi[i][j] = default_jacobi;
	}
	else
	{
		this->user_jacobi[i][j] = jac;
	}
}

//Print out header info to output file
void Dove::print_header()
{
	if (this->Output == nullptr)
		this->Output = fopen("output/DOVE_Result.txt", "w+");
	if (this->Output == nullptr)
	{
		system("mkdir output");
		this->Output = fopen("output/DOVE_Result.txt", "w+");
	}
	fprintf(this->Output,"\nIntegration type =\t");
	switch (this->int_type)
	{
		case IMPLICIT:
			fprintf(this->Output,"IMPLICIT\n");
			break;
			
		case EXPLICIT:
			fprintf(this->Output,"EXPLICIT\n");
			break;
			
		default:
			fprintf(this->Output,"IMPLICIT\n");
			break;
	}
	fprintf(this->Output,"Integration scheme =\t");
	switch (this->int_sub)
	{
		case BE:
			fprintf(this->Output,"Backward-Euler\n");
			break;
			
		case FE:
			fprintf(this->Output,"Forward-Euler\n");
			break;
			
		case CN:
			fprintf(this->Output,"Crank-Nicholson\n");
			break;
			
		case BDF2:
			fprintf(this->Output,"Backward-Differentiation-Formula-2\n");
			break;
			
		case RK4:
			fprintf(this->Output,"Runge-Kutta-4\n");
			break;
			
		case RKF:
			fprintf(this->Output,"Runge-Kutta-Fehlberg\n");
			break;
			
		default:
			break;
	}
	fprintf(this->Output,"Time");
	for (int i=0; i<this->num_func; i++)
		fprintf(this->Output,"\tu[%i]",i);
	fprintf(this->Output,"\n");
}

//Print new result
void Dove::print_newresult()
{
	fprintf(this->Output,"%.6g",this->time);
	for (int i=0; i<this->num_func; i++)
		fprintf(this->Output,"\t%.6g",this->unp1(i,0));
	fprintf(this->Output,"\n");
}

//Print result
void Dove::print_result()
{
	fprintf(this->Output,"%.6g",this->time);
	for (int i=0; i<this->num_func; i++)
		fprintf(this->Output,"\t%.6g",this->un(i,0));
	fprintf(this->Output,"\n");
}

//Return reference to un
Matrix<double>& Dove::getCurrentU()
{
	return this->un;
}

//Return reference to unm1
Matrix<double>& Dove::getOldU()
{
	return this->unm1;
}

//Return reference to unp1
Matrix<double>& Dove::getNewU()
{
	return this->unp1;
}

//Return pointer to user data
const void* Dove::getUserData()
{
	return this->user_data;
}

//Return number of functions
int Dove::getNumFunc()
{
	return this->num_func;
}

//Return dt
double Dove::getTimeStep()
{
	return this->dt;
}

//Return dt_old
double Dove::getTimeStepOld()
{
	return this->dt_old;
}

//Return end time
double Dove::getEndTime()
{
	return this->time_end;
}

//Return time
double Dove::getCurrentTime()
{
	return this->time;
}

//Return time old
double Dove::getOldTime()
{
	return this->time_old;
}

//Return older time
double Dove::getOlderTime()
{
	return this->time_older;
}

//Return dtmin
double Dove::getMinTimeStep()
{
	return this->dtmin;
}

//Return dtmax
double Dove::getMaxTimeStep()
{
	return this->dtmax;
}

//Return bool for convergence
bool Dove::hasConverged()
{
	return this->Converged;
}

//Return non-linear res
double Dove::getNonlinearResidual()
{
	switch (this->int_sub)
	{
		case BE:
			return this->newton_dat.nl_res;
			break;
			
		case FE:
			return this->tolerance*pow((1.0/0.84),-4.0);
			break;
			
		case CN:
			return this->newton_dat.nl_res;
			break;
			
		case BDF2:
			return this->newton_dat.nl_res;
			break;
			
		case RK4:
			return this->tolerance*pow((1.0/0.84),-4.0);
			break;
			
		default:
			return this->newton_dat.nl_res;
			break;
	}
}

//Return non-linear rel res
double Dove::getNonlinearRelativeRes()
{
	return this->newton_dat.nl_relres;
}

//Compute next time step
double Dove::ComputeTimeStep()
{
	double step = 0.0;
	if (this->time == 0.0)
		return this->dt;
	if (this->Converged == true)
	{
		switch (this->timestepper)
		{
			case CONSTANT:
				step = this->dt;
				break;
				
			case ADAPTIVE:
				step = 1.5 * this->dt;
				if (step >= this->dtmax)
					step = this->dtmax;
				break;
				
			case FEHLBERG:
				step = 0.84*pow((this->tolerance/this->getNonlinearResidual()),0.25)*this->dt;
				if (step >= this->dtmax)
					step = this->dtmax;
				if (step <= this->dtmin)
					step = this->dtmin;
				break;
				
			default:
				step = this->dt;
				break;
		}
	}
	else
	{
		switch (this->timestepper)
		{
			case CONSTANT:
				this->dtmax = this->dt;
				this->timestepper = ADAPTIVE;
				step = 0.5 * this->dt;
				if (step <= this->dtmin)
					step = this->dtmin;
				break;
				
			case ADAPTIVE:
				step = 0.5 * this->dt;
				if (step <= this->dtmin)
					step = this->dtmin;
				break;
				
			case FEHLBERG:
				step = 0.84*pow((this->tolerance/this->getNonlinearResidual()),0.25)*this->dt;
				if (step >= this->dtmax)
					step = this->dtmax;
				if (step <= this->dtmin)
					step = this->dtmin;
				break;
				
			default:
				step = this->dt;
				break;
		}
	}
	
	if (this->time_end < step + this->time_old)
		step = this->time_end - this->time_old;
	
	return step;
}

//Eval user function i
double Dove::Eval_Func(int i, const Matrix<double>& u, double t)
{
	return this->user_func(i,0)(i,u,t,this->user_data);
}

//Eval user time coefficient function i
double Dove::Eval_Coeff(int i, const Matrix<double>& u, double t)
{
	return this->user_coeff(i,0)(i,u,t,this->user_data);
}

//Eval user time coefficient function i
double Dove::Eval_Jacobi(int i, int j, const Matrix<double>& u, double t)
{
	std::map<int, double (*) (int i, int j, const Matrix<double> &u, double t, const void *data)>::iterator it = this->user_jacobi[i].find(j);
	if (it == this->user_jacobi[i].end())
	{
		return default_jacobi(i,j,u,t,this->user_data);
	}
	else
	{
		return it->second(i,j,u,t,this->user_data);
	}
}

//Function to solve a single timestep
int Dove::solve_timestep()
{
	int success = 0;
	if (this->int_type == IMPLICIT)
	{
		success = pjfnk(this->residual, this->precon, this->unp1, &this->newton_dat, this, this);
		this->Converged = this->newton_dat.Converged;
	}
	else
	{
		switch (this->int_sub)
		{
			case FE:
				success = this->solve_FE();
				break;
				
			case RK4:
				success = this->solve_RK4();
				break;
				
			case RKF:
				success = this->solve_RKF();
				break;
				
			default:
				success = this->solve_FE();
				break;
		}
	}
	
	//What to do on failure
	if (this->Converged == false)
	{
		if (this->dt > this->dtmin)
		{
			if (this->DoveOutput == true)
			{
				if (this->int_sub == RKF)
					std::cout << "Failed to converge: Residual(" << this->newton_dat.nl_res << ") > Tolerance(" << this->tolerance << ")\n";
				else
					std::cout << "Failed to converge: Residual(" << this->newton_dat.nl_res << ") > Tolerance(" << this->newton_dat.nl_tol_abs << ")\n";
			}
			if (this->DoveOutput == true)
				std::cout << "Retrying simulation with with new time step: dt_old(" << this->dt << ") --> dt_new(";
			this->dt = this->ComputeTimeStep();
			this->time = this->time_old + this->dt;
			if (this->DoveOutput == true)
				std::cout << this->dt << ") --> for time (" << this->time << ")\n\n";
			success = solve_timestep();
		}
		else
		{
			if (this->DoveOutput == true)
				std::cout << "Unable to further reduce time step. CRITICAL ERROR!!!\n";
			success = -1;
		}
		
	}
	return success;
}

//Update solution states
void Dove::update_states()
{
	this->unm1 = this->un;
	this->un = this->unp1;
	this->dt_old = this->dt;
	this->time_older = this->time_old;
	this->time_old = this->time;
}

//Update the time step
void Dove::update_timestep()
{
	this->dt = this->ComputeTimeStep();
	this->time = this->time_old + this->dt;
}

//Reset all states
void Dove::reset_all()
{
	this->time = 0.0;
	this->time_old = 0.0;
	this->time_older = 0.0;
	this->dt_old = 0.0;
	this->Converged = true;
}

//Function to solve all states and print output to file
int Dove::solve_all()
{
	int success = 0;
	this->reset_all();
	this->print_header();
	this->print_result();
	if (this->DoveOutput == true)
	{
		std::cout << "Dove Scheme: ";
		switch (this->int_sub)
		{
			case BE:
				std::cout << "Backwards-Euler method.";
				break;
				
			case FE:
				std::cout << "Forwards-Euler method.";
				break;
				
			case CN:
				std::cout << "Crank-Nicholson method. ";
				break;
				
			case BDF2:
				std::cout << "Backwards-Differentiation 2nd Order method.";
				break;
				
			case RK4:
				std::cout << "Runge-Kutta 4th Order method.";
				break;
				
			case RKF:
				std::cout << "Runge-Kutta-Fehlberg method.";
				break;
				
			default:
				std::cout << "Backwards-Euler method.";
				break;
				
		}
		std::cout << "\n------------------------------------------------------";
	}
	
	//Do-while loop
	do
	{
		this->update_timestep();
		if (this->DoveOutput == true)
		{
			std::cout << "\nSolving time (" << this->time << ") with time step (" << this->dt << "). Please wait...\n";
		}
		success = this->solve_timestep();
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		this->print_newresult();
		this->update_states();
	} while (this->time_end > (this->time+this->dtmin));
	if (this->DoveOutput == true)
		std::cout << "------------------------------------------------------\n\n";
	
	return success;
}

//Function to solve with Forward-Euler
int Dove::solve_FE()
{
	int success = 0;
	
	for (int i=0; i<this->num_func; i++)
	{
		double value = ( (this->Eval_Coeff(i, this->un, this->time_old)*this->un(i,0)) + (this->getTimeStep()*this->Eval_Func(i, this->un, this->time_old)) )/this->Eval_Coeff(i, this->un, this->time);
		this->unp1.edit(i, 0, value );
		
		if (isinf(value) || isnan(value))
		{
			this->Converged = false;
			return -1;
		}
	}
	this->Converged = true;
	return success;
}

//Function to solve with Runge-Kutta-4
int Dove::solve_RK4()
{
	int success = 0;
	Matrix<double> temp;
	temp.set_size(this->num_func, 1);
	for (int i=0; i<this->num_func; i++)
	{
		double k1,k2,k3,k4;
		temp = this->un;
		k1 = this->dt*this->Eval_Func(i, temp, this->time_old);
		temp(i,0) = this->un(i,0) + (k1/2.0);
		k2 = this->dt*this->Eval_Func(i, temp, this->time_old+(this->dt/2.0));
		temp(i,0) = this->un(i,0) + (k2/2.0);
		k3 = this->dt*this->Eval_Func(i, temp, this->time_old+(this->dt/2.0));
		temp(i,0) = this->un(i,0) + k3;
		k4 = this->dt*this->Eval_Func(i, temp, this->time);
		
		double value = ( (this->Eval_Coeff(i, this->un, this->time_old)*this->un(i,0)) + ((k1+(2.0*(k2+k3))+k4)/6.0) )/this->Eval_Coeff(i, this->un, this->time_old);
		
		this->unp1.edit(i, 0, value );
		
		if (isinf(value) || isnan(value))
		{
			this->Converged = false;
			return -1;
		}
	}
	this->Converged = true;
	return success;
}

//Function to solve with Runge-Kutta-Fehlberg
int Dove::solve_RKF()
{
	int success = 0;
	double res = 0.0;
	Matrix<double> temp;
	temp.set_size(this->num_func, 1);
	for (int i=0; i<this->num_func; i++)
	{
		double k1,k2,k3,k4,k5,k6;
		temp = this->un;
		k1 = this->dt*this->Eval_Func(i, temp, this->time_old);
		temp(i,0) = this->un(i,0) + (k1/4.0);
		k2 = this->dt*this->Eval_Func(i, temp, this->time_old+(this->dt/4.0));
		temp(i,0) = this->un(i,0) + (3.0*k1/32.0) + (9.0*k2/32.0);
		k3 = this->dt*this->Eval_Func(i, temp, this->time_old+(3.0*this->dt/8.0));
		temp(i,0) = this->un(i,0) + (1932.0*k1/2197.0) - (7200.0*k2/2197.0) + (7296.0*k3/2197.0);
		k4 = this->dt*this->Eval_Func(i, temp, this->time_old+(12.0*this->dt/13.0));
		temp(i,0) = this->un(i,0) + (439.0*k1/216.0) - (8.0*k2) + (3680.0*k3/513.0) - (845.0*k4/4104.0);
		k5 = this->dt*this->Eval_Func(i, temp, this->time);
		temp(i,0) = this->un(i,0) - (8.0*k1/27.0) + (2.0*k2) - (3544.0*k3/2565.0) + (1859.0*k4/4104.0) - (11.0*k5/40.0);
		k6 = this->dt*this->Eval_Func(i, temp, this->time_old+(this->dt/2.0));
		
		double value = ( (this->Eval_Coeff(i, this->un, this->time_old)*this->un(i,0)) + ( (25.0*k1/216.0)+(1408.0*k3/2565.0)+(2197.0*k4/4104.0)-(k5/5.0) ) )/this->Eval_Coeff(i, this->un, this->time_old);
		double value_til = ( (this->Eval_Coeff(i, this->un, this->time_old)*this->un(i,0)) + ( (16.0*k1/135.0)+(6656.0*k3/12825.0)+(28561.0*k4/56430.0)-(9.0*k5/50.0)+(2.0*k6/55.0) ) )/this->Eval_Coeff(i, this->un, this->time_old);
		
		this->unp1.edit(i, 0, value );
		if (fabs(value_til - value)/this->dt > res)
			res = fabs(value_til - value)/this->dt;
		
		if (isinf(value) || isnan(value))
		{
			this->Converged = false;
			return -1;
		}
	}
	this->newton_dat.nl_res = res;
	if (res > this->tolerance)
		this->Converged = false;
	else
		this->Converged = true;
	return success;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Dove Class Definitions
 */

//Function for implicit-BE method residual
int residual_BE(const Matrix<double> &u, Matrix<double> &Res, const void *data)
{
	int success = 0;
	Dove *dat = (Dove *) data;
	
	for (int i=0; i<dat->getNumFunc(); i++)
	{
		Res(i,0) = (dat->Eval_Coeff(i, u, dat->getCurrentTime())*u(i,0)) - (dat->Eval_Coeff(i, dat->getCurrentU(),dat->getOldTime())*dat->getCurrentU()(i,0)) - (dat->getTimeStep()*dat->Eval_Func(i, u,dat->getCurrentTime()));
	}
	
	return success;
}

//Function for implicit-CN method residual
int residual_CN(const Matrix<double> &u, Matrix<double> &Res, const void *data)
{
	int success = 0;
	Dove *dat = (Dove *) data;
	
	for (int i=0; i<dat->getNumFunc(); i++)
	{
		Res(i,0) = (dat->Eval_Coeff(i, u, dat->getCurrentTime())*u(i,0)) - (dat->Eval_Coeff(i, dat->getCurrentU(), dat->getOldTime())*dat->getCurrentU()(i,0)) - (0.5*dat->getTimeStep()*dat->Eval_Func(i, u, dat->getCurrentTime())) - (0.5*dat->getTimeStep()*dat->Eval_Func(i, dat->getCurrentU(), dat->getOldTime()));
	}
	
	return success;
}

//Function for implicit-BDF2 method residual
int residual_BDF2(const Matrix<double> &u, Matrix<double> &Res, const void *data)
{
	int success = 0;
	Dove *dat = (Dove *) data;
	
	double rn = 0.0;
	if (dat->getOldTime() > 0.0)
		rn = dat->getTimeStep()/dat->getTimeStepOld();
	
	double an, bn, cn;
	an = (1.0 + (2.0*rn)) / (1.0 + rn);
	bn = (1.0 + rn);
	cn = (rn*rn)/(1.0+rn);
	
	for (int i=0; i<dat->getNumFunc(); i++)
	{
		Res(i,0) = (an*dat->Eval_Coeff(i, u, dat->getCurrentTime())*u(i,0)) - (bn*dat->Eval_Coeff(i, dat->getCurrentU(), dat->getOldTime())*dat->getCurrentU()(i,0)) + (cn*dat->Eval_Coeff(i, dat->getOldU(), dat->getOlderTime())*dat->getOldU()(i,0)) - (dat->getTimeStep()*dat->Eval_Func(i, u, dat->getCurrentTime()));
	}
	
	return success;
}

/// Default  function
double default_func(int i, const Matrix<double> &u, double t, const void *data)
{
	return 0.0;
}

/// Default time coefficient function
double default_coeff(int i, const Matrix<double> &u, double t, const void *data)
{
	return 1.0;
}

/// Default Jacobian element function
double default_jacobi(int i, int j, const Matrix<double> &u, double t, const void *data)
{
	return 0.0;
}


// -------------------- Begin temporary testing --------------------------
double f0(int i, const Matrix<double> &x, double t, const void *res_data)
{
	return x(0,0) + 1;
}

double f1(int i, const Matrix<double> &x, double t, const void *res_data)
{
	return x(1,0) - x(0,0);
}

int test_matvec(const Matrix<double> &x, Matrix<double> &Mx, const void *data)
{
	Mx(0,0) = 5*x(0,0);
	return 0;
}

int test_res(const Matrix<double> &x, Matrix<double> &Mx, const void *data)
{
	Mx(0,0) = 5.0*x(0,0)*x(0,0) - 1.0;
	return 0;
}

double first_order_decay(int i, const Matrix<double> &u, double t, const void *data)
{
	return -u(i,0);
}

double nonlinear_first_order_decay(int i, const Matrix<double> &u, double t, const void *data)
{
	return u(0,0)*u(1,0);
}
// -------------------- End temporary testing --------------------------

//Test function
int DOVE_TESTS()
{
	int success = 0;
	
	/*  MISC Tests of Various Ideas
	Matrix<double (*) (int i, const Matrix<double> &x, double t, const void *res_data)> list_func;
	list_func.set_size(2, 1);
	list_func.edit(0, 0, f0);
	list_func.edit(1, 0, f1);
	
	Matrix<double> xtest;
	xtest.set_size(2, 1);
	xtest.edit(0, 0, 1);
	xtest.edit(1, 0, 2);
	xtest.Display("x");
	
	std::cout << "f0 = " << list_func(0,0)(0,xtest,0,NULL) << std::endl;
	std::cout << "f1 = " << list_func(1,0)(1,xtest,0,NULL) << std::endl;
	
	Dove test;
	FILE *testfile;
	testfile = fopen("output/Here_is_test.txt","w+");
	fprintf(testfile,"Here is a line of text.");
	test.set_outputfile(testfile);
	test.set_numfunc(2);
	test.registerFunction(0, f0);
	test.registerFunction(1, f1);
	test.registerCoeff(0, default_func);
	test.registerCoeff(1, default_func);
	test.set_timestep(1.0);
	test.set_initialcondition(0, 1);
	test.set_initialcondition(1, 0);
	test.getCurrentU().Display("un");
	
	test.print_header();
	test.solve_timestep();
	test.print_newresult();
	test.getNewU().Display("unp1"); //Working!!!
	
	Matrix<double> M(1,1), x(1,1), b(1,1);
	M(0,0) = 5;
	b(0,0) = 1;
	x.qrSolve(M, b);
	x.Display("x"); //Can use QR for 1-D system
	
	QR_DATA qrtest;
	success = QRsolve(test_matvec,b,&qrtest,NULL);
	qrtest.x.Display("another x");
	
	x.zeros();
	x.edit(0, 0, 0);
	PJFNK_DATA newtest;
	newtest.nl_maxit = 30;
	newtest.LineSearch = true;
	success = pjfnk(test_res, NULL, x, &newtest, NULL, NULL);
	newtest.x.Display("nonlin x"); //PJFNK now works with a single equation!!
	*/
	
	/**  ---------    Test 01: Various methods for First Order Decay (No Coupling) -------------- */
	Dove test01;
	FILE *file;
	file = fopen("output/DOVE_Tests.txt", "w+");
	if (file == nullptr)
	{
		system("mkdir output");
		file = fopen("output/DOVE_Tests.txt", "w+");
	}
	test01.set_outputfile(file);
	
	fprintf(file,"Test01: Single variable 1st Order Decay\n---------------------------------\ndu/dt = -u\n");
	
	test01.set_numfunc(1);
	test01.registerFunction(0, first_order_decay);
	test01.set_endtime(1.0);
	test01.set_timestepper(CONSTANT);
	test01.set_NonlinearOutput(false);
	test01.set_output(true);
	
	test01.set_initialcondition(0, 1);
	test01.set_timestep(0.05);
	test01.set_integrationtype(BE);
	test01.solve_all();
	
	test01.set_initialcondition(0, 1);
	test01.set_timestep(0.05);
	test01.set_integrationtype(FE);
	test01.solve_all();
	
	test01.set_initialcondition(0, 1);
	test01.set_timestep(0.05);
	test01.set_integrationtype(CN);
	test01.solve_all();
	
	test01.set_initialcondition(0, 1);
	test01.set_timestep(0.05);
	test01.set_integrationtype(BDF2);
	test01.solve_all();
	
	test01.set_initialcondition(0, 1);
	test01.set_timestep(0.05);
	test01.set_integrationtype(RK4);
	test01.solve_all();
	
	test01.set_initialcondition(0, 1);
	test01.set_timestep(0.05);
	test01.set_integrationtype(RKF);
	test01.solve_all();
	
	fprintf(file,"\n --------------- End of Test01 ---------------- \n\n");
	/**  ------------------------------    END Test01   ---------------------------------- */
	
	/**  ---------    Test 02: Various methods for Nonlinear Coupled ODEs-------------- */
	Dove test02;
	test02.set_outputfile(file);
	fprintf(file,"Test02: Two variable Nonlinear Decay\n---------------------------------\ndu1/dt = -u1\ndu2/dt = u1*u2\n");
	
	test02.set_numfunc(2);
	test02.registerFunction(0, first_order_decay);
	test02.registerFunction(1, nonlinear_first_order_decay);
	test02.set_endtime(4.0);
	test02.set_timestepper(ADAPTIVE);
	test02.set_timestepmax(0.2);
	test02.set_NonlinearOutput(true);
	test02.set_output(true);
	
	test02.set_initialcondition(0, 1);
	test02.set_initialcondition(1, 1);
	test02.set_timestep(0.01);
	test02.set_integrationtype(BE);
	test02.solve_all();
	
	test02.set_initialcondition(0, 1);
	test02.set_initialcondition(1, 1);
	test02.set_timestep(0.1);
	test02.set_integrationtype(FE);
	test02.solve_all();
	
	test02.set_initialcondition(0, 1);
	test02.set_initialcondition(1, 1);
	test02.set_timestep(0.1);
	test02.set_integrationtype(CN);
	test02.solve_all();
	
	test02.set_initialcondition(0, 1);
	test02.set_initialcondition(1, 1);
	test02.set_timestep(0.1);
	test02.set_integrationtype(BDF2);
	test02.solve_all();
	
	test02.set_initialcondition(0, 1);
	test02.set_initialcondition(1, 1);
	test02.set_timestep(0.1);
	test02.set_integrationtype(RK4);
	test02.solve_all();
	
	test02.set_initialcondition(0, 1);
	test02.set_initialcondition(1, 1);
	test02.set_timestep(0.1);
	test02.set_integrationtype(RKF);
	test02.solve_all();

	fprintf(file,"\n --------------- End of Test02 ---------------- \n\n");
	/**  ------------------------------    END Test02   ---------------------------------- */
	
	return success;
}
