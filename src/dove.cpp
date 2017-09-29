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
	dtmin = sqrt(DBL_EPSILON);
	dtmax = 100.0;
	int_type = IMPLICIT;
	int_sub = BE;
	timestepper = CONSTANT;
	Output = nullptr;
	num_func = 1;
	Converged = false;
	user_data = NULL;
	residual = residual_BE;
	precon = NULL;
}

//Default destructor
Dove::~Dove()
{
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
		this->user_jacobi.set_size(1, 1);
		this->num_func = 1;
	}
	else
	{
		this->un.set_size(i, 1);
		this->unp1.set_size(i, 1);
		this->unm1.set_size(i, 1);
		this->user_func.set_size(i, 1);
		this->user_coeff.set_size(i, 1);
		this->user_jacobi.set_size(i, i);
		this->num_func = i;
	}
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
			break;
			
		case CN:
			this->int_type = IMPLICIT;
			break;
			
		case BDF2:
			this->int_type = IMPLICIT;
			break;
			
		case RK4:
			this->int_type = EXPLICIT;
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
}

//Register user function
void Dove::registerFunction(int i, double (*func) (const Matrix<double> &u, const void *data) )
{
	if ((*func) == NULL)
	{
		mError(nullptr_func);
		this->user_func.edit(i, 0, default_jacobi);
	}
	else
		this->user_func.edit(i, 0, func);
}

//Register time coeff functions
void Dove::registerCoeff(int i, double (*coeff) (const Matrix<double> &u, const void *data) )
{
	if ((*coeff) == NULL)
		this->user_coeff.edit(i, 0, default_coeff);
	else
		this->user_coeff.edit(i, 0, coeff);
}

//Register jacobians
void Dove::registerJacobi(int i, int j, double (*jac) (const Matrix<double> &u, const void *data) )
{
	if ((*jac) == NULL)
	{
		if (i == j)
			this->user_jacobi.edit(i, j, default_coeff);
		else
			this->user_jacobi.edit(i, j, default_jacobi);
	}
	else
	{
		this->user_jacobi.edit(i, j, jac);
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
	fprintf(this->Output,"\nTime");
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

//Eval user function i
double Dove::Eval_Func(int i, const Matrix<double>& u)
{
	return this->user_func(i,0)(u,this->user_data);
}

//Eval user time coefficient function i
double Dove::Eval_Coeff(int i, const Matrix<double>& u)
{
	return this->user_coeff(i,0)(u,this->user_data);
}

//Eval user time coefficient function i
double Dove::Eval_Jacobi(int i, int j, const Matrix<double>& u)
{
	return this->user_jacobi(i,j)(u,this->user_data);
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
		
	}
	return success;
}

//Update solution states
void Dove::update_states()
{
	this->unm1 = this->un;
	this->un = this->unp1;
	this->dt_old = this->dt;
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
		Res(i,0) = (dat->Eval_Coeff(i, u)*u(i,0)) - (dat->Eval_Coeff(i, dat->getCurrentU())*dat->getCurrentU()(i,0)) - (dat->getTimeStep()*dat->Eval_Func(i, u));
	}
	
	return success;
}

/// Default time coefficient function
double default_coeff(const Matrix<double> &u, const void *data)
{
	return 1.0;
}

/// Default Jacobian element function
double default_jacobi(const Matrix<double> &u, const void *data)
{
	return 0.0;
}


// -------------------- Begin temporary testing --------------------------
double f0(const Matrix<double> &x, const void *res_data)
{
	return x(0,0) + 1;
}

double f1(const Matrix<double> &x, const void *res_data)
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
// -------------------- End temporary testing --------------------------

//Test function
int DOVE_TESTS()
{
	int success = 0;
	std::cout << "\nThis test is currently blank\n";
	
	//std::vector<double (*) (const Matrix<double> &x, const void *res_data)> list_func;
	//list_func.resize(2);
	//list_func[0] = f0;
	//list_func[1] = f1;
	Matrix<double (*) (const Matrix<double> &x, const void *res_data)> list_func;
	list_func.set_size(2, 1);
	list_func.edit(0, 0, f0);
	list_func.edit(1, 0, f1);
	
	Matrix<double> xtest;
	xtest.set_size(2, 1);
	xtest.edit(0, 0, 1);
	xtest.edit(1, 0, 2);
	xtest.Display("x");
	
	std::cout << "f0 = " << list_func(0,0)(xtest,NULL) << std::endl;
	std::cout << "f1 = " << list_func(1,0)(xtest,NULL) << std::endl;
	
	Dove test;
	FILE *testfile;
	testfile = fopen("output/Here_is_test.txt","w+");
	fprintf(testfile,"Here is a line of text.");
	test.set_outputfile(testfile);
	test.set_numfunc(2);
	test.registerFunction(0, f0);
	test.registerFunction(1, f1);
	test.registerCoeff(0, default_jacobi);
	test.registerCoeff(1, default_jacobi);
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
	
	
	return success;
}
