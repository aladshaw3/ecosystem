//----------------------------------------
//  Created by Austin Ladshaw on 04/11/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 *		The SANDBOX is a clean working space where new algorithms and ideas can be tested before
 *		attempting to have those ideas implemented in the FLOCK or SCHOOL. 
 *
 */

#include "sandbox.h"

//Test function number 01 for speciation problems
int Speciation_Test01_Function(const Matrix<double> &x, Matrix<double> &F, const void *res_data)
{
	int success = 0;
	Speciation_Test01_Data *dat = (Speciation_Test01_Data *) res_data;
	
	F(0,0) = dat->logKa1 + x(0,0) - x(3,0) - x(1,0);
	F(1,0) = dat->logKa2 + x(1,0) - x(3,0) - x(2,0);
	F(2,0) = 1.0 - ( (pow(10.0, x(0,0))+pow(10.0, x(1,0))+pow(10.0, x(2,0)) ) / dat->CT);
	F(3,0) = 1.0 - ( (pow(10.0, x(1,0))+(2.0*pow(10.0, x(2,0)))+pow(10.0, dat->logKw-x(3,0))) / (dat->NaT + pow(10.0, x(3,0))) );
	
	return success;
}

//Analytical Jacobian for Test function number 01
int Speciation_Test01_Jacobian(const Matrix<double> &x, Matrix<double> &J, const void *precon_data)
{
	int success = 0;
	Speciation_Test01_Data *dat = (Speciation_Test01_Data *) precon_data;
	
	J(0,0) = 1.0; J(0,1) = -1.0; J(0,2) = 0.0; J(0,3) = -1.0;
	J(1,0) = 0.0; J(1,1) = 1.0; J(1,2) = -1.0; J(1,3) = -1.0;
	
	J(2,0) = -pow(10.0, x(0,0))/dat->CT;
	J(2,1) = -pow(10.0, x(1,0))/dat->CT;
	J(2,2) = -pow(10.0, x(2,0))/dat->CT;
	J(2,3) = 0.0;
	
	J(3,0) = 0.0;
	J(3,1) = -pow(10.0, x(1,0))/(dat->NaT + pow(10.0, x(3,0)));
	J(3,2) = -(2.0*pow(10.0, x(2,0)))/(dat->NaT + pow(10.0, x(3,0)));
	J(3,3) = ( ( (pow(10.0, x(1,0))+(2.0*pow(10.0, x(2,0))))*pow(10.0, x(3,0)) ) + (dat->NaT*pow(10.0, dat->logKw-x(3,0))) + (2.0*pow(10.0, dat->logKw)) ) / pow(dat->NaT+pow(10.0, x(3,0)), 2.0);
	
	return success;
}

//Make initial guess for problem 1
int Speciation_Test01_Guess(const void *user_data)
{
	int success = 0;
	Speciation_Test01_Data *dat = (Speciation_Test01_Data *) user_data;
	
	dat->logC(0,0) = log10(0.05*dat->CT);
	dat->logC(1,0) = log10(0.99*dat->CT);
	dat->logC(2,0) = log10(0.05*dat->CT);
	dat->logC(3,0) = -7.0;
	
	return success;
}

//Simple MatVec product
int Speciation_Test01_MatVec(const Matrix<double> &x, Matrix<double> &Ax, const void *matvec_data)
{
	int success = 0;
	Speciation_Test01_Data *dat = (Speciation_Test01_Data *) matvec_data;
	
	Ax = dat->Jacobian*x;
	
	return success;
}

//Form a numerical jacobian matrix
int NumericalJacobian( int (*Func) (const Matrix<double> &x, Matrix<double> &F, const void *user_data),
					  const Matrix<double> &x, Matrix<double> &J, int Nx, int Nf, NUM_JAC_DATA *jac_dat,
					  const void *user_data)
{
	int success = 0;
	
	//Check input arguments for problems
	if ( (*Func) == NULL )
	{
		mError(nullptr_func);
		return -1;
	}
	if ( jac_dat == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	if (jac_dat->dxj.rows() != Nx)
	{
		jac_dat->dxj.set_size(Nx, 1);
	}
	if (J.rows() != Nf && J.columns() != Nx)
	{
		J.set_size(Nf, Nx);
	}
	if (jac_dat->Fx.rows() != Nf)
	{
		jac_dat->Fx.set_size(Nf, 1);
	}
	if (jac_dat->Fxp.rows() != Nf)
	{
		jac_dat->Fxp.set_size(Nf, 1);
	}
	if (jac_dat->eps < sqrt(DBL_EPSILON) || jac_dat->eps >= 1.0)
	{
		jac_dat->eps = sqrt(DBL_EPSILON);
	}
	
	//Form the first fuction evaluation
	success = (*Func) (x, jac_dat->Fx, user_data);
	if (success != 0) {mError(simulation_fail); return -1;}
	
	//Create a copy of x to change piecewise
	jac_dat->dxj = x;
	for (int j=0; j<Nx; j++)
	{
		//Change the jth variable
		jac_dat->dxj(j,0) = x(j,0) + jac_dat->eps;
		success = (*Func) (jac_dat->dxj, jac_dat->Fxp, user_data);
		if (success != 0) {mError(simulation_fail); return -1;}
		
		//Approximate each row of the first column of the jacobian
		for (int i=0; i<Nf; i++)
		{
			J(i,j) = (jac_dat->Fxp(i,0) - jac_dat->Fx(i,0)) / jac_dat->eps;
		}
		
		//Recover the jth variable before continuing
		jac_dat->dxj(j,0) = x(j,0);
	}
	
	return success;
}

//Run the sandbox tests
int RUN_SANDBOX()
{
	int success = 0;
	
	std::cout << "Start SANDBOX\n\n";
	
	Speciation_Test01_Data dat01;
	dat01.x.resize(dat01.N+1);
	dat01.x[0].Register("H2CO3 (aq)");
	dat01.x[1].Register("HCO3 - (aq)");
	dat01.x[2].Register("CO3 2- (aq)");
	dat01.x[3].Register("H + (aq)");
	dat01.x[4].Register("OH - (aq)");
	dat01.Jacobian.set_size(dat01.N, dat01.N);
	dat01.NumJac.set_size(dat01.N, dat01.N);
	dat01.logC.set_size(dat01.N, 1);
	dat01.C.set_size(dat01.N, 1);
	
	for (int i=0; i<dat01.x.size(); i++)
		dat01.x[i].DisplayInfo();
	
	//Calculate the constants for the problem
	
	//Make initial guess to solution
	success = Speciation_Test01_Guess((void *)&dat01);
	dat01.logC.Display("logCo");
	
	//Convert to mol/L
	success = Convert2Concentration(dat01.logC, dat01.C);
	dat01.C.Display("Co");
	
	//Test Jacobian formation
	success = Speciation_Test01_Jacobian(dat01.logC, dat01.Jacobian, (void *)&dat01);
	dat01.Jacobian.Display("Jac01");
	
	//Form a temp BC matrix to solve Jacobian iteratively
	Matrix<double> b;
	b.set_size(dat01.N, 1);
	b.edit(0, 0, 1);
	GMRESRP_DATA gmres_dat01;
	gmres_dat01.x.set_size(dat01.N, 1);
	gmres_dat01.x.edit(0, 0, 1);
	success = gmresRightPreconditioned(Speciation_Test01_MatVec, NULL, b, &gmres_dat01, (void *)&dat01, NULL);
	gmres_dat01.x.Display("x_GMRES_test");
	
	GMRESR_DATA gmresr_dat01;
	gmresr_dat01.gmres_restart = 3;
	gmresr_dat01.gcr_dat.x.set_size(dat01.N, 1);
	gmresr_dat01.gcr_dat.x.edit(0, 0, 1);
	success = gmresr(Speciation_Test01_MatVec, NULL, b, &gmresr_dat01, (void *)&dat01, NULL);
	gmresr_dat01.gcr_dat.x.Display("x_GMRESR_test");
	
	//Form a numerical jacobian and compare
	NUM_JAC_DATA jac_dat01;
	dat01.logC.Display("logCo");
	success = NumericalJacobian(Speciation_Test01_Function, dat01.logC, dat01.NumJac, dat01.N, dat01.N, &jac_dat01, (void *)&dat01);
	dat01.NumJac.Display("NumJac");
	
	//Check the jacobian
	Matrix<double> Check;
	Check = dat01.Jacobian - dat01.NumJac;
	Check.Display("Check");
	std::cout << "Check norm = " << Check.norm() << std::endl;
	
	//Run PJFNK optimization
	PJFNK_DATA jfnk_dat;
	//jfnk_dat.L_Output = true;
	jfnk_dat.linear_solver = GMRESR;	//Always use GMRESRP, GCR, or GMRESR for these problems!
	jfnk_dat.nl_maxit = 100;			//If you have a bad initial guess, up the max iterations and turn on linesearch
	//jfnk_dat.LineSearch = true;			//If you have a bad initial guess, up the max iterations and turn on linesearch
	//jfnk_dat.lin_tol = 1e-6;			//Slightly better convergence if taking near exact newton step
	//jfnk_dat.nl_tol_abs = 1e-10;		//All other defaults are fine
	//jfnk_dat.nl_tol_rel = 1e-10;		//All other defaults are fine
	success = pjfnk(Speciation_Test01_Function, NULL, dat01.logC, &jfnk_dat, (void *)&dat01, (void *)&dat01);
	dat01.logC.Display("logC_final");
	
	//Convert to mol/L
	success = Convert2Concentration(dat01.logC, dat01.C);
	dat01.C.Display("C_final");
	
	
	
	std::cout << "\nEnd SANDBOX\n\n";
	
	return success;
}