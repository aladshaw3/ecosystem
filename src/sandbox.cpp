/*!
 *  \file sandbox.h sandbox.cpp
 *	\brief Coding Test Area
 *
 *	\warning Functions and methods in this file are not meant to be used anywhere else.
 *  \author Austin Ladshaw
 *	\date 04/11/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
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

//Evaluation of a symmetric 1D set of orthonormal basis functions
double Bounded_1D_NormalBasis_Func(int i, double bound, double x)
{
	//Exact Eigenbasis
	return pow(2.0/bound, 0.5) * sin( (((double)i+1.0)*M_PI*x)/bound );
}

//Matrix elements for particle in a box with the above basis functions
double ParticleNBox_1D_H_MatrixElements_Bound1D_ONB(int i, int j, double bound, double mass, double h_bar)
{
	//Exact Eigenbasis
	if ( i == j)
		return (h_bar*h_bar)*(M_PI*M_PI)*((double)i)*((double)i)/mass/bound/bound/2.0;
	else
		return 0.0;
}

//Overlap elements
double ParticleNBox_1D_PHI_MatrixElements_Bound1D_ONB(int i, int j, double bound)
{
	//Exact Eigenbasis
	if (i == j)
		return 1.0;
	else
		return 0.0;
}

//Form the Hamiltonian matrix in the basis
int Form_H_Matrix_1DPIB_Bound1DBasis(PIB_1D_DATA *dat)
{
	int success = 0;
	
	for (int i=0; i<dat->m; i++)
	{
		for (int j=0; j<dat->m; j++)
		{
			dat->H.edit(i, j, ParticleNBox_1D_H_MatrixElements_Bound1D_ONB(i+1, j+1, dat->box_size, dat->mass, dat->h_bar));
		}
	}
	
	return success;
}

//Formation of overlap matrix
int Form_PHI_Matrix_1DPIB_Bound1DBasis(PIB_1D_DATA *dat)
{
	int success = 0;
	
	for (int i=0; i<dat->m; i++)
	{
		for (int j=0; j<dat->m; j++)
		{
			dat->PHI.edit(i, j, ParticleNBox_1D_PHI_MatrixElements_Bound1D_ONB(i+1, j+1, dat->box_size));
		}
	}
	
	return success;
}

//Evaluate residuals
int Eval_1DPIB_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	PIB_1D_DATA *dat = (PIB_1D_DATA *) data;
	
	double sum_sqs_i = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		double sum_sqs_j = 0.0;
		for (int j=0; j<dat->m; j++)
		{
			sum_sqs_j = sum_sqs_j + (x(i+1,0)*x(j+1,0)*dat->PHI(i,j));
		}
		sum_sqs_i	= sum_sqs_i + sum_sqs_j;
	}
	F.edit(0, 0, sum_sqs_i - 1.0);
	
	for (int i=0; i<dat->m; i++)
	{
		double sum_hij_cj = 0.0;
		for (int j=0; j<dat->m; j++)
		{
			sum_hij_cj = sum_hij_cj + (x(j+1,0) * ( dat->H(i,j) - (x(0,0)*dat->PHI(i,j)) ) );
		}
		F.edit(i+1, 0, sum_hij_cj);
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
	
	// ------------- Quantum Mechanics Example: Particle in a Box with Variational Method --------------------------
	
	std::cout << isEven(6) << "\t" << isEven(5) << std::endl;
	
	std::cout << "Solving the Schrodinger Equations for Particle in a Box, with Variational Method\n" << std::endl;
	
	PIB_1D_DATA pib_dat;
	pib_dat.m = 10;
	pib_dat.mass = 1.0;
	pib_dat.h_bar = 1.0;
	pib_dat.box_size = 1.0;     // i.e. bounds = 0 to 1
	pib_dat.N = pib_dat.m + 1;
	pib_dat.x.set_size(pib_dat.N, 1);
	pib_dat.c.set_size(pib_dat.m, 1);
	pib_dat.H.set_size(pib_dat.m, pib_dat.m);
	pib_dat.PHI.set_size(pib_dat.m, pib_dat.m);
	
	success = Form_H_Matrix_1DPIB_Bound1DBasis(&pib_dat);
	pib_dat.H.Display("H");
	
	success = Form_PHI_Matrix_1DPIB_Bound1DBasis(&pib_dat);
	pib_dat.PHI.Display("PHI");
	
	for (int i=0; i<pib_dat.x.rows(); i++)
	{
		if ( i == 0)
			pib_dat.x.edit(i, 0, 1.0);
		else
			pib_dat.x.edit(i, 0, 1.0 / (pib_dat.x.rows()-1.0));
	}
	/*
	pib_dat.x.edit(0, 0, 4.9348);
	pib_dat.x.edit(1, 0, 0.00127891);
	pib_dat.x.edit(2, 0, 2.53624706);
	pib_dat.x.edit(3, 0, 0.15174395);
	pib_dat.x.edit(4, 0, -2.9902389);
	pib_dat.x.edit(5, 0, -0.1474352);
	pib_dat.x.edit(6, 0, 1.59611461);
	pib_dat.x.edit(7, 0, 0.03853069);
	pib_dat.x.edit(8, 0, -0.8377256);
	pib_dat.x.edit(9, 0, 0.34744869);
	pib_dat.x.edit(10, 0, 0.0);
	*/
	Matrix<double> res(pib_dat.N,1);
	success = Eval_1DPIB_Residuals(pib_dat.x, res, &pib_dat);
	res.Display("F(x)");
	
	PJFNK_DATA QM_dat;
	QM_dat.linear_solver = GMRESRP;
	QM_dat.LineSearch = true;
	QM_dat.Bounce = false;
	QM_dat.gmresrp_dat.restart = pib_dat.N;
	QM_dat.nl_maxit = 100;
	QM_dat.nl_tol_abs = 1e-6;
	QM_dat.nl_tol_rel = 1e-6;
	QM_dat.lin_tol_abs = 1e-6;
	QM_dat.lin_tol_rel = 1e-6;
	//QM_dat.L_Output = true;
	
	success = pjfnk(Eval_1DPIB_Residuals, NULL, pib_dat.x, &QM_dat, (void *)&pib_dat, NULL);
	
	for (int i=1; i<pib_dat.m; i++)
	{
		pib_dat.c.edit(i-1, 0, pib_dat.x(i,0));
	}
	
	pib_dat.c.Display("c");
	
	std::cout << "Approximate Eo = " << pib_dat.x(0,0) << std::endl;
	std::cout << "Exact Eo = " << ((M_PI*M_PI)*(pib_dat.h_bar*pib_dat.h_bar))/(2.0*pib_dat.mass*pib_dat.box_size*pib_dat.box_size) << std::endl;
	
	std::cout << "\nThis demonstrates that the variational method will always approximate the lowest energy state of the system\n\n";
	
	
	// -------------------------------------- End Quantum Mechanics Example ----------------------------------------
	
	std::cout << "\nEnd SANDBOX\n\n";
	
	return success;
}