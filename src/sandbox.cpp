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

//Evaluation the polynomial basis function at x
double Eval_PolyBasisFunc(int i, double x)
{
	if (i < 1)
		return 1.0;
	else
		return pow(x, (double)i);
}

//Evaluation of the 1st derivative of the basis
double Eval_1stDerivative_PolyBasisFunc(int i, double x)
{
	if (i < 1)
		return 0.0;
	else if (i == 1)
		return 1.0;
	else
		return (double)i * pow(x, (double)i-1.0);
}

//Evaluation of 2nd derivative of basis
double Eval_2ndDerivative_PolyBasisFunc(int i, double x)
{
	if (i < 2)
		return 0.0;
	else if (i == 2)
		return 2.0;
	else
		return (double)i * ((double)i - 1.0) * pow(x, (double)i - 2.0);
}

//Evaluate the poly function
double Eval_ApproximatePolySolution(Matrix<double> &c, double x)
{
	double sum = 0.0;
	for (int i=0; i<c.rows(); i++)
	{
		sum = sum + (c(i,0)*Eval_PolyBasisFunc(i, x));
	}
	return sum;
}

//Gradient integration
double Gradient_Integral_PolyBasis(int i, int j, double lower, double upper)
{
	double Bij = 0.0;
	
	if (j == 0)
		return Bij;
	else
	{
		double exp = (double)(i+j);
		double pre = (double)(j);
		Bij = pre * ((pow(upper, exp)/exp) - (pow(lower, exp)/exp));
	}
	
	return Bij;
}

//Laplacian integration
double Laplacian_Integral_PolyBasis(int i, int j, double lower, double upper)
{
	double Aij = 0.0;
	
	if (j==0 || j==1)
		return Aij;
	else
	{
		double exp = (double)(i+j)-1.0;
		double pre = (double)(j*(j-1));
		Aij = pre * ((pow(upper, exp)/exp) - (pow(lower, exp)/exp));
	}
	
	return Aij;
}

//Overlap integrals
double Overlap_Integral_PolyBasis(int i, int j, double lower, double upper)
{
	double Oij = 0.0;
	
	double exp = (double)(i+j+1);
	Oij = ((pow(upper, exp)/exp) - (pow(lower, exp)/exp));
	
	return Oij;
}

//Residuals for Variational Polynomial Approximation method Test 1
int Eval_VPA_Test_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	VPA_Test_DATA *dat = (VPA_Test_DATA *) data;
	
	//First two residuals are derivative of lambda for constraints
	double res0 = 0.0, res1 = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		res0 = res0 + (x(i+2,0)*Eval_PolyBasisFunc(i, 0.0));
		res1 = res1 + (x(i+2,0)*Eval_1stDerivative_PolyBasisFunc(i, dat->L));
		//res1 = res1 + (x(i+2,0)*Eval_PolyBasisFunc(i, dat->L));
	}
	F.edit(0, 0, res0 - dat->uo);
	F.edit(1, 0, res1);
	
	double Lap_sum = 0.0, Over_sum = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		Lap_sum = 0.0;
		Over_sum = 0.0;
		for (int j=0; j<dat->m; j++)
		{
			Lap_sum = Lap_sum + (x(j+2,0)*Laplacian_Integral_PolyBasis(i, j, 0.0, dat->L));
			Over_sum = Over_sum + (x(j+2,0)*Overlap_Integral_PolyBasis(i, j, 0.0, dat->L));
		}
		F.edit(i+2, 0, Lap_sum - ((dat->k/dat->D)*Over_sum) + (x(0,0)*Eval_PolyBasisFunc(i, 0.0)) + (x(1,0)*Eval_1stDerivative_PolyBasisFunc(i, dat->L)));
	}
	
	return success;
}

//Evaluate Trig Basis at x
/*
	Basis: 
 
	phi_i(x) = 
				cos(n*x)	with n=i/2			for i=0,2,4,6...
				sin(n*x)	with n=(i+1)/2		for i=1,3,5,7...
 */
double Eval_TrigBasisFunc(int i, double x)
{
	double phi = 0.0;
	
	if (i == 0)
		phi = 1.0;
	else
		(i % 2 == 0) ? phi = cos((double)(i/2)*x): phi = sin((double)((i+1)/2)*x);
	
	return phi;
}

//Evaluation of the 1st derivative of the trig basis
double Eval_1stDerivative_TrigBasisFunc(int i, double x)
{
	double der = 0.0;
	if (i < 1)
		der = 0.0;
	else
		(i % 2 == 0) ? der = -(double)(i/2)*sin((double)(i/2)*x): der = (double)((i+1)/2)*cos((double)((i+1)/2)*x);
	
	return der;
}

//Evaluation of the 2nd derivative of the trig basis
double Eval_2ndDerivative_TrigBasisFunc(int i, double x)
{
	double der = 0.0;
	if (i < 1)
		der = 0.0;
	else
		(i % 2 == 0) ? der= -(double)(i/2)*(double)(i/2)*cos((double)(i/2)*x):der= -(double)((i+1)/2)*(double)((i+1)/2)*sin((double)((i+1)/2)*x);
	
	return der;
}

//Evaluate the trig function
double Eval_ApproximateTrigSolution(Matrix<double> &c, double x)
{
	double sum = 0.0;
	for (int i=0; i<c.rows(); i++)
	{
		sum = sum + (c(i,0)*Eval_TrigBasisFunc(i, x));
	}
	return sum;
}

//Gradient integration
double Gradient_Integral_TrigBasis(int i, int j, double lower, double upper)
{
	// Calculate n and m and determine function type
	int n = 0, m = 0;
	trig_type type_i, type_j;
	if (i > 0) (i % 2 == 0) ? n = (i/2): n = ((i+1)/2);
	if (j > 0) (j % 2 == 0) ? m = (j/2): m = ((j+1)/2);
	
	if (i > 0) (i % 2 == 0) ? type_i = COS: type_i = SIN;
	else type_i = CONST;
	if (j > 0) (j % 2 == 0) ? type_j = SIN: type_j = COS;
	else type_j = CONST;
	
	//Note: if type_j == CONST, then Bij = 0, otherwise Bij = Oij*coeff, coeff = +/-m
	double coeff = 0.0;
	switch (type_j)
	{
		case CONST:
			coeff = 0.0;
			break;
			
		case SIN:
			coeff = (double)m;
			break;
			
		case COS:
			coeff = -(double)m;
			break;
	}
	
	return coeff*TrigBasis_Integrals(type_i, type_j, n, m, lower, upper);
}

//Laplacian integration
double Laplacian_Integral_TrigBasis(int i, int j, double lower, double upper)
{
	// Calculate n and m and determine function type
	int n = 0, m = 0;
	trig_type type_i, type_j;
	if (i > 0) (i % 2 == 0) ? n = (i/2): n = ((i+1)/2);
	if (j > 0) (j % 2 == 0) ? m = (j/2): m = ((j+1)/2);
	
	if (i > 0) (i % 2 == 0) ? type_i = COS: type_i = SIN;
	else type_i = CONST;
	if (j > 0) (j % 2 == 0) ? type_j = COS: type_j = SIN;
	else type_j = CONST;
	
	//Note: if type_j == CONST, then Bij = 0, otherwise Bij = Oij*coeff, coeff = -m*m
	double coeff = 0.0;
	switch (type_j)
	{
		case CONST:
			coeff = 0.0;
			break;
			
		case SIN:
			coeff = -(double)(m*m);
			break;
			
		case COS:
			coeff = -(double)(m*m);
			break;
	}
	
	return coeff*TrigBasis_Integrals(type_i, type_j, n, m, lower, upper);
}

//Overlap integrals
double Overlap_Integral_TrigBasis(int i, int j, double lower, double upper)
{
	// Calculate n and m and determine function type
	int n = 0, m = 0;
	trig_type type_i, type_j;
	if (i > 0) (i % 2 == 0) ? n = (i/2): n = ((i+1)/2);
	if (j > 0) (j % 2 == 0) ? m = (j/2): m = ((j+1)/2);
	
	if (i > 0) (i % 2 == 0) ? type_i = COS: type_i = SIN;
	else type_i = CONST;
	if (j > 0) (j % 2 == 0) ? type_j = COS: type_j = SIN;
	else type_j = CONST;
	
	/*
	std::cout << n << "\t" << m << std::endl;
	switch (type_i)
	{
		case CONST:
			std::cout << "CONST\t";
			break;
			
		case COS:
			std::cout << "COS\t";
			break;
			
		case SIN:
			std::cout << "SIN\t";
			break;
	}
	switch (type_j)
	{
		case CONST:
			std::cout << "CONST\n";
			break;
			
		case COS:
			std::cout << "COS\n";
			break;
			
		case SIN:
			std::cout << "SIN\n";
			break;
	}
	*/
	return TrigBasis_Integrals(type_i, type_j, n, m, lower, upper);
}

//Function to compute a trig*trig integral from lower to upper bound in 1-D
double TrigBasis_Integrals(trig_type type_i, trig_type type_j, int n, int m, double lower, double upper)
{
	double Oij = 0.0;
	bool equal;
	if (n == m) equal = true;
	else equal = false;
	
	//Switch case for possible integrals
	switch (type_i)
	{
		case CONST:
		{
			switch (type_j)
			{
				case CONST:
					Oij = upper - lower;
					break;
					
				case SIN:
					Oij = -(cos(upper*(double)m) - cos(lower*(double)m))/(double)m;
					break;
					
				case COS:
					Oij = (sin(upper*(double)m) - sin(lower*(double)m))/(double)m;
					break;
			}
		}
			break;
			
		case SIN:
		{
			switch (type_j)
			{
				case CONST:
					Oij = -(cos(upper*(double)n) - cos(lower*(double)n))/(double)n;
					break;
					
				case SIN:
					if (equal == true)
						Oij = ((2.0*(double)n*(upper-lower)) + sin(2.0*lower*(double)n) - sin(2.0*upper*(double)n))/(4.0*(double)n);
					else
						Oij = ( -((double)n*sin(lower*(double)m)*cos(lower*(double)n)) + ((double)m*cos(lower*(double)m)*sin(lower*(double)n)) + ((double)n*sin(upper*(double)m)*cos(upper*(double)n)) - ((double)m*cos(upper*(double)m)*sin(upper*(double)n)) )/((double)(m*m) - (double)(n*n));
					break;
					
				case COS:
					if (equal == true)
						Oij = (cos(2.0*lower*(double)n) - cos(2.0*upper*(double)n))/(4.0*(double)n);
					else
						Oij = ( -((double)m*sin(lower*(double)m)*sin(lower*(double)n)) - ((double)n*cos(lower*(double)m)*cos(lower*(double)n)) + ((double)m*sin(upper*(double)m)*sin(upper*(double)n)) + ((double)n*cos(upper*(double)m)*cos(upper*(double)n)) )/((double)(m*m) - (double)(n*n));
					break;
			}
		}
			break;
			
		case COS:
		{
			switch (type_j)
			{
				case CONST:
					Oij = (sin(upper*(double)n) - sin(lower*(double)n))/(double)n;
					break;
					
				case SIN:
					if (equal == true)
						Oij = (cos(2.0*lower*(double)n) - cos(2.0*upper*(double)n))/(4.0*(double)n);
					else
						Oij = ( ((double)n*sin(lower*(double)m)*sin(lower*(double)n)) + ((double)m*cos(lower*(double)m)*cos(lower*(double)n)) - ((double)n*sin(upper*(double)m)*sin(upper*(double)n)) - ((double)m*cos(upper*(double)m)*cos(upper*(double)n)) )/((double)(m*m) - (double)(n*n));
					break;
					
				case COS:
					if (equal == true)
						Oij = ((2.0*(double)n*(upper-lower)) - sin(2.0*lower*(double)n) + sin(2.0*upper*(double)n))/(4.0*(double)n);
					else
						Oij = ( -((double)m*sin(lower*(double)m)*cos(lower*(double)n)) + ((double)n*cos(lower*(double)m)*sin(lower*(double)n)) + ((double)m*sin(upper*(double)m)*cos(upper*(double)n)) - ((double)n*cos(upper*(double)m)*sin(upper*(double)n)) )/((double)(m*m) - (double)(n*n));
					break;
			}
		}
			break;
	}
	
	return Oij;
}

//Residuals for Variational Polynomial Approximation method Test 2
int Eval_VPA_Test02_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	VPA_Test02_DATA *dat = (VPA_Test02_DATA *) data;
	
	//Note: x(0) = lambda_0 and x(1) = lambda_1
	
	//First two residuals are derivative of lambda for constraints
	double res0 = 0.0, res1 = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		res0 = res0 + (x(i+2,0)*Eval_PolyBasisFunc(i, 0.0));
		res1 = res1 + (x(i+2,0)*Eval_1stDerivative_PolyBasisFunc(i, dat->L));
	}
	F.edit(0, 0, res0 - dat->uo);
	F.edit(1, 0, res1);
	
	double Grad_sum_np1 = 0.0, Over_sum_np1 = 0.0, Over_sum_n = 0.0, Lap_sum_np1 = 0.0;
	double Grad_sum_n = 0.0, Lap_sum_n = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		Grad_sum_np1 = 0.0;
		Lap_sum_np1 = 0.0;
		Over_sum_np1 = 0.0;
		Over_sum_n = 0.0;
		Grad_sum_n = 0.0;
		Lap_sum_n = 0.0;
		for (int j=0; j<dat->m; j++)
		{
			Grad_sum_np1 = Grad_sum_np1 + (x(j+2,0)*Gradient_Integral_PolyBasis(i, j, 0.0, dat->L));
			Grad_sum_n = Grad_sum_n + (dat->xn(j+2,0)*Gradient_Integral_PolyBasis(i, j, 0.0, dat->L));
			Over_sum_np1 = Over_sum_np1 + (x(j+2,0)*Overlap_Integral_PolyBasis(i, j, 0.0, dat->L));
			Over_sum_n = Over_sum_n + (dat->xn(j+2,0)*Overlap_Integral_PolyBasis(i, j, 0.0, dat->L));
			Lap_sum_np1 = Lap_sum_np1 + (x(j+2,0)*Laplacian_Integral_PolyBasis(i, j, 0.0, dat->L));
			Lap_sum_n = Lap_sum_n + (dat->xn(j+2,0)*Laplacian_Integral_PolyBasis(i, j, 0.0, dat->L));
		}
		F.edit(i+2, 0, Over_sum_np1 - Over_sum_n + (dat->beta*dat->dt*dat->v*Grad_sum_np1) + ((1.0-dat->beta)*dat->dt*dat->v*Grad_sum_n) - (dat->beta*dat->dt*dat->D*Lap_sum_np1) - ((1.0-dat->beta)*dat->dt*dat->D*Lap_sum_n) + (x(0,0)*Eval_PolyBasisFunc(i, 0.0)) + (x(1,0)*Eval_1stDerivative_PolyBasisFunc(i, dat->L)));
	}

	
	return success;
}

//Residuals for Variational Polynomial Approximation method Test 2
int Eval_VPA_Test03_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	VPA_Test02_DATA *dat = (VPA_Test02_DATA *) data;
	
	//Note: x(0) = lambda_0 and x(1) = lambda_1
	
	//First two residuals are derivative of lambda for constraints
	double res0 = 0.0, res1 = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		res0 = res0 + (x(i+2,0)*Eval_TrigBasisFunc(i, 0.0));
		res1 = res1 + (x(i+2,0)*Eval_1stDerivative_TrigBasisFunc(i, dat->L));
	}
	F.edit(0, 0, res0 - dat->uo);
	F.edit(1, 0, res1);
	
	double Grad_sum_np1 = 0.0, Over_sum_np1 = 0.0, Over_sum_n = 0.0, Lap_sum_np1 = 0.0;
	double Grad_sum_n = 0.0, Lap_sum_n = 0.0;
	for (int i=0; i<dat->m; i++)
	{
		Grad_sum_np1 = 0.0;
		Lap_sum_np1 = 0.0;
		Over_sum_np1 = 0.0;
		Over_sum_n = 0.0;
		Grad_sum_n = 0.0;
		Lap_sum_n = 0.0;
		for (int j=0; j<dat->m; j++)
		{
			Grad_sum_np1 = Grad_sum_np1 + (x(j+2,0)*Gradient_Integral_TrigBasis(i, j, 0.0, dat->L));
			Grad_sum_n = Grad_sum_n + (dat->xn(j+2,0)*Gradient_Integral_TrigBasis(i, j, 0.0, dat->L));
			Over_sum_np1 = Over_sum_np1 + (x(j+2,0)*Overlap_Integral_TrigBasis(i, j, 0.0, dat->L));
			Over_sum_n = Over_sum_n + (dat->xn(j+2,0)*Overlap_Integral_TrigBasis(i, j, 0.0, dat->L));
			Lap_sum_np1 = Lap_sum_np1 + (x(j+2,0)*Laplacian_Integral_TrigBasis(i, j, 0.0, dat->L));
			Lap_sum_n = Lap_sum_n + (dat->xn(j+2,0)*Laplacian_Integral_TrigBasis(i, j, 0.0, dat->L));
			
		}
		F.edit(i+2, 0, Over_sum_np1 - Over_sum_n + (dat->beta*dat->dt*dat->v*Grad_sum_np1) + ((1.0-dat->beta)*dat->dt*dat->v*Grad_sum_n) - (dat->beta*dat->dt*dat->D*Lap_sum_np1) - ((1.0-dat->beta)*dat->dt*dat->D*Lap_sum_n) + (x(0,0)*Eval_TrigBasisFunc(i, 0.0)) + (x(1,0)*Eval_1stDerivative_TrigBasisFunc(i, dat->L)));
	}
	
	
	return success;
}

//Evaluation of 2D poly basis functions
double PolyBasis_2D(int i, int j, double x, double y)
{
	double xpart = 0.0, ypart = 0.0;
	if (i < 1)
		xpart = 1.0;
	else
		xpart = pow(x, (double)i);
	
	if (j < 1)
		ypart = 1.0;
	else
		ypart = pow(y, (double)j);
	return xpart*ypart;
}

//First derivative in x
double PolyBasis_2D_dx(int i, int j, double x, double y)
{
	return (double)i*PolyBasis_2D(i-1, j, x, y);
}

//Second derivative in x
double PolyBasis_2D_dx2(int i, int j, double x, double y)
{
	return (double)i*((double)i-1.0)*PolyBasis_2D(i-2, j, x, y);
}

//First derivative in y
double PolyBasis_2D_dy(int i, int j, double x, double y)
{
	return (double)j*PolyBasis_2D(i, j-1, x, y);
}

//Second derivative in y
double PolyBasis_2D_dy2(int i, int j, double x, double y)
{
	return (double)j*((double)j-1.0)*PolyBasis_2D(i, j-2, x, y);
}

//First derivative in x and y
double PolyBasis_2D_dxdy(int i, int j, double x, double y)
{
	return (double)i*(double)j*PolyBasis_2D(i-1, j-1, x, y);
}

//Evaluate the 2D polynomial approximation
double PolyBasis_2D_LinearComboAppox(Matrix<double> &c, double x, double y)
{
	double sum = 0.0;
	
	for (int i=0; i<c.rows(); i++)
	{
		for (int j=0; j<c.columns(); j++)
		{
			sum = sum + (c(i,j)*PolyBasis_2D(i,j, x, y));
		}
	}
	
	return sum;
}

//Evaluation of the 2D Laplacian integrals
double Laplacian_Integral_PolyBasis_2D(int i, int j, int l, int m, double x_low, double x_high, double y_low, double y_high)
{
	double Aij_lm = 0.0;
	double yx_exp = (double)j+m+1;
	double xx_exp = (double)i+l-1;
	double yy_exp = (double)j+m-1;
	double xy_exp = (double)i+l+1;
	
	double xpiece = 0.0;
	if (l == 0 || l == 1)
		xpiece = 0.0;
	else
		xpiece = (double)l*((double)l-1.0)*(pow(y_high, yx_exp) - pow(y_low, yx_exp))*(pow(x_high, xx_exp) - pow(x_low, xx_exp))/(yx_exp*xx_exp);
	
	double ypiece = 0.0;
	if (m == 0 || m == 1)
		ypiece = 0.0;
	else
		ypiece = (double)m*((double)m-1.0)*(pow(y_high, yy_exp) - pow(y_low, yy_exp))*(pow(x_high, xy_exp) - pow(x_low, xy_exp))/(yy_exp*xy_exp);
	
	Aij_lm = xpiece + ypiece;
	
	
	return Aij_lm;
}

//Evaluation of overlap integrals in 2D poly basis
double Overlap_Integral_PolyBasis_2D(int i, int j, int l, int m, double x_low, double x_high, double y_low, double y_high)
{
	double Oij_lm = 0.0;
	double y_exp = (double)j+m+1;
	double x_exp = (double)i+l+1;
	
	Oij_lm = (pow(y_high, y_exp)-pow(y_low, y_exp))*(pow(x_high, x_exp)-pow(x_low, x_exp))/(y_exp*x_exp);
	
	return Oij_lm;
}

int Eval_2D_VPA_TEST_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	/*
	VPA_2D_TEST_DATA *dat = (VPA_2D_TEST_DATA *) data;
	
	double Dij_sum = 0.0;
	
	int ij = 0;
	for (int i=0; i<dat->m; i++)
	{// i loop
		
		for (int j=0; j<dat->m; j++)
		{// j loop
			
			ij = (i*dat->m)+j+dat->bcs;
			double Alm_sum = 0.0;
			double Olm_sum = 0.0;
			
			int lm = 0;
			for (int l=0; l<dat->m; l++)
			{// l loop
				
				for (int m=0; m<dat->m; m++)
				{// m loop
					
					lm = (l*dat->m)+m+dat->bcs;
					Alm_sum = Alm_sum + (x(lm,0)*Laplacian_Integral_PolyBasis_2D(i, j, l, m, 0.0, dat->Lx, 0.0, dat->Ly));
					Olm_sum = Olm_sum + (x(lm,0)*Overlap_Integral_PolyBasis_2D(i, j, l, m, 0.0, dat->Lx, 0.0, dat->Ly));
					
				}// end m loop
				
			}// end l loop
			
			double l0 = PolyBasis_2D(i, j, 0, 0);
			
			Dij_sum = Dij_sum + (x(ij,0)*l0);
			
			F.edit(ij, 0, (dat->D*Alm_sum)-(dat->k*Olm_sum)+(x(0,0)*l0));
	 
		
	}// end i loop
	
	F.edit(0, 0, Dij_sum - dat->uo);
	
	*/
	return success;
}

//Example for QR
int ex1_mult(const Matrix<double>& x, Matrix<double> &Ax, const void *data)
{
	int success = 0;
	
	QR_EX1 *dat = (QR_EX1 *) data;
	
	Ax = dat->M * x;
	
	return success;
}

extern "C"
{
    int blah() {return RUN_SANDBOX();}
    double obj_func(double *list, int len, double *args)
    {
        double sum = 0.0;
        for (int i=0; i<len; i++)
            sum += args[i];
        return sum;
    }
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
	jfnk_dat.linear_solver = QR;	//Always use GMRESRP, GCR, or GMRESR for these problems!
	jfnk_dat.nl_maxit = 100;			//If you have a bad initial guess, up the max iterations and turn on linesearch
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
	
	
	// ----------------------------- Example of Varitational Polynomial Approximation ------------------------------
	
	std::cout << "Solve d^2/dx^2 (u) = (k/D) * u Approximately with a constrained variational method...\n\n";
	
	VPA_Test_DATA vpa_dat;
	vpa_dat.m = 3;
	vpa_dat.N = vpa_dat.m + 2;
	vpa_dat.L = 1.0;
	vpa_dat.k = 3.0;
	vpa_dat.D = 1.0;
	vpa_dat.uo = 1.0;
	vpa_dat.c.set_size(vpa_dat.m, 1);
	vpa_dat.x.set_size(vpa_dat.N, 1);
	Matrix<double> testb(vpa_dat.N,1);
	
	PJFNK_DATA vpa_newton;
	vpa_newton.linear_solver = QR;
	vpa_newton.gmresrp_dat.restart = vpa_dat.N;
	vpa_newton.LineSearch = true;
	vpa_newton.nl_tol_abs = 1e-6;
	vpa_newton.nl_tol_rel = 1e-6;
	
	success = pjfnk(Eval_VPA_Test_Residuals, NULL, vpa_dat.x, &vpa_newton, &vpa_dat, &vpa_dat);
	vpa_dat.x.Display("Variables");
	
	for (int i=0; i<vpa_dat.m; i++)
	{
		vpa_dat.c.edit(i, 0, vpa_dat.x(i+2,0));
	}
	//vpa_dat.c.edit(0, 0, 1);
	//vpa_dat.c.edit(1, 0, -1.35);
	//vpa_dat.c.edit(2, 0, 0.67);
	
	double x = 0.0;
	double dx = vpa_dat.L/20.0;
	Matrix<double> ux(21,1);
	for (int i=0; i<21; i++)
	{
		x = ((double)i*dx);
		ux.edit(i, 0, Eval_ApproximatePolySolution(vpa_dat.c,x));
	}
	ux.Display("u_x");
	
	std::cout << "Test was a HUGE SUCCESS!!! ^_^ \n\n";
	
	// --------------------------------------------- END VPA Example -----------------------------------------------
	
	// ----------------------------- Example of QR Solve ------------------------------
	QR_DATA qr_dat;
	QR_EX1 ex1;
	ex1.M.set_size(10, 10);
	ex1.M.tridiagonalFill(-1, 2, -1, false);
	ex1.M.Display("M");
	Matrix<double> bq;
	Matrix<double> xq;
	bq.set_size(10, 1);
	xq.set_size(10, 1);
	bq.edit(0, 0, 1.0);
	bq.Display("b");
	xq.ladshawSolve(ex1.M, bq);
	xq.Display("x_tri");
	
	xq.zeros();
	xq.qrSolve(ex1.M, bq);
	xq.Display("x_macaw");
	
	success = QRsolve(ex1_mult, bq, &qr_dat, (void *)&ex1);
	
	qr_dat.Ro.Display("R");
	qr_dat.x.Display("x_lark");
	
	Matrix<double> A1, b1, x1;
	A1.set_size(7, 7);
	b1.set_size(7, 1);
	x1.set_size(7, 1);
	b1(0,0) = 1.0;
	A1(0,0) = 1; A1(0,1) = 2; A1(0,2) = 5; A1(0,3) = 3; A1(0,4) = 6; A1(0,5) = 1; A1(0,6) = 0;
	A1(1,0) = 2; A1(1,1) = 6; A1(1,2) = 1; A1(1,3) = 0; A1(1,4) = 1; A1(1,5) = 3; A1(1,6) = 5;
	A1(2,0) = 5; A1(2,1) = 2; A1(2,2) = 3; A1(2,3) = 1; A1(2,4) = 1; A1(2,5) = 0; A1(2,6) = 6;
	A1(3,0) = 0; A1(3,1) = 6; A1(3,2) = 1; A1(3,3) = 1; A1(3,4) = 5; A1(3,5) = 3; A1(3,6) = 2;
	A1(4,0) = 3; A1(4,1) = 0; A1(4,2) = 6; A1(4,3) = 5; A1(4,4) = 1; A1(4,5) = 2; A1(4,6) = 1;
	A1(5,0) = 1; A1(5,1) = 1; A1(5,2) = 0; A1(5,3) = 3; A1(5,4) = 5; A1(5,5) = 6; A1(5,6) = 2;
	A1(6,0) = 0; A1(6,1) = 4; A1(6,2) = 5; A1(6,3) = 1; A1(6,4) = 2; A1(6,5) = 3; A1(6,6) = 0;
	
	x1.qrSolve(A1, b1);
	
	A1.Display("A");
	b1.Display("b");
	x1.Display("x");
	
	std::cout << "QR solve norm = " << (b1 - A1*x1).norm() << std::endl << std::endl;
	
	
	// ------------------------------------- END QR Solve Example -----------------------------------------
	
	// ----------------------------- Example of Gauss-Seidel ------------------------------
	Matrix<double> r1, U1, L1, s1;
	r1.set_size(200, 1);
	s1.set_size(200, 1);
	A1.set_size(200, 200);
	U1.set_size(200, 200);
	L1.set_size(200, 200);
	x1.set_size(200, 1);
	b1.set_size(200, 1);
	A1.tridiagonalFill(-1, 1.9999, -1, false);
	U1.tridiagonalFill(0, 0, -1, false);
	L1.tridiagonalFill(-1, 1.9999, 0, false);
	//A1.Display("A");
	x1.zeros();
	b1(0,0) = 1.0;
	r1 = (b1 - A1*x1);
	double norm = (b1 - A1*x1).norm();
	std::cout << "Gauss-Seidel Method" << std::endl;
	std::cout << 0 << "\t" << norm << std::endl;
	for (int i=0; i<10; i++)
	{
		s1.lowerTriangularSolve(A1, r1);
		x1 = s1+x1;
		r1 = (b1 - A1*x1);
		
		norm = (b1 - A1*x1).norm();
		std::cout << i+1 << "\t" << norm << std::endl;
	}
	std::cout << "\n\n";
	
	// ------------------------------------- END Gauss-Seidel Example -----------------------------------------
	
	// ----------------------------- Example 02 of Varitational Polynomial Approximation ------------------------------
	
	std::cout << "Solve {du/dt + v*du/dx = D*d^2u/dx^2} Approximately with a constrained variational method and implicit/CN integration...\n\n";
	
	VPA_Test02_DATA vpa_dat02;
	vpa_dat02.m = 11; //Polynomial order - this may be a bad basis set for advection problems
	vpa_dat02.N = vpa_dat02.m + 2;
	vpa_dat02.L = 1.0;
	vpa_dat02.D = 0.2;
	vpa_dat02.v = 5.0;
	vpa_dat02.dt = 0.01;
	vpa_dat02.uo = 1.0;
	vpa_dat02.beta = 1.0;
	vpa_dat02.cnp1.set_size(vpa_dat02.m, 1);
	vpa_dat02.xnp1.set_size(vpa_dat02.N, 1);
	
	//Initial Conditions
	vpa_dat02.cn.set_size(vpa_dat02.m, 1);
	vpa_dat02.xn.set_size(vpa_dat02.N, 1);
	vpa_dat02.cn.zeros();
	vpa_dat02.cnp1.zeros();
	vpa_dat02.xn.zeros();
	vpa_dat02.xnp1.zeros();
	
	PJFNK_DATA vpa_newton02;
	vpa_newton02.linear_solver = QR;
	vpa_newton02.LineSearch = true;
	vpa_newton02.nl_tol_abs = 1e-4;
	vpa_newton02.nl_tol_rel = 1e-6;
	vpa_newton02.NL_Output = false;
	
	double end_time = 1.0;
	double current_time = 0.0;
	dx = vpa_dat02.L/10.0;
	std::cout << "\t-------------------- x values ------------------------------ \n";
	std::cout << "Time";
	x = 0.0;
	for (int i=0; i<11; i++)
	{
		x = ((double)i*dx);
		std::cout << "\t" << x;
	}
	std::cout << "\n";
	do
	{
		x = 0.0;
		std::cout << current_time;
		for (int i=0; i<11; i++)
		{
			x = ((double)i*dx);
			std::cout  << "\t" << Eval_ApproximatePolySolution(vpa_dat02.cnp1,x);
			//std::cout  << "\t" << Eval_ApproximateTrigSolution(vpa_dat02.cnp1,x);
		}
		std::cout << "\n";
		
		success = pjfnk(Eval_VPA_Test02_Residuals, NULL, vpa_dat02.xnp1, &vpa_newton02, &vpa_dat02, &vpa_dat02);
		//success = pjfnk(Eval_VPA_Test03_Residuals, NULL, vpa_dat02.xnp1, &vpa_newton02, &vpa_dat02, &vpa_dat02);
	
		for (int i=0; i<vpa_dat02.m; i++)
		{
			vpa_dat02.cnp1.edit(i, 0, vpa_dat02.xnp1(i+2,0));
		}
		
		//Reset n = n+1 level
		vpa_dat02.cn = vpa_dat02.cnp1;
		vpa_dat02.xn = vpa_dat02.xnp1;
		vpa_dat02.xnp1.zeros(); //Note: these are zeroed out because the method seems more efficient this way
		current_time += vpa_dat02.dt;
	} while (current_time <= end_time+vpa_dat02.dt);
	
	
	// --------------------------------------------- END VPA Example 02 -----------------------------------------------
	std::cout << "\n\n";
	
	//std::cout << TrigBasis_Integrals(COS, SIN, 2, 1, 0, 1) << std::endl;
	//std::cout << Overlap_Integral_TrigBasis(4, 1, 0, 1) << std::endl;
	//std::cout << Gradient_Integral_TrigBasis(4, 1, 0, 1) << std::endl;
	//std::cout << Laplacian_Integral_TrigBasis(4, 1, 0, 1) << std::endl;
	
	std::cout << "\nEnd SANDBOX\n\n";
	
	return success;
}
