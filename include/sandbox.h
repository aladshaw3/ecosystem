//----------------------------------------
//  Created by Austin Ladshaw on 04/11/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef SANDBOX_HPP_
#define SANDBOX_HPP_

#include "flock.h"
#include "school.h"

typedef struct
{
	int N = 4;
	
	const double logKw = -14.0;
	const double logKa1 = -6.35;
	const double logKa2 = -10.33;
	double CT = 0.1786;
	double NaT = 0.1786;
	
	std::vector<Molecule> x;
	
	Matrix<double> Jacobian;
	Matrix<double> NumJac;
	Matrix<double> logC;
	Matrix<double> C;
	
}Speciation_Test01_Data;

typedef struct
{
	double eps = sqrt(DBL_EPSILON);
	Matrix<double> Fx;
	Matrix<double> Fxp;
	Matrix<double> dxj;
	
}NUM_JAC_DATA;

int Speciation_Test01_Function(const Matrix<double> &x, Matrix<double> &F, const void *res_data);

int Speciation_Test01_Jacobian(const Matrix<double> &x, Matrix<double> &J, const void *precon_data);

int Speciation_Test01_Guess(const void *user_data);

int Speciation_Test01_MatVec(const Matrix<double> &x, Matrix<double> &Ax, const void *matvec_data);

int NumericalJacobian( int (*Func) (const Matrix<double> &x, Matrix<double> &F, const void *user_data),
					   const Matrix<double> &x, Matrix<double> &J, int Nx, int Nf, NUM_JAC_DATA *jac_dat,
					   const void *user_data);

int RUN_SANDBOX();

#endif
