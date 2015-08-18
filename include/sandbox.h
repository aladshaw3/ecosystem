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
	
	Matrix Jacobian;
	Matrix NumJac;
	Matrix logC;
	Matrix C;
	
}Speciation_Test01_Data;

typedef struct
{
	double eps = sqrt(DBL_EPSILON);
	Matrix Fx;
	Matrix Fxp;
	Matrix dxj;
	
}NUM_JAC_DATA;

int Speciation_Test01_Function(const Matrix &x, Matrix &F, const void *res_data);

int Speciation_Test01_Jacobian(const Matrix &x, Matrix &J, const void *precon_data);

int Speciation_Test01_Guess(const void *user_data);

int Speciation_Test01_MatVec(const Matrix &x, Matrix &Ax, const void *matvec_data);

int NumericalJacobian( int (*Func) (const Matrix &x, Matrix &F, const void *user_data),
					   const Matrix &x, Matrix &J, int Nx, int Nf, NUM_JAC_DATA *jac_dat,
					   const void *user_data);

int RUN_SANDBOX();

#endif
