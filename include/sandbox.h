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

int Speciation_Test01_Function(const Matrix<double> &x, Matrix<double> &F, const void *res_data);

int Speciation_Test01_Jacobian(const Matrix<double> &x, Matrix<double> &J, const void *precon_data);

int Speciation_Test01_Guess(const void *user_data);

int Speciation_Test01_MatVec(const Matrix<double> &x, Matrix<double> &Ax, const void *matvec_data);

int RUN_SANDBOX();

#endif
