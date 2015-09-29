/*!
 *  \file sandbox.h sandbox.cpp
 *	\brief Coding Test Area
 *	\details This file contains a series of simple tests for routines used in other files
 *			and algorithms. Before any code or methods are used, they are tested here to
 *			make sure that they are useful. The tests in the sandbox are callable from the
 *			UI to make it easier to alter existing sandbox code and run tests on new proposed
 *			methods or algorithms.
 *
 *	\warning Functions and methods in this file are not meant to be used anywhere else.
 *  \author Austin Ladshaw
 *	\date 04/11/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#ifndef SANDBOX_HPP_
#define SANDBOX_HPP_

#include "flock.h"
#include "school.h"

/// \cond
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

/// \endcond

/// Function to run the methods implemented in the Sandbox
/** This function is callable from the UI and is used to observe results
	from the tests of newly developed algorithms. Edit header and source
	files here to test out your own routines or functions. Then you can
	run those functions by rebuilding the Ecosystem executable and running
	the sandbox tests. */
int RUN_SANDBOX();

#endif
