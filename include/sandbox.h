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

double Bounded_1D_NormalBasis_Func(int i, double bound, double x);

double ParticleNBox_1D_H_MatrixElements_Bound1D_ONB(int i, int j, double bound, double mass, double h_bar);

double ParticleNBox_1D_PHI_MatrixElements_Bound1D_ONB(int i, int j, double bound);

typedef struct
{
	int m;							//Number of basis functions
	int N;							//Size of H and number of total functions in the residual
	Matrix<double> c;				//Vector of coefficients for wavefunction
	Matrix<double> x;				//Solution Vector for non-linear formulation of problem x(0,0) = Eo and x(i,0) = c(i-1,0) for i>0
	double Eo;						//Approximation to the ground-state energy
	Matrix<double> H;				//Hamiltonian Matrix in the basis set
	Matrix<double> PHI;				//Matrix of overlap of basis set
	double mass;					//Mass of the particle
	double h_bar;					//Reduced Planck's Constant
	double box_size;				//Size of the box ( = 2*bounds )
	
}PIB_1D_DATA;

int Form_H_Matrix_1DPIB_Bound1DBasis(PIB_1D_DATA *dat);

int Form_PHI_Matrix_1DPIB_Bound1DBasis(PIB_1D_DATA *dat);

int Eval_1DPIB_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data);

double Eval_PolyBasisFunc(int i, double x);

double Eval_1stDerivative_PolyBasisFunc(int i, double x);

double Eval_2ndDerivative_PolyBasisFunc(int i, double x);

double Eval_ApproximatePolySolution(Matrix<double> &c, double x);

double Gradient_Integral_PolyBasis(int i, int j, double lower, double upper);

double Laplacian_Integral_PolyBasis(int i, int j, double lower, double upper);

double Overlap_Integral_PolyBasis(int i, int j, double lower, double upper);

typedef struct
{
	int m;						//Size of the basis
	int N;						//Size of the problem
	double D;					//Diffusivity parameter
	double L;					//Length of the domain
	double k;					//Reaction coefficient
	double uo;					//Boundary value
	Matrix<double> c;			//Coefficient Matrix (size = m)
	Matrix<double> x;			//Non-linear variable matrix (size = N) (x(0) = lambda_0 && x(1) = lambda_1)
} VPA_Test_DATA;

int Eval_VPA_Test_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data);

typedef struct
{
	int m;						//Size of the basis
	int N;						//Size of the problem
	double v;					//Velocity parameter
	double L;					//Length of the domain
	double D;					//Diffusivity parameter
	double uo;					//Boundary value
	double dt;					//Size of the time step
	Matrix<double> cnp1;		//Coefficient Matrix for n+1 time level (size = m)
	Matrix<double> cn;			//Coefficient Matrix for n time level (size = m)
	Matrix<double> xnp1;		//Non-linear variable matrix for n+1 level (size = N) (x(0) = lambda_0 && x(1) = lambda_1)
	Matrix<double> xn;			//Non-linear variable matrix for n level (size = N) (x(0) = lambda_0 && x(1) = lambda_1)
} VPA_Test02_DATA;

int Eval_VPA_Test02_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data);

double PolyBasis_2D(int i, int j, double x, double y);

double PolyBasis_2D_dx(int i, int j, double x, double y);

double PolyBasis_2D_dx2(int i, int j, double x, double y);

double PolyBasis_2D_dy(int i, int j, double x, double y);

double PolyBasis_2D_dy2(int i, int j, double x, double y);

double PolyBasis_2D_dxdy(int i, int j, double x, double y);

double PolyBasis_2D_LinearComboAppox(Matrix<double> &c, double x, double y);

double Laplacian_Integral_PolyBasis_2D(int i, int j, int l, int m, double x_low, double x_high, double y_low, double y_high);

double Overlap_Integral_PolyBasis_2D(int i, int j, int l, int m, double x_low, double x_high, double y_low, double y_high);

typedef struct
{
	int bcs;				//number of bcs
	int m;					//polynomial order
	int N;					//problem size
	double D;				//diffusivity
	double k;				//reaction
	double Lx;				//distance in x
	double Ly;				//distance in y
	double uo;				//boundary value at (0,0)
	Matrix<double> c;		//Coefficient matrix
	Matrix<double> x;		//variable matrix
	
} VPA_2D_TEST_DATA;

int Eval_2D_VPA_TEST_Residuals(const Matrix<double> &x, Matrix<double> &F, const void *data);

typedef struct
{
	Matrix<double> M;
	
} QR_EX1;

int ex1_mult(const Matrix<double>& x, Matrix<double> &Ax, const void *data);

/// \endcond

/// Function to run the methods implemented in the Sandbox
/** This function is callable from the UI and is used to observe results
	from the tests of newly developed algorithms. Edit header and source
	files here to test out your own routines or functions. Then you can
	run those functions by rebuilding the Ecosystem executable and running
	the sandbox tests. */
int RUN_SANDBOX();

#endif
