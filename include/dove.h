/*!
 *  \file dove.h
 *	\brief Dynamic Ode solver with Various Established methods
 *	\details This file creates objects and subroutines for solving systems of Ordinary Differential
 *			Equations using various established methods. The basic idea is that a user will create
 *			a function to calculate all the right-hand sides of a system of ODEs, then pass that
 *			function to the DOVE routine, which will seek a numerical solution to that system.
 *
 *			Methods for Integration
 *			-----------------------
 *			(None available - still under construction)
 *
 *	\note This kernel is still under construction.
 *
 *  \author Austin Ladshaw
 *	\date 09/25/2017
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2017, all
 *             rights reserved.
 */

#include "macaw.h"
#include "lark.h"
#include "yaml_wrapper.h"

#ifndef DOVE_HPP_
#define DOVE_HPP_

/// Enumeration for the list of valid time integration types
/** The only types that have been defined are for Implicit and Explicit methods.
	Sub-type enumeration is used to denote the specific methods.*/
typedef enum {IMPLICIT, EXPLICIT} integrate_type;

/// Enumeration for the list of valid time integration subtypes
/** Theses subtypes define the specific scheme to be used. The table below gives
	a brief description of each.
 
	\param BE Backwards-Euler: Standard implicit method.
	\param FE Forwards-Euler: Standard explicit method.
	\param CN Crank-Nicholson: Time averaged, 2nd order implicit scheme.
	\param BDF2 Backwards-Differentiation-Formula-2: 2nd order implicit method.
	\param RK4 Runge-Kutta-4: 4th order explicit method. */
typedef enum {BE, FE, CN, BDF2, RK4} integrate_subtype;

/// Enumeration for the list of valid time stepper types
/** Type of time stepper to be used by Dove.
 
	\param CONSTANT time stepper will use a constant dt value for all time steps.
	\param ADAPTIVE time stepper will adjust the time step according to simulation success.*/
typedef enum {CONSTANT, ADAPTIVE} timestep_type;

/// Dynamic ODE-solver with Various Established methods (DOVE) object
/** This class structure creates a C++ object that can be used to solve coupled systems of 
	Ordinary Differential Equations. A user will interface with this object by creating functions
	that evaluate the right-hand side of an ODE based on the given variable set. Those functions
	will collectively create the system to solve numerically using either explicit or implicit 
	methods. The choice of methods can be set by the user, or the object will default to the
	Backwards-Euler implicit method for stability. 
 
	User functions for the right-hand side are written as...
 
	du_i/dt = user_function_i(const Matrix<double> &u, const void *data_struct)
 
	In some cases, there is a need to include a time coefficient on the left-hand side of the rate
	expression. For those cases, the user may also provide a time coefficient function...
 
	user_time_coeff_i(const Matrix<double> &u, const void *data_struct) * du_i/dt = user_function_i(...)
 
	For most implicit problems, the ODE system must be solved iteratively using a Newton-style method. In 
	these cases, the user may also provide functions for Jacobian matrix elements...
	
	user_jacobi_element_i_j(const Matrix<double> &u, const void *data_struct)
 
	All of these above functions are to be put into Matrices inside of the Dove class object so that Dove
	will call those functions when it needs to be called. Data structures for all function calls are optional
	and are to be defined by the user to contain whatever parameter information is needed for their particular
	problem. */
class Dove
{
public:
	Dove();												///< Default constructor
	~Dove();											///< Default destructor
	
	void set_numfunc(int i);							///< Set the number of functions to solve and reserve necessary space
	void set_timestep(double d);						///< Set the value of the time step
	void set_timestepmin(double dmin);					///< Set the value of the minimum time step
	void set_timestepmax(double dmax);					///< Set the value of the maximum time step
	void set_endtime(double e);							///< Set the value of the end time
	void set_integrationtype(integrate_subtype type);	///< Set the type of integration scheme to use
	void set_timestepper(timestep_type type);			///< Set the time stepper scheme type
	void set_outputfile(FILE *file);					///< Set the output file for simulation results
	void set_userdata(const void *data);				///< Set the user defined data structure
	void set_initialcondition(int i, double ic);		///< Set the initial condition of variable i to value ic
	
	void registerFunction(int i, double (*func) (const Matrix<double> &u, const void *data) );		///< Register the ith user function
	void registerCoeff(int i, double (*coeff) (const Matrix<double> &u, const void *data) );		///< Register the ith time coeff function
	void registerJacobi(int i, int j, double (*jac) (const Matrix<double> &u, const void *data) );	///< Register the i-jth element of jacobian
	
	void print_header();								///< Function to print out a header to output file
	void print_newresult();								///< Function to print out the new result of n+1 time level
	
	Matrix<double>& getCurrentU();					///< Return reference to the n level solution
	Matrix<double>& getOldU();						///< Return reference to the n-1 level solution
	Matrix<double>& getNewU();						///< Return reference to the n+1 level solution
	const void *getUserData();						///< Return pointer to user data
	int getNumFunc();								///< Return the number of functions
	double getTimeStep();							///< Return the current time step
	double getTimeStepOld();						///< Return the old time step
	double getEndTime();							///< Return value of end time
	double getCurrentTime();						///< Return the value of current time
	double getMinTimeStep();						///< Return the value of the minimum time step
	double getMaxTimeStep();						///< Return the value of the maximum time step
	bool hasConverged();							///< Returns state of convergence
	
	double Eval_Func(int i, const Matrix<double>& u);	///< Evaluate user function i at given u matrix
	double Eval_Coeff(int i, const Matrix<double>& u);	///< Evaluate user time coefficient function i at given u matrix
	double Eval_Jacobi(int i, int j, const Matrix<double>& u);	///< Evaluate user jacobian function for (i,j) at given u matrix
	
	int solve_timestep();							///< Function to solve a single time step
	
	void update_states();							///< Function to update the stateful information
	
protected:
	Matrix<double> un;								///< Matrix for nth level solution vector
	Matrix<double> unp1;							///< Matrix for n+1 level solution vector
	Matrix<double> unm1;							///< Matrix for n-1 level solution vector
	double dt;										///< Time step between n and n+1 time levels
	double dt_old;									///< Time step between n and n-1 time levels
	double time_end;								///< Time on which to end the ODE simulations
	double time;									///< Value of current time
	double dtmin;									///< Minimum allowable time step
	double dtmax;									///< Maximum allowable time step
	integrate_type int_type;						///< Type of time integration to use
	integrate_subtype int_sub;						///< Subtype of time integration scheme to use
	timestep_type timestepper;						///< Type of time stepper to be used
	FILE *Output;									///< File to where simulation results will be place
	int num_func;									///< Number of functions in the system of ODEs
	bool Converged;									///< Boolean to hold information on whether or not last step converged
	
	Matrix<double (*) (const Matrix<double> &u, const void *data)> user_func;	///< Matrix object for user defined rate functions
	Matrix<double (*) (const Matrix<double> &u, const void *data)> user_coeff;	///< Matrix object for user defined time coefficients (optional)
	Matrix<double (*) (const Matrix<double> &u, const void *data)> user_jacobi;	///< Matrix object for user defined Jacobian elements (optional)
	const void *user_data;														///< Pointer for user defined data structure
	
	PJFNK_DATA newton_dat;							///< Data structure for the PJFNK iterative method
	/// Function pointer for the residual function of DOVE
	int (*residual) (const Matrix<double> &x, Matrix<double> &F, const void *res_data);
	/// Function pointer for the preconditioning operation of DOVE
	int (*precon) (const Matrix<double> &r, Matrix<double> &p, const void *precon_data);
	
private:
	
};

/// Residual function for implicit-BE method
/** This function will be passed to PJFNK as the residual function for the Dove object. In this function,
 DOVE will call the user defined rate functions to create a vector of residuals at the current iterate. That
 information will be passed into the pjfnk function (see lark.h) to iteratively solve the system of equations
 at a single time step.
 
 Res[i] = Rnp1[i]*unp1[i] - Rn[i]*un[i] - dt*func[i](unp1)   */
int residual_BE(const Matrix<double> &u, Matrix<double> &Res, const void *data);

/// Default time coefficient function
double default_coeff(const Matrix<double> &u, const void *data);

/// Default Jacobian element function
double default_jacobi(const Matrix<double> &u, const void *data);

/// Test function for DOVE kernel
/** This function sets up and solves a test problem for DOVE. It is callable from the UI. */
int DOVE_TESTS();

#endif /* DOVE_HPP_ */
