/*!
 *  \file magpie.h magpie.cpp
 *	\brief Multicomponent Adsorption Generalized Procedure for Isothermal Equilibria
 *	\details This file contains all functions and routines associated with predicting isothermal
 *			adsorption equilibria from only single component isotherm information. The basis of
 *			the model is the Adsorbed Solution Theory developed by Myers and Prausnitz (1965). 
 *			Added to that base model is a procedure by which we can predict the non-idealities
 *			present at the surface phase by solving a closed system of equations involving
 *			the activity model. 
 *
 *			For more details on this procedure, check out our publication in AIChE where we
 *			give a fully feature explaination of our Generalized Predictive Adsorbed Solution
 *			Theory (GPAST). 
 *
 *			Reference: Ladshaw, A., Yiacoumi, S., and Tsouris, C., "A generalized procedure for
 *				the prediction of multicomponent adsorption equilibria", AIChE J., vol. 61, No. 8,
 *				p. 2600-2610, 2015. 
 *
 *			MAGPIE represents a special case of the more general GPAST procedure, wherin the isotherm
 *			for each species is respresent by the GSTA isotherm (see gsta_opt.h) and the activity 
 *			model for non-ideality at the adsorbent surface is a Modified Spreading Pressure Dependent
 *			(MSPD) model. See the above paper reference for more details.
 *
 *  \author Austin Ladshaw
 *	\date 12/17/2013
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#ifndef MAGPIE_HPP_
#define MAGPIE_HPP_

#include "lmcurve.h"			  //Main include to use the lmfit solver library
#include <stdio.h>				  //Line to allow for printf functions
#include <math.h>                 //Line added to allow usage of the pow (e, x) function
#include <iostream>				  //Line to allow for read/write to the console using cpp functions
#include <fstream>				  //Line to allow for read/write to and from .txt files
#include <stdlib.h>				  //Line need to convert strings to doubles
#include <vector>				  //Line needed to use dynamic arrays called vectors
#include <time.h>				  //Line needed to display program runtime
#include <float.h>				  //Line to allow use of machine constants
#include <string>    			  //Line to allow use of c++ strings
#include "error.h"				  //Line to allow use of custom error reporting

#ifndef DBL_EPSILON
#define DBL_EPSILON    2.2204460492503131e-016 		///< Machine precision value used for approximating gradients
#endif

#ifndef	Z
#define Z 10.0	///< Surface coordination number used in the MSPD activity model
#endif

#ifndef	A
#define A 3.13E+09				///< Corresponding van der Waals standard area for our coordination number (cm^2/mol)
#endif

#ifndef	V
#define V 18.92					///< Corresponding van der Waals standard volume for our coordination number (cm^3/mol)
#endif

#ifndef	Po
#define Po 100.0				///< Standard State Pressure - Units: kPa
#endif

#ifndef	R
#define R 8.3144621				///< Gas Constant - Units: J/(K*mol) = kB * Na
#endif

#ifndef	Na
#define Na 6.0221413E+23		///< Avagadro's Number - Units: molecules/mol
#endif

#ifndef	kB
#define kB 1.3806488E-23		///< Boltzmann's Constant - Units: J/K
#endif

/// This macro replaces all instances of shapeFactor(#) with the following single line calculation
#ifndef shapeFactor
#define shapeFactor(v_i) ( ( (Z - 2) * v_i ) / ( Z * V ) ) + ( 2 / Z )
#endif

/// This macro calculates the natural log of the dimensionless isotherm parameter
#ifndef lnKo
#define lnKo(H,S,T)	-( H / ( R * T ) ) + ( S / R )
#endif

/// This macro calculates the Henry's Coefficient for the ith component
#ifndef He
#define He(qm,K1,m) ( qm * K1 ) / ( m * Po )
#endif

/// GSTA Data Structure
/** C-style object holding all parameter information associated with the Generalized 
	Statistical Thermodynamic Adsorption (GSTA) isotherm model. Each species in the gas phase
	will have one of these objects. */
typedef struct
{
	double qmax;						///< Theoretical maximum capacity of adsorbate-adsorbent pair (mol/kg)
	int m;								///< Number of parameters in the GSTA isotherm
	std::vector<double> dHo;			///< Enthalpies for each site (J/mol)
	std::vector<double> dSo;			///< Entropies for each site (J/(K*mol))
}GSTA_DATA;

/// MSPD Data Structure
/** C-Style object holding all parameter information associated with the Modified
	Spreading Pressure Dependent (SPD) activity model. Each species in the gas phase will have one
	of these objects. */
typedef struct
{
	double s;							///< Area shape factor
	double v;							///< van der Waals Volume (cm^3/mol)
	double eMax;						///< Maximum lateral interaction energy (J/mol)
	std::vector<double> eta;			///< Binary interaction parameter matrix (i,j)
	double gama;						///< Activity coefficient calculated from mSPD
}mSPD_DATA;

/// GPAST Data Structure
/** C-style object holding all parameter information associated with the Generalized
	Predictive Adsorbed Solution Theory (GPAST) system of equations. Each species in the 
	gas phase will have one of these objects. */
typedef struct
{
	double x;							///< Adsorbed mole fraction
	double y;							///< Gas phase mole fraction
	double He;							///< Henry's Coefficient (mol/kg/kPa)
	double q;							///< Amount adsorbed for each component (mol/kg)
	std::vector<double> gama_inf;		///< Infinite dilution activities
	double qo;							///< Pure component capacities (mol/kg)
	double PIo;							///< Pure component spreading pressures (mol/kg)
	std::vector<double> po;				///< Pure component reference state pressures (kPa)
	double poi;							///< Reference state pressures solved for using Recover eval GPAST
	bool present;						///< If true, then the component is present; if false, then the component is not present
}GPAST_DATA;

/// System Data Structure
/** C-style object holding all the data associated with the overall system to be modeled. */
typedef struct
{
	double T;							///< System Temperature (K)
	double PT;							///< Total Pressure (kPa)
	double qT;							///< Total Amount adsorbed (mol/kg)
	double PI;							///< Total Lumped Spreading Pressure (mol/kg)
	double pi;							///< Actual Spreading pressure (J/m^2)
	double As;							///< Specific surface area of adsorbent (m^2/kg)
	int N;								///< Total Number of Components
	int I,J,K;							///< Special indices used to keep track of sub-systems
	unsigned long int total_eval;		///< Counter to keep track of total number of non-linear steps
	double avg_norm;					///< Used to store all norms from evaluations then average at end of run
	double max_norm;					///< Used to store the maximum e.norm calculated from non-linear iterations
	int Sys;							///< Number of sub-systems to solve
	int Par;							///< Number of binary parameters to solve for
	bool Recover;						///< If Recover == false, standard GPAST using y's as knowns
	bool Carrier;						///< If there is an inert carrier gas, Carrier == true
	bool Ideal;							///< If the behavior of the system is determined to be ideal, then Ideal == true
	bool Output;						///< Boolean to suppress output if desired (true = display, false = no display
}SYSTEM_DATA;

/// MAGPIE Data Structure
/** C-style object holding all information necessary to run a MAGPIE simulation. This is the data
	structure that will be used in other sub-routines when a mixed gas adsorption simulation needs
	to be run.*/
typedef struct
{
	std::vector<GSTA_DATA> gsta_dat;
	std::vector<mSPD_DATA> mspd_dat;
	std::vector<GPAST_DATA> gpast_dat;
	SYSTEM_DATA sys_dat;
}MAGPIE_DATA;

/// Function computes the result of the GSTA isotherm for the ith species
/** This function just computes the result of the GSTA isotherm model for the ith
	species given the partial pressure po. 
 
	\param po partial pressure in kPa at which to evaluate the GSTA model
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double qo(double po, const void *data, int i);

/// Function computes the derivative of the GSTA model with respect to partial pressure
/** This function just computes the result of the derivative of GSTA isotherm model 
	for the ith species at the given the partial pressure p.
 
	\param p partial pressure in kPa at which to evaluate the GSTA model
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double dq_dp(double p, const void *data, int i);

/// Function computes the ratio between the adsorbed amount and partial pressure for the GSTA isotherm
/** This function just computes the ratio between the adsorbed amount q (mol/kg) and the 
	partial pressure p (kPa) at the given partial pressure. If p == 0, then this function
	returns the Henry's Law constant for the isotherm of the ith species.
 
	\param p partial pressure in kPa at which to evaluate the GSTA model
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double q_p(double p, const void *data, int i);

/// Function computes the spreading pressure integral of the ith species
/** This function uses an analytical solution to the spreading pressure integral with
	the GSTA isotherm to evaluate and return the value computed by that integral equation.
 
	\param po partial pressure in kPa at which to evaluate the lumped spreading pressure
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double PI(double po, const void *data, int i);

/// Function computes the heat of adsorption based on the ith species GSTA parameters
/** This function computes the isosteric heat of adsorption (J/mol) for the GSTA parameters
	of the ith species.
 
	\param po partial pressure in kPa at which to evaluate the heat of adsorption
	\param data void pointer to the MAGPIE_DATA data structure
	\param i index of the gas species for which the GSTA model is being evaluated*/
double Qst(double po, const void *data, int i);

double eMax(const void *data, int i);

double lnact_mSPD(const double *par, const void *data, int i, volatile double PI);

double grad_mSPD(const double *par, const void *data, int i);

double qT(const double *par, const void *data);

void initialGuess_mSPD(double *par, const void *data);

void eval_po_PI(const double *par, int m_dat, const void *data, double *fvec, int *info);

void eval_po_qo(const double *par, int m_dat, const void *data, double *fvec, int *info);

void eval_po(const double *par, int m_dat, const void *data, double *fvec, int *info);

void eval_eta(const double *par, int m_dat, const void *data, double *fvec, int *info);

void eval_GPAST(const double *par, int m_dat, const void *data, double *fvec, int *info);

int MAGPIE(const void *data);

int MAGPIE_SCENARIOS(const char *inputFileName, const char *sceneFileName);

#endif /* MAGPIE_HPP_ */
