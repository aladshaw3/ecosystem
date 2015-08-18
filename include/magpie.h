//----------------------------------------
//  Created by Austin Ladshaw on 12/17/13
//  Copyright (c) 2013
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

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

#ifndef DBL_EPSILON									//Used for gradient estimation
#define DBL_EPSILON    2.2204460492503131e-016 		//machine dependent (see float.h)
#endif

#ifndef	Z					//Define the coordination number
#define Z 10.0
#endif

#ifndef	A				//Define corresponding van der Waals standard area
#define A 3.13E+09				//sq.cm / mole
#endif

#ifndef	V				//Define corresponding van der Waals standard volume
#define V 18.92					//cu.cm / mole
#endif

#ifndef	Po				//Standard State Pressure
#define Po 100.0					//Units: kPa
#endif

#ifndef	R					//Gas Constant
#define R 8.3144621				//Units: J/(K*mol) = kB * Na
#endif

#ifndef	Na				//Avagadro's Number
#define Na 6.0221413E+23		//Units: molecules/mol
#endif

#ifndef	kB				//Boltzmann's Constant
#define kB 1.3806488E-23		//Units: J/K
#endif

//This macro replaces all instances of shapeFactor(#) with the following single line calculation
#ifndef shapeFactor
#define shapeFactor(v_i) ( ( (Z - 2) * v_i ) / ( Z * V ) ) + ( 2 / Z )
#endif

//This macro calculates the natural log of the dimensionless isotherm parameter
#ifndef lnKo
#define lnKo(H,S,T)	-( H / ( R * T ) ) + ( S / R )
#endif

//This macro calculates the Henry's Coefficient for the ith component
#ifndef He
#define He(qm,K1,m) ( qm * K1 ) / ( m * Po )
#endif

//GSTA Data Structure
typedef struct
{
	double qmax;						//Theoretical maximum capacity of adsorbate-adsorbent pair (mol/kg)
	int m;								//Number of parameters in the GSTA isotherm
	std::vector<double> dHo;			//Enthalpy (J/mol)
	std::vector<double> dSo;			//Entropy (J/(K*mol))
}GSTA_DATA;

//mSPD Data Structure
typedef struct
{
	double s;							//Area shape factor
	double v;							//van der Waals Volume (cm^3/mol)
	double eMax;							//Maximum lateral interaction energy (J/mol)
	std::vector<double> eta;			//Binary interaction parameter matrix (i,j)
	double gama;						//Activity calculated from mSPD
}mSPD_DATA;

//GPAST Data Structure
typedef struct
{
	double x;							//Adsorbed mole fraction
	double y;							//Gas phase mole fraction
	double He;							//Henry's Coefficient (mol/kg/kPa)
	double q;							//Amount adsorbed for each component (mol/kg)
	std::vector<double> gama_inf;		//Infinite dilution activities
	double qo;							//Pure component capacities
	double PIo;							//Pure component spreading pressures
	std::vector<double> po;				//Pure component reference state pressures
	double poi;							//Reference state pressures solved for using Recover eval GPAST
	bool present;						//If true, then the component is present; if false, then the component is not present
}GPAST_DATA;

//System Data Structure
typedef struct
{
	double T;							//System Temperature (K)
	double PT;							//Total Pressure (kPa)
	double qT;							//Total Amount adsorbed (mol/kg)
	double PI;							//Total Spreading Pressure (mol/kg)
	double pi;							//Spreading pressure (J/m^2)
	double As;							//Specific surface area of adsorbent (m^2/kg)
	int N;								//Total Number of Components (adsorbable components only)
	int I,J,K;							//Special indices used to keep track of sub-systems
	unsigned long int total_eval;		//Counter to keep track of total number of non-linear steps
	double avg_norm;					//Used to store all norms from evaluations then average at end of run
	double max_norm;					//Used to store the maximum e.norm calculated from non-linear iterations
	int Sys;							//Number of sub-systems to solve
	int Par;							//Number of binary parameters to solve for
	bool Recover;						//If Recover == false, standard GPAST using y's as knowns
	bool Carrier;						//If there is an inert carrier gas, Carrier == true
	bool Ideal;							//If the behavior of the system is determined to be ideal, then Ideal == true
	bool Output;						//Boolean to suppress output if desired (true = display, false = no display
}SYSTEM_DATA;

//MAGPIE Data Structure - Abstraction and Inheritance
typedef struct
{
	std::vector<GSTA_DATA> gsta_dat;
	std::vector<mSPD_DATA> mspd_dat;
	std::vector<GPAST_DATA> gpast_dat;
	SYSTEM_DATA sys_dat;
}MAGPIE_DATA;

double qo(double po, const void *data, int i);

double dq_dp(double p, const void *data, int i);

double q_p(double p, const void *data, int i);

double PI(double po, const void *data, int i);

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
