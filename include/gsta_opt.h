//----------------------------------------
//  Created by Austin Ladshaw on 12/17/13
//  Copyright (c) 2013
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef GSTA_OPT_HPP_
#define GSTA_OPT_HPP_

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
#include "error.h"

#ifndef	Po						//Standard State Pressure
#define Po 100.0					//Units: kPa
#endif

#ifndef	R						//Gas Constant
#define R 8.3144621				//Units: J/(K*mol) = kB * Na
#endif

#ifndef	Na						//Avagadro's Number
#define Na 6.0221413E+23		//Units: molecules/mol
#endif

//Data structure used in several functions and classes
typedef struct
{
	int total_eval;												//Keeps track of the total number of function evaluations
	int n_par;													//Number of parameters being optimized for
	double qmax;												//Maximum theoretical adsorption capacity (M/M) (0 if unknown)
 	int iso;													//Keeps isotherm that is currently being optimized
    std::vector<std::vector<double> > Fobj;						//Creates a dynamic array to store all Fobj values
    std::vector<std::vector<double> > q,P;						//Creates a dynamic array for q and P data pairs
    std::vector<std::vector<double> > best_par;					//Used to store the values of the parameters of best fit
    std::vector<std::vector<double> > Kno;						//Dimensionless parameters determined from best_par
    std::vector<std::vector<std::vector<double> > > all_pars;  	//Used to create a ragged array of all parameters
    std::vector<std::vector<double> > norms;					//Used to store the values of all the calculated norms
    std::vector<double> opt_qmax;								//If qmax is unknown, this vector holds it's optimized values
} GSTA_OPT_DATA;

void error();

int roundIt(double d);

int twoFifths(int m);

int orderMag(double x);

int minValue(std::vector<int> array);

int minIndex(std::vector<double> array);

int avgPar(std::vector<int> array);

double avgValue(std::vector<double> array);

double weightedAvg(double *enorm, double *x, int n);

double rSq(double *x, double *y, double slope, double vint, int m_dat);

bool isSmooth(double *par, void *data);

void orthoLinReg(double *x, double *y, double *par, int m_dat, int n_par);

void eduGuess(double *P, double *q, double *par, int k, int m_dat, void *data);

double gstaFunc( double p, const double *K, double qmax, int n_par);

double gstaObjFunc(double *t, double *y, double *par, int m_dat, void *data);

void eval_GSTA(const double *par, int m_dat, const void *data, double *fvec, int *info);

int gsta_optimize(const char* fileName);

#endif /* GSTA_OPT_HPP_ */
