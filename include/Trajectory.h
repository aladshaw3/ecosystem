/*!
 *  \file Trajectory.h Trajectory.cpp
 *	\brief Single Particle Trajectory Analysis for Magnetic Filtration
 *	\details Alex, Please provide details here... and elsewhere in the file.
 *  \author Alex Wiechert
 *	\date 08/25/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Alex Wiechert for PhD research in the area
 *             of environmental surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "macaw.h"
#include <random>
#include <chrono>

#ifndef TRAJECTORY_HPP_
#define TRAJECTORY_HPP_

typedef struct
{
	//Constants
	double mu_0 = 12.57e-7;						///< permeability of free space, H/m
	double rho_f = 1000.0;						///< Fluid density, Kg/m3
	double eta = 0.001;							//Dynamic viscosity, Kg/m-s
	double Hamaker = 1.3e-21;					//Hamaker constant for ferric oxide particles
	double Temp = 298;                          //Temperature, Kelvin
	double k = 1.38e-23;						//Boltzmann constant

	//Separator Parameters

	double Rs;									//separator radius, m
	double L;									//0.05separator length, m
	double porosity;							//0.965separator porosity
	double V_separator;							//volume of separator, m^3

	//System Parameters

	double a;										//Wire radius,m
	double V_wire;									//Total wire volume, m^3
	double L_wire;									//Total wire length, m
	double A_separator;
	double A_wire;
	double B0;										//Applied magnetic induction,T
	double H0;
	double Ms = 0.6;								//Saturated magnetization,T

	//Particle parameters

	double b;										//Particle radius,m
	double chi_p;							//Volume magnetic susceptibility, dimensionless
	double rho_p = 8000.0;							//Particle density, Kg/m3

	//Model parameters
	double Q_in;									//Volumetric fluid flow, m^3/sec
	double V0;										//Superficial velocity,m/s
	double Y_initial = 20.0;						//initial position in Y-axis, Y=y/a, dimensionless
	double dt;										//Time step size

	//Breakthrough parameters
	double M;
	double mp;										//Particle Mass, Kg
	
	//Brownian Variables
	double beta;									//Friction coefficient per unit mass
	double q_bar;
	double sigma_v; 								//Take square root to get actual value for sigma v
	double sigma_vz;
	double sigma_z; 								//Take square root for actual value of sigma z
	double sigma_n;
	double sigma_m;

	//Random Numbers
	double n_rand;
	double m_rand;
	double s_rand;
	double t_rand;

	Matrix<double> POL,H,dX,dY,Vr,Vt;
	Matrix<double> X, Y;
	Matrix<int> Cap;

}TRAJECTORY_DATA;

double Magnetic_R(const Matrix<double>& dX, const Matrix<double>& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a);

double Magnetic_T(const Matrix<double>& dX, const Matrix<double>& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a);

double Grav_R(const Matrix<double>& dX, int i, double b, double rho_p, double rho_f);

double Grav_T(const Matrix<double>& dX, int i, double b, double rho_p, double rho_f);

double Van_R(const Matrix<double>& dX, const Matrix<double>& dY, int i, double Hamaker, double b, double a);

double V_RAD (const Matrix<double>& dX, const Matrix<double>& dY, int i, double V0, double rho_f, double a, double eta);

double V_THETA (const Matrix<double>& dX, const Matrix<double>& dY, int i, double V0, double rho_f, double a, double eta);

double Brown_RAD (double n_rand, double m_rand, double sigma_n, double sigma_m);

double Brown_THETA (double s_rand, double t_rand, double sigma_n, double sigma_m);

int POLAR(Matrix<double>& POL, const Matrix<double>& dX, const Matrix<double>& dY, const void *data, int i);

double In_PVel_Rad (const Matrix<double>& POL);

double In_PVel_Theta (const Matrix<double>& POL);

int In_P_Velocity (const Matrix<double>& POL, Matrix<double>& Vr, Matrix<double>& Vt);

double PVel_Rad (const Matrix<double>& POL, const Matrix<double>& Vr,int i, double mp, double beta, double t, double sigma_v, double rand_n);

double PVel_Theta (const Matrix<double>& POL, const Matrix<double>& Vt,int i, double mp, double beta, double t, double sigma_v, double rand_s);

int P_Velocity (const Matrix<double>& POL, Matrix<double>& Vr, Matrix<double>& Vt, int i, const void *data);

double RADIAL_FORCE (const Matrix<double>& POL, const Matrix<double>& Vr, int i, double beta, double mp, double dt, double a);

double TANGENTIAL_FORCE (const Matrix<double>& POL, const Matrix<double>& Vt, const Matrix<double>& dY, int i, double beta, double mp, double dt, double a);

double Capture_Force (const Matrix<double>& POL, const Matrix<double>& Vr, int i, double beta, double mp, double dt, double a);

int CARTESIAN(const Matrix<double>& POL, const Matrix<double>& Vr, const Matrix<double>& Vt, Matrix<double>& H, const Matrix<double>& dY, int i, const void *data);

int DISPLACEMENT (Matrix<double>& dX, Matrix<double>& dY, const Matrix<double>& H, int i);

int LOCATION (const Matrix<double>& dY, const Matrix<double>& dX, Matrix<double>& X, Matrix<double>& Y, int i);

double Removal_Efficiency (double Sum_Cap, const void *data);

int Trajectory_SetupConstants(TRAJECTORY_DATA *dat);

int Number_Generator(TRAJECTORY_DATA *dat);

/// Run_Trajectory function
/** Function to run the Trajectory project. */
int Run_Trajectory();

#endif
