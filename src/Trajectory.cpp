/*!
 *  \file Trajectory.cpp Trajectory.h
 *	\brief Single Particle Trajectory Analysis for Magnetic Filtration
 *  \author Alex Wiechert
 *	\date 08/25/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Alex Wiechert for PhD research in the area
 *             of environmental surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "Trajectory.h"

//Force Calculations
double Magnetic_R(const Matrix<double>& dX, const Matrix<double>& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a)
{
	return -((4.0*M_PI*pow(b,3)*mu_0*chi_p*M*H0)/(3.0*a))*((M/(2.0*H0*pow(dY(i-1,0),5.0)))+(cos(2.0*dX(i-1,0))/pow(dY(i-1,0),3.0)));
}

double Magnetic_T(const Matrix<double>& dX, const Matrix<double>& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a)
{
	return -((4.0*M_PI*pow(b,3.0)*mu_0*chi_p*M*H0)/(3.0*a))*(sin(2.0*dX(i-1,0))/pow(dY(i-1,0),3.0));
}

double Grav_R(const Matrix<double>& dX, int i, double b, double rho_p, double rho_f)
{
	return -(4.0*M_PI*pow(b,3.0)/3.0)*(rho_p-rho_f)*9.81*sin(dX(i-1,0));
}

double Grav_T(const Matrix<double>& dX, int i, double b, double rho_p, double rho_f)
{
	return -(4.0*M_PI*pow(b,3.0)/3.0)*(rho_p-rho_f)*9.81*cos(dX(i-1,0));
}

double Van_R(const Matrix<double>& dX, const Matrix<double>& dY, int i, double Hamaker, double b, double a)
{
	return -2.0*Hamaker/(3.0*b*pow(((dY(i-1,0)*a-a-b)/b)+2.0,2.0)*pow((dY(i-1,0)*a-a-b)/b,2.0));
}

double V_RAD (const Matrix<double>& dX, const Matrix<double>& dY, int i, double V0, double rho_f, double a, double eta)
{
	return -V0*(log(dY(i-1,0))-0.5*(1-pow(1/dY(i-1,0),2.0)))*sin(dX(i-1,0))/(2.002-log(2.0*V0*rho_f*a/eta));
}

double V_THETA (const Matrix<double>& dX, const Matrix<double>& dY, int i, double V0, double rho_f, double a, double eta)
{
	return -V0*(log(dY(i-1,0))+0.5*(1-pow(1/dY(i-1,0),2.0)))*cos(dX(i-1,0))/(2.002-log(2.0*V0*rho_f*a/eta));
}

double Brown_RAD (double n_rand, double m_rand, double sigma_n, double sigma_m)
{
	return (n_rand*sigma_n)+(m_rand*sigma_m);
}

double Brown_THETA (double s_rand, double t_rand, double sigma_n, double sigma_m)
{
	return (s_rand*sigma_n)+(t_rand*sigma_m);
}

int POLAR(Matrix<double>& POL, const Matrix<double>& dX, const Matrix<double>& dY, const void *data, int i)
{
	TRAJECTORY_DATA *dat = (TRAJECTORY_DATA *) data;
	POL(0,0) = Magnetic_R(dX, dY, i, dat->b, dat->mu_0, dat->chi_p, dat-> M, dat->H0, dat->a);
	POL(1,0) = Grav_R(dX, i, dat->b, dat->rho_p, dat->rho_f);
	POL(2,0) = Van_R(dX, dY, i, dat->Hamaker, dat->b, dat->a);
	POL(3,0) = V_RAD(dX, dY, i, dat->V0, dat->rho_f, dat-> a, dat->eta);
	POL(4,0) = Brown_RAD (dat->n_rand, dat->m_rand, dat->sigma_n, dat->sigma_m);

	POL(5,0) = Magnetic_T(dX, dY, i, dat->b, dat->mu_0, dat->chi_p, dat-> M, dat->H0, dat->a);
	POL(6,0) = Grav_T(dX, i, dat->b, dat->rho_p, dat->rho_f);
	POL(7,0) = 0.0;//i;
	POL(8,0) = V_THETA(dX, dY, i, dat->V0, dat->rho_f, dat->a, dat->eta);
	POL(9,0) = Brown_THETA (dat->s_rand, dat->t_rand, dat->sigma_n, dat->sigma_m);

	return 0;
}

double In_PVel_Rad (const Matrix<double>& POL)
{
	return POL(3,0);
}

double In_PVel_Theta (const Matrix<double>& POL)
{
	return POL(8,0);
}

int In_P_Velocity (const Matrix<double>& POL, Matrix<double>& Vr, Matrix<double>& Vt)
{
	Vr(0,0) = In_PVel_Rad(POL);
	Vt(0,0) = In_PVel_Theta(POL);
	return 0;
}

double PVel_Rad (const Matrix<double>& POL, const Matrix<double>& Vr,int i, double mp, double beta, double dt, double sigma_v, double n_rand)
{
	return (Vr(i-2,0)*exp(-beta*dt))+(POL(3,0)+(POL(0,0)+POL(1,0)+POL(2,0))/(mp*beta))*(1-exp(-beta*dt))+(sigma_v*n_rand);
}

double PVel_Theta (const Matrix<double>& POL, const Matrix<double>& Vt,int i, double mp, double beta, double dt, double sigma_v, double s_rand)
{
	return (Vt(i-2,0)*exp(-beta*dt))+(POL(8,0)+((POL(5,0)+POL(6,0)))/(mp*beta))*(1-exp(-beta*dt))+(sigma_v*s_rand);
}

int P_Velocity (const Matrix<double>& POL, Matrix<double>& Vr, Matrix<double>& Vt, int i, const void *data)
{
	TRAJECTORY_DATA *dat = (TRAJECTORY_DATA *) data;
	Vr(i-1,0) = PVel_Rad(POL, Vr, i, dat->mp, dat-> beta, dat->dt, dat-> sigma_v, dat-> n_rand);
	Vt(i-1,0) = PVel_Theta(POL, Vt, i, dat->mp, dat-> beta, dat->dt, dat-> sigma_v, dat-> n_rand);
	return 0;
}

double RADIAL_FORCE (const Matrix<double>& POL, const Matrix<double>& Vr, int i, double beta, double mp, double dt, double a)
{
	return (((Vr(i-1,0)/beta)*(1-exp(-beta*dt)))+((POL(3,0)+(POL(0,0)+POL(1,0)+POL(2,0))/(mp*beta))*(dt-(1/beta)*(1-exp(-beta*dt))))+POL(4,0))/a;
}

double TANGENTIAL_FORCE (const Matrix<double>& POL, const Matrix<double>& Vt, const Matrix<double>& dY, int i, double beta, double mp, double dt, double a)
{
	return (((Vt(i-1,0)/beta)*(1-exp(-beta*dt)))+((POL(8,0)+(POL(5,0)+POL(6,0))/(mp*beta))*(dt-(1/beta)*(1-exp(-beta*dt))))+POL(9,0))/(a*dY(i-1,0));
}

double Capture_Force (const Matrix<double>& POL, const Matrix<double>& Vr, int i, double beta, double mp, double dt, double a)
{
	return (((Vr(i-1,0)/beta)*(1-exp(-beta*dt)))+((POL(3,0)+(POL(0,0)+POL(1,0)+POL(2,0))/(mp*beta))*(dt-(1/beta)*(1-exp(-beta*dt)))))/a;
}

int CARTESIAN(const Matrix<double>& POL, const Matrix<double>& Vr, const Matrix<double>& Vt, Matrix<double>& H, const Matrix<double>& dY, int i, const void *data)
{
	TRAJECTORY_DATA *dat = (TRAJECTORY_DATA *) data;
	H(0,0) = RADIAL_FORCE(POL, Vr, i, dat->beta, dat->mp, dat->dt, dat-> a);
	H(1,0) = TANGENTIAL_FORCE(POL, Vt, dY, i, dat->beta, dat->mp, dat->dt, dat->a);
	H(2,0) = Capture_Force(POL, Vr, i, dat->beta, dat->mp, dat->dt, dat-> a);

	return 0;
}

int DISPLACEMENT (Matrix<double>& dX, Matrix<double>& dY, const Matrix<double>& H, int i)
{
	dX(i,0) = dX(i-1,0)+H(1,0);
	dY(i,0) = dY(i-1,0)+H(0,0);
	return 0;
}


int LOCATION (const Matrix<double>& dY, const Matrix<double>& dX, Matrix<double>& X, Matrix<double>& Y, int i)
{
	X(i,0) = dY(i,0)*cos(dX(i,0));
	Y(i,0) = dY(i,0)*sin(dX(i,0));

	return 0;
}

double Removal_Efficiency (double Sum_Cap, const void *data)
{
	TRAJECTORY_DATA *dat = (TRAJECTORY_DATA *) data;
	return 1.0-exp((-0.12*(1-dat->porosity)*Sum_Cap*dat->L)/(58630.0*(500.0/500.0)*(41.0/41.0)*M_PI*dat->a));
}


//Set values of constants and other variable declared in the data structure
int Trajectory_SetupConstants(TRAJECTORY_DATA *dat)
{
	//Constants
	dat->mu_0 = 4.0*M_PI*pow(10.0, -7.0);						//Permeability of free space, H/m
	dat->rho_f = 1000.0;										//Fluid density, Kg/m3
	dat->eta = 0.001;											//Dynamic viscosity, Kg/m-s
	dat->Hamaker = 1.3e-21;										//Hamaker constant for ferric oxide particle
	
	//Missing some parameters here...
	
	//Separator Parameters
	
	dat->Rs = 0.002743;										//Separator radius, m
	dat->L = 0.0644;											//Separator length, m
	dat->porosity = 0.943;										//Separator porosity
	dat->V_separator=M_PI*pow(dat->Rs,2.0)*dat->L;				//Volume of separator, m^3
	
	//System Parameters
	
	dat->a = 3.9e-5;											//Wire radius,m
	dat->V_wire = (1-dat->porosity)*dat->V_separator;			//Total wire volume, m^3
	dat->L_wire = dat->V_wire/M_PI/pow(dat->a,2.0);				//Total wire length, m
	dat->A_separator = M_PI*pow(dat->Rs,2.0);
	dat->A_wire = dat->V_wire/dat->L_wire;
	dat->B0 = 1.1;												//Applied magnetic induction,T
	dat->H0 = dat->B0/dat->mu_0;
	dat->Ms = 0.6;												//Saturated magnetization,T
	
	//Particle parameters
	
	dat->b = 0.225e-6;											//Particle radius,m
	dat->chi_p = 0.00048;										//Volume magnetic susceptibility, dimensionless
	dat->rho_p = 5240.0;										//Particle density, Kg/m3
	
	//Some more missing info
	
	//Model parameters
	dat->Q_in = 9.5e-7;											//Volumetric fluid flow, m^3/sec
	dat->V0 = dat->Q_in/(dat->A_separator - dat->A_wire);		//Superficial velocity,m/s
	dat->Y_initial = 20.0;										//initial position in Y-axis, Y=y/a, dimensionless
	dat->dt = 2.0*dat->Y_initial*dat->a/(dat->V0*100.0);		//Time step
	
	//Breakthrough parameters
	dat->M = 2.0*dat->Ms/dat->mu_0;
	dat->mp = 4.0*3.1415*pow(dat->b,3.0)*dat->rho_p/3.0;		//Particle Mass, Kg
	
	//Brownian Variables
	dat-> beta = 6.0*M_PI*dat->eta*dat->b/(dat->mp);			//Friction coefficient per unit mass
	dat-> q_bar = dat->beta*dat->Temp*dat->k/dat->mp;
	dat-> sigma_v = (dat->q_bar/dat->beta)*(1-exp(-2.0*dat->beta*dat->dt)); 		//Take square root to get actual value for sigma v and z
	dat-> sigma_vz = (dat->q_bar/pow(dat->beta,2.0))*pow(1-exp(-dat->beta*dat->dt),2.0);
	dat-> sigma_z = (dat->q_bar/pow(dat->beta,3.0))*(2.0*dat->beta*dat->dt-3.0+4.0*exp(-dat->beta*dat->dt)-exp(-2.0*dat->beta*dat->dt));
	dat-> sigma_n = (dat->sigma_vz/pow(dat->sigma_v,0.5));
	dat-> sigma_m = pow((dat->sigma_z-(pow(dat->sigma_vz,2.0)/dat->sigma_v)),0.5);

	return 0;
}

//Generates Random Numbers
int Number_Generator(TRAJECTORY_DATA *dat)
{
	std::random_device generator;
	std::normal_distribution<double> distribution(0.0,1.0);
	dat-> n_rand = distribution(generator);
	dat-> m_rand = distribution(generator);
	dat-> s_rand = distribution(generator);
	dat-> t_rand = distribution(generator);

	return 0;
}

int Run_Trajectory()
{
	int success = 0;
	
	TRAJECTORY_DATA dat;
	
	success = Trajectory_SetupConstants(&dat);
	
	dat.POL.set_size(10,1);
	dat.Vr.set_size(400,1);
	dat.Vt.set_size(400,1);
	dat.H.set_size(3,1);
	dat.dX.set_size(400,1);
	dat.dY.set_size(400,1);
	dat.X.set_size(400, 1);
	dat.Y.set_size(400, 1);
	dat.Cap.set_size(500,41);
	
	bool Hit = false;
	for (int q = 0; q<41; q++)
	{
			dat.Y(0,0) = 20.0;
			dat.dX(0,0) = (89.94+q*0.003)/180.0*M_PI;				//FIXED negative starting values for x
			dat.dY(0,0) = dat.Y(0,0)/sin(dat.dX(0,0));

	for (int j = 0; j<500; j++)
	{
		for (int i = 1; i<399; i++)
		{
			int Random;
			Random = Number_Generator(&dat);

			int Force_Polar;
			Force_Polar = POLAR(dat.POL,dat.dX,dat.dY,(void *)&dat,i);

			if (i == 1)
			{
			int Initial_Particle_Velocity;
			Initial_Particle_Velocity = In_P_Velocity (dat.POL,dat.Vr,dat.Vt);
			}
			else
			{
			int Particle_Velocity;
			Particle_Velocity = P_Velocity(dat.POL,dat.Vr,dat.Vt,i,(void *)&dat);
			}

			int Force_Balance;
			Force_Balance = CARTESIAN(dat.POL,dat.Vr,dat.Vt,dat.H,dat.dY,i,(void *)&dat);

			int Next;
			Next = DISPLACEMENT(dat.dX,dat.dY,dat.H,i);
		
			if (dat.dY(i,0) <= (1.0+(0.225/39.0)) && dat.H(2,0) < 0.0)
			{
				Hit = true;
				break;
			}
			else
			{
				Hit = false;
			}

		}

		/*for (int i=0; i<300; i++)
		{
		
			int Final;
			Final = LOCATION (dat.dY,dat.dX,dat.X,dat.Y,i);
		
		}
		dat.X.Display("X Profile");
		dat.Y.Display("Y Profile");*/

		if (Hit == true)
		{
		dat.Cap(j,q) = 1;
		}
		else if(Hit == false)
		{
		dat.Cap(j,q) = 0;
		}
	}
	}
	
	dat.Cap.Display("Capture Profile");
	double Sum_Cap = dat.Cap.sum();
	double Output = Removal_Efficiency (Sum_Cap, (void *)&dat);
	std::cout<<"Capture:   "<<dat.Cap.sum()<<"   ";
	std::cout<<"Removal Efficiency:   "<<Output*100.0<<"%   ";
	return success;
}

