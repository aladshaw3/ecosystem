//  main.cpp
//  AdsorptionToolBox

#include "Trajectory.h"

//Force Calculations
double Magnetic_R(const Matrix<double>& dX, const Matrix<double>& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a)
{
	return -((4.0*3.1415*pow(b,3)*mu_0*chi_p*M*H0)/(3.0*a))*((M/(2.0*H0*pow(dY(i-1,0),5.0)))+(cos(2.0*dX(i-1,0))/pow(dY(i-1,0),3.0)));
}

double Magnetic_T(const Matrix<double>& dX, const Matrix<double>& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a)
{
	return -((4.0*3.1415*pow(b,3.0)*mu_0*chi_p*M*H0)/(3.0*a))*(sin(2.0*dX(i-1,0))/pow(dY(i-1,0),3.0));
}

double Grav_R(const Matrix<double>& dX, int i, double b, double rho_p, double rho_f)
{
	return -(4.0*3.1415*pow(b,3.0)/3.0)*(rho_p-rho_f)*9.81*sin(dX(i-1,0));
}

double Grav_T(const Matrix<double>& dX, int i, double b, double rho_p, double rho_f)
{
	return -(4.0*3.1415*pow(b,3.0)/3.0)*(rho_p-rho_f)*9.81*cos(dX(i-1,0));
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
	POL(7,0) = i;
	POL(8,0) = V_THETA(dX, dY, i, dat->V0, dat->rho_f, dat->a, dat->eta);
	POL(9,0) = Brown_THETA (dat->s_rand, dat->t_rand, dat->sigma_n, dat->sigma_m);

	return 0;
}

double RADIAL_FORCE (const Matrix<double>& POL, double eta, double b, double mp, double t, double a)
{
	return ((POL(0,0)+POL(1,0)+POL(2,0))*((t/(6.0*3.1415*eta*b))-((1-exp(-6.0*3.1415*eta*b*t/mp))/(pow(6.0*3.1415*eta*b,2.0)/mp)))+(POL(3,0)*t)+POL(4,0))/a;
}

double TANGENTIAL_FORCE (const Matrix<double>& POL, const Matrix<double>& dY, double eta, double b, double mp, double t, double a, int i)
{
	return (POL(5,0)+POL(6,0)*((t/(6.0*3.1415*eta*b))-((1-exp(-6.0*3.1415*eta*b*t/mp))/(pow(6.0*3.1415*eta*b,2.0)/mp)))+POL(8,0)*t+POL(9,0))/(a*dY(i-1,0));
}

int CARTESIAN(const Matrix<double>& POL, Matrix<double>& H, const Matrix<double>& dY, double i, const void *data)
{
	TRAJECTORY_DATA *dat = (TRAJECTORY_DATA *) data;
	H(0,0) = RADIAL_FORCE(POL, dat->eta, dat->b, dat->mp, dat->dt, dat-> a);
	H(1,0) = TANGENTIAL_FORCE(POL, dY, dat->eta, dat->b, dat->mp, dat->dt, dat->a, i);

	return 0;
}

int DISPLACEMENT (Matrix<double>& dX, Matrix<double>& dY, const Matrix<double>& H, int i)
{
	dX(i,0) = dX(i-1,0)+H(1,0);
	dY(i,0) = dY(i-1,0)+H(0,0);
	return 0;
}

int LOCATION (const Matrix<double>& dY, const Matrix<double>& dX, Matrix<double>& XX, Matrix<double>& YY, int i)
{
	XX(i,0) = dY(i,0)*cos(dX(i,0));
	YY(i,0) = dY(i,0)*sin(dX(i,0));

	return 0;
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
	
	dat->Rs = 0.0026925;										//Separator radius, m
	dat->L = 0.0611;											//Separator length, m
	dat->porosity = 0.9036;										//Separator porosity
	dat->V_separator=M_PI*pow(dat->Rs,2.0)*dat->L;				//Volume of separator, m^3
	
	//System Parameters
	
	dat->a = 3.3e-5;											//Wire radius,m
	dat->V_wire = (1-dat->porosity)*dat->V_separator;			//Total wire volume, m^3
	dat->L_wire = dat->V_wire/M_PI/pow(dat->a,2.0);				//Total wire length, m
	dat->A_separator = M_PI*pow(dat->Rs,2.0);
	dat->A_wire = dat->V_wire/dat->L_wire;
	dat->B0 = 1.1;												//Applied magnetic induction,T
	dat->H0 = dat->B0/dat->mu_0;
	dat->Ms = 0.6;												//Saturated magnetization,T
	
	//Particle parameters
	
	dat->b = 1.0e-6;											//Particle radius,m
	dat->chi_p = 3.87e-7;										//Volume magnetic susceptibility, dimensionless
	dat->rho_p = 8700.0;										//Particle density, Kg/m3
	
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
	dat-> sigma_z = (dat->q_bar/pow(dat->beta,3.0))*(2.0*dat->beta*dat->dt-3.0-4.0*exp(-dat->beta*dat->dt)-exp(-2.0*dat->beta*dat->dt));
	dat-> sigma_n = (dat->sigma_vz/pow(dat->sigma_v,0.5));
	dat-> sigma_m = (dat->sigma_z-pow(dat->sigma_vz,2.0))/dat->sigma_v;

	return 0;
}

//Generates Random Numbers
int Number_Generator(Matrix<double>& Temporary, TRAJECTORY_DATA *dat)
{
	unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::normal_distribution<double> distribution(0.0,1.0);
	dat-> n_rand = distribution(generator);
	dat-> m_rand = distribution(generator);
	dat-> s_rand = distribution(generator);
	dat-> t_rand = distribution(generator);

	Temporary(0,0) = dat->n_rand;
	Temporary(1,0) = dat->m_rand;
	Temporary(2,0) = dat->s_rand;
	Temporary(3,0) = dat->t_rand;

	return 0;
}

int Run_Trajectory()
{
	int success = 0;
	
	TRAJECTORY_DATA dat;
	
	success = Trajectory_SetupConstants(&dat);
	
	dat.POL.set_size(10,1);
	dat.H.set_size(2,1);
	dat.dX.set_size(300,1);
	dat.dY.set_size(300,1);
	dat.X.set_size(300, 1);
	dat.Y.set_size(300, 1);
	dat.Temporary.set_size(4,1);
	
	bool Hit = false;
	
	dat.dY(0,0) = 20.0;
	dat.dX(0,0) = 89.99/180.0*M_PI;			//FIXED negative starting values for x
	
	for (int i = 1; i<299; i++)
	{
		int Random;
		Random = Number_Generator(dat.Temporary,&dat);

		dat.Temporary.Display("Random Numbers");

		int Force_Polar;
		Force_Polar = POLAR(dat.POL,dat.dX,dat.dY,(void *)&dat,i);
		
		int Force_Balance;
		Force_Balance = CARTESIAN(dat.POL,dat.H,dat.dY,i,(void *)&dat);

		int Next;
		Next = DISPLACEMENT(dat.dX,dat.dY,dat.H,i);
		
		if (dat.dY(i,0) <= (1.0+(1.85/33.0)) && dat.H(0,0) < 0.0)
		{
			Hit = true;
			break;
		}
		else
		{
			Hit = false;
		}
		
	}
	
	for (int i=0; i<300; i++)
	{
		
		int Final;
		Final = LOCATION (dat.dY,dat.dX,dat.X,dat.Y,i);
		
	}
	
	//dat.X.Display("Distance X");
	//dat.Y.Display("Diatacne Y");

	return success;
}

