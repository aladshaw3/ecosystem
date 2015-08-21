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

int POLAR(Matrix<double>& POL, const Matrix<double>& dX, const Matrix<double>& dY, const void *data, int i)
{
	TRAJECTORY_DATA *dat = (TRAJECTORY_DATA *) data;
	POL(0,0) = Magnetic_R(dX, dY, i, dat->b, dat->mu_0, dat->chi_p, dat-> M, dat->H0, dat->a);
	POL(1,0) = Grav_R(dX, i, dat->b, dat->rho_p, dat->rho_f);
	POL(2,0) = Van_R(dX, dY, i, dat->Hamaker, dat->b, dat->a);
	POL(3,0) = V_RAD(dX, dY, i, dat->V0, dat->rho_f, dat-> a, dat->eta);

	POL(4,0) = Magnetic_T(dX, dY, i, dat->b, dat->mu_0, dat->chi_p, dat-> M, dat->H0, dat->a);
	POL(5,0) = Grav_T(dX, i, dat->b, dat->rho_p, dat->rho_f);
	POL(6,0) = i;
	POL(7,0) = V_THETA(dX, dY, i, dat->V0, dat->rho_f, dat->a, dat->eta);

	return 0;
}

double RADIAL_FORCE (const Matrix<double>& POL, double eta, double b, double mp, double t, double a)
{
	//return ((POL(0,0)+POL(1,0)+POL(2,0))*((t/(6.0*3.1415*eta*b))-((1-exp(-6.0*3.1415*eta*b*t/mp))/(pow(6.0*3.1415*eta*b,2.0)/mp)))+POL(3,0)*t)/a;
	return (POL(3,0)*t)/a;
}

double TANGENTIAL_FORCE (const Matrix<double>& POL, const Matrix<double>& dY, double eta, double b, double mp, double t, double a, int i)
{
	//return (POL(4,0)+POL(5,0)*((t/(6.0*3.1415*eta*b))-((1-exp(-6.0*3.1415*eta*b*t/mp))/(pow(6.0*3.1415*eta*b,2.0)/mp)))+POL(7,0)*t)/(a*dY(i-1,0));
	return (POL(7,0)*t)/(a*dY(i-1,0));
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
	dat->mu_0 = 4.0*M_PI*pow(10.0, -7.0);		//permeability of free space, H/m
	dat->rho_f = 1000.0;						//Fluid density, Kg/m3
	dat->eta = 0.001;							//Dynamic viscosity, Kg/m-s
	dat->Hamaker = 1.3e-21;					//Hamaker constant for ferric oxide particle
	
	//Separator Parameters
	
	dat->Rs = 0.0026925;								//separator radius, m
	dat->L = 0.0611;									//0.05separator length, m
	dat->porosity = 0.9036;								//0.965separator porosity
	dat->V_separator=M_PI*pow(dat->Rs,2.0)*dat->L;		//volume of separator, m^3
	
	//System Parameters
	
	dat->a = 3.3e-5;										//Wire radius,m
	dat->V_wire = (1-dat->porosity)*dat->V_separator;		//Total wire volume, m^3
	dat->L_wire = dat->V_wire/M_PI/pow(dat->a,2.0);		//Total wire length, m
	dat->A_separator = M_PI*pow(dat->Rs,2.0);
	dat->A_wire = dat->V_wire/dat->L_wire;
	dat->B0 = 1.1;								//Applied magnetic induction,T
	dat->H0 = dat->B0/dat->mu_0;
	dat->Ms = 0.6;								//Saturated magnetization,T
	
	//Particle parameters
	
	dat->b = 1.0e-6;							//Particle radius,m
	dat->chi_p = 3.87e-7;						//Volume magnetic susceptibility, dimensionless
	dat->rho_p = 8700.0;						//Particle density, Kg/m3
		
	//Model parameters
	dat->Q_in = 9.5e-7;											//Volumetric fluid flow, m^3/sec
	dat->V0 = dat->Q_in/(dat->A_separator - dat->A_wire);		//Superficial velocity,m/s
	dat->Y_initial = 20.0;										//initial position in Y-axis, Y=y/a, dimensionless
	dat->dt = 2.0*dat->Y_initial*dat->a/(dat->V0*100.0);		
	
	//Breakthrough parameters
	dat->M = 2.0*dat->Ms/dat->mu_0;			
	dat->mp = 4.0*3.1415*pow(dat->b,3.0)*dat->rho_p/3.0;	//Particle Mass, Kg 
		
	return 0;
}

int Run_Trajectory()
{
	int success = 0;
	
	TRAJECTORY_DATA dat;
	
	success = Trajectory_SetupConstants(&dat);
	
	dat.POL.set_size(8,1);
	dat.H.set_size(2,1);
	dat.dX.set_size(300,1);
	dat.dY.set_size(300,1);
	dat.X.set_size(300, 1);
	dat.Y.set_size(300, 1);
	
	bool Hit = false;
	
	dat.dY(0,0) = 20.0;
	dat.dX(0,0) = 89.99/180.0*M_PI;			//FIXED negative starting values for x
	
	for (int i = 1; i<299; i++)
	{
		int Force_Polar;
		Force_Polar = POLAR(dat.POL,dat.dX,dat.dY,(void *)&dat,i);
		
		//dat.POL.Display("Polar Force");
		
		int Force_Balance;
		Force_Balance = CARTESIAN(dat.POL,dat.H,dat.dY,i,(void *)&dat);
		
		//dat.H.Display("H-Matrix<double>");
		
		int Next;
		Next = DISPLACEMENT(dat.dX,dat.dY,dat.H,i);
		
		if (dat.dY(i,0) <= (1.0+(1.85/33.0)) && dat.POL(3,0) < 0.0)
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
	
	dat.X.Display("Distance X");
	dat.Y.Display("Diatacne Y");
	
	const int nrolls=10000;  // number of experiments
	const int nstars=100;    // maximum number of stars to distribute
	
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(5.0,2.0);
	
	int p[10]={};
	
	for (int i=0; i<nrolls; ++i)
	{
		double number = distribution(generator);
		if ((number>=0.0)&&(number<10.0)) ++p[int(number)];
	}
	
	std::cout << "normal_distribution (5.0,2.0):" << std::endl;
	
	for (int i=0; i<10; ++i)
	{
		std::cout << i << "-" << (i+1) << ": ";
		std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
	}

	
	return success;
}

