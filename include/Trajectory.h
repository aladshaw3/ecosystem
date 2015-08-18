#include "macaw.h"
#include <random>

typedef struct
{
	//Constants
	double mu_0 = 12.57e-7;						//permeability of free space, H/m
	double rho_f = 1000.0;						//Fluid density, Kg/m3
	double eta = 0.001;							//Dynamic viscosity, Kg/m-s
	double Hamaker = 1.3e-21;					//Hamaker constant for ferric oxide particles

	//Separator Parameters

	double Rs = 0.0026925;						//separator radius, m
	double L = 0.0611;							//0.05separator length, m
	double porosity = 0.9036;					//0.965separator porosity
	double V_separator;							//volume of separator, m^3

	//System Parameters

	double a = 33.0e-6;								//Wire radius,m
	double V_wire;									//Total wire volume, m^3
	double L_wire;									//Total wire length, m
	double A_separator;
	double A_wire;
	double B0 = 1.0;								//Applied magnetic induction,T
	double H0;
	double Ms = 0.6;								//Saturated magnetization,T

	//Particle parameters

	double b = 1.0e-6;								//Particle radius,m
	double chi_p = 3.87e-6;							//Volume magnetic susceptibility, dimensionless
	double rho_p = 8700.0;							//Particle density, Kg/m3

	//Model parameters
	double Q_in = 9.5e-7;							//Volumetric fluid flow, m^3/sec
	double V0;										//Superficial velocity,m/s
	double Y_initial = 20.0;						//initial position in Y-axis, Y=y/a, dimensionless
	double dt;										//Time step size

	//Breakthrough parameters
	double M;
	double mp;										//Particle Mass, Kg
	
	Matrix POL,H,dX,dY;
	Matrix X, Y;

}TRAJECTORY_DATA;

double Magnetic_R(const Matrix& dX, const Matrix& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a);

double Magnetic_T(const Matrix& dX, const Matrix& dY, int i, double b, double mu_0, double chi_p, double M, double H0, double a);

double Grav_R(const Matrix& dX, int i, double b, double rho_p, double rho_f);

double Grav_T(const Matrix& dX, int i, double b, double rho_p, double rho_f);

double Van_R(const Matrix& dX, const Matrix& dY, int i, double Hamaker, double b, double a);

double V_RAD (const Matrix& dX, const Matrix& dY, int i, double V0, double rho_f, double a, double eta);

double V_THETA (const Matrix& dX, const Matrix& dY, int i, double V0, double rho_f, double a, double eta);

int POLAR(Matrix& POL, const Matrix& dX, const Matrix& dY, const void *data, int i);

double RADIAL_FORCE (const Matrix& POL, double eta, double b, double mp, double t, double a);

double TANGENTIAL_FORCE (const Matrix& POL, const Matrix& dY, double eta, double b, double mp, double t, double a, int i);

int CARTESIAN(const Matrix& POL, Matrix& H, const Matrix& dY, double i, const void *data);

int DISPLACEMENT (Matrix& dX, Matrix& dY, const Matrix& H, int i);

int LOCATION (const Matrix& dY, const Matrix& dX, Matrix& XX, Matrix& YY, int i);

int Trajectory_SetupConstants(TRAJECTORY_DATA *dat);

int Run_Trajectory();