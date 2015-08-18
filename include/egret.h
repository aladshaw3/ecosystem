//----------------------------------------
//  Created by Austin Ladshaw on 1/29/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "macaw.h"

#ifndef EGRET_HPP_
#define EGRET_HPP_

#ifndef Rstd
#define Rstd 8.3144621						// J/K/mol (or) L*kPa/K/mol (Standard Units)
#endif

#ifndef RE3
#define RE3 8.3144621E+3					// cm^3*kPa/K/mol	(Convenient for density calculations)
#endif

#ifndef Po
#define Po 100.0							// Standard state pressure (kPa)
#endif

#ifndef Cstd
#define Cstd(p,T) ((p)/(Rstd*T))				// Calculation of concentration from partial pressure (Cstd = mol/L)
#endif										// Note: can also be used to calculate densities if p = PT * MW

#ifndef CE3
#define CE3(p,T) ((p)/(RE3*T))				// Calculation of concentration from partial pressure (CE3 = mol/cm^3)
#endif										// Note: can also be used to calculate densities if p = PT * MW

#ifndef Pstd
#define Pstd(c,T) ((c)*Rstd*T)				// Calculation of partial pressure from concentration (c = mol/L)
#endif

#ifndef PE3
#define PE3(c,T) ((c)*RE3*T)					// Calculation of partial pressure from concentration (c = mol/cm^3)
#endif

#ifndef Nu
#define Nu(mu,rho) ((mu)/(rho))					// Calculation of kinematic viscosity from dynamic vis. and density (cm^2/s)
#endif

#ifndef PSI
#define PSI(T) (0.873143 + (0.000072375*T)) // Calculation of temperature correction factor for dynamic vis.
#endif

#ifndef Dp_ij
#define Dp_ij(Dij,PT) ((PT*Dij)/Po)			// Calculation of the corrected binary diffusivity (cm^2/s)
#endif

#ifndef D_ij								// Calculation of binary diffusion based on MW, density, viscosity info (cm^2/s)
#define D_ij(MWi,MWj,rhoi,rhoj,mui,muj) ( (4.0 / sqrt(2.0)) * pow(((1/MWi)+(1/MWj)),0.5) ) / pow( (pow((pow((rhoi/(1.385*mui)),2.0)/MWi),0.25)+ pow((pow((rhoj/(1.385*muj)),2.0)/MWj),0.25)),2.0 )
#endif

#ifndef Mu									// Calculation of single species viscosity from Sutherland's Equ (g/cm/s)
#define Mu(muo,To,C,T) (muo * ((To + C)/(T + C)) * pow((T/To),1.5) )
#endif

#ifndef D_ii
#define D_ii(rhoi,mui) (1.385*mui/rhoi)		// Calculation of self-diffusivity (cm^2/s)
#endif

#ifndef ReNum
#define ReNum(u,L,nu) (u*L/nu)				// Calculation of the Reynold's Number (-)
#endif

#ifndef ScNum
#define ScNum(nu,D) (nu/D)					// Calculation of the Schmidt Number (-)
#endif

#ifndef FilmMTCoeff							// Calculation of film mass transfer coefficient (cm/s)
#define FilmMTCoeff(D,L,Re,Sc) ((D/L)*(2.0 + (1.1*pow(Re,0.6)*pow(Sc,0.3))))
#endif

typedef struct
{	
	//Constants
	double molecular_weight;				//Given: molecular weights (g/mol)
	double Sutherland_Temp;					//Given: Sutherland's Reference Temperature (K)
	double Sutherland_Const;				//Given: Sutherland's Constant (K)
	double Sutherland_Viscosity;			//Given: Sutherland's Reference Viscosity (g/cm/s)
	double specific_heat;					//Given: Specific heat of the gas (J/g/K)
	
	//Parameters
	double molecular_diffusion;				//Calculated: molecular diffusivities (cm^2/s)
	double dynamic_viscosity;				//Calculated: dynamic viscosities (g/cm/s)
	double density;							//Calculated: gas densities (g/cm^3) {use RE3}
	double Schmidt;							//Calculated: Value of the Schmidt number (-)
		
}PURE_GAS;

typedef struct
{
	//Constants
	int N;									//Given: Total number of gas species
	bool CheckMolefractions = true;			//Given: True = Check Molefractions for errors
	
	//Variables
	double total_pressure;					//Given: Total gas pressure (kPa)
	double gas_temperature;					//Given: Gas temperature (K)
	double velocity;						//Given: Gas phase velocity (cm/s)
	double char_length;						//Given: Characteristic Length (cm)
	std::vector<double> molefraction;		//Given: Gas molefractions of each species (-)
	
	//Parameters
	double total_density;					//Calculated: Total gas density (g/cm^3) {use RE3}
	double total_dyn_vis;					//Calculated: Total dynamic viscosity (g/cm/s)
	double kinematic_viscosity;				//Calculated: Kinematic viscosity (cm^2/s)
	double total_molecular_weight;			//Calculated: Total molecular weight (g/mol)
	double total_specific_heat;				//Calculated: Total specific heat (J/g/K)
	double Reynolds;						//Calculated: Value of the Reynold's number	(-)
	Matrix binary_diffusion;				//Calculated: Tensor matrix of binary gas diffusivities (cm^2/s)
	
	//All Species Info
	std::vector<PURE_GAS> species_dat;		//Vector of the pure gas info of all species
	
}MIXED_GAS;

int initialize_data(int N, MIXED_GAS *gas_dat);

int set_variables(double PT, double T, double us, double L, std::vector<double> &y, MIXED_GAS *gas_dat);

int calculate_properties(MIXED_GAS *gas_dat);

int EGRET_TESTS();

#endif
