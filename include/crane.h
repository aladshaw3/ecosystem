/*!
 *  \file crane.h
 *	\brief Cloud Rise After Nuclear Explosion
 *	\details This file creates objects and subroutines for solving the systems of equations for the mass,
 *			energy, and temperature of debris clouds caused by nuclear detonations. The equations solved
 *			and methods used come from the DEfense Land Fallout Interpretative Code (DELFIC) developed by
 *			U.S. DOD in the 1960s to 1970s. The original DELFIC software was written in Fortran77 and the
 *			source code is not available to the public. This software is a recreation of the Cloud Rise
 *			Module from DELFIC based on the reports made publically available. In this software, we are
 *			only interested in estimating the cloud rise and the shape of the nuclear debris cloud post-
 *			detonation of a nuclear weapon. This software does not perform any transport of the resulting
 *			fallout cloud of debris. Transport will be handled by a different code for modeling systems
 *			of PDEs. The cloud rise simulation performed here will become the initial conditions for a
 *			transport model that is to be developed later.
 *
 *			References for DELFIC
 *			------------------------------
 *			H.G. Normet, "DELFIC: Department of Defense Fallout Prediction System - Volume I - Fundamentals," 
 *				U.S. DOD, DNA-001-76-C-0010, DNA 5159F-1, December 1979.
 *
 *			H.G. Normet, "DELFIC: Department of Defense Fallout Prediction System - Volume II - User's Manual,"
 *				U.S. DOD, DNA-001-76-C-0010, DNA 5159F-2, December 1979.
 *
 *			H.G. Normet and S. Woolf, "Department of Defense Land Fallout Prediction System - Volume III - Cloud
 *				Rise," U.S. DOD, DASA01-69-C-0077, DASA-1800-III, September 1970.
 *
 *
 *  \author Austin Ladshaw
 *	\date 05/30/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of radioactive particle aggregation and transport. Copyright (c) 2018,
 *             all rights reserved.
 */

#include "dove.h"

#ifndef CRANE_HPP_
#define CRANE_HPP_

/// CRANE object to hold data and functions associated with Cloud Rise
/** This is a C++ object that contains data and functions associated with calculating debris cloud
	rise from nuclear detonations. This object must be passed to the Dove object and registered as
	the user defined data structure. Then, inside of residual functions developed for each non-linear
	variable, this object must be appropriately dereferenced so that its members and functions can
	be called appropriately. 
 
	\note Crane interfaces with Dove, but does not contain an instance of Dove. Also, Crane does not
	contain data for the non-linear variables since Dove holds the vectors of non-linear variables.*/
class Crane
{
public:
	Crane();															///< Default Constructor
	~Crane();															///< Default Destructor
	
protected:
	double eps;															///< Ratio of molecular wieghts for water-vapor and dry air (mol/g / mol/g)
	double grav;														///< Acceleration from gravity (m/s^2)
	double k;															///< Dimensionless initial cloud rise yield
	double k2;															///< Dimensionless power function yield
	double k3;															///< Dimensionless kinetic energy yield
	double k6;															///< Dimensionless wind shear factor
	double mu;															///< Dimensionless energy yield for given explosion
	double k_temp;														///< Dimensionless temperature factor for specific heat
	double apparent_temp;												///< Apparent temperature of cloud (K)
	double apparent_amb_temp;											///< Apparent ambient temperature at specific altitude (K)
	double q_x;															///< Ratio of apparent to actual temperature in the cloud
	double q_xe;														///< Ratio of apparent to actual ambient temperature
	double beta_prime;													///< Ratio of cloud gas density to total cloud density
	double xe;															///< Ratio of water-vapor to dry air for ambient air at given altitude
	double char_vel;													///< Characteristic velocity for cloud rise (m/s)
	double vert_rad;													///< Vertical radius of the cloud shape (m)
	double horz_rad;													///< Horizontal radius of the cloud shape (m)
	double virtual_mass;												///< Initial virtual mass of the cloud (kg)
	double gas_const;													///< Gas law constant for dry air (J/kg/K)
	double latent_heat;													///< Latent heat of evaporation of water (or ice) (J/kg)
	double sigma_turbulence;											///< Rate of dissipation of turbulent energy per mass (K/s)
	double mean_spec_heat;												///< Weighted mean specific heat of the cloud (J/kg/K)
	double actual_spec_heat;											///< Actual specific heat of the cloud mass (J/kg/K)
	double spec_heat_water;												///< Specific heat of water vapor (J/kg/K)
	double spec_heat_entrain;											///< Specific heat of entrained air (J/kg/K)
	double spec_heat_conds;												///< Specific heat of condensed matter (J/kg/K)
	double cloud_volume;												///< Volume of the debris cloud (m^3)
	double equil_temp;													///< Initial thermal equilibrium temperature (K)
	double total_mass_fallout_rate;										///< Overall rate of particle fallout (kg/s)
	double surf_area;													///< Surface area of the cloud (m^2)
	double shear_ratio;													///< (Surf_area/Vol)*mu*char_vel   Can be swapped to include shear (1/s)
	double shear_vel;													///< Magnitude of shear velocity impacting cloud shape (m/s)
	double part_density;												///< Average density of soil particles in cloud (kg/m^3)
	double adjusted_height;												///< Initial adjusted height of the cloud (m)
	double bomb_yield;													///< Size of nuclear bomb or explosion (kT)
	double force_factor;												///< Dimensionless factor of explosive force
	double det_height;													///< Height of detontation above mean sea level (m)
	double burst_height;												///< Height of detonation above ground level (m)
	double ground_height;												///< Height/altitude of the ground level above sea level (m)
	double energy_frac;													///< Fraction of energy available to heat the air
	double eccentricity;												///< Eccentricity of an oblate spheriod (default = 0.75)
	
	/// Need to double check units on variables below
	double air_density;													///< Density of air in cloud (kg/m^3)
	double air_viscosity;												///< Viscosity of air in cloud ...
	/// Note: Each particle size has its own slip factor and ND (Davies Number?)
	double slip_factor;													///< Slip factor for particle settling ...
	double ND;															///< Unitless number for particle settling analysis ...
	
	std::map<double, double> amb_temp;									///< Ambient Temperature (K) at various altitudes (m)
	std::map<double, double> atm_press;									///< Atmospheric Pressure (Pa) at various altitudes (m)
	std::map<double, double> rel_humid;									///< Relative Humidity (%) at various altitudes (m)
	/// Wind Velocities (m/s) at various altitudes (m)
	/** Velocities stored in x (0,0) and y (1,0) components at a given altitude */
	std::map<double, Matrix<double> > wind_vel;
	std::map<double, double> part_hist;									///< Normalized Histogram of Particle Distribution by size (um)
	std::map<double, double> settling_rate;								///< Particle settling rate (m/s) by particle size (um)
	double min_dia;														///< Minimum particle diameter for distribution (um)
	double max_dia;														///< Maximum particle diameter for distribution (um)
	double mean_dia;													///< Mean particle diameter for distribution (um)
	double std_dia;														///< Standard deviation for lognormal distribution of particle sizes
	int num_bins;														///< Number of desired size bins for particles
	
	/// List of Variables to be solved for by Dove
	/** These variables are stored internally by Dove,
		but will be placed into these parameters below
		after completing a solve for convenience. */
	double cloud_mass;													///< Total mass of the debris cloud (kg)
	double entrained_mass;												///< Mass of entrained materials in debris cloud (kg)
	double cloud_rise;													///< Rate of cloud rise through atmosphere (m/s)
	double cloud_height;												///< Height of the center of the coud above mean sea level (m)
	double x_water_vapor;												///< Mixing ratio for water vapor to dry air (kg/kg)
	double w_water_conds;												///< Mixing ratio for condensed water to dry air (kg/kg)
	double s_soil;														///< Mixing ratio for suspended soils to dry air (kg/kg)
	double temperature;													///< Temperature of the air in the cloud (K)
	double energy;														///< Mass-less Kinetic Energy in the cloud (m^2/s^2)
	std::map<double, double> part_conc;									///< Number of particles per volume of cloud (#/m^3) for given size (um)
	double current_time;												///< Current time since explosion (s)
	
private:
	
};

/// Test function for CRANE
/**  Test function is callable from the cli */
int CRANE_TESTS();
#endif /* CRANE_HPP_ */
