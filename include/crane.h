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
	double eps;															///< Ratio of molecular wieghts for water-vapor and dry air
	double grav;														///< Acceleration from gravity (m/s^2)
	double k2;															///< Dimensionless power function yield
	std::map<double, double> amb_temp;									///< Ambient Temperature (K) at various altitudes (m)
	std::map<double, double> atm_press;									///< Atmospheric Pressure (Pa) at various altitudes (m)
	std::map<double, double> rel_humid;									///< Relative Humidity (%) at various altitudes (m)
	/// Wind Velocities (m/s) at various altitudes (m)
	/** Velocities stored in x (0,0) and y (1,0) components at a given altitude */
	std::map<double, Matrix<double> > wind_vel;
	std::map<double, double> num_part;								///< Number of particles per cloud volume (#/m^3) at given particle size (m)
	
private:
	
};

/// Test function for CRANE
/**  Test function is callable from the cli */
int CRANE_TESTS();
#endif /* CRANE_HPP_ */
