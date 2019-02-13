/*!
 *  \file kea.h kea.cpp
 *	\brief Kernel for Estimating Activity-distribution
 *	\details This file contains a C++ object for determining the distribution of
 *			radioactivity particles and activity onto debris particles in specific
 *			size classes. It is directly coupled with FAIRY to determine yields of
 *			nuclides from a specific nuclear event, then will be integrated into
 *			CRANE to establish the distribution of nuclides in the debris cloud.
 *			For the sake of modularity, this kernel will not be coupled with CRANE
 *			and instead CRANE will integrate this kernel into its source. Thus, the
 *			activity-distribution in KEA will be determined from information that is
 *			anticipated to be passed to the functions and objects developed here.
 *			That will allow for independent testing of this kernel and allow for 
 *			changes to how the activity-distribution is determined to be made on
 *			the fly. 
 *
 *			References for Activity-Distribution
 *			------------------------------------
 *			E.C. Freiling, "Radionuclide Fractionation in Bomb Debris," Science,
 *				1991-1998, 1961.
 *
 *			J.T. McGahan, E.J. Kownaki, "Sensitivity of Fallout Predictions to 
 *				Initial Conditions and Model Assumptions," Defense Nuclear Agency,
 *				DNA-3439F, 1974.
 *
 *			H.G. Norment, "DELFIC: Department of Defense Fallout Prediction System:
 *				Volume I - Fundamentals," Defense Nuclear Agency, DNA-5159F, 1979.
 *
 *			R.C. Tompkins, "DELFIC: Department of Defense Fallout Prediction System:
 *				Volume V - Pacticle Activity," US Army Nuclear Defense Laboratory, 
 *				DASA-1800-V, 1968.
 *
 *			D.A. Hooper, V.J. Jodoin, "Revision of the DELFIC Particle Activity Module,"
 *				Oak Ridge National Laboratory, ORNL/TM-2010/220, 2010.
 *
 *			Reference for Induced-Soil Activity
 *			-----------------------------------
 *			T.H. Jones, "A Prediction System for the Neutron-Induced Activity Contribution
 *				to Fallout Exposure Rates," U.S. Naval Radiological Defense Laboratory,
 *				USNRDL-TR-1056, 1966.
 *
 *			Reference for Neutron Absorption and Scattering Cross Sections (in EEL)
 *			-----------------------------------------------------------------------
 *			V.F. Sears, "Neutron Scattering Lengths and Cross Sections," Neutron News, 3,
 *				26-37, 1992. 
 *
 *  \author Austin Ladshaw
 *	\date 02/07/2019
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay and transport. Copyright (c) 2018,
 *			   all rights reserved.
 */

#ifndef KEA_HPP_
#define KEA_HPP_

#include "fairy.h"
#include "mola.h"

/// Enumeration for the list of valid activity-size distribution methods
/** List of valid models for activity-size distributions.*/
typedef enum {freiling, freiling_tompkins, mod_freiling, mod_freiling_tompkins} asd_model;

/// Function to determine the activity-size distribution model type
asd_model activitymodel_choice(std::string &choice);

/// C++ Object for determining the activity-size distribution
/** This object inherits from FissionProducts and will determine the activity-size distributions
	of nuclides in a nuclear debris cloud. It will also be used to determine the induced activity
	in soil particles and from weapon material absorbtion of neutrons. While this kernel is developed
	independently from CRANE, it will be coupled with the size distributions and other parameters
	from CRANE. Then CRANE and KEA will be implemented together to fully describe the nuclear debris
	cloud post-detonation and to the time of cloud stabilization.
 */
class ActivityDistribution
{
public:
	ActivityDistribution();											///< Default constructor
	~ActivityDistribution();										///< Default destructor
	
protected:
	asd_model model_type;											///< Type of activity-size distribution model to use
    // capfis_ratio = No*(fc)_i
	double capfis_ratio;											///< Neutron capture-to-fission ratio for induced activity
	
	/// Below are all the parameters associated with the induced-soil-activity models
	double neutrons_emit;											///< Neutrons emitted per fission (No)
	double fusion_yield;											///< Fusion yield in kT (Wfu)
	double fission_yield;											///< Fission yield in kT (Wfis)
    double total_yield;                                             ///< Total weapon yield in kT (W)
	double casing_cap;												///< Weapon casing capture (Sigma)
	double casing_den;												///< Weapon casing material density in g/cm^3 (rho_c)
	double casing_thickness;										///< Weapon casing material thickness in cm (X)
	double casing_mw;												///< Weapon casing average molecular weight in g/mol (A)
	double casing_thermal;											///< Weapon casing average thermal neutron x-sec in barns (sigma_c)
	std::map<std::string, Molecule> casing_mat;						///< Weapon casing molecular composition
	std::map<std::string, double> casing_frac;						///< Weapon casing molefractions
    std::map<std::string, Isotope> weapon_mat;                      ///< Weapon molecular composition
    std::map<std::string, double> weapon_frac;                      ///< Weapon molefractions
	double burst_height;											///< Weapon burst height above ground (ft)
    
    /// Below are the parameters associated with the activity-size distributions
    std::map<double, FissionProducts> nuc_fractionation;            ///< Fractionation of nuclides with particle size (um)
    std::map<int, double> freiling_rat;                             ///< Freiling ratios for each mass number chain
	
private:
	
};

/// Test function for KEA
int KEA_TESTS();

#endif /* KEA_HPP_ */
