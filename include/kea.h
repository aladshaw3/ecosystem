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

/// Enumeration for the list of valid activity-size distribution methods
/** List of valid models for activity-size distributions.*/
typedef enum {freiling, freiling_tompkins, mod_freiling, mod_freiling_tompkins} asd_model;

/// Test function for KEA
int KEA_TESTS();

#endif /* KEA_HPP_ */
