/*!
 *  \file kea.h kea.cpp
 *	\brief Kernel for Estimating Activity-distribution
 *  \author Austin Ladshaw
 *	\date 02/07/2019
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay and transport. Copyright (c) 2018,
 *			   all rights reserved.
 */

#include "kea.h"
#include "crane.h"  ///< For testing purposes only.

//Model choice
asd_model activitymodel_choice(std::string &choice)
{
	asd_model type = freiling;
	
	std::string copy = choice;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	
	if (copy == "freiling")
		type = freiling;
	else if (copy == "freiling-tompkins")
		type = freiling_tompkins;
	else if (copy == "modified-freiling")
		type = mod_freiling;
	else if (copy == "modified-freiling-tompkins")
		type = mod_freiling_tompkins;
	else
		type = freiling;
	
	return type;
}

/*
 *								Start: ActivityDistribution Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Constructor
ActivityDistribution::ActivityDistribution()
{
	type = explosion;
	total_mass = 1.0;
	fiss_extent = 100.0;
	energy_level = 0.0;
	model_type = freiling;
}

//Destructor
ActivityDistribution::~ActivityDistribution()
{
	InitialMat.clear();
	MatFrac.clear();
	Yields.clear();
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: ActivityDistribution Class Definitions
 */

int KEA_TESTS()
{
	int success = 0;
	
	return success;
}
