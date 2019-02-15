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
	model_type = freiling;
	capfis_ratio = 0.0;
	neutrons_emit = 1.5;
	fusion_yield = 25.0;
	fission_yield = 25.0;
	total_yield = 50.0;
	casing_cap = 0.2;
	casing_den = 2600.0;
	casing_thickness = 10.0;
	casing_mw = 100.0;
	casing_thermal = 0.2;
	soil_thermal = 0.2;
	soil_scattering = 2.0;
	burst_height = 0.0;
	escape_fraction = 0.0;
	volatile_fraction = 0.0;
	soil_capture_fraction = 0.0;
}

//Destructor
ActivityDistribution::~ActivityDistribution()
{
	delete_fractionation();
	delete_capture_fractions();
	delete_casing_components();
}

void ActivityDistribution::set_model_type(asd_model type)
{
	this->model_type = type;
}

void ActivityDistribution::set_capfis_ratio(double val)
{
	this->capfis_ratio = val;
}

void ActivityDistribution::set_neutrons_emit(double val)
{
	this->neutrons_emit = val;
}

void ActivityDistribution::set_fusion_yield(double val)
{
	this->fusion_yield = val;
}

void ActivityDistribution::set_fission_yield(double val)
{
	this->fission_yield = val;
}

void ActivityDistribution::set_total_yield(double val)
{
	this->total_yield = val;
}

void ActivityDistribution::set_casing_cap(double val)
{
	this->casing_cap = val;
}

void ActivityDistribution::set_casing_den(double val)
{
	this->casing_den = val;
}

void ActivityDistribution::set_casing_thickness(double val)
{
	this->casing_thickness = val;
}

void ActivityDistribution::set_casing_mw(double val)
{
	this->casing_mw = val;
}

void ActivityDistribution::set_casing_thermal(double val)
{
	this->casing_thermal = val;
}

void ActivityDistribution::set_soil_thermal(double val)
{
	this->soil_thermal = val;
}

void ActivityDistribution::set_soil_scattering(double val)
{
	this->soil_scattering = val;
}

void ActivityDistribution::set_burst_height(double val)
{
	this->burst_height = val;
}

void ActivityDistribution::set_escape_fraction(double val)
{
	this->escape_fraction = val;
}

void ActivityDistribution::set_volatile_fraction(double val)
{
	this->volatile_fraction = val;
}

void ActivityDistribution::set_soil_capture_fraction(double val)
{
	this->soil_capture_fraction = val;
}

void ActivityDistribution::initialize_fractionation_map(std::map<double, double> & part_conc)
{
	
}

void ActivityDistribution::delete_casing_components()
{
	this->casing_atom_frac.clear();
	this->casing_mat.clear();
	this->casing_frac.clear();
}

void ActivityDistribution::add_casing_component(std::string name, double frac)
{
	
}

void ActivityDistribution::verify_casing_components()
{
	
}

void ActivityDistribution::delete_fractionation()
{
	this->nuc_fractionation.clear();
	this->freiling_rat.clear();
}

void ActivityDistribution::delete_capture_fractions()
{
	this->casing_capfrac.clear();
	this->soil_capfrac.clear();
	this->weapon_capfrac.clear();
}

asd_model ActivityDistribution::get_model_type()
{
	return this->model_type;
}

double ActivityDistribution::get_capfis_ratio()
{
	return this->capfis_ratio;
}

double ActivityDistribution::get_neutrons_emit()
{
	return this->neutrons_emit;
}

double ActivityDistribution::get_fusion_yield()
{
	return this->fusion_yield;
}

double ActivityDistribution::get_fission_yield()
{
	return this->fission_yield;
}

double ActivityDistribution::get_total_yield()
{
	return this->total_yield;
}

double ActivityDistribution::get_casing_cap()
{
	return this->casing_cap;
}

double ActivityDistribution::get_casing_den()
{
	return this->casing_den;
}

double ActivityDistribution::get_casing_thickness()
{
	return this->casing_thickness;
}

double ActivityDistribution::get_casing_mw()
{
	return this->casing_mw;
}

double ActivityDistribution::get_casing_thermal()
{
	return this->casing_thermal;
}

double ActivityDistribution::get_soil_thermal()
{
	return this->soil_thermal;
}

double ActivityDistribution::get_soil_scattering()
{
	return this->soil_scattering;
}

double ActivityDistribution::get_burst_height()
{
	return this->burst_height;
}

double ActivityDistribution::get_escape_fraction()
{
	return this->escape_fraction;
}

double ActivityDistribution::get_volatile_fraction()
{
	return this->volatile_fraction;
}

double ActivityDistribution::get_soil_capture_fraction()
{
	return this->soil_capture_fraction;
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
