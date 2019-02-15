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
	weapon_thermal = 0.2;
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

void ActivityDistribution::set_weapon_thermal(double val)
{
	this->weapon_thermal = val;
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
	if (frac < 0.0)
		frac = 0.0;
	if (frac > 1.0)
		frac = 1.0;
	this->casing_frac[name] = frac;
	Molecule temp;
	if (name == "Other")
		temp.Register(0, 0, 0, 0, false, false, "Solid", name, name, "Fe10C");
	else
		temp.Register(0, 0, 0, 0, false, false, "Solid", name, name, name);
	this->casing_mat[name] = temp;
}

void ActivityDistribution::verify_casing_components()
{
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	double sum = 0.0;
	int count = 0;
	for (it=this->casing_frac.begin(); it!=this->casing_frac.end(); it++)
	{
		sum += it->second;
		count++;
	}
	if (count == 0)
	{
		this->add_casing_component("Other", 1.0);
	}
	else
	{
		if ( sum > 1.0 )
		{
			for (it=this->casing_frac.begin(); it!=this->casing_frac.end(); it++)
			{
				it->second = it->second / sum;
			}
		}
		if ( sum < 1.0 )
		{
			double diff = 1.0 - sum;
			it = this->casing_frac.find("Other");
			if (it == this->casing_frac.end())
			{
				this->add_casing_component("Other", diff);
			}
			else
			{
				this->casing_frac["Other"] += diff;
			}
		}
	}
	
	//Create the soil atom map
	for (it=this->casing_frac.begin(); it!=this->casing_frac.end(); it++)
	{
		for (int i=0; i<this->casing_mat[it->first].getAtoms().size(); i++)
		{
			this->casing_atom_frac[this->casing_mat[it->first].getAtoms()[i].AtomSymbol()] = 0.0;
		}
	}
	//Fill the soil atom map
	for (it=this->casing_frac.begin(); it!=this->casing_frac.end(); it++)
	{
		for (int i=0; i<this->casing_mat[it->first].getAtoms().size(); i++)
		{
			this->casing_atom_frac[this->casing_mat[it->first].getAtoms()[i].AtomSymbol()] += it->second*(1.0/(double)this->casing_mat[it->first].getAtoms().size());
		}
	}
	
	for (it=this->casing_atom_frac.begin(); it!=this->casing_atom_frac.end(); it++)
	{
		Atom temp;
		temp.Register(it->first);
		this->casing_atom[it->first] = temp;
	}
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

double ActivityDistribution::get_weapon_thermal()
{
	return this->weapon_thermal;
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

void ActivityDistribution::compute_neutrons_emit(double fission, double fusion)
{
	this->set_neutrons_emit( 1.5 + (3.0*fusion/fission) );
}

void ActivityDistribution::compute_casing_mw()
{
	double mw = 0.0;
	
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=this->casing_frac.begin(); it!=this->casing_frac.end(); it++)
	{
		mw += it->second*this->casing_mat[it->first].MolarWeight();
	}
	this->set_casing_mw(mw);
}

void ActivityDistribution::compute_casing_thermal()
{
	double sum = 0.0;
	
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=this->casing_atom_frac.begin(); it!=this->casing_atom_frac.end(); it++)
	{
		sum += it->second*this->casing_atom[it->first].ThermalXSection();
	}
	this->set_casing_thermal(sum);
}

void ActivityDistribution::compute_soil_thermal(std::map<std::string, double> & soil_atom_frac, std::map<std::string, Atom> & soil_atom)
{
	double sum = 0.0;
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=soil_atom_frac.begin(); it!=soil_atom_frac.end(); it++)
	{
		sum += it->second*soil_atom[it->first].ThermalXSection();
	}
	this->set_soil_thermal(sum);
}

void ActivityDistribution::compute_soil_scattering(std::map<std::string, double> & soil_atom_frac, std::map<std::string, Atom> & soil_atom)
{
	double sum = 0.0;
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=soil_atom_frac.begin(); it!=soil_atom_frac.end(); it++)
	{
		sum += it->second*soil_atom[it->first].ScatterXSection();
	}
	this->set_soil_scattering(sum);
}

void ActivityDistribution::compute_weapon_thermal(FissionProducts & weapon)
{
	double sum = 0.0;
	for (int i=0; i<weapon.getWeaponMat().size(); i++)
	{
		sum += weapon.getWeaponFrac()[i]/100.0*weapon.getWeaponMat()[i].ThermalXSection();
	}
	this->set_weapon_thermal(sum);
}

void ActivityDistribution::compute_casing_capfrac()
{
	this->compute_casing_thermal();
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=this->casing_atom_frac.begin(); it!=this->casing_atom_frac.end(); it++)
	{
		this->casing_capfrac[it->first] = it->second*this->casing_atom[it->first].ThermalXSection()/this->get_casing_thermal();
		std::cout << it->first << "\t" << this->casing_capfrac[it->first] << std::endl;
	}
}

void ActivityDistribution::compute_soil_capfrac(std::map<std::string, double> & soil_atom_frac, std::map<std::string, Atom> & soil_atom)
{
	this->compute_soil_thermal(soil_atom_frac, soil_atom);
	this->compute_soil_scattering(soil_atom_frac, soil_atom);
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=soil_atom_frac.begin(); it!=soil_atom_frac.end(); it++)
	{
		this->soil_capfrac[it->first] = it->second*soil_atom[it->first].ThermalXSection()/this->get_soil_thermal();
		std::cout << it->first << "\t" << this->soil_capfrac[it->first] << std::endl;
	}
}

void ActivityDistribution::compute_weapon_capfrac(FissionProducts & weapon)
{
	this->compute_weapon_thermal(weapon);
	for (int i=0; i<weapon.getWeaponMat().size(); i++)
	{
		this->weapon_capfrac[weapon.getWeaponMat()[i].IsotopeName()] = 0.0;
	}
	for (int i=0; i<weapon.getWeaponMat().size(); i++)
	{
		this->weapon_capfrac[weapon.getWeaponMat()[i].IsotopeName()] += weapon.getWeaponFrac()[i]/100.0*weapon.getWeaponMat()[i].ThermalXSection()/this->get_weapon_thermal();
	}
	
	//Iterate through the map
	std::map<std::string,double>::iterator it;
	for (it=this->weapon_capfrac.begin(); it!=this->weapon_capfrac.end(); it++)
	{
		std::cout << it->first << "\t" << this->weapon_capfrac[it->first] << std::endl;
	}
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: ActivityDistribution Class Definitions
 */

int KEA_TESTS()
{
	int success = 0;
	double time;
	
	std::cout << "\nTesting of the KEA module for the individual functions\n";
	std::cout <<   "------------------------------------------------------\n\n";
	
	ActivityDistribution test;
	Crane cranetest;
	FissionProducts yeildtest;
	Dove dove;
	yaml_cpp_class nuc_data;
	time = clock();
	
	/// ----------------- Initialize FAIRY for the Weapon Components ----------------------------------
	nuc_data.executeYamlRead("database/NuclideLibrary.yml");
	yeildtest.loadNuclides(nuc_data);
	
	yeildtest.setTotalMass(9.7);
	yeildtest.setFissionExtent(99.0);
	yeildtest.setFissionType(explosion);
	yeildtest.setEnergyLevel(1000.0);
	yeildtest.timeSinceDetonation(3.0, 99.0); //Set yeild to 3 sec cutoff based on 99% conversion to daughters
	yeildtest.addIsotopeMaterial("U-235", 90.0);
	yeildtest.addIsotopeMaterial("U-238", 10.0);
	
	yeildtest.loadFissionProductYields();
	yeildtest.evaluateYields();
	//Don't form eigenvectors yet... sets ICs for FPY data
	
	/// ----------------- END Initialize FAIRY for the Weapon Components ----------------------------------
	
	FILE *file, *cloud;
	file = fopen("output/KEA_Tests.txt", "w+");
	cloud = fopen("output/KEA_Tests_CloudGrowth.txt", "w+");
	if (file == nullptr)
	{
		system("mkdir output");
		file = fopen("output/KEA_Tests.txt", "w+");
	}
	if (cloud == nullptr)
	{
		system("mkdir output");
		cloud = fopen("output/KEA_Tests_CloudGrowth.txt", "w+");
	}
	cranetest.set_CloudFile(cloud);
	
	//Weapon casing is typically steel with a lead or lead-bismuth alloy lining
	test.add_casing_component("Fe10C", 0.9);
	test.add_casing_component("PbBi", 0.1);
	test.verify_casing_components();
	
	test.set_fission_yield(25.0);
	test.set_fusion_yield(25.0);
	test.set_total_yield(test.get_fission_yield()+test.get_fusion_yield());
	test.compute_neutrons_emit(test.get_fission_yield(), test.get_fusion_yield());
	
	test.compute_casing_mw();
	//test.compute_casing_thermal();
	
	/// ----------------- Initialize CRANE for the Soil Components ----------------------------------
	
	//V. Jodoin Test Case from 1994 Thesis
	double W = 50.0; //50 kT
	double hb = 0.0;// 0 m
	double gz = 500.0; //500 m
	
	int bins = 10;
	bool includeShear = true;
	bool isTight = true;
	
	cranetest.add_soil_component("SiO2", 0.75);
	cranetest.add_soil_component("CaO", 0.25);
	cranetest.establish_initial_conditions(dove, W, gz, hb, bins, includeShear, isTight);
	
	std::cout << "Shear Velocity          = \t";
	if (includeShear == true)
	std::cout << "True\n";
	else
	std::cout << "False\n";
	std::cout << "Bomb Yield (kT)         =\t" << W << std::endl;
	std::cout << "Burst Height (m)        =\t" << hb << std::endl;
	std::cout << "Ground Altitude (m)     =\t" << gz << std::endl;
	std::cout << "Initial Time (s)        =\t" << cranetest.get_current_time() << std::endl;
	std::cout << "Number of air parcels   =\t" << cranetest.return_parcel_alt_top().rows() << std::endl;
	std::cout << "Number of particle bins =\t" << cranetest.return_parcel_alt_top().columns() << std::endl;
	cranetest.display_part_hist();
	std::cout << "Soil Solid. Temp. (K)   =\t" << cranetest.get_solidification_temp() << std::endl;
	std::cout << "Soil Vapor. Temp. (K)   =\t" << cranetest.get_vaporization_temp() << std::endl;
	std::cout << "Vaporized Soil (kg)     =\t" << cranetest.get_initial_soil_vapor() << std::endl;
	cranetest.display_soil_characteristics();
	std::cout << "\n";
	
	/// ----------------- END Initialize CRANE for the Soil Components ----------------------------------
	
	//test.compute_soil_thermal(cranetest.get_soil_atom_frac(), cranetest.get_soil_atom());
	//test.compute_soil_scattering(cranetest.get_soil_atom_frac(), cranetest.get_soil_atom());
	//test.compute_weapon_thermal(yeildtest);
	
	test.compute_casing_capfrac();
	test.compute_soil_capfrac(cranetest.get_soil_atom_frac(), cranetest.get_soil_atom());
	test.compute_weapon_capfrac(yeildtest);
	
	if (file!= nullptr)
		fclose(file);
	if (cloud!=nullptr)
		fclose(cloud);
	
	time = clock() - time;
	std::cout << "\nTest Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
