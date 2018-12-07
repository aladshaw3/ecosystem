/*!
 *  \file fairy.h fairy.cpp
 *	\brief Fission-products from Atomic Incident and their Respective Yields
 *  \author Austin Ladshaw
 *	\date 12/07/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay and transport. Copyright (c) 2018,
 *				all rights reserved.
 */

#include "fairy.h"

//Pick fission type
fiss_type fisstype_choice(std::string &choice)
{
	fiss_type type = neutron;
	
	std::string copy = choice;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	
	if (copy == "neutron" || copy == "neutron-induced")
		type = neutron;
	else if (copy == "spontaneous")
		type = spontaneous;
	else if (copy == "explosion" || copy == "bomb" || copy == "blast")
		type = explosion;
	else
		type = neutron;
	
	return type;
}


/*
 *								Start: FissionProducts Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
FissionProducts::FissionProducts()
{
	type = neutron;
	total_mass = 1.0;
	fiss_extent = 100.0;
	energy_level = 0.0;
}

//Default destructor
FissionProducts::~FissionProducts()
{
	InitialMat.clear();
	MatFrac.clear();
}

//Display information
void FissionProducts::DisplayInfo()
{
	std::cout << "Type of Fission Event:\t";
	switch (this->type)
	{
		case neutron:
			std::cout << "Neutron-Induced\n";
			break;
			
		case spontaneous:
			std::cout << "Spontaneous\n";
			break;
			
		case explosion:
			std::cout << "Explosion\n";
			break;
			
		default:
			std::cout << "Neutron-Induced\n";
			break;
	}
	std::cout << "Neutron Source Energy (eV) =\t" << this->energy_level << std::endl;
	std::cout << "Extent of Fission (%)  =\t" << this->fiss_extent << std::endl;
	std::cout << "Mass of Materials (kg) = \t" << this->total_mass << std::endl;
	std::cout << "\tMaterial\tPercentage\n";
	std::cout << "\t--------\t----------\n";
	for (int i=0; i<this->InitialMat.size(); i++)
	{
		std::cout << "\t" << this->InitialMat[i].IsotopeName() << "\t\t" << this->MatFrac[i] << std::endl;
	}
}

//Load library
void FissionProducts::loadNuclides(yaml_cpp_class &data)
{
	this->DecayChain::loadNuclides(data);
}

//Unload library
void FissionProducts::unloadNuclides()
{
	this->DecayChain::unloadNuclides();
}

//Load FPY library
void FissionProducts::loadFissionProductYields()
{
	switch (this->type)
	{
		case neutron:
			this->fpy_data.executeYamlRead("database/NeutronFissionProductYields.yml");
			break;
			
		case spontaneous:
			this->fpy_data.executeYamlRead("database/SpontaneousFissionProductYields.yml");
			break;
			
		case explosion:
			this->fpy_data.executeYamlRead("database/NeutronFissionProductYields.yml");
			break;
			
		default:
			this->fpy_data.executeYamlRead("database/NeutronFissionProductYields.yml");
			break;
	}
}

//Unload FPY library
void FissionProducts::unloadFissionProductYields()
{
	this->fpy_data.DeleteContents();
}

//Set type
void FissionProducts::setFissionType(fiss_type opt)
{
	this->type = opt;
}

//Set mass
void FissionProducts::setTotalMass(double mass)
{
	if (mass <= 0.0)
		mass = 1.0;
	this->total_mass = mass;
}

//Set fission extent
void FissionProducts::setFissionExtent(double per)
{
	if (per > 100.)
		per = 100;
	if (per <= 0.0)
		per = 0.0;
	this->fiss_extent = per;
}

//Set energy
void FissionProducts::setEnergyLevel(double el)
{
	if (el < 0.0)
		el = 0.0;
	this->energy_level = el;
}

//Add isotope
void FissionProducts::addIsotopeMaterial(std::string iso, double percent)
{
	Isotope temp;
	temp.loadNuclides(*this->nuclides);
	temp.registerIsotope(iso);
	this->InitialMat.push_back(temp);
	if (percent > 100.)
		percent = 100;
	if (percent <= 0.0)
		percent = 0.0;
	this->MatFrac.push_back(percent);
}

//Add isotope
void FissionProducts::addIsotopeMaterial(std::string sym, int iso_num, double percent)
{
	Isotope temp;
	temp.loadNuclides(*this->nuclides);
	temp.registerIsotope(sym, iso_num);
	this->InitialMat.push_back(temp);
	if (percent > 100.)
		percent = 100;
	if (percent <= 0.0)
		percent = 0.0;
	this->MatFrac.push_back(percent);
}

//Add isotope
void FissionProducts::addIsotopeMaterial(int atom_num, int iso_num, double percent)
{
	Isotope temp;
	temp.loadNuclides(*this->nuclides);
	temp.registerIsotope(atom_num, iso_num);
	this->InitialMat.push_back(temp);
	if (percent > 100.)
		percent = 100;
	if (percent <= 0.0)
		percent = 0.0;
	this->MatFrac.push_back(percent);
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: FissionProducts Class Definitions
 */

/// Test function for FAIRY
int FAIRY_TESTS()
{
	//Initializations
	int success = 0;
	double time;
	yaml_cpp_class nuc_data;
	std::cout << "\nRunning FAIRY tests...\n";
	
	//Declarations
	time = clock();
	nuc_data.executeYamlRead("database/NuclideLibrary.yml");
	FissionProducts test;
	test.loadNuclides(nuc_data);
	test.loadFissionProductYields();
	
	test.setTotalMass(1.0);
	test.setFissionExtent(100);
	test.setFissionType(neutron);
	test.addIsotopeMaterial("U-235", 90);
	test.addIsotopeMaterial("U-238", 10);
	
	test.DisplayInfo();
	
	test.unloadNuclides();
	test.unloadFissionProductYields();
	time = clock() - time;
	std::cout << "\nRuntime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
