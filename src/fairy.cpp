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

//set threshold
void FissionProducts::setThreshold(double val)
{
	if (val < 0.0)
		val = 0.0;
	this->DecayChain::setThreshold(val);
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

//Check mass error
void FissionProducts::checkFractions()
{
	double sum = 0.0;
	for (int i=0; i<this->MatFrac.size(); i++)
	{
		sum += this->MatFrac[i];
	}
	if (sum <= 0.0) {mError(invalid_molefraction); return;}
	if (sum-100.0 > 0.0)
	{
		double correction = 100.0/sum;
		for (int i=0; i<this->MatFrac.size(); i++)
		{
			this->MatFrac[i] = this->MatFrac[i]*correction;
		}
	}
}

//Evaluate the yields from data
int FissionProducts::evaluateYields()
{
	int success = 0;
	if (this->type == explosion)
	{
		std::string high_key;
		for (int i=0; i<this->InitialMat.size(); i++)
		{
			//Calculate moles of isotope and find header
			double moles;
			moles = this->total_mass*this->MatFrac[i]*1000.0/this->InitialMat[i].AtomicWeight();
			int levels = 0;
			try
			{
				levels = (int)this->fpy_data.getYamlWrapper()(this->InitialMat[i].IsotopeName()).getHeadMap().size();
			}
			catch (std::out_of_range)
			{
				mError(key_not_found);
				std::cout << this->InitialMat[i].IsotopeName() << " does not exist in yield data...\n\n";
				return -1;
			}
			
			//Loop through energy levels to find high key value
			double high = 0.0;
			for (auto &x: this->fpy_data.getYamlWrapper()(this->InitialMat[i].IsotopeName()).getHeadMap())
			{
				try
				{
					double energy = x.second["Energy"].getDouble();
					if (energy >= high)
					{
						high = energy;
						high_key = x.first;
					}
				}
				catch (std::out_of_range)
				{
					mError(key_not_found);
					std::cout << "Energy key does not exist in yield data...\n\n";
					return -1;
				}
			}
			
			//Read in yields for all isotopes from high_key
			try
			{
				levels = (int)this->fpy_data.getYamlWrapper()(this->InitialMat[i].IsotopeName())(high_key).getSubMap().size();
			}
			catch (std::out_of_range)
			{
				mError(key_not_found);
				std::cout << "Error in determining the High_Key value...\n\n";
				return -1;
			}
			
			//Loop through all isotopes produced in fission
			double yield = 0.0;
			double sum = 0.0;
			for (auto &x: this->fpy_data.getYamlWrapper()(this->InitialMat[i].IsotopeName())(high_key).getSubMap())
			{
				try
				{
					yield = x.second["Yield"].getDouble();
					success = this->registerInitialNuclide(x.first, yield*moles*this->fiss_extent/100.0);
					if (success != 0)
					{
						success = 0;
						sum += yield;
					}
				}
				catch (std::out_of_range)
				{
					mError(key_not_found);
					std::cout << "Missing Yield key value pair...\n\n";
					return -1;
				}
			}
			success = this->registerInitialNuclide(this->InitialMat[i].IsotopeName(), moles*(100.0-this->fiss_extent)/100.0);
			if (sum >= 1e-6)
			{
				mError(invalid_molefraction);
				std::cout << "Sum of independent yields contains too much error...\n\n";
				return -1;
			}
		}
	}
	else
	{
		//Need to do linear interpolation for neutron energy sources
	}
	this->createChains();
	this->formEigenvectors();
	//this->verifyEigenSoln();
	
	//this->DisplayList();
	return success;
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
	
	test.setTotalMass(9.7);
	test.setFissionExtent(99);
	test.setFissionType(explosion);
	test.setEnergyLevel(1000.0);
	test.setThreshold(3.0*log(0.5)/log(1-0.99)); //Set yeild to 3 sec cutoff based on 99% conversion to daughters
	test.addIsotopeMaterial("U-235", 90);
	test.addIsotopeMaterial("U-238", 10);
	test.checkFractions();
	
	test.DisplayInfo();
	
	test.loadFissionProductYields();
	test.evaluateYields();
	
	test.unloadNuclides();					//May be redundant
	test.unloadFissionProductYields();		//May be redundant
	time = clock() - time;
	std::cout << "\nRuntime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
