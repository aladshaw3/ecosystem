/*!
 *  \file ibis.h ibis.cpp
 *	\brief Implicit Branching Isotope System
 *  \author Austin Ladshaw
 *	\date 09/04/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay. Copyright (c) 2018, all
 *             rights reserved.
 */

#include "ibis.h"

/// Function to convert from a starting unit and value to and ending unit and value (returns converted value)
double time_conversion(time_units end_unit, double start_value, time_units start_unit)
{
	double end_value = 0.0;
	
	switch (end_unit)
	{
		case seconds:
			switch (start_unit)
			{
				case seconds:
					end_value = start_value;
					break;
					
				case minutes:
					end_value = start_value * 60.0;
					break;
					
				case hours:
					end_value = start_value * 60.0 * 60.0;
					
				case days:
					end_value = start_value * 60.0 * 60.0 * 24.0;
					
				case years:
					end_value = start_value * 60.0 * 60.0 * 24.0 * 365.25;
					
				default:
					end_value = start_value;
					break;
			}
			break;
			
		case minutes:
			switch (start_unit)
		{
			case seconds:
				end_value = start_value / 60.0;
				break;
				
			case minutes:
				end_value = start_value;
				break;
				
			case hours:
				end_value = start_value * 60.0;
				
			case days:
				end_value = start_value * 60.0 * 24.0;
				
			case years:
				end_value = start_value * 60.0 * 24.0 * 365.25;
				
			default:
				end_value = start_value;
				break;
		}
			break;
			
		case hours:
			switch (start_unit)
		{
			case seconds:
				end_value = start_value / 60.0 / 60.0;
				break;
				
			case minutes:
				end_value = start_value / 60.0;
				break;
				
			case hours:
				end_value = start_value;
				
			case days:
				end_value = start_value * 24.0;
				
			case years:
				end_value = start_value * 24.0 * 365.25;
				
			default:
				end_value = start_value;
				break;
		}
			break;
			
		case days:
			switch (start_unit)
		{
			case seconds:
				end_value = start_value / 60.0 / 60.0 / 24.0;
				break;
				
			case minutes:
				end_value = start_value / 60.0 / 24.0;
				break;
				
			case hours:
				end_value = start_value / 24.0;
				
			case days:
				end_value = start_value;
				
			case years:
				end_value = start_value * 365.25;
				
			default:
				end_value = start_value;
				break;
		}
			break;
			
		case years:
			switch (start_unit)
		{
			case seconds:
				end_value = start_value / 60.0 / 60.0 / 24.0 / 365.25;
				break;
				
			case minutes:
				end_value = start_value / 60.0 / 24.0 / 365.25;
				break;
				
			case hours:
				end_value = start_value / 24.0 / 365.25;
				
			case days:
				end_value = start_value / 365.25;
				
			case years:
				end_value = start_value;
				
			default:
				end_value = start_value;
				break;
		}
			break;
			
		default:
			end_value = start_value;
			break;
	}
	
	return end_value;
}

/*
 *								Start: Ibis Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
Isotope::Isotope() : Atom(0)
{
	IsoName = "No Name";
	decay_rate = 0.0;
	half_life = INFINITY;
	hl_units = years;
	isotope_number = 0;
}

//Default destructor
Isotope::~Isotope()
{
	decay_modes.clear();
	branch_ratios.clear();
}

//Register via atomic number and isotop number
void Isotope::registerIsotope(int atom_num, int iso_num)
{
	this->Register(atom_num);
	this->isotope_number = iso_num;
	if (this->Neutrons() + atom_num != iso_num)
		this->editAtomicWeight(this->AtomicWeight() + (iso_num - atom_num));
	this->editNeutrons(iso_num - atom_num);
	
	this->setConstants();
	this->computeDecayRate();
}

//Return isotope number
int Isotope::IsotopeNumber()
{
	return this->isotope_number;
}

//Return decay rate
double Isotope::DecayRate()
{
	return this->decay_rate;
}

//Return half-life
double Isotope::HalfLife(time_units units)
{
	return time_conversion(units, this->half_life, this->hl_units);
}

//Return name
std::string Isotope::IsotopeName()
{
	return this->IsoName;
}

//Set the decay information based on a registered atomic number and isotope number
void Isotope::setConstants()
{
	// n (neutron)
	if (this->AtomicNumber() == 0)
	{
		//Register decay modes and constants
		if (this->isotope_number == 1)
		{
			this->hl_units = seconds;
			this->half_life = 613.9;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}

		//Unregistered isotope
		else
		{
			mError(invalid_isotope);
		}
	}
	
	// H atom
	else if (this->AtomicNumber() == 1)
	{
		//Register decay modes and constants
		if (this->isotope_number == 1)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 2)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 3)
		{
			this->hl_units = years;
			this->half_life = 12.32;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		
		//Unregistered isotope
		else
		{
			mError(invalid_isotope);
		}
	}
	
	// He atom
	else if (this->AtomicNumber() == 2)
	{
		//Register decay modes and constants
		if (this->isotope_number == 2)
		{
			this->hl_units = seconds;
			this->half_life = 1E-9;
			this->decay_modes.resize(2);
			this->branch_ratios.resize(2);
			
			this->decay_modes[0] = prot_emiss;
			this->branch_ratios[0] = 99.99/100.0;
			
			this->decay_modes[1] = beta_plus;
			this->branch_ratios[1] = 0.01/100.0;
		}
		else if (this->isotope_number == 3)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 4)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 5)
		{
			this->hl_units = seconds;
			this->half_life = 7.04E-22;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = neut_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 6)
		{
			this->hl_units = seconds;
			this->half_life = 0.807;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 7)
		{
			this->hl_units = seconds;
			this->half_life = 3.04E-21;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = neut_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 8)
		{
			this->hl_units = seconds;
			this->half_life = 0.119;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 9)
		{
			this->hl_units = seconds;
			this->half_life = 7E-21;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = neut_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		
		//Unregistered isotope
		else
		{
			mError(invalid_isotope);
		}
	}

	// Li atom
	if (this->AtomicNumber() == 3)
	{
		//Register decay modes and constants
		if (this->isotope_number == 4)
		{
			this->hl_units = seconds;
			this->half_life = 91E-24;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = prot_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 5)
		{
			this->hl_units = seconds;
			this->half_life = 3.71E-22;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = prot_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 6)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 7)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 8)
		{
			this->hl_units = seconds;
			this->half_life = 0.84;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 9)
		{
			this->hl_units = seconds;
			this->half_life = 0.178;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 10)
		{
			this->hl_units = seconds;
			this->half_life = 2E-21;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = neut_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 11)
		{
			this->hl_units = seconds;
			this->half_life = 8.75E-3;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		
		//Unregistered isotope
		else
		{
			mError(invalid_isotope);
		}
	}
	
	// Be Atom
	if (this->AtomicNumber() == 4)
	{
		//Register decay modes and constants
		if (this->isotope_number == 7)
		{
			this->hl_units = days;
			this->half_life = 53.22;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = elec_cap;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 8)
		{
			this->hl_units = seconds;
			this->half_life = 8.19E-17;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = alpha;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 9)
		{
			this->hl_units = years;
			this->half_life = INFINITY;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = stable;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 10)
		{
			this->hl_units = years;
			this->half_life = 1.39E6;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 11)
		{
			this->hl_units = seconds;
			this->half_life = 13.81;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 12)
		{
			this->hl_units = seconds;
			this->half_life = 0.02149;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = beta_minus;
			this->branch_ratios[0] = 100.0/100.0;
		}
		else if (this->isotope_number == 13)
		{
			this->hl_units = seconds;
			this->half_life = 2.7E-21;
			this->decay_modes.resize(1);
			this->branch_ratios.resize(1);
			
			this->decay_modes[0] = neut_emiss;
			this->branch_ratios[0] = 100.0/100.0;
		}
		
		//Unregistered isotope
		else
		{
			mError(invalid_isotope);
		}
	}
	
	//Unregistered Atom
	else
	{
		//No action
		mError(invalid_atom);
	}
}

//Compute decay rate
void Isotope::computeDecayRate()
{
	
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Ibis Class Definitions
 */

//Test function
int IBIS_TESTS()
{
	int success = 0;
	
	return success;
}
