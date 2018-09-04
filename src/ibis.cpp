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

/*
 *	-------------------------------------------------------------------------------------
 *								End: Ibis Class Definitions
 */
