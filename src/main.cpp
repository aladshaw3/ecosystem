/*!
 *  \file main.cpp
 *	\brief Main Function
 *	\details User input provided at time of execution is used to call the ui functions
 *  \author Austin Ladshaw
 *	\version 0.0 beta
 *	\date 08/25/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "ui.h"

int main(int argc, const char * argv[])
{	
	int success = 0;
	
	/// Command Line Interface: To Override, comment out this line and replace with your own run function
	success = run_executable(argc, argv);
	
	//Space below for overrides of the run_executable function
	//success = MACAW_TESTS();
	//success = SHARK_SCENARIO("shark_input.yml");
	//success = SHARK_TESTS();
	//success = SKUA_SCENARIOS("Scenario_Kr_Xe.txt", "AdsorbentProperties_FMOFZn.txt", "AllComponentProperties_N2_O2_Ar_CO2_Kr_Xe.txt", "AllAdsorbateProperties_Kr_Xe.txt");
	
	std::cout << "Exit Code:\t" << success << std::endl;
	return success;
}

