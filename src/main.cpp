//----------------------------------------
//  Created by Austin Ladshaw on 08/25/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "ui.h"

int main(int argc, const char * argv[])
{	
	int success = 0;
	
	//Command Line Interface: To Override, comment out this line and replace with your own run function
	success = run_executable(argc, argv);
	
	//Space below for overrides of the run_executable function
	//success = MACAW_TESTS();
	//success = SHARK_SCENARIO("shark_input.yml");
	//success = SHARK_TESTS();
	//success = SKUA_SCENARIOS("Scenario_Kr_Xe.txt", "AdsorbentProperties_FMOFZn.txt", "AllComponentProperties_N2_O2_Ar_CO2_Kr_Xe.txt", "AllAdsorbateProperties_Kr_Xe.txt");
	
	std::cout << "Exit Code:\t" << success << std::endl;
	return success;
}

