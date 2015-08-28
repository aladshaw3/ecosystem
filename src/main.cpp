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
	
	std::cout << "Exit Code:\t" << success << std::endl;
	return success;
}

