//----------------------------------------
//  Created by Austin Ladshaw on 08/25/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include <fstream>
#include <string>
#include <iostream>
#include "error.h"

#ifndef UI_HPP_
#define UI_HPP_

typedef enum
{
	TESTS,
	EXECUTABLES
} option;

int run_executable(int argc, const char * argv[]);

#endif
