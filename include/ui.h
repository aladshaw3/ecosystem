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
#include "yaml_wrapper.h"

#ifndef UI_HPP_
#define UI_HPP_

std::string allLower(const std::string input);

bool exit(const std::string input);

void invalid_input(int count, int max);

int run_executable(int argc, const char * argv[]);

#endif
