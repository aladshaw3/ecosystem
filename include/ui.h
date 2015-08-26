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
#include "flock.h"
#include "school.h"
#include "sandbox.h"
#include "Trajectory.h"

#ifndef UI_HPP_
#define UI_HPP_

typedef enum {TEST, EXECUTE, EXIT, CONTINUE,
				dogfish, eel, egret, finch, lark,
				macaw, mola, monkfish, sandbox,
				scopsowl, shark, skua} valid_options;

typedef struct
{
	ValueTypePair value_type;					//Data pair for input, tells what the input is and it's type
	std::string user_input;						//What is read in from the console at any point
	int count = 0;								//Number of times a questing has been asked
	int max = 3;								//Maximum allowable recursions of a question
	int option;									//Current option choosen by the user
	
}UI_DATA;

std::string allLower(const std::string input);

bool exit(const std::string input);

int invalid_input(int count, int max);

bool valid_input_main(UI_DATA *ui_dat);

bool valid_input_tests(UI_DATA *ui_dat);

int run_executable(int argc, const char * argv[]);

#endif