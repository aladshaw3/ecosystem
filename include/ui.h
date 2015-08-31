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

#define ECO_VERSION "0.0 alpha"
#define ECO_EXECUTABLE "eco0"

typedef enum {TEST, EXECUTE, EXIT, CONTINUE, HELP,
	
				dogfish, eel, egret, finch, lark,
				macaw, mola, monkfish, sandbox,
	
				scopsowl, shark, skua, gsta_opt, magpie,
				scops_opt, skua_opt, trajectory} valid_options;

typedef struct
{
	ValueTypePair value_type;					//Data pair for input, tells what the input is and it's type
	std::vector<std::string> user_input;		//What is read in from the console at any point
	std::vector<std::string> input_files;		//A vector of input file names and directories given by user
	std::string path;							//Path to where input files are located
	int count = 0;								//Number of times a questing has been asked
	int max = 3;								//Maximum allowable recursions of a question
	int option;									//Current option choosen by the user
	bool Path = false;							//True if user gives path as an option
	bool Files = false;							//True if user gives input files as an option
	bool MissingArg = true;						//True if an input argument is missing; False if everything is ok
	bool BasicUI = true;						//True if using Basic UI; False if using Advanced UI
	int argc;									//Number of console arguments given on input
	const char * argv[];						//Actual console arguments given at execution
	
}UI_DATA;

std::string allLower(const std::string &input);

bool exit(const std::string &input);

bool help(const std::string &input);

bool version(const std::string &input);

bool test(const std::string &input);

bool exec(const std::string &input);

bool path(const std::string &input);

bool input(const std::string &input);

bool valid_test_string(const std::string &input, UI_DATA *ui_dat);

bool valid_exec_string(const std::string &input, UI_DATA *ui_dat);

int number_files(UI_DATA *ui_dat);

bool valid_addon_options(UI_DATA *ui_dat);

void exec_option(UI_DATA *ui_dat);

void display_help(UI_DATA *ui_dat);

void display_version(UI_DATA *ui_dat);

int invalid_input(int count, int max);

bool valid_input_main(UI_DATA *ui_dat);

bool valid_input_tests(UI_DATA *ui_dat);

bool valid_input_execute(UI_DATA *ui_dat);

int test_loop(UI_DATA *ui_dat);

int exec_loop(UI_DATA *ui_dat);

int run_test(UI_DATA *ui_dat);

int run_exec(UI_DATA *ui_dat);

int run_executable(int argc, const char * argv[]);

#endif
