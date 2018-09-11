/*!
 *  \file ui.h ui.cpp
 *	\brief User Interface for Ecosystem
 *	\details These routines define how the user will interface with the software
 *  \author Austin Ladshaw
 *	\date 08/25/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all 
 *             rights reserved.
 */

/*!	\mainpage Introduction
 *
 *	\section copyright Copyright Statement
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 *
 *	\section info General Information
 *
 *	The source code contained within the ecosystem project was designed as a
 *	standalone tool set for performing numerical modeling and data analyses
 *	associated with adsorption phenomena in both gaseous and aqueous systems.
 *	Many of the lower level tools are general enough to be applied to any 
 *	system you desire to be modeled. Such algorithms included are Krylov
 *	subspace methods for linear systems and a Jacobian-Free Newton-Krylov
 *	method for non-linear systems. There is also a templated matrix object
 *	for generic matrix construction and modification. For specific information
 *	on each individual kernel, navigate through the class and file indices or
 *	table of contents.
 *
 *	\warning Many of these algorithms may still be under development. This library is 
 *	distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
 *	even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

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

/// Macro expansion for executable current version number
#define ECO_VERSION "1.0.0"
/// Macro expansion for executable current name
#define ECO_EXECUTABLE "eco"

/// Valid options available upon execution of the code.
/** Enumeration of valid options for executing the ecosystem code.
    More options become available as the code updates. Some options that
    appear here may not be viewable in the "help" screen of the executable.
    Those options are hidden, but are still valid entries.*/
typedef enum {TEST, EXECUTE, EXIT, CONTINUE, HELP,
	
				dogfish, eel, egret, finch, lark,
				macaw, mola, monkfish, sandbox,
	
				scopsowl, shark, skua, gsta_opt, magpie,
				scops_opt, skua_opt, trajectory, dove,
				crow, mesh, crane, ibis } valid_options;

/// Data structure holding the UI arguments.
/** C-Style object for interfacing with users request upon execution of
    the program. User input is stored in objects below and a series of 
	booleans is used to determine how and what to execute. */
typedef struct
{
	ValueTypePair value_type;					///< Data pair for input, tells what the input is and it's type
	std::vector<std::string> user_input;		///< What is read in from the console at any point
	std::vector<std::string> input_files;		///< A vector of input file names and directories given by user
	std::string path;							///< Path to where input files are located
	int count = 0;								///< Number of times a questing has been asked
	int max = 3;								///< Maximum allowable recursions of a question
	int option;									///< Current option choosen by the user
	bool Path = false;							///< True if user gives path as an option
	bool Files = false;							///< True if user gives input files as an option
	bool MissingArg = true;						///< True if an input argument is missing; False if everything is ok
	bool BasicUI = true;						///< True if using Basic UI; False if using Advanced UI
	int argc;									///< Number of console arguments given on input
	const char * argv[];						///< Actual console arguments given at execution
	
}UI_DATA;

/// Function to display help for Advanced User Interface
/** The Advanved User Interface help screen is accessed by including
	run option -h or --help when executing the program from command line. */
void aui_help();

/// Function to display help for Basic User Interface
/** The Basic User Interface help screen is accessed by running the executable,
	then typing "help" at any point during the console prompts. Exception to this
	occurs when the console prompts you to provide input files for your choosen routine.
	In this circumstance, the executable always assumes that what the user types in will
	be an input file. */
void bui_help();

/// Function returns true if user requests exit
/** This function will check the input string for "exit" or "quit" and 
	terminate the executable. Only checked if using the Basic User
	Interface. 
 
	\param input input string user gives to the console */
bool exit(const std::string &input);

/// Function returns trun if the user requests help
/** This function will check the input string for "help", "-h", or "--help"
	and will tell the executable to display the help menu. The help menu that
	gets displayed depends on how the executable was run to begin with. 
 
	\param input input string user gives to the console*/
bool help(const std::string &input);

/// Function returns true if user requests to know the executable version
/** This function will check the input string for "version", "-v", or "--version"
	and will tell the executable to display version information about the executable.
 
	\param input input string user gives to the console */
bool version(const std::string &input);

/// Function returns true if user requests to run a test
/** This function will check the input string for "-t" or "--test" and
	determine whether or not the user requests to run an ecosystem test function.
 
	\param input input string user gives to the console*/
bool test(const std::string &input);

/// Function returns true if the user requests to run a simulation/executable
/** This function will check the input string for "-e" or "--execute" and 
	determine whether or not the user requests to run an ecosystem executable function. 
 
	\param input input string the user gives to the console*/
bool exec(const std::string &input);

/// Function returns true if the user indicates that input files share a common path
/** This function will check the input string for "-p" or "--path" and
	determine whether or not the user will give a common path to all input files
	needed for the specified simulation. Only used in Advanced User Interface.
 
	\param input input string the user gives to the console*/
bool path(const std::string &input);

/// Function returns true if the user indicates that the next arguments are input files
/** This function will check the input string for "-i" or "--input" and
	determine whether or not the user's next arguments are input files for 
	a specific simulation. Only used in Advanced User Interface.
 
	\param input input string the user gives to the console*/
bool input(const std::string &input);

/// Function returns true if the user gave a valid test option
/** This function will check the input string given by the user and determine whether
	that string denotes a valid test. Then, it will mark the option variable in ui_dat
	with the appropriate option from the valid_options enum. 
	
	\param input input string the user gives to the console
	\param ui_dat pointer to the data structure for the ui object*/
bool valid_test_string(const std::string &input, UI_DATA *ui_dat);

/// Function returns true if the user gave a valid execution option
/** This function will check the input string given by the user and determine whether
	that string denotes a valid execution option. Then, it will mark the option variable in ui_dat
	with the appropriate option from the valid_options enum.
	
	\param input input string the user gives to the console
	\param ui_dat pointer to the data structure for the ui object*/
bool valid_exec_string(const std::string &input, UI_DATA *ui_dat);

/// Function returns the number of expected input files for the user's run option
/** This function will check the option variable in the ui_dat structure to determine
	the number of input files that is expected to be given. Running different executable
	functions in ecosystem may require various number of input files.
 
	\param ui_dat pointer to the data structure for the ui object*/
int number_files(UI_DATA *ui_dat);

/// Function returns true if the user has choosen a valid additional runtime option
/** This function will check all additional input options in the user_input variable
	of ui_dat to determine if the user requests any additional options during runtime.
	Valid additional options are -p or --path and -i or --input. 
 
	\param ui_dat pointer to the data structure for the ui object*/
bool valid_addon_options(UI_DATA *ui_dat);

/// Function to call the appropriate help menu based on type of interface
/** This function looks at the ui_dat structure and the user's OS files to determine 
	what help menu to display and how to display it. There are two different types of
	help menus that can be displayed: (i) Advanced Help and (ii) Basic Help. Additionally,
	this function checks the OS file system for the existence of installed help files. If it
	finds those files, then it instructs the command terminal to read the contents of those
	files with the "less" command. Otherwise, it will just print the appropriate help menu to
	the console window.
 
	\param ui_dat pointer to the data structure for the ui object*/
void display_help(UI_DATA *ui_dat);

/// Function to display ecosystem version information to the console
/** This function will check the ui_dat structure to see which type of interface the 
	user is using, then print out the version information for the executable being run.
 
	\param ui_dat pointer to the data structure for the ui object*/
void display_version(UI_DATA *ui_dat);

/// Function returns a CONTINUE or EXIT when invalid input is given
/** This function looks at the current count and the max iterations and determines
	whether or not to force the executable to terminate. If the user provides too many
	incorrect options during the Basic User Interface, then the executable will force quit.
 
	\param count number of times the user has provided a bad option
	\param max maximum allowable bad options before force quit*/
int invalid_input(int count, int max);

/// Function returns true if user gave valid input in Basic UI
/** This function is only called if the user is running the Basic UI. It checks the
	given console argument stored in user_input of ui_dat for a valid option. If no
	valid option is given, then this function returns false. 
 
	\param ui_dat pointer to the data structure for the ui object*/
bool valid_input_main(UI_DATA *ui_dat);

/// Function returns true if user gave a valid test function to run
/** This function checks the user_input argument of ui_dat for a valid test option. If
	no valid test was given, then this function returns false.
 
	\param ui_dat pointer to the data structure for the ui object*/
bool valid_input_tests(UI_DATA *ui_dat);

/// Function returns true if user gave a valid executable function to run
/** This function checks the user_input argument of ui_dat for a valid executable option. If
	no valid executable was given, then this function returns false.
 
	\param ui_dat pointer to the data structure for the ui object*/
bool valid_input_execute(UI_DATA *ui_dat);

/// Function that loops the Basic UI until a valid test option was selected
/** This function loops the Basic UI menu for running a test until a valid test is
	selected by the user. If a valid test is not selected, and the maximum number of
	loops has been reached, then this function will cause the program to force quit. 
 
	\param ui_dat pointer to the data structure for the ui object*/
int test_loop(UI_DATA *ui_dat);

/// Function that loops the Basic UI until a valid executable option was selected
/** This function loops the Basic UI menu for running an executable until a valid executable is
	selected by the user. If a valid executable is not selected, and the maximum number of
	loops has been reached, then this function will cause the program to force quit.
 
	\param ui_dat pointer to the data structure for the ui object*/
int exec_loop(UI_DATA *ui_dat);

/// Function will call the user requested test function
/** This function checks the option variable of the ui_dat structure and runs the 
	corresponding test function.
 
	\param ui_dat pointer to the data structure for the ui object*/
int run_test(UI_DATA *ui_dat);

/// Function will call the user requested executable function
/** This function checks the option variable of the ui_dat structure and runs the
	corresponding executable function.
 
	\param ui_dat pointer to the data structure for the ui object*/
int run_exec(UI_DATA *ui_dat);

/// Function called by the main and runs both user interfaces for the program
/** This function is called in the main.cpp file and passes the console
	arguments given at run time.
 
	\param argc number of arguments provided by the user at the time of execution
	\param argv list of C-strings that was provided by the user at the time of execution*/
int run_executable(int argc, const char * argv[]);

#endif
