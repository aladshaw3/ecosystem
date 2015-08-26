//----------------------------------------
//  Created by Austin Ladshaw on 08/25/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "ui.h"

//Convert input to all lower case
std::string allLower(const std::string input)
{
	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	return copy;
}

//Check user string input for keyword "exit"
bool exit(const std::string input)
{
	bool exit = false;
	std::string copy = allLower(input);
	if (copy == "exit" || copy == "quit")
		exit = true;
	return exit;
}

//Check to see if user requests help
bool help(const std::string input)
{
	bool help = false;
	std::string copy = allLower(input);
	if (copy == "help")
		help = true;
	return help;
}

//Display a help message
void display_help(int argc)
{
	if (argc == 1)
	{
		std::cout << "This is the help message for basic ui\n";
	}
	else
	{
		std::cout << "This is the help message for console commands\n";
	}
	std::cout << std::endl;
}

//Message to display for invalid input
int invalid_input(int count, int max)
{
	std::cout << "Invalid selection!\n";
	if (count < max-1)
	{
		std::cout << "Please try again...\n\n";
		return CONTINUE;
	}
	else
	{
		std::cout << "Max attempts reached! Force Quit...\n\n";
		return EXIT;
	}
}

//Messages to display for option 1
bool valid_input_main(UI_DATA *ui_dat)
{
	bool valid_input = false;
	ui_dat->option = EXIT;
	
	std::cout << "What would you like to run?\n---------------------------\n";
	std::cout << "(" << 1 << ") TESTS    (" << 2 << ") EXECUTABLES\n\nChoice: ";
	std::cin >> ui_dat->user_input;
	std::cout << std::endl;
	ui_dat->value_type.editValue(ui_dat->user_input);
	ui_dat->value_type.findType();
	
	//Check for user to request exiting
	if (exit(ui_dat->value_type.getValue()) == true)
	{
		std::cout << "Exiting program...\n\n";
		ui_dat->option = EXIT;
		return true;
	}
	if (help(ui_dat->value_type.getValue()) == true)
	{
		display_help(ui_dat->argc);
		ui_dat->option = HELP;
		return true;
	}
	
	if (ui_dat->value_type.getType() == INT)
	{
		if (ui_dat->value_type.getInt() == 1)
		{
			ui_dat->option = TEST;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 2)
		{
			ui_dat->option = EXECUTE;
			valid_input = true;
		}
		else
		{
			valid_input = false;
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
		}
	}
	else if (ui_dat->value_type.getType() == STRING)
	{
		if (allLower(ui_dat->value_type.getString()) == "tests")
		{
			ui_dat->option = TEST;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "executables")
		{
			ui_dat->option = EXECUTE;
			valid_input = true;
		}
		else
		{
			valid_input = false;
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
		}
	}
	else
	{
		valid_input = false;
		ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
	}

	
	return valid_input;
}

//Ask for test options
bool valid_input_tests(UI_DATA *ui_dat)
{
	bool valid_input = false;
	ui_dat->option = EXIT;
	
	std::cout << "Choose a test to run from the list below\n----------------------------------------\n\n";
	std::cout << "(1)  DOGFISH     (2)  EEL          (3)  EGRET\n";
	std::cout << "(4)  FINCH       (5)  LARK         (6)  MACAW\n";
	std::cout << "(7)  MOLA        (8)  MONKFISH     (9)  SANDBOX\n";
	std::cout << "(10) SCOPSOWL    (11) SHARK        (12) SKUA\n";
	std::cout << "\nChoice: ";
	std::cin >> ui_dat->user_input;
	std::cout << std::endl;
	ui_dat->value_type.editValue(ui_dat->user_input);
	ui_dat->value_type.findType();
	
	//Check for user to request exiting
	if (exit(ui_dat->value_type.getValue()) == true)
	{
		std::cout << "Exiting program...\n\n";
		ui_dat->option = EXIT;
		return true;
	}
	if (help(ui_dat->value_type.getValue()) == true)
	{
		display_help(ui_dat->argc);
		ui_dat->option = HELP;
		return true;
	}
	
	if (ui_dat->value_type.getType() == INT)
	{
		if (ui_dat->value_type.getInt() == 1)
		{
			ui_dat->option = dogfish;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 2)
		{
			ui_dat->option = eel;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 3)
		{
			ui_dat->option = egret;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 4)
		{
			ui_dat->option = finch;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 5)
		{
			ui_dat->option = lark;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 6)
		{
			ui_dat->option = macaw;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 7)
		{
			ui_dat->option = mola;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 8)
		{
			ui_dat->option = monkfish;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 9)
		{
			ui_dat->option = sandbox;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 10)
		{
			ui_dat->option = scopsowl;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 11)
		{
			ui_dat->option = shark;
			valid_input = true;
		}
		else if (ui_dat->value_type.getInt() == 12)
		{
			ui_dat->option = skua;
			valid_input = true;
		}
		else
		{
			valid_input = false;
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
		}
	}
	else if (ui_dat->value_type.getType() == STRING)
	{
		if (allLower(ui_dat->value_type.getString()) == "dogfish")
		{
			ui_dat->option = dogfish;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "eel")
		{
			ui_dat->option = eel;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "finch")
		{
			ui_dat->option = finch;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "lark")
		{
			ui_dat->option = lark;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "macaw")
		{
			ui_dat->option = macaw;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "mola")
		{
			ui_dat->option = mola;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "monkfish")
		{
			ui_dat->option = monkfish;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "sandbox")
		{
			ui_dat->option = sandbox;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "scopsowl")
		{
			ui_dat->option = scopsowl;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "shark")
		{
			ui_dat->option = shark;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "skua")
		{
			ui_dat->option = skua;
			valid_input = true;
		}
		else
		{
			valid_input = false;
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
		}
	}
	else
	{
		valid_input = false;
		ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
	}
	
	return valid_input;
}

//Check executable options
bool valid_input_execute(UI_DATA *ui_dat)
{
	bool valid_input = false;
	ui_dat->option = EXIT;
	
	std::cout << "Choose a simulation to run from the list below\n-----------------------------------------------\n\n";
	std::cout << "(1)  GSTA_OPT      (2)  MAGPIE   (3)  SCOPSOWL\n";
	std::cout << "(4)  SCOPSOWL_OPT  (5)  SKUA     (6)  SKUA_OPT\n";
	std::cout << "\nChoice: ";
	std::cin >> ui_dat->user_input;
	std::cout << std::endl;
	ui_dat->value_type.editValue(ui_dat->user_input);
	ui_dat->value_type.findType();
	
	//Check for user to request exiting
	if (exit(ui_dat->value_type.getValue()) == true)
	{
		std::cout << "Exiting program...\n\n";
		ui_dat->option = EXIT;
		return true;
	}
	if (help(ui_dat->value_type.getValue()) == true)
	{
		display_help(ui_dat->argc);
		ui_dat->option = HELP;
		return true;
	}
	
	if (ui_dat->value_type.getType() == INT)
	{
		switch (ui_dat->value_type.getInt())
		{
			case 1:
			{
				ui_dat->option = gsta_opt;
				valid_input = true;
				break;
			}
				
			case 2:
			{
				ui_dat->option = magpie;
				valid_input = true;
				break;
			}
				
			case 3:
			{
				ui_dat->option = scopsowl;
				valid_input = true;
				break;
			}
				
			case 4:
			{
				ui_dat->option = scops_opt;
				valid_input = true;
				break;
			}
				
			case 5:
			{
				ui_dat->option = skua;
				valid_input = true;
				break;
			}
				
			case 6:
			{
				ui_dat->option = skua_opt;
				valid_input = true;
				break;
			}
				
			default:
			{
				valid_input = false;
				ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
				break;
			}
		}
	}
	else if (ui_dat->value_type.getType() == STRING)
	{
		if (allLower(ui_dat->value_type.getString()) == "gsta_opt")
		{
			ui_dat->option = gsta_opt;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "magpie")
		{
			ui_dat->option = magpie;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "scopsowl")
		{
			ui_dat->option = scopsowl;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "scopsowl_opt")
		{
			ui_dat->option = scops_opt;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "skua")
		{
			ui_dat->option = skua;
			valid_input = true;
		}
		else if (allLower(ui_dat->value_type.getString()) == "skua_opt")
		{
			ui_dat->option = skua_opt;
			valid_input = true;
		}
		else
		{
			valid_input = false;
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
		}
	}
	else
	{
		valid_input = false;
		ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
	}
	
	return valid_input;
}

//Run the test case declared by user
int run_test(UI_DATA *ui_dat)
{
	int success = 0;
	
	switch (ui_dat->option)
	{
		case dogfish:
			success = DOGFISH_TESTS();
			break;
			
		case eel:
			success = EEL_TESTS();
			break;
			
		case egret:
			success = EGRET_TESTS();
			break;
			
		case finch:
			success = FINCH_TESTS();
			break;
			
		case lark:
			success = LARK_TESTS();
			break;
			
		case macaw:
			success = MACAW_TESTS();
			break;
			
		case mola:
			success = MOLA_TESTS();
			break;
			
		case monkfish:
			success = MONKFISH_TESTS();
			break;
			
		case sandbox:
			success = RUN_SANDBOX();
			break;
			
		case scopsowl:
			success = SCOPSOWL_TESTS();
			break;
			
		case shark:
			success = SHARK_TESTS();
			break;
			
		case skua:
			success = SKUA_TESTS();
			break;
			
		default:
		{
			mError(opt_no_support);
			return -1;
			break;
		}
	}
	
	return success;
}

//Run the execution option choosen by the user
int run_exec(UI_DATA *ui_dat)
{
	int success = 0;
	
	std::cout << "WARNING! All simulation runs will require properly formatted input files.\n";
	std::cout << "User is required to provide full path and extension to each necessary file.\n";
	std::cout << "(Example: this/is/path/to/input/files.txt)\n\n";
	
	switch (ui_dat->option)
	{
		case scopsowl:
		{
			ui_dat->input_files.resize(4);
			std::cout << "SCOPSOWL needs 4 input files containing scenario, adsorbent, component, and adsorbate information.\n";
			std::cout << "Please provide each necessary file and full path for the requested file...\n\n";
			std::cout << "Scenario file: ";
			std::cin >> ui_dat->input_files[0];
			std::cout << std::endl;
			std::cout << "Adsorbent Properties file: ";
			std::cin >> ui_dat->input_files[1];
			std::cout << std::endl;
			std::cout << "Component Properties file: ";
			std::cin >> ui_dat->input_files[2];
			std::cout << std::endl;
			std::cout << "Adsorbate Properties file: ";
			std::cin >> ui_dat->input_files[3];
			std::cout << std::endl;
			
			success = SCOPSOWL_SCENARIOS(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str());
			break;
		}
			
		case skua:
		{
			ui_dat->input_files.resize(4);
			std::cout << "SKUA needs 4 input files containing scenario, adsorbent, component, and adsorbate information.\n";
			std::cout << "Please provide each necessary file and full path for the requested file...\n\n";
			std::cout << "Scenario file: ";
			std::cin >> ui_dat->input_files[0];
			std::cout << std::endl;
			std::cout << "Adsorbent Properties file: ";
			std::cin >> ui_dat->input_files[1];
			std::cout << std::endl;
			std::cout << "Component Properties file: ";
			std::cin >> ui_dat->input_files[2];
			std::cout << std::endl;
			std::cout << "Adsorbate Properties file: ";
			std::cin >> ui_dat->input_files[3];
			std::cout << std::endl;
			
			success = SKUA_SCENARIOS(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str());
			break;
		}
			
		case gsta_opt:
		{
			ui_dat->input_files.resize(1);
			std::cout << "GSTA_OPT needs 1 input file containing the equilibrium adsorption data.\n";
			std::cout << "Please provide the file and full path...\n\n";
			std::cout << "GSTA input: ";
			std::cin >> ui_dat->input_files[0];
			std::cout << std::endl;
			success = gsta_optimize(ui_dat->input_files[0].c_str());
			break;
		}
			
		case magpie:
		{
			ui_dat->input_files.resize(2);
			std::cout << "MAGPIE needs 2 input files containing gsta isotherm equilibrium parameters and scenario information.\n";
			std::cout << "Please provide each necessary file and full path for the requested file...\n\n";
			std::cout << "Isotherm Parameters file: ";
			std::cin >> ui_dat->input_files[0];
			std::cout << std::endl;
			std::cout << "Scenario file: ";
			std::cin >> ui_dat->input_files[1];
			std::cout << std::endl;
			
			success = MAGPIE_SCENARIOS(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str());
			break;
		}
			
		case scops_opt:
		{
			ui_dat->input_files.resize(5);
			std::cout << "SCOPSOWL_OPT needs 5 input files containing scenario, adsorbent, component, and adsorbate information,\n";
			std::cout << "as well as the data file containing the experimental data to compare the model against.\n";
			std::cout << "Please provide each necessary file and full path for the requested file...\n\n";
			std::cout << "Scenario file: ";
			std::cin >> ui_dat->input_files[0];
			std::cout << std::endl;
			std::cout << "Adsorbent Properties file: ";
			std::cin >> ui_dat->input_files[1];
			std::cout << std::endl;
			std::cout << "Component Properties file: ";
			std::cin >> ui_dat->input_files[2];
			std::cout << std::endl;
			std::cout << "Adsorbate Properties file: ";
			std::cin >> ui_dat->input_files[3];
			std::cout << std::endl;
			std::cout << "Experimental Data file: ";
			std::cin >> ui_dat->input_files[4];
			std::cout << std::endl;
			
			success = SCOPSOWL_OPTIMIZE(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str(), ui_dat->input_files[4].c_str());
			break;
		}
			
		case skua_opt:
		{
			ui_dat->input_files.resize(5);
			std::cout << "SKUA_OPT needs 5 input files containing scenario, adsorbent, component, and adsorbate information,\n";
			std::cout << "as well as the data file containing the experimental data to compare the model against.\n";
			std::cout << "Please provide each necessary file and full path for the requested file...\n\n";
			std::cout << "Scenario file: ";
			std::cin >> ui_dat->input_files[0];
			std::cout << std::endl;
			std::cout << "Adsorbent Properties file: ";
			std::cin >> ui_dat->input_files[1];
			std::cout << std::endl;
			std::cout << "Component Properties file: ";
			std::cin >> ui_dat->input_files[2];
			std::cout << std::endl;
			std::cout << "Adsorbate Properties file: ";
			std::cin >> ui_dat->input_files[3];
			std::cout << std::endl;
			std::cout << "Experimental Data file: ";
			std::cin >> ui_dat->input_files[4];
			std::cout << std::endl;
			
			success = SKUA_OPTIMIZE(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str(), ui_dat->input_files[4].c_str());
			break;
		}
			
		default:
		{
			mError(opt_no_support);
			return -1;
			break;
		}
	}
	
	return success;
}

//Run the executable based on the input arguments
int run_executable(int argc, const char * argv[])
{
	int success = 0;
	UI_DATA ui_dat;
	ui_dat.argc = argc;
	*ui_dat.argv = *argv;
	
	//Run a text based menu in command line
	if (argc == 1)
	{
		bool valid_input = false;
		do
		{
			valid_input = valid_input_main(&ui_dat);
			if (ui_dat.option == EXIT)
				return 0;
			if (ui_dat.option == HELP)
			{
				ui_dat.count--;
				valid_input = false;
			}
			ui_dat.count++;
		} while (valid_input == false && ui_dat.count < ui_dat.max);
		ui_dat.count = 0;
		valid_input = false;
		
		if (ui_dat.option == TEST)
		{
			ui_dat.exec = false;
			do
			{
				valid_input = valid_input_tests(&ui_dat);
				if (ui_dat.option == EXIT)
					return 0;
				if (ui_dat.option == HELP)
				{
					ui_dat.count--;
					valid_input = false;
				}
				ui_dat.count++;
			} while (valid_input == false && ui_dat.count < ui_dat.max);
			
			success = run_test(&ui_dat);
		}
		else if (ui_dat.option == EXECUTE)
		{
			ui_dat.exec = true;
			do
			{
				valid_input = valid_input_execute(&ui_dat);
				if (ui_dat.option == EXIT)
					return 0;
				if (ui_dat.option == HELP)
				{
					ui_dat.count--;
					valid_input = false;
				}
				ui_dat.count++;
			} while (valid_input == false && ui_dat.count < ui_dat.max);
			
			success = run_exec(&ui_dat);
		}
		else
		{
			mError(opt_no_support);
			return -1;
		}
	}
	else
	{
		mError(opt_no_support);
		return -1;
	}
	
	return success;
}