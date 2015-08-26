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
	std::cout << "(" << 1 << ") TESTS\t(" << 2 << ") EXECUTABLES\n\nChoice: ";
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
	std::cout << "(1)  DOGFISH \t\t(2)  EEL      \t\t(3)  EGRET\n";
	std::cout << "(4)  FINCH   \t\t(5)  LARK     \t\t(6)  MACAW\n";
	std::cout << "(7)  MOLA    \t\t(8)  MONKFISH \t\t(9)  SANDBOX\n";
	std::cout << "(10) SCOPSOWL\t\t(11) SHARK    \t\t(12) SKUA\n";
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

//Run the executable based on the input arguments
int run_executable(int argc, const char * argv[])
{
	int success = 0;
	UI_DATA ui_dat;
	
	//Run a text based menu in command line
	if (argc == 1)
	{
		bool valid_input = false;
		do
		{
			valid_input = valid_input_main(&ui_dat);
			if (ui_dat.option == EXIT)
				return 0;
			
			ui_dat.count++;
		} while (valid_input == false && ui_dat.count < ui_dat.max);
		ui_dat.count = 0;
		valid_input = false;
		
		if (ui_dat.option == TEST)
		{
			do
			{
				valid_input = valid_input_tests(&ui_dat);
				if (ui_dat.option == EXIT)
					return 0;
				
			} while (valid_input == false && ui_dat.count < ui_dat.max);
			
			switch (ui_dat.option)
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
		}
		else if (ui_dat.option == EXECUTE)
		{
			mError(opt_no_support);
			return -1;
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