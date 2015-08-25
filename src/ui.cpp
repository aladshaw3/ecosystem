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
	if (copy == "exit")
		exit = true;
	return exit;
}

//Message to display for invalid input
void invalid_input(int count, int max)
{
	std::cout << "Invalid selection...\n\n";
	if (count < max-1)
	{
		std::cout << "Please try again...\n\n";
	}
	else
	{
		std::cout << "Max attempts reached! Force Quit...\n\n";
	}
}

//Run the executable based on the input arguments
int run_executable(int argc, const char * argv[])
{
	int success = 0;
	
	//Run a text based menu in command line
	if (argc == 1)
	{
		ValueTypePair value_type;
		std::string user_input;
		int count = 0;
		int max = 3;
		bool valid_input = false;
		do
		{
			std::cout << "What would you like to run?\n---------------------------\n";
			std::cout << "(" << 1 << ") TESTS\t(" << 2 << ") EXECUTABLES\n\nChoice: ";
			std::cin >> user_input;
			value_type.editValue(user_input);
			value_type.findType();
			
			if (exit(value_type.getValue()) == true)
			{
				std::cout << "Exiting program...\n\n";
				return 0;
			}
			
			if (value_type.getType() == INT)
			{
				if (value_type.getInt() == 1)
				{
					valid_input = true;
				}
				else if (value_type.getInt() == 2)
				{
					valid_input = true;
				}
				else
				{
					valid_input = false;
					invalid_input(count, max);
				}
			}
			else if (value_type.getType() == STRING)
			{
				if (allLower(value_type.getString()) == "tests")
				{
					valid_input = true;
				}
				else if (allLower(value_type.getString()) == "executables")
				{
					valid_input = true;
				}
				else
				{
					valid_input = false;
					invalid_input(count, max);
				}
			}
			else
			{
				valid_input = false;
				invalid_input(count, max);
			}
			count++;
		} while (valid_input == false && count < max);
	}
	
	return success;
}