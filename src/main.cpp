//
//  main.cpp
//  AdsorptionToolBox
//
//  Created by aladshaw3 on 2/23/15.
//  Copyright (c) 2015 Georgia Institute of Technology. All rights reserved.
//

#include "ui.h"

int main(int argc, const char * argv[])
{	
	int success = 0;
	
	//------------------------------Command Line Interface------------------------
	
	success = run_executable(argc, argv);
	
	//std::cout << argv[0] << std::endl;	//Name of executable with path
	/*
	if (argc > 1)						//Next array of arguments
		std::cout << argv[1] << std::endl;
	if (argc > 2)						//Next array of arguments
		std::cout << argv[2] << std::endl;
	
	//Assume argv[2] is a file to open
	if (argc > 2)
	{
		// We assume argv[2] is a filename to open
		std::ifstream the_file ( argv[2] );
		// Always check to see if file opening succeeded
		if ( !the_file.is_open() )
			std::cout<<"Could not open file\n";
		else
		{
			char x;
			// the_file.get ( x ) returns false if the end of the file
			//  is reached or an error occurs
			while ( the_file.get ( x ) )
				std::cout << x << std::endl;
			
			//Note this will get the file even if said file is in another location 
		}
	}
	*/
	
	std::cout << "Exit Code:\t" << success << std::endl;
	return success;
}

