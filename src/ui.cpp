//----------------------------------------
//  Created by Austin Ladshaw on 08/25/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "ui.h"

//Convert input to all lower case
std::string allLower(const std::string &input)
{
	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	return copy;
}

//Check user string input for keyword "exit"
bool exit(const std::string &input)
{
	bool exit = false;
	std::string copy = allLower(input);
	if (copy == "exit" || copy == "quit")
		exit = true;
	return exit;
}

//Check to see if user requests help
bool help(const std::string &input)
{
	bool help = false;
	std::string copy = allLower(input);
	if (copy == "help" || copy == "-h" || copy == "--help")
		help = true;
	return help;
}

//Check to see if user requests version
bool version(const std::string &input)
{
	bool version = false;
	std::string copy = allLower(input);
	if (copy == "-v" || copy == "--version" || copy == "version")
		version = true;
	return version;
}

//Check to see if user wants to run a test
bool test(const std::string &input)
{
	bool test = false;
	std::string copy = allLower(input);
	if (copy == "-t" || copy == "--test")
		test = true;
	return test;
}

//Check to see if user wants to run an executable
bool exec(const std::string &input)
{
	bool exec = false;
	std::string copy = allLower(input);
	if (copy == "-e" || copy == "--execute")
		exec = true;
	return exec;
}

//Check to see if user wants to specify a single input path
bool path(const std::string &input)
{
	bool path = false;
	std::string copy = allLower(input);
	if (copy == "-p" || copy == "--path")
		path = true;
	return path;
}

//Check to see if user is giving input files
bool input(const std::string &input)
{
	bool file = false;
	std::string copy = allLower(input);
	if (copy == "-i" || copy == "--input")
		file = true;
	return file;
}

//Display to console the current test info
void current_tests()
{
	std::cout << "           CURRENTLY AVAILABLE TESTS           \n-----------------------------------------------\n\n";
	
	std::cout << "(1) DOGFISH (Diffusion Object Governing Fiber Interior Sorption History)\n\n";
	puts("\tThis test runs a simple example of the intraparticle mass transfer uptake of aqueous ions into cylindrical adsorbent fibers. Currently, there is no executable for this set of algorithms.\n\n");
	
	std::cout << "(2) EEL (Easy-access Element Library)\n\n";
	puts("\tThis test runs a series of checks on our digital atom library to ensure that all objects are operating as they should. There is no output file associated with this test and these algorithms are not used directly by the user, but are called by other algorithms in the library.\n\n");
	
	std::cout << "(3) EGRET (Estimation of Gas-phase pRopErTies)\n\n";
	puts("\tThis test runs a series of checks on our implementations of kinetic gas theory to predict various gas phase properties from the basic molecular information of each molecule in a gas. Properties calculated include binary diffusivities, molecular diffusivities, and film mass transfer coefficients for each individual gas species, as well as determining the overall gas viscosity, density, and heat capacity. There is no output file associated with this test and these algorithms are not used directly by the user, but are called by other algorithms in the library.\n\n");
	
	std::cout << "(4) FINCH (Flux-limiting Implicit Non-oscillatory Conservative High-resolution scheme)\n\n";
	puts("\tThis test runs an example calculation for a 1-D PDE representing a conservation law. Our algorithms are based on a MUSCL scheme for high accuracy solutions to PDEs involving highly advective components. There is an output file associated with this test. However, these algorithms are not used directly by the user, but are called by other algorithms in the library when PDE solutions are needed.\n\n");
	
	std::cout << "(5) LARK (Linear Algebra Residual Kernels)\n\n";
	puts("\tThis test runs a series of checks on our implementations of various iterative solvers to systems of equations. The available solvers include multiple Krylov Subpace techniques, such as PCG, GMRES, CGS, BiCGSTAB, and GMRESR, as well as two different non-linear solvers: Picard's Method and the Jacobian-Free Newton-Krylov method. We also have an implementation of the Arnoldi Iteration to produce a full orthonormal basis from any non-singular matrix operator. There is no output file associated with this test and these algorithms are not used directly by the user, but are called by other algorithms in the library.\n\n");
	
	std::cout << "(6) MACAW (MAtrix CAlculation Workspace)\n\n";
	puts("\tThis test runs a series of checks on our Matrix template object. This object is use extensively throughout the entire library. Therefore, it is critical that it runs correctly. If errors are reported in other simulations or tests, be sure to run the MACAW tests and check for any error messages. If this runs without error, then all is well with these sub-routines and objects. here is no output file associated with this test and these algorithms are not used directly by the user, but are called by other algorithms in the library.\n\n");
	
	std::cout << "(7) MOLA (Molecule Object Library from Atoms)\n\n";
	puts("\tThis test runs a series of checks on our Molecule objects. These objects are built from the EEL atoms and allow for registration of new and/or existing molecules. If you ever get a message say that a molecule is not registerd in the library. First, check to make sure you used the proper alias/name for the molecule, then look through the mola.cpp file to see if that molecule has not yet been added to the record. This is a growing digital library so do not expect it to contain every molecule you want. There is no output file associated with this test and these algorithms are not used directly by the user, but are called by other algorithms in the library.\n\n");
	
	std::cout << "(8) MONKFISH (Multi-fiber wOven Nest Kernel For Interparticle Sorption History)\n\n";
	puts("\tThis test runs a simple example of the interparticle mass transfer uptake of aqueous ions into a woven conglomeration of cylindrical adsorbent fibers. Currently, there is no executable for this set of algorithms. NOTE: THIS IS A PLACE HOLDER! THIS TEST DOES NOTHING RIGHT NOW!\n\n");
	
	std::cout << "(9) SANDBOX - this one has no fancy acronym :-( \n\n";
	puts("\tThis test runs the sandbox executable. The sandbox is just an application where we store different temporary algorithms, methods, and functions before applying them in the rest of the library. Feel free to modify these source files (sandbox.h and sandbox.cpp) with any of your own methods and tests. Please note that after modifying any source file you must run make and make install in the primary directory of the ecosystem project folder before any changes that you make will be implemented\n\n");
	
	std::cout << "(10) SCOPSOWL (Simultaneously Coupled Objects for Pore and Surface diffusion Operations With Linear systems)\n\n";
	puts("\tThis test runs an example problem for gasoues, multi-species adsorption via a pore and surface diffusion mechanism through a bi-porous, spherical adsorbent pellet. The SCOPSOWL object is coupled to both the SKUA and MAGPIE objects, which are responsible for resolving surface diffusion and multi-species adsorption, respectively. There is an executable available for the user to interface with to run simulations and optimizations.\n\n");
	
	std::cout << "(11) SHARK (Speciation-object Hierarchy for Aqueous Reaction Kinetics)\n\n";
	puts("\tThis test runs an example problem for aqeous adsorption in a multi-species solution. Adsorption is represented as a metal-ligand complexation reaction. All other species in solution are at a pseudo-steady-state and are resolved in a series of speciation reactions coupled with overall mass balances on all major sub-species in solution. Currently, there is no executable for this set of algorithms.\n\n");
	
	std::cout << "(12) SKUA (Surface Kinetics for Uptake by Adsorption)\n\n";
	puts("\tThis test runs an example problem for gasoues, multi-species adsorption via a surface diffusion mechanism spherical adsorbent pellet. The SKUA object is coupled to the MAGPIE object, which resolves the multi-species adsorption equilibria between gas and solid phases. There is an executable available for the user to interface with to run simulations and optimizations.\n\n");
	
}

//Display info about the currently available executable routines
void current_execs()
{
	std::cout << "           CURRENTLY AVAILABLE EXECUTABLES           \n---------------------------------------------------------\n\n";
	
	std::cout << "(1) GSTA_OPT (Optimization Routine for GSTA analysis of adsorption equilibrium data)\n\n";
	puts("\tThis algorithm requires a single input file containing the raw adsorption data, as well as some other basis information about the system the data represents. It will produce a series of output files giving the full analysis of the data and the optimium equilibrium parameters associated with the Generalized Statistical Thermodynamic Adsorption (GSTA) isotherm.\n\n");
	
	std::cout << "(2) MAGPIE (Multicomponent Adsorption Generalized Procedure for Isothermal Equilibria)\n\n";
	puts("\tThis algorithm requires two input files: (i) a file containing the GSTA isotherm parameters for each adsorbing species in the system and (ii) a file detailing all the scenarios you wish to simulate. It will produce a single output file showing the results of the scenario simulations requested for adsorption capacity at various temperature, pressures, and gas compositions.\n\n");
	
	std::cout << "(3) SCOPSOWL (Simultaneously Coupled Objects for Pore and Surface diffusion Operations With Linear systems)\n\n";
	puts("\tThis algorithm requires four input files: (i) a scenario file detailing the system parameters and the species of interest, (ii) an adsorbent properties file giving information on the type of adsorbent used, (iii) a component properties file that gives basic information on specifc properies of each gaseous species, and (iv) an adsorbate properties file detailing the type of surface diffusion equation to use, the parameters of surface diffusion, and the isotherm parameters for the GSTA isotherm.\n\n");
	
	std::cout << "(4) SCOPSOWL_OPT (Optimization scheme for analysis of kinetic uptake data with the SCOPSOWL model)\n\n";
	puts("\tThis algorithm requires five input files: (i) a scenario file detailing some system parameters and the species of interest, (ii) an adsorbent properties file giving information on the type of adsorbent used, (iii) a component properties file that gives basic information on specifc properies of each gaseous species, (iv) an adsorbate properties file detailing the type of surface diffusion equation to use, the parameters of surface diffusion, and the isotherm parameters for the GSTA isotherm, and (v) a data file containing the actual adsorption time-series data for the model to be compared against. PLEASE NOTE, that the structure of these input files vary compared to running a standard SCOPSOWL simulation.\n\n");
	
	std::cout << "(5) SKUA (Surface Kinetics for Uptake by Adsorption)\n\n";
	puts("\tThis algorithm requires four input files: (i) a scenario file detailing the system parameters and the species of interest, (ii) an adsorbent properties file giving information on the type of adsorbent used, (iii) a component properties file that gives basic information on specifc properies of each gaseous species, and (iv) an adsorbate properties file detailing the type of surface diffusion equation to use, the parameters of surface diffusion, and the isotherm parameters for the GSTA isotherm.\n\n");
	
	std::cout << "(6) SKUA_OPT (Optimization scheme for analysis of kinetic uptake data with the SKUA model)\n\n";
	puts("\tThis algorithm requires five input files: (i) a scenario file detailing some system parameters and the species of interest, (ii) an adsorbent properties file giving information on the type of adsorbent used, (iii) a component properties file that gives basic information on specifc properies of each gaseous species, (iv) an adsorbate properties file detailing the type of surface diffusion equation to use, the parameters of surface diffusion, and the isotherm parameters for the GSTA isotherm, and (v) a data file containing the actual adsorption time-series data for the model to be compared against. PLEASE NOTE, that the structure of these input files vary compared to running a standard SKUA simulation.\n\n");
	
	puts("ADDITIONAL NOTES\n-----------------\n(i) Details on how each input file must be structured can be viewed the instructions file under doc sub-directory of the ecosystem project folder\n(ii) INSTRUCTIONS ARE CURRENTLY UNDER CONSTRUCTION AND ARE LIKELY INCOMPLETE!!!\n\n");
}

//Determine the users Test Option
bool valid_test_string(const std::string &input, UI_DATA *ui_dat)
{
	bool valid_input = false;
	
	if (allLower(input) == "dogfish")
	{
		ui_dat->option = dogfish;
		valid_input = true;
	}
	else if (allLower(input) == "eel")
	{
		ui_dat->option = eel;
		valid_input = true;
	}
	else if (allLower(input) == "finch")
	{
		ui_dat->option = finch;
		valid_input = true;
	}
	else if (allLower(input) == "lark")
	{
		ui_dat->option = lark;
		valid_input = true;
	}
	else if (allLower(input) == "macaw")
	{
		ui_dat->option = macaw;
		valid_input = true;
	}
	else if (allLower(input) == "mola")
	{
		ui_dat->option = mola;
		valid_input = true;
	}
	else if (allLower(input) == "monkfish")
	{
		ui_dat->option = monkfish;
		valid_input = true;
	}
	else if (allLower(input) == "sandbox")
	{
		ui_dat->option = sandbox;
		valid_input = true;
	}
	else if (allLower(input) == "scopsowl")
	{
		ui_dat->option = scopsowl;
		valid_input = true;
	}
	else if (allLower(input) == "shark")
	{
		ui_dat->option = shark;
		valid_input = true;
	}
	else if (allLower(input) == "skua")
	{
		ui_dat->option = skua;
		valid_input = true;
	}
	else if (allLower(input) == "trajectory")
	{
		ui_dat->option = trajectory;
		valid_input = true;
	}
	else
	{
		valid_input = false;
		if (ui_dat->BasicUI == true)
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
	}

	return valid_input;
}

//Check the string for a valid option
bool valid_exec_string(const std::string &input, UI_DATA *ui_dat)
{
	bool valid_input = false;
	
	if (allLower(input) == "gsta_opt")
	{
		ui_dat->option = gsta_opt;
		valid_input = true;
	}
	else if (allLower(input) == "magpie")
	{
		ui_dat->option = magpie;
		valid_input = true;
	}
	else if (allLower(input) == "scopsowl")
	{
		ui_dat->option = scopsowl;
		valid_input = true;
	}
	else if (allLower(input) == "scopsowl_opt")
	{
		ui_dat->option = scops_opt;
		valid_input = true;
	}
	else if (allLower(input) == "skua")
	{
		ui_dat->option = skua;
		valid_input = true;
	}
	else if (allLower(input) == "skua_opt")
	{
		ui_dat->option = skua_opt;
		valid_input = true;
	}
	else
	{
		valid_input = false;
		if (ui_dat->BasicUI == true)
			ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
	}
	
	return valid_input;
}

//Return the required number of files based on the executable option choosen
int number_files(UI_DATA *ui_dat)
{
	int num = 0;
	
	switch (ui_dat->option)
	{
		case scopsowl:
		{
			num = 4;
			break;
		}
			
		case skua:
		{
			num = 4;
			break;
		}
			
		case gsta_opt:
		{
			num = 1;
			break;
		}
			
		case magpie:
		{
			num = 2;
			break;
		}
			
		case scops_opt:
		{
			num = 5;
			break;
		}
			
		case skua_opt:
		{
			num = 5;
			break;
		}
			
		default:
		{
			mError(opt_no_support);
			return -1;
			break;
		}
	}
	
	return num;
}

//Check for valid additional options requested
bool valid_addon_options(UI_DATA *ui_dat)
{
	bool valid_input = false;
	int remander = (int) ui_dat->user_input.size();
	int files = 0;
	
	//Read the next option stored at offset 3
	ui_dat->Path = path(ui_dat->user_input[3]);
	ui_dat->Files = input(ui_dat->user_input[3]);
	remander = remander - 4;
	
	if (ui_dat->Path == false && ui_dat->Files == false)
		return false;
	
	//Read in files 
	if (ui_dat->Files == true)
	{
		files = number_files(ui_dat);
		if (files != remander)
		{
			ui_dat->MissingArg = true;
			return true;
		}
		valid_input = true;
		ui_dat->MissingArg = false;
		ui_dat->input_files.resize(files);
		for (int i=0; i<files; i++)
		{
			ui_dat->input_files[i] = ui_dat->user_input[i+4];
		}
	}
	//Read in path then input files
	else
	{
		ui_dat->path = ui_dat->user_input[4];
		remander--;
		
		ui_dat->Files = input(ui_dat->user_input[5]);
		remander--;
		
		if (ui_dat->Files == false)
			return false;
		
		files = number_files(ui_dat);
		if (files != remander)
		{
			ui_dat->MissingArg = true;
			return true;
		}
		valid_input = true;
		ui_dat->MissingArg = false;
		ui_dat->input_files.resize(files);
		for (int i=0; i<files; i++)
		{
			ui_dat->input_files[i] = ui_dat->path + ui_dat->user_input[i+6];
		}
	}
	
	return valid_input;
}

//Display a help message
void display_help(UI_DATA *ui_dat)
{
	if (ui_dat->argc == 1 || ui_dat->BasicUI == true)
	{
		std::cout << "The following message provides information on running this software with the basic ui...\n\n";
		std::cout << "           BASIC USER INTERFACE (BUI)          \n-----------------------------------------------\n";
		std::cout << "\tInitial BUI options (TEST and EXECUTABLES) allow the user to choose to run either a library test function\n";
		std::cout << "or a simulation based on the currently available, problem specific algorithms. The ouput from any test or\n";
		std::cout << "executable will be placed into a sub-directory named output. If no such directory exists, then one will be\n";
		std::cout << "created. You can then navigate to this directory to view output from the software.\n\n";
		
		current_tests();
		current_execs();
	}
	else
	{
		std::cout << "\nThe following message provides information on running this software with the advanced ui...\n\n";
		std::cout << "           ADVANCED USER INTERFACE (AUI)          \n----------------------------------------------------\n";
	}
	std::cout << std::endl;
}

//Display version info
void display_version(UI_DATA *ui_dat)
{
	if (ui_dat->argc == 1 || ui_dat->BasicUI == true)
	{
		std::cout << "This is the Basic UI for the " << ECO_VERSION << " version of the " << ECO_EXECUTABLE << " executable\n\n";
	}
	else
	{
		std::cout << "\nSpecifications\n---------------\nExecutable: " << ECO_EXECUTABLE << "\nVersion: " << ECO_VERSION << "\n\n";
	}
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
	std::cin >> ui_dat->user_input[0];
	std::cout << std::endl;
	ui_dat->value_type.editValue(ui_dat->user_input[0]);
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
		display_help(ui_dat);
		ui_dat->option = HELP;
		return true;
	}
	if (version(ui_dat->value_type.getValue()) == true)
	{
		display_version(ui_dat);
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
	std::cin >> ui_dat->user_input[0];
	std::cout << std::endl;
	ui_dat->value_type.editValue(ui_dat->user_input[0]);
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
		display_help(ui_dat);
		ui_dat->option = HELP;
		return true;
	}
	if (version(ui_dat->value_type.getValue()) == true)
	{
		display_version(ui_dat);
		ui_dat->option = HELP;
		return true;
	}
	
	if (ui_dat->value_type.getType() == INT)
	{
		
		switch (ui_dat->value_type.getInt())
		{
			case 1:
			{
				ui_dat->option = dogfish;
				valid_input = true;
				break;
			}
				
			case 2:
			{
				ui_dat->option = eel;
				valid_input = true;
				break;
			}
				
			case 3:
			{
				ui_dat->option = egret;
				valid_input = true;
				break;
			}
				
			case 4:
			{
				ui_dat->option = finch;
				valid_input = true;
				break;
			}
				
			case 5:
			{
				ui_dat->option = lark;
				valid_input = true;
				break;
			}
				
			case 6:
			{
				ui_dat->option = macaw;
				valid_input = true;
				break;
			}
				
			case 7:
			{
				ui_dat->option = mola;
				valid_input = true;
				break;
			}
				
			case 8:
			{
				ui_dat->option = monkfish;
				valid_input = true;
				break;
			}
				
			case 9:
			{
				ui_dat->option = sandbox;
				valid_input = true;
				break;
			}
				
			case 10:
			{
				ui_dat->option = scopsowl;
				valid_input = true;
				break;
			}
				
			case 11:
			{
				ui_dat->option = shark;
				valid_input = true;
				break;
			}
				
			case 12:
			{
				ui_dat->option = skua;
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
		valid_input = valid_test_string(ui_dat->value_type.getString(), ui_dat);
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
	std::cin >> ui_dat->user_input[0];
	std::cout << std::endl;
	ui_dat->value_type.editValue(ui_dat->user_input[0]);
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
		display_help(ui_dat);
		ui_dat->option = HELP;
		return true;
	}
	if (version(ui_dat->value_type.getValue()) == true)
	{
		display_version(ui_dat);
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
		valid_input = valid_exec_string(ui_dat->value_type.getString(), ui_dat);
	}
	else
	{
		valid_input = false;
		ui_dat->option = invalid_input(ui_dat->count, ui_dat->max);
	}
	
	return valid_input;
}

//Loop through basis ui test menu until success
int test_loop(UI_DATA *ui_dat)
{
	int success = 0;
	bool valid_input = false;
	do
	{
		valid_input = valid_input_tests(ui_dat);
		if (ui_dat->option == EXIT)
			return 0;
		if (ui_dat->option == HELP)
		{
			ui_dat->count--;
			valid_input = false;
		}
		ui_dat->count++;
	} while (valid_input == false && ui_dat->count < ui_dat->max);
	
	success = run_test(ui_dat);
	
	return success;
}

//Run through executable loop till valid choice is made
int exec_loop(UI_DATA *ui_dat)
{
	int success = 0;
	bool valid_input = false;
	do
	{
		valid_input = valid_input_execute(ui_dat);
		if (ui_dat->option == EXIT)
			return 0;
		if (ui_dat->option == HELP)
		{
			ui_dat->count--;
			valid_input = false;
		}
		ui_dat->count++;
	} while (valid_input == false && ui_dat->count < ui_dat->max);
	
	success = run_exec(ui_dat);
	
	return success;
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
			
		case trajectory:
			success = Run_Trajectory();
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
	
	if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
	{
		std::cout << "WARNING! All simulation runs will require properly formatted input files.\n";
		std::cout << "User is required to provide full path and extension to each necessary file.\n";
		std::cout << "(Example: this/is/path/to/input/files.txt)\n\n";
	}
	
	switch (ui_dat->option)
	{
		case scopsowl:
		{
			if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
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
			}
			
			success = SCOPSOWL_SCENARIOS(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str());
			break;
		}
			
		case skua:
		{
			if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
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
			}
			
			success = SKUA_SCENARIOS(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str());
			break;
		}
			
		case gsta_opt:
		{
			if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
			{
				ui_dat->input_files.resize(1);
				std::cout << "GSTA_OPT needs 1 input file containing the equilibrium adsorption data.\n";
				std::cout << "Please provide the file and full path...\n\n";
				std::cout << "GSTA input: ";
				std::cin >> ui_dat->input_files[0];
				std::cout << std::endl;
			}
			
			success = gsta_optimize(ui_dat->input_files[0].c_str());
			break;
		}
			
		case magpie:
		{
			if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
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
			}
			
			success = MAGPIE_SCENARIOS(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str());
			break;
		}
			
		case scops_opt:
		{
			if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
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
			}
			
			success = SCOPSOWL_OPTIMIZE(ui_dat->input_files[0].c_str(), ui_dat->input_files[1].c_str(), ui_dat->input_files[2].c_str(), ui_dat->input_files[3].c_str(), ui_dat->input_files[4].c_str());
			break;
		}
			
		case skua_opt:
		{
			if (ui_dat->argc == 1 || ui_dat->MissingArg == true)
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
			}
			
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
		ui_dat.user_input.resize(1);
		ui_dat.MissingArg = false;
		ui_dat.BasicUI = true;
		ui_dat.Path = false;
		ui_dat.Files = false;
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
			success = test_loop(&ui_dat);
		}
		else if (ui_dat.option == EXECUTE)
		{
			success = exec_loop(&ui_dat);
		}
		else
		{
			mError(opt_no_support);
			return -1;
		}
	}//END Text based menu for Basic UI
	
	//Run the Advanced UI
	else
	{
		ui_dat.BasicUI = false;
		ui_dat.MissingArg = true;
		ui_dat.Path = false;
		ui_dat.Files = false;
		bool valid_input = false;
		
		//Store all console commands given
		for (int i=0; i<argc; i++)
		{
			ui_dat.user_input.push_back(argv[i]);
		}
		//NOTE: The first command will always be the excutable, so we can ignore or override it later
		
		//Check to see if the user requests help
		if (help(ui_dat.user_input[1]) == true)
		{
			display_help(&ui_dat);
			return 0;
		}
		
		//Check to see if user asks for version info
		if (version(ui_dat.user_input[1]) == true)
		{
			display_version(&ui_dat);
			return 0;
		}
		
		//Check to see if user requests running tests (-t) or executables (-e)
		if (test(ui_dat.user_input[1]) == true)
		{
			if (argc > 2)
			{
				ui_dat.MissingArg = false;
				valid_input = valid_test_string(ui_dat.user_input[2], &ui_dat);
				if (valid_input == true)
				{
					success = run_test(&ui_dat);
				}
				else
				{
					mError(invalid_console_input);
					display_help(&ui_dat);
					return -1;
				}
			}
			else
			{
				ui_dat.BasicUI = true;
				success = test_loop(&ui_dat);
			}
		}
		else if (exec(ui_dat.user_input[1]) == true)
		{
			if (argc > 2)
			{
				valid_input = valid_exec_string(ui_dat.user_input[2], &ui_dat);
				if (valid_input == true)
				{
					if (argc == 3)
					{
						ui_dat.MissingArg = true;
						success = run_exec(&ui_dat);
					}
					else if (argc == 4)
					{
						mError(invalid_console_input);
						display_help(&ui_dat);
						return -1;
					}
					else
					{
						valid_input = valid_addon_options(&ui_dat);
						if (valid_input == true)
						{
							success = run_exec(&ui_dat);
						}
						else
						{
							mError(invalid_console_input);
							display_help(&ui_dat);
							return -1;
						}
					}
				}
				else
				{
					mError(invalid_console_input);
					display_help(&ui_dat);
					return -1;
				}
			}
			else
			{
				ui_dat.BasicUI = true;
				success = exec_loop(&ui_dat);
			}
		}
		else
		{
			mError(invalid_console_input);
			display_help(&ui_dat);
			return -1;
		}
		
	}//END Advanced UI
	
	return success;
}