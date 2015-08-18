//----------------------------------------
//  Created by Austin Ladshaw on 1/2/14
//  Copyright (c) 2014
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 * error.C
 *
 *      v1.0.0
 *
 *      0.0.1 - This class file will be expanded upon as need be for more specific error reporting by all programs
 *      		which are under development. Currently, serves to give only a generic error message.
 *
 *      1.0.0 - Error messages are flagged by an integer. All current messages are those previously defined within
 *      		each class file.
 */

#include "error.h"

//Error message to be displayed when exceptions are caught or problems occur
void error(int flag)
{
	if (flag == 1)
		std::cout << "\nError!\n\nInput File does not exist!" << std::endl;
	else if (flag == 2)
		std::cout << "\nError!\n\nIndexing Mistake! Check input file for errors..." << std::endl;
	else if (flag == 3)
		std::cout << "\nError!\n\nNot used in reverse evaluations!" << std::endl;
	else if (flag == 4)
		std::cout << "\nError!\n\nSimulation Failed!" << std::endl;
	else if (flag == 5)
		std::cout << "\nError!\n\nInvalid number of components!" << std::endl;
	else if (flag == 6)
		std::cout << "\nError!\n\nInvalid Boolean Selection!" << std::endl;
	else if (flag == 7)
		std::cout << "\nError!\n\nInvalid Mole Fraction!" << std::endl;
	else if (flag == 8)
		std::cout << "\nError!\n\nInvalid Sum of Gas Mole Fractions!" << std::endl;
	else if (flag == 9)
		std::cout << "\nError!\n\nInvalid Sum of Solid Mole Fractions!" << std::endl;
	else if (flag == 10)
		std::cout << "\nError!\n\nScenario simulation failed!" << std::endl;
	else if (flag == 11)
		std::cout << "\nError!\n\nIndex Out of Bounds of Matrix!" << std::endl;
	else if (flag == 12)
		std::cout << "\nError!\n\nNon-Square Matrix!" << std::endl;
	else if (flag == 13)
		std::cout << "\nError!\n\nDimesional Mis-Match for Matrices!" << std::endl;
	else if (flag == 14)
		std::cout << "\nError!\n\nEmpty or Non-Initialized Matrix!" << std::endl;
	else if (flag == 15)
		std::cout << "\nError!\n\nOption not currently supported..." << std::endl;
	else if (flag == 16)
		std::cout << "\nError!\n\nFractions should be between 0 and 1!" << std::endl;
	else if (flag == 17)
		std::cout << "\nWarning!\n\nOrthogonallity Check Failed!" << std::endl;
	else if (flag == 18)
		std::cout << "\nWarning!\n\nInstability in Coefficient Matrix!" << std::endl;
	else if (flag == 19)
		std::cout << "\nWarning!\n\nNothing to Simulate! No Diffusion Can Occur!" << std::endl;
	else if (flag == 20)
		std::cout << "\nWarning!\n\nNegative Mass Encountered!" << std::endl;
	else if (flag == 21)
		std::cout << "\nError!\n\nNegative Time Encountered!" << std::endl;
	else if (flag == 22)
		std::cout << "\nError!\n\nMatrix and Vector Sizes Do Not Match!" << std::endl;
  	else if (flag == 23)
    	std::cout << "\nError!\n\nMatrix argument and storage matrix are the same!" << std::endl;
    else if (flag == 24)
        std::cout << "\nWarning!\n\nMatrix is Singular or Close to Singular!" << std::endl;
    else if (flag == 25)
        std::cout << "\nError!\n\nMatrix is too small for this function!" << std::endl;
    else if (flag == 26)
        std::cout << "\nError!\n\nCan not build size 0 or smaller matrix!" << std::endl;
  	else if (flag == 27)
    	std::cout << "\nError!\n\nGiven NULL pointer as function!" << std::endl;
	else if (flag == 28)
    	std::cout << "\nError!\n\nNorms must be non-negative scalars!" << std::endl;
	else if (flag == 29)
		std::cout << "\nError!\n\nOut of Bounds of Real Vector!" << std::endl;
	else if (flag == 30)
		std::cout << "\nError!\n\nZero vector cannot be evaluated on!" << std::endl;
	else if (flag == 31)
		std::cout << "\nError!\n\nOut of Bounds of Real Tensor!" << std::endl;
	else if (flag == 32)
		std::cout << "\nError!\n\nCannot construct edge from same two points!" << std::endl;
	else if (flag == 33)
		std::cout << "\nError!\n\nNULL POINTER ERROR!!!" << std::endl;
	else if (flag == 34)
		std::cout << "\nError!\n\nInvalid or Unsupported Atom!" << std::endl;
	else if (flag == 35)
		std::cout << "\nError!\n\nInvalid number of protons!" << std::endl;
	else if (flag == 36)
		std::cout << "\nError!\n\nInvalid number of neutrons!" << std::endl;
	else if (flag == 37)
		std::cout << "\nError!\n\nInvalid number of electrons!" << std::endl;
	else if (flag == 38)
		std::cout << "\nError!\n\nInvalid number of valence electrons!" << std::endl;
	else if (flag == 39)
		std::cout << "\nError!\n\nUnexpected value during string parsing!" << std::endl;
	else if (flag == 40)
		std::cout << "\nError!\n\nMolecule formula not registered or not found!" << std::endl;
	else if (flag == 41)
		std::cout << "\nError!\n\nReaction rates ratios do not match equilibrium constant!" << std::endl;
	else if (flag == 42)
		std::cout << "\nError!\n\nSpecies not found in list! Check names for errors." << std::endl;
	else if (flag == 43)
		std::cout << "\nError!\n\nDuplicate variable found in list! Check names for errors." << std::endl;
	else if (flag == 44)
		std::cout << "\nError!\n\nMissing Information! Cannot complete requested tasks." << std::endl;
	else if (flag == 45)
		std::cout << "\nError!\n\nInvalid Type Specification!" << std::endl;
	else if (flag == 46)
		std::cout << "\nError!\n\nKey does not exist in the hash table!" << std::endl;
	else if (flag == 47)
		std::cout << "\nError!\n\nAnchor Alias Pair does not exist!" << std::endl;
	else
		std::cout << "\nUndefined Error!!!" << std::endl;
}
