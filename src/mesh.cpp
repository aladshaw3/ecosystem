/*!
 *  \file mesh.h
 *	\brief Mesh Objects and Associated Sub-Objects
 *	\details This kernel allows for creation of mesh objects from constitutient sub-objects
 *			such as Nodes and Elements (Line, Surface, or Volume). Mesh objects can be used
 *			to establish physical-chemcial simulations in multi-dimensional space.
 *
 *			\warning This kernel is still under active development. Use with caution!
 *
 *
 *  \author Austin Ladshaw
 *	\date 04/09/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2018, all
 *             rights reserved.
 */

#include "mesh.h"

// Test function for MESH kernel
int MESH_TESTS()
{
	// --- Initializations ----
	int success = 0;
	double time = clock();
	std::cout << "\nStart Running Tests on Mesh Objects...\n";
	
	//---Exit Messages and cleanup---
	time = clock() - time;
	std::cout << "\nMESH Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	return success;
}
