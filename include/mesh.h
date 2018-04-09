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

#ifndef MESH_HPP_
#define MESH_HPP_

#include <stdio.h>				//Line to allow cout functionality
#include <math.h>               //Line added to allow usage of the pow (e, x) function
#include <iostream>				//Line to allow for read/write to the console using cpp functions
#include <fstream>				//Line to allow for read/write to and from .txt files
#include <stdlib.h>				//Line need to convert strings to doubles
#include <vector>				//Line needed to use dynamic arrays called vectors
#include <time.h>				//Line needed to display program runtime
#include <float.h>				//Line to allow use of machine precision constants
#include <string>				//Line to allow use of strings as a data type
#include <exception>            //Line to allow use of try-catch statements
#include "error.h"

/// Enumeration for the list of valid node types
/** The only types that have been defined are for Boudary and Interior nodes.*/
typedef enum {BOUNDARY, INTERIOR} node_type;

/// Node object
/** This class structure creates a C++ object for a node in a mesh. The node will have
	an identifying ID number and a sub_type (either BOUNDARY or INTERIOR) to add in identification
	of different aspects the node has in the overall mesh. All nodes are considered to be points
	in 3D space, whether or not the final mesh is 3D. As such, every node will have a coordinate
	vector (x, y, z) associated with it. */
class Node
{
public:
	Node();												///< Default constructor
	~Node();											///< Default destructor
	
protected:
	
private:
	
	
};

/// Test function for MESH kernel
/** This function runs tests on the mesh objects. */
int MESH_TESTS();

#endif /* MESH_HPP_ */
