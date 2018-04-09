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




#endif /* MESH_HPP_ */
