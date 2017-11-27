/*!
 *  \file crow.h
 *	\brief Coupled Reaction Object Workspace
 *	\details This file creates objects and subroutines for setting up and solving systems
 *			of reaction driven equations using the DOVE (see dove.h) solver. It combines
 *			a generalized description of chemical reaction mathematics with a comprehensive
 *			input file framework to allow systems of equations to be developed on the fly
 *			and solved with reasonable accuracy and efficiency. 
 *
 *			Mathematical description of generic reaction:
 *			---------------------------------------------
 *			a_i*du_i/dt = k1*Product(j,u_j^v_j) - k2*Product(l,u_l^v_l)
 *
 *			where i,j,l are indices of variables, k1,k2 are reaction constants, a is a time
 *			coefficient, and Product(i,arg) is the product of all args in i.
 *
 *			\warning This kernel is still under active development. Use with caution!
 *
 *
 *  \author Austin Ladshaw
 *	\date 11/27/2017
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2017, all
 *             rights reserved.
 */

#include "dove.h"

#ifndef crow_h
#define crow_h

//Run CROW scenario
int CROW_SCENARIO(const char *yaml_input);

//Run the CROW test
int CROW_TESTS();

#endif /* crow_h */
