/*!
 *  \file dove.h
 *	\brief Dynamic Ode solver with Various Established methods
 *	\details This file creates objects and subroutines for solving systems of Ordinary Differential
 *			Equations using various established methods. The basic idea is that a user will create
 *			a function to calculate all the right-hand sides of a system of ODEs, then pass that
 *			function to the DOVE routine, which will seek a numerical solution to that system.
 *
 *			Methods for Integration
 *			-----------------------
 *			(None available - still under construction)
 *
 *	\note This kernel is still under construction.
 *
 *  \author Austin Ladshaw
 *	\date 09/25/2017
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2017, all
 *             rights reserved.
 */

#include "macaw.h"
#include "lark.h"
#include "yaml_wrapper.h"

#ifndef DOVE_HPP_
#define DOVE_HPP_

/// Test function for DOVE kernel
/** This function sets up and solves a test problem for DOVE. It is callable from the UI. */
int DOVE_TESTS();

#endif /* DOVE_HPP_ */
