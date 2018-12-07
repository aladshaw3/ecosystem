/*!
 *  \file fairy.h fairy.cpp
 *	\brief Fission-products from Atomic Incident and their Respective Yields
 *	\details This file contains a C++ object for determining fission products and
 *			their yields from some nuclear event based on: (i) type of fission,
 *			either neutron-induced or spontaneous, (ii) energy level of neutron
 *			source or bomb yield, (iii) extent of fission, and (iv) initial mass
 *			and composition of fuel or bomb materials. Data for fission products
 *			comes from ENDF-6 data libraries that were read with a python script
 *			and output into a yaml format (see 'scripts/fission-product-yields').
 *
 *  \author Austin Ladshaw
 *	\date 12/07/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay and transport. Copyright (c) 2018, 
 *				all rights reserved.
 */

#ifndef FAIRY_HPP_
#define FAIRY_HPP_

#include "ibis.h"


/// Test function for FAIRY
int FAIRY_TESTS();

#endif /* FAIRY_HPP_ */
