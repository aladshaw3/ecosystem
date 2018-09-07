/*!
 *  \file ibis.h ibis.cpp
 *	\brief Implicit Branching Isotope System
 *	\details This file contains a C++ object for creating and utilizing isotopes 
 *			in a branching isotope decay chain. The object inherits from Class Atom
 *			(see eel.h), which creates individual atoms from names, symbols, or 
 *			atomic numbers, then adds the ability to decay those atoms to different
 *			isotopes through various decay modes. Added to the Atom object are parameters
 *			for decay rates (half-lifes), branching ratios, and decay modes (alpha, beta,
 *			etc). The intent of this system is to determine how radioactive decay occurs
 *			when given a sinle or multiple starting isotopes.
 *
 *  \author Austin Ladshaw
 *	\date 09/04/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay. Copyright (c) 2018, all
 *             rights reserved.
 */

#ifndef IBIS_HPP_
#define IBIS_HPP_

#include "eel.h"


/// Enumeration for the list of valid decay modes
/** List of valid types of radioactive decay. The type of decay dictates
	how the Isotope object will transform the given isotope. */
typedef enum {alpha, beta_plus, beta_minus, stable, elec_cap, neut_emiss, prot_emiss} decay_mode;

/// Enumeration for the list of valid units of half-life
/** List of valid units for half-lifes for better readability of code.*/
typedef enum {seconds, minutes, hours, days, years} time_units;

/// Function to convert from a starting unit and value to and ending unit and value (returns converted value)
double time_conversion(time_units end_unit, double start_value, time_units start_unit);

/// Isotope object to hold information and provide decay operations
/** This is the C++ object to store information and functions associated with the
	decay of radioactive isotopes. It herits from the Atom object and extends that
	object to include information such as decay constants, branching ratios, and
	decay modes. It can be used to determine branching paths and setup systems of
	equations involving decay products. */
class Isotope : Atom
{
public:
	Isotope();											///< Default constructor
	~Isotope();											///< Default destructor
	
	void registerIsotope(std::string isotope_name);		///< Register an isotope given the isotope name
	void registerIsotope(std::string symbol, int iso);	///< Register an isotope given an atomic symbol (e.g., H) and isotope number (e.g., 2)
	void registerIsotope(int atom_num, int iso_num);	///< Register an isotope given both an atomic and isotope number
	
	int IsotopeNumber();								///< Return the isotope number of the atom
	double DecayRate();									///< Return the decay rate of the isotope
	double HalfLife(time_units units);					///< Return the half-life in the given units
	std::string IsotopeName();							///< Return the name of the isotope
	
protected:
	std::string IsoName;								///< Name of the isotope (e.g., H-2)
	std::vector<decay_mode> decay_modes;				///< List of decay modes the given isotope can undergo
	std::vector<double> branch_ratios;					///< Branching ratios for each possible decay mode
	
	double decay_rate;									///< Rate of decay for the given isotope (1/s)
	double half_life;									///< Half-life of the isotope (in hl_units)
	time_units hl_units;								///< Units given for the half-life
	int isotope_number;									///< isotope number for the object
	
	void setConstants();								///< Set the decay_modes, branch_ratios, and other info based on atomic and isotope numbers
	void computeDecayRate();							///< Compute the decay rate (in 1/s) based on the half-life
	
private:
	
};

///< Function to test the implementation of Isotope 
int IBIS_TESTS();

#endif /* IBIS_HPP_ */
