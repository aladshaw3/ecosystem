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
#include "yaml_wrapper.h"


/// Enumeration for the list of valid decay modes
/** List of valid types of radioactive decay. The type of decay dictates
	how the Isotope object will transform the given isotope. */
typedef enum {alpha, beta_min, beta_plus, stable, spon_fiss, iso_trans, neutron_em, proton_em,
				beta_min_neutron_em, beta_plus_neutron_em, beta_plus_alpha, beta_plus_beta_plus,
				beta_min_beta_min, beta_min_2neutron_em, beta_min_alpha, proton_em_proton_em,
				neutron_em_neutron_em, beta_min_3neutron_em, beta_min_4neutron_em, beta_plus_2proton_em,
				beta_plus_3proton_em, specific_isotope, beta_plus_proton_em, undefined} decay_mode;

/// Enumeration for the list of valid units of half-life
/** List of valid units for half-lifes for better readability of code.*/
typedef enum {seconds, minutes, hours, days, years} time_units;

/// Function to convert from a starting unit and value to and ending unit and value (returns converted value)
double time_conversion(time_units end_unit, double start_value, time_units start_unit);

/// Function to determine what type of time units are used based on given string argument
time_units timeunits_choice(std::string &choice);

/// Function to determine the decay mode based on given string
decay_mode decaymode_choice(std::string &choice);

/// Function to return a string for the name of a decay mode given the enum
std::string decaymode_string(decay_mode mode);

/// Function to return a string for the units
std::string timeunits_string(time_units units);

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
	
	void loadNuclides(yaml_cpp_class &data);			///< Function to load the nuclide library into the pointer
	void unloadNuclides();								///< Delete the pointer to nuclide library to free space
	void registerIsotope(std::string isotope_name);		///< Register an isotope given the isotope name
	void registerIsotope(std::string symbol, int iso);	///< Register an isotope given an atomic symbol (e.g., H) and isotope number (e.g., 2)
	void registerIsotope(int atom_num, int iso_num);	///< Register an isotope given both an atomic and isotope number
	
	void DisplayInfo();									///< Print out isotope information to the console
	
	int IsotopeNumber();								///< Return the isotope number of the atom
	double DecayRate();									///< Return the decay rate of the isotope
	double HalfLife(time_units units);					///< Return the half-life in the given units
	time_units HalfLifeUnits();							///< Return the half-life units
	std::string IsotopeName();							///< Return the name of the isotope
	bool isStable();									///< Return stability condition
	bool isIsomericState();								///< Return isomeric condition
	int DecayModes();									///< Return the number of decay modes
	
	decay_mode DecayMode(int i);						///< Return the ith decay mode
	double BranchFraction(int i);						///< Return the ith branch fraction
	std::string ParticleEmitted(int i);					///< Return the name of the particle emitted for the ith decay mode
	int NumberParticlesEmitted(int i);					///< Return the number of particles that get emitted
	std::string Daughter(int i);						///< Return the name of the daughter isotope 
	
protected:
	std::string IsoName;								///< Name of the isotope (e.g., H-2)
	std::vector<decay_mode> decay_modes;				///< List of decay modes the given isotope can undergo
	std::vector<double> branch_ratios;					///< Branching ratios for each possible decay mode
	std::vector<std::string> particle_emitted;			///< Name of the particle(s) ejected during spec_iso decay
	std::vector<int> num_particles;						///< Numbers of particles emitted during decay mode
	std::vector<std::string> daughter;					///< Name of the daughter isotope formed
	
	double decay_rate;									///< Rate of decay for the given isotope (1/s)
	double half_life;									///< Half-life of the isotope (in hl_units)
	time_units hl_units;								///< Units given for the half-life
	int isotope_number;									///< isotope number for the object
	bool Stable;										///< Boolean is True if isotope is stable
	bool IsomericState;									///< Boolean is True if isotope is in an isomeric state
	
	yaml_cpp_class *nuclides;							///< Pointer to a yaml object storing the digital library of all nuclides
	
	void setConstants();								///< Set the decay_modes, branch_ratios, and other info based on load library
	void computeDecayRate();							///< Compute the decay rate (in 1/s) based on the half-life
	YamlWrapper& getNuclideLibrary();					///< Return reference to the nuclide library
	
private:
	
};

///< Function to test the implementation of Isotope 
int IBIS_TESTS();

#endif /* IBIS_HPP_ */
