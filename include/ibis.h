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
#include "dove.h"
#include <unordered_map>


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
	void clearChain();									///< Delete the chain for this isotope to free space
	void registerIsotope(std::string isotope_name);		///< Register an isotope given the isotope name (e.g., H-2)
	void registerIsotope(std::string symbol, int iso);	///< Register an isotope given an atomic symbol (e.g., H) and isotope number (e.g., 2)
	void registerIsotope(int atom_num, int iso_num);	///< Register an isotope given both an atomic and isotope number (e.g., H-2 = 1, 2)
	
	void DisplayInfo();									///< Print out isotope information to the console
	void DisplayChain();								///< Print out chain information to the console
	void createChain();									///< Function to create and fill in the chain of nuclides for this starting isotope
	
	void setInitialCondition(double ic);				///< Set the value for the initial condition of this nuclide
	void setConcentration(double c);					///< Set the concentration value for the nuclide
	
	int IsotopeNumber();								///< Return the isotope number of the atom
	double DecayRate();									///< Return the decay rate of the isotope
	double HalfLife(time_units units);					///< Return the half-life in the given units
	time_units HalfLifeUnits();							///< Return the half-life units
	std::string IsotopeName();							///< Return the name of the isotope
	bool isStable();									///< Return stability condition
	bool isIsomericState();								///< Return isomeric condition
	int DecayModes();									///< Return the number of decay modes
	double getInitialCondition();						///< Return the value of the initial condition
	double getConcentration();							///< Return the concentration value of the nuclide 
	
	decay_mode DecayMode(int i);						///< Return the ith decay mode
	double BranchFraction(int i);						///< Return the ith branch fraction
	std::string ParticleEmitted(int i);					///< Return the name of the particle emitted for the ith decay mode
	int NumberParticlesEmitted(int i);					///< Return the number of particles that get emitted
	std::string Daughter(int i);						///< Return the name of the daughter isotope
	
	/// Return a list of indices of the decay modes that this daughter isotope is formed from given the parent isotope's name
	/** This function will iterate through the decay modes for the parent isotope we are investigating and return a
		list of indices that represent the modes of decay that form this isotope from the given parent. If the mode
		does not form this isotope, then it will have an index of -1. The size of the vector on return will equal
		the size of the list of decay modes for the parent. If the parent is invalid, then it will return a vector
		of size 0. If the parent does not form this isotope, then all indices will be -1. */
	std::vector<int> DaughterIndices(std::string parent);
	
	/// Return a list of indices of the decay modes that this particle emission is formed from given the parent isotope's name
	/** This function will iterate through the decay modes for the parent isotope we are investigating and return a
		list of indices that represent the modes of decay that form this isotope from the given parent. If the mode
		does not form this isotope, then it will have an index of -1. The size of the vector on return will equal
		the size of the list of decay modes for the parent. If the parent is invalid, then it will return a vector
		of size 0. If the parent does not form this isotope, then all indices will be -1. */
	std::vector<int> EmissionIndices(std::string parent);
	
protected:
	std::string IsoName;								///< Name of the isotope (e.g., H-2)
	std::vector<decay_mode> decay_modes;				///< List of decay modes the given isotope can undergo
	std::vector<double> branch_ratios;					///< Branching ratios for each possible decay mode
	std::vector<std::string> particle_emitted;			///< Name of the particle(s) ejected during spec_iso decay
	std::vector<int> num_particles;						///< Numbers of particles emitted during decay mode
	std::vector<std::string> daughter;					///< Name of the daughter isotope formed
	
	/// List of all parent and daughter/particles pairs in a chain given this isotope
	/** Stores the list of parent and daughter pairs. If the parent produces emitted particles, those go
		in this list as well. If the parent produces more than one daughter, all those daughters go in the 
		list as well. 
	 
		chain[i] = contains lists of parent daughter pairs at the ith level of the chain
		chain[i][j] = contains the jth pair at the ith level
		chain[i][j].first = contains the name of the parent
		chain[i][j].second = contains the name of the daughter/particle emitted
		
		Example:
		--------
	 
		Levels:							A					A = Principal Parent
		-------						   / \
		level 0:					  B   C					B = 1st Daughter (pair A/B), C = 1st particle emission (pair A/C)
									 / \   \
		level 1:                    D   E   F               D = 1st Daugher of 1st Daughter, E = 2nd Daughter of 1st Daughter, F = Stable Daughter
								  / |   | \
		level 2:                 G  H   I  J				So on and so forth...
	 
	 */
	std::vector< std::vector< std::pair<std::string,std::string> > > chain;
	
	double decay_rate;									///< Rate of decay for the given isotope (1/s)
	double half_life;									///< Half-life of the isotope (in hl_units)
	time_units hl_units;								///< Units given for the half-life
	int isotope_number;									///< isotope number for the object
	bool Stable;										///< Boolean is True if isotope is stable
	bool IsomericState;									///< Boolean is True if isotope is in an isomeric state
	double initial_condition;							///< Value to hold initial condition for this nuclide
	double concentration;								///< Value to hold concentration after a point in time
	
	yaml_cpp_class *nuclides;							///< Pointer to a yaml object storing the digital library of all nuclides
	
	void setConstants();								///< Set the decay_modes, branch_ratios, and other info based on load library
	void computeDecayRate();							///< Compute the decay rate (in 1/s) based on the half-life
	int addPairs(int i, std::string parent);			///< Function to add parent/daughter pairs given the parent's name and the current level
	YamlWrapper& getNuclideLibrary();					///< Return reference to the nuclide library
	
private:
	
};

/// DecayChain object to hold and store a set of unique isotopes in a branched decay chain
/** C++ style object that will contain a list of unique nuclides that can be used to numerically
	solve a decay chain system. User will provide a list of initial nuclides present and this
	object will then use that information to build the list of all possible nuclides that can
	be formed. One of the key features of this object is that the nuclide list will be unique
	(i.e., no duplicate nuclides) so that we can iterate through that list to apply the branch
	fractions and decay constants to solve the fractionation as a function of time.
 
	First, create chains of isotope from the initial nuclides given. Then, reduce chains to remove
	redundant isotopes and order the isotopes in final_nuc list to ensure that only the highest
	mass number isotopes are listed first. Lastly, fill out the parents and branches lists to
	allow for direct access to all necessary decay information to speed up computation and
	evaluation of the chains. */
class DecayChain
{
public:
	DecayChain();											///< Default constructor
	~DecayChain();											///< Default destructor
	
	void DisplayList();										///< Display list of nuclides to console
	void DisplayInfo();										///< Display final nuclides and their parents and branches
	void DisplayMap();										///< Display the coefficient map to the console
	
	void loadNuclides(yaml_cpp_class &data);				///< Function to load the nuclide library into the pointer
	void unloadNuclides();									///< Delete the pointer to nuclide library to free space
	
	void registerInitialNuclide(std::string isotope_name);	///< Register an initial nuclide by name (e.g., H-2)
	void registerInitialNuclide(std::string symb, int iso);	///< Register an initial nuclide by symbol (e.g., H) and isotope number (e.g., 2)
	void registerInitialNuclide(int atom_num, int iso_num);	///< Register an initial nuclide by atomic and mass numbers (e.g., H-2 = 1, 2)
	
	void createChains();									///< Function to create unique list of final nuclides from decay chains of initial
	void formEigenvectors();								///< Function to produce eigenvectors from coefficient matrix
	void verifyEigenSoln();									///< Function will verify that the eigenvectors and eigenvalues are correct 
	
	int getNumberNuclides();								///< Return the number of nuclides in the decay chain
	std::vector<int>& getParentList(int i);					///< Return the vector list of parents for the ith isotope in the nuclide list
	std::vector<int>& getBranchList(int i, int j);			///< Return the vector list of branch fractions for the jth parent of the ith nuclide
	Isotope& getIsotope(int i);								///< Return the ith isotope in the nuclide list
	Matrix<double>& getEigenvectors();						///< Return the eigenvectors matrix
	Matrix<double>& getInverseEigenvectors();				///< Return the eigenvectors inverse matrix
	
protected:
	
	void roughInsertSort(Isotope iso);						///< Insert an isotope to the initial nuclide list and sort according to isotope number
	void finalSort();										///< Sort the list of nuclides after creating the chains
	
	/// Return a sorted list of isotopes given reference to another list
	/** All isotopes in the list given will have the same isotope number, so we are 
		sorting the list based on partents and daughters. */
	std::vector<Isotope> sameIsoNumSort(std::vector<Isotope> &list);
	
	void fillOutBranchData();								///< Function to fill out all branching data and parent indices
	
	void fillOutCoefMap();									///< Function to fill out all coefficients in the map
	
private:
	std::vector<Isotope> nuc_list;							///< List of (ith) nuclides that make up the decay chain
	
	/// List of the indices that each isotope has from the final_nuc list
	/** List of indices for all parents of an isotope in the list of all isotopes. 
		
		parents[i] = list of parent indices for the ith isotope
		parents[i][j] = jth parent index of the ith isotope */
	std::vector< std::vector<int> > parents;
	
	/// List of the indices for all branches of a parent to an isotope in final_nuc list
	/** List of indices for direct access to necessary information to develop the chain scheme.
	 
		branches[i] = list of a list of all branch indices (by parent) for ith isotope
		branches[i][j] = list of all branches that contribute to formation of ith isotope by jth parent
		branches[i][j][k] = kth branch index for the jth parent that forms the ith isotope */
	std::vector< std::vector< std::vector<int> > > branches;
	
	std::vector< std::map<int, double> > CoefMap;		///< Coefficient Map for matrix representing the ODE system
	
	Matrix<double> Eigs;								///< Matrix of eigenvectors for the coefficient matrix
	Matrix<double> invEigs;								///< Inverse matrix of eigenvectors for the coefficient matrix
	
	yaml_cpp_class *nuclides;							///< Pointer to a yaml object storing the digital library of all nuclides
	
};

///< Function to test the implementation of Isotope 
int IBIS_TESTS();

#endif /* IBIS_HPP_ */
