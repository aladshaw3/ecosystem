/*!
 *  \file mola.h mola.cpp
 *	\brief Molecule Object Library from Atoms
 *	\details This file contains a C++ Class for creating Molecule objects from the
 *			Atom objects that were defined in eel.h. Molecules can be created and
 *			registered from basic information or can be registered from a growing
 *			list of pre-registered molecules that are accessible by name/formula.
 *
 *			Registered Molecules are are known and defined prior to runtime. They
 *			have a charge, energy characteristics, phase, name, and formula that they
 *			are recongized by. The formula is used to create the atoms that they are
 *			made from. If some information is incomplete, it must be specified as to
 *			what information is missing (i.e. denote whether the formation energies are
 *			known).
 *
 *			Formation energies are used to determine stability/dissociation/acidity
 *			equilibrium constants during runtime. If the formation energies are unknown,
 *			then the equilibrium constants must be given to a reaction object on when it
 *			is initialized.
 *
 *			The molecule formula's are given as strings which are parsed in the constructor
 *			to determine what atoms from the EEL files will be registered and used. Note,
 *			you will be able to build molecules from an input file, but the library molecules
 *			here are ready to be used in applications and require no more input other that
 *			the molecule's formula.
 *
 *			List of Currently Registered Molecules
 *			--------------------------------------
 *			CO3 2- (aq) \n
 *			Cl - (aq) \n
 *			H2O (l) \n
 *			H + (aq) \n
 *			H2CO3 (aq) \n
 *			HCO3 - (aq) \n
 *			HNO3 (aq) \n
 *			HCl (aq) \n
 *			NaHCO3 (aq) \n
 *			NaCO3 - (aq) \n
 *			Na + (aq) \n
 *			NaCl (aq) \n
 *			NaOH (aq) \n
 *			NO3 - (aq) \n
 *			OH - (aq) \n
 *			UO2 2+ (aq) \n
 *			UO2NO3 + (aq) \n
 *			UO2(NO3)2 (aq) \n
 *			UO2OH + (aq) \n
 *			UO2(OH)2 (aq) \n
 *			UO2(OH)3 - (aq) \n
 *			UO2(OH)4 2- (aq) \n
 *			(UO2)2OH 3+ (aq) \n
 *			(UO2)2(OH)2 2+ (aq) \n
 *			(UO2)3(OH)4 2+ (aq) \n
 *			(UO2)3(OH)5 + (aq) \n
 *			(UO2)3(OH)7 - (aq) \n
 *			(UO2)4(OH)7 + (aq) \n
 *			UO2CO3 (aq) \n
 *			UO2(CO3)2 2- (aq) \n
 *			UO2(CO3)3 4- (aq) \n
 *
 *			Those registered molecules follow a strict naming convention by which they
 *			can be recognized (see below)...
 *
 *			Naming Convention
 *			-----------------
 *			Plus (+) and minus (-) charges are denoted by the numeric value of the charge
 *			followed by a + or - sign, respectively ( e.g. UO2(CO3)3 4- (aq) )
 *
 *			The phase is always denoted last and will be marked as (l) for
 *			liquid, (s) for solid, (aq) for aqueous, and (g) for gas (see above).
 *
 *			When registering a molecule that is not in the library, you must
 *			also provide a linear formula during construction or registration.
 *			This is needed so that the string parsing is easier to handle when the
 *			molecule subsequently registers the necessary atoms. (e.g. UO2(CO3)3 =
 *			UO2C3O9 or UO11C3).
 *
 *  \author Austin Ladshaw
 *	\date 02/24/2014
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#ifndef MOLA_HPP_
#define MOLA_HPP_

#include <ctype.h>
#include "eel.h"

typedef enum {SOLID, LIQUID, AQUEOUS, GAS, PLASMA, OTHER} valid_phase;

/// C++ Molecule Object built from Atom Objects (click Molecule to go to function definitions)
/** C++ Class Object that stores information and certain operations associated with molecules.
	Registered molecules are built up from their respective atoms so that the molecule can keep
	track of information such as molecular weigth and oxidation states. Primarily, this object
	is used in conjunction with shark.h to formulate the system of equations necessary for solving
	speciation type problems in aqueous systems. However, this object is generalized enough to be
	of use in RedOx calculations, reaction formulation, and molecular transformations.
 
	All information for a molecule should be initialized prior to performing operations with or 
	on the object. There are several molecules already defined for construction by the formulas
	listed at the top of this section. */
class Molecule
{
public:
	Molecule();								///< Default Constructor (builds an empty molecule object)
	~Molecule();							///< Default Destructor (clears out memory)
	
	/// Construct any molecule from the available information
	/** This constructor will build a user defined custom molecule. 
	 
		\param charge the ionic charge of the molecule
		\param enthalpy the standard formation enthalpy of the molecule (J/mol)
		\param entropy the standard formation entropy of the molecule (J/K/mol)
		\param energy the standard Gibb's Free Energy of formation of the molecule (J/mol)
		\param HS boolean to be set to true if enthalpy and entropy were given
		\param G boolean to be set to true if the energy was given
		\param Phase string denoting molecule's phase (i.e., Liquid, Aqueous, Gas, Solid)
		\param Name string denoting the common name of the molecule (i.e., H2O -> Water)
		\param Formula string denoting the formula by which the molecule is referened (i.e., Cl - (aq))
		\param lin_formula string denoting all the atoms in the molecule (i.e., UO2(OH)2 -> UO4H2)*/
	Molecule(int charge,
			 double enthalpy,
			 double entropy,
			 double energy,
			 bool HS,
			 bool G,
			 std::string Phase,
			 std::string Name,
			 std::string Formula,
			 std::string lin_formula);
	
	/// Function to register this molecule from the available information
	/** This function will build a user defined custom molecule.
	 
		\param charge the ionic charge of the molecule
		\param enthalpy the standard formation enthalpy of the molecule (J/mol)
		\param entropy the standard formation entropy of the molecule (J/K/mol)
		\param energy the standard Gibb's Free Energy of formation of the molecule (J/mol)
		\param HS boolean to be set to true if enthalpy and entropy were given
		\param G boolean to be set to true if the energy was given
		\param Phase string denoting molecule's phase (i.e., Liquid, Aqueous, Gas, Solid)
		\param Name string denoting the common name of the molecule (i.e., H2O -> Water)
		\param Formula string denoting the formula by which the molecule is referened (i.e., Cl - (aq))
		\param lin_formula string denoting all the atoms in the molecule (i.e., UO2(OH)2 -> UO4H2)*/
	void Register(int charge,
				  double enthalpy,
				  double entropy,
				  double energy,
				  bool HS,
				  bool G,
				  std::string Phase,
				  std::string Name,
				  std::string Formula,
				  std::string lin_formula);
	
	/// Function to register this molecule based on the given formula (if formula is in library)
	/** This function will create this molecule object from the given formula, but only if that
		formula is already registered in the library. See the top of this class section for a 
		list of all currently registered formulas. */
	
	/** \note The formula is checked against a known set of molecules inside of the registration function
	*	If the formula is unknown, an error will print to the screen. Unknown molecules should be registered
	*	using the full registration function from above. The library can only be added to by a going in and
	*	editing the source code of the mola.cpp file. However, this is a relatively simple task.
	*/
	void Register(std::string formula);
	
	void setFormula(std::string form);				///< Sets the formula for a molecule
	void recalculateMolarWeight();					///< Forces molecule to recalculate its molar weight
	void setMolarWeigth(double mw);					///< Set the molar weight of species to a constant
	void editCharge(int c);							///< Change the ionic charge of a molecule
	
	/// Change oxidation state of one of the given atoms (always first match found)
	/** This function will search the list of Atoms that make up the Molecule for the given
		atomic Symbol. It will change the oxidation state of the first found matching atom with
		the given state. */
	void editOneOxidationState(int state,
							   std::string Symbol);
	
	/// Change oxidation state of all of the given atoms
	/** This function will search the list of Atoms that make up the Molecule for the given
		atomic Symbol. It will change the oxidation state of all found matching atoms with
		the given state. */
	void editAllOxidationStates(int state,
								std::string Symbol);
	
	/// Function to calculate the average oxidation state of the atoms
	/** This function search the atoms in the molecule for the matching atomic Symbol. It
		then looks at all oxidation states of that atom in the molecule and then sets all 
		the oxidation states of that atom to the average value calculated. */
	void calculateAvgOxiState(std::string Symbol);
	
	void editEnthalpy(double enthalpy);				///< Edit the molecules formation enthalpy (J/mol)
	void editEntropy(double entropy);				///< Edit the molecules formation entropy (J/K/mol)
	
	/// Edit both formation enthalpy and entropy
	/** This function will change or set the values for formation enthalpy (J/mol) and
		formation entropy (J/K/mol) based on the given values.
		
		\param H formation enthalpy (J/mol)
		\param S formation entropy (J/K/mol)*/
	void editHS(double H, double S);
	
	void editEnergy(double energy);					///< Edit Gibb's formation energy
	
	void removeOneAtom(std::string Symbol);			///< Removes one atom of the symbol given (always the first atom found)
	void removeAllAtoms(std::string Symbol);		///< Removes all atoms of the symbol given
	
	int Charge();							///< Return the charge of the molecule
	double MolarWeight();					///< Return the molar weight of the molecule
	bool HaveHS();							///< Returns true if enthalpy and entropy are known
	bool HaveEnergy();						///< Returns true if the Gibb's energy is known
	bool isRegistered();					///< Returns true if the molecule has been registered
	double Enthalpy();						///< Return the formation enthalpy of the molecule
	double Entropy();						///< Return the formation entropy of the molecule
	double Energy();						///< Return the Gibb's formation energy of the molecule
	std::string MoleculeName();				///< Return the common name of the molecule
	std::string MolecularFormula();			///< Return the molecular formula of the molecule
	std::string MoleculePhase();			///< Return the phase of the molecule
	int MoleculePhaseID();					///< Return the enum phase ID of the molecule
	
	void DisplayInfo();						///< Function to display molecule information
	
protected:
	int charge;								///< Ionic charge of the molecule - specified
	double molar_weight;					///< Molar weight of the molecule (g/mol) - determined from atoms or specified
	double formation_enthalpy;				///< Enthalpy of formation of the molecule (J/mol) - constant
	double formation_entropy;				///< Entropy of formation of the molecule (J/K/mol) - constant
	double formation_energy;				///< Gibb's energy of formation (J/mol) - given
	
	std::string Phase;						///< Phase of the molecule (i.e. Solid, Liquid, Aqueous, Gas...)
	int PhaseID;							///< Phase ID of the molecule (from the enum)
	std::vector<Atom> atoms;				///< Atoms which make up the molecule - based on Formula
	
private:
	std::string Name;						///< Name of the Molecule - Common Name (i.e. H2O = Water)
	std::string Formula;					///< Formula for the molecule - specified (i.e. H2O)
	bool haveG;								///< True = given Gibb's energy of formation
	bool haveHS;							///< True = give enthalpy and entropy of formation
	bool registered;						///< True = the object was registered
	
};

/// Function to run the MOLA tests
/** This function is callable from the UI and is used to run several
	algorithm tests for the Molecule objects. This test should never
	report any errors. */
int MOLA_TESTS();

#endif
