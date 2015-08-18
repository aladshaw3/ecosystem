//----------------------------------------
//  Created by Austin Ladshaw on 02/24/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

//	MOLA = Molecule Object Library from Atoms

#ifndef MOLA_HPP_
#define MOLA_HPP_

#include <ctype.h>
#include "eel.h"

class Molecule : Atom
{
public:
	Molecule();								//Default Constructor
	~Molecule();							//Default Destructor
	Molecule(int charge,
			 double enthalpy,
			 double entropy,
			 double energy,
			 bool HS,
			 bool G,
			 std::string Phase,
			 std::string Name,
			 std::string Formula,
			 std::string lin_formula);		//Construct any molecule from all available information
	
	void Register(int charge,
				  double enthalpy,
				  double entropy,
				  double energy,
				  bool HS,
				  bool G,
				  std::string Phase,
				  std::string Name,
				  std::string Formula,
				  std::string lin_formula);	//Register any molecule from all available information
	
	void Register(std::string formula);		//Register a known molecule from its common formula (built from library)
	
	//------------------------- NOTES ON LIBRARY REGISTRATION -----------------------
	/*		The formula is checked against a known set of molecules inside of the registration function
	*	If the formula is unknown, an error will print to the screen. Unknown molecules should be registered
	*	using the full registration function from above. The library can only be added to by a going in and
	*	editing the source code of the mola.cpp file. However, this is a relatively simple task.
	*/
	
	void recalculateMolarWeight();					//Forces molecule to recalculate its molar weight
	void setMolarWeigth(double mw);					//Set the molar weight of species to a constant
	void editCharge(int c);							//Change the ionic charge of a molecule
	void editOneOxidationState(int state,
							   std::string Symbol);	//Change oxidation state of one of the given atoms (always first atom)
	void editAllOxidationStates(int state,
								std::string Symbol);//Change the oxidation state of all the given atoms in the molecule
	void calculateAvgOxiState(std::string Symbol);	//Calculates average oxidation state of all given atoms based
													//on current oxidation states of other atoms and the molecule's charge
	void editEnthalpy(double enthalpy);				//Edit the molecules formation enthalpy
	void editEntropy(double entropy);				//Edit the molecules formation entropy
	void editHS(double H, double S);				//Edit both formation enthalpy and entropy
	void editEnergy(double energy);					//Edit Gibb's formation energy
	
	void removeOneAtom(std::string Symbol);			//Removes one atom of the symbol given (always the first atom)
	void removeAllAtoms(std::string Symbol);		//Removes all atoms of the symbol given
	
	int Charge();							//Return the charge of the molecule
	double MolarWeight();					//Return the molar weight of the molecule
	bool HaveHS();							//Returns true if enthalpy and entropy are known
	bool HaveEnergy();						//Returns true if the Gibb's energy is known
	bool isRegistered();					//Returns true if the molecule has been registered
	double Enthalpy();						//Return the formation enthalpy of the molecule
	double Entropy();						//Return the formation entropy of the molecule
	double Energy();						//Return the Gibb's formation energy of the molecule
	std::string MoleculeName();				//Return the common name of the molecule
	std::string MolecularFormula();			//Return the molecular formula of the molecule
	std::string MoleculePhase();			//Return the phase of the molecule
	
	void DisplayInfo();						//Function to display molecule information
	
protected:
	int charge;								//Ionic charge of the molecule - specified
	double molar_weight;					//Molar weight of the molecule (g/mol) - determined from atoms
	double formation_enthalpy;				//Enthalpy of formation of the molecule (J/mol) - constant
	double formation_entropy;				//Entropy of formation of the molecule (J/K/mol) - constant
	double formation_energy;				//Gibb's energy of formation (J/mol) - given
	
	std::string Phase;						//Phase of the molecule (i.e. solid, liquid, aqueous, gas...) - may change
	std::vector<Atom> atoms;				//Atoms which make up the molecule - based on Formula (constant)
	
private:
	std::string Name;						//Name of the Molecule - Common Name (i.e. H2O = Water)
	std::string Formula;					//Formula for the molecule - specified (i.e. H2O)
	bool haveG;								//True = given Gibb's energy of formation
	bool haveHS;							//True = give enthalpy and entropy of formation
	bool registered;						//True = the object was registered
	
};

//Naming convention: Formula \space charge \space phase
/*
 *						Details on naming convention
 *					-------------------------------------
 *
 *		+ and - charges are denoted by the numeric value of the charge
 *	followed by a + or - sign, respectively ( e.g. UO2(CO3)3 4- (aq) )
 *
 *		The phase is always denoted last and will be marked as (l) for
 *	liquid, (s) for solid, (aq) for aqueous, and (g) for gas (see above).
 *
 *		When registering a molecule that is not in the library, you must 
 *	also provide a linear formula during construction or registration.
 *	This is needed so that the string parsing is easier to handle when the
 *	molecule subsequently registers the necessary atoms. (e.g. UO2(CO3)3 = 
 *	UO2C3O9 or UO11C3).
 *
 *
 *						Notes on Library Molecules
 *					---------------------------------
 *
 *		These are molecules that are known and defined prior to runtime. They
 *	have a charge, energy characteristics, phase, name, and formula that they
 *	are recongized by. The formula is used to create the atoms that they are 
 *	made from. If some information is incomplete, it must be specified as to
 *	what information is missing (i.e. denote whether the formation energies are
 *	known).
 *
 *		Formation energies are used to determine stability/dissociation/acidity
 *	equilibrium constants during runtime. If the formation energies are unknown,
 *	then the equilibrium constants must be given to a reaction object on when it
 *	is initialized.
 *
 *		The molecule formula's are given as strings which are parsed in the constructor
 *	to determine what atoms from the EEL files will be registered and used. Note,
 *	you will be able to build molecules from an input file, but the library molecules
 *	here are ready to be used in applications and require no more input other that
 *	the molecule's formula. 
 */

int MOLA_TESTS();

#endif
