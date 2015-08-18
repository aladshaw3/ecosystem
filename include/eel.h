//----------------------------------------
//  Created by Austin Ladshaw on 02/23/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

//	EEL = Easy-access Element Library

#ifndef EEL_HPP_
#define EEL_HPP_

#include <stdio.h>				//Line to allow cout functionality
#include <math.h>               //Line added to allow usage of the pow (e, x) function
#include <iostream>				//Line to allow for read/write to the console using cpp functions
#include <fstream>				//Line to allow for read/write to and from .txt files
#include <stdlib.h>				//Line need to convert strings to doubles
#include <vector>				//Line needed to use dynamic arrays called vectors
#include <time.h>				//Line needed to display program runtime
#include <float.h>				//Line to allow use of machine precision constants
#include <string>				//Line to allow use of strings as a data type
#include "error.h"				//Line to allow use of the custom error file

class Atom
{
public:
	Atom();									//Default Constructor
	~Atom();								//Default Destructor
	Atom(std::string Name);					//Constructor by Atom Name
	Atom(int number);						//Constructor by Atomic number
	
	void Register(std::string Symbol);		//Register an atom object by symbol
	void Register(int number);				//Register an atom object by number
	
	void editAtomicWeight(double AW);		//Manually changes the atomic weight
	void editOxidationState(int state);		//Manually changes the oxidation state
	void editProtons(int proton);			//Manually changes the number of protons
	void editNeutrons(int neutron);			//Manually changes the number of neutrons
	void editElectrons(int electron);		//Manually changes the number of electrons
	void editValence(int val);				//Manually changes the number of valence electrons
	
	void removeProton();					//Manually removes 1 proton and adjusts weight
	void removeNeutron();					//Manually removes 1 neutron and adjusts weight
	void removeElectron();					//Manually removes 1 electron from valence
	
	double AtomicWeight();					//Returns the current atomic weight (g/mol)
	int OxidationState();					//Returns the current oxidation state
	int Protons();							//Returns the current number of protons
	int Neutrons();							//Returns the current number of neutrons
	int Electrons();						//Returns the current number of electrons
	int BondingElectrons();					//Returns the number of electrons available for bonding
	
	std::string AtomName();					//Returns the name of the atom
	std::string AtomSymbol();				//Returns the symbol of the atom
	std::string AtomCategory();				//Returns the category of the atom
	std::string AtomState();				//Returns the state of the atom
	int AtomicNumber();						//Returns the atomic number of the atom
	
	void DisplayInfo();						//Displays Atom information to console
	
protected:
	double atomic_weight;
	int oxidation_state;
	int protons;
	int neutrons;
	int electrons;
	int valence_e;
	
private:
	std::string Name;
	std::string Symbol;
	std::string Category;
	std::string NaturalState;
	int atomic_number;
	
};

class PeriodicTable : Atom
{
public:
	PeriodicTable();								//Default Constructor - Build Perodic Table
	~PeriodicTable();								//Default Destructor - Destroy the table
	PeriodicTable(int *n, int N);					//Construct a partial table from a list of atomic numbers
	PeriodicTable(std::vector<std::string> &Symbol);//Construct a partial table from a vector of atom symbols
	PeriodicTable(std::vector<int> &n);				//Construct a partial table from a vector of atomic numbers
	
	void DisplayTable();							//Displays the periodic table via symbols
	
protected:
	std::vector<Atom> Table;
	
private:
	int number_elements;
	
};

int EEL_TESTS();

#endif
