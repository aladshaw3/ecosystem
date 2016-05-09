/*!
 *  \file eel.cpp eel.h
 *	\brief Easy-access Element Library
 *  \author Austin Ladshaw
 *	\date 02/23/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */
#include "eel.h"

//Default Constructor
Atom::Atom()
:
atomic_weight(0.0),
oxidation_state(0),
protons(0),
neutrons(0),
electrons(0),
valence_e(0),
Name("No Name"),
Symbol("N/A"),
Category("N/A"),
NaturalState("N/A"),
atomic_number(0),
atomic_radii(0)
{
}

//Default Destructor
Atom::~Atom()
{
}

//Constructor by Name
Atom::Atom(std::string name)
{
	if (name.compare("Hydrogen") == 0)
	{
		atomic_weight = 1.0081;
		oxidation_state = 1;
		protons = 1;
		neutrons = 0;
		electrons = 1;
		valence_e = 1;
		Name = "Hydrogen";
		Symbol = "H";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 1;
		atomic_radii = 1.2;
	}
	else if (name.compare("Helium") == 0)
	{
		atomic_weight = 4.0026022;
		oxidation_state = 0;
		protons = 2;
		neutrons = 2;
		electrons = 2;
		valence_e = 2;
		Name = "Helium";
		Symbol = "He";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 2;
		atomic_radii = 1.4;
	}
	else if (name.compare("Lithium") == 0)
	{
		atomic_weight = 6.941;
		oxidation_state = 1;
		protons = 3;
		neutrons = 4;
		electrons = 3;
		valence_e = 1;
		Name = "Lithium";
		Symbol = "Li";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 3;
		atomic_radii = 1.82;
	}
	else if (name.compare("Beryllium") == 0)
	{
		atomic_weight = 9.01218315;
		oxidation_state = 2;
		protons = 4;
		neutrons = 5;
		electrons = 4;
		valence_e = 2;
		Name = "Beryllium";
		Symbol = "Be";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 4;
		atomic_radii = 1.53;
	}
	else if (name.compare("Boron") == 0)
	{
		atomic_weight = 10.811;
		oxidation_state = 3;
		protons = 5;
		neutrons = 6;
		electrons = 5;
		valence_e = 3;
		Name = "Boron";
		Symbol = "B";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 5;
		atomic_radii = 1.92;
	}
	else if (name.compare("Carbon") == 0)
	{
		atomic_weight = 12.0111;
		oxidation_state = 4;
		protons = 6;
		neutrons = 6;
		electrons = 6;
		valence_e = 4;
		Name = "Carbon";
		Symbol = "C";
		Category = "Polyatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 6;
		atomic_radii = 1.70;
	}
	else if (name.compare("Nitrogen") == 0)
	{
		atomic_weight = 14.0071;
		oxidation_state = -3;
		protons = 7;
		neutrons = 7;
		electrons = 7;
		valence_e = 5;
		Name = "Nitrogen";
		Symbol = "N";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 7;
		atomic_radii = 1.55;
	}
	else if (name.compare("Oxygen") == 0)
	{
		atomic_weight = 15.9994;
		oxidation_state = -2;
		protons = 8;
		neutrons = 8;
		electrons = 8;
		valence_e = 6;
		Name = "Oxygen";
		Symbol = "O";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 8;
		atomic_radii = 1.52;
	}
	else if (name.compare("Fluorine") == 0)
	{
		atomic_weight = 18.9984031636;
		oxidation_state = -1;
		protons = 9;
		neutrons = 10;
		electrons = 9;
		valence_e = 7;
		Name = "Fluorine";
		Symbol = "F";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 9;
		atomic_radii = 1.35;
	}
	else if (name.compare("Neon") == 0)
	{
		atomic_weight = 20.17976;
		oxidation_state = 0;
		protons = 10;
		neutrons = 10;
		electrons = 10;
		valence_e = 8;
		Name = "Neon";
		Symbol = "Ne";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 10;
		atomic_radii = 1.54;
	}
	else if (name.compare("Sodium") == 0)
	{
		atomic_weight = 22.989769282;
		oxidation_state = 1;
		protons = 11;
		neutrons = 12;
		electrons = 11;
		valence_e = 1;
		Name = "Sodium";
		Symbol = "Na";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 11;
		atomic_radii = 2.27;
	}
	else if (name.compare("Magnesium") == 0)
	{
		atomic_weight = 24.3051;
		oxidation_state = 2;
		protons = 12;
		neutrons = 12;
		electrons = 12;
		valence_e = 2;
		Name = "Magnesium";
		Symbol = "Mg";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 12;
		atomic_radii = 1.73;
	}
	else if (name.compare("Aluminium") == 0)
	{
		atomic_weight = 26.98153857;
		oxidation_state = 3;
		protons = 13;
		neutrons = 14;
		electrons = 13;
		valence_e = 3;
		Name = "Aluminium";
		Symbol = "Al";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 13;
		atomic_radii = 1.84;
	}
	else if (name.compare("Silicon") == 0)
	{
		atomic_weight = 28.0851;
		oxidation_state = 4;
		protons = 14;
		neutrons = 14;
		electrons = 14;
		valence_e = 4;
		Name = "Silicon";
		Symbol = "Si";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 14;
		atomic_radii = 2.10;
	}
	else if (name.compare("Phosphorus") == 0)
	{
		atomic_weight = 30.9737619985;
		oxidation_state = 5;
		protons = 15;
		neutrons = 16;
		electrons = 15;
		valence_e = 5;
		Name = "Phosphorus";
		Symbol = "P";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 15;
		atomic_radii = 1.80;
	}
	else if (name.compare("Sulfur") == 0)
	{
		atomic_weight = 32.061;
		oxidation_state = 6;
		protons = 16;
		neutrons = 16;
		electrons = 16;
		valence_e = 6;
		Name = "Sulfur";
		Symbol = "S";
		Category = "Polyatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 16;
		atomic_radii = 1.80;
	}
	else if (name.compare("Chlorine") == 0)
	{
		atomic_weight = 35.451;
		oxidation_state = -1;
		protons = 17;
		neutrons = 18;
		electrons = 17;
		valence_e = 7;
		Name = "Chlorine";
		Symbol = "Cl";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 17;
		atomic_radii = 1.75;
	}
	else if (name.compare("Argon") == 0)
	{
		atomic_weight = 39.9481;
		oxidation_state = 0;
		protons = 18;
		neutrons = 22;
		electrons = 18;
		valence_e = 8;
		Name = "Argon";
		Symbol = "Ar";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 18;
		atomic_radii = 1.88;
	}
	else if (name.compare("Potassium") == 0)
	{
		atomic_weight = 39.09831;
		oxidation_state = 1;
		protons = 19;
		neutrons = 20;
		electrons = 19;
		valence_e = 1;
		Name = "Potassium";
		Symbol = "K";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 19;
		atomic_radii = 2.75;
	}
	else if (name.compare("Calcium") == 0)
	{
		atomic_weight = 40.0784;
		oxidation_state = 2;
		protons = 20;
		neutrons = 20;
		electrons = 20;
		valence_e = 2;
		Name = "Calcium";
		Symbol = "Ca";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 20;
		atomic_radii = 2.31;
	}
	else if (name.compare("Scandium") == 0)
	{
		atomic_weight = 44.9559085;
		oxidation_state = 3;
		protons = 21;
		neutrons = 24;
		electrons = 21;
		valence_e = 3;
		Name = "Scandium";
		Symbol = "Sc";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 21;
		atomic_radii = 2.11;
	}
	else if (name.compare("Titanium") == 0)
	{
		atomic_weight = 47.8671;
		oxidation_state = 4;
		protons = 22;
		neutrons = 26;
		electrons = 22;
		valence_e = 4;
		Name = "Titanium";
		Symbol = "Ti";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 22;
		atomic_radii = 1.60;
	}
	else if (name.compare("Vanadium") == 0)
	{
		atomic_weight = 50.94151;
		oxidation_state = 5;
		protons = 23;
		neutrons = 28;
		electrons = 23;
		valence_e = 5;
		Name = "Vanadium";
		Symbol = "V";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 23;
		atomic_radii = 1.53;
	}
	else if (name.compare("Chromium") == 0)
	{
		atomic_weight = 51.99616;
		oxidation_state = 6;
		protons = 24;
		neutrons = 28;
		electrons = 24;
		valence_e = 6;
		Name = "Chromium";
		Symbol = "Cr";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 24;
		atomic_radii = 1.39;
	}
	else if (name.compare("Manganese") == 0)
	{
		atomic_weight = 54.9380443;
		oxidation_state = 7;
		protons = 25;
		neutrons = 30;
		electrons = 25;
		valence_e = 7;
		Name = "Manganese";
		Symbol = "Mn";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 25;
		atomic_radii = 1.39;
	}
	else if (name.compare("Iron") == 0)
	{
		atomic_weight = 55.8452;
		oxidation_state = 3;
		protons = 26;
		neutrons = 30;
		electrons = 26;
		valence_e = 8;
		Name = "Iron";
		Symbol = "Fe";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 26;
		atomic_radii = 1.32;
	}
	else if (name.compare("Cobalt") == 0)
	{
		atomic_weight = 58.9331944;
		oxidation_state = 3;
		protons = 27;
		neutrons = 32;
		electrons = 27;
		valence_e = 9;
		Name = "Cobalt";
		Symbol = "Co";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 27;
		atomic_radii = 1.26;
	}
	else if (name.compare("Nickel") == 0)
	{
		atomic_weight = 58.69344;
		oxidation_state = 2;
		protons = 28;
		neutrons = 31;
		electrons = 28;
		valence_e = 10;
		Name = "Nickel";
		Symbol = "Ni";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 28;
		atomic_radii = 1.63;
	}
	else if (name.compare("Copper") == 0)
	{
		atomic_weight = 63.5463;
		oxidation_state = 2;
		protons = 29;
		neutrons = 35;
		electrons = 29;
		valence_e = 11;
		Name = "Copper";
		Symbol = "Cu";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 29;
		atomic_radii = 1.40;
	}
	else if (name.compare("Zinc") == 0)
	{
		atomic_weight = 65.382;
		oxidation_state = 2;
		protons = 30;
		neutrons = 35;
		electrons = 30;
		valence_e = 12;
		Name = "Zinc";
		Symbol = "Zn";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 30;
		atomic_radii = 1.39;
	}
	else if (name.compare("Gallium") == 0)
	{
		atomic_weight = 69.7231;
		oxidation_state = 3;
		protons = 31;
		neutrons = 39;
		electrons = 31;
		valence_e = 3;
		Name = "Gallium";
		Symbol = "Ga";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 31;
		atomic_radii = 1.87;
	}
	else if (name.compare("Germanium") == 0)
	{
		atomic_weight = 72.6308;
		oxidation_state = 4;
		protons = 32;
		neutrons = 41;
		electrons = 32;
		valence_e = 4;
		Name = "Germanium";
		Symbol = "Ge";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 32;
		atomic_radii = 2.11;
	}
	else if (name.compare("Arsenic") == 0)
	{
		atomic_weight = 74.9215956;
		oxidation_state = 5;
		protons = 33;
		neutrons = 42;
		electrons = 33;
		valence_e = 5;
		Name = "Arsenic";
		Symbol = "As";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 33;
		atomic_radii = 1.85;
	}
	else if (name.compare("Selenium") == 0)
	{
		atomic_weight = 78.9718;
		oxidation_state = 6;
		protons = 34;
		neutrons = 45;
		electrons = 34;
		valence_e = 6;
		Name = "Selenium";
		Symbol = "Se";
		Category = "Polyatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 34;
		atomic_radii = 1.90;
	}
	else if (name.compare("Bromine") == 0)
	{
		atomic_weight = 79.9041;
		oxidation_state = -1;
		protons = 35;
		neutrons = 45;
		electrons = 35;
		valence_e = 7;
		Name = "Bromine";
		Symbol = "Br";
		Category = "Diatomic Non-metal";
		NaturalState = "Liquid";
		atomic_number = 35;
		atomic_radii = 1.85;
	}
	else if (name.compare("Krypton") == 0)
	{
		atomic_weight = 83.798;
		oxidation_state = 0;
		protons = 36;
		neutrons = 48;
		electrons = 36;
		valence_e = 8;
		Name = "Krypton";
		Symbol = "Kr";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 36;
		atomic_radii = 2.02;
	}
	else if (name.compare("Rubidium") == 0)
	{
		atomic_weight = 85.46783;
		oxidation_state = 1;
		protons = 37;
		neutrons = 48;
		electrons = 37;
		valence_e = 1;
		Name = "Rubidium";
		Symbol = "Rb";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 37;
		atomic_radii = 3.03;
	}
	else if (name.compare("Strontium") == 0)
	{
		atomic_weight = 87.621;
		oxidation_state = 2;
		protons = 38;
		neutrons = 50;
		electrons = 38;
		valence_e = 2;
		Name = "Strontium";
		Symbol = "Sr";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 38;
		atomic_radii = 2.49;
	}
	else if (name.compare("Yttrium") == 0)
	{
		atomic_weight = 88.905842;
		oxidation_state = 3;
		protons = 39;
		neutrons = 50;
		electrons = 39;
		valence_e = 3;
		Name = "Yttrium";
		Symbol = "Y";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 39;
		atomic_radii = 1.90;
	}
	else if (name.compare("Zirconium") == 0)
	{
		atomic_weight = 91.2242;
		oxidation_state = 4;
		protons = 40;
		neutrons = 51;
		electrons = 40;
		valence_e = 4;
		Name = "Zirconium";
		Symbol = "Zr";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 40;
		atomic_radii = 1.75;
	}
	else if (name.compare("Niobium") == 0)
	{
		atomic_weight = 92.906372;
		oxidation_state = 5;
		protons = 41;
		neutrons = 52;
		electrons = 41;
		valence_e = 5;
		Name = "Niobium";
		Symbol = "Nb";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 41;
		atomic_radii = 1.64;
	}
	else if (name.compare("Molybdenum") == 0)
	{
		atomic_weight = 95.951;
		oxidation_state = 6;
		protons = 42;
		neutrons = 54;
		electrons = 42;
		valence_e = 6;
		Name = "Molybdenum";
		Symbol = "Mo";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 42;
		atomic_radii = 1.54;
	}
	else if (name.compare("Technetium") == 0)
	{
		atomic_weight = 98.0;
		oxidation_state = 7;
		protons = 43;
		neutrons = 55;
		electrons = 43;
		valence_e = 7;
		Name = "Technetium";
		Symbol = "Tc";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 43;
		atomic_radii = 1.47;
	}
	else if (name.compare("Ruthenium") == 0)
	{
		atomic_weight = 101.072;
		oxidation_state = 4;
		protons = 44;
		neutrons = 57;
		electrons = 44;
		valence_e = 8;
		Name = "Ruthenium";
		Symbol = "Ru";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 44;
		atomic_radii = 1.46;
	}
	else if (name.compare("Rhodium") == 0)
	{
		atomic_weight = 102.905502;
		oxidation_state = 3;
		protons = 45;
		neutrons = 58;
		electrons = 45;
		valence_e = 9;
		Name = "Rhodium";
		Symbol = "Rh";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 45;
		atomic_radii = 1.42;
	}
	else if (name.compare("Palladium") == 0)
	{
		atomic_weight = 106.421;
		oxidation_state = 4;
		protons = 46;
		neutrons = 60;
		electrons = 46;
		valence_e = 10;
		Name = "Palladium";
		Symbol = "Pd";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 46;
		atomic_radii = 1.63;
	}
	else if (name.compare("Silver") == 0)
	{
		atomic_weight = 107.86822;
		oxidation_state = 1;
		protons = 47;
		neutrons = 61;
		electrons = 47;
		valence_e = 11;
		Name = "Silver";
		Symbol = "Ag";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 47;
		atomic_radii = 1.72;
	}
	else if (name.compare("Cadmium") == 0)
	{
		atomic_weight = 112.4144;
		oxidation_state = 2;
		protons = 48;
		neutrons = 64;
		electrons = 48;
		valence_e = 12;
		Name = "Cadmium";
		Symbol = "Cd";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 48;
		atomic_radii = 1.58;
	}
	else if (name.compare("Indium") == 0)
	{
		atomic_weight = 114.8181;
		oxidation_state = 3;
		protons = 49;
		neutrons = 66;
		electrons = 49;
		valence_e = 3;
		Name = "Indium";
		Symbol = "In";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 49;
		atomic_radii = 1.93;
	}
	else if (name.compare("Tin") == 0)
	{
		atomic_weight = 118.7107;
		oxidation_state = 4;
		protons = 50;
		neutrons = 69;
		electrons = 50;
		valence_e = 4;
		Name = "Tin";
		Symbol = "Sn";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 50;
		atomic_radii = 2.17;
	}
	else if (name.compare("Antimony") == 0)
	{
		atomic_weight = 121.7601;
		oxidation_state = 5;
		protons = 51;
		neutrons = 71;
		electrons = 51;
		valence_e = 5;
		Name = "Antimony";
		Symbol = "Sb";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 51;
		atomic_radii = 2.06;
	}
	else if (name.compare("Tellurium") == 0)
	{
		atomic_weight = 127.603;
		oxidation_state = 6;
		protons = 52;
		neutrons = 76;
		electrons = 52;
		valence_e = 6;
		Name = "Tellurium";
		Symbol = "Te";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 52;
		atomic_radii = 2.06;
	}
	else if (name.compare("Iodine") == 0)
	{
		atomic_weight = 126.904473;
		oxidation_state = -1;
		protons = 53;
		neutrons = 74;
		electrons = 53;
		valence_e = 7;
		Name = "Iodine";
		Symbol = "I";
		Category = "Diatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 53;
		atomic_radii = 1.98;
	}
	else if (name.compare("Xenon") == 0)
	{
		atomic_weight = 131.2936;
		oxidation_state = 0;
		protons = 54;
		neutrons = 77;
		electrons = 54;
		valence_e = 8;
		Name = "Xenon";
		Symbol = "Xe";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 54;
		atomic_radii = 2.16;
	}
	else if (name.compare("Caesium") == 0)
	{
		atomic_weight = 132.905451966;
		oxidation_state = 1;
		protons = 55;
		neutrons = 78;
		electrons = 55;
		valence_e = 1;
		Name = "Caesium";
		Symbol = "Cs";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 55;
		atomic_radii = 3.43;
	}
	else if (name.compare("Barium") == 0)
	{
		atomic_weight = 137.3277;
		oxidation_state = 2;
		protons = 56;
		neutrons = 81;
		electrons = 56;
		valence_e = 2;
		Name = "Barium";
		Symbol = "Ba";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 56;
		atomic_radii = 2.68;
	}
	else if (name.compare("Lanthanum") == 0)
	{
		atomic_weight = 138.90547;
		oxidation_state = 3;
		protons = 57;
		neutrons = 82;
		electrons = 57;
		valence_e = 3;
		Name = "Lanthanum";
		Symbol = "La";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 57;
		atomic_radii = 2.07;
	}
	else if (name.compare("Cerium") == 0)
	{
		atomic_weight = 140.116;
		oxidation_state = 4;
		protons = 58;
		neutrons = 82;
		electrons = 58;
		valence_e = 4;
		Name = "Cerium";
		Symbol = "Ce";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 58;
		atomic_radii = 2.04;
	}
	else if (name.compare("Praseodymium") == 0)
	{
		atomic_weight = 140.907662;
		oxidation_state = 4;
		protons = 59;
		neutrons = 82;
		electrons = 59;
		valence_e = 5;
		Name = "Praseodymium";
		Symbol = "Pr";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 59;
		atomic_radii = 2.03;
	}
	else if (name.compare("Neodymium") == 0)
	{
		atomic_weight = 144.242;
		oxidation_state = 3;
		protons = 60;
		neutrons = 84;
		electrons = 60;
		valence_e = 6;
		Name = "Neodymium";
		Symbol = "Nd";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 60;
		atomic_radii = 2.01;
	}
	else if (name.compare("Promethium") == 0)
	{
		atomic_weight = 145.0;
		oxidation_state = 3;
		protons = 61;
		neutrons = 84;
		electrons = 61;
		valence_e = 7;
		Name = "Promethium";
		Symbol = "Pm";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 61;
		atomic_radii = 1.99;
	}
	else if (name.compare("Samarium") == 0)
	{
		atomic_weight = 150.362;
		oxidation_state = 3;
		protons = 62;
		neutrons = 88;
		electrons = 62;
		valence_e = 8;
		Name = "Samarium";
		Symbol = "Sm";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 62;
		atomic_radii = 1.98;
	}
	else if (name.compare("Europium") == 0)
	{
		atomic_weight = 151.964;
		oxidation_state = 3;
		protons = 63;
		neutrons = 89;
		electrons = 63;
		valence_e = 9;
		Name = "Europium";
		Symbol = "Eu";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 63;
		atomic_radii = 1.98;
	}
	else if (name.compare("Gadolinium") == 0)
	{
		atomic_weight = 157.253;
		oxidation_state = 3;
		protons = 64;
		neutrons = 93;
		electrons = 64;
		valence_e = 10;
		Name = "Gadolinium";
		Symbol = "Gd";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 64;
		atomic_radii = 1.96;
	}
	else if (name.compare("Terbium") == 0)
	{
		atomic_weight = 158.92535;
		oxidation_state = 3;
		protons = 65;
		neutrons = 94;
		electrons = 65;
		valence_e = 11;
		Name = "Terbium";
		Symbol = "Tb";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 65;
		atomic_radii = 1.94;
	}
	else if (name.compare("Dysprosium") == 0)
	{
		atomic_weight = 162.5001;
		oxidation_state = 3;
		protons = 66;
		neutrons = 97;
		electrons = 66;
		valence_e = 12;
		Name = "Dysprosium";
		Symbol = "Dy";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 66;
		atomic_radii = 1.92;
	}
	else if (name.compare("Holmium") == 0)
	{
		atomic_weight = 164.930332;
		oxidation_state = 3;
		protons = 67;
		neutrons = 98;
		electrons = 67;
		valence_e = 13;
		Name = "Holmium";
		Symbol = "Ho";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 67;
		atomic_radii = 1.92;
	}
	else if (name.compare("Erbium") == 0)
	{
		atomic_weight = 167.259;
		oxidation_state = 3;
		protons = 68;
		neutrons = 99;
		electrons = 68;
		valence_e = 14;
		Name = "Erbium";
		Symbol = "Er";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 68;
		atomic_radii = 1.89;
	}
	else if (name.compare("Thulium") == 0)
	{
		atomic_weight = 168.934222;
		oxidation_state = 3;
		protons = 69;
		neutrons = 100;
		electrons = 69;
		valence_e = 15;
		Name = "Thulium";
		Symbol = "Tm";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 69;
		atomic_radii = 1.90;
	}
	else if (name.compare("Ytterbium") == 0)
	{
		atomic_weight = 173.0545;
		oxidation_state = 3;
		protons = 70;
		neutrons = 103;
		electrons = 70;
		valence_e = 16;
		Name = "Ytterbium";
		Symbol = "Yb";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 70;
		atomic_radii = 1.87;
	}
	else if (name.compare("Lutetium") == 0)
	{
		atomic_weight = 174.96684;
		oxidation_state = 3;
		protons = 71;
		neutrons = 104;
		electrons = 71;
		valence_e = 3;
		Name = "Lutetium";
		Symbol = "Lu";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 71;
		atomic_radii = 1.87;
	}
	else if (name.compare("Hafnium") == 0)
	{
		atomic_weight = 178.492;
		oxidation_state = 4;
		protons = 72;
		neutrons = 106;
		electrons = 72;
		valence_e = 4;
		Name = "Hafnium";
		Symbol = "Hf";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 72;
		atomic_radii = 1.75;
	}
	else if (name.compare("Tantalum") == 0)
	{
		atomic_weight = 180.947882;
		oxidation_state = 5;
		protons = 73;
		neutrons = 108;
		electrons = 73;
		valence_e = 5;
		Name = "Tantalum";
		Symbol = "Ta";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 73;
		atomic_radii = 1.70;
	}
	else if (name.compare("Tungsten") == 0)
	{
		atomic_weight = 183.841;
		oxidation_state = 6;
		protons = 74;
		neutrons = 110;
		electrons = 74;
		valence_e = 6;
		Name = "Tungsten";
		Symbol = "W";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 74;
		atomic_radii = 1.62;
	}
	else if (name.compare("Rhenium") == 0)
	{
		atomic_weight = 186.2071;
		oxidation_state = 7;
		protons = 75;
		neutrons = 111;
		electrons = 75;
		valence_e = 7;
		Name = "Rhenium";
		Symbol = "Re";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 75;
		atomic_radii = 1.51;
	}
	else if (name.compare("Osmium") == 0)
	{
		atomic_weight = 190.233;
		oxidation_state = 4;
		protons = 76;
		neutrons = 114;
		electrons = 76;
		valence_e = 8;
		Name = "Osmium";
		Symbol = "Os";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 76;
		atomic_radii = 1.44;
	}
	else if (name.compare("Iridium") == 0)
	{
		atomic_weight = 192.2173;
		oxidation_state = 4;
		protons = 77;
		neutrons = 115;
		electrons = 77;
		valence_e = 9;
		Name = "Iridium";
		Symbol = "Ir";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 77;
		atomic_radii = 1.41;
	}
	else if (name.compare("Platinum") == 0)
	{
		atomic_weight = 195.0849;
		oxidation_state = 4;
		protons = 78;
		neutrons = 117;
		electrons = 78;
		valence_e = 10;
		Name = "Platinum";
		Symbol = "Pt";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 78;
		atomic_radii = 1.75;
	}
	else if (name.compare("Gold") == 0)
	{
		atomic_weight = 196.9665694;
		oxidation_state = 3;
		protons = 79;
		neutrons = 118;
		electrons = 79;
		valence_e = 11;
		Name = "Gold";
		Symbol = "Au";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 79;
		atomic_radii = 1.66;
	}
	else if (name.compare("Mercury") == 0)
	{
		atomic_weight = 200.5923;
		oxidation_state = 2;
		protons = 80;
		neutrons = 121;
		electrons = 80;
		valence_e = 12;
		Name = "Mercury";
		Symbol = "Hg";
		Category = "Transition Metal";
		NaturalState = "Liquid";
		atomic_number = 80;
		atomic_radii = 1.55;
	}
	else if (name.compare("Thallium") == 0)
	{
		atomic_weight = 204.381;
		oxidation_state = 1;
		protons = 81;
		neutrons = 123;
		electrons = 81;
		valence_e = 3;
		Name = "Thallium";
		Symbol = "Tl";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 81;
		atomic_radii = 1.96;
	}
	else if (name.compare("Lead") == 0)
	{
		atomic_weight = 207.21;
		oxidation_state = 2;
		protons = 82;
		neutrons = 125;
		electrons = 82;
		valence_e = 4;
		Name = "Lead";
		Symbol = "Pb";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 82;
		atomic_radii = 2.02;
	}
	else if (name.compare("Bismuth") == 0)
	{
		atomic_weight = 208.980401;
		oxidation_state = 3;
		protons = 83;
		neutrons = 126;
		electrons = 83;
		valence_e = 5;
		Name = "Bismuth";
		Symbol = "Bi";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 83;
		atomic_radii = 2.07;
	}
	else if (name.compare("Polonium") == 0)
	{
		atomic_weight = 209.0;
		oxidation_state = 4;
		protons = 84;
		neutrons = 125;
		electrons = 84;
		valence_e = 6;
		Name = "Polonium";
		Symbol = "Po";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 84;
		atomic_radii = 1.97;
	}
	else if (name.compare("Astatine") == 0)
	{
		atomic_weight = 210.0;
		oxidation_state = -1;
		protons = 85;
		neutrons = 125;
		electrons = 85;
		valence_e = 7;
		Name = "Astatine";
		Symbol = "At";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 85;
		atomic_radii = 2.02;
	}
	else if (name.compare("Radon") == 0)
	{
		atomic_weight = 222.0;
		oxidation_state = 0;
		protons = 86;
		neutrons = 136;
		electrons = 86;
		valence_e = 8;
		Name = "Radon";
		Symbol = "Rn";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 86;
		atomic_radii = 2.20;
	}
	else if (name.compare("Francium") == 0)
	{
		atomic_weight = 223.0;
		oxidation_state = 1;
		protons = 87;
		neutrons = 136;
		electrons = 87;
		valence_e = 1;
		Name = "Francium";
		Symbol = "Fr";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 87;
		atomic_radii = 3.48;
	}
	else if (name.compare("Radium") == 0)
	{
		atomic_weight = 226.0;
		oxidation_state = 2;
		protons = 88;
		neutrons = 138;
		electrons = 88;
		valence_e = 2;
		Name = "Radium";
		Symbol = "Ra";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 88;
		atomic_radii = 2.83;
	}
	else if (name.compare("Actinium") == 0)
	{
		atomic_weight = 227.0;
		oxidation_state = 3;
		protons = 89;
		neutrons = 138;
		electrons = 89;
		valence_e = 3;
		Name = "Actinium";
		Symbol = "Ac";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 89;
		atomic_radii = 2.15;
	}
	else if (name.compare("Thorium") == 0)
	{
		atomic_weight = 232.03774;
		oxidation_state = 4;
		protons = 90;
		neutrons = 142;
		electrons = 90;
		valence_e = 4;
		Name = "Thorium";
		Symbol = "Th";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 90;
		atomic_radii = 2.06;
	}
	else if (name.compare("Protactinium") == 0)
	{
		atomic_weight = 231.03588;
		oxidation_state = 5;
		protons = 91;
		neutrons = 140;
		electrons = 91;
		valence_e = 5;
		Name = "Protactinium";
		Symbol = "Pa";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 91;
		atomic_radii = 2.00;
	}
	else if (name.compare("Uranium") == 0)
	{
		atomic_weight = 238.028913;
		oxidation_state = 6;
		protons = 92;
		neutrons = 146;
		electrons = 92;
		valence_e = 6;
		Name = "Uranium";
		Symbol = "U";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 92;
		atomic_radii = 1.86;
	}
	else if (name.compare("Neptunium") == 0)
	{
		atomic_weight = 237.0;
		oxidation_state = 5;
		protons = 93;
		neutrons = 144;
		electrons = 93;
		valence_e = 7;
		Name = "Neptunium";
		Symbol = "Np";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 93;
		atomic_radii = 1.90;
	}
	else if (name.compare("Plutonium") == 0)
	{
		atomic_weight = 244.0;
		oxidation_state = 4;
		protons = 94;
		neutrons = 150;
		electrons = 94;
		valence_e = 8;
		Name = "Plutonium";
		Symbol = "Pu";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 94;
		atomic_radii = 1.87;
	}
	else if (name.compare("Americium") == 0)
	{
		atomic_weight = 243.0;
		oxidation_state = 3;
		protons = 95;
		neutrons = 148;
		electrons = 95;
		valence_e = 9;
		Name = "Americium";
		Symbol = "Am";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 95;
		atomic_radii = 1.80;
	}
	else if (name.compare("Curium") == 0)
	{
		atomic_weight = 247.0;
		oxidation_state = 3;
		protons = 96;
		neutrons = 151;
		electrons = 96;
		valence_e = 10;
		Name = "Curium";
		Symbol = "Cm";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 96;
		atomic_radii = 1.69;
	}
	else if (name.compare("Berkelium") == 0)
	{
		atomic_weight = 247.0;
		oxidation_state = 3;
		protons = 97;
		neutrons = 150;
		electrons = 97;
		valence_e = 11;
		Name = "Berkelium";
		Symbol = "Bk";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 97;
		atomic_radii = 1.70;
	}
	else if (name.compare("Californium") == 0)
	{
		atomic_weight = 251.0;
		oxidation_state = 3;
		protons = 98;
		neutrons = 153;
		electrons = 98;
		valence_e = 12;
		Name = "Californium";
		Symbol = "Cf";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 98;
		atomic_radii = 1.70;
	}
	else if (name.compare("Einsteinium") == 0)
	{
		atomic_weight = 252.0;
		oxidation_state = 3;
		protons = 99;
		neutrons = 153;
		electrons = 99;
		valence_e = 13;
		Name = "Einsteinium";
		Symbol = "Es";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 99;
		atomic_radii = 1.70;
	}
	else if (name.compare("Fermium") == 0)
	{
		atomic_weight = 257.0;
		oxidation_state = 3;
		protons = 100;
		neutrons = 157;
		electrons = 100;
		valence_e = 14;
		Name = "Fermium";
		Symbol = "Fm";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 100;
		atomic_radii = 1.70;
	}
	else if (name.compare("Mendelevium") == 0)
	{
		atomic_weight = 258.0;
		oxidation_state = 3;
		protons = 101;
		neutrons = 157;
		electrons = 101;
		valence_e = 15;
		Name = "Mendelevium";
		Symbol = "Md";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 101;
		atomic_radii = 1.70;
	}
	else if (name.compare("Nobelium") == 0)
	{
		atomic_weight = 259.0;
		oxidation_state = 2;
		protons = 102;
		neutrons = 157;
		electrons = 102;
		valence_e = 16;
		Name = "Nobelium";
		Symbol = "No";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 102;
		atomic_radii = 1.70;
	}
	else if (name.compare("Lawrencium") == 0)
	{
		atomic_weight = 266.0;
		oxidation_state = 3;
		protons = 103;
		neutrons = 159;
		electrons = 103;
		valence_e = 3;
		Name = "Lawrencium";
		Symbol = "Lr";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 103;
		atomic_radii = 1.70;
	}
	else if (name.compare("Rutherfordium") == 0)
	{
		atomic_weight = 267.0;
		oxidation_state = 4;
		protons = 104;
		neutrons = 157;
		electrons = 104;
		valence_e = 4;
		Name = "Rutherfordium";
		Symbol = "Rf";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 104;
		atomic_radii = 1.57;
	}
	else if (name.compare("Dubnium") == 0)
	{
		atomic_weight = 268.0;
		oxidation_state = 5;
		protons = 105;
		neutrons = 157;
		electrons = 105;
		valence_e = 5;
		Name = "Dubnium";
		Symbol = "Db";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 105;
		atomic_radii = 1.49;
	}
	else if (name.compare("Seaborgium") == 0)
	{
		atomic_weight = 269.0;
		oxidation_state = 6;
		protons = 106;
		neutrons = 157;
		electrons = 106;
		valence_e = 6;
		Name = "Seaborgium";
		Symbol = "Sg";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 106;
		atomic_radii = 1.43;
	}
	else if (name.compare("Bohrium") == 0)
	{
		atomic_weight = 270.0;
		oxidation_state = 7;
		protons = 107;
		neutrons = 155;
		electrons = 107;
		valence_e = 7;
		Name = "Bohrium";
		Symbol = "Bh";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 107;
		atomic_radii = 1.41;
	}
	else if (name.compare("Hassium") == 0)
	{
		atomic_weight = 269.0;
		oxidation_state = 8;
		protons = 108;
		neutrons = 157;
		electrons = 108;
		valence_e = 8;
		Name = "Hassium";
		Symbol = "Hs";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 108;
		atomic_radii = 1.34;
	}
	else if (name.compare("Meitnerium") == 0)
	{
		atomic_weight = 278.0;
		oxidation_state = 6;
		protons = 109;
		neutrons = 157;
		electrons = 109;
		valence_e = 9;
		Name = "Meitnerium";
		Symbol = "Mt";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 109;
		atomic_radii = 1.29;
	}
	else if (name.compare("Darmstadium") == 0)
	{
		atomic_weight = 281.0;
		oxidation_state = 8;
		protons = 110;
		neutrons = 171;
		electrons = 110;
		valence_e = 10;
		Name = "Darmstadium";
		Symbol = "Ds";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 110;
		atomic_radii = 1.28;
	}
	else if (name.compare("Roentgenium") == 0)
	{
		atomic_weight = 281.0;
		oxidation_state = 3;
		protons = 111;
		neutrons = 170;
		electrons = 111;
		valence_e = 11;
		Name = "Roentgenium";
		Symbol = "Rg";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 111;
		atomic_radii = 1.21;
	}
	else if (name.compare("Copernicium") == 0)
	{
		atomic_weight = 285.0;
		oxidation_state = 4;
		protons = 112;
		neutrons = 173;
		electrons = 112;
		valence_e = 12;
		Name = "Copernicium";
		Symbol = "Cn";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 112;
		atomic_radii = 1.22;
	}
	else if (name.compare("Ununtrium") == 0)
	{
		atomic_weight = 286.0;
		oxidation_state = 1;
		protons = 113;
		neutrons = 173;
		electrons = 113;
		valence_e = 3;
		Name = "Ununtrium";
		Symbol = "Uut";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 113;
		atomic_radii = 1.76;
	}
	else if (name.compare("Flerovium") == 0)
	{
		atomic_weight = 289.0;
		oxidation_state = 2;
		protons = 114;
		neutrons = 175;
		electrons = 114;
		valence_e = 4;
		Name = "Flerovium";
		Symbol = "Fl";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 114;
		atomic_radii = 1.74;
	}
	else if (name.compare("Ununpentium") == 0)
	{
		atomic_weight = 289.0;
		oxidation_state = 1;
		protons = 115;
		neutrons = 174;
		electrons = 115;
		valence_e = 5;
		Name = "Ununpentium";
		Symbol = "Uup";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 115;
		atomic_radii = 1.57;
	}
	else if (name.compare("Livermorium") == 0)
	{
		atomic_weight = 293.0;
		oxidation_state = 2;
		protons = 116;
		neutrons = 177;
		electrons = 116;
		valence_e = 6;
		Name = "Livermorium";
		Symbol = "Lv";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 116;
		atomic_radii = 1.64;
	}
	else if (name.compare("Ununseptium") == 0)
	{
		atomic_weight = 294.0;
		oxidation_state = 1;
		protons = 117;
		neutrons = 177;
		electrons = 117;
		valence_e = 7;
		Name = "Ununseptium";
		Symbol = "Uus";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 117;
		atomic_radii = 1.57;
	}
	else if (name.compare("Ununoctium") == 0)
	{
		atomic_weight = 294.0;
		oxidation_state = 0;
		protons = 118;
		neutrons = 176;
		electrons = 118;
		valence_e = 8;
		Name = "Ununoctium";
		Symbol = "Uuo";
		Category = "Nobel Gas";
		NaturalState = "Solid";
		atomic_number = 118;
		atomic_radii = 1.57;
	}
	else
	{
		atomic_weight = 0.0;
		oxidation_state = 0;
		protons = 0;
		neutrons = 0;
		electrons = 0;
		Name = "No Name";
		Symbol = "N/A";
		Category = "N/A";
		NaturalState = "N/A";
		atomic_number = 0;
	}
}

//Register an instance of an Atom by atomic number
Atom::Atom(int number)
{
	if (number == 1)
	{
		atomic_weight = 1.0081;
		oxidation_state = 1;
		protons = 1;
		neutrons = 0;
		electrons = 1;
		valence_e = 1;
		Name = "Hydrogen";
		Symbol = "H";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 1;
		atomic_radii = 1.2;
	}
	else if (number == 2)
	{
		atomic_weight = 4.0026022;
		oxidation_state = 0;
		protons = 2;
		neutrons = 2;
		electrons = 2;
		valence_e = 2;
		Name = "Helium";
		Symbol = "He";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 2;
		atomic_radii = 1.4;
	}
	else if (number == 3)
	{
		atomic_weight = 6.941;
		oxidation_state = 1;
		protons = 3;
		neutrons = 4;
		electrons = 3;
		valence_e = 1;
		Name = "Lithium";
		Symbol = "Li";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 3;
		atomic_radii = 1.82;
	}
	else if (number == 4)
	{
		atomic_weight = 9.01218315;
		oxidation_state = 2;
		protons = 4;
		neutrons = 5;
		electrons = 4;
		valence_e = 2;
		Name = "Beryllium";
		Symbol = "Be";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 4;
		atomic_radii = 1.53;
	}
	else if (number == 5)
	{
		atomic_weight = 10.811;
		oxidation_state = 3;
		protons = 5;
		neutrons = 6;
		electrons = 5;
		valence_e = 3;
		Name = "Boron";
		Symbol = "B";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 5;
		atomic_radii = 1.92;
	}
	else if (number == 6)
	{
		atomic_weight = 12.0111;
		oxidation_state = 4;
		protons = 6;
		neutrons = 6;
		electrons = 6;
		valence_e = 4;
		Name = "Carbon";
		Symbol = "C";
		Category = "Polyatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 6;
		atomic_radii = 1.70;
	}
	else if (number == 7)
	{
		atomic_weight = 14.0071;
		oxidation_state = -3;
		protons = 7;
		neutrons = 7;
		electrons = 7;
		valence_e = 5;
		Name = "Nitrogen";
		Symbol = "N";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 7;
		atomic_radii = 1.55;
	}
	else if (number == 8)
	{
		atomic_weight = 15.9994;
		oxidation_state = -2;
		protons = 8;
		neutrons = 8;
		electrons = 8;
		valence_e = 6;
		Name = "Oxygen";
		Symbol = "O";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 8;
		atomic_radii = 1.52;
	}
	else if (number == 9)
	{
		atomic_weight = 18.9984031636;
		oxidation_state = -1;
		protons = 9;
		neutrons = 10;
		electrons = 9;
		valence_e = 7;
		Name = "Fluorine";
		Symbol = "F";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 9;
		atomic_radii = 1.35;
	}
	else if (number == 10)
	{
		atomic_weight = 20.17976;
		oxidation_state = 0;
		protons = 10;
		neutrons = 10;
		electrons = 10;
		valence_e = 8;
		Name = "Neon";
		Symbol = "Ne";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 10;
		atomic_radii = 1.54;
	}
	else if (number == 11)
	{
		atomic_weight = 22.989769282;
		oxidation_state = 1;
		protons = 11;
		neutrons = 12;
		electrons = 11;
		valence_e = 1;
		Name = "Sodium";
		Symbol = "Na";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 11;
		atomic_radii = 2.27;
	}
	else if (number == 12)
	{
		atomic_weight = 24.3051;
		oxidation_state = 2;
		protons = 12;
		neutrons = 12;
		electrons = 12;
		valence_e = 2;
		Name = "Magnesium";
		Symbol = "Mg";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 12;
		atomic_radii = 1.73;
	}
	else if (number == 13)
	{
		atomic_weight = 26.98153857;
		oxidation_state = 3;
		protons = 13;
		neutrons = 14;
		electrons = 13;
		valence_e = 3;
		Name = "Aluminium";
		Symbol = "Al";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 13;
		atomic_radii = 1.84;
	}
	else if (number == 14)
	{
		atomic_weight = 28.0851;
		oxidation_state = 4;
		protons = 14;
		neutrons = 14;
		electrons = 14;
		valence_e = 4;
		Name = "Silicon";
		Symbol = "Si";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 14;
		atomic_radii = 2.10;
	}
	else if (number == 15)
	{
		atomic_weight = 30.9737619985;
		oxidation_state = 5;
		protons = 15;
		neutrons = 16;
		electrons = 15;
		valence_e = 5;
		Name = "Phosphorus";
		Symbol = "P";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 15;
		atomic_radii = 1.80;
	}
	else if (number == 16)
	{
		atomic_weight = 32.061;
		oxidation_state = 6;
		protons = 16;
		neutrons = 16;
		electrons = 16;
		valence_e = 6;
		Name = "Sulfur";
		Symbol = "S";
		Category = "Polyatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 16;
		atomic_radii = 1.80;
	}
	else if (number == 17)
	{
		atomic_weight = 35.451;
		oxidation_state = -1;
		protons = 17;
		neutrons = 18;
		electrons = 17;
		valence_e = 7;
		Name = "Chlorine";
		Symbol = "Cl";
		Category = "Diatomic Non-metal";
		NaturalState = "Gas";
		atomic_number = 17;
		atomic_radii = 1.75;
	}
	else if (number == 18)
	{
		atomic_weight = 39.9481;
		oxidation_state = 0;
		protons = 18;
		neutrons = 22;
		electrons = 18;
		valence_e = 8;
		Name = "Argon";
		Symbol = "Ar";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 18;
		atomic_radii = 1.88;
	}
	else if (number == 19)
	{
		atomic_weight = 39.09831;
		oxidation_state = 1;
		protons = 19;
		neutrons = 20;
		electrons = 19;
		valence_e = 1;
		Name = "Potassium";
		Symbol = "K";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 19;
		atomic_radii = 2.75;
	}
	else if (number == 20)
	{
		atomic_weight = 40.0784;
		oxidation_state = 2;
		protons = 20;
		neutrons = 20;
		electrons = 20;
		valence_e = 2;
		Name = "Calcium";
		Symbol = "Ca";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 20;
		atomic_radii = 2.31;
	}
	else if (number == 21)
	{
		atomic_weight = 44.9559085;
		oxidation_state = 3;
		protons = 21;
		neutrons = 24;
		electrons = 21;
		valence_e = 3;
		Name = "Scandium";
		Symbol = "Sc";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 21;
		atomic_radii = 2.11;
	}
	else if (number == 22)
	{
		atomic_weight = 47.8671;
		oxidation_state = 4;
		protons = 22;
		neutrons = 26;
		electrons = 22;
		valence_e = 4;
		Name = "Titanium";
		Symbol = "Ti";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 22;
		atomic_radii = 1.60;
	}
	else if (number == 23)
	{
		atomic_weight = 50.94151;
		oxidation_state = 5;
		protons = 23;
		neutrons = 28;
		electrons = 23;
		valence_e = 5;
		Name = "Vanadium";
		Symbol = "V";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 23;
		atomic_radii = 1.53;
	}
	else if (number == 24)
	{
		atomic_weight = 51.99616;
		oxidation_state = 6;
		protons = 24;
		neutrons = 28;
		electrons = 24;
		valence_e = 6;
		Name = "Chromium";
		Symbol = "Cr";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 24;
		atomic_radii = 1.39;
	}
	else if (number == 25)
	{
		atomic_weight = 54.9380443;
		oxidation_state = 7;
		protons = 25;
		neutrons = 30;
		electrons = 25;
		valence_e = 7;
		Name = "Manganese";
		Symbol = "Mn";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 25;
		atomic_radii = 1.39;
	}
	else if (number == 26)
	{
		atomic_weight = 55.8452;
		oxidation_state = 3;
		protons = 26;
		neutrons = 30;
		electrons = 26;
		valence_e = 8;
		Name = "Iron";
		Symbol = "Fe";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 26;
		atomic_radii = 1.32;
	}
	else if (number == 27)
	{
		atomic_weight = 58.9331944;
		oxidation_state = 3;
		protons = 27;
		neutrons = 32;
		electrons = 27;
		valence_e = 9;
		Name = "Cobalt";
		Symbol = "Co";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 27;
		atomic_radii = 1.26;
	}
	else if (number == 28)
	{
		atomic_weight = 58.69344;
		oxidation_state = 2;
		protons = 28;
		neutrons = 31;
		electrons = 28;
		valence_e = 10;
		Name = "Nickel";
		Symbol = "Ni";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 28;
		atomic_radii = 1.63;
	}
	else if (number == 29)
	{
		atomic_weight = 63.5463;
		oxidation_state = 2;
		protons = 29;
		neutrons = 35;
		electrons = 29;
		valence_e = 11;
		Name = "Copper";
		Symbol = "Cu";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 29;
		atomic_radii = 1.40;
	}
	else if (number == 30)
	{
		atomic_weight = 65.382;
		oxidation_state = 2;
		protons = 30;
		neutrons = 35;
		electrons = 30;
		valence_e = 12;
		Name = "Zinc";
		Symbol = "Zn";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 30;
		atomic_radii = 1.39;
	}
	else if (number == 31)
	{
		atomic_weight = 69.7231;
		oxidation_state = 3;
		protons = 31;
		neutrons = 39;
		electrons = 31;
		valence_e = 3;
		Name = "Gallium";
		Symbol = "Ga";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 31;
		atomic_radii = 1.87;
	}
	else if (number == 32)
	{
		atomic_weight = 72.6308;
		oxidation_state = 4;
		protons = 32;
		neutrons = 41;
		electrons = 32;
		valence_e = 4;
		Name = "Germanium";
		Symbol = "Ge";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 32;
		atomic_radii = 2.11;
	}
	else if (number == 33)
	{
		atomic_weight = 74.9215956;
		oxidation_state = 5;
		protons = 33;
		neutrons = 42;
		electrons = 33;
		valence_e = 5;
		Name = "Arsenic";
		Symbol = "As";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 33;
		atomic_radii = 1.85;
	}
	else if (number == 34)
	{
		atomic_weight = 78.9718;
		oxidation_state = 6;
		protons = 34;
		neutrons = 45;
		electrons = 34;
		valence_e = 6;
		Name = "Selenium";
		Symbol = "Se";
		Category = "Polyatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 34;
		atomic_radii = 1.90;
	}
	else if (number == 35)
	{
		atomic_weight = 79.9041;
		oxidation_state = -1;
		protons = 35;
		neutrons = 45;
		electrons = 35;
		valence_e = 7;
		Name = "Bromine";
		Symbol = "Br";
		Category = "Diatomic Non-metal";
		NaturalState = "Liquid";
		atomic_number = 35;
		atomic_radii = 1.85;
	}
	else if (number == 36)
	{
		atomic_weight = 83.798;
		oxidation_state = 0;
		protons = 36;
		neutrons = 48;
		electrons = 36;
		valence_e = 8;
		Name = "Krypton";
		Symbol = "Kr";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 36;
		atomic_radii = 2.02;
	}
	else if (number == 37)
	{
		atomic_weight = 85.46783;
		oxidation_state = 1;
		protons = 37;
		neutrons = 48;
		electrons = 37;
		valence_e = 1;
		Name = "Rubidium";
		Symbol = "Rb";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 37;
		atomic_radii = 3.03;
	}
	else if (number == 38)
	{
		atomic_weight = 87.621;
		oxidation_state = 2;
		protons = 38;
		neutrons = 50;
		electrons = 38;
		valence_e = 2;
		Name = "Strontium";
		Symbol = "Sr";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 38;
		atomic_radii = 2.49;
	}
	else if (number == 39)
	{
		atomic_weight = 88.905842;
		oxidation_state = 3;
		protons = 39;
		neutrons = 50;
		electrons = 39;
		valence_e = 3;
		Name = "Yttrium";
		Symbol = "Y";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 39;
		atomic_radii = 1.90;
	}
	else if (number == 40)
	{
		atomic_weight = 91.2242;
		oxidation_state = 4;
		protons = 40;
		neutrons = 51;
		electrons = 40;
		valence_e = 4;
		Name = "Zirconium";
		Symbol = "Zr";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 40;
		atomic_radii = 1.75;
	}
	else if (number == 41)
	{
		atomic_weight = 92.906372;
		oxidation_state = 5;
		protons = 41;
		neutrons = 52;
		electrons = 41;
		valence_e = 5;
		Name = "Niobium";
		Symbol = "Nb";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 41;
		atomic_radii = 1.64;
	}
	else if (number == 42)
	{
		atomic_weight = 95.951;
		oxidation_state = 6;
		protons = 42;
		neutrons = 54;
		electrons = 42;
		valence_e = 6;
		Name = "Molybdenum";
		Symbol = "Mo";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 42;
		atomic_radii = 1.54;
	}
	else if (number == 43)
	{
		atomic_weight = 98.0;
		oxidation_state = 7;
		protons = 43;
		neutrons = 55;
		electrons = 43;
		valence_e = 7;
		Name = "Technetium";
		Symbol = "Tc";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 43;
		atomic_radii = 1.47;
	}
	else if (number == 44)
	{
		atomic_weight = 101.072;
		oxidation_state = 4;
		protons = 44;
		neutrons = 57;
		electrons = 44;
		valence_e = 8;
		Name = "Ruthenium";
		Symbol = "Ru";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 44;
		atomic_radii = 1.46;
	}
	else if (number == 45)
	{
		atomic_weight = 102.905502;
		oxidation_state = 3;
		protons = 45;
		neutrons = 58;
		electrons = 45;
		valence_e = 9;
		Name = "Rhodium";
		Symbol = "Rh";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 45;
		atomic_radii = 1.42;
	}
	else if (number == 46)
	{
		atomic_weight = 106.421;
		oxidation_state = 4;
		protons = 46;
		neutrons = 60;
		electrons = 46;
		valence_e = 10;
		Name = "Palladium";
		Symbol = "Pd";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 46;
		atomic_radii = 1.63;
	}
	else if (number == 47)
	{
		atomic_weight = 107.86822;
		oxidation_state = 1;
		protons = 47;
		neutrons = 61;
		electrons = 47;
		valence_e = 11;
		Name = "Silver";
		Symbol = "Ag";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 47;
		atomic_radii = 1.72;
	}
	else if (number == 48)
	{
		atomic_weight = 112.4144;
		oxidation_state = 2;
		protons = 48;
		neutrons = 64;
		electrons = 48;
		valence_e = 12;
		Name = "Cadmium";
		Symbol = "Cd";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 48;
		atomic_radii = 1.58;
	}
	else if (number == 49)
	{
		atomic_weight = 114.8181;
		oxidation_state = 3;
		protons = 49;
		neutrons = 66;
		electrons = 49;
		valence_e = 3;
		Name = "Indium";
		Symbol = "In";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 49;
		atomic_radii = 1.93;
	}
	else if (number == 50)
	{
		atomic_weight = 118.7107;
		oxidation_state = 4;
		protons = 50;
		neutrons = 69;
		electrons = 50;
		valence_e = 4;
		Name = "Tin";
		Symbol = "Sn";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 50;
		atomic_radii = 2.17;
	}
	else if (number == 51)
	{
		atomic_weight = 121.7601;
		oxidation_state = 5;
		protons = 51;
		neutrons = 71;
		electrons = 51;
		valence_e = 5;
		Name = "Antimony";
		Symbol = "Sb";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 51;
		atomic_radii = 2.06;
	}
	else if (number == 52)
	{
		atomic_weight = 127.603;
		oxidation_state = 6;
		protons = 52;
		neutrons = 76;
		electrons = 52;
		valence_e = 6;
		Name = "Tellurium";
		Symbol = "Te";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 52;
		atomic_radii = 2.06;
	}
	else if (number == 53)
	{
		atomic_weight = 126.904473;
		oxidation_state = -1;
		protons = 53;
		neutrons = 74;
		electrons = 53;
		valence_e = 7;
		Name = "Iodine";
		Symbol = "I";
		Category = "Diatomic Non-metal";
		NaturalState = "Solid";
		atomic_number = 53;
		atomic_radii = 1.98;
	}
	else if (number == 54)
	{
		atomic_weight = 131.2936;
		oxidation_state = 0;
		protons = 54;
		neutrons = 77;
		electrons = 54;
		valence_e = 8;
		Name = "Xenon";
		Symbol = "Xe";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 54;
		atomic_radii = 2.16;
	}
	else if (number == 55)
	{
		atomic_weight = 132.905451966;
		oxidation_state = 1;
		protons = 55;
		neutrons = 78;
		electrons = 55;
		valence_e = 1;
		Name = "Caesium";
		Symbol = "Cs";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 55;
		atomic_radii = 3.43;
	}
	else if (number == 56)
	{
		atomic_weight = 137.3277;
		oxidation_state = 2;
		protons = 56;
		neutrons = 81;
		electrons = 56;
		valence_e = 2;
		Name = "Barium";
		Symbol = "Ba";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 56;
		atomic_radii = 2.68;
	}
	else if (number == 57)
	{
		atomic_weight = 138.90547;
		oxidation_state = 3;
		protons = 57;
		neutrons = 82;
		electrons = 57;
		valence_e = 3;
		Name = "Lanthanum";
		Symbol = "La";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 57;
		atomic_radii = 2.07;
	}
	else if (number == 58)
	{
		atomic_weight = 140.116;
		oxidation_state = 4;
		protons = 58;
		neutrons = 82;
		electrons = 58;
		valence_e = 4;
		Name = "Cerium";
		Symbol = "Ce";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 58;
		atomic_radii = 2.04;
	}
	else if (number == 59)
	{
		atomic_weight = 140.907662;
		oxidation_state = 4;
		protons = 59;
		neutrons = 82;
		electrons = 59;
		valence_e = 5;
		Name = "Praseodymium";
		Symbol = "Pr";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 59;
		atomic_radii = 2.03;
	}
	else if (number == 60)
	{
		atomic_weight = 144.242;
		oxidation_state = 3;
		protons = 60;
		neutrons = 84;
		electrons = 60;
		valence_e = 6;
		Name = "Neodymium";
		Symbol = "Nd";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 60;
		atomic_radii = 2.01;
	}
	else if (number == 61)
	{
		atomic_weight = 145.0;
		oxidation_state = 3;
		protons = 61;
		neutrons = 84;
		electrons = 61;
		valence_e = 7;
		Name = "Promethium";
		Symbol = "Pm";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 61;
		atomic_radii = 1.99;
	}
	else if (number == 62)
	{
		atomic_weight = 150.362;
		oxidation_state = 3;
		protons = 62;
		neutrons = 88;
		electrons = 62;
		valence_e = 8;
		Name = "Samarium";
		Symbol = "Sm";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 62;
		atomic_radii = 1.98;
	}
	else if (number == 63)
	{
		atomic_weight = 151.964;
		oxidation_state = 3;
		protons = 63;
		neutrons = 89;
		electrons = 63;
		valence_e = 9;
		Name = "Europium";
		Symbol = "Eu";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 63;
		atomic_radii = 1.98;
	}
	else if (number == 64)
	{
		atomic_weight = 157.253;
		oxidation_state = 3;
		protons = 64;
		neutrons = 93;
		electrons = 64;
		valence_e = 10;
		Name = "Gadolinium";
		Symbol = "Gd";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 64;
		atomic_radii = 1.96;
	}
	else if (number == 65)
	{
		atomic_weight = 158.92535;
		oxidation_state = 3;
		protons = 65;
		neutrons = 94;
		electrons = 65;
		valence_e = 11;
		Name = "Terbium";
		Symbol = "Tb";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 65;
		atomic_radii = 1.94;
	}
	else if (number == 66)
	{
		atomic_weight = 162.5001;
		oxidation_state = 3;
		protons = 66;
		neutrons = 97;
		electrons = 66;
		valence_e = 12;
		Name = "Dysprosium";
		Symbol = "Dy";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 66;
		atomic_radii = 1.92;
	}
	else if (number == 67)
	{
		atomic_weight = 164.930332;
		oxidation_state = 3;
		protons = 67;
		neutrons = 98;
		electrons = 67;
		valence_e = 13;
		Name = "Holmium";
		Symbol = "Ho";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 67;
		atomic_radii = 1.92;
	}
	else if (number == 68)
	{
		atomic_weight = 167.259;
		oxidation_state = 3;
		protons = 68;
		neutrons = 99;
		electrons = 68;
		valence_e = 14;
		Name = "Erbium";
		Symbol = "Er";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 68;
		atomic_radii = 1.89;
	}
	else if (number == 69)
	{
		atomic_weight = 168.934222;
		oxidation_state = 3;
		protons = 69;
		neutrons = 100;
		electrons = 69;
		valence_e = 15;
		Name = "Thulium";
		Symbol = "Tm";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 69;
		atomic_radii = 1.90;
	}
	else if (number == 70)
	{
		atomic_weight = 173.0545;
		oxidation_state = 3;
		protons = 70;
		neutrons = 103;
		electrons = 70;
		valence_e = 16;
		Name = "Ytterbium";
		Symbol = "Yb";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 70;
		atomic_radii = 1.87;
	}
	else if (number == 71)
	{
		atomic_weight = 174.96684;
		oxidation_state = 3;
		protons = 71;
		neutrons = 104;
		electrons = 71;
		valence_e = 3;
		Name = "Lutetium";
		Symbol = "Lu";
		Category = "Lanthanide";
		NaturalState = "Solid";
		atomic_number = 71;
		atomic_radii = 1.87;
	}
	else if (number == 72)
	{
		atomic_weight = 178.492;
		oxidation_state = 4;
		protons = 72;
		neutrons = 106;
		electrons = 72;
		valence_e = 4;
		Name = "Hafnium";
		Symbol = "Hf";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 72;
		atomic_radii = 1.75;
	}
	else if (number == 73)
	{
		atomic_weight = 180.947882;
		oxidation_state = 5;
		protons = 73;
		neutrons = 108;
		electrons = 73;
		valence_e = 5;
		Name = "Tantalum";
		Symbol = "Ta";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 73;
		atomic_radii = 1.70;
	}
	else if (number == 74)
	{
		atomic_weight = 183.841;
		oxidation_state = 6;
		protons = 74;
		neutrons = 110;
		electrons = 74;
		valence_e = 6;
		Name = "Tungsten";
		Symbol = "W";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 74;
		atomic_radii = 1.62;
	}
	else if (number == 75)
	{
		atomic_weight = 186.2071;
		oxidation_state = 7;
		protons = 75;
		neutrons = 111;
		electrons = 75;
		valence_e = 7;
		Name = "Rhenium";
		Symbol = "Re";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 75;
		atomic_radii = 1.51;
	}
	else if (number == 76)
	{
		atomic_weight = 190.233;
		oxidation_state = 4;
		protons = 76;
		neutrons = 114;
		electrons = 76;
		valence_e = 8;
		Name = "Osmium";
		Symbol = "Os";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 76;
		atomic_radii = 1.44;
	}
	else if (number == 77)
	{
		atomic_weight = 192.2173;
		oxidation_state = 4;
		protons = 77;
		neutrons = 115;
		electrons = 77;
		valence_e = 9;
		Name = "Iridium";
		Symbol = "Ir";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 77;
		atomic_radii = 1.41;
	}
	else if (number == 78)
	{
		atomic_weight = 195.0849;
		oxidation_state = 4;
		protons = 78;
		neutrons = 117;
		electrons = 78;
		valence_e = 10;
		Name = "Platinum";
		Symbol = "Pt";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 78;
		atomic_radii = 1.75;
	}
	else if (number == 79)
	{
		atomic_weight = 196.9665694;
		oxidation_state = 3;
		protons = 79;
		neutrons = 118;
		electrons = 79;
		valence_e = 11;
		Name = "Gold";
		Symbol = "Au";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 79;
		atomic_radii = 1.66;
	}
	else if (number == 80)
	{
		atomic_weight = 200.5923;
		oxidation_state = 2;
		protons = 80;
		neutrons = 121;
		electrons = 80;
		valence_e = 12;
		Name = "Mercury";
		Symbol = "Hg";
		Category = "Transition Metal";
		NaturalState = "Liquid";
		atomic_number = 80;
		atomic_radii = 1.55;
	}
	else if (number == 81)
	{
		atomic_weight = 204.381;
		oxidation_state = 1;
		protons = 81;
		neutrons = 123;
		electrons = 81;
		valence_e = 3;
		Name = "Thallium";
		Symbol = "Tl";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 81;
		atomic_radii = 1.96;
	}
	else if (number == 82)
	{
		atomic_weight = 207.21;
		oxidation_state = 2;
		protons = 82;
		neutrons = 125;
		electrons = 82;
		valence_e = 4;
		Name = "Lead";
		Symbol = "Pb";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 82;
		atomic_radii = 2.02;
	}
	else if (number == 83)
	{
		atomic_weight = 208.980401;
		oxidation_state = 3;
		protons = 83;
		neutrons = 126;
		electrons = 83;
		valence_e = 5;
		Name = "Bismuth";
		Symbol = "Bi";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 83;
		atomic_radii = 2.07;
	}
	else if (number == 84)
	{
		atomic_weight = 209.0;
		oxidation_state = 4;
		protons = 84;
		neutrons = 125;
		electrons = 84;
		valence_e = 6;
		Name = "Polonium";
		Symbol = "Po";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 84;
		atomic_radii = 1.97;
	}
	else if (number == 85)
	{
		atomic_weight = 210.0;
		oxidation_state = -1;
		protons = 85;
		neutrons = 125;
		electrons = 85;
		valence_e = 7;
		Name = "Astatine";
		Symbol = "At";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 85;
		atomic_radii = 2.02;
	}
	else if (number == 86)
	{
		atomic_weight = 222.0;
		oxidation_state = 0;
		protons = 86;
		neutrons = 136;
		electrons = 86;
		valence_e = 8;
		Name = "Radon";
		Symbol = "Rn";
		Category = "Nobel Gas";
		NaturalState = "Gas";
		atomic_number = 86;
		atomic_radii = 2.20;
	}
	else if (number == 87)
	{
		atomic_weight = 223.0;
		oxidation_state = 1;
		protons = 87;
		neutrons = 136;
		electrons = 87;
		valence_e = 1;
		Name = "Francium";
		Symbol = "Fr";
		Category = "Alkali Metal";
		NaturalState = "Solid";
		atomic_number = 87;
		atomic_radii = 3.48;
	}
	else if (number == 88)
	{
		atomic_weight = 226.0;
		oxidation_state = 2;
		protons = 88;
		neutrons = 138;
		electrons = 88;
		valence_e = 2;
		Name = "Radium";
		Symbol = "Ra";
		Category = "Alkaline Earth Metal";
		NaturalState = "Solid";
		atomic_number = 88;
		atomic_radii = 2.83;
	}
	else if (number == 89)
	{
		atomic_weight = 227.0;
		oxidation_state = 3;
		protons = 89;
		neutrons = 138;
		electrons = 89;
		valence_e = 3;
		Name = "Actinium";
		Symbol = "Ac";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 89;
		atomic_radii = 2.15;
	}
	else if (number == 90)
	{
		atomic_weight = 232.03774;
		oxidation_state = 4;
		protons = 90;
		neutrons = 142;
		electrons = 90;
		valence_e = 4;
		Name = "Thorium";
		Symbol = "Th";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 90;
		atomic_radii = 2.06;
	}
	else if (number == 91)
	{
		atomic_weight = 231.03588;
		oxidation_state = 5;
		protons = 91;
		neutrons = 140;
		electrons = 91;
		valence_e = 5;
		Name = "Protactinium";
		Symbol = "Pa";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 91;
		atomic_radii = 2.00;
	}
	else if (number == 92)
	{
		atomic_weight = 238.028913;
		oxidation_state = 6;
		protons = 92;
		neutrons = 146;
		electrons = 92;
		valence_e = 6;
		Name = "Uranium";
		Symbol = "U";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 92;
		atomic_radii = 1.86;
	}
	else if (number == 93)
	{
		atomic_weight = 237.0;
		oxidation_state = 5;
		protons = 93;
		neutrons = 144;
		electrons = 93;
		valence_e = 7;
		Name = "Neptunium";
		Symbol = "Np";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 93;
		atomic_radii = 1.90;
	}
	else if (number == 94)
	{
		atomic_weight = 244.0;
		oxidation_state = 4;
		protons = 94;
		neutrons = 150;
		electrons = 94;
		valence_e = 8;
		Name = "Plutonium";
		Symbol = "Pu";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 94;
		atomic_radii = 1.87;
	}
	else if (number == 95)
	{
		atomic_weight = 243.0;
		oxidation_state = 3;
		protons = 95;
		neutrons = 148;
		electrons = 95;
		valence_e = 9;
		Name = "Americium";
		Symbol = "Am";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 95;
		atomic_radii = 1.80;
	}
	else if (number == 96)
	{
		atomic_weight = 247.0;
		oxidation_state = 3;
		protons = 96;
		neutrons = 151;
		electrons = 96;
		valence_e = 10;
		Name = "Curium";
		Symbol = "Cm";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 96;
		atomic_radii = 1.69;
	}
	else if (number == 97)
	{
		atomic_weight = 247.0;
		oxidation_state = 3;
		protons = 97;
		neutrons = 150;
		electrons = 97;
		valence_e = 11;
		Name = "Berkelium";
		Symbol = "Bk";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 97;
		atomic_radii = 1.70;
	}
	else if (number == 98)
	{
		atomic_weight = 251.0;
		oxidation_state = 3;
		protons = 98;
		neutrons = 153;
		electrons = 98;
		valence_e = 12;
		Name = "Californium";
		Symbol = "Cf";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 98;
		atomic_radii = 1.70;
	}
	else if (number == 99)
	{
		atomic_weight = 252.0;
		oxidation_state = 3;
		protons = 99;
		neutrons = 153;
		electrons = 99;
		valence_e = 13;
		Name = "Einsteinium";
		Symbol = "Es";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 99;
		atomic_radii = 1.70;
	}
	else if (number == 100)
	{
		atomic_weight = 257.0;
		oxidation_state = 3;
		protons = 100;
		neutrons = 157;
		electrons = 100;
		valence_e = 14;
		Name = "Fermium";
		Symbol = "Fm";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 100;
		atomic_radii = 1.70;
	}
	else if (number == 101)
	{
		atomic_weight = 258.0;
		oxidation_state = 3;
		protons = 101;
		neutrons = 157;
		electrons = 101;
		valence_e = 15;
		Name = "Mendelevium";
		Symbol = "Md";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 101;
		atomic_radii = 1.70;
	}
	else if (number == 102)
	{
		atomic_weight = 259.0;
		oxidation_state = 2;
		protons = 102;
		neutrons = 157;
		electrons = 102;
		valence_e = 16;
		Name = "Nobelium";
		Symbol = "No";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 102;
		atomic_radii = 1.70;
	}
	else if (number == 103)
	{
		atomic_weight = 266.0;
		oxidation_state = 3;
		protons = 103;
		neutrons = 159;
		electrons = 103;
		valence_e = 3;
		Name = "Lawrencium";
		Symbol = "Lr";
		Category = "Actinide";
		NaturalState = "Solid";
		atomic_number = 103;
		atomic_radii = 1.70;
	}
	else if (number == 104)
	{
		atomic_weight = 267.0;
		oxidation_state = 4;
		protons = 104;
		neutrons = 157;
		electrons = 104;
		valence_e = 4;
		Name = "Rutherfordium";
		Symbol = "Rf";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 104;
		atomic_radii = 1.57;
	}
	else if (number == 105)
	{
		atomic_weight = 268.0;
		oxidation_state = 5;
		protons = 105;
		neutrons = 157;
		electrons = 105;
		valence_e = 5;
		Name = "Dubnium";
		Symbol = "Db";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 105;
		atomic_radii = 1.49;
	}
	else if (number == 106)
	{
		atomic_weight = 269.0;
		oxidation_state = 6;
		protons = 106;
		neutrons = 157;
		electrons = 106;
		valence_e = 6;
		Name = "Seaborgium";
		Symbol = "Sg";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 106;
		atomic_radii = 1.43;
	}
	else if (number == 107)
	{
		atomic_weight = 270.0;
		oxidation_state = 7;
		protons = 107;
		neutrons = 155;
		electrons = 107;
		valence_e = 7;
		Name = "Bohrium";
		Symbol = "Bh";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 107;
		atomic_radii = 1.41;
	}
	else if (number == 108)
	{
		atomic_weight = 269.0;
		oxidation_state = 8;
		protons = 108;
		neutrons = 157;
		electrons = 108;
		valence_e = 8;
		Name = "Hassium";
		Symbol = "Hs";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 108;
		atomic_radii = 1.34;
	}
	else if (number == 109)
	{
		atomic_weight = 278.0;
		oxidation_state = 6;
		protons = 109;
		neutrons = 157;
		electrons = 109;
		valence_e = 9;
		Name = "Meitnerium";
		Symbol = "Mt";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 109;
		atomic_radii = 1.29;
	}
	else if (number == 110)
	{
		atomic_weight = 281.0;
		oxidation_state = 8;
		protons = 110;
		neutrons = 171;
		electrons = 110;
		valence_e = 10;
		Name = "Darmstadium";
		Symbol = "Ds";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 110;
		atomic_radii = 1.28;
	}
	else if (number == 111)
	{
		atomic_weight = 281.0;
		oxidation_state = 3;
		protons = 111;
		neutrons = 170;
		electrons = 111;
		valence_e = 11;
		Name = "Roentgenium";
		Symbol = "Rg";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 111;
		atomic_radii = 1.21;
	}
	else if (number == 112)
	{
		atomic_weight = 285.0;
		oxidation_state = 4;
		protons = 112;
		neutrons = 173;
		electrons = 112;
		valence_e = 12;
		Name = "Copernicium";
		Symbol = "Cn";
		Category = "Transition Metal";
		NaturalState = "Solid";
		atomic_number = 112;
		atomic_radii = 1.22;
	}
	else if (number == 113)
	{
		atomic_weight = 286.0;
		oxidation_state = 1;
		protons = 113;
		neutrons = 173;
		electrons = 113;
		valence_e = 3;
		Name = "Ununtrium";
		Symbol = "Uut";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 113;
		atomic_radii = 1.76;
	}
	else if (number == 114)
	{
		atomic_weight = 289.0;
		oxidation_state = 2;
		protons = 114;
		neutrons = 175;
		electrons = 114;
		valence_e = 4;
		Name = "Flerovium";
		Symbol = "Fl";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 114;
		atomic_radii = 1.74;
	}
	else if (number == 115)
	{
		atomic_weight = 289.0;
		oxidation_state = 1;
		protons = 115;
		neutrons = 174;
		electrons = 115;
		valence_e = 5;
		Name = "Ununpentium";
		Symbol = "Uup";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 115;
		atomic_radii = 1.57;
	}
	else if (number == 116)
	{
		atomic_weight = 293.0;
		oxidation_state = 2;
		protons = 116;
		neutrons = 177;
		electrons = 116;
		valence_e = 6;
		Name = "Livermorium";
		Symbol = "Lv";
		Category = "Post-Transition Metal";
		NaturalState = "Solid";
		atomic_number = 116;
		atomic_radii = 1.64;
	}
	else if (number == 117)
	{
		atomic_weight = 294.0;
		oxidation_state = 1;
		protons = 117;
		neutrons = 177;
		electrons = 117;
		valence_e = 7;
		Name = "Ununseptium";
		Symbol = "Uus";
		Category = "Metalloid";
		NaturalState = "Solid";
		atomic_number = 117;
		atomic_radii = 1.57;
	}
	else if (number == 118)
	{
		atomic_weight = 294.0;
		oxidation_state = 0;
		protons = 118;
		neutrons = 176;
		electrons = 118;
		valence_e = 7;
		Name = "Ununoctium";
		Symbol = "Uuo";
		Category = "Nobel Gas";
		NaturalState = "Solid";
		atomic_number = 118;
		atomic_radii = 1.57;
	}
	else
	{
		atomic_weight = 0.0;
		oxidation_state = 0;
		protons = 0;
		neutrons = 0;
		electrons = 0;
		Name = "No Name";
		Symbol = "N/A";
		Category = "N/A";
		NaturalState = "N/A";
		atomic_number = 0;
	}
}

//Register an instance of an Atom object by it's symbol
void Atom::Register(std::string Symbol)
{
	if (Symbol.compare("H") == 0)
	{
		this->atomic_weight = 1.0081;
		this->oxidation_state = 1;
		this->protons = 1;
		this->neutrons = 0;
		this->electrons = 1;
		this->valence_e = 1;
		this->Name = "Hydrogen";
		this->Symbol= "H";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 1;
		this->atomic_radii = 1.2;
	}
	else if (Symbol.compare("He") == 0)
	{
		this->atomic_weight = 4.0026022;
		this->oxidation_state = 0;
		this->protons = 2;
		this->neutrons = 2;
		this->electrons = 2;
		this->valence_e = 2;
		this->Name = "Helium";
		this->Symbol= "He";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 2;
		this->atomic_radii = 1.4;
	}
	else if (Symbol.compare("Li") == 0)
	{
		this->atomic_weight = 6.941;
		this->oxidation_state = 1;
		this->protons = 3;
		this->neutrons = 4;
		this->electrons = 3;
		this->valence_e = 1;
		this->Name = "Lithium";
		this->Symbol= "Li";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 3;
		this->atomic_radii = 1.82;
	}
	else if (Symbol.compare("Be") == 0)
	{
		this->atomic_weight = 9.01218315;
		this->oxidation_state = 2;
		this->protons = 4;
		this->neutrons = 5;
		this->electrons = 4;
		this->valence_e = 2;
		this->Name = "Beryllium";
		this->Symbol= "Be";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 4;
		this->atomic_radii = 1.53;
	}
	else if (Symbol.compare("B") == 0)
	{
		this->atomic_weight = 10.811;
		this->oxidation_state = 3;
		this->protons = 5;
		this->neutrons = 6;
		this->electrons = 5;
		this->valence_e = 3;
		this->Name = "Boron";
		this->Symbol= "B";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 5;
		this->atomic_radii = 1.92;
	}
	else if (Symbol.compare("C") == 0)
	{
		this->atomic_weight = 12.0111;
		this->oxidation_state = 4;
		this->protons = 6;
		this->neutrons = 6;
		this->electrons = 6;
		this->valence_e = 4;
		this->Name = "Carbon";
		this->Symbol= "C";
		this->Category = "Polyatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 6;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("N") == 0)
	{
		this->atomic_weight = 14.0071;
		this->oxidation_state = -3;
		this->protons = 7;
		this->neutrons = 7;
		this->electrons = 7;
		this->valence_e = 5;
		this->Name = "Nitrogen";
		this->Symbol= "N";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 7;
		this->atomic_radii = 1.55;
	}
	else if (Symbol.compare("O") == 0)
	{
		this->atomic_weight = 15.9994;
		this->oxidation_state = -2;
		this->protons = 8;
		this->neutrons = 8;
		this->electrons = 8;
		this->valence_e = 6;
		this->Name = "Oxygen";
		this->Symbol= "O";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 8;
		this->atomic_radii = 1.52;
	}
	else if (Symbol.compare("F") == 0)
	{
		this->atomic_weight = 18.9984031636;
		this->oxidation_state = -1;
		this->protons = 9;
		this->neutrons = 10;
		this->electrons = 9;
		this->valence_e = 7;
		this->Name = "Fluorine";
		this->Symbol= "F";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 9;
		this->atomic_radii = 1.35;
	}
	else if (Symbol.compare("Ne") == 0)
	{
		this->atomic_weight = 20.17976;
		this->oxidation_state = 0;
		this->protons = 10;
		this->neutrons = 10;
		this->electrons = 10;
		this->valence_e = 8;
		this->Name = "Neon";
		this->Symbol= "Ne";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 10;
		this->atomic_radii = 1.54;
	}
	else if (Symbol.compare("Na") == 0)
	{
		this->atomic_weight = 22.989769282;
		this->oxidation_state = 1;
		this->protons = 11;
		this->neutrons = 12;
		this->electrons = 11;
		this->valence_e = 1;
		this->Name = "Sodium";
		this->Symbol= "Na";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 11;
		this->atomic_radii = 2.27;
	}
	else if (Symbol.compare("Mg") == 0)
	{
		this->atomic_weight = 24.3051;
		this->oxidation_state = 2;
		this->protons = 12;
		this->neutrons = 12;
		this->electrons = 12;
		this->valence_e = 2;
		this->Name = "Magnesium";
		this->Symbol= "Mg";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 12;
		this->atomic_radii = 1.73;
	}
	else if (Symbol.compare("Al") == 0)
	{
		this->atomic_weight = 26.98153857;
		this->oxidation_state = 3;
		this->protons = 13;
		this->neutrons = 14;
		this->electrons = 13;
		this->valence_e = 3;
		this->Name = "Aluminium";
		this->Symbol= "Al";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 13;
		this->atomic_radii = 1.84;
	}
	else if (Symbol.compare("Si") == 0)
	{
		this->atomic_weight = 28.0851;
		this->oxidation_state = 4;
		this->protons = 14;
		this->neutrons = 14;
		this->electrons = 14;
		this->valence_e = 4;
		this->Name = "Silicon";
		this->Symbol= "Si";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 14;
		this->atomic_radii = 2.10;
	}
	else if (Symbol.compare("P") == 0)
	{
		this->atomic_weight = 30.9737619985;
		this->oxidation_state = 5;
		this->protons = 15;
		this->neutrons = 16;
		this->electrons = 15;
		this->valence_e = 5;
		this->Name = "Phosphorus";
		this->Symbol= "P";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 15;
		this->atomic_radii = 1.80;
	}
	else if (Symbol.compare("S") == 0)
	{
		this->atomic_weight = 32.061;
		this->oxidation_state = 6;
		this->protons = 16;
		this->neutrons = 16;
		this->electrons = 16;
		this->valence_e = 6;
		this->Name = "Sulfur";
		this->Symbol= "S";
		this->Category = "Polyatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 16;
		this->atomic_radii = 1.80;
	}
	else if (Symbol.compare("Cl") == 0)
	{
		this->atomic_weight = 35.451;
		this->oxidation_state = -1;
		this->protons = 17;
		this->neutrons = 18;
		this->electrons = 17;
		this->valence_e = 7;
		this->Name = "Chlorine";
		this->Symbol= "Cl";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 17;
		this->atomic_radii = 1.75;
	}
	else if (Symbol.compare("Ar") == 0)
	{
		this->atomic_weight = 39.9481;
		this->oxidation_state = 0;
		this->protons = 18;
		this->neutrons = 22;
		this->electrons = 18;
		this->valence_e = 8;
		this->Name = "Argon";
		this->Symbol= "Ar";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 18;
		this->atomic_radii = 1.88;
	}
	else if (Symbol.compare("K") == 0)
	{
		this->atomic_weight = 39.09831;
		this->oxidation_state = 1;
		this->protons = 19;
		this->neutrons = 20;
		this->electrons = 19;
		this->valence_e = 1;
		this->Name = "Potassium";
		this->Symbol= "K";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 19;
		this->atomic_radii = 2.75;
	}
	else if (Symbol.compare("Ca") == 0)
	{
		this->atomic_weight = 40.0784;
		this->oxidation_state = 2;
		this->protons = 20;
		this->neutrons = 20;
		this->electrons = 20;
		this->valence_e = 2;
		this->Name = "Calcium";
		this->Symbol= "Ca";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 20;
		this->atomic_radii = 2.31;
	}
	else if (Symbol.compare("Sc") == 0)
	{
		this->atomic_weight = 44.9559085;
		this->oxidation_state = 3;
		this->protons = 21;
		this->neutrons = 24;
		this->electrons = 21;
		this->valence_e = 3;
		this->Name = "Scandium";
		this->Symbol= "Sc";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 21;
		this->atomic_radii = 2.11;
	}
	else if (Symbol.compare("Ti") == 0)
	{
		this->atomic_weight = 47.8671;
		this->oxidation_state = 4;
		this->protons = 22;
		this->neutrons = 26;
		this->electrons = 22;
		this->valence_e = 4;
		this->Name = "Titanium";
		this->Symbol= "Ti";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 22;
		this->atomic_radii = 1.60;
	}
	else if (Symbol.compare("V") == 0)
	{
		this->atomic_weight = 50.94151;
		this->oxidation_state = 5;
		this->protons = 23;
		this->neutrons = 28;
		this->electrons = 23;
		this->valence_e = 5;
		this->Name = "Vanadium";
		this->Symbol= "V";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 23;
		this->atomic_radii = 1.53;
	}
	else if (Symbol.compare("Cr") == 0)
	{
		this->atomic_weight = 51.99616;
		this->oxidation_state = 6;
		this->protons = 24;
		this->neutrons = 28;
		this->electrons = 24;
		this->valence_e = 6;
		this->Name = "Chromium";
		this->Symbol= "Cr";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 24;
		this->atomic_radii = 1.39;
	}
	else if (Symbol.compare("Mn") == 0)
	{
		this->atomic_weight = 54.9380443;
		this->oxidation_state = 7;
		this->protons = 25;
		this->neutrons = 30;
		this->electrons = 25;
		this->valence_e = 7;
		this->Name = "Manganese";
		this->Symbol= "Mn";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 25;
		this->atomic_radii = 1.39;
	}
	else if (Symbol.compare("Fe") == 0)
	{
		this->atomic_weight = 55.8452;
		this->oxidation_state = 3;
		this->protons = 26;
		this->neutrons = 30;
		this->electrons = 26;
		this->valence_e = 8;
		this->Name = "Iron";
		this->Symbol= "Fe";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 26;
		this->atomic_radii = 1.32;
	}
	else if (Symbol.compare("Co") == 0)
	{
		this->atomic_weight = 58.9331944;
		this->oxidation_state = 3;
		this->protons = 27;
		this->neutrons = 32;
		this->electrons = 27;
		this->valence_e = 9;
		this->Name = "Cobalt";
		this->Symbol= "Co";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 27;
		this->atomic_radii = 1.26;
	}
	else if (Symbol.compare("Ni") == 0)
	{
		this->atomic_weight = 58.69344;
		this->oxidation_state = 2;
		this->protons = 28;
		this->neutrons = 31;
		this->electrons = 28;
		this->valence_e = 10;
		this->Name = "Nickel";
		this->Symbol= "Ni";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 28;
		this->atomic_radii = 1.63;
	}
	else if (Symbol.compare("Cu") == 0)
	{
		this->atomic_weight = 63.5463;
		this->oxidation_state = 2;
		this->protons = 29;
		this->neutrons = 35;
		this->electrons = 29;
		this->valence_e = 11;
		this->Name = "Copper";
		this->Symbol= "Cu";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 29;
		this->atomic_radii = 1.40;
	}
	else if (Symbol.compare("Zn") == 0)
	{
		this->atomic_weight = 65.382;
		this->oxidation_state = 2;
		this->protons = 30;
		this->neutrons = 35;
		this->electrons = 30;
		this->valence_e = 12;
		this->Name = "Zinc";
		this->Symbol= "Zn";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 30;
		this->atomic_radii = 1.39;
	}
	else if (Symbol.compare("Ga") == 0)
	{
		this->atomic_weight = 69.7231;
		this->oxidation_state = 3;
		this->protons = 31;
		this->neutrons = 39;
		this->electrons = 31;
		this->valence_e = 3;
		this->Name = "Gallium";
		this->Symbol= "Ga";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 31;
		this->atomic_radii = 1.87;
	}
	else if (Symbol.compare("Ge") == 0)
	{
		this->atomic_weight = 72.6308;
		this->oxidation_state = 4;
		this->protons = 32;
		this->neutrons = 41;
		this->electrons = 32;
		this->valence_e = 4;
		this->Name = "Germanium";
		this->Symbol= "Ge";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 32;
		this->atomic_radii = 2.11;
	}
	else if (Symbol.compare("As") == 0)
	{
		this->atomic_weight = 74.9215956;
		this->oxidation_state = 5;
		this->protons = 33;
		this->neutrons = 42;
		this->electrons = 33;
		this->valence_e = 5;
		this->Name = "Arsenic";
		this->Symbol= "As";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 33;
		this->atomic_radii = 1.85;
	}
	else if (Symbol.compare("Se") == 0)
	{
		this->atomic_weight = 78.9718;
		this->oxidation_state = 6;
		this->protons = 34;
		this->neutrons = 45;
		this->electrons = 34;
		this->valence_e = 6;
		this->Name = "Selenium";
		this->Symbol= "Se";
		this->Category = "Polyatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 34;
		this->atomic_radii = 1.90;
	}
	else if (Symbol.compare("Br") == 0)
	{
		this->atomic_weight = 79.9041;
		this->oxidation_state = -1;
		this->protons = 35;
		this->neutrons = 45;
		this->electrons = 35;
		this->valence_e = 7;
		this->Name = "Bromine";
		this->Symbol= "Br";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Liquid";
		this->atomic_number = 35;
		this->atomic_radii = 1.85;
	}
	else if (Symbol.compare("Kr") == 0)
	{
		this->atomic_weight = 83.798;
		this->oxidation_state = 0;
		this->protons = 36;
		this->neutrons = 48;
		this->electrons = 36;
		this->valence_e = 8;
		this->Name = "Krypton";
		this->Symbol= "Kr";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 36;
		this->atomic_radii = 2.02;
	}
	else if (Symbol.compare("Rb") == 0)
	{
		this->atomic_weight = 85.46783;
		this->oxidation_state = 1;
		this->protons = 37;
		this->neutrons = 48;
		this->electrons = 37;
		this->valence_e = 1;
		this->Name = "Rubidium";
		this->Symbol= "Rb";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 37;
		this->atomic_radii = 3.03;
	}
	else if (Symbol.compare("Sr") == 0)
	{
		this->atomic_weight = 87.621;
		this->oxidation_state = 2;
		this->protons = 38;
		this->neutrons = 50;
		this->electrons = 38;
		this->valence_e = 2;
		this->Name = "Strontium";
		this->Symbol= "Sr";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 38;
		this->atomic_radii = 2.49;
	}
	else if (Symbol.compare("Y") == 0)
	{
		this->atomic_weight = 88.905842;
		this->oxidation_state = 3;
		this->protons = 39;
		this->neutrons = 50;
		this->electrons = 39;
		this->valence_e = 3;
		this->Name = "Yttrium";
		this->Symbol= "Y";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 39;
		this->atomic_radii = 1.90;
	}
	else if (Symbol.compare("Zr") == 0)
	{
		this->atomic_weight = 91.2242;
		this->oxidation_state = 4;
		this->protons = 40;
		this->neutrons = 51;
		this->electrons = 40;
		this->valence_e = 4;
		this->Name = "Zirconium";
		this->Symbol= "Zr";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 40;
		this->atomic_radii = 1.75;
	}
	else if (Symbol.compare("Nb") == 0)
	{
		this->atomic_weight = 92.906372;
		this->oxidation_state = 5;
		this->protons = 41;
		this->neutrons = 52;
		this->electrons = 41;
		this->valence_e = 5;
		this->Name = "Niobium";
		this->Symbol= "Nb";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 41;
		this->atomic_radii = 1.64;
	}
	else if (Symbol.compare("Mo") == 0)
	{
		this->atomic_weight = 95.951;
		this->oxidation_state = 6;
		this->protons = 42;
		this->neutrons = 54;
		this->electrons = 42;
		this->valence_e = 6;
		this->Name = "Molybdenum";
		this->Symbol= "Mo";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 42;
		this->atomic_radii = 1.54;
	}
	else if (Symbol.compare("Tc") == 0)
	{
		this->atomic_weight = 98.0;
		this->oxidation_state = 7;
		this->protons = 43;
		this->neutrons = 55;
		this->electrons = 43;
		this->valence_e = 7;
		this->Name = "Technetium";
		this->Symbol= "Tc";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 43;
		this->atomic_radii = 1.47;
	}
	else if (Symbol.compare("Ru") == 0)
	{
		this->atomic_weight = 101.072;
		this->oxidation_state = 4;
		this->protons = 44;
		this->neutrons = 57;
		this->electrons = 44;
		this->valence_e = 8;
		this->Name = "Ruthenium";
		this->Symbol= "Ru";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 44;
		this->atomic_radii = 1.46;
	}
	else if (Symbol.compare("Rh") == 0)
	{
		this->atomic_weight = 102.905502;
		this->oxidation_state = 3;
		this->protons = 45;
		this->neutrons = 58;
		this->electrons = 45;
		this->valence_e = 9;
		this->Name = "Rhodium";
		this->Symbol= "Rh";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 45;
		this->atomic_radii = 1.42;
	}
	else if (Symbol.compare("Pd") == 0)
	{
		this->atomic_weight = 106.421;
		this->oxidation_state = 4;
		this->protons = 46;
		this->neutrons = 60;
		this->electrons = 46;
		this->valence_e = 10;
		this->Name = "Palladium";
		this->Symbol= "Pd";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 46;
		this->atomic_radii = 1.63;
	}
	else if (Symbol.compare("Ag") == 0)
	{
		this->atomic_weight = 107.86822;
		this->oxidation_state = 1;
		this->protons = 47;
		this->neutrons = 61;
		this->electrons = 47;
		this->valence_e = 11;
		this->Name = "Silver";
		this->Symbol= "Ag";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 47;
		this->atomic_radii = 1.72;
	}
	else if (Symbol.compare("Cd") == 0)
	{
		this->atomic_weight = 112.4144;
		this->oxidation_state = 2;
		this->protons = 48;
		this->neutrons = 64;
		this->electrons = 48;
		this->valence_e = 12;
		this->Name = "Cadmium";
		this->Symbol= "Cd";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 48;
		this->atomic_radii = 1.58;
	}
	else if (Symbol.compare("In") == 0)
	{
		this->atomic_weight = 114.8181;
		this->oxidation_state = 3;
		this->protons = 49;
		this->neutrons = 66;
		this->electrons = 49;
		this->valence_e = 3;
		this->Name = "Indium";
		this->Symbol= "In";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 49;
		this->atomic_radii = 1.93;
	}
	else if (Symbol.compare("Sn") == 0)
	{
		this->atomic_weight = 118.7107;
		this->oxidation_state = 4;
		this->protons = 50;
		this->neutrons = 69;
		this->electrons = 50;
		this->valence_e = 4;
		this->Name = "Tin";
		this->Symbol= "Sn";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 50;
		this->atomic_radii = 2.17;
	}
	else if (Symbol.compare("Sb") == 0)
	{
		this->atomic_weight = 121.7601;
		this->oxidation_state = 5;
		this->protons = 51;
		this->neutrons = 71;
		this->electrons = 51;
		this->valence_e = 5;
		this->Name = "Antimony";
		this->Symbol= "Sb";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 51;
		this->atomic_radii = 2.06;
	}
	else if (Symbol.compare("Te") == 0)
	{
		this->atomic_weight = 127.603;
		this->oxidation_state = 6;
		this->protons = 52;
		this->neutrons = 76;
		this->electrons = 52;
		this->valence_e = 6;
		this->Name = "Tellurium";
		this->Symbol= "Te";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 52;
		this->atomic_radii = 2.06;
	}
	else if (Symbol.compare("I") == 0)
	{
		this->atomic_weight = 126.904473;
		this->oxidation_state = -1;
		this->protons = 53;
		this->neutrons = 74;
		this->electrons = 53;
		this->valence_e = 7;
		this->Name = "Iodine";
		this->Symbol= "I";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 53;
		this->atomic_radii = 1.98;
	}
	else if (Symbol.compare("Xe") == 0)
	{
		this->atomic_weight = 131.2936;
		this->oxidation_state = 0;
		this->protons = 54;
		this->neutrons = 77;
		this->electrons = 54;
		this->valence_e = 8;
		this->Name = "Xenon";
		this->Symbol= "Xe";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 54;
		this->atomic_radii = 2.16;
	}
	else if (Symbol.compare("Cs") == 0)
	{
		this->atomic_weight = 132.905451966;
		this->oxidation_state = 1;
		this->protons = 55;
		this->neutrons = 78;
		this->electrons = 55;
		this->valence_e = 1;
		this->Name = "Caesium";
		this->Symbol= "Cs";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 55;
		this->atomic_radii = 3.43;
	}
	else if (Symbol.compare("Ba") == 0)
	{
		this->atomic_weight = 137.3277;
		this->oxidation_state = 2;
		this->protons = 56;
		this->neutrons = 81;
		this->electrons = 56;
		this->valence_e = 2;
		this->Name = "Barium";
		this->Symbol= "Ba";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 56;
		this->atomic_radii = 2.68;
	}
	else if (Symbol.compare("La") == 0)
	{
		this->atomic_weight = 138.90547;
		this->oxidation_state = 3;
		this->protons = 57;
		this->neutrons = 82;
		this->electrons = 57;
		this->valence_e = 3;
		this->Name = "Lanthanum";
		this->Symbol= "La";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 57;
		this->atomic_radii = 2.07;
	}
	else if (Symbol.compare("Ce") == 0)
	{
		this->atomic_weight = 140.116;
		this->oxidation_state = 4;
		this->protons = 58;
		this->neutrons = 82;
		this->electrons = 58;
		this->valence_e = 4;
		this->Name = "Cerium";
		this->Symbol= "Ce";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 58;
		this->atomic_radii = 2.04;
	}
	else if (Symbol.compare("Pr") == 0)
	{
		this->atomic_weight = 140.907662;
		this->oxidation_state = 4;
		this->protons = 59;
		this->neutrons = 82;
		this->electrons = 59;
		this->valence_e = 5;
		this->Name = "Praseodymium";
		this->Symbol= "Pr";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 59;
		this->atomic_radii = 2.03;
	}
	else if (Symbol.compare("Nd") == 0)
	{
		this->atomic_weight = 144.242;
		this->oxidation_state = 3;
		this->protons = 60;
		this->neutrons = 84;
		this->electrons = 60;
		this->valence_e = 6;
		this->Name = "Neodymium";
		this->Symbol= "Nd";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 60;
		this->atomic_radii = 2.01;
	}
	else if (Symbol.compare("Pm") == 0)
	{
		this->atomic_weight = 145.0;
		this->oxidation_state = 3;
		this->protons = 61;
		this->neutrons = 84;
		this->electrons = 61;
		this->valence_e = 7;
		this->Name = "Promethium";
		this->Symbol= "Pm";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 61;
		this->atomic_radii = 1.99;
	}
	else if (Symbol.compare("Sm") == 0)
	{
		this->atomic_weight = 150.362;
		this->oxidation_state = 3;
		this->protons = 62;
		this->neutrons = 88;
		this->electrons = 62;
		this->valence_e = 8;
		this->Name = "Samarium";
		this->Symbol= "Sm";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 62;
		this->atomic_radii = 1.98;
	}
	else if (Symbol.compare("Eu") == 0)
	{
		this->atomic_weight = 151.964;
		this->oxidation_state = 3;
		this->protons = 63;
		this->neutrons = 89;
		this->electrons = 63;
		this->valence_e = 9;
		this->Name = "Europium";
		this->Symbol= "Eu";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 63;
		this->atomic_radii = 1.98;
	}
	else if (Symbol.compare("Gd") == 0)
	{
		this->atomic_weight = 157.253;
		this->oxidation_state = 3;
		this->protons = 64;
		this->neutrons = 93;
		this->electrons = 64;
		this->valence_e = 10;
		this->Name = "Gadolinium";
		this->Symbol= "Gd";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 64;
		this->atomic_radii = 1.96;
	}
	else if (Symbol.compare("Tb") == 0)
	{
		this->atomic_weight = 158.92535;
		this->oxidation_state = 3;
		this->protons = 65;
		this->neutrons = 94;
		this->electrons = 65;
		this->valence_e = 11;
		this->Name = "Terbium";
		this->Symbol= "Tb";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 65;
		this->atomic_radii = 1.94;
	}
	else if (Symbol.compare("Dy") == 0)
	{
		this->atomic_weight = 162.5001;
		this->oxidation_state = 3;
		this->protons = 66;
		this->neutrons = 97;
		this->electrons = 66;
		this->valence_e = 12;
		this->Name = "Dysprosium";
		this->Symbol= "Dy";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 66;
		this->atomic_radii = 1.92;
	}
	else if (Symbol.compare("Ho") == 0)
	{
		this->atomic_weight = 164.930332;
		this->oxidation_state = 3;
		this->protons = 67;
		this->neutrons = 98;
		this->electrons = 67;
		this->valence_e = 13;
		this->Name = "Holmium";
		this->Symbol= "Ho";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 67;
		this->atomic_radii = 1.92;
	}
	else if (Symbol.compare("Er") == 0)
	{
		this->atomic_weight = 167.259;
		this->oxidation_state = 3;
		this->protons = 68;
		this->neutrons = 99;
		this->electrons = 68;
		this->valence_e = 14;
		this->Name = "Erbium";
		this->Symbol= "Er";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 68;
		this->atomic_radii = 1.89;
	}
	else if (Symbol.compare("Tm") == 0)
	{
		this->atomic_weight = 168.934222;
		this->oxidation_state = 3;
		this->protons = 69;
		this->neutrons = 100;
		this->electrons = 69;
		this->valence_e = 15;
		this->Name = "Thulium";
		this->Symbol= "Tm";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 69;
		this->atomic_radii = 1.90;
	}
	else if (Symbol.compare("Yb") == 0)
	{
		this->atomic_weight = 173.0545;
		this->oxidation_state = 3;
		this->protons = 70;
		this->neutrons = 103;
		this->electrons = 70;
		this->valence_e = 16;
		this->Name = "Ytterbium";
		this->Symbol= "Yb";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 70;
		this->atomic_radii = 1.87;
	}
	else if (Symbol.compare("Lu") == 0)
	{
		this->atomic_weight = 174.96684;
		this->oxidation_state = 3;
		this->protons = 71;
		this->neutrons = 104;
		this->electrons = 71;
		this->valence_e = 3;
		this->Name = "Lutetium";
		this->Symbol= "Lu";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 71;
		this->atomic_radii = 1.87;
	}
	else if (Symbol.compare("Hf") == 0)
	{
		this->atomic_weight = 178.492;
		this->oxidation_state = 4;
		this->protons = 72;
		this->neutrons = 106;
		this->electrons = 72;
		this->valence_e = 4;
		this->Name = "Hafnium";
		this->Symbol= "Hf";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 72;
		this->atomic_radii = 1.75;
	}
	else if (Symbol.compare("Ta") == 0)
	{
		this->atomic_weight = 180.947882;
		this->oxidation_state = 5;
		this->protons = 73;
		this->neutrons = 108;
		this->electrons = 73;
		this->valence_e = 5;
		this->Name = "Tantalum";
		this->Symbol= "Ta";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 73;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("W") == 0)
	{
		this->atomic_weight = 183.841;
		this->oxidation_state = 6;
		this->protons = 74;
		this->neutrons = 110;
		this->electrons = 74;
		this->valence_e = 6;
		this->Name = "Tungsten";
		this->Symbol= "W";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 74;
		this->atomic_radii = 1.62;
	}
	else if (Symbol.compare("Re") == 0)
	{
		this->atomic_weight = 186.2071;
		this->oxidation_state = 7;
		this->protons = 75;
		this->neutrons = 111;
		this->electrons = 75;
		this->valence_e = 7;
		this->Name = "Rhenium";
		this->Symbol= "Re";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 75;
		this->atomic_radii = 1.51;
	}
	else if (Symbol.compare("Os") == 0)
	{
		this->atomic_weight = 190.233;
		this->oxidation_state = 4;
		this->protons = 76;
		this->neutrons = 114;
		this->electrons = 76;
		this->valence_e = 8;
		this->Name = "Osmium";
		this->Symbol= "Os";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 76;
		this->atomic_radii = 1.44;
	}
	else if (Symbol.compare("Ir") == 0)
	{
		this->atomic_weight = 192.2173;
		this->oxidation_state = 4;
		this->protons = 77;
		this->neutrons = 115;
		this->electrons = 77;
		this->valence_e = 9;
		this->Name = "Iridium";
		this->Symbol= "Ir";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 77;
		this->atomic_radii = 1.41;
	}
	else if (Symbol.compare("Pt") == 0)
	{
		this->atomic_weight = 195.0849;
		this->oxidation_state = 4;
		this->protons = 78;
		this->neutrons = 117;
		this->electrons = 78;
		this->valence_e = 10;
		this->Name = "Platinum";
		this->Symbol= "Pt";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 78;
		this->atomic_radii = 1.75;
	}
	else if (Symbol.compare("Au") == 0)
	{
		this->atomic_weight = 196.9665694;
		this->oxidation_state = 3;
		this->protons = 79;
		this->neutrons = 118;
		this->electrons = 79;
		this->valence_e = 11;
		this->Name = "Gold";
		this->Symbol= "Au";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 79;
		this->atomic_radii = 1.66;
	}
	else if (Symbol.compare("Hg") == 0)
	{
		this->atomic_weight = 200.5923;
		this->oxidation_state = 2;
		this->protons = 80;
		this->neutrons = 121;
		this->electrons = 80;
		this->valence_e = 12;
		this->Name = "Mercury";
		this->Symbol= "Hg";
		this->Category = "Transition Metal";
		this->NaturalState = "Liquid";
		this->atomic_number = 80;
		this->atomic_radii = 1.55;
	}
	else if (Symbol.compare("Tl") == 0)
	{
		this->atomic_weight = 204.381;
		this->oxidation_state = 1;
		this->protons = 81;
		this->neutrons = 123;
		this->electrons = 81;
		this->valence_e = 3;
		this->Name = "Thallium";
		this->Symbol= "Tl";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 81;
		this->atomic_radii = 1.96;
	}
	else if (Symbol.compare("Pb") == 0)
	{
		this->atomic_weight = 207.21;
		this->oxidation_state = 2;
		this->protons = 82;
		this->neutrons = 125;
		this->electrons = 82;
		this->valence_e = 4;
		this->Name = "Lead";
		this->Symbol= "Pb";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 82;
		this->atomic_radii = 2.02;
	}
	else if (Symbol.compare("Bi") == 0)
	{
		this->atomic_weight = 208.980401;
		this->oxidation_state = 3;
		this->protons = 83;
		this->neutrons = 126;
		this->electrons = 83;
		this->valence_e = 5;
		this->Name = "Bismuth";
		this->Symbol= "Bi";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 83;
		this->atomic_radii = 2.07;
	}
	else if (Symbol.compare("Po") == 0)
	{
		this->atomic_weight = 209.0;
		this->oxidation_state = 4;
		this->protons = 84;
		this->neutrons = 125;
		this->electrons = 84;
		this->valence_e = 6;
		this->Name = "Polonium";
		this->Symbol= "Po";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 84;
		this->atomic_radii = 1.97;
	}
	else if (Symbol.compare("At") == 0)
	{
		this->atomic_weight = 210.0;
		this->oxidation_state = -1;
		this->protons = 85;
		this->neutrons = 125;
		this->electrons = 85;
		this->valence_e = 7;
		this->Name = "Astatine";
		this->Symbol= "At";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 85;
		this->atomic_radii = 2.02;
	}
	else if (Symbol.compare("Rn") == 0)
	{
		this->atomic_weight = 222.0;
		this->oxidation_state = 0;
		this->protons = 86;
		this->neutrons = 136;
		this->electrons = 86;
		this->valence_e = 8;
		this->Name = "Radon";
		this->Symbol= "Rn";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 86;
		this->atomic_radii = 2.20;
	}
	else if (Symbol.compare("Fr") == 0)
	{
		this->atomic_weight = 223.0;
		this->oxidation_state = 1;
		this->protons = 87;
		this->neutrons = 136;
		this->electrons = 87;
		this->valence_e = 1;
		this->Name = "Francium";
		this->Symbol= "Fr";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 87;
		this->atomic_radii = 3.48;
	}
	else if (Symbol.compare("Ra") == 0)
	{
		this->atomic_weight = 226.0;
		this->oxidation_state = 2;
		this->protons = 88;
		this->neutrons = 138;
		this->electrons = 88;
		this->valence_e = 2;
		this->Name = "Radium";
		this->Symbol= "Ra";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 88;
		this->atomic_radii = 2.83;
	}
	else if (Symbol.compare("Ac") == 0)
	{
		this->atomic_weight = 227.0;
		this->oxidation_state = 3;
		this->protons = 89;
		this->neutrons = 138;
		this->electrons = 89;
		this->valence_e = 3;
		this->Name = "Actinium";
		this->Symbol= "Ac";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 89;
		this->atomic_radii = 2.15;
	}
	else if (Symbol.compare("Th") == 0)
	{
		this->atomic_weight = 232.03774;
		this->oxidation_state = 4;
		this->protons = 90;
		this->neutrons = 142;
		this->electrons = 90;
		this->valence_e = 4;
		this->Name = "Thorium";
		this->Symbol= "Th";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 90;
		this->atomic_radii = 2.06;
	}
	else if (Symbol.compare("Pa") == 0)
	{
		this->atomic_weight = 231.03588;
		this->oxidation_state = 5;
		this->protons = 91;
		this->neutrons = 140;
		this->electrons = 91;
		this->valence_e = 5;
		this->Name = "Protactinium";
		this->Symbol= "Pa";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 91;
		this->atomic_radii = 2.00;
	}
	else if (Symbol.compare("U") == 0)
	{
		this->atomic_weight = 238.028913;
		this->oxidation_state = 6;
		this->protons = 92;
		this->neutrons = 146;
		this->electrons = 92;
		this->valence_e = 6;
		this->Name = "Uranium";
		this->Symbol= "U";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 92;
		this->atomic_radii = 1.86;
	}
	else if (Symbol.compare("Np") == 0)
	{
		this->atomic_weight = 237.0;
		this->oxidation_state = 5;
		this->protons = 93;
		this->neutrons = 144;
		this->electrons = 93;
		this->valence_e = 7;
		this->Name = "Neptunium";
		this->Symbol= "Np";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 93;
		this->atomic_radii = 1.90;
	}
	else if (Symbol.compare("Pu") == 0)
	{
		this->atomic_weight = 244.0;
		this->oxidation_state = 4;
		this->protons = 94;
		this->neutrons = 150;
		this->electrons = 94;
		this->valence_e = 8;
		this->Name = "Plutonium";
		this->Symbol= "Pu";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 94;
		this->atomic_radii = 1.87;
	}
	else if (Symbol.compare("Am") == 0)
	{
		this->atomic_weight = 243.0;
		this->oxidation_state = 3;
		this->protons = 95;
		this->neutrons = 148;
		this->electrons = 95;
		this->valence_e = 9;
		this->Name = "Americium";
		this->Symbol= "Am";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 95;
		this->atomic_radii = 1.80;
	}
	else if (Symbol.compare("Cm") == 0)
	{
		this->atomic_weight = 247.0;
		this->oxidation_state = 3;
		this->protons = 96;
		this->neutrons = 151;
		this->electrons = 96;
		this->valence_e = 10;
		this->Name = "Curium";
		this->Symbol= "Cm";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 96;
		this->atomic_radii = 1.69;
	}
	else if (Symbol.compare("Bk") == 0)
	{
		this->atomic_weight = 247.0;
		this->oxidation_state = 3;
		this->protons = 97;
		this->neutrons = 150;
		this->electrons = 97;
		this->valence_e = 11;
		this->Name = "Berkelium";
		this->Symbol= "Bk";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 97;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("Cf") == 0)
	{
		this->atomic_weight = 251.0;
		this->oxidation_state = 3;
		this->protons = 98;
		this->neutrons = 153;
		this->electrons = 98;
		this->valence_e = 12;
		this->Name = "Californium";
		this->Symbol= "Cf";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 98;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("Es") == 0)
	{
		this->atomic_weight = 252.0;
		this->oxidation_state = 3;
		this->protons = 99;
		this->neutrons = 153;
		this->electrons = 99;
		this->valence_e = 13;
		this->Name = "Einsteinium";
		this->Symbol= "Es";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 99;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("Fm") == 0)
	{
		this->atomic_weight = 257.0;
		this->oxidation_state = 3;
		this->protons = 100;
		this->neutrons = 157;
		this->electrons = 100;
		this->valence_e = 14;
		this->Name = "Fermium";
		this->Symbol= "Fm";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 100;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("Md") == 0)
	{
		this->atomic_weight = 258.0;
		this->oxidation_state = 3;
		this->protons = 101;
		this->neutrons = 157;
		this->electrons = 101;
		this->valence_e = 15;
		this->Name = "Mendelevium";
		this->Symbol= "Md";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 101;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("No") == 0)
	{
		this->atomic_weight = 259.0;
		this->oxidation_state = 2;
		this->protons = 102;
		this->neutrons = 157;
		this->electrons = 102;
		this->valence_e = 16;
		this->Name = "Nobelium";
		this->Symbol= "No";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 102;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("Lr") == 0)
	{
		this->atomic_weight = 266.0;
		this->oxidation_state = 3;
		this->protons = 103;
		this->neutrons = 159;
		this->electrons = 103;
		this->valence_e = 3;
		this->Name = "Lawrencium";
		this->Symbol= "Lr";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 103;
		this->atomic_radii = 1.70;
	}
	else if (Symbol.compare("Rf") == 0)
	{
		this->atomic_weight = 267.0;
		this->oxidation_state = 4;
		this->protons = 104;
		this->neutrons = 157;
		this->electrons = 104;
		this->valence_e = 4;
		this->Name = "Rutherfordium";
		this->Symbol= "Rf";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 104;
		this->atomic_radii = 1.57;
	}
	else if (Symbol.compare("Db") == 0)
	{
		this->atomic_weight = 268.0;
		this->oxidation_state = 5;
		this->protons = 105;
		this->neutrons = 157;
		this->electrons = 105;
		this->valence_e = 5;
		this->Name = "Dubnium";
		this->Symbol= "Db";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 105;
		this->atomic_radii = 1.49;
	}
	else if (Symbol.compare("Sg") == 0)
	{
		this->atomic_weight = 269.0;
		this->oxidation_state = 6;
		this->protons = 106;
		this->neutrons = 157;
		this->electrons = 106;
		this->valence_e = 6;
		this->Name = "Seaborgium";
		this->Symbol= "Sg";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 106;
		this->atomic_radii = 1.43;
	}
	else if (Symbol.compare("Bh") == 0)
	{
		this->atomic_weight = 270.0;
		this->oxidation_state = 7;
		this->protons = 107;
		this->neutrons = 155;
		this->electrons = 107;
		this->valence_e = 7;
		this->Name = "Bohrium";
		this->Symbol= "Bh";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 107;
		this->atomic_radii = 1.41;
	}
	else if (Symbol.compare("Hs") == 0)
	{
		this->atomic_weight = 269.0;
		this->oxidation_state = 8;
		this->protons = 108;
		this->neutrons = 157;
		this->electrons = 108;
		this->valence_e = 8;
		this->Name = "Hassium";
		this->Symbol= "Hs";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 108;
		this->atomic_radii = 1.34;
	}
	else if (Symbol.compare("Mt") == 0)
	{
		this->atomic_weight = 278.0;
		this->oxidation_state = 6;
		this->protons = 109;
		this->neutrons = 157;
		this->electrons = 109;
		this->valence_e = 9;
		this->Name = "Meitnerium";
		this->Symbol= "Mt";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 109;
		this->atomic_radii = 1.29;
	}
	else if (Symbol.compare("Ds") == 0)
	{
		this->atomic_weight = 281.0;
		this->oxidation_state = 8;
		this->protons = 110;
		this->neutrons = 171;
		this->electrons = 110;
		this->valence_e = 10;
		this->Name = "Darmstadium";
		this->Symbol= "Ds";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 110;
		this->atomic_radii = 1.28;
	}
	else if (Symbol.compare("Rg") == 0)
	{
		this->atomic_weight = 281.0;
		this->oxidation_state = 3;
		this->protons = 111;
		this->neutrons = 170;
		this->electrons = 111;
		this->valence_e = 11;
		this->Name = "Roentgenium";
		this->Symbol= "Rg";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 111;
		this->atomic_radii = 1.21;
	}
	else if (Symbol.compare("Cn") == 0)
	{
		this->atomic_weight = 285.0;
		this->oxidation_state = 4;
		this->protons = 112;
		this->neutrons = 173;
		this->electrons = 112;
		this->valence_e = 12;
		this->Name = "Copernicium";
		this->Symbol= "Cn";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 112;
		this->atomic_radii = 1.22;
	}
	else if (Symbol.compare("Uut") == 0)
	{
		this->atomic_weight = 286.0;
		this->oxidation_state = 1;
		this->protons = 113;
		this->neutrons = 173;
		this->electrons = 113;
		this->valence_e = 3;
		this->Name = "Ununtrium";
		this->Symbol= "Uut";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 113;
		this->atomic_radii = 1.76;
	}
	else if (Symbol.compare("Fl") == 0)
	{
		this->atomic_weight = 289.0;
		this->oxidation_state = 2;
		this->protons = 114;
		this->neutrons = 175;
		this->electrons = 114;
		this->valence_e = 4;
		this->Name = "Flerovium";
		this->Symbol= "Fl";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 114;
		this->atomic_radii = 1.74;
	}
	else if (Symbol.compare("Uup") == 0)
	{
		this->atomic_weight = 289.0;
		this->oxidation_state = 1;
		this->protons = 115;
		this->neutrons = 174;
		this->electrons = 115;
		this->valence_e = 5;
		this->Name = "Ununpentium";
		this->Symbol= "Uup";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 115;
		this->atomic_radii = 1.57;
	}
	else if (Symbol.compare("Lv") == 0)
	{
		this->atomic_weight = 293.0;
		this->oxidation_state = 2;
		this->protons = 116;
		this->neutrons = 177;
		this->electrons = 116;
		this->valence_e = 6;
		this->Name = "Livermorium";
		this->Symbol= "Lv";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 116;
		this->atomic_radii = 1.64;
	}
	else if (Symbol.compare("Uus") == 0)
	{
		this->atomic_weight = 294.0;
		this->oxidation_state = 1;
		this->protons = 117;
		this->neutrons = 177;
		this->electrons = 117;
		this->valence_e = 7;
		this->Name = "Ununseptium";
		this->Symbol= "Uus";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 117;
		this->atomic_radii = 1.57;
	}
	else if (Symbol.compare("Uuo") == 0)
	{
		this->atomic_weight = 294.0;
		this->oxidation_state = 0;
		this->protons = 118;
		this->neutrons = 176;
		this->electrons = 118;
		this->valence_e = 8;
		this->Name = "Ununoctium";
		this->Symbol= "Uuo";
		this->Category = "Nobel Gas";
		this->NaturalState = "Solid";
		this->atomic_number = 118;
		this->atomic_radii = 1.57;
	}
	else
	{
		std::cout << "Given Symbol: " << Symbol << std::endl;
		std::cout << "Symbol registration is case sensitive!" << std::endl;
		mError(invalid_atom);
		this->atomic_weight = 0.0;
		this->oxidation_state = 0;
		this->protons = 0;
		this->neutrons = 0;
		this->electrons = 0;
		this->Name = "No Name";
		this->Symbol= "N/A";
		this->Category = "N/A";
		this->NaturalState = "N/A";
		this->atomic_number = 0;
	}
}

//Register a instance of an Atom object by it's symbol
void Atom::Register(int number)
{
	if (number == 1)
	{
		this->atomic_weight = 1.0081;
		this->oxidation_state = 1;
		this->protons = 1;
		this->neutrons = 0;
		this->electrons = 1;
		this->valence_e = 1;
		this->Name = "Hydrogen";
		this->Symbol= "H";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 1;
		this->atomic_radii = 1.2;
	}
	else if (number == 2)
	{
		this->atomic_weight = 4.0026022;
		this->oxidation_state = 0;
		this->protons = 2;
		this->neutrons = 2;
		this->electrons = 2;
		this->valence_e = 2;
		this->Name = "Helium";
		this->Symbol= "He";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 2;
		this->atomic_radii = 1.4;
	}
	else if (number == 3)
	{
		this->atomic_weight = 6.941;
		this->oxidation_state = 1;
		this->protons = 3;
		this->neutrons = 4;
		this->electrons = 3;
		this->valence_e = 1;
		this->Name = "Lithium";
		this->Symbol= "Li";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 3;
		this->atomic_radii = 1.82;
	}
	else if (number == 4)
	{
		this->atomic_weight = 9.01218315;
		this->oxidation_state = 2;
		this->protons = 4;
		this->neutrons = 5;
		this->electrons = 4;
		this->valence_e = 2;
		this->Name = "Beryllium";
		this->Symbol= "Be";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 4;
		this->atomic_radii = 1.53;
	}
	else if (number == 5)
	{
		this->atomic_weight = 10.811;
		this->oxidation_state = 3;
		this->protons = 5;
		this->neutrons = 6;
		this->electrons = 5;
		this->valence_e = 3;
		this->Name = "Boron";
		this->Symbol= "B";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 5;
		this->atomic_radii = 1.92;
	}
	else if (number == 6)
	{
		this->atomic_weight = 12.0111;
		this->oxidation_state = 4;
		this->protons = 6;
		this->neutrons = 6;
		this->electrons = 6;
		this->valence_e = 4;
		this->Name = "Carbon";
		this->Symbol= "C";
		this->Category = "Polyatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 6;
		this->atomic_radii = 1.70;
	}
	else if (number == 7)
	{
		this->atomic_weight = 14.0071;
		this->oxidation_state = -3;
		this->protons = 7;
		this->neutrons = 7;
		this->electrons = 7;
		this->valence_e = 5;
		this->Name = "Nitrogen";
		this->Symbol= "N";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 7;
		this->atomic_radii = 1.55;
	}
	else if (number == 8)
	{
		this->atomic_weight = 15.9994;
		this->oxidation_state = -2;
		this->protons = 8;
		this->neutrons = 8;
		this->electrons = 8;
		this->valence_e = 6;
		this->Name = "Oxygen";
		this->Symbol= "O";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 8;
		this->atomic_radii = 1.52;
	}
	else if (number == 9)
	{
		this->atomic_weight = 18.9984031636;
		this->oxidation_state = -1;
		this->protons = 9;
		this->neutrons = 10;
		this->electrons = 9;
		this->valence_e = 7;
		this->Name = "Fluorine";
		this->Symbol= "F";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 9;
		this->atomic_radii = 1.35;
	}
	else if (number == 10)
	{
		this->atomic_weight = 20.17976;
		this->oxidation_state = 0;
		this->protons = 10;
		this->neutrons = 10;
		this->electrons = 10;
		this->valence_e = 8;
		this->Name = "Neon";
		this->Symbol= "Ne";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 10;
		this->atomic_radii = 1.54;
	}
	else if (number == 11)
	{
		this->atomic_weight = 22.989769282;
		this->oxidation_state = 1;
		this->protons = 11;
		this->neutrons = 12;
		this->electrons = 11;
		this->valence_e = 1;
		this->Name = "Sodium";
		this->Symbol= "Na";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 11;
		this->atomic_radii = 2.27;
	}
	else if (number == 12)
	{
		this->atomic_weight = 24.3051;
		this->oxidation_state = 2;
		this->protons = 12;
		this->neutrons = 12;
		this->electrons = 12;
		this->valence_e = 2;
		this->Name = "Magnesium";
		this->Symbol= "Mg";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 12;
		this->atomic_radii = 1.73;
	}
	else if (number == 13)
	{
		this->atomic_weight = 26.98153857;
		this->oxidation_state = 3;
		this->protons = 13;
		this->neutrons = 14;
		this->electrons = 13;
		this->valence_e = 3;
		this->Name = "Aluminium";
		this->Symbol= "Al";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 13;
		this->atomic_radii = 1.84;
	}
	else if (number == 14)
	{
		this->atomic_weight = 28.0851;
		this->oxidation_state = 4;
		this->protons = 14;
		this->neutrons = 14;
		this->electrons = 14;
		this->valence_e = 4;
		this->Name = "Silicon";
		this->Symbol= "Si";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 14;
		this->atomic_radii = 2.10;
	}
	else if (number == 15)
	{
		this->atomic_weight = 30.9737619985;
		this->oxidation_state = 5;
		this->protons = 15;
		this->neutrons = 16;
		this->electrons = 15;
		this->valence_e = 5;
		this->Name = "Phosphorus";
		this->Symbol= "P";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 15;
		this->atomic_radii = 1.80;
	}
	else if (number == 16)
	{
		this->atomic_weight = 32.061;
		this->oxidation_state = 6;
		this->protons = 16;
		this->neutrons = 16;
		this->electrons = 16;
		this->valence_e = 6;
		this->Name = "Sulfur";
		this->Symbol= "S";
		this->Category = "Polyatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 16;
		this->atomic_radii = 1.80;
	}
	else if (number == 17)
	{
		this->atomic_weight = 35.451;
		this->oxidation_state = -1;
		this->protons = 17;
		this->neutrons = 18;
		this->electrons = 17;
		this->valence_e = 7;
		this->Name = "Chlorine";
		this->Symbol= "Cl";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 17;
		this->atomic_radii = 1.75;
	}
	else if (number == 18)
	{
		this->atomic_weight = 39.9481;
		this->oxidation_state = 0;
		this->protons = 18;
		this->neutrons = 22;
		this->electrons = 18;
		this->valence_e = 8;
		this->Name = "Argon";
		this->Symbol= "Ar";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 18;
		this->atomic_radii = 1.88;
	}
	else if (number == 19)
	{
		this->atomic_weight = 39.09831;
		this->oxidation_state = 1;
		this->protons = 19;
		this->neutrons = 20;
		this->electrons = 19;
		this->valence_e = 1;
		this->Name = "Potassium";
		this->Symbol= "K";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 19;
		this->atomic_radii = 2.75;
	}
	else if (number == 20)
	{
		this->atomic_weight = 40.0784;
		this->oxidation_state = 2;
		this->protons = 20;
		this->neutrons = 20;
		this->electrons = 20;
		this->valence_e = 2;
		this->Name = "Calcium";
		this->Symbol= "Ca";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 20;
		this->atomic_radii = 2.31;
	}
	else if (number == 21)
	{
		this->atomic_weight = 44.9559085;
		this->oxidation_state = 3;
		this->protons = 21;
		this->neutrons = 24;
		this->electrons = 21;
		this->valence_e = 3;
		this->Name = "Scandium";
		this->Symbol= "Sc";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 21;
		this->atomic_radii = 2.11;
	}
	else if (number == 22)
	{
		this->atomic_weight = 47.8671;
		this->oxidation_state = 4;
		this->protons = 22;
		this->neutrons = 26;
		this->electrons = 22;
		this->valence_e = 4;
		this->Name = "Titanium";
		this->Symbol= "Ti";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 22;
		this->atomic_radii = 1.60;
	}
	else if (number == 23)
	{
		this->atomic_weight = 50.94151;
		this->oxidation_state = 5;
		this->protons = 23;
		this->neutrons = 28;
		this->electrons = 23;
		this->valence_e = 5;
		this->Name = "Vanadium";
		this->Symbol= "V";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 23;
		this->atomic_radii = 1.53;
	}
	else if (number == 24)
	{
		this->atomic_weight = 51.99616;
		this->oxidation_state = 6;
		this->protons = 24;
		this->neutrons = 28;
		this->electrons = 24;
		this->valence_e = 6;
		this->Name = "Chromium";
		this->Symbol= "Cr";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 24;
		this->atomic_radii = 1.39;
	}
	else if (number == 25)
	{
		this->atomic_weight = 54.9380443;
		this->oxidation_state = 7;
		this->protons = 25;
		this->neutrons = 30;
		this->electrons = 25;
		this->valence_e = 7;
		this->Name = "Manganese";
		this->Symbol= "Mn";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 25;
		this->atomic_radii = 1.39;
	}
	else if (number == 26)
	{
		this->atomic_weight = 55.8452;
		this->oxidation_state = 3;
		this->protons = 26;
		this->neutrons = 30;
		this->electrons = 26;
		this->valence_e = 8;
		this->Name = "Iron";
		this->Symbol= "Fe";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 26;
		this->atomic_radii = 1.32;
	}
	else if (number == 27)
	{
		this->atomic_weight = 58.9331944;
		this->oxidation_state = 3;
		this->protons = 27;
		this->neutrons = 32;
		this->electrons = 27;
		this->valence_e = 9;
		this->Name = "Cobalt";
		this->Symbol= "Co";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 27;
		this->atomic_radii = 1.26;
	}
	else if (number == 28)
	{
		this->atomic_weight = 58.69344;
		this->oxidation_state = 2;
		this->protons = 28;
		this->neutrons = 31;
		this->electrons = 28;
		this->valence_e = 10;
		this->Name = "Nickel";
		this->Symbol= "Ni";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 28;
		this->atomic_radii = 1.63;
	}
	else if (number == 29)
	{
		this->atomic_weight = 63.5463;
		this->oxidation_state = 2;
		this->protons = 29;
		this->neutrons = 35;
		this->electrons = 29;
		this->valence_e = 11;
		this->Name = "Copper";
		this->Symbol= "Cu";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 29;
		this->atomic_radii = 1.40;
	}
	else if (number == 30)
	{
		this->atomic_weight = 65.382;
		this->oxidation_state = 2;
		this->protons = 30;
		this->neutrons = 35;
		this->electrons = 30;
		this->valence_e = 12;
		this->Name = "Zinc";
		this->Symbol= "Zn";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 30;
		this->atomic_radii = 1.39;
	}
	else if (number == 31)
	{
		this->atomic_weight = 69.7231;
		this->oxidation_state = 3;
		this->protons = 31;
		this->neutrons = 39;
		this->electrons = 31;
		this->valence_e = 3;
		this->Name = "Gallium";
		this->Symbol= "Ga";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 31;
		this->atomic_radii = 1.87;
	}
	else if (number == 32)
	{
		this->atomic_weight = 72.6308;
		this->oxidation_state = 4;
		this->protons = 32;
		this->neutrons = 41;
		this->electrons = 32;
		this->valence_e = 4;
		this->Name = "Germanium";
		this->Symbol= "Ge";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 32;
		this->atomic_radii = 2.11;
	}
	else if (number == 33)
	{
		this->atomic_weight = 74.9215956;
		this->oxidation_state = 5;
		this->protons = 33;
		this->neutrons = 42;
		this->electrons = 33;
		this->valence_e = 5;
		this->Name = "Arsenic";
		this->Symbol= "As";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 33;
		this->atomic_radii = 1.85;
	}
	else if (number == 34)
	{
		this->atomic_weight = 78.9718;
		this->oxidation_state = 6;
		this->protons = 34;
		this->neutrons = 45;
		this->electrons = 34;
		this->valence_e = 6;
		this->Name = "Selenium";
		this->Symbol= "Se";
		this->Category = "Polyatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 34;
		this->atomic_radii = 1.90;
	}
	else if (number == 35)
	{
		this->atomic_weight = 79.9041;
		this->oxidation_state = -1;
		this->protons = 35;
		this->neutrons = 45;
		this->electrons = 35;
		this->valence_e = 7;
		this->Name = "Bromine";
		this->Symbol= "Br";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Liquid";
		this->atomic_number = 35;
		this->atomic_radii = 1.85;
	}
	else if (number == 36)
	{
		this->atomic_weight = 83.798;
		this->oxidation_state = 0;
		this->protons = 36;
		this->neutrons = 48;
		this->electrons = 36;
		this->valence_e = 8;
		this->Name = "Krypton";
		this->Symbol= "Kr";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 36;
		this->atomic_radii = 2.02;
	}
	else if (number == 37)
	{
		this->atomic_weight = 85.46783;
		this->oxidation_state = 1;
		this->protons = 37;
		this->neutrons = 48;
		this->electrons = 37;
		this->valence_e = 1;
		this->Name = "Rubidium";
		this->Symbol= "Rb";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 37;
		this->atomic_radii = 3.03;
	}
	else if (number == 38)
	{
		this->atomic_weight = 87.621;
		this->oxidation_state = 2;
		this->protons = 38;
		this->neutrons = 50;
		this->electrons = 38;
		this->valence_e = 2;
		this->Name = "Strontium";
		this->Symbol= "Sr";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 38;
		this->atomic_radii = 2.49;
	}
	else if (number == 39)
	{
		this->atomic_weight = 88.905842;
		this->oxidation_state = 3;
		this->protons = 39;
		this->neutrons = 50;
		this->electrons = 39;
		this->valence_e = 3;
		this->Name = "Yttrium";
		this->Symbol= "Y";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 39;
		this->atomic_radii = 1.90;
	}
	else if (number == 40)
	{
		this->atomic_weight = 91.2242;
		this->oxidation_state = 4;
		this->protons = 40;
		this->neutrons = 51;
		this->electrons = 40;
		this->valence_e = 4;
		this->Name = "Zirconium";
		this->Symbol= "Zr";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 40;
		this->atomic_radii = 1.75;
	}
	else if (number == 41)
	{
		this->atomic_weight = 92.906372;
		this->oxidation_state = 5;
		this->protons = 41;
		this->neutrons = 52;
		this->electrons = 41;
		this->valence_e = 5;
		this->Name = "Niobium";
		this->Symbol= "Nb";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 41;
		this->atomic_radii = 1.64;
	}
	else if (number == 42)
	{
		this->atomic_weight = 95.951;
		this->oxidation_state = 6;
		this->protons = 42;
		this->neutrons = 54;
		this->electrons = 42;
		this->valence_e = 6;
		this->Name = "Molybdenum";
		this->Symbol= "Mo";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 42;
		this->atomic_radii = 1.54;
	}
	else if (number == 43)
	{
		this->atomic_weight = 98.0;
		this->oxidation_state = 7;
		this->protons = 43;
		this->neutrons = 55;
		this->electrons = 43;
		this->valence_e = 7;
		this->Name = "Technetium";
		this->Symbol= "Tc";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 43;
		this->atomic_radii = 1.47;
	}
	else if (number == 44)
	{
		this->atomic_weight = 101.072;
		this->oxidation_state = 4;
		this->protons = 44;
		this->neutrons = 57;
		this->electrons = 44;
		this->valence_e = 8;
		this->Name = "Ruthenium";
		this->Symbol= "Ru";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 44;
		this->atomic_radii = 1.46;
	}
	else if (number == 45)
	{
		this->atomic_weight = 102.905502;
		this->oxidation_state = 3;
		this->protons = 45;
		this->neutrons = 58;
		this->electrons = 45;
		this->valence_e = 9;
		this->Name = "Rhodium";
		this->Symbol= "Rh";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 45;
		this->atomic_radii = 1.42;
	}
	else if (number == 46)
	{
		this->atomic_weight = 106.421;
		this->oxidation_state = 4;
		this->protons = 46;
		this->neutrons = 60;
		this->electrons = 46;
		this->valence_e = 10;
		this->Name = "Palladium";
		this->Symbol= "Pd";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 46;
		this->atomic_radii = 1.63;
	}
	else if (number == 47)
	{
		this->atomic_weight = 107.86822;
		this->oxidation_state = 1;
		this->protons = 47;
		this->neutrons = 61;
		this->electrons = 47;
		this->valence_e = 11;
		this->Name = "Silver";
		this->Symbol= "Ag";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 47;
		this->atomic_radii = 1.72;
	}
	else if (number == 48)
	{
		this->atomic_weight = 112.4144;
		this->oxidation_state = 2;
		this->protons = 48;
		this->neutrons = 64;
		this->electrons = 48;
		this->valence_e = 12;
		this->Name = "Cadmium";
		this->Symbol= "Cd";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 48;
		this->atomic_radii = 1.58;
	}
	else if (number == 49)
	{
		this->atomic_weight = 114.8181;
		this->oxidation_state = 3;
		this->protons = 49;
		this->neutrons = 66;
		this->electrons = 49;
		this->valence_e = 3;
		this->Name = "Indium";
		this->Symbol= "In";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 49;
		this->atomic_radii = 1.93;
	}
	else if (number == 50)
	{
		this->atomic_weight = 118.7107;
		this->oxidation_state = 4;
		this->protons = 50;
		this->neutrons = 69;
		this->electrons = 50;
		this->valence_e = 4;
		this->Name = "Tin";
		this->Symbol= "Sn";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 50;
		this->atomic_radii = 2.17;
	}
	else if (number == 51)
	{
		this->atomic_weight = 121.7601;
		this->oxidation_state = 5;
		this->protons = 51;
		this->neutrons = 71;
		this->electrons = 51;
		this->valence_e = 5;
		this->Name = "Antimony";
		this->Symbol= "Sb";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 51;
		this->atomic_radii = 2.06;
	}
	else if (number == 52)
	{
		this->atomic_weight = 127.603;
		this->oxidation_state = 6;
		this->protons = 52;
		this->neutrons = 76;
		this->electrons = 52;
		this->valence_e = 6;
		this->Name = "Tellurium";
		this->Symbol= "Te";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 52;
		this->atomic_radii = 2.06;
	}
	else if (number == 53)
	{
		this->atomic_weight = 126.904473;
		this->oxidation_state = -1;
		this->protons = 53;
		this->neutrons = 74;
		this->electrons = 53;
		this->valence_e = 7;
		this->Name = "Iodine";
		this->Symbol= "I";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Solid";
		this->atomic_number = 53;
		this->atomic_radii = 1.98;
	}
	else if (number == 54)
	{
		this->atomic_weight = 131.2936;
		this->oxidation_state = 0;
		this->protons = 54;
		this->neutrons = 77;
		this->electrons = 54;
		this->valence_e = 8;
		this->Name = "Xenon";
		this->Symbol= "Xe";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 54;
		this->atomic_radii = 2.16;
	}
	else if (number == 55)
	{
		this->atomic_weight = 132.905451966;
		this->oxidation_state = 1;
		this->protons = 55;
		this->neutrons = 78;
		this->electrons = 55;
		this->valence_e = 1;
		this->Name = "Caesium";
		this->Symbol= "Cs";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 55;
		this->atomic_radii = 3.43;
	}
	else if (number == 56)
	{
		this->atomic_weight = 137.3277;
		this->oxidation_state = 2;
		this->protons = 56;
		this->neutrons = 81;
		this->electrons = 56;
		this->valence_e = 2;
		this->Name = "Barium";
		this->Symbol= "Ba";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 56;
		this->atomic_radii = 2.68;
	}
	else if (number == 57)
	{
		this->atomic_weight = 138.90547;
		this->oxidation_state = 3;
		this->protons = 57;
		this->neutrons = 82;
		this->electrons = 57;
		this->valence_e = 3;
		this->Name = "Lanthanum";
		this->Symbol= "La";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 57;
		this->atomic_radii = 2.07;
	}
	else if (number == 58)
	{
		this->atomic_weight = 140.116;
		this->oxidation_state = 4;
		this->protons = 58;
		this->neutrons = 82;
		this->electrons = 58;
		this->valence_e = 4;
		this->Name = "Cerium";
		this->Symbol= "Ce";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 58;
		this->atomic_radii = 2.04;
	}
	else if (number == 59)
	{
		this->atomic_weight = 140.907662;
		this->oxidation_state = 4;
		this->protons = 59;
		this->neutrons = 82;
		this->electrons = 59;
		this->valence_e = 5;
		this->Name = "Praseodymium";
		this->Symbol= "Pr";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 59;
		this->atomic_radii = 2.03;
	}
	else if (number == 60)
	{
		this->atomic_weight = 144.242;
		this->oxidation_state = 3;
		this->protons = 60;
		this->neutrons = 84;
		this->electrons = 60;
		this->valence_e = 6;
		this->Name = "Neodymium";
		this->Symbol= "Nd";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 60;
		this->atomic_radii = 2.01;
	}
	else if (number == 61)
	{
		this->atomic_weight = 145.0;
		this->oxidation_state = 3;
		this->protons = 61;
		this->neutrons = 84;
		this->electrons = 61;
		this->valence_e = 7;
		this->Name = "Promethium";
		this->Symbol= "Pm";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 61;
		this->atomic_radii = 1.99;
	}
	else if (number == 62)
	{
		this->atomic_weight = 150.362;
		this->oxidation_state = 3;
		this->protons = 62;
		this->neutrons = 88;
		this->electrons = 62;
		this->valence_e = 8;
		this->Name = "Samarium";
		this->Symbol= "Sm";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 62;
		this->atomic_radii = 1.98;
	}
	else if (number == 63)
	{
		this->atomic_weight = 151.964;
		this->oxidation_state = 3;
		this->protons = 63;
		this->neutrons = 89;
		this->electrons = 63;
		this->valence_e = 9;
		this->Name = "Europium";
		this->Symbol= "Eu";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 63;
		this->atomic_radii = 1.98;
	}
	else if (number == 64)
	{
		this->atomic_weight = 157.253;
		this->oxidation_state = 3;
		this->protons = 64;
		this->neutrons = 93;
		this->electrons = 64;
		this->valence_e = 10;
		this->Name = "Gadolinium";
		this->Symbol= "Gd";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 64;
		this->atomic_radii = 1.96;
	}
	else if (number == 65)
	{
		this->atomic_weight = 158.92535;
		this->oxidation_state = 3;
		this->protons = 65;
		this->neutrons = 94;
		this->electrons = 65;
		this->valence_e = 11;
		this->Name = "Terbium";
		this->Symbol= "Tb";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 65;
		this->atomic_radii = 1.94;
	}
	else if (number == 66)
	{
		this->atomic_weight = 162.5001;
		this->oxidation_state = 3;
		this->protons = 66;
		this->neutrons = 97;
		this->electrons = 66;
		this->valence_e = 12;
		this->Name = "Dysprosium";
		this->Symbol= "Dy";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 66;
		this->atomic_radii = 1.92;
	}
	else if (number == 67)
	{
		this->atomic_weight = 164.930332;
		this->oxidation_state = 3;
		this->protons = 67;
		this->neutrons = 98;
		this->electrons = 67;
		this->valence_e = 13;
		this->Name = "Holmium";
		this->Symbol= "Ho";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 67;
		this->atomic_radii = 1.92;
	}
	else if (number == 68)
	{
		this->atomic_weight = 167.259;
		this->oxidation_state = 3;
		this->protons = 68;
		this->neutrons = 99;
		this->electrons = 68;
		this->valence_e = 14;
		this->Name = "Erbium";
		this->Symbol= "Er";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 68;
		this->atomic_radii = 1.89;
	}
	else if (number == 69)
	{
		this->atomic_weight = 168.934222;
		this->oxidation_state = 3;
		this->protons = 69;
		this->neutrons = 100;
		this->electrons = 69;
		this->valence_e = 15;
		this->Name = "Thulium";
		this->Symbol= "Tm";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 69;
		this->atomic_radii = 1.90;
	}
	else if (number == 70)
	{
		this->atomic_weight = 173.0545;
		this->oxidation_state = 3;
		this->protons = 70;
		this->neutrons = 103;
		this->electrons = 70;
		this->valence_e = 16;
		this->Name = "Ytterbium";
		this->Symbol= "Yb";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 70;
		this->atomic_radii = 1.87;
	}
	else if (number == 71)
	{
		this->atomic_weight = 174.96684;
		this->oxidation_state = 3;
		this->protons = 71;
		this->neutrons = 104;
		this->electrons = 71;
		this->valence_e = 3;
		this->Name = "Lutetium";
		this->Symbol= "Lu";
		this->Category = "Lanthanide";
		this->NaturalState = "Solid";
		this->atomic_number = 71;
		this->atomic_radii = 1.87;
	}
	else if (number == 72)
	{
		this->atomic_weight = 178.492;
		this->oxidation_state = 4;
		this->protons = 72;
		this->neutrons = 106;
		this->electrons = 72;
		this->valence_e = 4;
		this->Name = "Hafnium";
		this->Symbol= "Hf";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 72;
		this->atomic_radii = 1.75;
	}
	else if (number == 73)
	{
		this->atomic_weight = 180.947882;
		this->oxidation_state = 5;
		this->protons = 73;
		this->neutrons = 108;
		this->electrons = 73;
		this->valence_e = 5;
		this->Name = "Tantalum";
		this->Symbol= "Ta";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 73;
		this->atomic_radii = 1.70;
	}
	else if (number == 74)
	{
		this->atomic_weight = 183.841;
		this->oxidation_state = 6;
		this->protons = 74;
		this->neutrons = 110;
		this->electrons = 74;
		this->valence_e = 6;
		this->Name = "Tungsten";
		this->Symbol= "W";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 74;
		this->atomic_radii = 1.62;
	}
	else if (number == 75)
	{
		this->atomic_weight = 186.2071;
		this->oxidation_state = 7;
		this->protons = 75;
		this->neutrons = 111;
		this->electrons = 75;
		this->valence_e = 7;
		this->Name = "Rhenium";
		this->Symbol= "Re";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 75;
		this->atomic_radii = 1.51;
	}
	else if (number == 76)
	{
		this->atomic_weight = 190.233;
		this->oxidation_state = 4;
		this->protons = 76;
		this->neutrons = 114;
		this->electrons = 76;
		this->valence_e = 8;
		this->Name = "Osmium";
		this->Symbol= "Os";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 76;
		this->atomic_radii = 1.44;
	}
	else if (number == 77)
	{
		this->atomic_weight = 192.2173;
		this->oxidation_state = 4;
		this->protons = 77;
		this->neutrons = 115;
		this->electrons = 77;
		this->valence_e = 9;
		this->Name = "Iridium";
		this->Symbol= "Ir";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 77;
		this->atomic_radii = 1.41;
	}
	else if (number == 78)
	{
		this->atomic_weight = 195.0849;
		this->oxidation_state = 4;
		this->protons = 78;
		this->neutrons = 117;
		this->electrons = 78;
		this->valence_e = 10;
		this->Name = "Platinum";
		this->Symbol= "Pt";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 78;
		this->atomic_radii = 1.75;
	}
	else if (number == 79)
	{
		this->atomic_weight = 196.9665694;
		this->oxidation_state = 3;
		this->protons = 79;
		this->neutrons = 118;
		this->electrons = 79;
		this->valence_e = 11;
		this->Name = "Gold";
		this->Symbol= "Au";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 79;
		this->atomic_radii = 1.66;
	}
	else if (number == 80)
	{
		this->atomic_weight = 200.5923;
		this->oxidation_state = 2;
		this->protons = 80;
		this->neutrons = 121;
		this->electrons = 80;
		this->valence_e = 12;
		this->Name = "Mercury";
		this->Symbol= "Hg";
		this->Category = "Transition Metal";
		this->NaturalState = "Liquid";
		this->atomic_number = 80;
		this->atomic_radii = 1.55;
	}
	else if (number == 81)
	{
		this->atomic_weight = 204.381;
		this->oxidation_state = 1;
		this->protons = 81;
		this->neutrons = 123;
		this->electrons = 81;
		this->valence_e = 3;
		this->Name = "Thallium";
		this->Symbol= "Tl";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 81;
		this->atomic_radii = 1.96;
	}
	else if (number == 82)
	{
		this->atomic_weight = 207.21;
		this->oxidation_state = 2;
		this->protons = 82;
		this->neutrons = 125;
		this->electrons = 82;
		this->valence_e = 4;
		this->Name = "Lead";
		this->Symbol= "Pb";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 82;
		this->atomic_radii = 2.02;
	}
	else if (number == 83)
	{
		this->atomic_weight = 208.980401;
		this->oxidation_state = 3;
		this->protons = 83;
		this->neutrons = 126;
		this->electrons = 83;
		this->valence_e = 5;
		this->Name = "Bismuth";
		this->Symbol= "Bi";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 83;
		this->atomic_radii = 2.07;
	}
	else if (number == 84)
	{
		this->atomic_weight = 209.0;
		this->oxidation_state = 4;
		this->protons = 84;
		this->neutrons = 125;
		this->electrons = 84;
		this->valence_e = 6;
		this->Name = "Polonium";
		this->Symbol= "Po";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 84;
		this->atomic_radii = 1.97;
	}
	else if (number == 85)
	{
		this->atomic_weight = 210.0;
		this->oxidation_state = -1;
		this->protons = 85;
		this->neutrons = 125;
		this->electrons = 85;
		this->valence_e = 7;
		this->Name = "Astatine";
		this->Symbol= "At";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 85;
		this->atomic_radii = 2.02;
	}
	else if (number == 86)
	{
		this->atomic_weight = 222.0;
		this->oxidation_state = 0;
		this->protons = 86;
		this->neutrons = 136;
		this->electrons = 86;
		this->valence_e = 8;
		this->Name = "Radon";
		this->Symbol= "Rn";
		this->Category = "Nobel Gas";
		this->NaturalState = "Gas";
		this->atomic_number = 86;
		this->atomic_radii = 2.20;
	}
	else if (number == 87)
	{
		this->atomic_weight = 223.0;
		this->oxidation_state = 1;
		this->protons = 87;
		this->neutrons = 136;
		this->electrons = 87;
		this->valence_e = 1;
		this->Name = "Francium";
		this->Symbol= "Fr";
		this->Category = "Alkali Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 87;
		this->atomic_radii = 3.48;
	}
	else if (number == 88)
	{
		this->atomic_weight = 226.0;
		this->oxidation_state = 2;
		this->protons = 88;
		this->neutrons = 138;
		this->electrons = 88;
		this->valence_e = 2;
		this->Name = "Radium";
		this->Symbol= "Ra";
		this->Category = "Alkaline Earth Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 88;
		this->atomic_radii = 2.83;
	}
	else if (number == 89)
	{
		this->atomic_weight = 227.0;
		this->oxidation_state = 3;
		this->protons = 89;
		this->neutrons = 138;
		this->electrons = 89;
		this->valence_e = 3;
		this->Name = "Actinium";
		this->Symbol= "Ac";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 89;
		this->atomic_radii = 2.15;
	}
	else if (number == 90)
	{
		this->atomic_weight = 232.03774;
		this->oxidation_state = 4;
		this->protons = 90;
		this->neutrons = 142;
		this->electrons = 90;
		this->valence_e = 4;
		this->Name = "Thorium";
		this->Symbol= "Th";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 90;
		this->atomic_radii = 2.06;
	}
	else if (number == 91)
	{
		this->atomic_weight = 231.03588;
		this->oxidation_state = 5;
		this->protons = 91;
		this->neutrons = 140;
		this->electrons = 91;
		this->valence_e = 5;
		this->Name = "Protactinium";
		this->Symbol= "Pa";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 91;
		this->atomic_radii = 2.00;
	}
	else if (number == 92)
	{
		this->atomic_weight = 238.028913;
		this->oxidation_state = 6;
		this->protons = 92;
		this->neutrons = 146;
		this->electrons = 92;
		this->valence_e = 6;
		this->Name = "Uranium";
		this->Symbol= "U";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 92;
		this->atomic_radii = 1.86;
	}
	else if (number == 93)
	{
		this->atomic_weight = 237.0;
		this->oxidation_state = 5;
		this->protons = 93;
		this->neutrons = 144;
		this->electrons = 93;
		this->valence_e = 7;
		this->Name = "Neptunium";
		this->Symbol= "Np";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 93;
		this->atomic_radii = 1.90;
	}
	else if (number == 94)
	{
		this->atomic_weight = 244.0;
		this->oxidation_state = 4;
		this->protons = 94;
		this->neutrons = 150;
		this->electrons = 94;
		this->valence_e = 8;
		this->Name = "Plutonium";
		this->Symbol= "Pu";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 94;
		this->atomic_radii = 1.87;
	}
	else if (number == 95)
	{
		this->atomic_weight = 243.0;
		this->oxidation_state = 3;
		this->protons = 95;
		this->neutrons = 148;
		this->electrons = 95;
		this->valence_e = 9;
		this->Name = "Americium";
		this->Symbol= "Am";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 95;
		this->atomic_radii = 1.80;
	}
	else if (number == 96)
	{
		this->atomic_weight = 247.0;
		this->oxidation_state = 3;
		this->protons = 96;
		this->neutrons = 151;
		this->electrons = 96;
		this->valence_e = 10;
		this->Name = "Curium";
		this->Symbol= "Cm";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 96;
		this->atomic_radii = 1.69;
	}
	else if (number == 97)
	{
		this->atomic_weight = 247.0;
		this->oxidation_state = 3;
		this->protons = 97;
		this->neutrons = 150;
		this->electrons = 97;
		this->valence_e = 11;
		this->Name = "Berkelium";
		this->Symbol= "Bk";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 97;
		this->atomic_radii = 1.70;
	}
	else if (number == 98)
	{
		this->atomic_weight = 251.0;
		this->oxidation_state = 3;
		this->protons = 98;
		this->neutrons = 153;
		this->electrons = 98;
		this->valence_e = 12;
		this->Name = "Californium";
		this->Symbol= "Cf";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 98;
		this->atomic_radii = 1.70;
	}
	else if (number == 99)
	{
		this->atomic_weight = 252.0;
		this->oxidation_state = 3;
		this->protons = 99;
		this->neutrons = 153;
		this->electrons = 99;
		this->valence_e = 13;
		this->Name = "Einsteinium";
		this->Symbol= "Es";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 99;
		this->atomic_radii = 1.70;
	}
	else if (number == 100)
	{
		this->atomic_weight = 257.0;
		this->oxidation_state = 3;
		this->protons = 100;
		this->neutrons = 157;
		this->electrons = 100;
		this->valence_e = 14;
		this->Name = "Fermium";
		this->Symbol= "Fm";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 100;
		this->atomic_radii = 1.70;
	}
	else if (number == 101)
	{
		this->atomic_weight = 258.0;
		this->oxidation_state = 3;
		this->protons = 101;
		this->neutrons = 157;
		this->electrons = 101;
		this->valence_e = 15;
		this->Name = "Mendelevium";
		this->Symbol= "Md";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 101;
		this->atomic_radii = 1.70;
	}
	else if (number == 102)
	{
		this->atomic_weight = 259.0;
		this->oxidation_state = 2;
		this->protons = 102;
		this->neutrons = 157;
		this->electrons = 102;
		this->valence_e = 16;
		this->Name = "Nobelium";
		this->Symbol= "No";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 102;
		this->atomic_radii = 1.70;
	}
	else if (number == 103)
	{
		this->atomic_weight = 266.0;
		this->oxidation_state = 3;
		this->protons = 103;
		this->neutrons = 159;
		this->electrons = 103;
		this->valence_e = 3;
		this->Name = "Lawrencium";
		this->Symbol= "Lr";
		this->Category = "Actinide";
		this->NaturalState = "Solid";
		this->atomic_number = 103;
		this->atomic_radii = 1.70;
	}
	else if (number == 104)
	{
		this->atomic_weight = 267.0;
		this->oxidation_state = 4;
		this->protons = 104;
		this->neutrons = 157;
		this->electrons = 104;
		this->valence_e = 4;
		this->Name = "Rutherfordium";
		this->Symbol= "Rf";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 104;
		this->atomic_radii = 1.57;
	}
	else if (number == 105)
	{
		this->atomic_weight = 268.0;
		this->oxidation_state = 5;
		this->protons = 105;
		this->neutrons = 157;
		this->electrons = 105;
		this->valence_e = 5;
		this->Name = "Dubnium";
		this->Symbol= "Db";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 105;
		this->atomic_radii = 1.49;
	}
	else if (number == 106)
	{
		this->atomic_weight = 269.0;
		this->oxidation_state = 6;
		this->protons = 106;
		this->neutrons = 157;
		this->electrons = 106;
		this->valence_e = 6;
		this->Name = "Seaborgium";
		this->Symbol= "Sg";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 106;
		this->atomic_radii = 1.43;
	}
	else if (number == 107)
	{
		this->atomic_weight = 270.0;
		this->oxidation_state = 7;
		this->protons = 107;
		this->neutrons = 155;
		this->electrons = 107;
		this->valence_e = 7;
		this->Name = "Bohrium";
		this->Symbol= "Bh";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 107;
		this->atomic_radii = 1.41;
	}
	else if (number == 108)
	{
		this->atomic_weight = 269.0;
		this->oxidation_state = 8;
		this->protons = 108;
		this->neutrons = 157;
		this->electrons = 108;
		this->valence_e = 8;
		this->Name = "Hassium";
		this->Symbol= "Hs";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 108;
		this->atomic_radii = 1.34;
	}
	else if (number == 109)
	{
		this->atomic_weight = 278.0;
		this->oxidation_state = 6;
		this->protons = 109;
		this->neutrons = 157;
		this->electrons = 109;
		this->valence_e = 9;
		this->Name = "Meitnerium";
		this->Symbol= "Mt";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 109;
		this->atomic_radii = 1.29;
	}
	else if (number == 110)
	{
		this->atomic_weight = 281.0;
		this->oxidation_state = 8;
		this->protons = 110;
		this->neutrons = 171;
		this->electrons = 110;
		this->valence_e = 10;
		this->Name = "Darmstadium";
		this->Symbol= "Ds";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 110;
		this->atomic_radii = 1.28;
	}
	else if (number == 111)
	{
		this->atomic_weight = 281.0;
		this->oxidation_state = 3;
		this->protons = 111;
		this->neutrons = 170;
		this->electrons = 111;
		this->valence_e = 11;
		this->Name = "Roentgenium";
		this->Symbol= "Rg";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 111;
		this->atomic_radii = 1.21;
	}
	else if (number == 112)
	{
		this->atomic_weight = 285.0;
		this->oxidation_state = 4;
		this->protons = 112;
		this->neutrons = 173;
		this->electrons = 112;
		this->valence_e = 12;
		this->Name = "Copernicium";
		this->Symbol= "Cn";
		this->Category = "Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 112;
		this->atomic_radii = 1.22;
	}
	else if (number == 113)
	{
		this->atomic_weight = 286.0;
		this->oxidation_state = 1;
		this->protons = 113;
		this->neutrons = 173;
		this->electrons = 113;
		this->valence_e = 3;
		this->Name = "Ununtrium";
		this->Symbol= "Uut";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 113;
		this->atomic_radii = 1.76;
	}
	else if (number == 114)
	{
		this->atomic_weight = 289.0;
		this->oxidation_state = 2;
		this->protons = 114;
		this->neutrons = 175;
		this->electrons = 114;
		this->valence_e = 4;
		this->Name = "Flerovium";
		this->Symbol= "Fl";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 114;
		this->atomic_radii = 1.74;
	}
	else if (number == 115)
	{
		this->atomic_weight = 289.0;
		this->oxidation_state = 1;
		this->protons = 115;
		this->neutrons = 174;
		this->electrons = 115;
		this->valence_e = 5;
		this->Name = "Ununpentium";
		this->Symbol= "Uup";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 115;
		this->atomic_radii = 1.57;
	}
	else if (number == 116)
	{
		this->atomic_weight = 293.0;
		this->oxidation_state = 2;
		this->protons = 116;
		this->neutrons = 177;
		this->electrons = 116;
		this->valence_e = 6;
		this->Name = "Livermorium";
		this->Symbol= "Lv";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 116;
		this->atomic_radii = 1.64;
	}
	else if (number == 117)
	{
		this->atomic_weight = 294.0;
		this->oxidation_state = 1;
		this->protons = 117;
		this->neutrons = 177;
		this->electrons = 117;
		this->valence_e = 7;
		this->Name = "Ununseptium";
		this->Symbol= "Uus";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 117;
		this->atomic_radii = 1.57;
	}
	else if (number == 118)
	{
		this->atomic_weight = 294.0;
		this->oxidation_state = 0;
		this->protons = 118;
		this->neutrons = 176;
		this->electrons = 118;
		this->valence_e = 8;
		this->Name = "Ununoctium";
		this->Symbol= "Uuo";
		this->Category = "Nobel Gas";
		this->NaturalState = "Solid";
		this->atomic_number = 118;
		this->atomic_radii = 1.57;
	}
	else
	{
		std::cout << "Given Number: " << number << "\tValid Option: 1-118" << std::endl;
		mError(invalid_atom);
		this->atomic_weight = 0.0;
		this->oxidation_state = 0;
		this->protons = 0;
		this->neutrons = 0;
		this->electrons = 0;
		this->Name = "No Name";
		this->Symbol= "N/A";
		this->Category = "N/A";
		this->NaturalState = "N/A";
		this->atomic_number = 0;
	}
}

//Edit the atomic weight of the atom
void Atom::editAtomicWeight(double AW)
{
	if (AW <= 0.0) {mError(negative_mass); return;}
	this->atomic_weight = AW;
}

//Edit the oxidation state of the atom
void Atom::editOxidationState(int state)
{
	this->oxidation_state = state;
}

//Edit the number of protons in the atom
void Atom::editProtons(int proton)
{
	if (proton <= 0) {mError(invalid_proton); return;}
	this->protons = proton;
}

//Edit the number of neutrons in the atom
void Atom::editNeutrons(int neutron)
{
	if (neutron < 0) {mError(invalid_neutron); return;}
	this->neutrons = neutron;
}

//Edit the number of electrons in the atom
void Atom::editElectrons(int electron)
{
	if (electron < 0) {mError(invalid_electron); return;}
	this->electrons = electron;
}

//Edit the number of electrons in the valence shells
void Atom::editValence(int val)
{
	//Min Valence = 0; Max Valence = 12
	if (this->Category == "Transition Metal")
	{
		if (val > 12 || val < 0)
		{
			mError(invalid_valence);
			return;
		}
	}
	//Min Valence = 0; Max Valence = 16
	else if (this->Category == "Actinide" || this->Category == "Lanthanide")
	{
		if (val > 16 || val < 0)
		{
			mError(invalid_valence);
			return;
		}
	}
	//Min Valence = 0; Max Valence = 2
	else if (this->atomic_number == 1 || this->atomic_number == 2)
	{
		if (val > 2 || val < 0)
		{
			mError(invalid_valence);
			return;
		}
	}
	//Min Valence = 0; Max Valence = 8
	else
	{
		if (val > 8 || val < 0)
		{
			mError(invalid_valence);
			return;
		}
	}
	this->valence_e = val;
}

//Edit the atomic radii
void Atom::editRadii(double r)
{
	if (r <= 0)
		this->atomic_radii = 1.2;
	else
		this->atomic_radii = r;
}

//Remove a proton
void Atom::removeProton()
{
	if (this->protons <= 1)
	{
		mError(invalid_proton);
		return;
	}
	int neutrons = this->neutrons;
	int electrons = this->electrons;
	double atom_weight = this->atomic_weight - 1.0;
	int valence = this->valence_e;
	this->Register(this->atomic_number-1);
	this->editNeutrons(neutrons);
	this->editElectrons(electrons);
	this->editAtomicWeight(atom_weight);
	this->editValence(valence);
}

//Remove a neutron
void Atom::removeNeutron()
{
	if (this->neutrons < 1)
	{
		mError(invalid_neutron);
		return;
	}
	this->neutrons--;
	double atom_weight = this->atomic_weight - 1.0;
	this->editAtomicWeight(atom_weight);
}

//Remove an electron
void Atom::removeElectron()
{
	if (this->electrons < 1)
	{
		mError(invalid_electron);
		return;
	}
	this->electrons--;
	if (this->valence_e == 1)
		this->valence_e = 8;
	else
		this->valence_e--;
}

//Return the atomic weight
double Atom::AtomicWeight()
{
	return this->atomic_weight;
}

//Return the oxidation state
int Atom::OxidationState()
{
	return this->oxidation_state;
}

//Return the number of protons
int Atom::Protons()
{
	return this->protons;
}

//Return the number of neutrons
int Atom::Neutrons()
{
	return this->neutrons;
}

//Return the number of electrons
int Atom::Electrons()
{
	return this->electrons;
}

//Return the number of bonding electrons
int Atom::BondingElectrons()
{
	return this->valence_e;
}

//Return the van der Waals radii
double Atom::AtomicRadii()
{
	return this->atomic_radii;
}

//Return the name of the atom
std::string Atom::AtomName()
{
	return this->Name;
}

//Return the atom symbol
std::string Atom::AtomSymbol()
{
	return this->Symbol;
}

//Return the atom category
std::string Atom::AtomCategory()
{
	return this->Category;
}

//Return the atom state
std::string Atom::AtomState()
{
	return this->NaturalState;
}

//Return teh atomic number
int Atom::AtomicNumber()
{
	return this->atomic_number;
}

//Display atom info to console
void Atom::DisplayInfo()
{
	std::cout << this->Name << " (" << this->Symbol << "): \n";
	std::cout << "Atomic Number: " << this->atomic_number << "\tAtomic Weight (g/mol): " << this->atomic_weight << std::endl;
	std::cout << "Natural State: " << this->NaturalState << "\tCategory: " << this->Category << std::endl;
	std::cout << "Current protons = " << this->protons << ", neutrons = " << this->neutrons << ", electrons = " << this->electrons << std::endl;
	std::cout << "Current Oxidation State: ";
	if (this->oxidation_state > 0)
		std::cout << "+" << this->oxidation_state;
	else
		std::cout << this->oxidation_state;
	std::cout << "\tValence Electrons: " << this->valence_e << std::endl;
	std::cout << "van der Waals Radius (Angstroms): " << this->atomic_radii << std::endl;
	std::cout << std::endl;
}


//Default Constructor for Full Digital Periodic Table
PeriodicTable::PeriodicTable()
:
Table(118),
number_elements(118)
{
	for (int i=0; i<number_elements; i++)
	{
		Table[i].Register(i+1);
	}
}

//Default Destructor
PeriodicTable::~PeriodicTable()
{
	
}

//Construct a Partial Table from an array of atomic numbers
PeriodicTable::PeriodicTable(int *n, int N)
:
Table(N),
number_elements(N)
{
	for (int i=0; i<N; i++)
		Table[i].Register(n[i]);
}

//Construct a Partial Table from a vector of atomic symbols
PeriodicTable::PeriodicTable(std::vector<std::string> &Symbol)
:
Table(Symbol.size()),
number_elements((int)Symbol.size())
{
	for (int i=0; i<number_elements; i++)
		Table[i].Register(Symbol[i]);
}

//Construct a Partial Table from a vector of atomic numbers
PeriodicTable::PeriodicTable(std::vector<int> &n)
:
Table(n.size()),
number_elements((int)n.size())
{
	for (int i=0; i<number_elements; i++)
		Table[i].Register(n[i]);
}

//Display the periodic table
void PeriodicTable::DisplayTable()
{
	//Display the full Periodic Table
	if (this->number_elements == 118)
	{
		std::cout << "----------------------------- Full Periodic Table -----------------------------\n";
		std::cout << "\t" << Table[0].AtomSymbol();
		for (int i=0; i<18; i++)
			std::cout << "\t";
		std::cout << Table[1].AtomSymbol() << std::endl;
		std::cout << "\t" << Table[2].AtomSymbol() << "\t" << Table[3].AtomSymbol();
		for (int i=0; i<11; i++)
			std::cout << "\t";
		for (int i=0; i<6; i++)
			std::cout << "\t" << Table[i+4].AtomSymbol();
		std::cout << std::endl;
		std::cout << "\t" << Table[10].AtomSymbol() << "\t" << Table[11].AtomSymbol();
		for (int i=0; i<11; i++)
			std::cout << "\t";
		for (int i=0; i<6; i++)
			std::cout << "\t" << Table[i+12].AtomSymbol();
		std::cout << std::endl;
		std::cout << "\t" << Table[18].AtomSymbol() << "\t" << Table[19].AtomSymbol() << "\t";
		for (int i=0; i<16; i++)
			std::cout << "\t" << Table[i+20].AtomSymbol();
		std::cout << std::endl;
		std::cout << "\t" << Table[36].AtomSymbol() << "\t" << Table[37].AtomSymbol() << "\t";
		for (int i=0; i<16; i++)
			std::cout << "\t" << Table[i+38].AtomSymbol();
		std::cout << std::endl;
		std::cout << "\t" << Table[54].AtomSymbol() << "\t" << Table[55].AtomSymbol() << "\t*";
		for (int i=0; i<16; i++)
			std::cout << "\t" << Table[i+56+14].AtomSymbol();
		std::cout << std::endl;
		std::cout << "\t" << Table[86].AtomSymbol() << "\t" << Table[87].AtomSymbol() << "\t**";
		for (int i=0; i<16; i++)
			std::cout << "\t" << Table[i+88+14].AtomSymbol();
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "\t\t*";
		for (int i=0; i<14; i++)
			std::cout << "\t" << Table[i+56].AtomSymbol();
		std::cout << std::endl;
		std::cout << "\t\t**";
		for (int i=0; i<14; i++)
			std::cout << "\t" << Table[i+88].AtomSymbol();
		std::cout << std::endl;
	}
	//Partial Table listed
	else
	{
		std::cout << "---- Partial List of Registered Atoms -----\n";
		std::cout << "\tIndex\tSymbol\tAtomic Number\n";
		for (int i=0; i<this->number_elements; i++)
		{
			std::cout << "\t" << i << "\t" << Table[i].AtomSymbol() << "\t" << Table[i].AtomicNumber() << std::endl;
		}
	}
	
	std::cout << std::endl;
}

//Function to run tests on the classes
int EEL_TESTS()
{
	int success = 0;
	double time = clock();
	
	//--------------------Atom Tests---------------------------
	std::string name1 = "Hydrogen";
	Atom particle1(name1), particle2("Ununoctium");
	particle1.DisplayInfo();
	particle2.DisplayInfo();
	
	Atom regTest1, regTest2;
	std::string sym1 = "K";
	regTest2.Register(sym1);
	regTest1.Register("Ne");
	regTest1.DisplayInfo();
	regTest2.DisplayInfo();
	
	Atom invalid;
	invalid.Register("n");
	invalid.Register("Na");
	invalid.editValence(9);
	
	invalid.DisplayInfo();
	invalid.removeProton();
	invalid.removeNeutron();
	invalid.removeElectron();
	invalid.DisplayInfo();
	
	std::cout << particle1.AtomName() << std::endl << std::endl;
	//----------------End Atom Tests----------------------------
	
	//----------------Periodic Table Tests----------------------
	PeriodicTable pt;
	pt.DisplayTable();
	
	int N = 5;
	int nums[N];
	nums[0] = 2;
	nums[1] = 12;
	nums[2] = 1;
	nums[3] = 100;
	nums[4] = 35;
	
	PeriodicTable partial(nums,N);
	partial.DisplayTable();
	
	std::vector<std::string> Sym;
	Sym.resize(N);
	Sym[0] = "Ne";
	Sym[1] = "Ar";
	Sym[2] = "U";
	Sym[3] = "H";
	Sym[4] = "Be";
	
	PeriodicTable part2(Sym);
	part2.DisplayTable();
	
	std::vector<int> NUMS;
	NUMS.resize(N);
	for (int i=0; i<N; i++)
		NUMS[i] = nums[i];
	PeriodicTable part3(NUMS);
	part3.DisplayTable();
	//----------------End Periodic Table Tests------------------
	
	//-------------- Test Alex's Contribution to EEL -----------------------------------
	Atom a(1);
	a.DisplayInfo();
	
	for (int i=0; i<10; i++)
	{
		//Note: Because of scope, these atoms only exist inside this loop
		Atom b(i+1);
		b.DisplayInfo();
	}
	for (int i=22; i<54; i++)
	{
		//Note: Because of scope, these atoms only exist inside this loop
		Atom b(i+1);
		b.DisplayInfo();
	}
	
	//All tests passed
	//-------------- Alex Tests --------------------------------------------------------
	
	time = clock() - time;
	std::cout << "Runtime (s): " << (time/CLOCKS_PER_SEC) << std::endl;
	return success;
}
