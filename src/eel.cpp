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
atomic_radii(0),
melting_point(0),
boiling_point(0),
thermal_x_sec(0)
{
}

//Default Destructor
Atom::~Atom()
{
}

//Register an instance of an Atom object by it's symbol
void Atom::Register(std::string Symbol)
{
	if (Symbol.compare("n") == 0)
	{
		this->atomic_weight = 1.0081;
		this->oxidation_state = 0;
		this->protons = 0;
		this->neutrons = 1;
		this->electrons = 0;
		this->valence_e = 0;
		this->Name = "neutron";
		this->Symbol= "n";
		this->Category = "Diatomic Non-metal";
		this->NaturalState = "Gas";
		this->atomic_number = 0;
		this->atomic_radii = 1.2;
		this->melting_point = 13.99;
		this->boiling_point = 20.271;
		this->thermal_x_sec = 0.0;
	}
	else if (Symbol.compare("H") == 0)
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
		this->melting_point = 13.99;
		this->boiling_point = 20.271;
		this->thermal_x_sec = 0.3326;
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
		this->melting_point = 0.95;
		this->boiling_point = 4.222;
		this->thermal_x_sec = 0.007;
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
		this->melting_point = 453.65;
		this->boiling_point = 1603;
		this->thermal_x_sec = 70.5;
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
		this->melting_point = 1560;
		this->boiling_point = 2742;
		this->thermal_x_sec = 0.0092;
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
		this->melting_point = 2349;
		this->boiling_point = 4200;
		this->thermal_x_sec = 767.0;
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
		this->melting_point = 3915;
		this->boiling_point = 3915;
		this->thermal_x_sec = 0.0035;
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
		this->melting_point = 63.15;
		this->boiling_point = 77.355;
		this->thermal_x_sec = 1.91;
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
		this->melting_point = 54.36;
		this->boiling_point = 90.188;
		this->thermal_x_sec = 0.00019;
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
		this->melting_point = 53.48;
		this->boiling_point = 85.03;
		this->thermal_x_sec = 0.0096;
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
		this->melting_point = 24.56;
		this->boiling_point = 27.104;
		this->thermal_x_sec = 0.04;
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
		this->melting_point = 370.944;
		this->boiling_point = 1156.09;
		this->thermal_x_sec = 0.53;
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
		this->melting_point = 923;
		this->boiling_point = 1363;
		this->thermal_x_sec = 0.063;
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
		this->melting_point = 933.47;
		this->boiling_point = 2743;
		this->thermal_x_sec = 0.232;
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
		this->melting_point = 1687;
		this->boiling_point = 3538;
		this->thermal_x_sec = 0.171;
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
		this->melting_point = 317.3;
		this->boiling_point = 553.7;
		this->thermal_x_sec = 0.172;
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
		this->melting_point = 388.36;
		this->boiling_point = 717.8;
		this->thermal_x_sec = 0.53;
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
		this->melting_point = 171.6;
		this->boiling_point = 239.11;
		this->thermal_x_sec = 35.5;
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
		this->melting_point = 83.81;
		this->boiling_point = 87.302;
		this->thermal_x_sec = 0.675;
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
		this->melting_point = 336.7;
		this->boiling_point = 1032;
		this->thermal_x_sec = 2.1;
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
		this->melting_point = 1115;
		this->boiling_point = 1757;
		this->thermal_x_sec = 0.43;
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
		this->melting_point = 1814;
		this->boiling_point = 3109;
		this->thermal_x_sec = 27.2;
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
		this->melting_point = 1941;
		this->boiling_point = 3034;
		this->thermal_x_sec = 6.09;
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
		this->melting_point = 2183;
		this->boiling_point = 3680;
		this->thermal_x_sec = 5.08;
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
		this->melting_point = 2180;
		this->boiling_point = 2944;
		this->thermal_x_sec = 3.1;
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
		this->melting_point = 1519;
		this->boiling_point = 2334;
		this->thermal_x_sec = 13.3;
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
		this->melting_point = 1811;
		this->boiling_point = 3134;
		this->thermal_x_sec = 2.56;
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
		this->melting_point = 1768;
		this->boiling_point = 2723;
		this->thermal_x_sec = 37.2;
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
		this->melting_point = 1728;
		this->boiling_point = 3003;
		this->thermal_x_sec = 4.49;
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
		this->melting_point = 1357.77;
		this->boiling_point = 2835;
		this->thermal_x_sec = 3.78;
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
		this->melting_point = 692.68;
		this->boiling_point = 1180;
		this->thermal_x_sec = 1.11;
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
		this->melting_point = 302.9146;
		this->boiling_point = 2673;
		this->thermal_x_sec = 2.9;
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
		this->melting_point = 1211.40;
		this->boiling_point = 3106;
		this->thermal_x_sec = 2.2;
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
		this->melting_point = 887;
		this->boiling_point = 887;
		this->thermal_x_sec = 4.3;
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
		this->melting_point = 494;
		this->boiling_point = 958;
		this->thermal_x_sec = 11.7;
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
		this->melting_point = 265.8;
		this->boiling_point = 332;
		this->thermal_x_sec = 6.8;
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
		this->melting_point = 115.78;
		this->boiling_point = 119.93;
		this->thermal_x_sec = 25.0;
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
		this->melting_point = 312.45;
		this->boiling_point = 961;
		this->thermal_x_sec = 0.38;
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
		this->melting_point = 1050;
		this->boiling_point = 1650;
		this->thermal_x_sec = 1.28;
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
		this->melting_point = 1799;
		this->boiling_point = 3203;
		this->thermal_x_sec = 1.28;
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
		this->melting_point = 2128;
		this->boiling_point = 4650;
		this->thermal_x_sec = 0.184;
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
		this->melting_point = 2750;
		this->boiling_point = 5017;
		this->thermal_x_sec = 1.15;
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
		this->melting_point = 2896;
		this->boiling_point = 4912;
		this->thermal_x_sec = 2.6;
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
		this->melting_point = 2430;
		this->boiling_point = 4538;
		this->thermal_x_sec = 20.0;
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
		this->melting_point = 2607;
		this->boiling_point = 4423;
		this->thermal_x_sec = 2.56;
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
		this->melting_point = 2237;
		this->boiling_point = 3968;
		this->thermal_x_sec = 144.8;
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
		this->melting_point = 1828.05;
		this->boiling_point = 3236;
		this->thermal_x_sec = 6.9;
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
		this->melting_point = 1234.93;
		this->boiling_point = 2435;
		this->thermal_x_sec = 63.6;
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
		this->melting_point = 594.22;
		this->boiling_point = 1040;
		this->thermal_x_sec = 2450.0;
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
		this->melting_point = 429.7485;
		this->boiling_point = 2345;
		this->thermal_x_sec = 194;
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
		this->melting_point = 505.08;
		this->boiling_point = 2875;
		this->thermal_x_sec = 0.626;
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
		this->melting_point = 903.78;
		this->boiling_point = 1908;
		this->thermal_x_sec = 4.91;
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
		this->melting_point = 722.66;
		this->boiling_point = 1261;
		this->thermal_x_sec = 4.7;
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
		this->melting_point = 386.85;
		this->boiling_point = 457.4;
		this->thermal_x_sec = 6.2;
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
		this->melting_point = 161.40;
		this->boiling_point = 165.051;
		this->thermal_x_sec = 23.9;
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
		this->melting_point = 301.7;
		this->boiling_point = 944;
		this->thermal_x_sec = 29;
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
		this->melting_point = 1000;
		this->boiling_point = 2118;
		this->thermal_x_sec = 1.3;
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
		this->melting_point = 1193;
		this->boiling_point = 3737;
		this->thermal_x_sec = 8.98;
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
		this->melting_point = 1068;
		this->boiling_point = 3716;
		this->thermal_x_sec = 0.6;
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
		this->melting_point = 1208;
		this->boiling_point = 3403;
		this->thermal_x_sec = 11.5;
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
		this->melting_point = 1297;
		this->boiling_point = 3347;
		this->thermal_x_sec = 49;
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
		this->melting_point = 1315;
		this->boiling_point = 3273;
		this->thermal_x_sec = 168.4;
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
		this->melting_point = 1345;
		this->boiling_point = 2173;
		this->thermal_x_sec = 5922;
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
		this->melting_point = 1099;
		this->boiling_point = 1802;
		this->thermal_x_sec = 4600;
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
		this->melting_point = 1585;
		this->boiling_point = 3273;
		this->thermal_x_sec = 49000;
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
		this->melting_point = 1629;
		this->boiling_point = 3396;
		this->thermal_x_sec = 23.4;
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
		this->melting_point = 1680;
		this->boiling_point = 2840;
		this->thermal_x_sec = 920;
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
		this->melting_point = 1734;
		this->boiling_point = 2873;
		this->thermal_x_sec = 65;
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
		this->melting_point = 1802;
		this->boiling_point = 3141;
		this->thermal_x_sec = 160;
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
		this->melting_point = 1818;
		this->boiling_point = 2223;
		this->thermal_x_sec = 100;
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
		this->melting_point = 1097;
		this->boiling_point = 1469;
		this->thermal_x_sec = 34.8;
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
		this->melting_point = 1925;
		this->boiling_point = 3675;
		this->thermal_x_sec = 84;
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
		this->melting_point = 2506;
		this->boiling_point = 4876;
		this->thermal_x_sec = 104;
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
		this->melting_point = 3290;
		this->boiling_point = 5463;
		this->thermal_x_sec = 20.6;
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
		this->melting_point = 3695;
		this->boiling_point = 6203;
		this->thermal_x_sec = 18.3;
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
		this->melting_point = 3459;
		this->boiling_point = 5903;
		this->thermal_x_sec = 89.7;
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
		this->melting_point = 3306;
		this->boiling_point = 5285;
		this->thermal_x_sec = 15;
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
		this->melting_point = 2719;
		this->boiling_point = 4403;
		this->thermal_x_sec = 425;
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
		this->melting_point = 2041.4;
		this->boiling_point = 4098;
		this->thermal_x_sec = 0.96;
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
		this->melting_point = 1337.33;
		this->boiling_point = 3243;
		this->thermal_x_sec = 98.7;
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
		this->melting_point = 234.321;
		this->boiling_point = 629.88;
		this->thermal_x_sec = 374;
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
		this->melting_point = 577;
		this->boiling_point = 1746;
		this->thermal_x_sec = 3.43;
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
		this->melting_point = 600.61;
		this->boiling_point = 2022;
		this->thermal_x_sec = 0.171;
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
		this->melting_point = 544.7;
		this->boiling_point = 1837;
		this->thermal_x_sec = 0.034;
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
		this->melting_point = 527;
		this->boiling_point = 1235;
		this->thermal_x_sec = 0.03;
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
		this->melting_point = 575;
		this->boiling_point = 610;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 202;
		this->boiling_point = 211.5;
		this->thermal_x_sec = 0.72;
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
		this->melting_point = 300;
		this->boiling_point = 950;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 973;
		this->boiling_point = 2010;
		this->thermal_x_sec = 12.8;
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
		this->melting_point = 1500;
		this->boiling_point = 3500;
		this->thermal_x_sec = 515;
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
		this->melting_point = 2023;
		this->boiling_point = 5061;
		this->thermal_x_sec = 7.37;
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
		this->melting_point = 1841;
		this->boiling_point = 4300;
		this->thermal_x_sec = 200.6;
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
		this->melting_point = 1405.3;
		this->boiling_point = 4404;
		this->thermal_x_sec = 7.57;
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
		this->melting_point = 912;
		this->boiling_point = 4447;
		this->thermal_x_sec = 180;
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
		this->melting_point = 912.5;
		this->boiling_point = 3505;
		this->thermal_x_sec = 1017.3;
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
		this->melting_point = 1449;
		this->boiling_point = 2880;
		this->thermal_x_sec = 75.3;
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
		this->melting_point = 1613;
		this->boiling_point = 3383;
		this->thermal_x_sec = 79;
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
		this->melting_point = 1259;
		this->boiling_point = 2900;
		this->thermal_x_sec = 710;
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
		this->melting_point = 1173;
		this->boiling_point = 1743;
		this->thermal_x_sec = 2900;
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
		this->melting_point = 1133;
		this->boiling_point = 1269;
		this->thermal_x_sec = 160;
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
		this->melting_point = 1800;
		this->boiling_point = 3000;
		this->thermal_x_sec = 5800;
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
		this->melting_point = 1100;
		this->boiling_point = 2500;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 1100;
		this->boiling_point = 2500;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 1900;
		this->boiling_point = 3100;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 2400;
		this->boiling_point = 5800;
		this->thermal_x_sec = 0.0;
	}
	else if (Symbol.compare("Nh") == 0)
	{
		this->atomic_weight = 286.0;
		this->oxidation_state = 1;
		this->protons = 113;
		this->neutrons = 173;
		this->electrons = 113;
		this->valence_e = 3;
		this->Name = "Nihonium";
		this->Symbol= "Nh";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 113;
		this->atomic_radii = 1.76;
		this->melting_point = 700;
		this->boiling_point = 1430;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 200;
		this->boiling_point = 210;
		this->thermal_x_sec = 0.0;
	}
	else if (Symbol.compare("Mc") == 0)
	{
		this->atomic_weight = 289.0;
		this->oxidation_state = 1;
		this->protons = 115;
		this->neutrons = 174;
		this->electrons = 115;
		this->valence_e = 5;
		this->Name = "Moscovium";
		this->Symbol= "Mc";
		this->Category = "Post-Transition Metal";
		this->NaturalState = "Solid";
		this->atomic_number = 115;
		this->atomic_radii = 1.57;
		this->melting_point = 670;
		this->boiling_point = 1400;
		this->thermal_x_sec = 0.0;
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
		this->melting_point = 740;
		this->boiling_point = 1085;
		this->thermal_x_sec = 0.0;
	}
	else if (Symbol.compare("Ts") == 0)
	{
		this->atomic_weight = 294.0;
		this->oxidation_state = 1;
		this->protons = 117;
		this->neutrons = 177;
		this->electrons = 117;
		this->valence_e = 7;
		this->Name = "Tennessine";
		this->Symbol= "Ts";
		this->Category = "Metalloid";
		this->NaturalState = "Solid";
		this->atomic_number = 117;
		this->atomic_radii = 1.57;
		this->melting_point = 723;
		this->boiling_point = 883;
		this->thermal_x_sec = 0.0;
	}
	else if (Symbol.compare("Og") == 0)
	{
		this->atomic_weight = 294.0;
		this->oxidation_state = 0;
		this->protons = 118;
		this->neutrons = 176;
		this->electrons = 118;
		this->valence_e = 8;
		this->Name = "Oganesson";
		this->Symbol= "Og";
		this->Category = "Nobel Gas";
		this->NaturalState = "Solid";
		this->atomic_number = 118;
		this->atomic_radii = 1.57;
		this->melting_point = 300;
		this->boiling_point = 350;
		this->thermal_x_sec = 0.0;
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
		this->atomic_radii = 1.20;
	}
}

//Register a instance of an Atom object by it's symbol
void Atom::Register(int number)
{
	if (number == 0)
	{
		this->Register("n");
	}
	else if (number == 1)
	{
		this->Register("H");
	}
	else if (number == 2)
	{
		this->Register("He");
	}
	else if (number == 3)
	{
		this->Register("Li");
	}
	else if (number == 4)
	{
		this->Register("Be");
	}
	else if (number == 5)
	{
		this->Register("B");
	}
	else if (number == 6)
	{
		this->Register("C");
	}
	else if (number == 7)
	{
		this->Register("N");
	}
	else if (number == 8)
	{
		this->Register("O");
	}
	else if (number == 9)
	{
		this->Register("F");
	}
	else if (number == 10)
	{
		this->Register("Ne");
	}
	else if (number == 11)
	{
		this->Register("Na");
	}
	else if (number == 12)
	{
		this->Register("Mg");
	}
	else if (number == 13)
	{
		this->Register("Al");
	}
	else if (number == 14)
	{
		this->Register("Si");
	}
	else if (number == 15)
	{
		this->Register("P");
	}
	else if (number == 16)
	{
		this->Register("S");
	}
	else if (number == 17)
	{
		this->Register("Cl");
	}
	else if (number == 18)
	{
		this->Register("Ar");
	}
	else if (number == 19)
	{
		this->Register("K");
	}
	else if (number == 20)
	{
		this->Register("Ca");
	}
	else if (number == 21)
	{
		this->Register("Sc");
	}
	else if (number == 22)
	{
		this->Register("Ti");
	}
	else if (number == 23)
	{
		this->Register("V");
	}
	else if (number == 24)
	{
		this->Register("Cr");
	}
	else if (number == 25)
	{
		this->Register("Mn");
	}
	else if (number == 26)
	{
		this->Register("Fe");
	}
	else if (number == 27)
	{
		this->Register("Co");
	}
	else if (number == 28)
	{
		this->Register("Ni");
	}
	else if (number == 29)
	{
		this->Register("Cu");
	}
	else if (number == 30)
	{
		this->Register("Zn");
	}
	else if (number == 31)
	{
		this->Register("Ga");
	}
	else if (number == 32)
	{
		this->Register("Ge");
	}
	else if (number == 33)
	{
		this->Register("As");
	}
	else if (number == 34)
	{
		this->Register("Se");
	}
	else if (number == 35)
	{
		this->Register("Br");
	}
	else if (number == 36)
	{
		this->Register("Kr");
	}
	else if (number == 37)
	{
		this->Register("Rb");
	}
	else if (number == 38)
	{
		this->Register("Sr");
	}
	else if (number == 39)
	{
		this->Register("Y");
	}
	else if (number == 40)
	{
		this->Register("Zr");
	}
	else if (number == 41)
	{
		this->Register("Nb");
	}
	else if (number == 42)
	{
		this->Register("Mo");
	}
	else if (number == 43)
	{
		this->Register("Tc");
	}
	else if (number == 44)
	{
		this->Register("Ru");
	}
	else if (number == 45)
	{
		this->Register("Rh");
	}
	else if (number == 46)
	{
		this->Register("Pd");
	}
	else if (number == 47)
	{
		this->Register("Ag");
	}
	else if (number == 48)
	{
		this->Register("Cd");
	}
	else if (number == 49)
	{
		this->Register("In");
	}
	else if (number == 50)
	{
		this->Register("Sn");
	}
	else if (number == 51)
	{
		this->Register("Sb");
	}
	else if (number == 52)
	{
		this->Register("Te");
	}
	else if (number == 53)
	{
		this->Register("I");
	}
	else if (number == 54)
	{
		this->Register("Xe");
	}
	else if (number == 55)
	{
		this->Register("Cs");
	}
	else if (number == 56)
	{
		this->Register("Ba");
	}
	else if (number == 57)
	{
		this->Register("La");
	}
	else if (number == 58)
	{
		this->Register("Ce");
	}
	else if (number == 59)
	{
		this->Register("Pr");
	}
	else if (number == 60)
	{
		this->Register("Nd");
	}
	else if (number == 61)
	{
		this->Register("Pm");
	}
	else if (number == 62)
	{
		this->Register("Sm");
	}
	else if (number == 63)
	{
		this->Register("Eu");
	}
	else if (number == 64)
	{
		this->Register("Gd");
	}
	else if (number == 65)
	{
		this->Register("Tb");
	}
	else if (number == 66)
	{
		this->Register("Dy");
	}
	else if (number == 67)
	{
		this->Register("Ho");
	}
	else if (number == 68)
	{
		this->Register("Er");
	}
	else if (number == 69)
	{
		this->Register("Tm");
	}
	else if (number == 70)
	{
		this->Register("Yb");
	}
	else if (number == 71)
	{
		this->Register("Lu");
	}
	else if (number == 72)
	{
		this->Register("Hf");
	}
	else if (number == 73)
	{
		this->Register("Ta");
	}
	else if (number == 74)
	{
		this->Register("W");
	}
	else if (number == 75)
	{
		this->Register("Re");
	}
	else if (number == 76)
	{
		this->Register("Os");
	}
	else if (number == 77)
	{
		this->Register("Ir");
	}
	else if (number == 78)
	{
		this->Register("Pt");
	}
	else if (number == 79)
	{
		this->Register("Au");
	}
	else if (number == 80)
	{
		this->Register("Hg");
	}
	else if (number == 81)
	{
		this->Register("Tl");
	}
	else if (number == 82)
	{
		this->Register("Pb");
	}
	else if (number == 83)
	{
		this->Register("Bi");
	}
	else if (number == 84)
	{
		this->Register("Po");
	}
	else if (number == 85)
	{
		this->Register("At");
	}
	else if (number == 86)
	{
		this->Register("Rn");
	}
	else if (number == 87)
	{
		this->Register("Fr");
	}
	else if (number == 88)
	{
		this->Register("Ra");
	}
	else if (number == 89)
	{
		this->Register("Ac");
	}
	else if (number == 90)
	{
		this->Register("Th");
	}
	else if (number == 91)
	{
		this->Register("Pa");
	}
	else if (number == 92)
	{
		this->Register("U");
	}
	else if (number == 93)
	{
		this->Register("Np");
	}
	else if (number == 94)
	{
		this->Register("Pu");
	}
	else if (number == 95)
	{
		this->Register("Am");
	}
	else if (number == 96)
	{
		this->Register("Cm");
	}
	else if (number == 97)
	{
		this->Register("Bk");
	}
	else if (number == 98)
	{
		this->Register("Cf");
	}
	else if (number == 99)
	{
		this->Register("Es");
	}
	else if (number == 100)
	{
		this->Register("Fm");
	}
	else if (number == 101)
	{
		this->Register("Md");
	}
	else if (number == 102)
	{
		this->Register("No");
	}
	else if (number == 103)
	{
		this->Register("Lr");
	}
	else if (number == 104)
	{
		this->Register("Rf");
	}
	else if (number == 105)
	{
		this->Register("Db");
	}
	else if (number == 106)
	{
		this->Register("Sg");
	}
	else if (number == 107)
	{
		this->Register("Bh");
	}
	else if (number == 108)
	{
		this->Register("Hs");
	}
	else if (number == 109)
	{
		this->Register("Mt");
	}
	else if (number == 110)
	{
		this->Register("Ds");
	}
	else if (number == 111)
	{
		this->Register("Rg");
	}
	else if (number == 112)
	{
		this->Register("Cn");
	}
	else if (number == 113)
	{
		this->Register("Nh");
	}
	else if (number == 114)
	{
		this->Register("Fl");
	}
	else if (number == 115)
	{
		this->Register("Mc");
	}
	else if (number == 116)
	{
		this->Register("Lv");
	}
	else if (number == 117)
	{
		this->Register("Ts");
	}
	else if (number == 118)
	{
		this->Register("Og");
	}
	else
	{
		std::cout << "Given Number: " << number << "\tValid Option: 0-118" << std::endl;
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
		this->atomic_radii = 1.20;
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

//Edit melting point
void Atom::editMeltingPoint(double val)
{
	this->melting_point = val;
}

//Edit boiling point
void Atom::editBoilingPoint(double val)
{
	this->boiling_point = val;
}

//Edit thermal
void Atom::editThermalXSection(double val)
{
	this->thermal_x_sec = val;
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

//Return melting point
double Atom::MeltingPoint()
{
	return this->melting_point;
}

//Return boiling point
double Atom::BoilingPoint()
{
	return this->boiling_point;
}

//Return thermal
double Atom::ThermalXSection()
{
	return this->thermal_x_sec;
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
	//Atom particle1(name1), particle2("Oganesson");
	Atom particle1, particle2;
	particle1.Register("H");
	particle1.Register("Og");
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
	
	//All tests passed
	//-------------- Alex Tests --------------------------------------------------------
	
	time = clock() - time;
	std::cout << "Runtime (s): " << (time/CLOCKS_PER_SEC) << std::endl;
	return success;
}
