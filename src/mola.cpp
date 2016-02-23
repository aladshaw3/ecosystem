/*!
 *  \file mola.cpp mola.h
 *	\brief Molecule Object Library from Atoms
 *  \author Austin Ladshaw
 *	\date 02/24/2014
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "mola.h"

//Default constructor
Molecule::Molecule()
:
atoms(0)
{
	registered = false;
	haveG = false;
	haveHS = false;
	Name = "Unregistered";
	Formula = Name;
	Phase = Name;
	charge = 0;
	molar_weight = 0;
	formation_energy = 0;
	formation_enthalpy = 0;
	formation_entropy = 0;
	PhaseID = OTHER;
}

//Default destructor
Molecule::~Molecule()
{
	
}

//Forming a molecule from its information
Molecule::Molecule(int c,
				   double enthalpy,
				   double entropy,
				   double energy,
				   bool HS,
				   bool G,
				   std::string phase,
				   std::string name,
				   std::string formula,
				   std::string lin_form)
{
	charge = c;
	formation_enthalpy = enthalpy;
	formation_entropy = entropy;
	formation_energy = energy;
	haveHS = HS;
	haveG = G;
	Phase = phase;
	for (int i=0; i<phase.size(); i++)
		phase[i] = tolower(phase[i]);
	if (phase == "aqueous")
		PhaseID = AQUEOUS;
	else if (phase == "solid")
		PhaseID = SOLID;
	else if (phase == "liquid")
		PhaseID = LIQUID;
	else if (phase == "gas")
		PhaseID = GAS;
	else if (phase == "plasma")
		PhaseID = PLASMA;
	else
		PhaseID = OTHER;
	Name = name;
	Formula = formula;
	molar_weight = 0.0;
	registered = true;
	Atom temp;
	std::string c_atom_r, c_atom;
	char c_char;
	int c_num = 1;
	for (long int i=lin_form.length()-1; i>=0; i--)
	{
		c_char = lin_form[i];
		if (isdigit(c_char))
		{
			if (i==0) {mError(string_parse_error); break;}
			c_num = atoi(&c_char);
			int x = 1;
			int a_num;
			char a_char;
			for (long int j=i-1; j>=0; j--)
			{
				a_char = lin_form[j];
				if (!isdigit(a_char)) {break;}
				else
				{
					a_num = atoi(&a_char);
					c_num = c_num + a_num*pow(10,x);
					x++;
					i--;
				}
			}
		}
		else
		{
			c_atom_r+=c_char;
			if (isupper(c_char))
			{
				long int r = c_atom_r.length()-1;
				for (int n=0; n<c_atom_r.length(); n++)
				{
					c_atom.push_back(c_atom_r.at(r));
					r--;
				}
				temp.Register(c_atom);
				for (int j=0; j<c_num; j++)
					atoms.push_back(temp);
				c_num = 1;
				c_atom.clear();
				c_atom_r.clear();
			}
			else{/*No Action*/}
		}
	}
	for (int i=0; i<atoms.size(); i++)
		molar_weight+=atoms[i].AtomicWeight();
}

//Register a molecule given all the necessary information
void Molecule::Register(int c,
						double enthalpy,
						double entropy,
						double energy,
						bool HS,
						bool G,
						std::string phase,
						std::string name,
						std::string formula,
						std::string lin_form)
{
	this->charge = c;
	this->formation_enthalpy = enthalpy;
	this->formation_entropy = entropy;
	this->formation_energy = energy;
	this->haveHS = HS;
	this->haveG = G;
	this->Phase = phase;
	for (int i=0; i<phase.size(); i++)
		phase[i] = tolower(phase[i]);
	if (phase == "aqueous")
		this->PhaseID = AQUEOUS;
	else if (phase == "solid")
		this->PhaseID = SOLID;
	else if (phase == "liquid")
		this->PhaseID = LIQUID;
	else if (phase == "gas")
		this->PhaseID = GAS;
	else if (phase == "plasma")
		this->PhaseID = PLASMA;
	else
		this->PhaseID = OTHER;
	this->Name = name;
	this->Formula = formula;
	this->molar_weight = 0.0;
	this->atoms.clear();
	this->registered = true;
	if (lin_form == "" || lin_form == " ")
		return;
	Atom temp;
	std::string c_atom_r, c_atom;
	char c_char;
	int c_num = 1;
	for (long int i=lin_form.length()-1; i>=0; i--)
	{
		c_char = lin_form[i];
		if (isdigit(c_char))
		{
			if (i==0) {mError(string_parse_error); return;}
			c_num = atoi(&c_char);
			int x = 1;
			int a_num;
			char a_char;
			for (long int j=i-1; j>=0; j--)
			{
				a_char = lin_form[j];
				if (!isdigit(a_char)) {break;}
				else
				{
					a_num = atoi(&a_char);
					c_num = c_num + a_num*pow(10,x);
					x++;
					i--;
				}
			}
		}
		else
		{
			c_atom_r+=c_char;
			if (isupper(c_char))
			{
				long int r = c_atom_r.length()-1;
				for (int n=0; n<c_atom_r.length(); n++)
				{
					c_atom.push_back(c_atom_r.at(r));
					r--;
				}
				temp.Register(c_atom);
				for (int j=0; j<c_num; j++)
					this->atoms.push_back(temp);
				c_num = 1;
				c_atom.clear();
				c_atom_r.clear();
			}
			else{/*No Action*/}
		}
	}
	for (int i=0; i<this->atoms.size(); i++)
		this->molar_weight+=this->atoms[i].AtomicWeight();
}

//Register a molecule based on a hard-coded library of known molecules
void Molecule::Register(std::string formula)
{
	if (formula.length() == 0) {mError(unregistered_name); return;}
	char first = formula[0];
	
	//Check to see if the second character needs to be read instead
	if (first == '(')
		first = formula[1];
	
	//Searching is done in alphabetical order by first atom in formula for greater efficiency
	if (first == 'A')
	{
		//List of molecules starting with A
		if (formula == "Ag (s)")
		{
			this->Register(0, 0.0, 42.6, 0.0, true, true, "Solid", "Silver", formula, "Ag");
		}
		else if (formula == "Ag + (aq)")
		{
			this->Register(1, 105600.0, 73.4, 77120.0, true, true, "Aqueous", "Silver", formula, "Ag");
		}
		else if (formula == "AgBr (s)")
		{
			this->Register(0, -100600.0, 107.0, -96900.0, true, true, "Solid", "Silver-bromide", formula, "AgBr");
		}
		else if (formula == "AgCl (s)")
		{
			this->Register(0, -127100.0, 96.0, -109800.0, true, true, "Solid", "Silver-chloride", formula, "AgCl");
		}
		else if (formula == "AgI (s)")
		{
			this->Register(0, -61840.0, 115.0, -66200.0, true, true, "Solid", "Silver-iodide", formula, "AgI");
		}
		else if (formula == "Ag2S (s)")
		{
			this->Register(0, -29400.0, 14.0, -40700.0, true, true, "Solid", "DiSilver-sulfide", formula, "Ag2S");
		}
		else if (formula == "AgOH (aq)")
		{
			this->Register(0, 0.0, 0.0, -92000.0, false, true, "Aqueous", "Silver-hydroxide", formula, "AgOH");
		}
		else if (formula == "Ag(OH)2 - (aq)")
		{
			this->Register(-1, 0.0, 0.0, -260200.0, false, true, "Aqueous", "Silver-dihydroxide", formula, "AgO2H2");
		}
		else if (formula == "AgCl (aq)")
		{
			this->Register(0, -72800.0, 154.0, -72800.0, true, true, "Aqueous", "Silver-chloride", formula, "AgCl");
		}
		else if (formula == "AgCl2 - (aq)")
		{
			this->Register(-1, -245200.0, 231.0, -215500.0, true, true, "Aqueous", "Silver-dichloride", formula, "AgCl2");
		}
		else if (formula == "Al (s)")
		{
			this->Register(0, 0.0, 28.3, 0.0, true, true, "Solid", "Aluminum", formula, "Al");
		}
		else if (formula == "Al 3+ (aq)")
		{
			this->Register(3, -531000.0, -308.0, -489400.0, true, true, "Aqueous", "Aluminum", formula, "Al");
		}
		else if (formula == "AlOH 2+ (aq)")
		{
			this->Register(2, 0.0, 0.0, -698000.0, false, true, "Aqueous", "Aluminum-hydroxide", formula, "AlOH");
		}
		else if (formula == "Al(OH)2 + (aq)")
		{
			this->Register(1, 0.0, 0.0, -911000.0, false, true, "Aqueous", "Aluminum-dihydroxide", formula, "AlO2H2");
		}
		else if (formula == "Al(OH)3 (aq)")
		{
			this->Register(0, 0.0, 0.0, -1115000.0, false, true, "Aqueous", "Aluminum-trihydroxide", formula, "AlO3H3");
		}
		else if (formula == "Al(OH)4 - (aq)")
		{
			this->Register(-1, 0.0, 0.0, -1325000.0, false, true, "Aqueous", "Aluminum-tetrahydroxide", formula, "AlO4H4");
		}
		else if (formula == "Al2O3 (s)")
		{
			this->Register(0, -1676000.0, 50.9, -1582000.0, true, true, "Solid", "Corundum", formula, "Al2O3");
		}
		else if (formula == "AlOOH (s)")
		{
			this->Register(0, -1000000.0, 17.8, -922000.0, true, true, "Solid", "Boehmite", formula, "AlOOH");
		}
		else if (formula == "Al(OH)3 (s)")
		{
			this->Register(0, -1293000.0, 68.4, -1155000.0, true, true, "Solid", "Gibbsite", formula, "AlO3H3");
		}
		else if (formula == "Al2Si2(OH)4 (s)")
		{
			this->Register(0, -4120000.0, 203.0, -3799000.0, true, true, "Solid", "Kaolinite", formula, "Al2Si2O4H4");
		}
		else if (formula == "As (s)")
		{
			this->Register(0, 0.0, 35.1, 0.0, true, true, "Solid", "Arsenic", formula, "As");
		}
		else if (formula == "AsO4 3- (aq)")
		{
			this->Register(-3, -870300.0, -145.0, -636000.0, true, true, "Aqueous", "Arsenate", formula, "AsO4");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'B')
	{
		//List of molecules starting with B
		if (formula == "Ba 2+ (aq)")
		{
			this->Register(2, -537600.0, 9.6, -560700.0, true, true, "Aqueous", "Barium", formula, "Ba");
		}
		else if (formula == "BaSO4 (s)")
		{
			this->Register(0, -1473000.0, 132.0, -1362000.0, true, true, "Solid", "Barite", formula, "BaSO4");
		}
		else if (formula == "BaCO3 (s)")
		{
			this->Register(0, -1211000.0, 112.0, -1132000.0, true, true, "Solid", "Witherite", formula, "BaCO3");
		}
		else if (formula == "Be 2+ (aq)")
		{
			this->Register(2, -382000.0, -130.0, -380000.0, true, true, "Aqueous", "Beryllium", formula, "Be");
		}
		else if (formula == "Be(OH)2 (s)")
		{
			this->Register(0, -902000.0, 51.9, -815000.0, true, true, "Solid", "Beryllium-dihydroxide", formula, "BeO2H2");
		}
		else if (formula == "Be3(OH)3 3+ (aq)")
		{
			this->Register(3, 0.0, 0.0, -1802000.0, false, true, "Aqueous", "TriBeryllium-trihydroxide", formula, "Be3O3H3");
		}
		else if (formula == "B(OH)4 - (aq)")
		{
			this->Register(-1, -1344000.0, 102.0, -1153300.0, true, true, "Aqueous", "Tetrahydroxy-borate", formula, "BO4H4");
		}
		else if (formula == "Br2 (l)")
		{
			this->Register(0, 0.0, 152.0, 0.0, true, true, "Liquid", "Bromide", formula, "Br2");
		}
		else if (formula == "Br2 (aq)")
		{
			this->Register(0, -2590.0, 130.5, 3930.0, true, true, "Aqueous", "Bromide", formula, "Br2");
		}
		else if (formula == "Br - (aq)")
		{
			this->Register(0, -121500.0, 82.4, -104000.0, true, true, "Aqueous", "Bromine", formula, "Br");
		}
		else if (formula == "BrO - (aq)")
		{
			this->Register(0, -94100.0, 42.0, -33500.0, true, true, "Aqueous", "Hypobromite", formula, "BrO");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'C')
	{
		//List of molecules starting with C
		if (formula == "CO3 2- (aq)")
		{
			this->Register(-2, -677100.0, -56.9, -527900.0, true, true, "Aqueous", "Carbonate", "CO3 2- (aq)", "CO3");
		}
		else if (formula == "Cl - (aq)")
		{
			this->Register(-1, -167200.0, 56.5, -131300.0, true, true, "Aqueous", "Chlorine", formula, "Cl");
		}
		else if (formula == "CaCl2 (aq)")
		{
			this->Register(0, -795420.0, 108.4, -816050.0, true, true, "Aqueous", "Calcium-chloride", formula, "CaCl2");
		}
		else if (formula == "CaAl2Si2O8 (s)")
		{
			this->Register(0, -4243000, 199.0, -4017300.0, true, true, "Solid", "Anorthite", formula, "CaAl2Si2O8");
		}
		else if (formula == "C (s)")
		{
			this->Register(0, 0.0, 152.0, 0.0, true, true, "Solid", "Graphite", formula, "C");
		}
		else if (formula == "CO2 (g)")
		{
			this->Register(0, -393500.0, 213.6, -394370.0, true, true, "Gas", "Carbon-dioxide", formula, "CO2");
		}
		else if (formula == "CH4 (g)")
		{
			this->Register(0, -74800.0, 186.0, -50750.0, true, true, "Gas", "Methane", formula, "CH4");
		}
		else if (formula == "CH4 (aq)")
		{
			this->Register(0, -89040.0, 83.7, -34390.0, true, true, "Aqueous", "Methane", formula, "CH4");
		}
		else if (formula == "CH3OH (aq)")
		{
			this->Register(0, -245900.0, 133.0, -175400.0, true, true, "Aqueous", "Methanol", formula, "CH3OH");
		}
		else if (formula == "CN - (aq)")
		{
			this->Register(-1, 150600.0, 94.1, 172400.0, true, true, "Aqueous", "Cyanide", formula, "CN");
		}
		else if (formula == "CH3COOH (aq)")
		{
			this->Register(0, -485800.0, 179.0, -396600.0, true, true, "Aqueous", "Acetic-Acid", formula, "CH3COOH");
		}
		else if (formula == "CH3COO - (aq)")
		{
			this->Register(-1, -486000.0, 86.6, -369400.0, true, true, "Aqueous", "Acetate", formula, "CH3COO");
		}
		else if (formula == "C2H5OH (aq)")
		{
			this->Register(0, -288300.0, 149.0, -181800.0, true, true, "Aqueous", "Ethanol", formula, "C2H5OH");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'D')
	{
		//List of molecules starting with D
		mError(unregistered_name); return;
	}
	else if (first == 'E')
	{
		//List of molecules starting with E
		mError(unregistered_name); return;
	}
	else if (first == 'F')
	{
		//List of molecules starting with F
		mError(unregistered_name); return;
	}
	else if (first == 'G')
	{
		//List of molecules starting with G
		mError(unregistered_name); return;
	}
	else if (first == 'H')
	{
		//List of molecules starting with H
		if (formula == "H2O (l)")
		{
			this->Register(0,-285830.0,69.95,-23780.0,true,true,"Liquid","Water","H2O (l)","H2O");
		}
		else if (formula == "H + (aq)")
		{
			this->Register(1,0.0,0.0,0.0,true,true,"Aqueous","Proton","H + (aq)","H");
		}
		else if (formula == "H2CO3 (aq)")
		{
			this->Register(0, -699700.0, 187.0, -623200.0, true, true, "Aqueous", "Carbonic-Acid", "H2CO3 (aq)", "H2CO3");
		}
		else if (formula == "HCO3 - (aq)")
		{
			this->Register(-1, -692000.0, 91.2, -586800.0, true, true, "Aqueous", "Bicarbonate", "HCO3 - (aq)", "HCO3");
		}
		else if (formula == "HNO3 (aq)")
		{
			this->Register(0, -207300.0, 146.0, -111300.0, true, true, "Aqueous", "Nitric-Acid", formula, "HNO3");
		}
		else if (formula == "HCl (aq)")
		{
			this->Register(0,-167160.0,56.48,-131300.0,true,true,"Aqueous","Hydrochloric-Acid","HCl (aq)","HCl");
		}
		else if (formula == "H3AsO4 (aq)")
		{
			this->Register(0, -898700.0, 206.0, -766000.0, true, true, "Aqueous", "Arsenic-Acid", formula, "H3AsO4");
		}
		else if (formula == "H2AsO4 - (aq)")
		{
			this->Register(-1, -904500.0, 117.0, -748500.0, true, true, "Aqueous", "Dihydrogen-arsenate", formula, "H2AsO4");
		}
		else if (formula == "HAsO4 2- (aq)")
		{
			this->Register(-2, -898700.0, 3.8, -707100.0, true, true, "Aqueous", "Hydrogen-arsenate", formula, "HAsO4");
		}
		else if (formula == "H2AsO3 - (aq)")
		{
			this->Register(-1, 0.0, 0.0, -587400.0, false, true, "Aqueous", "Dihydrogen-arsenite", formula, "H2AsO3");
		}
		else if (formula == "H3BO3 (aq)")
		{
			this->Register(0, -1072000.0, 162.0, -968700.0, true, true, "Aqueous", "Boric-Acid", formula, "H3BO3");
		}
		else if (formula == "HBrO (aq)")
		{
			this->Register(0, -113000.0, 147.0, -82200.0, true, true, "Aqueous", "Hypobromous-Acid", formula, "HBrO");
		}
		else if (formula == "HCOOH (aq)")
		{
			this->Register(0, -425400.0, 163.0, -372300.0, true, true, "Aqueous", "Formic-Acid", formula, "HCOOH");
		}
		else if (formula == "HCOO - (aq)")
		{
			this->Register(-1, -425600.0, 92.0, -351000.0, true, true, "Aqueous", "Formate", formula, "HCOO");
		}
		else if (formula == "HCN (aq)")
		{
			this->Register(0, 107100.0, 124.6, 119700.0, true, true, "Aqueous", "Hydrogen-cyanide", formula, "HCN");
		}
		else
		{
			mError(unregistered_name); return;
		}
		
	}
	else if (first == 'I')
	{
		//List of molecules starting with I
		mError(unregistered_name); return;
	}
	else if (first == 'J')
	{
		//List of molecules starting with J
		mError(unregistered_name); return;
	}
	else if (first == 'K')
	{
		//List of molecules starting with K
		if (formula == "KAl3Si3O10(OH)2 (s)")
		{
			this->Register(0, 0.0, 0.0, -1341000.0, false, true, "Solid", "Muscovite", formula, "KAl3Si3O12H2");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'L')
	{
		//List of molecules starting with L
		mError(unregistered_name); return;
	}
	else if (first == 'M')
	{
		//List of molecules starting with M
		if (formula == "Mg(OH)2 (aq)")
		{
			this->Register(0, -744700.0, 64.0, -769400.0, true, true, "Aqueous", "Magnesium-hydroxide", formula, "MgO2H2");
		}
		else if (formula == "Mg5Al2Si3O10(OH)8 (s)")
		{
			this->Register(0, 0.0, 0.0, -1962000.0, false, true, "Solid", "Chlorite", formula, "Mg5Al2Si3O18H8");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'N')
	{
		//List of molecules starting with N
		if (formula == "NaHCO3 (aq)")
		{
			this->Register(0,-945530.0,100.0,-847328.0,true,true,"Aqueous","Sodium-Bicarbonate",formula,"NaHCO3");
		}
		else if (formula == "NaCO3 - (aq)")
		{
			this->Register(-1,-937550.0,-41.85,-797049.0,true,true,"Aqueous","Sodium-Carbonate",formula,"NaCO3");
		}
		else if (formula == "Na2CO3 (aq)")
		{
			this->Register(0,0.0,0.0,-1051600.0,false,true,"Aqueous","DiSodium-Carbonate",formula,"Na2CO3");
		}
		else if (formula == "Na + (aq)")
		{
			this->Register(1,-240100.0,59.0,-261900.0,true,true,"Aqueous","Sodium",formula,"Na");
		}
		else if (formula == "NaCl (aq)")
		{
			this->Register(0,-407300.0,115.5,-393170.0,true,true,"Aqueous","Sodium-Chloride",formula,"NaCl");
		}
		else if (formula == "NaOH (aq)")
		{
			this->Register(0,-470110.0,48.1,-419200.0,true,true,"Aqueous","Sodium-Hydroxide",formula,"NaOH");
		}
		else if (formula == "NaNO3 (aq)")
		{
			this->Register(0,0.0,0.0,-373210.0,false,true,"Aqueous","Sodium-Nitrate",formula,"NaNO3");
		}
		else if (formula == "NO3 - (aq)")
		{
			this->Register(-1,-207300.0,146.4,-111300.0,true,true,"Aqueous","Nitrate",formula,"NO3");
		}
		else if (formula == "NH3 (aq)")
		{
			this->Register(0,0.0,193.0,-26570.0,false,true,"Aqueous","Ammonia",formula,"NH3");
		}
		else if (formula == "NaAlSiO3O8 (s)")
		{
			this->Register(0, -3935100, -749.7, -3711700.0, true, true, "Solid", "Albite", formula, "NaAlSiO3O8");
		}
		else if (formula == "NH2CH2COOH (aq)")
		{
			this->Register(0, -514000.0, 158.0, -370800.0, true, true, "Aqueous", "Glycine", formula, "NH2CH2COOH");
		}
		else if (formula == "NH2CH2COO - (aq)")
		{
			this->Register(-1, -469800.0, 119.0, -315000.0, true, true, "Aqueous", "Glycinate", formula, "NH2CH2COO");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'O')
	{
		//List of molecules starting with O
		if (formula == "OH - (aq)")
		{
			this->Register(-1,-230000.0,-10.75,-157300.0,true,true,"Aqueous","Hydroxide","OH - (aq)","OH");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'P')
	{
		//List of molecules starting with P
		mError(unregistered_name); return;
	}
	else if (first == 'Q')
	{
		//List of molecules starting with Q
		mError(unregistered_name); return;
	}
	else if (first == 'R')
	{
		//List of molecules starting with R
		mError(unregistered_name); return;
	}
	else if (first == 'S')
	{
		//List of molecules starting with S
		mError(unregistered_name); return;
	}
	else if (first == 'T')
	{
		//List of molecules starting with T
		mError(unregistered_name); return;
	}
	else if (first == 'U')
	{
		//List of molecules starting with U
		if (formula == "UO2 2+ (aq)")
		{
			//Note: Information is incomplete
			this->Register(2,-1019000.0,-98.2,-952551.0,true,true,"Aqueous","Uranyl","UO2 2+ (aq)","UO2");
		}
		else if (formula == "UO2NO3 + (aq)")
		{
			//Note: Information is incomplete
			this->Register(1,0.0,0.0,-1065557.0,false,true,"Aqueous","Uranyl-nitrate",formula,"UO2NO3");//-----------
		}
		else if (formula == "UO2(NO3)2 (aq)")
		{
			//Note: Information is incomplete
			this->Register(0,0.0,0.0,-1105803.0,false,true,"Aqueous","Uranyl-dinitrate",formula,"UO2N2O6");
		}
		else if (formula == "UO2OH + (aq)")
		{
			//Note: Information is incomplete
			this->Register(1,-1261371.0,17.0,-1159724.0,true,true,"Aqueous","Uranyl-hydroxide",formula,"UO2OH");
		}
		else if (formula == "UO2(OH)2 (aq)")
		{
			//Note: Information is incomplete
			this->Register(0,0,0,-1357479.0,false,true,"Aqueous","Uranyl-dihydroxide",formula,"UO2O2H2");
		}
		else if (formula == "UO2(OH)3 - (aq)")
		{
			//Note: Information is incomplete
			this->Register(-1,0.0,0.0,-908303.0,false,true,"Aqueous","Uranyl-trihydroxide",formula,"UO2O3H3");
		}
		else if (formula == "UO2(OH)4 2- (aq)")
		{
			//Note: Information is incomplete
			this->Register(-2,0.0,0.0,-1716171.0,false,true,"Aqueous","Uranyl-tetrahydroxide",formula,"UO2O4H4");
		}
		else if (formula == "(UO2)2OH 3+ (aq)")
		{
			//Note: Information is incomplete
			this->Register(3,0.0,0.0,-2126830.0,false,true,"Aqueous","Diuranyl-hydroxide",formula,"U2O4OH");
		}
		else if (formula == "(UO2)2(OH)2 2+ (aq)")
		{
			//Note: Information is incomplete
			this->Register(2,-2572065.0,-38.0,-2347303.0,true,true,"Aqueous","Diuranyl-dihydroxide",formula,"U2O4O2H2");
		}
		else if (formula == "(UO2)3(OH)4 2+ (aq)")
		{
			//Note: Information is incomplete
			this->Register(2,0.0,0.0,-3738288.0,false,true,"Aqueous","Triuranyl-tetrahydroxide",formula,"U3O6O4H4");
		}
		else if (formula == "(UO2)3(OH)5 + (aq)")
		{
			//Note: Information is incomplete
			this->Register(1,-4389086.0,83.0,-3954594.0,true,true,"Aqueous","Triuranyl-pentahydroxide",formula,"U3O6O5H5");
		}
		else if (formula == "(UO2)3(OH)7 - (aq)")
		{
			//Note: Information is incomplete
			this->Register(-1,0,0,-4333835.0,false,true,"Aqueous","Triuranyl-heptahydroxide",formula,"U3O6O7H7");
		}
		else if (formula == "(UO2)4(OH)7 + (aq)")
		{
			//Note: Information is incomplete
			this->Register(1,0,0,-5345179.0,false,true,"Aqueous","Tetrauranyl-heptahydroxide",formula,"U4O8O7H7");
		}
		else if (formula == "UO2CO3 (aq)")
		{
			//Note: Information is incomplete
			this->Register(0,-1689230.0,58.870,-1537188.0,true,true,"Aqueous","Uranyl-carbonate",formula,"UO2CO3");
		}
		else if (formula == "UO2(CO3)2 2- (aq)")
		{
			//Note: Information is incomplete
			this->Register(-2,-2350960.0,181.846,-2103161.0,true,true,"Aqueous","Uranyl-dicarbonate",formula,"UO2C2O6");
		}
		else if (formula == "UO2(CO3)3 4- (aq)")
		{
			//Note: Information is incomplete
			this->Register(-4,-3083890.0,38.446,-2660914.0,true,true,"Aqueous","Uranyl-tricarbonate","UO2(CO3)3 4- (aq)","UO2C3O9");
		}
		else
		{
			mError(unregistered_name); return;
		}
	}
	else if (first == 'V')
	{
		//List of molecules starting with V
		mError(unregistered_name); return;
	}
	else if (first == 'W')
	{
		//List of molecules starting with W
		mError(unregistered_name); return;
	}
	else if (first == 'X')
	{
		//List of molecules starting with X
		mError(unregistered_name); return;
	}
	else if (first == 'Y')
	{
		//List of molecules starting with Y
		mError(unregistered_name); return;
	}
	else if (first == 'Z')
	{
		//List of molecules starting with Z
		mError(unregistered_name); return;
	}
	else {mError(unregistered_name); return;}
	this->registered = true;
}

//Set the formula for a given molecule (useful for working with unregistered molecules)
void Molecule::setFormula(std::string form)
{
	this->Formula = form;
}

//Function forces recalculation of molar weight based on current atomic makeup
void Molecule::recalculateMolarWeight()
{
	this->molar_weight = 0.0;
	for (int i=0; i<this->atoms.size(); i++)
	{
		this->molar_weight+=this->atoms[i].AtomicWeight();
	}
}

//Set the molar weight of the molecule to a constant
void Molecule::setMolarWeigth(double mw)
{
	this->molar_weight = mw;
}

//Function to change the ionic charge of a molecule
void Molecule::editCharge(int c)
{
	//Note: this is not particularly useful as we don't know which atom lost or gained electrons
	this->charge = c;
}

//Edit the oxidation state of one atom
void Molecule::editOneOxidationState(int state, std::string Symbol)
{
	bool changed = false;
	for (int i=0; i<this->atoms.size(); i++)
	{
		if (Symbol == this->atoms[i].AtomSymbol())
		{
			this->atoms[i].editOxidationState(state);
			changed = true;
			break;
		}
	}
	if (changed == false) {mError(invalid_atom);}
}

//Edit all oxidation states of a given atom
void Molecule::editAllOxidationStates(int state, std::string Symbol)
{
	bool changed = false;
	for (int i=0; i<this->atoms.size(); i++)
	{
		if (Symbol == this->atoms[i].AtomSymbol())
		{
			this->atoms[i].editOxidationState(state);
			changed = true;
		}
	}
	if (changed == false) {mError(invalid_atom);}
}

//Calculate the average oxidation state of the given atom based on charge and other atoms oxidation states
void Molecule::calculateAvgOxiState(std::string Symbol)
{
	//WARNING! This is an incomplete function! (see below)
	bool changed = false;
	int sum = 0;
	int count = 0;
	std::vector<int> indices;
	for (int i=0; i<this->atoms.size(); i++)
	{
		if (Symbol != this->atoms[i].AtomSymbol())
		{
			sum+=this->atoms[i].OxidationState();
		}
		else
		{
			changed = true;
			count++;
			indices.push_back(i);
		}
	}
	if (changed == false) {mError(invalid_atom); return;}
	int remaining = this->charge - sum;
	remaining = remaining / count;
	for (int i=0; i<indices.size(); i++)
	{
		//Note: Accuracy of this function may be bad if (remaining / count) is not a true int
		this->atoms[indices[i]].editOxidationState(remaining);
	}
}

//Edit the formation enthalpy of the molecule
void Molecule::editEnthalpy(double enthalpy)
{
	this->formation_enthalpy = enthalpy;
}

//Edit the formation enthalpy of the molecule
void Molecule::editEntropy(double entropy)
{
	this->formation_entropy = entropy;
}

//Edit both enthalpy and entropy
void Molecule::editHS(double H, double S)
{
	this->formation_enthalpy = H;
	this->formation_entropy = S;
	this->haveHS = true;
}

//Edit the Gibb's formation energy
void Molecule::editEnergy(double energy)
{
	this->formation_energy = energy;
	this->haveG = true;
}

//Function to remove one atom from the molecule (does not affect formula or name)
void Molecule::removeOneAtom(std::string Symbol)
{
	bool changed = false;
	for (int i=0; i<this->atoms.size(); i++)
	{
		if (Symbol == this->atoms[i].AtomSymbol())
		{
			this->atoms.erase(this->atoms.begin()+i);
			changed = true;
			break;
		}
	}
	if (changed == false) {mError(invalid_atom); return;}
	this->recalculateMolarWeight();
}

//Function to remove all atoms of symbol given (does not affect formula or name)
void Molecule::removeAllAtoms(std::string Symbol)
{
	bool changed = false;
	for (int i=0; i<this->atoms.size(); i++)
	{
		if (Symbol == this->atoms[i].AtomSymbol())
		{
			this->atoms.erase(this->atoms.begin()+i);
			i--;
			changed = true;
		}
	}
	if (changed == false) {mError(invalid_atom); return;}
	this->recalculateMolarWeight();
}

//Return the current ionic charge
int Molecule::Charge()
{
	return this->charge;
}

//Return the current molecular weight
double Molecule::MolarWeight()
{
	this->recalculateMolarWeight();
	return this->molar_weight;
}

//Return whether or not enthalpy and entropy are known
bool Molecule::HaveHS()
{
	return this->haveHS;
}

//Return whether or not Gibb's energy is known
bool Molecule::HaveEnergy()
{
	return this->haveG;
}

//Return whether or not the molecule was registered
bool Molecule::isRegistered()
{
	return this->registered;
}

//Return the enthalpy
double Molecule::Enthalpy()
{
	return this->formation_enthalpy;
}

//Return the entropy
double Molecule::Entropy()
{
	return this->formation_entropy;
}

//Return the Gibb's energy
double Molecule::Energy()
{
	return this->formation_energy;
}

//Return molecule common name
std::string Molecule::MoleculeName()
{
	return this->Name;
}

//Return the molecular formula
std::string Molecule::MolecularFormula()
{
	return this->Formula;
}

//Return the molecule's phase
std::string Molecule::MoleculePhase()
{
	return this->Phase;
}

//Return the enum flag of the phase
int Molecule::MoleculePhaseID()
{
	return this->PhaseID;
}

//Display molecule information
void Molecule::DisplayInfo()
{	
	std::cout << "\nCommon Name: " << this->Name << "\tFormula: " << this->Formula << std::endl;
	std::cout << "-------------------------------------------------\n";
	std::cout << "Molar Weight (g/mol): " << this->molar_weight << "\tIonic Charge: " << this->charge << std::endl;
	std::cout << "Phase: " << this->Phase << "\nSTP Formation Energy (J/mol): ";
	if (this->haveG == true) std::cout << this->formation_energy << "\n";
	else std::cout << "Unknown\n";
	std::cout << "STP Formation Enthalpy (J/mol): ";
	if (this->haveHS) std::cout << this->formation_enthalpy << "\t";
	else std::cout << "Unknown\t";
	std::cout << "\nSTP Formation Entropy (J/K/mol): ";
	if (this->haveHS) std::cout << this->formation_entropy << std::endl;
	else std::cout << "Unknown\n\n";
}

//Testing for MOLA
int MOLA_TESTS()
{
	int success = 0;
	double time = clock();
	
	//------------- Testing of the Molecule Object---------------------------------------
	Molecule H2O(0,-285830.0,69.95,-23780.0,true,true,"Liquid","Water","H2O (l)","H2O");
	Molecule HHe(0,0,0,0,false,false,"N/A","N/A","H10He22","H10He22");
	
	HHe = H2O;
	HHe.Register(0,0,0,0,false,false,"N/A","N/A","H10He22","H10He22");
	
	H2O.DisplayInfo();
	HHe.DisplayInfo();
	
	Molecule NewH2O;
	
	NewH2O.Register("H2O (l)"); //Example of registering molecule information based on just the formula
	NewH2O.DisplayInfo();
	
	Molecule H,OH;
	H.Register("H + (aq)");		//Note: Registration must follow the standard naming convention
	OH.Register("OH - (aq)");	//		formula \space charge \space phase
	H.DisplayInfo();
	OH.DisplayInfo();
	
	Molecule Uranyltricarb;
	Uranyltricarb.Register("UO2(CO3)3 4- (aq)");
	Uranyltricarb.editOneOxidationState(4, "U"); //Set the oxidation state of U to 4+
	Uranyltricarb.editAllOxidationStates(-2, "O");//Set all oxidation states of O to -2
	Uranyltricarb.calculateAvgOxiState("C");		//Set all oxidation states of C based on prior info
	Uranyltricarb.editEnergy(555);					//Edits the formation energy (note: this is not correct info)
	Uranyltricarb.removeOneAtom("U");				//Removes the U atom from the molecule
	std::cout << "\n\nTesting of atom removals and other edits to the molecule\n";
	Uranyltricarb.removeAllAtoms("O");				//Removes all O atoms from the molecule
	Uranyltricarb.DisplayInfo();
	//-------------- END molecule testing ----------------------------------------------
	
	Molecule NaHCO3;
	NaHCO3.Register("NaHCO3 (aq)");
	NaHCO3.DisplayInfo();
	
	Molecule Mus;
	Mus.Register("KAl3Si3O10(OH)2 (s)");
	Mus.DisplayInfo();
	
	time = clock() - time;
	std::cout << "\nRuntime (s): " << (time/CLOCKS_PER_SEC) << std::endl;
	return success;
}