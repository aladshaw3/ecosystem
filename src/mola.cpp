//----------------------------------------
//  Created by Austin Ladshaw on 02/24/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#include "mola.h"

//Default constructor
Molecule::Molecule()
:
Atom(0)
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
		std::cout << formula << " Not Found!\n";
		mError(unregistered_name); return;
	}
	else if (first == 'B')
	{
		//List of molecules starting with B
		mError(unregistered_name); return;
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
		mError(unregistered_name); return;
	}
	else if (first == 'L')
	{
		//List of molecules starting with L
		mError(unregistered_name); return;
	}
	else if (first == 'M')
	{
		//List of molecules starting with M
		mError(unregistered_name); return;
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
		else if (formula == "NO3 - (aq)")
		{
			this->Register(-1,-207300.0,146.4,-111300.0,true,true,"Aqueous","Nitrate",formula,"NO3");
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
	
	time = clock() - time;
	std::cout << "\nRuntime (s): " << (time/CLOCKS_PER_SEC) << std::endl;
	return success;
}