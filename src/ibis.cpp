/*!
 *  \file ibis.h ibis.cpp
 *	\brief Implicit Branching Isotope System
 *  \author Austin Ladshaw
 *	\date 09/04/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay. Copyright (c) 2018, all
 *             rights reserved.
 */

#include "ibis.h"

/// Function to convert from a starting unit and value to and ending unit and value (returns converted value)
double time_conversion(time_units end_unit, double start_value, time_units start_unit)
{
	double end_value = 0.0;
	
	switch (end_unit)
	{
		case seconds:
			switch (start_unit)
			{
				case seconds:
					end_value = start_value;
					break;
					
				case minutes:
					end_value = start_value * 60.0;
					break;
					
				case hours:
					end_value = start_value * 60.0 * 60.0;
					break;
					
				case days:
					end_value = start_value * 60.0 * 60.0 * 24.0;
					break;
					
				case years:
					end_value = start_value * 60.0 * 60.0 * 24.0 * 365.25;
					break;
					
				default:
					end_value = start_value;
					break;
			}
			break;
			
		case minutes:
			switch (start_unit)
			{
			case seconds:
				end_value = start_value / 60.0;
				break;
				
			case minutes:
				end_value = start_value;
				break;
				
			case hours:
				end_value = start_value * 60.0;
				break;
				
			case days:
				end_value = start_value * 60.0 * 24.0;
				break;
				
			case years:
				end_value = start_value * 60.0 * 24.0 * 365.25;
				break;
				
			default:
				end_value = start_value;
				break;
			}
			break;
			
		case hours:
			switch (start_unit)
			{
			case seconds:
				end_value = start_value / 60.0 / 60.0;
				break;
				
			case minutes:
				end_value = start_value / 60.0;
				break;
				
			case hours:
				end_value = start_value;
				break;
				
			case days:
				end_value = start_value * 24.0;
				break;
				
			case years:
				end_value = start_value * 24.0 * 365.25;
				break;
				
			default:
				end_value = start_value;
				break;
			}
			break;
			
		case days:
			switch (start_unit)
			{
			case seconds:
				end_value = start_value / 60.0 / 60.0 / 24.0;
				break;
				
			case minutes:
				end_value = start_value / 60.0 / 24.0;
				break;
				
			case hours:
				end_value = start_value / 24.0;
				break;
				
			case days:
				end_value = start_value;
				break;
				
			case years:
				end_value = start_value * 365.25;
				break;
				
			default:
				end_value = start_value;
				break;
			}
			break;
			
		case years:
			switch (start_unit)
			{
			case seconds:
				end_value = start_value / 60.0 / 60.0 / 24.0 / 365.25;
				break;
				
			case minutes:
				end_value = start_value / 60.0 / 24.0 / 365.25;
				break;
				
			case hours:
				end_value = start_value / 24.0 / 365.25;
				break;
				
			case days:
				end_value = start_value / 365.25;
				break;
				
			case years:
				end_value = start_value;
				break;
				
			default:
				end_value = start_value;
				break;
			}
			break;
			
		default:
			end_value = start_value;
			break;
	}
	
	return end_value;
}

//Pick out time units
time_units timeunits_choice(std::string &choice)
{
	time_units units = seconds;
	
	std::string copy = choice;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	
	if (copy == "seconds" || copy == "s")
		units = seconds;
	else if (copy == "minutes" || copy == "min")
		units = minutes;
	else if (copy == "hours" || copy == "hr" || copy == "h")
		units = hours;
	else if (copy == "days" || copy == "d")
		units = days;
	else if (copy == "years" || copy == "yr" || copy == "y")
		units = years;
	else
		units = seconds;
	
	return units;
}

//Pick decay mode
decay_mode decaymode_choice(std::string &choice)
{
	decay_mode mode = stable;
	
	std::string copy = choice;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	
	if (copy == "stable")
		mode = stable;
	else if (copy == "alpha")
		mode = alpha;
	else if (copy == "spontaneous-fission")
		mode = spon_fiss;
	else if (copy == "beta+")
		mode = beta_plus;
	else if (copy == "beta-")
		mode = beta_min;
	else if (copy == "isomeric-transition")
		mode = iso_trans;
	else if (copy == "neutron-emission")
		mode = neutron_em;
	else if (copy == "beta-/neutron-emission")
		mode = beta_min_neutron_em;
	else if (copy == "beta+/proton-emission")
		mode = beta_plus_proton_em;
	else if (copy == "proton-emission")
		mode = proton_em;
	else if (copy == "beta+/alpha")
		mode = beta_plus_alpha;
	else if (copy == "beta+/beta+")
		mode = beta_plus_beta_plus;
	else if (copy == "beta-/beta-")
		mode = beta_min_beta_min;
	else if (copy == "beta-/neutron-emission/neutron-emission")
		mode = beta_min_2neutron_em;
	else if (copy == "beta-/alpha")
		mode = beta_min_alpha;
	else if (copy == "proton-emission/proton-emission")
		mode = proton_em_proton_em;
	else if (copy == "neutron-emission/neutron-emission")
		mode = neutron_em_neutron_em;
	else if (copy == "beta-/neutron-emission/neutron-emission/neutron-emission")
		mode = beta_min_3neutron_em;
	else if (copy == "beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission")
		mode = beta_min_4neutron_em;
	else if (copy == "beta+/proton-emission/proton-emission")
		mode = beta_plus_2proton_em;
	else if (copy == "beta+/proton-emission/proton-emission/proton-emission")
		mode = beta_plus_3proton_em;
	else if (copy == "specific-isotope")
		mode = specific_isotope;
	else if (copy == "undefined")
		mode = undefined;
	else
		mode = stable;
	
	return mode;
}

// Return name of mode
std::string decaymode_string(decay_mode mode)
{
	std::string name = "stable";
	
	switch (mode)
	{
		case stable:
			name = "stable";
			break;
			
		case alpha:
			name = "alpha";
			break;
			
		case spon_fiss:
			name = "spontaneous-fission";
			break;
			
		case beta_plus:
			name = "beta+";
			break;
			
		case beta_min:
			name = "beta-";
			break;
			
		case iso_trans:
			name = "isomeric-transition";
			break;
			
		case neutron_em:
			name = "neutron-emission";
			break;
			
		case beta_min_neutron_em:
			name = "beta-/neutron-emission";
			break;
			
		case beta_plus_proton_em:
			name = "beta+/proton-emission";
			break;
			
		case proton_em:
			name = "proton-emission";
			break;
			
		case beta_plus_alpha:
			name = "beta+/alpha";
			break;
			
		case beta_plus_beta_plus:
			name = "beta+/beta+";
			break;
			
		case beta_min_beta_min:
			name = "beta-/beta-";
			break;
			
		case beta_min_2neutron_em:
			name = "beta-/neutron-emission/neutron-emission";
			break;
			
		case beta_min_alpha:
			name = "beta-/alpha";
			break;
			
		case proton_em_proton_em:
			name = "proton-emission/proton-emission";
			break;
			
		case neutron_em_neutron_em:
			name = "neutron-emission/neutron-emission";
			break;
			
		case beta_min_3neutron_em:
			name = "beta-/neutron-emission/neutron-emission/neutron-emission";
			break;
			
		case beta_min_4neutron_em:
			name = "beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission";
			break;
			
		case beta_plus_2proton_em:
			name = "beta+/proton-emission/proton-emission";
			break;
			
		case beta_plus_3proton_em:
			name = "beta+/proton-emission/proton-emission/proton-emission";
			break;
			
		case specific_isotope:
			name = "specific-isotope";
			break;
			
		case undefined:
			name = "undefined";
			break;
			
		default:
			break;
	}
	
	return name;
}

// String of units
std::string timeunits_string(time_units units)
{
	std::string name = "seconds";
	
	switch (units)
	{
		case seconds:
			name = "seconds";
			break;
			
		case minutes:
			name = "minutes";
			break;
			
		case hours:
			name = "hours";
			break;
			
		case days:
			name = "days";
			break;
			
		case years:
			name = "years";
			break;
			
		default:
			break;
	}
	
	return name;
}

/*
 *								Start: Isotope Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
Isotope::Isotope() : Atom(0)
{
	IsoName = "No Name";
	decay_rate = 0.0;
	half_life = INFINITY;
	hl_units = years;
	isotope_number = 0;
	nuclides = nullptr;
	Stable = false;
	IsomericState = false;
	initial_condition = 0.0;
	concentration = 0.0;
}

//Default destructor
Isotope::~Isotope()
{
	decay_modes.clear();
	branch_ratios.clear();
	particle_emitted.clear();
	num_particles.clear();
	daughter.clear();
	chain.clear();
	//nuclides->DeleteContents();  //Must not delete contents on destructor or this invalidates usage of temporary isotopes
}

//Load library
void Isotope::loadNuclides(yaml_cpp_class &data)
{
	this->nuclides = &data;
}

//Unload library
void Isotope::unloadNuclides()
{
	this->nuclides->DeleteContents();
}

//Clear chain
void Isotope::clearChain()
{
	this->chain.clear();
}

//Register via isotope name (e.g., H-2)
void Isotope::registerIsotope(std::string isotope_name)
{
	char *str;
	char iso[256];
	strcpy(iso, isotope_name.c_str());
	str = strtok(iso, "-");
	std::string sym;
	int iso_num;
	
	int i=0;
	while (str != NULL)
	{
		if (i == 0)
			sym = str;
		if (i == 1)
			iso_num = atoi(str);
		str = strtok(NULL, "-");
		i++;
	}
	
	this->registerIsotope(sym, iso_num);
}

//Register via atomic symbol and isotope number
void Isotope::registerIsotope(std::string sym, int iso)
{
	this->Register(sym);
	this->isotope_number = iso;
	this->editNeutrons(iso - this->AtomicNumber());
	this->setConstants();
	this->computeDecayRate();
}

//Register via atomic number and isotope number
void Isotope::registerIsotope(int atom_num, int iso_num)
{
	this->Register(atom_num);
	this->isotope_number = iso_num;
	this->editNeutrons(iso_num - atom_num);
	this->setConstants();
	this->computeDecayRate();
}

//Display isotope information
void Isotope::DisplayInfo()
{
	std::cout << std::endl;
	
	std::cout << "                    Isotope: " << this->IsotopeName() << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Atomic Number: " << this->AtomicNumber() << std::endl;
	std::cout << "Mass Number: " << this->IsotopeNumber() << std::endl;
	std::cout << "Weight (g/mol): " << this->AtomicWeight() << std::endl;
	std::cout << "Isomeric: ";
	if (this->isIsomericState() == true)
		std::cout << "True\n";
	else
		std::cout << "False\n";
	std::cout << "Stable: ";
	if (this->isStable() == true)
		std::cout << "True\n";
	else
		std::cout << "False\n";
	
	if (this->isStable() == false)
	{
		std::cout << "Half-life (" << timeunits_string(this->HalfLifeUnits()) << "): " << this->HalfLife(this->HalfLifeUnits()) << std::endl;
		std::cout << "Decay Rate (s): " << this->DecayRate() << std::endl;
		std::cout << "Decay Modes:\n";
		std::cout << "------------\n";
		std::cout << "\tMode : Fraction : Emission : Number Emitted : Daughter \n";
		std::cout << "\t------------------------------------------------------\n";
		for (int i=0; i<this->DecayModes(); i++)
		{
			std::cout << "\t" << decaymode_string(this->DecayMode(i)) << " : " << this->BranchFraction(i);
			std::cout << " : " << this->ParticleEmitted(i) << " : " << this->NumberParticlesEmitted(i) << " : " << this->Daughter(i);
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
	std::cout << std::endl;
}

//Display chain
void Isotope::DisplayChain()
{
	//Level loop
	for (int i=0; i<this->chain.size(); i++)
	{
		std::cout << "Level " << i << ":\n";
		std::cout << "------------\n";
		//Daughters loop
		for (int j=0; j<this->chain[i].size(); j++)
		{
			std::cout << this->chain[i][j].first << " ---> " << this->chain[i][j].second << std::endl;
		}
		std::cout << std::endl;
	}
}

//Create branching decay chain
void Isotope::createChain()
{
	bool all_stable = true;
	int success = 0;
	std::vector< std::pair<std::string,std::string> > temp;
	int i = 0;
	do
	{
		all_stable = true;
		try
		{
			if (i == 0)
			{
				all_stable = this->isStable();
				
				if (all_stable == true) break;
				else
				{
					this->chain.push_back(temp);
					success = this->addPairs(i, this->IsotopeName());
					if (success != 0) {mError(read_error); all_stable = true; break;}
				}
			}
			else
			{
				for (int j=0; j<this->chain[i-1].size(); j++)
				{
					if (this->chain[i-1][j].second == "stable")
					{
						this->chain[i-1].erase(this->chain[i-1].begin()+j);
						all_stable = true;
					}
					else
						all_stable = this->getNuclideLibrary()(this->chain[i-1][j].second)["stable"].getBool();
					
					if (all_stable == false) break;
				}
				if (all_stable == true)
				{
					//No Action
				}
				else
				{
					this->chain.push_back(temp);
					
					for (int j=0; j<this->chain[i-1].size(); j++)
					{
						if (this->chain[i-1][j].second == "stable")
						{
							this->chain[i-1].erase(this->chain[i-1].begin()+j);
						}
						else
						{
							if (this->getNuclideLibrary()(this->chain[i-1][j].second)["stable"].getBool() == true)
							{
								//Only push_back is daughter for next level in chain is not stable
							}
							else
							{
								success = this->addPairs(i, this->chain[i-1][j].second);
							}
						}
						if (success != 0) {mError(read_error); all_stable = true; break;}
					}
				}
			}
		}
		catch (std::out_of_range)
		{
			mError(invalid_isotope);
			all_stable = true;
		}
		
		i++;
	} while (all_stable == false);
}

//set initial cond
void Isotope::setInitialCondition(double ic)
{
	if (ic < 0.0)
		ic = 0.0;
	this->initial_condition = ic;
}

//set concentration
void Isotope::setConcentration(double c)
{
	this->concentration = c;
}

//Return isotope number
int Isotope::IsotopeNumber()
{
	return this->isotope_number;
}

//Return decay rate
double Isotope::DecayRate()
{
	return this->decay_rate;
}

//Return half-life
double Isotope::HalfLife(time_units units)
{
	return time_conversion(units, this->half_life, this->hl_units);
}

//Return units
time_units Isotope::HalfLifeUnits()
{
	return this->hl_units;
}

//Return name
std::string Isotope::IsotopeName()
{
	return this->IsoName;
}

//Return stability
bool Isotope::isStable()
{
	return this->Stable;
}

//Return isomeric state
bool Isotope::isIsomericState()
{
	return this->IsomericState;
}

//Return number of decay modes
int Isotope::DecayModes()
{
	return (int)this->decay_modes.size();
}

//Return initial cond
double Isotope::getInitialCondition()
{
	return this->initial_condition;
}

//Return concentration
double Isotope::getConcentration()
{
	return this->concentration;
}

//return decay mode
decay_mode Isotope::DecayMode(int i)
{
	if (i < 0 || i >= this->DecayModes())
		return stable;
	else
		return this->decay_modes[i];
		
}

//return branch fraction
double Isotope::BranchFraction(int i)
{
	if (i < 0 || i >= this->DecayModes())
		return 0.0;
	else
		return this->branch_ratios[i];
	
}

//return isotope emission
std::string Isotope::ParticleEmitted(int i)
{
	if (i < 0 || i >= this->DecayModes())
		return "None";
	else
		return this->particle_emitted[i];
	
}

//return the number of particles emitted
int Isotope::NumberParticlesEmitted(int i)
{
	if (i < 0 || i >= this->DecayModes())
		return 0;
	else
		return this->num_particles[i];
}

//return name of daughter
std::string Isotope::Daughter(int i)
{
	if (i < 0 || i >= this->DecayModes())
		return "None";
	else
		return this->daughter[i];
	
}

//return list of decay mode indices that form this isotope
std::vector<int> Isotope::DaughterIndices(std::string parent)
{
	std::vector<int> indices;
	
	//Read from the library to find the isotope
	try
	{
		//Loop through the decay_modes header
		int i = 0;
		for (auto &x: this->getNuclideLibrary()(parent)("decay_modes").getSubMap())
		{
			std::string dau = this->getNuclideLibrary()(parent)("decay_modes")(x.first)["daughter"].getString();
			
			if (dau == this->IsotopeName())
				indices.push_back(i);
			else
				indices.push_back(-1);
			
			i++;
		}
	}
	catch (std::out_of_range)
	{
		//Parent not found
		indices.clear();
	}
	
	return indices;
}

//return list of decay mode indices that form this isotope
std::vector<int> Isotope::EmissionIndices(std::string parent)
{
	std::vector<int> indices;
	
	//Read from the library to find the isotope
	try
	{
		//Loop through the decay_modes header
		int i = 0;
		for (auto &x: this->getNuclideLibrary()(parent)("decay_modes").getSubMap())
		{
			std::string part = this->getNuclideLibrary()(parent)("decay_modes")(x.first)["part_emitted"].getString();
			
			if (part == this->IsotopeName())
				indices.push_back(i);
			else
				indices.push_back(-1);
			
			i++;
		}
	}
	catch (std::out_of_range)
	{
		//Parent not found
		indices.clear();
	}
	
	return indices;
}

//Set the decay information based on a registered atomic number and isotope number
void Isotope::setConstants()
{
	//Set the name of the isotope
	char buff[10];
	sprintf(buff, "-%i", this->IsotopeNumber());
	std::string name = this->Atom::AtomSymbol();
	name.append(buff);
	this->IsoName = name;
	
	this->decay_modes.clear();
	this->branch_ratios.clear();
	this->particle_emitted.clear();
	this->num_particles.clear();
	this->daughter.clear();
	
	//Read from the library to find the isotope
	try
	{
		this->editAtomicWeight( this->getNuclideLibrary()(this->IsotopeName())["atom_weight"].getDouble() );
		this->IsomericState = this->getNuclideLibrary()(this->IsotopeName())["isomeric"].getBool();
		this->Stable = this->getNuclideLibrary()(this->IsotopeName())["stable"].getBool();
		if (this->Stable == false)
			this->half_life = this->getNuclideLibrary()(this->IsotopeName())["half_life"].getDouble();
		else
			this->half_life = INFINITY;
		std::string read_units = this->getNuclideLibrary()(this->IsotopeName())["hl_units"].getString();
		this->hl_units = timeunits_choice( read_units );
		
		//Loop through the decay_modes header
		int i = 0;
		for (auto &x: this->getNuclideLibrary()(this->IsotopeName())("decay_modes").getSubMap())
		{
			std::string read_decay = this->getNuclideLibrary()(this->IsotopeName())("decay_modes")(x.first)["type"].getString();
			this->decay_modes.push_back( decaymode_choice(read_decay) );
			this->branch_ratios.push_back(this->getNuclideLibrary()(this->IsotopeName())("decay_modes")(x.first)["branch_frac"].getDouble());
			this->particle_emitted.push_back( this->getNuclideLibrary()(this->IsotopeName())("decay_modes")(x.first)["part_emitted"].getString() );
			this->num_particles.push_back(this->getNuclideLibrary()(this->IsotopeName())("decay_modes")(x.first)["num_parts"].getInt());
			this->daughter.push_back( this->getNuclideLibrary()(this->IsotopeName())("decay_modes")(x.first)["daughter"].getString() );
			i++;
		}
	}
	//If isotope is not found, assume it is stable and set default values
	catch (std::out_of_range)
	{
		mError(invalid_isotope);
		std::cout << std::endl << this->IsoName << std::endl;
		std::cout << "Setting some default values...\n\n";
		
		this->editAtomicWeight(this->IsotopeNumber());
		this->IsomericState = false;
		this->Stable = true;
		this->half_life = INFINITY;
		this->hl_units = years;
		
		this->decay_modes.push_back(stable);
		this->branch_ratios.push_back(0.0);
		this->particle_emitted.push_back("None");
		this->num_particles.push_back(0);
		this->daughter.push_back("None");
	}
}

//Compute decay rate
void Isotope::computeDecayRate()
{
	double hl_sec = time_conversion(seconds, this->half_life, this->hl_units);
	this->decay_rate = log(2.0)/hl_sec;
}

//Append pairs to end of vectors
int Isotope::addPairs(int i, std::string parent)
{
	int success = 0;
	std::pair<std::string, std::string> temp = std::make_pair(parent, "stable");
	
	try
	{
		if (this->getNuclideLibrary()(parent)["stable"].getBool() == true)
		{
			this->chain[i].push_back(temp);
		}
		else
		{
			//Loop through the decay_modes header
			int m = 0;
			for (auto &x: this->getNuclideLibrary()(parent)("decay_modes").getSubMap())
			{
				std::string read_decay = this->getNuclideLibrary()(parent)("decay_modes")(x.first)["type"].getString();
				decay_mode read_mode = decaymode_choice(read_decay);
				std::string dau = this->getNuclideLibrary()(parent)("decay_modes")(x.first)["daughter"].getString();
				std::string part = this->getNuclideLibrary()(parent)("decay_modes")(x.first)["part_emitted"].getString();
				if (read_mode == stable || read_mode == undefined || read_mode == iso_trans)
				{
					//If stable, don't push back
				}
				else
				{
					temp = std::make_pair(parent, part);
					if (part != "None")
					{
						this->chain[i].push_back(temp);
					}
					temp = std::make_pair(parent, dau);
					this->chain[i].push_back(temp);
				}
				m++;
			}
		}
	}
	catch (std::out_of_range)
	{
		mError(invalid_isotope);
		return -1;
	}
	
	return success;
}

//Return library
YamlWrapper& Isotope::getNuclideLibrary()
{
	return this->nuclides->getYamlWrapper();
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Isotope Class Definitions
 */

/*
 *								Start: DecayChain Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
DecayChain::DecayChain()
{
	
}

//Default destructor
DecayChain::~DecayChain()
{
	nuc_list.clear();
	nuc_map.clear();
	stable_list.clear();
	stable_map.clear();
	parents.clear();
	branches.clear();
	CoefMap.clear();
	stable_parents.clear();
	stable_branches.clear();
	stable_CoefMap.clear();
}

//Display list
void DecayChain::DisplayList()
{
	std::cout << "List of Nuclides:\n";
	std::cout << "-----------------\n";
	
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		std::cout << this->nuc_list[i].IsotopeName() << std::endl;
	}
	std::cout << std::endl;
}

//Display stable list
void DecayChain::DisplayStableList()
{
	std::cout << "List of Stable Nuclides:\n";
	std::cout << "------------------------\n";
	
	for (int i=0; i<this->stable_list.size(); i++)
	{
		std::cout << this->stable_list[i].IsotopeName() << std::endl;
	}
	std::cout << std::endl;
}

//Display all stable nuclide and decay chain information
void DecayChain::DisplayStableInfo()
{
	std::cout << "List of Stable Nuclide Information:\n";
	std::cout << "-----------------------------------\n";
	std::cout << "-----------------------------------\n";
	
	for (int i=0; i<this->stable_list.size(); i++)
	{
		std::cout << "Stable Nuc Index: " << i << "\tName: " << this->getStableIsotope(i).IsotopeName() << std::endl;
		for (int j=0; j<this->getStableParentList(i).size(); j++)
		{
			if (j == 0) std::cout << "----------------------- List of Parents ------------------------\n";
			std::cout << "\t" << this->getIsotope( this->getStableParentList(i)[j] ).IsotopeName() << "\tDecay Const = " << this->getIsotope( this->getStableParentList(i)[j] ).DecayRate() << "\tFraction(s): ";
			for (int k=0; k<this->getStableBranchList(i, j).size(); k++)
			{
				std::cout << this->getIsotope( this->getStableParentList(i)[j] ).BranchFraction( this->getStableBranchList(i, j)[k] ) << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//Display all nuclide and decay chain information
void DecayChain::DisplayInfo()
{
	std::cout << "List of Nuclide Information:\n";
	std::cout << "----------------------------\n";
	std::cout << "----------------------------\n";
	
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		std::cout << "Nuclide Index: " << i << "\tName: " << this->getIsotope(i).IsotopeName() << "\tDecay Const = " << this->getIsotope(i).DecayRate() << std::endl;
		for (int j=0; j<this->getParentList(i).size(); j++)
		{
			if (j == 0) std::cout << "----------------------- List of Parents ------------------------\n";
			std::cout << "\t" << this->getIsotope( this->getParentList(i)[j] ).IsotopeName() << "\tDecay Const = " << this->getIsotope( this->getParentList(i)[j] ).DecayRate()	<< "\tFraction(s): ";
			for (int k=0; k<this->getBranchList(i, j).size(); k++)
			{
				std::cout << this->getIsotope( this->getParentList(i)[j] ).BranchFraction( this->getBranchList(i, j)[k] ) << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//Display map
void DecayChain::DisplayMap()
{
	std::cout << "Coefficient Map:\n";
	std::cout << "----------------\n";
	
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		for (int j=0; j<this->nuc_list.size(); j++)
		{
			try
			{
				double val = this->CoefMap[i].at(j);
				if (val < 0 && i!=j)
				{
					std::cout << "Error!\n";
					break;
				}
				std::cout << "X ";
			}
			catch (std::out_of_range)
			{
				std::cout << "O ";
			}
		}
		std::cout << std::endl;
	}
}

//Load library
void DecayChain::loadNuclides(yaml_cpp_class &data)
{
	this->nuclides = &data;
}

//Unload library
void DecayChain::unloadNuclides()
{
	this->nuclides->DeleteContents();
}

//Register initial nuclide
void DecayChain::registerInitialNuclide(std::string isotope_name)
{
	Isotope temp;
	temp.loadNuclides(*this->nuclides);
	temp.registerIsotope(isotope_name);
	this->roughInsertSort(temp);
}

//Register initial nuclide
void DecayChain::registerInitialNuclide(std::string symb, int iso)
{
	Isotope temp;
	temp.loadNuclides(*this->nuclides);
	temp.registerIsotope(symb, iso);
	this->roughInsertSort(temp);
}

//Register initial nuclide
void DecayChain::registerInitialNuclide(int atom_num, int iso_num)
{
	Isotope temp;
	temp.loadNuclides(*this->nuclides);
	temp.registerIsotope(atom_num, iso_num);
	this->roughInsertSort(temp);
}

//Create the decay chains and list of final nuclides
void DecayChain::createChains()
{
	std::string name;
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		for (int j=0; j<this->nuc_list[i].DecayModes(); j++)
		{
			name = this->nuc_list[i].Daughter(j);
			if (name != "None" && name != "none")
				this->registerInitialNuclide(name);
			name = this->nuc_list[i].ParticleEmitted(j);
			if (name != "None" && name != "none")
				this->registerInitialNuclide(name);
		}//End decay_mode loop
		
	}//End nuc_list loop
	
	this->finalSort();
	this->parents.resize(this->nuc_list.size());
	this->branches.resize(this->nuc_list.size());
	this->stable_parents.resize(this->stable_list.size());
	this->stable_branches.resize(this->stable_list.size());
	this->fillOutBranchData();
	this->CoefMap.resize(this->nuc_list.size());
	this->stable_CoefMap.resize(this->stable_list.size());
	this->fillOutCoefMap();
	
	//Unload library to save memory
	this->unloadNuclides();
	
}

//Form eigenvectors
void DecayChain::formEigenvectors()
{
	if (this->CoefMap.size() != this->nuc_list.size())
	{
		mError(empty_matrix);
		std::cout << "Call the 'createChains() function first...\n";
		return;
	}
	
	this->Eigs.set_size((int)this->CoefMap.size(), (int)this->CoefMap.size());
	this->invEigs.set_size((int)this->CoefMap.size(), (int)this->CoefMap.size());
	
	//Loop for the kth eigenvector
	for (int k=0; k<this->nuc_list.size(); k++)
	{
		double eigenvalue = -this->nuc_list[k].DecayRate();
		double sum = 0.0;
		this->Eigs.edit(k, k, 1.0);
		
		//Loop over the ith elements of the kth eigenvector
		for (int i=k+1; i<this->nuc_list.size(); i++)
		{
			sum = 0.0;
			//Iterate through the map (loop over j columns)
			for (std::map<int,double>::iterator jt=this->CoefMap[i].begin(); jt!=this->CoefMap[i].end(); jt++)
			{
				sum = sum + this->Eigs(jt->first,k)*jt->second;
			}
			double diff = -this->nuc_list[i].DecayRate()-eigenvalue;
			if (diff == 0.0)
			{
				mError(unstable_matrix);
				std::cout << "Non-unique eigenvalues breakdown formation of eigenvectors!\n";
				std::cout << "Problem Isotopes: " << this->nuc_list[k].IsotopeName() << "\t" << this->nuc_list[i].IsotopeName() << std::endl;
				//Fix with small adjustment to the kth eigenvalue
				diff = -this->nuc_list[i].DecayRate()-(eigenvalue+sqrt(DBL_EPSILON));
			}
			this->Eigs.edit(i, k, -sum/diff);
		}
		
	}
	
	//Loop for the kth eigenvector
	for (int k=0; k<this->nuc_list.size(); k++)
	{
		this->invEigs.edit(k, k, 1.0);
		double sum = 0.0;
		//Loop over the ith elements of the kth eigenvector
		for (int i=k+1; i<this->nuc_list.size(); i++)
		{
			sum = 0.0;
			//Loop over the jth columns
			for (int j=0; j<i; j++)
			{
				sum = sum + this->invEigs(j,k)*this->Eigs(i,j);
			}
			this->invEigs.edit(i, k, -sum);
		}
	}
}

//verfiy eigenvectors and eigenvalues
void DecayChain::verifyEigenSoln()
{
	if (this->Eigs.rows() != this->Eigs.columns() && this->Eigs.rows() != this->nuc_list.size())
	{
		mError(empty_matrix);
		std::cout << "Call the 'formEigenvectors() function first...\n";
		return;
	}
	
	std::cout << "Verifying Eigen Solution...\n";
	std::cout << "---------------------------\n";
	
	Matrix<double> zero((int)this->CoefMap.size(),1);
	
	//Loop over kth eigenvector and eigenvalue pair
	for (int k=0; k<this->nuc_list.size(); k++)
	{
		double eigenvalue = -this->nuc_list[k].DecayRate();
		double sum = 0.0;
		
		//Loop over the ith elements of the kth eigenvector
		for (int i=0; i<this->nuc_list.size(); i++)
		{
			sum = 0.0;
			
			//Iterate through the map (loop over j columns)
			for (std::map<int,double>::iterator jt=this->CoefMap[i].begin(); jt!=this->CoefMap[i].end(); jt++)
			{
				sum = sum + this->Eigs(jt->first,k)*jt->second;
			}
			//Use relative error for comparison because machine precicision is based on relative error
			zero.edit(i, 0, (sum - eigenvalue*this->Eigs(i,k))/sum);
		}
		
		double error = zero.norm();
		if (error > MIN_TOL)
		{
			std::cout << "Eigen Solution NOT within tolerance (" << MIN_TOL << ") at " << k << "-th eigen pair!\n";
			std::cout << "\tNorm = " << error << std::endl;
		}
		else
		{
			std::cout << "Eigen Solution is within tolerance (" << MIN_TOL << ") at " << k << "-th eigen pair!\n";
		}
	}
}

//Estimate fractionation
void DecayChain::calculateFractionation(double t)
{
	if (this->Eigs.rows() != this->Eigs.columns() && this->Eigs.rows() != this->nuc_list.size())
	{
		mError(empty_matrix);
		std::cout << "Call the 'formEigenvectors() function first...\n";
		return;
	}
	
	//Loop over all ith isotopes
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		double sum_outer = 0.0;
		//Loop over all j columns of the eigenvector matrices
		for (int j=0; j<=i; j++)
		{
			double sum_inner = 0.0;
			//Loop over all k eigenvalues
			for (int k=j; k<=i; k++)
			{
				sum_inner = sum_inner + this->Eigs(i,k)*this->invEigs(k,j)*exp(-this->nuc_list[k].DecayRate()*t);
			}
			sum_outer = sum_outer + sum_inner*this->nuc_list[j].getInitialCondition();
		}
		this->nuc_list[i].setConcentration(sum_outer);
	}
	
	//Loop over all ith stable isotopes
	for (int i=0; i<this->stable_list.size(); i++)
	{
		//Iterate through the map (loop over j columns)
		double j_sum = 0.0;
		for (std::map<int,double>::iterator jt=this->stable_CoefMap[i].begin(); jt!=this->stable_CoefMap[i].end(); jt++)
		{
			int j = jt->first;
			
			//Loop over k
			double k_sum = 0.0;
			for (int k=0; k<=j; k++)
			{
				//Loop over l
				double l_sum = 0.0;
				for (int l=k; l<=j; l++)
				{
					l_sum = l_sum + ( this->Eigs(j,l)*this->invEigs(l,k)*(1.0-exp(-this->nuc_list[l].DecayRate()*t))/this->nuc_list[l].DecayRate() );
				}
				k_sum = k_sum + this->nuc_list[k].getInitialCondition()*l_sum;
			}
			j_sum = j_sum + k_sum*jt->second;
		}
		this->stable_list[i].setConcentration( this->stable_list[i].getInitialCondition() + j_sum );
	}
}

//Return unstable frac at t
double DecayChain::returnUnstableFractionation(int i, double t)
{
	if (i >= (int)this->nuc_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	
	double sum_outer = 0.0;
	//Loop over all j columns of the eigenvector matrices
	for (int j=0; j<=i; j++)
	{
		double sum_inner = 0.0;
		//Loop over all k eigenvalues
		for (int k=j; k<=i; k++)
		{
			sum_inner = sum_inner + this->Eigs(i,k)*this->invEigs(k,j)*exp(-this->nuc_list[k].DecayRate()*t);
		}
		sum_outer = sum_outer + sum_inner*this->nuc_list[j].getInitialCondition();
	}
	this->nuc_list[i].setConcentration(sum_outer);
	
	return this->nuc_list[i].getConcentration();
}

//Return stable frac at t
double DecayChain::returnStableFractionation(int i, double t)
{
	if (i >= (int)this->stable_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	
	//Iterate through the map (loop over j columns)
	double j_sum = 0.0;
	for (std::map<int,double>::iterator jt=this->stable_CoefMap[i].begin(); jt!=this->stable_CoefMap[i].end(); jt++)
	{
		int j = jt->first;
		
		//Loop over k
		double k_sum = 0.0;
		for (int k=0; k<=j; k++)
		{
			//Loop over l
			double l_sum = 0.0;
			for (int l=k; l<=j; l++)
			{
				l_sum = l_sum + ( this->Eigs(j,l)*this->invEigs(l,k)*(1.0-exp(-this->nuc_list[l].DecayRate()*t))/this->nuc_list[l].DecayRate() );
			}
			k_sum = k_sum + this->nuc_list[k].getInitialCondition()*l_sum;
		}
		j_sum = j_sum + k_sum*jt->second;
	}
	this->stable_list[i].setConcentration( this->stable_list[i].getInitialCondition() + j_sum );
	
	return this->stable_list[i].getConcentration();
}

//Return fractionation of given isotope
double DecayChain::returnFractionation(std::string iso_name, double t)
{
	try
	{
		return this->returnUnstableFractionation(this->nuc_map.at(iso_name), t);
	}
	catch (std::out_of_range)
	{
		//Test for stable isotope
	}
	try
	{
		return this->returnStableFractionation(this->stable_map.at(iso_name), t);
	}
	catch (std::out_of_range)
	{
		mError(out_of_bounds);
		return 0.0;
	}
}

//Return num nuc
int DecayChain::getNumberNuclides()
{
	return (int)this->nuc_list.size();
}

//Return num stable nuc
int DecayChain::getNumberStableNuclides()
{
	return (int)this->stable_list.size();
}

//Return unstable index
int DecayChain::getIsotopeIndex(std::string iso_name)
{
	try
	{
		return  this->nuc_map.at(iso_name);
	}
	catch (std::out_of_range)
	{
		mError(out_of_bounds);
		return 0;
	}
}

//Return stable index
int DecayChain::getStableIsotopeIndex(std::string iso_name)
{
	try
	{
		return  this->stable_map.at(iso_name);
	}
	catch (std::out_of_range)
	{
		mError(out_of_bounds);
		return 0;
	}
}

//Return parents
std::vector<int>& DecayChain::getParentList(int i)
{
	if (i >= (int)this->nuc_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	return this->parents[i];
}

//Return stable parents
std::vector<int>& DecayChain::getStableParentList(int i)
{
	if (i >= (int)this->stable_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	return this->stable_parents[i];
}

//Return branches
std::vector<int>& DecayChain::getBranchList(int i, int j)
{
	if (i >= (int)this->nuc_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	if (j >= (int)this->parents[i].size() || j < 0)
	{
		j = 0;
		mError(out_of_bounds);
	}
	return this->branches[i][j];
}

//Return stable branches
std::vector<int>& DecayChain::getStableBranchList(int i, int j)
{
	if (i >= (int)this->stable_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	if (j >= (int)this->stable_parents[i].size() || j < 0)
	{
		j = 0;
		mError(out_of_bounds);
	}
	return this->stable_branches[i][j];
}

//Return isotope
Isotope& DecayChain::getIsotope(int i)
{
	if (i >= (int)this->nuc_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	return this->nuc_list[i];
}

//Return stable isotope
Isotope& DecayChain::getStableIsotope(int i)
{
	if (i >= (int)this->stable_list.size() || i < 0)
	{
		i = 0;
		mError(out_of_bounds);
	}
	return this->stable_list[i];
}

//Return stable or unstable isotope
Isotope& DecayChain::getIsotope(std::string iso_name)
{
	try
	{
		return this->nuc_list[ this->nuc_map.at(iso_name) ];
	}
	catch (std::out_of_range)
	{
		//Test for stable isotope
	}
	try
	{
		return this->stable_list[ this->stable_map.at(iso_name) ];
	}
	catch (std::out_of_range)
	{
		mError(out_of_bounds);
		return this->nuc_list[0];
	}
}

//Return eigenvectors
Matrix<double>& DecayChain::getEigenvectors()
{
	return this->Eigs;
}

//Return inverse eigenvectors
Matrix<double>& DecayChain::getInverseEigenvectors()
{
	return this->invEigs;
}

//Insert an isotope to the initial list
void DecayChain::roughInsertSort(Isotope iso)
{
	//If iso is stable, do not insert
	if (iso.isStable() == true)
	{
		Isotope pivot;
		pivot = iso;
		int i = 0;
		for (i=0; i<this->stable_list.size(); i++)
		{
			//Check is temp == ith stable nuclide
			if (iso.IsotopeName() == this->stable_list[i].IsotopeName())
				return;	//Don't add a redundant isotope
			
			//Check temp vs ith stable nuclide
			if (iso.IsotopeNumber() > this->stable_list[i].IsotopeNumber())
			{
				//Replace ith nuclide with temp and push all other nuclides downward
				pivot = this->stable_list[i];
				this->stable_list[i] = iso;
				iso = pivot;
			}
		}
		this->stable_list.push_back(iso);
		return;
	}
	
	Isotope pivot;
	pivot = iso;
	int i = 0;
	for (i=0; i<this->nuc_list.size(); i++)
	{
		//Check is temp == ith nuclide
		if (iso.IsotopeName() == this->nuc_list[i].IsotopeName())
			return;	//Don't add a redundant isotope
		
		//Check temp vs ith nuclide
		if (iso.IsotopeNumber() > this->nuc_list[i].IsotopeNumber())
		{
			//Replace ith nuclide with temp and push all other nuclides downward
			pivot = this->nuc_list[i];
			this->nuc_list[i] = iso;
			iso = pivot;
		}
	}
	this->nuc_list.push_back(iso);
}

//Perform final sort
void DecayChain::finalSort()
{
	std::vector<int> common_iso;
	std::vector<Isotope> iso_list;
	std::vector<Isotope> sorted_list;
	
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		common_iso.push_back(i);
		iso_list.push_back(this->nuc_list[i]);
		//Look ahead
		for (int j=i+1; j<this->nuc_list.size(); j++)
		{
			if (this->nuc_list[i].IsotopeNumber() == this->nuc_list[j].IsotopeNumber())
			{
				common_iso.push_back(j);
				iso_list.push_back(this->nuc_list[j]);
			}
			else
				break;
		}//End look ahead
		
		//Check common_iso size
		if (common_iso.size() < 2)
		{
			common_iso.clear();
			iso_list.clear();
		}
		else
		{
			sorted_list = this->sameIsoNumSort(iso_list);
			for (int j=0; j<common_iso.size(); j++)
				this->nuc_list[common_iso[j]] = sorted_list[j];
			
			i = i + (int)common_iso.size() - 1;
			common_iso.clear();
			iso_list.clear();
		}
		
	}//End nuc_list loop
}

std::vector<Isotope> DecayChain::sameIsoNumSort(std::vector<Isotope> &list)
{
	std::vector<Isotope> sorted;
	int size = (int)list.size();
	bool daughter = false;
	
	//Loop through all isotopes
	for (int i=0; i<size; i++)
	{
		//Loop through the current list to check j
		for (int j=0; j<list.size(); j++)
		{
			daughter = false;
			
			//Loop through list to check k against j
			for (int k=0; k<list.size(); k++)
			{
				//Check only if we are looking at a different isotope
				if (k!=j && daughter == false)
				{
					//check daughters of k against j
					for (int d=0; d<list[k].DecayModes(); d++)
					{
						if (list[k].Daughter(d) == list[j].IsotopeName())
						{
							daughter = true;
							break;
						}
					}
				}
				
				if (daughter == true)
					break;
			}
			
			if (daughter == false)
			{
				sorted.push_back(list[j]);
				list.erase(list.begin()+j);
				break;
			}
			
		}
		
	}
	
	return sorted;
}

// Fill out all data
void DecayChain::fillOutBranchData()
{
	std::vector<int> temp;
	//Loop over all i nuclides (daughters)
	for (int i=(int)this->nuc_list.size()-1; i>=0; i--)
	{
		//Loop over all j nuclides for j<i
		for (int j=i-1; j>=0; j--)
		{
			bool hasPushed = false;
			//Check daughters and particles emitted from j
			for (int d=0; d<this->nuc_list[j].DecayModes(); d++)
			{
				if (this->nuc_list[j].Daughter(d) == this->nuc_list[i].IsotopeName() ||
					this->nuc_list[j].ParticleEmitted(d) == this->nuc_list[i].IsotopeName())
				{
					if (hasPushed == false)
					{
						this->parents[i].push_back(j);
						hasPushed = true;
					}
					temp.push_back(d);
					
				}
			}
			if (hasPushed == true)
				this->branches[i].push_back(temp);
			temp.clear();
		}
	}
	
	temp.clear();
	//Loop over all i stable nuclides (daughters)
	for (int i=0; i<(int)this->stable_list.size(); i++)
	{
		//Loop over all j unstable nuclides (potential parents)
		for (int j=0; j<(int)this->nuc_list.size(); j++)
		{
			bool hasPushed = false;
			//Check daughters and particles emitted from j
			for (int d=0; d<this->nuc_list[j].DecayModes(); d++)
			{
				if (this->nuc_list[j].Daughter(d) == this->stable_list[i].IsotopeName() ||
					this->nuc_list[j].ParticleEmitted(d) == this->stable_list[i].IsotopeName())
				{
					if (hasPushed == false)
					{
						this->stable_parents[i].push_back(j);
						hasPushed = true;
					}
					temp.push_back(d);
					
				}
			}
			if (hasPushed == true)
				this->stable_branches[i].push_back(temp);
			temp.clear();
		}
	}
}

//Fill out coef map
void DecayChain::fillOutCoefMap()
{
	//Loop over all nuclides (daughters)
	for (int i=0; i<this->nuc_list.size(); i++)
	{
		this->nuc_map[this->nuc_list[i].IsotopeName()] = i;
		
		//Element [i][i] == -DecayConst
		this->CoefMap[i][i] = -this->nuc_list[i].DecayRate();
		
		//Loop over all parents
		for (int j=0; j<this->getParentList(i).size(); j++)
		{
			int J = this->getParentList(i)[j]; //Index in 'Matrix'
			double coef = 0.0;
			double decay = this->getIsotope( this->getParentList(i)[j] ).DecayRate();
			
			//Loop over branches of parent that forms daughter
			for (int k=0; k<this->getBranchList(i, j).size(); k++)
			{
				coef = coef + decay*this->getIsotope( this->getParentList(i)[j] ).BranchFraction( this->getBranchList(i, j)[k] );
			}
			
			//Fill out off-diags
			this->CoefMap[i][J] = coef;
		}
	}
	
	//Loop over all stable nuclides (daughters)
	for (int i=0; i<this->stable_list.size(); i++)
	{
		this->stable_map[this->stable_list[i].IsotopeName()] = i;
		
		//Loop over all parents
		for (int j=0; j<this->getStableParentList(i).size(); j++)
		{
			int J = this->getStableParentList(i)[j]; //Index in 'Matrix'
			double coef = 0.0;
			double decay = this->getIsotope( this->getStableParentList(i)[j] ).DecayRate();
			
			//Loop over branches of parent that forms daughter
			for (int k=0; k<this->getStableBranchList(i, j).size(); k++)
			{
				coef = coef + decay*this->getIsotope( this->getStableParentList(i)[j] ).BranchFraction( this->getStableBranchList(i, j)[k] );
			}
			
			//Fill out map
			this->stable_CoefMap[i][J] = coef;
		}
	}
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: DecayChain Class Definitions
 */

//Test function
int IBIS_TESTS()
{
	int success = 0;
	yaml_cpp_class nuc_data;
	
	//Read in the library (uses ~ 7.3 MB to hold)
	nuc_data.executeYamlRead("database/NuclideLibrary.yml");
	
	//Test decay chain object
	DecayChain test;
	test.loadNuclides(nuc_data);
	
	//test.registerInitialNuclide("Ba-114");
	//test.registerInitialNuclide("U-235");
	//test.registerInitialNuclide("U-238");
	//test.registerInitialNuclide("U-235");		//Not added to list because it is redundant
	//test.registerInitialNuclide("H-1");		//Not added to list because it is stable
	//test.registerInitialNuclide("O-20");
	//test.registerInitialNuclide("F-20");
	//test.registerInitialNuclide("Na-20");
	//test.registerInitialNuclide("N-20");
	//test.registerInitialNuclide("O-19");
	//test.registerInitialNuclide("N-19");
	//test.registerInitialNuclide("Th-234");
	//test.registerInitialNuclide("He-5");
	//test.registerInitialNuclide("Be-8");
	//test.registerInitialNuclide("Rn-222");
	//test.registerInitialNuclide("n-1");
	
	//test.registerInitialNuclide("Te-132");
	//test.registerInitialNuclide("Xe-132");
	
	//test.registerInitialNuclide("Cf-252");
	//test.registerInitialNuclide("H-3");
	//test.registerInitialNuclide("Tm-172");
	
	test.registerInitialNuclide("Cd-121");
	test.registerInitialNuclide("O-20");
	
	test.createChains();					//Creates list of nuclides and sorts the list from parent to daughter
	test.DisplayInfo();
	test.DisplayStableInfo();
	test.formEigenvectors();				//Mandatory before solving
	test.verifyEigenSoln();					//Completely optional
	
	//test.getIsotope( "Te-132" ).setInitialCondition(100.0);
	
	//test.getIsotope("Cf-252").setInitialCondition(1.0);
	//test.getIsotope("Tm-172").setInitialCondition(2.0);
	//test.getIsotope("H-3").setInitialCondition(0.5);
	
	test.getIsotope("Cd-121").setInitialCondition(1.0);
	test.getIsotope("O-20").setInitialCondition(2.0);
	
	std::cout << "Time(s)";
	for (int i=0; i<test.getNumberNuclides(); i++)
	{
		std::cout << "\t" << test.getIsotope(i).IsotopeName();
	}
	for (int i=0; i<test.getNumberStableNuclides(); i++)
	{
		std::cout << "\t" << test.getStableIsotope(i).IsotopeName();
	}
	std::cout << std::endl;
	double time=0;
	for (int i=0; i<=40; i++)
	{
		time = (double)(i)*10.0;
		test.calculateFractionation(time);
		std::cout << time;
		for (int j=0; j<test.getNumberNuclides(); j++)
		{
			std::cout << "\t" << test.getIsotope(j).getConcentration();
			//std::cout << "\t" << test.returnUnstableFractionation(j, time);
		}
		for (int j=0; j<test.getNumberStableNuclides(); j++)
		{
			std::cout << "\t" << test.getStableIsotope(j).getConcentration();
			//std::cout << "\t" << test.returnStableFractionation(j, time);
		}
		std::cout << std::endl;
		
	}
	
	//std::cout << test.returnFractionation("Ne-20", time) << std::endl;
	
	//Clear the library when no longer needed (redundant)
	test.unloadNuclides();
	
	return success;
}
