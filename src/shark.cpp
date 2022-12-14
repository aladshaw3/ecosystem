/*!
 *  \file shark.cpp shark.h
 *	\brief Speciation-object Hierarchy for Adsorption Reactions and Kinetics
 *
 *  \author Austin Ladshaw
 *	\date 05/27/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "shark.h"

/*
 *								Start: MasterSpeciesList
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
MasterSpeciesList::MasterSpeciesList()
:
species(1)
{
	size = 1;
	residual_alkalinity = 0.0;
}

//Default destructor
MasterSpeciesList::~MasterSpeciesList()
{
	species.clear();
}

//Copy constructor
MasterSpeciesList::MasterSpeciesList(MasterSpeciesList const &msl)
:
species(msl.size)
{
	species = msl.species;
	size = msl.size;
	residual_alkalinity = msl.residual_alkalinity;
}


//Equals operator overload
MasterSpeciesList& MasterSpeciesList::operator=(const MasterSpeciesList &msl)
{
	if (this == &msl)
	{
		return *this;
	}
	else
	{
		this->size = msl.size;
		this->species = msl.species;
		this->residual_alkalinity = msl.residual_alkalinity;
		return *this;
	}
}

//Initialize the size of the list
void MasterSpeciesList::set_list_size(int i)
{
	this->species.clear();
	this->species.resize(i);
	this->size = i;
}

//Set the ith species in the list
void MasterSpeciesList::set_species(int i, std::string formula)
{
	//Check args
	if (i >= this->size)
	{
		mError(out_of_bounds);
		return;
	}

	this->species[i].Register(formula);
}

//Set the ith species to custom object
void MasterSpeciesList::set_species(int i, int charge, double enthalpy, double entropy, double energy,
									bool HS, bool G, std::string Phase, std::string Name, std::string Formula, std::string lin_formula)
{
	//Check args
	if (i >= this->size)
	{
		mError(out_of_bounds);
		return;
	}
	this->species[i].Register(charge, enthalpy, entropy, energy, HS, G, Phase, Name, Formula, lin_formula);
}

//Function to display object info
void MasterSpeciesList::DisplayInfo(int i)
{
	//Check args
	if (i >= this->size)
	{
		mError(out_of_bounds);
		return;
	}
	this->species[i].DisplayInfo();
}

//Function to display all species info
void MasterSpeciesList::DisplayAll()
{
	for (int i=0; i<this->size; i++)
	{
		this->DisplayInfo(i);
		std::cout << "\n";
	}
}

//Function to display the given concentrations of all species
void MasterSpeciesList::DisplayConcentrations(Matrix<double> &C)
{
	if (this->list_size() != C.rows())
	{
		mError(dim_mis_match);
		return;
	}
	else
	{
		std::cout << "\nConcentrations of Species:\n---------------------------------\n";
		for (int i=0; i<this->list_size(); i++)
		{
			std::cout << "[ " << this->speciesName(i) << " ] =\t" << C(i,0) << std::endl;
		}
		std::cout << "---------------------------------\n";
	}
}

//Function to set the alkalinity of solution to constant
void MasterSpeciesList::set_alkalinity(double alk)
{
	this->residual_alkalinity = alk;
}

//Function to return list size
int MasterSpeciesList::list_size()
{
	return this->size;
}

//Function to return a reference to a listed molecule
Molecule& MasterSpeciesList::get_species(int i)
{
	if (i >= size)
	{
		mError(invalid_species);
		return species[0];
	}
	else
	{
		return species[i];
	}
}

//Function to return the index that a species is located in
int MasterSpeciesList::get_index(std::string name)
{
	//Function returns an invalid index if name is not in the list
	int index = -1;
	for (int i=0; i<this->species.size(); i++)
	{
		if (this->species[i].MolecularFormula() == name)
		{
			index = i;
			break;
		}
	}
	return index;
}

//Fetch and return charge of ith species in list
double MasterSpeciesList::charge(int i)
{
	double charge = 0.0;

	//Check args
	if (i >= this->size)
	{
		mError(out_of_bounds);
		return 0.0;
	}

	charge = (double) this->species[i].Charge();

	return charge;
}

//Function to fetch the alkalinity of the solution
double MasterSpeciesList::alkalinity()
{
	return this->residual_alkalinity;
}

//Function to return the name/formula of the ith species in the list
std::string MasterSpeciesList::speciesName(int i)
{
	//Check args
	if (i >= this->size)
	{
		mError(out_of_bounds);
		return "Error";
	}
	else
	{
		return this->species[i].MolecularFormula();
	}
}

//Function to evaluate charge balance residual from given concentrations and activities
double MasterSpeciesList::Eval_ChargeResidual(const Matrix<double> &x)
{
	double res = 0.0;

	//Loop for all species in list
	for (int i=0; i<this->list_size(); i++)
	{
		if (this->get_species(i).MoleculePhaseID() == AQUEOUS || this->get_species(i).MoleculePhaseID() == LIQUID)
			res = res + (this->charge(i) * pow(10.0, x(i,0)));
	}
	res = res + this->alkalinity();

	if (isnan(res) || isinf(res))
		res = sqrt(DBL_MAX)/this->list_size();

	return res;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: MasterSpeciesList
 */

/*
 *								Start: Reaction
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
Reaction::Reaction()
:
Stoichiometric(1)
{
	Equilibrium = 0;
	enthalpy = 0.0;
	entropy = 0.0;
	energy = 0.0;
	CanCalcHS = false;
	CanCalcG = false;
	HaveHS = false;
	HaveG = false;
	HaveEquil = false;
}

//Default destructor
Reaction::~Reaction()
{
	Stoichiometric.clear();
}

//Function to initialize the object from the master species list
void Reaction::Initialize_Object(MasterSpeciesList &List)
{
	this->List = &List;
	this->Stoichiometric.resize(this->List->list_size());
	for (int i=0; i<this->Stoichiometric.size(); i++)
		this->Stoichiometric[i] = 0.0;
}

//Function to display all objects in the list
void Reaction::Display_Info()
{

	//Loop over all species to display reactants
	bool first = true;
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (this->Get_Stoichiometric(i) < 0.0)
		{
			if (i == 0 || first == true)
			{
				if (fabs(this->Get_Stoichiometric(i)) > 1.0)
					std::cout << fabs(this->Get_Stoichiometric(i)) << " x { " << this->List->get_species(i).MolecularFormula() << " }";
				else
					std::cout << "{ " << this->List->get_species(i).MolecularFormula() << " }";
				first = false;
			}
			else
			{
				if (fabs(this->Get_Stoichiometric(i)) > 1.0)
					std::cout << " + " << fabs(this->Get_Stoichiometric(i)) << " x { " << this->List->get_species(i).MolecularFormula() << " }";
				else
					std::cout << " + { " << this->List->get_species(i).MolecularFormula() << " }";
			}

		}
		if (i == this->List->list_size() - 1)
			std::cout << " = ";
	}
	//Loop over all species to display products
	first = true;
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (this->Get_Stoichiometric(i) > 0.0)
		{
			if (i == 0 || first == true)
			{
				if (fabs(this->Get_Stoichiometric(i)) > 1.0)
					std::cout << fabs(this->Get_Stoichiometric(i)) << " x { " << this->List->get_species(i).MolecularFormula() << " }";
				else
					std::cout << "{ " << this->List->get_species(i).MolecularFormula() << " }";
				first = false;
			}
			else
			{
				if (fabs(this->Get_Stoichiometric(i)) > 1.0)
					std::cout << " + " << fabs(this->Get_Stoichiometric(i)) << " x { " << this->List->get_species(i).MolecularFormula() << " }";
				else
					std::cout << " + { " << this->List->get_species(i).MolecularFormula() << " }";
			}
		}
		if (i == this->List->list_size() - 1)
			std::cout << "\t:\tlogK = " << this->Get_Equilibrium() << "\n\n";
	}
}

//Set the ith value of the stoichiometry list
void Reaction::Set_Stoichiometric (int i, double v)
{
	//Check args
	if (i >= this->Stoichiometric.size())
	{
		mError(out_of_bounds);
		return;
	}

	this->Stoichiometric[i] = v;
}

//Set the equilibrium value for the reaction
void Reaction::Set_Equilibrium (double v)
{
	this->Equilibrium = v;
	this->HaveEquil = true;
}

//Set the Enthalpy of the reaction
void Reaction::Set_Enthalpy(double H)
{
	this->enthalpy = H;
}

//Set the Entropy of the reaction
void Reaction::Set_Entropy(double S)
{
	this->entropy = S;
}

//Set the Enthalpy and Entropy
void Reaction::Set_EnthalpyANDEntropy(double H, double S)
{
	this->enthalpy = H;
	this->entropy = S;
	this->HaveHS = true;
	this->HaveEquil = true;
}

//Set the Gibb's Energy of the reaction
void Reaction::Set_Energy(double G)
{
	this->energy = G;
	this->HaveG = true;
	this->HaveEquil = true;
}

//Function to check MasterList Reference for species energy info
void Reaction::checkSpeciesEnergies()
{
	bool HS = true, G = true;

	for (int i=0; i<this->List->list_size(); i++)
	{
		//Check if species is involved in reaction
		if (this->Stoichiometric[i] != 0.0)
		{
			//Check if species i has HS
			if (HS == true)
			{
				HS = this->List->get_species(i).HaveHS();
			}
			if (HS == false)
				break;
		}
	}
	this->CanCalcHS = HS;

	for (int i=0; i<this->List->list_size(); i++)
	{
		//Check if species is involved in reaction
		if (this->Stoichiometric[i] != 0.0)
		{
			//Check if species i has G
			if (G == true)
			{
				G = this->List->get_species(i).HaveEnergy();
			}
			if (G == false)
				break;
		}
	}
	this->CanCalcG = G;

	if (this->CanCalcHS == true || this->CanCalcG == true)
		this->HaveEquil = true;
}

//Function to calculate the energies of the reaction
void Reaction::calculateEnergies()
{
	if (this->CanCalcHS == true)
	{
		this->enthalpy = 0.0;
		this->entropy = 0.0;
	}
	if (this->CanCalcG == true)
	{
		this->energy = 0.0;
	}

	for (int i=0; i<this->Stoichiometric.size(); i++)
	{
		if (this->CanCalcHS == true)
		{
			this->enthalpy = this->enthalpy + (this->Stoichiometric[i] * this->List->get_species(i).Enthalpy());
			this->entropy = this->entropy + (this->Stoichiometric[i] * this->List->get_species(i).Entropy());
		}
		if (this->CanCalcG == true)
		{
			this->energy = this->energy + (this->Stoichiometric[i] * this->List->get_species(i).Energy());
		}
	}
}

//Function to calculate equilibrium constant
void Reaction::calculateEquilibrium(double T)
{
	if (this->CanCalcHS == true || this->HaveHS == true)
	{
		this->energy = this->enthalpy - (T * this->entropy);
		this->Equilibrium = log10(exp(-this->energy/(Rstd*T)));
		return;
	}
	else if (this->CanCalcG == true || this->HaveG == true)
	{
		this->Equilibrium = log10(exp(-this->energy/(Rstd*T)));
		return;
	}
	else if (this->HaveEquil == true)
	{
		this->energy = -Rstd*T*log(pow(10.0,this->Equilibrium));
		return;
	}
	else
	{
		mError(missing_information);
		return;
	}
}

//Function to return whether or not equilibrium is known
bool Reaction::haveEquilibrium()
{
	return this->HaveEquil;
}

//Get the stoichiometric constant
double Reaction::Get_Stoichiometric (int i)
{
	//Check args
	if (i >= this->Stoichiometric.size())
	{
		mError(out_of_bounds);
		return 0.0;
	}
	return this->Stoichiometric[i];
}

//Get the equilibrium constant
double Reaction::Get_Equilibrium ()
{
	return this->Equilibrium;
}

//Get the Enthalpy
double Reaction::Get_Enthalpy()
{
	return this->enthalpy;
}

//Get the Entropy
double Reaction::Get_Entropy()
{
	return this->entropy;
}

//Get the Energy
double Reaction::Get_Energy()
{
	return this->energy;
}

//Evaluate the reaction object residual from variable x and activity gama
double Reaction::Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama)
{
	double res = 0.0;

	//Loop for all species in master list
	for (int i=0; i<this->List->list_size(); i++)
	{
		res = res + ( this->Stoichiometric[i] * (log10(gama(i,0)) + x(i,0)) );
	}
	res = res - this->Get_Equilibrium();

	return res;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Reaction
 */

/*
 *								Start: MassBalance
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
MassBalance::MassBalance()
:
Delta(1)
{
	TotalConcentration = 0.0;
	Type = BATCH;
	volume = 1.0;
	flow_rate = 1.0;
	xsec_area = 1.0;
	dt = 0.1;
	InitialConcentration = 0.0;
	InletConcentration = 0.0;
	SteadyState = true;
	ZeroInitialSolids = false;
}

//Default destructor
MassBalance::~MassBalance()
{
	Delta.clear();
}

//Function to initialize the object from the master species list
void MassBalance::Initialize_Object(MasterSpeciesList &List)
{
	this->List = &List;
	this->Delta.resize(this->List->list_size());
	for (int i=0; i<this->Delta.size(); i++)
		this->Delta[i] = 0.0;
}

//Function to display all objects in the list
void MassBalance::Display_Info()
{
	//Loop through all species to display mass balance
	std::cout << this->Get_Name() << "\n";
	bool first = true;
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (i == 0 || first == true)
		{
			if (this->Delta[i] != 0.0)
			{
				if (this->Delta[i] > 1.0)
					std::cout << this->Delta[i] << " x [ " << this->List->get_species(i).MolecularFormula() << " ]";
				else
					std::cout << "[ " << this->List->get_species(i).MolecularFormula() << " ]";
				first = false;
			}
		}
		else
		{
			if (this->Delta[i] != 0.0)
			{
				if (this->Delta[i] > 1.0)
					std::cout << " + " << this->Delta[i] << " x [ " << this->List->get_species(i).MolecularFormula() << " ]";
				else
					std::cout << " + [ " << this->List->get_species(i).MolecularFormula() << " ]";
			}
		}
	}
	std::cout << " = " << this->Get_TotalConcentration() << "\n\n";
}

//Set the ith value of the delta list
void MassBalance::Set_Delta (int i, double v)
{
	//Check args
	if (i >= this->Delta.size())
	{
		mError(out_of_bounds);
		return;
	}

	this->Delta[i] = v;
}

//Set the total concentration
void MassBalance::Set_TotalConcentration (double v)
{
	if (v <= 0.0)
		v = DBL_MIN;
	this->TotalConcentration = v;
}

//Set the reactor type
void MassBalance::Set_Type(int type)
{
	switch (type)
	{
		case BATCH:
			this->Type = BATCH;
			break;

		case CSTR:
			this->Type = CSTR;
			break;

		case PFR:
			this->Type = PFR;
			break;

		default:
			this->Type = BATCH;
			break;
	}
}

//Set the volume of the reactor
void MassBalance::Set_Volume(double v)
{
	if (v <= 0.0)
		v = DBL_MIN;
	this->volume = v;
}

//Set the flow rate of the reactor
void MassBalance::Set_FlowRate(double v)
{
	if (v <= 0.0)
		v = DBL_MIN;
	this->flow_rate = v;
}

//Set the cross sectional area of the reactor
void MassBalance::Set_Area(double v)
{
	if (v <= 0.0)
		v = DBL_MIN;
	this->xsec_area	= v;
}

//Set the time step of the reactor
void MassBalance::Set_TimeStep(double v)
{
	if (v <= 0.0)
		v = 0.1;
	this->dt = v;
}

//Set the initial concentration
void MassBalance::Set_InitialConcentration (double v)
{
	if (v <= 0.0)
		v = DBL_MIN;
	this->InitialConcentration = v;
}

//Set the inlet concentration
void MassBalance::Set_InletConcentration (double v)
{
	if (v <= 0.0)
		v = DBL_MIN;
	this->InletConcentration = v;
}

//Set the steady-state option
void MassBalance::Set_SteadyState(bool ss)
{
	this->SteadyState = ss;
}

//Set the initial solids option
void MassBalance::Set_ZeroInitialSolids(bool ss)
{
	this->ZeroInitialSolids = ss;
}

//Set the name of the mass balance
void MassBalance::Set_Name(std::string name)
{
	this->Name = name;
}

//Get the delta constant
double MassBalance::Get_Delta (int i)
{
	double d = 0.0;

	//Check args
	if (i >= this->Delta.size())
	{
		mError(out_of_bounds);
		return 0.0;
	}

	d = this->Delta[i];
	return d;
}

//Sum up the Deltas and return the total
double MassBalance::Sum_Delta()
{
	double sum = 0.0;
	for (int i=0; i<this->Delta.size(); i++)
		sum = sum + this->Delta[i];
	return sum;
}

//Get the total concentration
double MassBalance::Get_TotalConcentration ()
{
	if (this->TotalConcentration <= 0.0)
		return DBL_MIN;
	else
		return this->TotalConcentration;
}

//Get the reactor type
int MassBalance::Get_Type()
{
	return this->Type;
}

//Get the reactor volume
double MassBalance::Get_Volume()
{
	return this->volume;
}

//Get the reactor flow rate
double MassBalance::Get_FlowRate()
{
	return this->flow_rate;
}

//Get the reactor cross sectional area
double MassBalance::Get_Area()
{
	return this->xsec_area;
}

//Get the time step
double MassBalance::Get_TimeStep()
{
	return this->dt;
}

//Get the initial concentration
double MassBalance::Get_InitialConcentration()
{
	return this->InitialConcentration;
}

//Get the inlet concentration
double MassBalance::Get_InletConcentration()
{
	return this->InletConcentration;
}

//Get the steady-state value
bool MassBalance::isSteadyState()
{
	return this->SteadyState;
}

//Get the initial solids value
bool MassBalance::isZeroInitialSolids()
{
	return this->ZeroInitialSolids;
}

//Get the name of the Mass Balance
std::string MassBalance::Get_Name()
{
	return this->Name;
}

//Function to evaluate the mass balance residual given the concentration matrix x
double MassBalance::Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old)
{
	double res = 0.0;

	if (this->Type == BATCH)
	{
		//Loop for all species in list
		for (int i=0; i<this->List->list_size(); i++)
		{
			res = res + ( this->Delta[i] * pow(10.0, x_new(i,0)) );
		}

		if (this->Get_TotalConcentration() <= DBL_MIN)
			res = (res - this->Get_TotalConcentration())/pow(DBL_EPSILON, 2.0);
		else
			res = (res / this->Get_TotalConcentration()) - 1.0;

		if (isnan(res) || isinf(res))
			res = sqrt(DBL_MAX)/this->List->list_size();

		return res;
	}
	else if (this->Type == CSTR)
	{
		double CT = 0.0;
		double CT_old = 0.0;
		double ST = 0.0;
		double ST_old = 0.0;
		double theta = this->Get_FlowRate() / this->Get_Volume();
		if (this->isSteadyState() == true)
		{
			//Loop for all species in list
			for (int i=0; i<this->List->list_size(); i++)
			{
				if (this->List->get_species(i).MoleculePhaseID() == AQUEOUS || this->List->get_species(i).MoleculePhaseID() == LIQUID)
					CT = CT + ( this->Delta[i] * pow(10.0, x_new(i,0)) );

			}

			if (this->Get_InletConcentration() <= DBL_MIN)
				res = (CT - this->Get_InletConcentration())/pow(DBL_EPSILON, 2.0);
			else
				res = (CT / this->Get_InletConcentration()) - 1.0;

			if (isnan(res) || isinf(res))
				res = sqrt(DBL_MAX)/this->List->list_size();

			return res;
		}
		else
		{
			//Loop for all species in list
			for (int i=0; i<this->List->list_size(); i++)
			{
				if (this->List->get_species(i).MoleculePhaseID() == AQUEOUS || this->List->get_species(i).MoleculePhaseID() == LIQUID)
				{
					CT = CT + ( this->Delta[i] * pow(10.0, x_new(i,0)) );
					CT_old = CT_old + ( this->Delta[i] * pow(10.0, x_old(i,0)) );
				}
				else if (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED)
				{
					ST = ST + ( this->Delta[i] * pow(10.0, x_new(i,0)) );
					ST_old = ST_old + ( this->Delta[i] * pow(10.0, x_old(i,0)) );
				}
				else
				{
					//Nothing
				}
			}

			if (this->Get_InletConcentration() <= DBL_MIN)
				res = ((CT_old + (theta*this->Get_TimeStep()*this->Get_InletConcentration()) - ST + ST_old) - (CT + (theta*this->Get_TimeStep()*CT))) /(this->Get_InitialConcentration()*10.0) ;
			else
				res = ((CT_old + (theta*this->Get_TimeStep()*this->Get_InletConcentration()) - ST + ST_old) - (CT + (theta*this->Get_TimeStep()*CT))) /(this->Get_InletConcentration()*10.0) ;

			if (isnan(res) || isinf(res))
				res = sqrt(DBL_MAX)/this->List->list_size();

			return res;
		}
	}
	else if (this->Type == PFR)
	{
		//----------Default to Batch--------------

		//Loop for all species in list
		for (int i=0; i<this->List->list_size(); i++)
		{
			res = res + ( this->Delta[i] * pow(10.0, x_new(i,0)) );
		}

		if (this->Get_TotalConcentration() <= DBL_MIN)
			res = (res - this->Get_TotalConcentration())/pow(DBL_EPSILON, 2.0);
		else
			res = (res / this->Get_TotalConcentration()) - 1.0;

		if (isnan(res) || isinf(res))
			res = sqrt(DBL_MAX)/this->List->list_size();

		return res;
	}
	else
	{
		//--------------Default to Batch---------------

		//Loop for all species in list
		for (int i=0; i<this->List->list_size(); i++)
		{
			res = res + ( this->Delta[i] * pow(10.0, x_new(i,0)) );
		}

		if (this->Get_TotalConcentration() <= DBL_MIN)
			res = (res - this->Get_TotalConcentration())/pow(DBL_EPSILON, 2.0);
		else
			res = (res / this->Get_TotalConcentration()) - 1.0;

		if (isnan(res) || isinf(res))
			res = sqrt(DBL_MAX)/this->List->list_size();

		return res;
	}

}

//Function to evaluate the mass balance residual given the concentration matrix x
double MassBalance::Eval_IC_Residual(const Matrix<double> &x)
{
	double res = 0.0;

	if (this->Type == BATCH)
	{
		//Loop for all species in list
		for (int i=0; i<this->List->list_size(); i++)
		{
			res = res + ( this->Delta[i] * pow(10.0, x(i,0)) );
		}

		if (this->Get_TotalConcentration() <= DBL_MIN)
			res = (res - this->Get_TotalConcentration())/pow(DBL_EPSILON, 2.0);
		else
			res = (res / this->Get_TotalConcentration()) - 1.0;

		if (isnan(res) || isinf(res))
			res = sqrt(DBL_MAX)/this->List->list_size();
	}
	else
	{

		//Loop for all species in list
		for (int i=0; i<this->List->list_size(); i++)
		{
			if (this->isZeroInitialSolids() == false)
				res = res + ( this->Delta[i] * pow(10.0, x(i,0)) );
			else
			{
				if (this->List->get_species(i).MoleculePhaseID() == AQUEOUS || this->List->get_species(i).MoleculePhaseID() == LIQUID)
					res = res + ( this->Delta[i] * pow(10.0, x(i,0)) );
			}
		}

		if (this->Get_InitialConcentration() <= DBL_MIN)
			res = (res - this->Get_InitialConcentration())/pow(DBL_EPSILON, 2.0);
		else
			res = (res / this->Get_InitialConcentration()) - 1.0;

		if (isnan(res) || isinf(res))
			res = sqrt(DBL_MAX)/this->List->list_size();
	}

	return res;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: MassBalance
 */

/*
 *								Start: UnsteadyReaction
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
UnsteadyReaction::UnsteadyReaction()
{
	Stoichiometric.resize(1);
	Equilibrium = 0;
	initial_value = 0;
	max_value = 0;
	forward_rate = 0.0;
	forward_ref_rate = 0.0;
	reverse_ref_rate = 0.0;
	reverse_rate = 0.0;
	time_step = 0.0;
	enthalpy = 0.0;
	entropy = 0.0;
	energy = 0.0;
	activation_energy = 0.0;
	temperature_affinity = 0.0;
	CanCalcHS = false;
	CanCalcG = false;
	HaveHS = false;
	HaveG = false;
	HaveEquil = false;
	HaveForward = false;
	HaveReverse = false;
	HaveForRef = false;
	HaveRevRef = false;
}

//Default destructor
UnsteadyReaction::~UnsteadyReaction()
{
	Stoichiometric.clear();
}

//Function to initialize the object from the master species list
void UnsteadyReaction::Initialize_Object(MasterSpeciesList &List)
{
	this->List = &List;
	this->Stoichiometric.resize(this->List->list_size());
	for (int i=0; i<this->Stoichiometric.size(); i++)
		this->Stoichiometric[i] = 0.0;
}

//Function to display all objects in the list
void UnsteadyReaction::Display_Info()
{
	std::cout << "d { " << this->List->get_species(this->species_index).MolecularFormula() << " } / dt\t:\tk_forward = " << this->forward_rate << "\t:\tk_reverse = " << this->reverse_rate << "\n";
	this->Reaction::Display_Info();
}

//Set the species index by an integer value
void UnsteadyReaction::Set_Species_Index(int i)
{
	this->species_index = i;
}

//Set the species index by the name of the species in master list
void UnsteadyReaction::Set_Species_Index(std::string formula)
{
	bool foundit = false;
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (formula == this->List->speciesName(i))
		{
			this->species_index = i;
			foundit = true;
			break;
		}
	}

	//Check for error
	if (foundit == false)
	{
		mError(invalid_species);
		this->species_index = 0;
	}
}

//Set the ith value of the stoichiometry list
void UnsteadyReaction::Set_Stoichiometric (int i, double v)
{
	//Check args
	if (i >= this->Stoichiometric.size())
	{
		mError(out_of_bounds);
		return;
	}

	this->Stoichiometric[i] = v;
}

//Set the equilibrium value for the reaction
void UnsteadyReaction::Set_Equilibrium (double v)
{
	this->Equilibrium = v;
	this->HaveEquil = true;
}

//Set the Enthalpy of the reaction
void UnsteadyReaction::Set_Enthalpy(double H)
{
	this->enthalpy = H;
}

//Set the Entropy of the reaction
void UnsteadyReaction::Set_Entropy(double S)
{
	this->entropy = S;
}

//Set both the Enthalpy and Entropy
void UnsteadyReaction::Set_EnthalpyANDEntropy(double H, double S)
{
	this->enthalpy = H;
	this->entropy = S;
	this->HaveHS = true;
	this->HaveEquil = true;
}

//Set both the Gibb's Energy for the reaction
void UnsteadyReaction::Set_Energy(double G)
{
	this->energy = G;
	this->HaveG = true;
	this->HaveEquil = true;
}

//Set the initial value of the unsteady variable
void UnsteadyReaction::Set_InitialValue(double ic)
{
	if (ic <= 0.0)
		ic = DBL_MIN;
	this->initial_value = ic;
}

//Set the maximum value of the variable in mol/L
void UnsteadyReaction::Set_MaximumValue(double max)
{
	if (max >= DBL_MAX)
		max = DBL_MAX;
	this->max_value = max;
}

//Set the forward rate
void UnsteadyReaction::Set_Forward(double forward)
{
	this->forward_rate = forward;
	if (this->HaveEquil == false)
	{
		mError(rxn_rate_error);
		return;
	}
	this->reverse_rate = this->forward_rate / pow(10.0,this->Equilibrium);
	this->HaveForward = true;
	this->HaveReverse = false;
}

//Set the reverse rate
void UnsteadyReaction::Set_Reverse(double reverse)
{
	this->reverse_rate = reverse;
	if (this->HaveEquil == false)
	{
		mError(rxn_rate_error);
		return;
	}
	this->forward_rate = this->reverse_rate * pow(10.0,this->Equilibrium);
	this->HaveForward = false;
	this->HaveReverse = true;
}

//Set the forward reference rate
void UnsteadyReaction::Set_ForwardRef(double Fref)
{
	this->forward_ref_rate = Fref;
	this->HaveForRef = true;
	this->HaveRevRef = false;
}

//Set the reverse reference rate
void UnsteadyReaction::Set_ReverseRef(double Rref)
{
	this->reverse_ref_rate = Rref;
	this->HaveForRef = false;
	this->HaveRevRef = true;
}

//Set the activation energy for the reaction
void UnsteadyReaction::Set_ActivationEnergy(double E)
{
	if (E < 0.0)
		E = 0.0;
	this->activation_energy = E;
}

//Set the temperature affinity parameter
void UnsteadyReaction::Set_Affinity(double b)
{
	this->temperature_affinity = b;
}

//Set time step
void UnsteadyReaction::Set_TimeStep(double dt)
{
	this->time_step = dt;
}

//Check the energies of species in master list
void UnsteadyReaction::checkSpeciesEnergies()
{
	this->Reaction::checkSpeciesEnergies();
}

//Calculate the energies for the reaction
void UnsteadyReaction::calculateEnergies()
{
	this->Reaction::calculateEnergies();
}

//Calculate the equilibrium constant
void UnsteadyReaction::calculateEquilibrium(double T)
{
	this->Reaction::calculateEquilibrium(T);
}

//Calculate the rate constant for the reaction
void UnsteadyReaction::calculateRate(double T)
{
	this->calculateEquilibrium(T);
	if (this->HaveForRef == true)
	{
		this->Set_Forward(this->forward_ref_rate * pow(T, this->temperature_affinity) * exp(-this->activation_energy/(Rstd*T)));
		return;
	}
	else if (this->HaveRevRef == true)
	{
		this->Set_Reverse(this->reverse_ref_rate * pow(T, this->temperature_affinity) * exp(-this->activation_energy/(Rstd*T)));
		return;
	}
	else if (this->HaveForward == true)
	{
		this->reverse_rate = this->forward_rate / pow(10.0,this->Equilibrium);
		return;
	}
	else if (this->HaveReverse == true)
	{
		this->forward_rate = this->reverse_rate * pow(10.0,this->Equilibrium);
		return;
	}
	else
	{
		mError(missing_information);
		return;
	}
}

//Return the equilibrium status
bool UnsteadyReaction::haveEquilibrium()
{
	return this->HaveEquil;
}

//Return the reaction rate status
bool UnsteadyReaction::haveRate()
{
	if ((this->HaveReverse == true || this->HaveForward == true) && this->HaveEquil == true)
	{
		return true;
	}
	else if ((this->HaveForRef == true || this->HaveRevRef == true) && this->HaveEquil == true)
	{
		return true;
	}
	else
		return false;
}

//Return true if have forward reference
bool UnsteadyReaction::haveForwardRef()
{
	return this->HaveForRef;
}

//Return true if have reverse reference
bool UnsteadyReaction::haveReverseRef()
{
	return this->HaveRevRef;
}

//Return true if have forward rate
bool UnsteadyReaction::haveForward()
{
	return this->HaveForward;
}

//Return true if have reverse rate
bool UnsteadyReaction::haveReverse()
{
	return this->HaveReverse;
}

//Get the species index for the Unsteady species
int UnsteadyReaction::Get_Species_Index()
{
	return this->species_index;
}

//Get the stoichiometric constant
double UnsteadyReaction::Get_Stoichiometric (int i)
{
	double Sto = 0.0;

	//Check args
	if (i >= this->Stoichiometric.size())
	{
		mError(out_of_bounds);
		return 0.0;
	}

	Sto = this->Stoichiometric[i];
	return Sto;
}

//Get the equilibrium constant
double UnsteadyReaction::Get_Equilibrium ()
{
	return this->Equilibrium;
}

//Fetch the Enthalpy
double UnsteadyReaction::Get_Enthalpy()
{
	return this->enthalpy;
}

//Fetch the Entropy
double UnsteadyReaction::Get_Entropy()
{
	return this->entropy;
}

//Fetch the Energy
double UnsteadyReaction::Get_Energy()
{
	return this->energy;
}

//Get the intial value of the unsteady variable
double UnsteadyReaction::Get_InitialValue()
{
	return this->initial_value;
}

//Get the maximum concentration value
double UnsteadyReaction::Get_MaximumValue()
{
	return this->max_value;
}

//Get forward rate
double UnsteadyReaction::Get_Forward()
{
	return this->forward_rate;
}

//Get reverse rate
double UnsteadyReaction::Get_Reverse()
{
	return this->reverse_rate;
}

//Get the forward referenece rate
double UnsteadyReaction::Get_ForwardRef()
{
	return this->forward_ref_rate;
}

//Get the reverse reference rate
double UnsteadyReaction::Get_ReverseRef()
{
	return this->reverse_ref_rate;
}

//Get the activation energy
double UnsteadyReaction::Get_ActivationEnergy()
{
	return this->activation_energy;
}

//Get the temperature affinity
double UnsteadyReaction::Get_Affinity()
{
	return this->temperature_affinity;
}

//Get time step
double UnsteadyReaction::Get_TimeStep()
{
	return this->time_step;
}

//Function to calculate the reaction rate from given concentration and activities
double UnsteadyReaction::Eval_ReactionRate(const Matrix<double> &x, const Matrix<double> &gama)
{
	double R = 0.0;

	//Loop over all species in list
	double reactants = 0.0, products = 0.0;
	bool first_prod = true, first_reac = true;
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (this->Stoichiometric[i] > 0.0)
		{
			if (first_prod == true)
			{
				products = ( pow(gama(i,0),fabs(this->Stoichiometric[i])) * pow(10.0,(fabs(this->Stoichiometric[i])*x(i,0)) ) );
				first_prod = false;
			}
			else
				products = products * ( pow(gama(i,0),fabs(this->Stoichiometric[i])) * pow(10.0,(fabs(this->Stoichiometric[i])*x(i,0)) ) );
		}
		else if (this->Stoichiometric[i] < 0.0)
		{
			if (first_reac == true)
			{
				reactants = ( pow(gama(i,0),fabs(this->Stoichiometric[i])) * pow(10.0,(fabs(this->Stoichiometric[i])*x(i,0)) ) );
				first_reac = false;
			}
			else
				reactants = reactants * ( pow(gama(i,0),fabs(this->Stoichiometric[i])) * pow(10.0,(fabs(this->Stoichiometric[i])*x(i,0)) ) );
		}
		else
		{
			//No action
		}
	}
	R = fabs(this->Get_Stoichiometric(this->Get_Species_Index())) * ((this->Get_Forward() * reactants) - (this->Get_Reverse() * products));

	return R;
}

//Function to calculate the unsteady reaction residual contribution
double UnsteadyReaction::Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old, const Matrix<double> &gama_new, const Matrix<double> &gama_old)
{
	double res = 0.0, rate;
	double step, log_step;

	//Take full implicit time step
	rate = this->Eval_ReactionRate(x_new, gama_new);
	step = ( (gama_old(this->species_index,0) * pow(10.0, x_old(this->species_index,0))) + (this->time_step * rate) );

	if (step <= 0.0)
		step = DBL_MIN;
	log_step = log10(step);
	if (log_step >= log10(this->Get_MaximumValue()))
		res = this->Reaction::Eval_Residual(x_new, gama_new);
	else
		res = log10(gama_new(this->species_index,0)) + x_new(this->species_index,0) - log_step;

	if (isnan(res) || isinf(res))
		res = sqrt(DBL_MAX)/this->List->list_size();

	return res;
}

//Call the steady-state version of the residual function
double UnsteadyReaction::Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama)
{
	return this->Reaction::Eval_Residual(x, gama);
}

//Calculate the residual for the initial condtions
double UnsteadyReaction::Eval_IC_Residual(const Matrix<double> &x)
{
	double ic = 0.0;
	if (this->Get_InitialValue() <= 0.0)
		ic = -DBL_MAX;
	else
		ic = log10(this->Get_InitialValue());

	return x(this->Get_Species_Index(),0) - ic;
}

//Function to estimate the new concentration by an explicit step
double UnsteadyReaction::Explicit_Eval(const Matrix<double> &x, const Matrix<double> &gama)
{
	double conc = pow(10.0,x(this->Get_Species_Index(),0))+(this->Get_TimeStep()*this->Eval_ReactionRate(x, gama)/gama(this->Get_Species_Index(),0));
	if (conc <= 0.0)
		conc = DBL_MIN;
	if (conc > this->Get_MaximumValue())
		conc = this->Get_MaximumValue() * 0.90;
	return conc;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: UnsteadyReaction
 */

/*
 *								Start: AdsorptionReaction
 *	-------------------------------------------------------------------------------------
 */

//Default Constructor for Adsorption Reaction
AdsorptionReaction::AdsorptionReaction()
:
ads_rxn(0),
area_factors(0),
volume_factors(0),
adsorb_index(0),
molar_factor(0),
activities()
{
	surface_activity = (*ideal_solution);
	activity_data = nullptr;
	act_fun = IDEAL_ADS;
	specific_area = 1.0;
	specific_molality = 1.0;
	surface_charge = 0.0;
	total_mass = 0.0;
	total_volume = 1.0;
	charge_density = 0.0;
	ionic_strength = 0.0;
	num_rxns = 0;
	AreaBasis = true;
	IncludeSurfCharge = true;
	adsorbent_name = "AX";
}

//Default destructor for Adsorption Reaction
AdsorptionReaction::~AdsorptionReaction()
{
	ads_rxn.clear();
	area_factors.clear();
	volume_factors.clear();
	adsorb_index.clear();
}

//Initialization of the object
void AdsorptionReaction::Initialize_Object(MasterSpeciesList &List, int n)
{
	this->List = &List;
	this->num_rxns = n;
	this->area_factors.resize(this->List->list_size());
	this->volume_factors.resize(this->List->list_size());
	this->activities.set_size(this->List->list_size(),1);
	for (int i=0; i<this->List->list_size(); i++)
	{
		this->area_factors[i] = 0.0;
		this->volume_factors[i] = 0.0;
		this->activities(i,0) = 1.0;
	}

	this->ads_rxn.resize(this->num_rxns);
	this->adsorb_index.resize(this->num_rxns);
	this->aqueous_index.resize(this->num_rxns);
	this->molar_factor.resize(this->num_rxns);
	for (int i=0; i<this->num_rxns; i++)
	{
		this->ads_rxn[i].Initialize_Object(List);
		this->adsorb_index[i] = -1;
		this->aqueous_index[i] = -1;
		this->molar_factor[i] = 1.0;
	}
}

//Display info about the AdsorptionReaction Object (PLACE HOLDER)
void AdsorptionReaction::Display_Info()
{
	for (int i=0; i<this->ads_rxn.size(); i++)
		this->ads_rxn[i].Display_Info();
}

//Find and set the adsorption indices
int AdsorptionReaction::setAdsorbIndices()
{
	int success = 0;
	for (int i=0; i<this->num_rxns; i++)
	{
		for (int n=0; n<this->List->list_size(); n++)
		{
			if (this->ads_rxn[i].Get_Stoichiometric(n) != 0.0 && (this->List->get_species(n).MoleculePhaseID() == SOLID || this->List->get_species(n).MoleculePhaseID() == ADSORBED))
			{
				this->adsorb_index[i] = n;
				break;
			}
		}
		double normal_factor = this->ads_rxn[i].Get_Stoichiometric(this->adsorb_index[i]);
		for (int n=0; n<this->List->list_size(); n++)
		{
			this->ads_rxn[i].Set_Stoichiometric(n, this->ads_rxn[i].Get_Stoichiometric(n)/normal_factor);
		}
		if (this->adsorb_index[i] < 0)
		{
			mError(invalid_species);
			return -1;
		}
	}
	return success;
}

//Function to check for indexing errors in aqueous indices
int AdsorptionReaction::checkAqueousIndices()
{
	int success = 0;
	for (int i=0; i<this->num_rxns; i++)
	{
		if (this->aqueous_index[i] >= this->List->list_size() || this->aqueous_index[i] < 0)
		{
			mError(indexing_error);
			return -1;
		}
		if (this->getReaction(i).Get_Stoichiometric(this->aqueous_index[i]) > 0.0)
		{
			if (this->getReaction(i).Get_Stoichiometric(this->adsorb_index[i]) > 0.0)
			{
				mError(indexing_error);
				return -1;
			}
		}
		if (this->getReaction(i).Get_Stoichiometric(this->aqueous_index[i]) < 0.0)
		{
			if (this->getReaction(i).Get_Stoichiometric(this->adsorb_index[i]) < 0.0)
			{
				mError(indexing_error);
				return -1;
			}
		}
	}
	return success;
}

//Set the activity function and data structure
void AdsorptionReaction::setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
											  const void *act_data)
{
	if ( (*act) == NULL )
	{
		this->surface_activity = (*ideal_solution);
	}
	else
	{
		this->surface_activity = (*act);
	}
	if ( (act_data) == NULL	)
	{
		this->activity_data = this;
	}
	else
	{
		this->activity_data = act_data;
	}
}

//Set the aqueous species index
void AdsorptionReaction::setAqueousIndex(int rxn, int sp)
{
	if (rxn >= this->num_rxns || rxn < 0)
	{
		mError(out_of_bounds);
		return;
	}
	if (sp >= this->List->list_size() || sp < 0)
	{
		mError(out_of_bounds);
		return;
	}
	this->aqueous_index[rxn] = sp;
}

//Automatically set the aqueous species index
int AdsorptionReaction::setAqueousIndexAuto()
{
	int success = 0;
	for (int i=0; i<this->num_rxns; i++)
	{
		if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) > 0.0)
		{
			for (int n=0; n<this->List->list_size(); n++)
			{
				if (this->ads_rxn[i].Get_Stoichiometric(n) < 0.0 && this->List->get_species(n).MoleculePhaseID() == AQUEOUS	)
				{
					this->setAqueousIndex(i, n);
					break;
				}
			}
		}
		else if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) < 0.0)
		{
			for (int n=0; n<this->List->list_size(); n++)
			{
				if (this->ads_rxn[i].Get_Stoichiometric(n) > 0.0 && this->List->get_species(n).MoleculePhaseID() == AQUEOUS	)
				{
					this->setAqueousIndex(i, n);
					break;
				}
			}
		}
		else
		{
			mError(invalid_species);
			return -1;
		}
	}
	return success;
}

//Set the enum for surface activity
void AdsorptionReaction::setActivityEnum(int act)
{
	switch (act)
	{
		case IDEAL_ADS:
			this->act_fun = IDEAL_ADS;
			break;

		case FLORY_HUGGINS:
			this->act_fun = FLORY_HUGGINS;
			break;

		case UNIQUAC_ACT:
			this->act_fun = UNIQUAC_ACT;
			break;

		default:
			this->act_fun = IDEAL_ADS;
			break;
	}
}

//Set the molar factor for the given reaction
void AdsorptionReaction::setMolarFactor(int rxn_i, double m)
{
	if (rxn_i < 0 || rxn_i >= this->num_rxns)
	{
		mError(out_of_bounds);
		return;
	}
	this->molar_factor[rxn_i] = m;
}

//Set the ith volume factor for the species
void AdsorptionReaction::setVolumeFactor(int i, double v)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		return;
	}
	if (v < 0)
		v = DBL_MIN;
	this->volume_factors[i] = v;
}

//Set the ith area factor for the species
void AdsorptionReaction::setAreaFactor(int i, double a)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		return;
	}
	if (a < 0)
		a = DBL_MIN;
	this->area_factors[i] = a;
}

//Set the specific area for the adsorbent
void AdsorptionReaction::setSpecificArea(double a)
{
	if (a < 0)
		a = DBL_MIN;
	this->specific_area = a;
}

//Set the specific molality for the adsorbent
void AdsorptionReaction::setSpecificMolality(double a)
{
	if (a < 0)
		a = DBL_MIN;
	this->specific_molality = a;
}

//Set the surface charge
void AdsorptionReaction::setSurfaceCharge(double c)
{
	this->surface_charge = c;
}

//Set the total mass of adsorbent in system
void AdsorptionReaction::setTotalMass(double m)
{
	if (m < 0)
		m = DBL_MIN;
	this->total_mass = m;
}

//Set the total volume of the system
void AdsorptionReaction::setTotalVolume(double v)
{
	if (v < 0)
		v = DBL_MIN;
	this->total_volume = v;
}

//Directly set the area basis boolean
void AdsorptionReaction::setAreaBasisBool(bool opt)
{
	this->AreaBasis = opt;
}

//Directly set the surface charging boolean
void AdsorptionReaction::setSurfaceChargeBool(bool opt)
{
	this->IncludeSurfCharge = opt;
}

//Set the area basis from the string argument
void AdsorptionReaction::setBasis(std::string option)
{
	option = allLower(option);
	if (option == "area" || option == "surface")
	{
		this->setAreaBasisBool(true);
	}
	else if (option == "ligand" || option == "molar" || option == "mole")
	{
		this->setAreaBasisBool(false);
	}
	else
	{
		mError(invalid_type);
		this->setAreaBasisBool(true);
	}
}

//Set the name of the adsorbent
void AdsorptionReaction::setAdsorbentName(std::string name)
{
	this->adsorbent_name = name;
}

//Function to set the charge density to a given value
void AdsorptionReaction::setChargeDensityValue(double a)
{
	this->charge_density = a;
}

//Function to set ionic strength to a given value
void AdsorptionReaction::setIonicStrengthValue(double a)
{
	this->ionic_strength = a;
}

//Function to set the activity matrix to a specific matrix
void AdsorptionReaction::setActivities(Matrix<double> &x)
{
	if (x.rows() != this->activities.rows() || x.columns() != this->activities.columns())
	{
		mError(dim_mis_match);
		return;
	}

	this->activities = x;
}

//Modify the deltas in the given mass balance
void AdsorptionReaction::modifyDeltas(MassBalance &mbo)
{
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED)
		{
			mbo.Set_Delta(i, (mbo.Get_Delta(i) * this->getBulkDensity()) );
		}
	}
}

//Calculation and initialization of the area factors
void AdsorptionReaction::calculateAreaFactors()
{
	for (int r=0; r<this->getNumberRxns(); r++)
	{
		for (int i=0; i<this->List->list_size(); i++)
		{

			if (this->volume_factors[i] == 0.0 || (this->List->get_species(i).MoleculePhaseID() != SOLID && this->List->get_species(i).MoleculePhaseID() != ADSORBED))
			{
				this->area_factors[i] = 0.0;
			}
			else if (this->volume_factors[i] == 0.0 && (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED))
			{
				this->area_factors[i] = 4.0 * M_PI * pow((3.0/(4.0*M_PI))*(4.33/Na), (2.0/3.0)) / 10000.0 * Na;
			}
			else
			{
				if (this->AreaBasis == false)
				{
					this->volume_factors[i] = this->volume_factors[i] * this->molar_factor[r];
				}
				this->area_factors[i] = 4.0 * M_PI * pow((3.0/(4.0*M_PI))*(this->volume_factors[i]/Na), (2.0/3.0)) / 10000.0 * Na;
			}
		}
	}
}

//Calculation of equilibria parameters from temperature
void AdsorptionReaction::calculateEquilibria(double T)
{
	for (int i=0; i<this->num_rxns; i++)
		this->getReaction(i).calculateEquilibrium(T);
}

//Set the charge density
void AdsorptionReaction::setChargeDensity(const Matrix<double> &x)
{
	this->charge_density = this->calculateSurfaceChargeDensity(x);
}

//Set the Ionic strength of the solution
void AdsorptionReaction::setIonicStrength(const Matrix<double> &x)
{
	this->ionic_strength = calculate_ionic_strength(x, *this->List);
}

//Call the activity model passing the logx concentrations of species
int AdsorptionReaction::callSurfaceActivity(const Matrix<double> &x)
{
	int success = this->surface_activity(x,this->activities,this->activity_data);
	return success;
}

//Calculation of the surface charge density
double AdsorptionReaction::calculateSurfaceChargeDensity( const Matrix<double> &x)
{
	double sigma = 0.0;
	double sum = 0.0;

	for (int i=0; i<this->getNumberRxns(); i++)
	{
		sum = sum + (this->List->charge(this->getAdsorbIndex(i)) * pow(10.0, x(this->getAdsorbIndex(i),0)));
	}
	sum = sum + (this->getSurfaceCharge() * this->getSpecificMolality() * this->calculateActiveFraction(x));

	sigma = sum * (Faraday/this->getSpecificArea());

	return sigma;
}

//Calculation of the active fraction of the surface area
double AdsorptionReaction::calculateActiveFraction(const Matrix<double> &x)
{
	double phi = 0.0, sum = 0.0;

	if (this->isAreaBasis() == true)
	{
		for (int i=0; i<this->List->list_size(); i++)
		{
			sum = sum + (this->getAreaFactor(i) * pow(10.0, x(i,0)));
		}
		if (sum > this->specific_area)
			sum = this->specific_area;
		phi = 1.0 - (sum / (this->specific_area));
		if (phi < DBL_MIN)
			phi = DBL_MIN;
	}
	else
	{
		for (int i=0; i<this->num_rxns; i++)
		{
			sum = sum + (this->getMolarFactor(i) * pow(10.0, x(this->getAdsorbIndex(i),0)));
		}
		if (sum > this->specific_molality)
			sum = this->specific_molality;
		phi = 1.0 - (sum / (this->specific_molality));
		if (phi < DBL_MIN)
			phi = DBL_MIN;
	}

	return phi;
}

//Calculation of maximum capacity for adsorption (may vary with concentrations)
double AdsorptionReaction::calculateLangmuirMaxCapacity(int i)
{
	if (i < 0 || i >= this->num_rxns)
	{
		mError(out_of_bounds);
		return 0.0;
	}
	if (this->isAreaBasis() == true)
		return (this->specific_area / (this->getAreaFactor(this->getAdsorbIndex(i))));
	else
		return (this->specific_molality / this->molar_factor[i]);
}

//Calculation of the effective Langmuir parameter for the ith adsorption reaction
double AdsorptionReaction::calculateLangmuirEquParam(const Matrix<double> &x, const Matrix<double> &gama, int i)
{
	if (i < 0 || i >= this->num_rxns)
	{
		mError(out_of_bounds);
		return 1.0;
	}
	double K = pow(10.0,this->getReaction(i).Get_Equilibrium());
	double prod = 1.0, react = 1.0;

	for (int n=0; n<this->List->list_size(); n++)
	{
		if (n != this->getAqueousIndex(i) && this->List->get_species(n).MoleculePhaseID() != SOLID && this->List->get_species(n).MoleculePhaseID() != ADSORBED)
		{
			if (this->getReaction(i).Get_Stoichiometric(n) > 0.0)
			{
				//Products
				prod = prod * ( pow(gama(n,0),fabs(this->getReaction(i).Get_Stoichiometric(n))) * pow(10.0, fabs(this->getReaction(i).Get_Stoichiometric(n))*x(n,0)));
			}
			else if (this->getReaction(i).Get_Stoichiometric(n) < 0.0)
			{
				//Reactants
				react = react * ( pow(gama(n,0),fabs(this->getReaction(i).Get_Stoichiometric(n))) * pow(10.0, fabs(this->getReaction(i).Get_Stoichiometric(n))*x(n,0)));
			}
			else
			{
				//No Action
			}
		}
	}

	if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) > 0.0)
	{
		if (this->isAreaBasis() == true)
			K = K * this->getAreaFactor(this->getAdsorbIndex(i)) * pow(gama(this->getAqueousIndex(i),0),fabs(this->getReaction(i).Get_Stoichiometric(this->getAqueousIndex(i))) ) * (react/prod) / this->getActivity(this->getAdsorbIndex(i));
		else
			K = K * this->getMolarFactor(i) * pow(gama(this->getAqueousIndex(i),0),fabs(this->getReaction(i).Get_Stoichiometric(this->getAqueousIndex(i))) ) * (react/prod) / this->getActivity(this->getAdsorbIndex(i));
	}
	else
	{
		if (this->isAreaBasis() == true)
			K = 1.0 / (K * this->getAreaFactor(this->getAdsorbIndex(i)) * pow(gama(this->getAqueousIndex(i),0),fabs(this->getReaction(i).Get_Stoichiometric(this->getAqueousIndex(i))) ) * (react/prod) / this->getActivity(this->getAdsorbIndex(i)) );
		else
			K = 1.0 / (K * this->getMolarFactor(i) * pow(gama(this->getAqueousIndex(i),0),fabs(this->getReaction(i).Get_Stoichiometric(this->getAqueousIndex(i))) ) * (react/prod) / this->getActivity(this->getAdsorbIndex(i)) );
	}

	return K;
}

//Calculation of the Langmuir Adsorption for the given reaction
double AdsorptionReaction::calculateLangmuirAdsorption(const Matrix<double> &x, const Matrix<double> &gama, int i)
{
	double q = 0.0;
	double qmax, top, bot = 1.0;
	qmax = this->calculateLangmuirMaxCapacity(i);
	top = this->calculateLangmuirEquParam(x, gama, i) * pow(10.0, fabs(this->getReaction(i).Get_Stoichiometric(this->getAqueousIndex(i)))*x(this->getAqueousIndex(i),0));
	for (int n=0; n<this->num_rxns; n++)
	{
		bot = bot + ( this->calculateLangmuirEquParam(x, gama, n) * pow(10.0, fabs(this->getReaction(n).Get_Stoichiometric(this->getAqueousIndex(n)))*x(this->getAqueousIndex(n),0)));
	}
	q = (qmax * (top / bot));
	return q;
}

//Approximation of the electric potential of the surface
double AdsorptionReaction::calculatePsi(double sigma, double T, double I, double rel_epsilon)
{
	double psi = 0.0;

	I = I * 1000.0;//First, convert ionic strength to mol/m^3
	double coeff = sigma / sqrt(8.0*AbsPerm(rel_epsilon)*Rstd*T*I);
	double inv = log(coeff + sqrt((coeff*coeff)+1.0));
	psi = (2.0 * kB * T * inv) / e;

	return psi; //Units: V
}

//Calculation of charge exchange in given reaction
double AdsorptionReaction::calculateAqueousChargeExchange(int i)
{
	double n = 0.0;

	for (int j = 0; j<this->List->list_size(); j++)
	{
		if (this->List->get_species(j).MoleculePhaseID() == AQUEOUS || this->List->get_species(j).MoleculePhaseID() == LIQUID)
		{
			n = n + (this->getReaction(i).Get_Stoichiometric(j)*this->List->charge(j));
		}
	}

	return -n;
}

//Calculation of equilibrium correct term
double AdsorptionReaction::calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int i)
{
	if (this->includeSurfaceCharge() == true)
		return -(this->calculateAqueousChargeExchange(i)*e*calculatePsi(sigma, T, I, rel_epsilon))/(kB*T);
	else
		return 0.0;
}

//Calculation of residual of ith reaction for the solver to work on
double AdsorptionReaction::Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i)
{
	double res = this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) * ( log10(this->getActivity(this->getAdsorbIndex(i))) + x(this->getAdsorbIndex(i),0) );

	if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) > 0.0)
	{
		if (this->isAreaBasis() == true)
			res = res - log10(this->getSpecificArea()*this->calculateActiveFraction(x));
		else
			res = res - (this->getMolarFactor(i) * log10(this->getSpecificMolality()*this->calculateActiveFraction(x)));
	}
	else
	{
		if (this->isAreaBasis() == true)
			res = res + log10(this->getSpecificArea()*this->calculateActiveFraction(x));
		else
			res = res + (this->getMolarFactor(i) * log10(this->getSpecificMolality()*this->calculateActiveFraction(x)));
	}

	for (int n=0; n<this->List->list_size(); n++)
	{
		if (this->List->get_species(n).MoleculePhaseID() == AQUEOUS || this->List->get_species(n).MoleculePhaseID() == LIQUID)
			res = res + ( this->getReaction(i).Get_Stoichiometric(n)*( log10(gama(n,0))+x(n,0) ) );
	}

	double logK = this->getReaction(i).Get_Equilibrium();
	logK = logK + ((this->calculateEquilibriumCorrection(this->getChargeDensity(), T, this->getIonicStrength(), rel_perm, i))/log(10.0));

	res = res - logK;

	return res;
}

//Return reference to the ith reaction object
Reaction& AdsorptionReaction::getReaction(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->ads_rxn[i];
}

//Return the ith molar factor
double AdsorptionReaction::getMolarFactor(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->molar_factor[i];
}

//Return the ith volume factor
double AdsorptionReaction::getVolumeFactor(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->volume_factors[i];
}

//Return the ith area factor
double AdsorptionReaction::getAreaFactor(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->area_factors[i];
}

//Return the ith activity coefficient
double AdsorptionReaction::getActivity(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->activities(i,0);
}

//Return the specific molality
double AdsorptionReaction::getSpecificMolality()
{
	return this->specific_molality;
}

//Return the specific area
double AdsorptionReaction::getSpecificArea()
{
	return this->specific_area;
}

//Return the surface charge
double AdsorptionReaction::getSurfaceCharge()
{
	return this->surface_charge;
}

//Calculate and return the bulk density of the adsorbent in the system
double AdsorptionReaction::getBulkDensity()
{
	return this->total_mass / this->total_volume;
}

//Return the total mass
double AdsorptionReaction::getTotalMass()
{
	return this->total_mass;
}

//Return the total system volume
double AdsorptionReaction::getTotalVolume()
{
	return this->total_volume;
}

//Return charge density
double AdsorptionReaction::getChargeDensity()
{
	return this->charge_density;
}

//Return the ionic strength
double AdsorptionReaction::getIonicStrength()
{
	return this->ionic_strength;
}

//Return the number of reactions
int AdsorptionReaction::getNumberRxns()
{
	return this->num_rxns;
}

//Return index of the adsorbed species
int AdsorptionReaction::getAdsorbIndex(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	if (this->adsorb_index[i] >= this->List->list_size() || this->adsorb_index[i] < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	return this->adsorb_index[i];
}

//Return the index of the primary aqueous species
int AdsorptionReaction::getAqueousIndex(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	if (this->aqueous_index[i] >= this->List->list_size() || this->aqueous_index[i] < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	return this->aqueous_index[i];
}

//Return the activity enum
int AdsorptionReaction::getActivityEnum()
{
	return this->act_fun;
}

//Return true is in Area basis
bool AdsorptionReaction::isAreaBasis()
{
	return this->AreaBasis;
}

//Return state of surface charge inclusion
bool AdsorptionReaction::includeSurfaceCharge()
{
	return this->IncludeSurfCharge;
}

//Return the name of the adsorbent
std::string AdsorptionReaction::getAdsorbentName()
{
	return this->adsorbent_name;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: AdsorptionReaction
 */

/*
 *								Start: UnsteadyAdsorption
 *	-------------------------------------------------------------------------------------
 */
//Default Constructor
UnsteadyAdsorption::UnsteadyAdsorption()
{
	ads_rxn.resize(0);
	area_factors.resize(0);
	volume_factors.resize(0);
	adsorb_index.resize(0);
	molar_factor.resize(0);
	surface_activity = (*ideal_solution);
	activity_data = nullptr;
	act_fun = IDEAL_ADS;
	specific_area = 1.0;
	specific_molality = 1.0;
	surface_charge = 0.0;
	total_mass = 0.0;
	total_volume = 1.0;
	charge_density = 0.0;
	ionic_strength = 0.0;
	num_rxns = 0;
	AreaBasis = true;
	IncludeSurfCharge = true;
	adsorbent_name = "AX";
}

//Default Destructor
UnsteadyAdsorption::~UnsteadyAdsorption()
{
	ads_rxn.clear();
	area_factors.clear();
	volume_factors.clear();
	adsorb_index.clear();
}

//Initialize memory for unsteady adsorption
void UnsteadyAdsorption::Initialize_Object(MasterSpeciesList &List, int n)
{
	this->List = &List;
	this->num_rxns = n;
	this->area_factors.resize(this->List->list_size());
	this->volume_factors.resize(this->List->list_size());
	this->activities.set_size(this->List->list_size(),1);
	this->activities_old.set_size(this->List->list_size(),1);
	for (int i=0; i<this->List->list_size(); i++)
	{
		this->area_factors[i] = 0.0;
		this->volume_factors[i] = 0.0;
		this->activities(i,0) = 1.0;
	}

	this->ads_rxn.resize(this->num_rxns);
	this->adsorb_index.resize(this->num_rxns);
	this->aqueous_index.resize(this->num_rxns);
	this->molar_factor.resize(this->num_rxns);
	for (int i=0; i<this->num_rxns; i++)
	{
		this->ads_rxn[i].Initialize_Object(List);
		this->adsorb_index[i] = -1;
		this->aqueous_index[i] = -1;
		this->molar_factor[i] = 1.0;
	}
}

//Display object information
void UnsteadyAdsorption::Display_Info()
{
	for (int i=0; i<this->ads_rxn.size(); i++)
		this->ads_rxn[i].Display_Info();
}

//Modify the deltas in the given mass balance
void UnsteadyAdsorption::modifyDeltas(MassBalance &mbo)
{
	this->AdsorptionReaction::modifyDeltas(mbo);
}

//Find and set the adsorption indices
int UnsteadyAdsorption::setAdsorbIndices()
{
	int success = 0;
	for (int i=0; i<this->num_rxns; i++)
	{
		for (int n=0; n<this->List->list_size(); n++)
		{
			if (this->ads_rxn[i].Get_Stoichiometric(n) != 0.0 && (this->List->get_species(n).MoleculePhaseID() == SOLID || this->List->get_species(n).MoleculePhaseID() == ADSORBED))
			{
				this->adsorb_index[i] = n;
				this->ads_rxn[i].Set_Species_Index(n);
				break;
			}
		}
		double normal_factor = this->ads_rxn[i].Get_Stoichiometric(this->adsorb_index[i]);
		for (int n=0; n<this->List->list_size(); n++)
		{
			this->ads_rxn[i].Set_Stoichiometric(n, this->ads_rxn[i].Get_Stoichiometric(n)/normal_factor);
		}
		if (this->adsorb_index[i] < 0)
		{
			mError(invalid_species);
			return -1;
		}
	}
	return success;
}

//Function to check for indexing errors in aqueous indices
int UnsteadyAdsorption::checkAqueousIndices()
{
	int success = 0;
	for (int i=0; i<this->num_rxns; i++)
	{
		if (this->aqueous_index[i] >= this->List->list_size() || this->aqueous_index[i] < 0)
		{
			mError(indexing_error);
			return -1;
		}
		if (this->getReaction(i).Get_Stoichiometric(this->aqueous_index[i]) > 0.0)
		{
			if (this->getReaction(i).Get_Stoichiometric(this->adsorb_index[i]) > 0.0)
			{
				mError(indexing_error);
				return -1;
			}
		}
		if (this->getReaction(i).Get_Stoichiometric(this->aqueous_index[i]) < 0.0)
		{
			if (this->getReaction(i).Get_Stoichiometric(this->adsorb_index[i]) < 0.0)
			{
				mError(indexing_error);
				return -1;
			}
		}
	}
	return success;
}

//Set the activity function and data structure
void UnsteadyAdsorption::setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
											  const void *act_data)
{
	this->AdsorptionReaction::setActivityModelInfo(act,act_data);
}

//Set the aqueous species index
void UnsteadyAdsorption::setAqueousIndex(int rxn_i, int species_i)
{
	this->AdsorptionReaction::setAqueousIndex(rxn_i, species_i);
}

//Automatically set the aqueous species index
int UnsteadyAdsorption::setAqueousIndexAuto()
{
	int success = 0;
	for (int i=0; i<this->num_rxns; i++)
	{
		if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) > 0.0)
		{
			for (int n=0; n<this->List->list_size(); n++)
			{
				if (this->ads_rxn[i].Get_Stoichiometric(n) < 0.0 && this->List->get_species(n).MoleculePhaseID() == AQUEOUS	)
				{
					this->setAqueousIndex(i, n);
					break;
				}
			}
		}
		else if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) < 0.0)
		{
			for (int n=0; n<this->List->list_size(); n++)
			{
				if (this->ads_rxn[i].Get_Stoichiometric(n) > 0.0 && this->List->get_species(n).MoleculePhaseID() == AQUEOUS	)
				{
					this->setAqueousIndex(i, n);
					break;
				}
			}
		}
		else
		{
			mError(invalid_species);
			return -1;
		}
	}
	return success;
}

//Set the enum for surface activity
void UnsteadyAdsorption::setActivityEnum(int act)
{
	switch (act)
	{
		case IDEAL_ADS:
			this->act_fun = IDEAL_ADS;
			break;

		case FLORY_HUGGINS:
			this->act_fun = FLORY_HUGGINS;
			break;

		case UNIQUAC_ACT:
			this->act_fun = UNIQUAC_ACT;
			break;

		default:
			this->act_fun = IDEAL_ADS;
			break;
	}
}

//Set the molar factor for the given reaction
void UnsteadyAdsorption::setMolarFactor(int rxn_i, double m)
{
	this->AdsorptionReaction::setMolarFactor(rxn_i,m);
}

//Set the ith volume factor for the species
void UnsteadyAdsorption::setVolumeFactor(int i, double v)
{
	this->AdsorptionReaction::setVolumeFactor(i,v);
}

//Set the ith area factor for the species
void UnsteadyAdsorption::setAreaFactor(int i, double a)
{
	this->AdsorptionReaction::setAreaFactor(i,a);
}

//Set the specific area for the adsorbent
void UnsteadyAdsorption::setSpecificArea(double a)
{
	this->AdsorptionReaction::setSpecificArea(a);
}

//Set the specific molality for the adsorbent
void UnsteadyAdsorption::setSpecificMolality(double a)
{
	this->AdsorptionReaction::setSpecificMolality(a);
}

//Set the surface charge
void UnsteadyAdsorption::setSurfaceCharge(double c)
{
	this->AdsorptionReaction::setSurfaceCharge(c);
}

//Set the total mass of adsorbent in system
void UnsteadyAdsorption::setTotalMass(double m)
{
	this->AdsorptionReaction::setTotalMass(m);
}

//Set the total volume of the system
void UnsteadyAdsorption::setTotalVolume(double v)
{
	this->AdsorptionReaction::setTotalVolume(v);
}

//Directly set the area basis boolean
void UnsteadyAdsorption::setAreaBasisBool(bool opt)
{
	this->AdsorptionReaction::setAreaBasisBool(opt);
}

//Directly set the surface charging boolean
void UnsteadyAdsorption::setSurfaceChargeBool(bool opt)
{
	this->AdsorptionReaction::setSurfaceChargeBool(opt);
}

//Set the area basis from the string argument
void UnsteadyAdsorption::setBasis(std::string option)
{
	this->AdsorptionReaction::setBasis(option);
}

//Set the name of the adsorbent
void UnsteadyAdsorption::setAdsorbentName(std::string name)
{
	this->AdsorptionReaction::setAdsorbentName(name);
}

//Update activities for next time step
void UnsteadyAdsorption::updateActivities()
{
	this->activities_old = this->activities;
}

//Return reference to the corresponding reaction
UnsteadyReaction& UnsteadyAdsorption::getReaction(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->ads_rxn[i];
}

//Return the ith molar factor
double UnsteadyAdsorption::getMolarFactor(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->molar_factor[i];
}

//Return the ith volume factor
double UnsteadyAdsorption::getVolumeFactor(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->volume_factors[i];
}

//Calculation and initialization of the area factors
void UnsteadyAdsorption::calculateAreaFactors()
{
	this->AdsorptionReaction::calculateAreaFactors();
}

//Calculate equilibria parameters in each reaction
void UnsteadyAdsorption::calculateEquilibria(double T)
{
	for (int i=0; i<this->num_rxns; i++)
		this->getReaction(i).calculateEquilibrium(T);
}

//Calculate reaction rate parameters in each reaction
void UnsteadyAdsorption::calculateRates(double T)
{
	for (int i=0; i<this->num_rxns; i++)
		this->getReaction(i).calculateRate(T);
}

//Set the charge density
void UnsteadyAdsorption::setChargeDensity(const Matrix<double> &x)
{
	this->charge_density = this->calculateSurfaceChargeDensity(x);
}

//Set the Ionic strength of the solution
void UnsteadyAdsorption::setIonicStrength(const Matrix<double> &x)
{
	this->AdsorptionReaction::setIonicStrength(x);
}

//Call the activity model passing the logx concentrations of species
int UnsteadyAdsorption::callSurfaceActivity(const Matrix<double> &x)
{
	int success = this->surface_activity(x,this->activities,this->activity_data);
	return success;
}

//Calculation of the active fraction of the surface area
double UnsteadyAdsorption::calculateActiveFraction(const Matrix<double> &x)
{
	double phi = 0.0, sum = 0.0;

	if (this->isAreaBasis() == true)
	{
		for (int i=0; i<this->List->list_size(); i++)
		{
			sum = sum + (this->getAreaFactor(i) * pow(10.0, x(i,0)));
		}
		if (sum > this->specific_area)
			sum = this->specific_area;
		phi = 1.0 - (sum / this->specific_area);
		if (phi < DBL_MIN)
			phi = DBL_MIN;
	}
	else
	{
		for (int i=0; i<this->num_rxns; i++)
		{
			sum = sum + (this->getMolarFactor(i) * pow(10.0, x(this->getAdsorbIndex(i),0)));
		}
		if (sum > this->specific_molality)
			sum = this->specific_molality;
		phi = 1.0 - (sum / this->specific_molality);
		if (phi < DBL_MIN)
			phi = DBL_MIN;
	}

	return phi;
}

//Calculation of the surface charge density
double UnsteadyAdsorption::calculateSurfaceChargeDensity( const Matrix<double> &x)
{
	double sigma = 0.0;
	double sum = 0.0;

	for (int i=0; i<this->getNumberRxns(); i++)
	{
		sum = sum + (this->List->charge(this->getAdsorbIndex(i)) * pow(10.0, x(this->getAdsorbIndex(i),0)));
	}
	sum = sum + (this->getSurfaceCharge() * this->getSpecificMolality() * this->calculateActiveFraction(x));

	sigma = sum * (Faraday/this->getSpecificArea());

	return sigma;
}

//Approximation of the electric potential of the surface
double UnsteadyAdsorption::calculatePsi(double sigma, double T, double I, double rel_epsilon)
{
	return this->AdsorptionReaction::calculatePsi(sigma,T,I,rel_epsilon);
}

//Calculation of charge exchange in given reaction
double UnsteadyAdsorption::calculateAqueousChargeExchange(int i)
{
	double n = 0.0;

	for (int j = 0; j<this->List->list_size(); j++)
	{
		if (this->List->get_species(j).MoleculePhaseID() == AQUEOUS || this->List->get_species(j).MoleculePhaseID() == LIQUID)
		{
			n = n + (this->getReaction(i).Get_Stoichiometric(j)*this->List->charge(j));
		}
	}

	return -n;
}

//Calculation of equilibrium correct term
double UnsteadyAdsorption::calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int i)
{
	if (this->includeSurfaceCharge() == true)
		return -(this->calculateAqueousChargeExchange(i)*e*calculatePsi(sigma, T, I, rel_epsilon))/(kB*T);
	else
		return 0.0;
}

//Calculation of residual of ith reaction for the solver to work on
double UnsteadyAdsorption::Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i)
{
	double res = this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) * ( log10(this->getActivity(this->getAdsorbIndex(i))) + x(this->getAdsorbIndex(i),0) );

	if (this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) > 0.0)
	{
		if (this->isAreaBasis() == true)
			res = res - log10(this->getSpecificArea()*this->calculateActiveFraction(x));
		else
			res = res - (this->getMolarFactor(i) * log10(this->getSpecificMolality()*this->calculateActiveFraction(x)));
	}
	else
	{
		if (this->isAreaBasis() == true)
			res = res + log10(this->getSpecificArea()*this->calculateActiveFraction(x));
		else
			res = res + (this->getMolarFactor(i) * log10(this->getSpecificMolality()*this->calculateActiveFraction(x)));
	}

	for (int n=0; n<this->List->list_size(); n++)
	{
		if (this->List->get_species(n).MoleculePhaseID() == AQUEOUS || this->List->get_species(n).MoleculePhaseID() == LIQUID)
			res = res + ( this->getReaction(i).Get_Stoichiometric(n)*( log10(gama(n,0))+x(n,0) ) );
	}

	double logK = this->getReaction(i).Get_Equilibrium();
	logK = logK + ((this->calculateEquilibriumCorrection(this->getChargeDensity(), T, this->getIonicStrength(), rel_perm, i))/log(10.0));
	res = res - logK;

	return res;
}

//Unsteady Adsorption residual
double UnsteadyAdsorption::Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old, const Matrix<double> &gama_new, const Matrix<double> &gama_old, double T, double rel_perm, int i)
{
	double res = 0.0, rate;
	double step, log_step;

	//Take full implicit time step
	rate = this->Eval_ReactionRate(x_new, gama_new, T, rel_perm, i);
	step = ( (this->getOldActivity(this->getReaction(i).Get_Species_Index()) * pow(10.0, x_old(this->getReaction(i).Get_Species_Index(),0))) + (this->getReaction(i).Get_TimeStep() * rate) );

	if (step <= 0.0)
		step = DBL_MIN;
	log_step = log10(step);
	if (log_step >= log10(this->getReaction(i).Get_MaximumValue()))
		res = this->UnsteadyAdsorption::Eval_Residual(x_new, gama_new, T, rel_perm, i);
	else
		res = log10(this->getActivity(this->getReaction(i).Get_Species_Index())) + x_new(this->getReaction(i).Get_Species_Index(),0) - log_step;

	if (isnan(res) || isinf(res))
		res = sqrt(DBL_MAX)/this->List->list_size();

	return res;
}

//Evaluation of the reaction rate
double UnsteadyAdsorption::Eval_ReactionRate(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int n)
{
	double R = 0.0;

	//Loop over all species in list
	double reactants = 0.0, products = 0.0;
	bool first_prod = true, first_reac = true;

	this->getReaction(n).Set_Equilibrium( this->getReaction(n).Get_Equilibrium() + ((this->calculateEquilibriumCorrection(this->getChargeDensity(), T, this->getIonicStrength(), rel_perm, n))/log(10.0)) );

	if (this->getReaction(n).haveForwardRef() == true)
	{
		this->getReaction(n).Set_Forward(this->getReaction(n).Get_ForwardRef() * pow(T, this->getReaction(n).Get_Affinity()) * exp(-this->getReaction(n).Get_ActivationEnergy()/(Rstd*T)));
	}
	else if (this->getReaction(n).haveReverseRef() == true)
	{
		this->getReaction(n).Set_Reverse(this->getReaction(n).Get_ReverseRef() * pow(T, this->getReaction(n).Get_Affinity()) * exp(-this->getReaction(n).Get_ActivationEnergy()/(Rstd*T)));
	}
	else if (this->getReaction(n).haveForward() == true)
	{
		this->getReaction(n).Set_Reverse( this->getReaction(n).Get_Forward() / pow(10.0,this->getReaction(n).Get_Equilibrium()) );
	}
	else if (this->getReaction(n).haveReverse() == true)
	{
		this->getReaction(n).Set_Forward( this->getReaction(n).Get_Reverse() * pow(10.0,this->getReaction(n).Get_Equilibrium()) );
	}
	else
	{
		//No Action
	}

	for (int i=0; i<this->List->list_size(); i++)
	{
		if (this->getReaction(n).Get_Stoichiometric(i) > 0.0)
		{
			if (first_prod == true)
			{
				if (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED)
					products = ( pow(this->getActivity(i),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
				else
					products = ( pow(gama(i,0),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
				first_prod = false;
			}
			else
			{
				if (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED)
					products = products * ( pow(this->getActivity(i),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
				else
					products = products * ( pow(gama(i,0),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
			}
		}
		else if (this->getReaction(n).Get_Stoichiometric(i) < 0.0)
		{
			if (first_reac == true)
			{
				if (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED)
					reactants = ( pow(this->getActivity(i),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
				else
					reactants = ( pow(gama(i,0),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
				first_reac = false;
			}
			else
			{
				if (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED)
					reactants = reactants * ( pow(this->getActivity(i),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
				else
					reactants = reactants * ( pow(gama(i,0),fabs(this->getReaction(n).Get_Stoichiometric(i))) * pow(10.0,(fabs(this->getReaction(n).Get_Stoichiometric(i))*x(i,0)) ) );
			}
		}
		else
		{
			//No action
		}
	}

	if (this->isAreaBasis() == true)
		reactants = reactants * this->getSpecificArea()*this->calculateActiveFraction(x);
	else
		reactants = reactants * pow(this->getSpecificMolality()*this->calculateActiveFraction(x),this->getMolarFactor(n));

	R = fabs(this->getReaction(n).Get_Stoichiometric(this->getReaction(n).Get_Species_Index())) * ((this->getReaction(n).Get_Forward() * reactants) - (this->getReaction(n).Get_Reverse() * products));

	return R;
}

//Calculate the residual for the initial condtions
double UnsteadyAdsorption::Eval_IC_Residual(const Matrix<double> &x, int i)
{
	double ic = 0.0;
	if (this->getReaction(i).Get_InitialValue() <= 0.0)
		ic = -DBL_MAX;
	else
		ic = log10(this->getReaction(i).Get_InitialValue());

	return x(this->getReaction(i).Get_Species_Index(),0) - ic;
}

//Function to estimate the new concentration by an explicit step
double UnsteadyAdsorption::Explicit_Eval(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i)
{
	double conc = pow(10.0,x(this->getReaction(i).Get_Species_Index(),0))+(this->getReaction(i).Get_TimeStep()*this->Eval_ReactionRate(x, gama, T, rel_perm,i))/this->getActivity(this->getReaction(i).Get_Species_Index());
	if (conc <= 0.0)
		conc = DBL_MIN;
	if (conc > this->getReaction(i).Get_MaximumValue())
		conc = this->getReaction(i).Get_MaximumValue() * 0.90;
	return conc;
}

//Return the ith area factor
double UnsteadyAdsorption::getAreaFactor(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->area_factors[i];
}

//Return the ith activity coefficient
double UnsteadyAdsorption::getActivity(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->activities(i,0);
}

//Return the ith old activity coefficient
double UnsteadyAdsorption::getOldActivity(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->activities_old(i,0);
}

//Return the specific molality
double UnsteadyAdsorption::getSpecificMolality()
{
	return this->specific_molality;
}

//Return the specific area
double UnsteadyAdsorption::getSpecificArea()
{
	return this->specific_area;
}

//Return the surface charge
double UnsteadyAdsorption::getSurfaceCharge()
{
	return this->surface_charge;
}

//Calculate and return the bulk density of the adsorbent in the system
double UnsteadyAdsorption::getBulkDensity()
{
	return this->total_mass / this->total_volume;
}

//Return the total mass
double UnsteadyAdsorption::getTotalMass()
{
	return this->total_mass;
}

//Return the total system volume
double UnsteadyAdsorption::getTotalVolume()
{
	return this->total_volume;
}

//Return charge density
double UnsteadyAdsorption::getChargeDensity()
{
	return this->charge_density;
}

//Return the ionic strength
double UnsteadyAdsorption::getIonicStrength()
{
	return this->ionic_strength;
}

//Return the number of reactions
int UnsteadyAdsorption::getNumberRxns()
{
	return this->num_rxns;
}

//Return index of the adsorbed species
int UnsteadyAdsorption::getAdsorbIndex(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	if (this->adsorb_index[i] >= this->List->list_size() || this->adsorb_index[i] < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	return this->adsorb_index[i];
}

//Return the index of the primary aqueous species
int UnsteadyAdsorption::getAqueousIndex(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	if (this->aqueous_index[i] >= this->List->list_size() || this->aqueous_index[i] < 0)
	{
		mError(out_of_bounds);
		return 0;
	}
	return this->aqueous_index[i];
}

//Return the activity enum
int UnsteadyAdsorption::getActivityEnum()
{
	return this->act_fun;
}

//Return true is in Area basis
bool UnsteadyAdsorption::isAreaBasis()
{
	return this->AreaBasis;
}

//Return state of surface charge inclusion
bool UnsteadyAdsorption::includeSurfaceCharge()
{
	return this->IncludeSurfCharge;
}

//Return the name of the adsorbent
std::string UnsteadyAdsorption::getAdsorbentName()
{
	return this->adsorbent_name;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: UnsteadyAdsorption
 */

/*
 *								Start: MultiligandAdsorption
 *	-------------------------------------------------------------------------------------
 */
//Default Constructor
MultiligandAdsorption::MultiligandAdsorption()
:
ligand_obj(0),
activities()
{
	num_ligands = 0;
	adsorbent_name = "Adsorbent A";

	surface_activity = (*ideal_solution);
	activity_data = nullptr;
	specific_area = 1.0;
	total_mass = 0.0;
	total_volume = 1.0;
	charge_density = 0.0;
	ionic_strength = 0.0;
	electric_potential = 0.0;
	IncludeSurfCharge = true;

}

//Default Destructor
MultiligandAdsorption::~MultiligandAdsorption()
{
	ligand_obj.clear();
}

//Initialize the list of adsorption objects
void MultiligandAdsorption::Initialize_Object(MasterSpeciesList &List, int ligand, std::vector<int> n)
{
	this->List = &List;
	this->num_ligands = ligand;
	this->ligand_obj.resize(ligand);
	this->activities.set_size(this->List->list_size(), 1);

	if (n.size() != ligand)
	{
		this->num_ligands = 0;
		mError(invalid_size);
		return;
	}

	for (int l=0; l<this->num_ligands; l++)
	{
		this->getAdsorptionObject(l).Initialize_Object(List, n[l]);
		this->getAdsorptionObject(l).setAreaBasisBool(false);
	}
}

//Modify the deltas for the mass balance object
void MultiligandAdsorption::modifyDeltas(MassBalance &mbo)
{
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getAdsorptionObject(l).setTotalMass(this->getTotalMass());
		this->getAdsorptionObject(l).setTotalVolume(this->getTotalVolume());
	}
	this->getAdsorptionObject(0).modifyDeltas(mbo);
}

//Set the adsorb indices for all reactions and ligands
int MultiligandAdsorption::setAdsorbIndices()
{
	int success = 0;

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		success = this->getAdsorptionObject(l).setAdsorbIndices();
		if (success == -1)
			return success;
	}

	return success;
}

//Check all aqueous indices from all reactions and all ligands for errors
int MultiligandAdsorption::checkAqueousIndices()
{
	int success = 0;

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		success = this->getAdsorptionObject(l).checkAqueousIndices();
		if (success == -1)
			return success;
	}

	return success;
}

//Set the activity function and data structure
void MultiligandAdsorption::setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data), const void *act_data)
{
	if ( (*act) == NULL )
	{
		this->surface_activity = (*ideal_solution);
	}
	else
	{
		this->surface_activity = (*act);
	}
	if ( (act_data) == NULL	)
	{
		this->activity_data = this;
	}
	else
	{
		this->activity_data = act_data;
	}

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getAdsorptionObject(l).setActivityModelInfo(NULL,NULL);
	}
}

//Set the aqueous indices for each adsorption object
int MultiligandAdsorption::setAqueousIndexAuto()
{
	int success = 0;

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		success = this->getAdsorptionObject(l).setAqueousIndexAuto();
		if (success == -1)
			return success;
	}

	return success;
}


//Set the enum for surface activity
void MultiligandAdsorption::setActivityEnum(int act)
{
	switch (act)
	{
		case IDEAL_ADS:
			this->act_fun = IDEAL_ADS;
			break;

		case FLORY_HUGGINS:
			this->act_fun = FLORY_HUGGINS;
			break;

		case UNIQUAC_ACT:
			this->act_fun = UNIQUAC_ACT;
			break;

		default:
			this->act_fun = IDEAL_ADS;
			break;
	}

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getAdsorptionObject(l).setActivityEnum(this->act_fun);
	}
}

//Set the molar factor
void MultiligandAdsorption::setMolarFactor(int ligand, int rxn, double m)
{
	if (ligand >= this->getNumberLigands() || ligand < 0)
	{
		mError(out_of_bounds);
		ligand = 0;
	}

	this->getAdsorptionObject(ligand).setMolarFactor(rxn, m);
}

//Set all volume factors
void MultiligandAdsorption::setVolumeFactor(int i, double v)
{
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getAdsorptionObject(l).setVolumeFactor(i, v);
	}
}

//Set all area factors
void MultiligandAdsorption::setAreaFactor(int i, double a)
{
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getAdsorptionObject(l).setAreaFactor(i, a);
	}
}

//Set the molality of the given ligand
void MultiligandAdsorption::setSpecificMolality(int ligand, double a)
{
	if (ligand >= this->getNumberLigands() || ligand < 0)
	{
		mError(out_of_bounds);
		ligand = 0;
	}

	this->getAdsorptionObject(ligand).setSpecificMolality(a);
}

//Set the surface charge of the given ligand
void MultiligandAdsorption::setSurfaceCharge(int ligand, double a)
{
	if (ligand >= this->getNumberLigands() || ligand < 0)
	{
		mError(out_of_bounds);
		ligand = 0;
	}

	this->getAdsorptionObject(ligand).setSurfaceCharge(a);
}

//Set the adsorbent name
void MultiligandAdsorption::setAdsorbentName(std::string name)
{
	this->adsorbent_name = name;
}

//Set the name of the ith ligand
void MultiligandAdsorption::setLigandName(int ligand, std::string name)
{
	this->getAdsorptionObject(ligand).setAdsorbentName(name);
}

//Set the specific area of the adsorbent
void MultiligandAdsorption::setSpecificArea(double area)
{
	if (area < 0)
		area = DBL_MIN;
	this->specific_area = area;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setSpecificArea(area);
}

//Set the total mass
void MultiligandAdsorption::setTotalMass(double mass)
{
	if (mass < 0)
		mass = DBL_MIN;
	this->total_mass = mass;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setTotalMass(mass);
}

//Set the total volume of the system
void MultiligandAdsorption::setTotalVolume(double volume)
{
	if (volume < 0)
		volume = DBL_MIN;
	this->total_volume = volume;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setTotalVolume(volume);
}

//Directly set the surface charging boolean
void MultiligandAdsorption::setSurfaceChargeBool(bool opt)
{
	this->IncludeSurfCharge = opt;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setSurfaceChargeBool(opt);
}

//Function to set the electic potential of the surface
void MultiligandAdsorption::setElectricPotential(double a)
{
	this->electric_potential = a;
}

//Function to calculate the area factors used by adsorption models
void MultiligandAdsorption::calculateAreaFactors()
{
	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).calculateAreaFactors();
}

//Function to calculate all equilibria parameters for all reactions
void MultiligandAdsorption::calculateEquilibria(double T)
{
	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).calculateEquilibria(T);
}

//Function to calculate and set the appropriate charge density
void MultiligandAdsorption::setChargeDensity(const Matrix<double> &x)
{

	this->charge_density = 0.0;
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->charge_density = this->charge_density + this->getAdsorptionObject(l).calculateSurfaceChargeDensity(x);
	}

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setChargeDensityValue(this->charge_density);

}

//Function to set the ionic strength parameter
void MultiligandAdsorption::setIonicStrength(const Matrix<double> &x)
{
	this->ionic_strength = calculate_ionic_strength(x, *this->List);

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setIonicStrengthValue(this->ionic_strength);
}

//Function to call the surface activity model
int MultiligandAdsorption::callSurfaceActivity(const Matrix<double> &x)
{
	int success = this->surface_activity(x,this->activities,this->activity_data);

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getAdsorptionObject(l).setActivities(this->activities);

	return success;
}

//Function to calculate the electric surface potential
void MultiligandAdsorption::calculateElecticPotential(double sigma, double T, double I, double rel_epsilon)
{
	I = I * 1000.0;//First, convert ionic strength to mol/m^3
	double coeff = sigma / sqrt(8.0*AbsPerm(rel_epsilon)*Rstd*T*I);
	double inv = log(coeff + sqrt((coeff*coeff)+1.0));
	this->electric_potential = (2.0 * kB * T * inv) / e;
}

//Calculation of equilibrium correct term
double MultiligandAdsorption::calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int rxn, int ligand)
{
	this->calculateElecticPotential(sigma, T, I, rel_epsilon);

	if (this->includeSurfaceCharge() == true)
		return -(this->getAdsorptionObject(ligand).calculateAqueousChargeExchange(rxn)*e*this->getElectricPotential())/(kB*T);
	else
		return 0.0;
}

//Calculation of residual
double MultiligandAdsorption::Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int rxn, int ligand)
{
	double res = this->getAdsorptionObject(ligand).getReaction(rxn).Get_Stoichiometric(this->getAdsorptionObject(ligand).getAdsorbIndex(rxn)) * ( log10(this->getActivity(this->getAdsorptionObject(ligand).getAdsorbIndex(rxn))) + x(this->getAdsorptionObject(ligand).getAdsorbIndex(rxn),0) );

	if (this->getAdsorptionObject(ligand).getReaction(rxn).Get_Stoichiometric(this->getAdsorptionObject(ligand).getAdsorbIndex(rxn)) > 0.0)
	{
		res = res - (this->getAdsorptionObject(ligand).getMolarFactor(rxn) * log10(this->getAdsorptionObject(ligand).getSpecificMolality()*this->getAdsorptionObject(ligand).calculateActiveFraction(x)));
	}
	else
	{
		res = res + (this->getAdsorptionObject(ligand).getMolarFactor(rxn) * log10(this->getAdsorptionObject(ligand).getSpecificMolality()*this->getAdsorptionObject(ligand).calculateActiveFraction(x)));
	}

	for (int n=0; n<this->List->list_size(); n++)
	{
		if (this->List->get_species(n).MoleculePhaseID() == AQUEOUS || this->List->get_species(n).MoleculePhaseID() == LIQUID)
			res = res + ( this->getAdsorptionObject(ligand).getReaction(rxn).Get_Stoichiometric(n)*( log10(gama(n,0))+x(n,0) ) );
	}

	double logK = this->getAdsorptionObject(ligand).getReaction(rxn).Get_Equilibrium();
	logK = logK + ((this->calculateEquilibriumCorrection(this->getChargeDensity(), T, this->getIonicStrength(), rel_perm, rxn, ligand))/log(10.0));

	res = res - logK;

	return res;
}

//Return reference to the corresponding adsorption object
AdsorptionReaction& MultiligandAdsorption::getAdsorptionObject(int l)
{
	if (l >= this->num_ligands || l < 0)
	{
		mError(out_of_bounds);
		l = 0;
	}
	return this->ligand_obj[l];
}

//Return number of ligands on the surface
int MultiligandAdsorption::getNumberLigands()
{
	return this->num_ligands;
}

//Return the activity enum
int MultiligandAdsorption::getActivityEnum()
{
	return this->act_fun;
}

//Return the value of activity for the ith species
double MultiligandAdsorption::getActivity(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->activities(i,0);
}

//Return the specific area of the adsorbent
double MultiligandAdsorption::getSpecificArea()
{
	return this->specific_area;
}

//Return the bulk density
double MultiligandAdsorption::getBulkDensity()
{
	return this->total_mass / this->total_volume;
}

//Return the total mass of adsorbent
double MultiligandAdsorption::getTotalMass()
{
	return this->total_mass;
}

//Return the total volume of the system
double MultiligandAdsorption::getTotalVolume()
{
	return this->total_volume;
}

//Return the charge density
double MultiligandAdsorption::getChargeDensity()
{
	return this->charge_density;
}

//Return ionic strength
double MultiligandAdsorption::getIonicStrength()
{
	return this->ionic_strength;
}

//Return the electric potential
double MultiligandAdsorption::getElectricPotential()
{
	return this->electric_potential;
}

//Return the value of the boolean for including surface charge
bool MultiligandAdsorption::includeSurfaceCharge()
{
	return this->IncludeSurfCharge;
}

//Return the name of the ligand at i
std::string MultiligandAdsorption::getLigandName(int l)
{
	return this->getAdsorptionObject(l).getAdsorbentName();
}

//Return the name of the adsorbent
std::string MultiligandAdsorption::getAdsorbentName()
{
	return this->adsorbent_name;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: MultiligandAdsorption
 */

/*
 *								Start: ChemisorptionReaction
 *	-------------------------------------------------------------------------------------
 */
//Default Constructor for Chemisorption Reaction
ChemisorptionReaction::ChemisorptionReaction()
{
	ads_rxn.resize(0);
	area_factors.resize(0);
	volume_factors.resize(0);
	adsorb_index.resize(0);
	surface_activity = (*ideal_solution);
	activity_data = nullptr;
	act_fun = IDEAL_ADS;
	specific_area = 1.0;
	specific_molality = 1.0;
	surface_charge = 0.0;
	total_mass = 0.0;
	total_volume = 1.0;
	charge_density = 0.0;
	ionic_strength = 0.0;
	num_rxns = 0;
	AreaBasis = false;
	IncludeSurfCharge = true;
	adsorbent_name = "HA";
	ligand_index = -1;
	Delta.resize(0);
}

//Default destructor for Chemisorption Reaction
ChemisorptionReaction::~ChemisorptionReaction()
{
	ads_rxn.clear();
	area_factors.clear();
	volume_factors.clear();
	adsorb_index.clear();
	Delta.clear();
}

//Initialize the object
void ChemisorptionReaction::Initialize_Object(MasterSpeciesList &List, int n)
{
	this->List = &List;
	this->num_rxns = n;
	this->area_factors.resize(this->List->list_size());
	this->volume_factors.resize(this->List->list_size());
	this->activities.set_size(this->List->list_size(),1);
	this->Delta.resize(this->List->list_size());
	for (int i=0; i<this->List->list_size(); i++)
	{
		this->area_factors[i] = 0.0;
		this->volume_factors[i] = 0.0;
		this->activities(i,0) = 1.0;
		this->Delta[i] = 0.0;
	}

	this->ads_rxn.resize(this->num_rxns);
	this->adsorb_index.resize(this->num_rxns);
	for (int i=0; i<this->num_rxns; i++)
	{
		this->ads_rxn[i].Initialize_Object(List);
		this->adsorb_index[i] = -1;
	}

}

//Display information about the object to the console
void ChemisorptionReaction::Display_Info()
{
	std::cout << "Ligand Name = " << this->getAdsorbentName() << std::endl;
	std::cout << "Site Balance (mol/kg) = " << this->getSpecificMolality() << " = ";
	//Loop through for site balance
	bool first = true;
	for (int i=0; i<this->List->list_size(); i++)
	{
		if (i == 0 || first == true)
		{
			if (this->getDelta(i) != 0.0)
			{
				if (this->getDelta(i) > 1.0)
					std::cout << this->getDelta(i) << " x [ " << this->List->get_species(i).MolecularFormula() << " ]";
				else
					std::cout << "[ " << this->List->get_species(i).MolecularFormula() << " ]";
				first = false;
			}
		}
		else
		{
			if (this->getDelta(i) != 0.0)
			{
				if (this->getDelta(i) > 1.0)
					std::cout << " + " << this->getDelta(i) << " x [ " << this->List->get_species(i).MolecularFormula() << " ]";
				else
					std::cout << " + [ " << this->List->get_species(i).MolecularFormula() << " ]";
			}
		}
	}

	std::cout << "\nReactions Involved...\n";
	for (int i=0; i<this->ads_rxn.size(); i++)
		this->getReaction(i).Display_Info();
}

//Modify the deltas of the system mass balance
void ChemisorptionReaction::modifyMBEdeltas(MassBalance &mbo)
{
	this->AdsorptionReaction::modifyDeltas(mbo);
}

//Set the adsorbed species index automatically
int ChemisorptionReaction::setAdsorbIndices()
{
	int success = 0;

	//Loop through all reactions to find solid/adsorbed species that are NOT the named ligand
	for (int i=0; i<this->getNumberRxns(); i++)
	{
		this->adsorb_index[i] = -1;
		for (int j=0; j<this->List->list_size(); j++)
		{
			//Check to see if jth species is involved with the ith reaction
			if (this->getReaction(i).Get_Stoichiometric(j) != 0.0)
			{
				//Check to see if species is solid or adsorbed
				if (this->List->get_species(j).MoleculePhaseID() == SOLID || this->List->get_species(j).MoleculePhaseID() == ADSORBED)
				{
					//Check to see if species is NOT the named ligand
					if (this->List->get_species(j).MolecularFormula() != this->getAdsorbentName())
					{
						this->adsorb_index[i] = j;
					}
				}
			}
		}
		if (this->adsorb_index[i] == -1)
		{
			mError(invalid_species);
			return -1;
		}
	}

	return success;
}

//Set the ligand species index automatically
int ChemisorptionReaction::setLigandIndex()
{
	int success = 0;

	//Loop through all reactions to find solid/adsorbed species that IS the named ligand
	for (int i=0; i<this->getNumberRxns(); i++)
	{
		for (int j=0; j<this->List->list_size(); j++)
		{
			//Check to see if jth species is involved with the ith reaction
			if (this->getReaction(i).Get_Stoichiometric(j) != 0.0)
			{
				//Check to see if species is solid or adsorbed
				if (this->List->get_species(j).MoleculePhaseID() == SOLID || this->List->get_species(j).MoleculePhaseID() == ADSORBED)
				{
					//Check to see if species IS the named ligand
					if (this->List->speciesName(j) == this->getAdsorbentName())
					{
						if (i == 0)
							this->ligand_index = j;
						else
						{
							if (j != this->ligand_index)
							{
								this->ligand_index = -1;
								mError(invalid_species);
								return -1;
							}
						}
					}
				}
			}
		}
		if (this->ligand_index == -1)
		{
			mError(invalid_species);
			return -1;
		}
	}
	this->setDelta(this->ligand_index, 1.0);

	return success;
}

//Set the deltas of the site balance automatically (given ligand index and stoicheometry)
int ChemisorptionReaction::setDeltas()
{
	int success = 0;

	for (int i=0; i<this->getNumberRxns(); i++)
	{
		this->setDelta( this->getAdsorbIndex(i), fabs(this->getReaction(i).Get_Stoichiometric(this->getLigandIndex())) );
	}

	return success;
}

//Set the activity model
void ChemisorptionReaction::setActivityModelInfo(int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
											  const void *act_data)
{
	this->AdsorptionReaction::setActivityModelInfo(act, act_data);
}

//Set the activity enum choice
void ChemisorptionReaction::setActivityEnum(int act)
{
	this->AdsorptionReaction::setActivityEnum(act);
}

//Set the delta value for the site balance
void ChemisorptionReaction::setDelta(int i, double v)
{
	//Check args
	if (i >= this->Delta.size() || i < 0)
	{
		mError(out_of_bounds);
		return;
	}

	this->Delta[i] = v;
}

//Set the volume factor of the given species
void ChemisorptionReaction::setVolumeFactor(int i, double v)
{
	this->AdsorptionReaction::setVolumeFactor(i, v);
}

//Set the area factor of the given species
void ChemisorptionReaction::setAreaFactor(int i, double a)
{
	this->AdsorptionReaction::setAreaFactor(i, a);
}

//Set the specific area of the adsorbent
void ChemisorptionReaction::setSpecificArea(double a)
{
	this->AdsorptionReaction::setSpecificArea(a);
}

//Set the specific molality of the adsorbent
void ChemisorptionReaction::setSpecificMolality(double a)
{
	this->AdsorptionReaction::setSpecificMolality(a);
}

//Set the total mass of the adsorbent
void ChemisorptionReaction::setTotalMass(double m)
{
	this->AdsorptionReaction::setTotalMass(m);
}

//Set the total volume of the system
void ChemisorptionReaction::setTotalVolume(double v)
{
	this->AdsorptionReaction::setTotalVolume(v);
}

//Set the surface charge boolean arg
void ChemisorptionReaction::setSurfaceChargeBool(bool opt)
{
	this->AdsorptionReaction::setSurfaceChargeBool(opt);
}

//Set the name of the ligand/adsorbent
void ChemisorptionReaction::setAdsorbentName(std::string name)
{
	this->AdsorptionReaction::setAdsorbentName(name);
}

//Set the value of charge density
void ChemisorptionReaction::setChargeDensityValue(double a)
{
	this->AdsorptionReaction::setChargeDensityValue(a);
}

//Set the value of ionic strength
void ChemisorptionReaction::setIonicStrengthValue(double a)
{
	this->AdsorptionReaction::setIonicStrengthValue(a);
}

//Set the activities to the given matrix
void ChemisorptionReaction::setActivities(Matrix<double> &x)
{
	this->AdsorptionReaction::setActivities(x);
}

//Calculate the area factors for adsorption
void ChemisorptionReaction::calculateAreaFactors()
{
	for (int r=0; r<this->getNumberRxns(); r++)
	{
		for (int i=0; i<this->List->list_size(); i++)
		{

			if (this->volume_factors[i] == 0.0 || (this->List->get_species(i).MoleculePhaseID() != SOLID && this->List->get_species(i).MoleculePhaseID() != ADSORBED))
			{
				this->area_factors[i] = 0.0;
			}
			else if (this->volume_factors[i] == 0.0 && (this->List->get_species(i).MoleculePhaseID() == SOLID || this->List->get_species(i).MoleculePhaseID() == ADSORBED))
			{
				this->area_factors[i] = 4.0 * M_PI * pow((3.0/(4.0*M_PI))*(4.33/Na), (2.0/3.0)) / 10000.0 * Na;
			}
			else
			{
				if (this->AreaBasis == false)
				{
					this->volume_factors[i] = this->volume_factors[i] * this->getReaction(r).Get_Stoichiometric(this->getLigandIndex());
				}
				this->area_factors[i] = 4.0 * M_PI * pow((3.0/(4.0*M_PI))*(this->volume_factors[i]/Na), (2.0/3.0)) / 10000.0 * Na;
			}
		}
	}

}

//Calculate the reaction equilibria constants based on temperature
void ChemisorptionReaction::calculateEquilibria(double T)
{
	for (int i=0; i<this->num_rxns; i++)
		this->getReaction(i).calculateEquilibrium(T);
}

//Set the value of surface charge density based on the calculation function
void ChemisorptionReaction::setChargeDensity(const Matrix<double> &x)
{
	this->charge_density = this->calculateSurfaceChargeDensity(x);
}

//Set the value of ionic strength based on the calculation function
void ChemisorptionReaction::setIonicStrength(const Matrix<double> &x)
{
	this->AdsorptionReaction::setIonicStrength(x);
}

//Function call to the surface activity function
int ChemisorptionReaction::callSurfaceActivity(const Matrix<double> &x)
{
	return this->AdsorptionReaction::callSurfaceActivity(x);
}

//Return the calculation of the surface charge density
double ChemisorptionReaction::calculateSurfaceChargeDensity(const Matrix<double> &x)
{
	double sigma = 0.0;
	double sum = 0.0;

	for (int i=0; i<this->getNumberRxns(); i++)
	{
		sum = sum + (this->List->charge(this->getAdsorbIndex(i)) * pow(10.0, x(this->getAdsorbIndex(i),0)));
	}
	sum = sum + (this->List->charge(this->getLigandIndex()) * pow(10.0, x(this->getLigandIndex(),0)));
	sigma = sum * (Faraday/this->getSpecificArea());

	return sigma;
}

//Return the calculation of the electric surface potential
double ChemisorptionReaction::calculateElecticPotential(double sigma, double T, double I, double rel_epsilon)
{
	return this->AdsorptionReaction::calculatePsi(sigma, T, I, rel_epsilon);
}

//Return the calculation of the net charge exchange term for the ith reaction
double ChemisorptionReaction::calculateAqueousChargeExchange(int i)
{
	double n = 0.0;

	for (int j = 0; j<this->List->list_size(); j++)
	{
		if (this->List->get_species(j).MoleculePhaseID() == AQUEOUS || this->List->get_species(j).MoleculePhaseID() == LIQUID)
		{
			n = n + (this->getReaction(i).Get_Stoichiometric(j)*this->List->charge(j));
		}
	}

	return -n;
}

//Return the calculation of the equilibrium correction term of the ith reaction
double ChemisorptionReaction::calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int i)
{
	if (this->includeSurfaceCharge() == true)
		return -(this->calculateAqueousChargeExchange(i)*e*calculateElecticPotential(sigma, T, I, rel_epsilon))/(kB*T);
	else
		return 0.0;
}

//Return the residual contributed by the ith reaction
double ChemisorptionReaction::Eval_RxnResidual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i)
{
	double res = 0.0;

	//Residual contributions from all aqueous species
	for (int j=0; j<this->List->list_size(); j++)
	{
		if (this->List->get_species(j).MoleculePhaseID() == AQUEOUS || this->List->get_species(j).MoleculePhaseID() == LIQUID || this->List->get_species(j).MoleculePhaseID() == GAS)
			res = res + ( this->getReaction(i).Get_Stoichiometric(j)*( log10(gama(j,0))+x(j,0) ) );
	}

	//Add residual contribution from the adsorbed species and ligand species
	res = res + this->getReaction(i).Get_Stoichiometric(this->getAdsorbIndex(i)) * ( log10(this->getActivity(this->getAdsorbIndex(i))) +  x(this->getAdsorbIndex(i),0));
	res = res + this->getReaction(i).Get_Stoichiometric(this->getLigandIndex()) * x(this->getLigandIndex(),0);

	//Subtract residual contribution from the binding constant
	double logK = this->getReaction(i).Get_Equilibrium();
	logK = logK + ((this->calculateEquilibriumCorrection(this->getChargeDensity(), T, this->getIonicStrength(), rel_perm, i))/log(10.0));
	res = res - logK;

	return res;
}

//Return the residual contributed the the site balance
double ChemisorptionReaction::Eval_SiteBalanceResidual(const Matrix<double> &x)
{
	double res = 0.0;

	//Loop for all species in list
	for (int i=0; i<this->List->list_size(); i++)
	{
		res = res + ( this->getDelta(i) * pow(10.0, x(i,0)) );
	}

	if (this->getSpecificMolality() <= DBL_MIN)
		res = (res - this->getSpecificMolality())/pow(DBL_EPSILON, 2.0);
	else
		res = (res / this->getSpecificMolality()) - 1.0;

	if (isnan(res) || isinf(res))
		res = sqrt(DBL_MAX)/this->List->list_size();

	return res;
}

//Return reference to the ith reaction object
Reaction& ChemisorptionReaction::getReaction(int i)
{
	if (i >= this->num_rxns || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->ads_rxn[i];
}

//Return the delta of the ith species
double ChemisorptionReaction::getDelta(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->Delta[i];
}

//Return the volume factor of the ith species
double ChemisorptionReaction::getVolumeFactor(int i)
{
	return this->AdsorptionReaction::getVolumeFactor(i);
}

//Return the area factor of the ith species
double ChemisorptionReaction::getAreaFactor(int i)
{
	return this->AdsorptionReaction::getAreaFactor(i);
}

//Return the activity coefficient of the ith species
double ChemisorptionReaction::getActivity(int i)
{
	return this->AdsorptionReaction::getActivity(i);
}

//Return the specific area of the adsorbent
double ChemisorptionReaction::getSpecificArea()
{
	return this->AdsorptionReaction::getSpecificArea();
}

//Return the specific area of the adsorbent
double ChemisorptionReaction::getSpecificMolality()
{
	return this->AdsorptionReaction::getSpecificMolality();
}

//Return the bulk density of the adsorbent in the system
double ChemisorptionReaction::getBulkDensity()
{
	return this->AdsorptionReaction::getBulkDensity();
}

//Return the total mass of the adsorbent
double ChemisorptionReaction::getTotalMass()
{
	return this->AdsorptionReaction::getTotalMass();
}

//Return the total volume of the system
double ChemisorptionReaction::getTotalVolume()
{
	return this->AdsorptionReaction::getTotalVolume();
}

//Return the charge density of the adsorbent
double ChemisorptionReaction::getChargeDensity()
{
	return this->AdsorptionReaction::getChargeDensity();
}

//Return the ionic strength of the system
double ChemisorptionReaction::getIonicStrength()
{
	return this->AdsorptionReaction::getIonicStrength();
}

//Return the number of reactions involved
int ChemisorptionReaction::getNumberRxns()
{
	return this->AdsorptionReaction::getNumberRxns();
}

//Return the adsorbed species index of the ith reaction
int ChemisorptionReaction::getAdsorbIndex(int i)
{
	return this->AdsorptionReaction::getAdsorbIndex(i);
}

//Return the ligand species index of the ith reaction
int ChemisorptionReaction::getLigandIndex()
{
	return this->ligand_index;
}

//Return the activity enum value
int ChemisorptionReaction::getActivityEnum()
{
	return this->AdsorptionReaction::getActivityEnum();
}

//Return the boolean for inclusion of surface charging
bool ChemisorptionReaction::includeSurfaceCharge()
{
	return this->AdsorptionReaction::includeSurfaceCharge();
}

//Return the name of the ligand/adsorbent
std::string ChemisorptionReaction::getAdsorbentName()
{
	return this->AdsorptionReaction::getAdsorbentName();
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: ChemisorptionReaction
 */

/*
 *								Start: MultiligandChemisorption
 *	-------------------------------------------------------------------------------------
 */
//Default constructor
MultiligandChemisorption::MultiligandChemisorption()
:
ligand_obj(0),
activities()
{
	num_ligands = 0;
	adsorbent_name = "Adsorbent A";

	surface_activity = (*ideal_solution);
	activity_data = nullptr;
	specific_area = 1.0;
	total_mass = 0.0;
	total_volume = 1.0;
	charge_density = 0.0;
	ionic_strength = 0.0;
	electric_potential = 0.0;
	IncludeSurfCharge = true;

}

//Default destructor
MultiligandChemisorption::~MultiligandChemisorption()
{
	ligand_obj.clear();
}

//Initialization of the object
void MultiligandChemisorption::Initialize_Object(MasterSpeciesList &List, int ligand, std::vector<int> n)
{
	this->List = &List;
	this->num_ligands = ligand;
	this->ligand_obj.resize(ligand);
	this->activities.set_size(this->List->list_size(), 1);

	if (n.size() != ligand)
	{
		this->num_ligands = 0;
		mError(invalid_size);
		return;
	}

	for (int l=0; l<this->num_ligands; l++)
	{
		this->getChemisorptionObject(l).Initialize_Object(List, n[l]);
	}

}

//Display info function
void MultiligandChemisorption::Display_Info()
{
	std::cout << "Adsorbent Name: " << this->adsorbent_name << std::endl;
	std::cout << "---------------------------------\n\n";
	for (int i=0; i<this->num_ligands; i++)
	{
		this->getChemisorptionObject(i).Display_Info();
	}
}

//Modify the deltas in the mass balance object given
void MultiligandChemisorption::modifyMBEdeltas(MassBalance &mbo)
{
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getChemisorptionObject(l).setTotalMass(this->getTotalMass());
		this->getChemisorptionObject(l).setTotalVolume(this->getTotalVolume());
	}
	this->getChemisorptionObject(0).modifyMBEdeltas(mbo);
}

//Automatically set the adsorbed species indices for all chemisorption objects and reactions
int MultiligandChemisorption::setAdsorbIndices()
{
	int success = 0;

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		success = this->getChemisorptionObject(l).setAdsorbIndices();
		if (success == -1)
			return success;
	}

	return success;
}

//Automatically set the ligand species indices for all chemisorption objects
int MultiligandChemisorption::setLigandIndices()
{
	int success = 0;

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		success = this->getChemisorptionObject(l).setLigandIndex();
		if (success == -1)
			return success;
	}

	return success;
}

//Automatically set the deltas for the site balances of all chemisorption objects
int MultiligandChemisorption::setDeltas()
{
	int success = 0;

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		success = this->getChemisorptionObject(l).setDeltas();
		if (success == -1)
			return success;
	}

	return success;
}

//Set the information necessary for the activity model
void MultiligandChemisorption::setActivityModelInfo(int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
													const void *act_data)
{
	if ( (*act) == NULL )
	{
		this->surface_activity = (*ideal_solution);
	}
	else
	{
		this->surface_activity = (*act);
	}
	if ( (act_data) == NULL	)
	{
		this->activity_data = this;
	}
	else
	{
		this->activity_data = act_data;
	}

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getChemisorptionObject(l).setActivityModelInfo(NULL,NULL);
	}
}

//Set the activity enum for the object
void MultiligandChemisorption::setActivityEnum(int act)
{
	switch (act)
	{
		case IDEAL_ADS:
			this->act_fun = IDEAL_ADS;
			break;

		case FLORY_HUGGINS:
			this->act_fun = FLORY_HUGGINS;
			break;

		case UNIQUAC_ACT:
			this->act_fun = UNIQUAC_ACT;
			break;

		default:
			this->act_fun = IDEAL_ADS;
			break;
	}

	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getChemisorptionObject(l).setActivityEnum(this->act_fun);
	}
}

//Set the volume factor for the ith species
void MultiligandChemisorption::setVolumeFactor(int i, double v)
{
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getChemisorptionObject(l).setVolumeFactor(i, v);
	}
}

//Set the area factor for the ith species
void MultiligandChemisorption::setAreaFactor(int i, double a)
{
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->getChemisorptionObject(l).setAreaFactor(i, a);
	}
}

//Set the specific molality of the given ligand
void MultiligandChemisorption::setSpecificMolality(int ligand, double a)
{
	if (ligand >= this->getNumberLigands() || ligand < 0)
	{
		mError(out_of_bounds);
		ligand = 0;
	}

	this->getChemisorptionObject(ligand).setSpecificMolality(a);
}

//Set the adsorbent name
void MultiligandChemisorption::setAdsorbentName(std::string name)
{
	this->adsorbent_name = name;
}

//Set the name of the given ligand
void MultiligandChemisorption::setLigandName(int ligand, std::string name)
{
	this->getChemisorptionObject(ligand).setAdsorbentName(name);
}

//Set the specific area of the adsorbent
void MultiligandChemisorption::setSpecificArea(double area)
{
	if (area < 0)
		area = DBL_MIN;
	this->specific_area = area;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setSpecificArea(area);
}

//Set the total mass of the adsorbent in the system
void MultiligandChemisorption::setTotalMass(double mass)
{
	if (mass < 0)
		mass = DBL_MIN;
	this->total_mass = mass;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setTotalMass(mass);
}

//Set the total volume of the system
void MultiligandChemisorption::setTotalVolume(double volume)
{
	if (volume < 0)
		volume = DBL_MIN;
	this->total_volume = volume;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setTotalVolume(volume);
}

//Set the boolean value for surface charging
void MultiligandChemisorption::setSurfaceChargeBool(bool opt)
{
	this->IncludeSurfCharge = opt;

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setSurfaceChargeBool(opt);
}

//Set the electric surface potential of the adsorbent to a given value
void MultiligandChemisorption::setElectricPotential(double a)
{
	this->electric_potential = a;
}

//Calculate the area factors needed by the chemisorption object
void MultiligandChemisorption::calculateAreaFactors()
{
	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).calculateAreaFactors();
}

//Calculate the equilibria constants for all reactions with all ligands given the temperature
void MultiligandChemisorption::calculateEquilibria(double T)
{
	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).calculateEquilibria(T);
}

//Calculation of the charge density of the adsorbent given the vector of non-linear system variables
void MultiligandChemisorption::setChargeDensity(const Matrix<double> &x)
{
	this->charge_density = 0.0;
	for (int l=0; l<this->getNumberLigands(); l++)
	{
		this->charge_density = this->charge_density + this->getChemisorptionObject(l).calculateSurfaceChargeDensity(x);
	}

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setChargeDensityValue(this->charge_density);
}

//Calculation of the ionic strength of the system given the vector of non-linear system variables
void MultiligandChemisorption::setIonicStrength(const Matrix<double> &x)
{
	this->ionic_strength = calculate_ionic_strength(x, *this->List);

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setIonicStrengthValue(this->ionic_strength);
}

//Call the surface activity function for the object
int MultiligandChemisorption::callSurfaceActivity(const Matrix<double> &x)
{
	int success = this->surface_activity(x,this->activities,this->activity_data);

	for (int l=0; l<this->getNumberLigands(); l++)
		this->getChemisorptionObject(l).setActivities(this->activities);

	return success;
}

//Calculation of the electric surface potential for the adsorbent
void MultiligandChemisorption::calculateElecticPotential(double sigma, double T, double I, double rel_epsilon)
{
	I = I * 1000.0;//First, convert ionic strength to mol/m^3
	double coeff = sigma / sqrt(8.0*AbsPerm(rel_epsilon)*Rstd*T*I);
	double inv = log(coeff + sqrt((coeff*coeff)+1.0));
	this->electric_potential = (2.0 * kB * T * inv) / e;
}

//Calculation of the equilibrium constant's correction factor
double MultiligandChemisorption::calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int rxn, int ligand)
{
	this->calculateElecticPotential(sigma, T, I, rel_epsilon);

	if (this->includeSurfaceCharge() == true)
		return -(this->getChemisorptionObject(ligand).calculateAqueousChargeExchange(rxn)*e*this->getElectricPotential())/(kB*T);
	else
		return 0.0;
}

//Evaluation of a reaction residual for the given rxn index and ligand index
double MultiligandChemisorption::Eval_RxnResidual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int rxn, int ligand)
{
	double res = 0.0;

	//Residual constributions from all aqueous species
	for (int j=0; j<this->List->list_size(); j++)
	{
		if (this->List->get_species(j).MoleculePhaseID() == AQUEOUS || this->List->get_species(j).MoleculePhaseID() == LIQUID || this->List->get_species(j).MoleculePhaseID() == GAS)
			res = res + ( this->getChemisorptionObject(ligand).getReaction(rxn).Get_Stoichiometric(j)*( log10(gama(j,0))+x(j,0) ) );
	}

	//Add residual contribution from the adsorbed species and ligand species
	res = res + this->getChemisorptionObject(ligand).getReaction(rxn).Get_Stoichiometric(this->getChemisorptionObject(ligand).getAdsorbIndex(rxn)) * ( log10(this->getActivity(this->getChemisorptionObject(ligand).getAdsorbIndex(rxn))) +  x(this->getChemisorptionObject(ligand).getAdsorbIndex(rxn),0));
	res = res + this->getChemisorptionObject(ligand).getReaction(rxn).Get_Stoichiometric(this->getChemisorptionObject(ligand).getLigandIndex()) * x(this->getChemisorptionObject(ligand).getLigandIndex(),0);

	//Subtract residual contribution from the binding constant
	double logK = this->getChemisorptionObject(ligand).getReaction(rxn).Get_Equilibrium();
	logK = logK + ((this->calculateEquilibriumCorrection(this->getChargeDensity(), T, this->getIonicStrength(), rel_perm, rxn, ligand))/log(10.0));
	res = res - logK;

	return res;
}

//Evaluation of the site balance residual for the given ligand
double MultiligandChemisorption::Eval_SiteBalanceResidual(const Matrix<double> &x, int ligand)
{
	return this->getChemisorptionObject(ligand).Eval_SiteBalanceResidual(x);
}

//Get the chemisorption object at the index
ChemisorptionReaction& MultiligandChemisorption::getChemisorptionObject(int ligand)
{
	if (ligand >= this->num_ligands || ligand < 0)
	{
		mError(out_of_bounds);
		ligand = 0;
	}
	return this->ligand_obj[ligand];
}

//Get the number of ligands for the object
int MultiligandChemisorption::getNumberLigands()
{
	return this->num_ligands;
}

//Get the activity enum for the object
int MultiligandChemisorption::getActivityEnum()
{
	return this->act_fun;
}

//Get the activity coefficient for the ith species
double MultiligandChemisorption::getActivity(int i)
{
	if (i >= this->List->list_size() || i < 0)
	{
		mError(out_of_bounds);
		i = 0;
	}
	return this->activities(i,0);
}

//Get the specific area of the adsorbent
double MultiligandChemisorption::getSpecificArea()
{
	return this->specific_area;
}

//Get the bulk density of the adsorbent in solution
double MultiligandChemisorption::getBulkDensity()
{
	return this->total_mass/this->total_volume;
}

//Get the total mass of adsorbent in the system
double MultiligandChemisorption::getTotalMass()
{
	return this->total_mass;
}

//Get the total volume of the system
double MultiligandChemisorption::getTotalVolume()
{
	return this->total_volume;
}

//Get the charge density for the adsorbent
double MultiligandChemisorption::getChargeDensity()
{
	return this->charge_density;
}

//Get the ionic strength of the solution
double MultiligandChemisorption::getIonicStrength()
{
	return this->ionic_strength;
}

//Get the electric potential of the adsorbent surface
double MultiligandChemisorption::getElectricPotential()
{
	return this->electric_potential;
}

//Return boolean to determine whether or not to include surface charge
bool MultiligandChemisorption::includeSurfaceCharge()
{
	return this->IncludeSurfCharge;
}

//Get the name of the ligand at the given index
std::string MultiligandChemisorption::getLigandName(int ligand)
{
	return this->getChemisorptionObject(ligand).getAdsorbentName();
}

//Get the name of the adsorbent
std::string MultiligandChemisorption::getAdsorbentName()
{
	return this->adsorbent_name;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: MultiligandChemisorption
 */

//Print SHARK info to output file
void print2file_shark_info(SHARK_DATA *shark_dat)
{
	if (shark_dat->File_Output == false)
		return;
	if (shark_dat->OutputFile == NULL)
	{
		mError(nullptr_error);
		return;
	}

	//Loop for all Reactions
	fprintf(shark_dat->OutputFile, "---------------SHARK SIMULATION: SOLUTION HEADER-----------------\n\n");
	fprintf(shark_dat->OutputFile, "Steady-State Reactions\n-----------------------------------------------------\n");
	for (int j=0; j<shark_dat->num_ssr; j++)
	{
		fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->ReactionList[j].Get_Equilibrium());
		bool first = true;
		for (int i=0; i<shark_dat->numvar; i++)
		{
			if (shark_dat->ReactionList[j].Get_Stoichiometric(i) < 0.0)
			{
				if (i == 0 || first == true)
				{
					if (fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					first = false;
				}
				else
				{
					if (fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
				}
			}
			if (i == shark_dat->numvar-1)
				fprintf(shark_dat->OutputFile, " = ");
		}

		first = true;
		for (int i=0; i<shark_dat->numvar; i++)
		{
			if (shark_dat->ReactionList[j].Get_Stoichiometric(i) > 0.0)
			{
				if (i == 0 || first == true)
				{
					if (fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					first = false;
				}
				else
				{
					if (fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->ReactionList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
				}
			}
			if (i == shark_dat->numvar-1)
				fprintf(shark_dat->OutputFile, "\n\n");
		}
	}

	//Loop for all Unsteady Reactions
	if (shark_dat->UnsteadyList.size() > 0)
		fprintf(shark_dat->OutputFile, "Unsteady Reactions\n-----------------------------------------------------\n");
	for (int j=0; j<shark_dat->num_usr; j++)
	{
		if (shark_dat->UnsteadyList[j].Get_Stoichiometric(shark_dat->UnsteadyList[j].Get_Species_Index()) != 1.0)
			fprintf(shark_dat->OutputFile, "(1 / %g) x ", fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(shark_dat->UnsteadyList[j].Get_Species_Index())));
		fprintf(shark_dat->OutputFile, "d { %s } / dt\nk_reverse =\t%.6g\t:\tk_forward =\t%.6g\n", shark_dat->MasterList.get_species(shark_dat->UnsteadyList[j].Get_Species_Index()).MolecularFormula().c_str(), shark_dat->UnsteadyList[j].Get_Reverse(), shark_dat->UnsteadyList[j].Get_Forward());
		fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->UnsteadyList[j].Get_Equilibrium());
		bool first = true;
		for (int i=0; i<shark_dat->numvar; i++)
		{
			if (shark_dat->UnsteadyList[j].Get_Stoichiometric(i) < 0.0)
			{
				if (i == 0 || first == true)
				{
					if (fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					first = false;
				}
				else
				{
					if (fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
				}
			}
			if (i == shark_dat->numvar-1)
				fprintf(shark_dat->OutputFile, "  <-- reverse : forward -->  ");
		}

		first = true;
		for (int i=0; i<shark_dat->numvar; i++)
		{
			if (shark_dat->UnsteadyList[j].Get_Stoichiometric(i) > 0.0)
			{
				if (i == 0 || first == true)
				{
					if (fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					first = false;
				}
				else
				{
					if (fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)) > 1.0)
						fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->UnsteadyList[j].Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
				}
			}
			if (i == shark_dat->numvar-1)
				fprintf(shark_dat->OutputFile, "\n\n");
		}
	}

	//Loop for all AdsorptionReaction Objects (steady-state)
	if (shark_dat->AdsorptionList.size() > 0)
	{
		fprintf(shark_dat->OutputFile, "Steady-State Adsorption Reactions\n-----------------------------------------------------\n");
		fprintf(shark_dat->OutputFile, "Number of Different Adsorbents = %i\n\n", (int) shark_dat->AdsorptionList.size());
		for (int n=0; n<shark_dat->AdsorptionList.size(); n++)
		{
			fprintf(shark_dat->OutputFile, "Reactions for Adsorbent %i \n-----------------------------------------------------\n", n+1);
			if (shark_dat->AdsorptionList[n].isAreaBasis() == true)
			{
				fprintf(shark_dat->OutputFile, "Adsorption Basis =\tArea\n");
				fprintf(shark_dat->OutputFile, "Specific Area (m^2/kg) = \t%.6g\n", shark_dat->AdsorptionList[n].getSpecificArea());
			}
			else
			{
				fprintf(shark_dat->OutputFile, "Adsorption Basis =\tMolar\n");
				fprintf(shark_dat->OutputFile, "Specific Molality (mol/kg) = \t%.6g\n", shark_dat->AdsorptionList[n].getSpecificMolality());
			}
			fprintf(shark_dat->OutputFile, "Total Mass (kg) = \t%.6g\n", shark_dat->AdsorptionList[n].getTotalMass());
			if (shark_dat->AdsorptionList[n].includeSurfaceCharge() == false)
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t FALSE\n");
			else
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t TRUE\n");
			int surf_act = shark_dat->AdsorptionList[n].getActivityEnum();
			switch (surf_act)
			{
				case IDEAL_ADS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;

				case FLORY_HUGGINS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t FLORY-HUGGINS\n");
					break;

				case UNIQUAC_ACT:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t UNIQUAC\n");
					break;

				default:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;
			}
			for (int j=0; j<shark_dat->AdsorptionList[n].getNumberRxns(); j++)
			{
				fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->AdsorptionList[n].getReaction(j).Get_Equilibrium());
				bool first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i) < 0.0)
					{
						if (i == 0 || first == true)
						{
							if (fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
						else
						{
							if (fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
					if (i == shark_dat->numvar-1)
					{
						if (shark_dat->AdsorptionList[n].isAreaBasis() == true || shark_dat->AdsorptionList[n].getMolarFactor(j) == 1.0)
							fprintf(shark_dat->OutputFile, " + { %s }_%i = ",shark_dat->AdsorptionList[n].getAdsorbentName().c_str(),n+1);
						else
							fprintf(shark_dat->OutputFile, " + %g x { %s }_%i = ",shark_dat->AdsorptionList[n].getMolarFactor(j),shark_dat->AdsorptionList[n].getAdsorbentName().c_str(),n+1);
					}
				}

				first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i) > 0.0)
					{
						if (i == 0 || first == true)
						{
							if (fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
						else
						{
							if (fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->AdsorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
					if (i == shark_dat->numvar-1)
						fprintf(shark_dat->OutputFile, "\n\n");
				}

			}//END RXN LOOP
		}
	}

	//Loop for all UnsteadyAdsorption Objects
	if (shark_dat->UnsteadyAdsList.size() > 0)
	{
		fprintf(shark_dat->OutputFile, "Unsteady Adsorption Reactions\n-----------------------------------------------------\n");
		fprintf(shark_dat->OutputFile, "Number of Different Adsorbents = %i\n\n", (int) shark_dat->UnsteadyAdsList.size());
		for (int n=0; n<shark_dat->UnsteadyAdsList.size(); n++)
		{
			fprintf(shark_dat->OutputFile, "Reactions for Adsorbent %i \n-----------------------------------------------------\n", n+1);
			if (shark_dat->UnsteadyAdsList[n].isAreaBasis() == true)
			{
				fprintf(shark_dat->OutputFile, "Adsorption Basis =\tArea\n");
				fprintf(shark_dat->OutputFile, "Specific Area (m^2/kg) = \t%.6g\n", shark_dat->UnsteadyAdsList[n].getSpecificArea());
			}
			else
			{
				fprintf(shark_dat->OutputFile, "Adsorption Basis =\tMolar\n");
				fprintf(shark_dat->OutputFile, "Specific Molality (mol/kg) = \t%.6g\n", shark_dat->UnsteadyAdsList[n].getSpecificMolality());
			}
			fprintf(shark_dat->OutputFile, "Total Mass (kg) = \t%.6g\n", shark_dat->UnsteadyAdsList[n].getTotalMass());
			if (shark_dat->UnsteadyAdsList[n].includeSurfaceCharge() == false)
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t FALSE\n");
			else
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t TRUE\n");
			int surf_act = shark_dat->UnsteadyAdsList[n].getActivityEnum();
			switch (surf_act)
			{
				case IDEAL_ADS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;

				case FLORY_HUGGINS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t FLORY-HUGGINS\n");
					break;

				case UNIQUAC_ACT:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t UNIQUAC\n");
					break;

				default:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;
			}
			for (int j=0; j<shark_dat->UnsteadyAdsList[n].getNumberRxns(); j++)
			{
				if (shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Species_Index()) != 1.0)
					fprintf(shark_dat->OutputFile, "(1 / %g) x ", fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Species_Index())));
				fprintf(shark_dat->OutputFile, "d { %s } / dt\nk_reverse =\t%.6g\t:\tk_forward =\t%.6g\n", shark_dat->MasterList.get_species(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Species_Index()).MolecularFormula().c_str(), shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Reverse(), shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Forward());

				fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Equilibrium());
				bool first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i) < 0.0)
					{
						if (i == 0 || first == true)
						{
							if (fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
						else
						{
							if (fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
					if (i == shark_dat->numvar-1)
					{
						if (shark_dat->UnsteadyAdsList[n].isAreaBasis() == true || shark_dat->UnsteadyAdsList[n].getMolarFactor(j) == 1.0)
							fprintf(shark_dat->OutputFile, " + { %s }_%i  <-- reverse : forward -->  ",shark_dat->UnsteadyAdsList[n].getAdsorbentName().c_str(),n+1);
						else
							fprintf(shark_dat->OutputFile, " + %g x { %s }_%i  <-- reverse : forward -->  ",shark_dat->UnsteadyAdsList[n].getMolarFactor(j),shark_dat->UnsteadyAdsList[n].getAdsorbentName().c_str(),n+1);
					}
				}

				first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i) > 0.0)
					{
						if (i == 0 || first == true)
						{
							if (fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
						else
						{
							if (fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->UnsteadyAdsList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
					if (i == shark_dat->numvar-1)
						fprintf(shark_dat->OutputFile, "\n\n");
				}

			}//END RXN LOOP
		}
	}//END Unsteady Ads Output


	//Print out all multisite adsorbent objects
	if (shark_dat->MultiAdsList.size() > 0)
	{
		fprintf(shark_dat->OutputFile, "Steady-State Multi-ligand Adsorbent Objects\n-----------------------------------------------------\n");
		fprintf(shark_dat->OutputFile, "Number of Different Adsorbents = %i\n\n", (int) shark_dat->MultiAdsList.size());

		for (int k=0; k<shark_dat->MultiAdsList.size(); k++)
		{
			fprintf(shark_dat->OutputFile, "Adsorbent:\t %s \n-----------------------------------------------------\n", shark_dat->MultiAdsList[k].getAdsorbentName().c_str());
			fprintf(shark_dat->OutputFile, "Specific Area (m^2/kg) = \t%.6g\n", shark_dat->MultiAdsList[k].getSpecificArea());
			fprintf(shark_dat->OutputFile, "Total Mass (kg) = \t%.6g\n", shark_dat->MultiAdsList[k].getTotalMass());
			if (shark_dat->MultiAdsList[k].includeSurfaceCharge() == false)
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t FALSE\n");
			else
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t TRUE\n");
			int surf_act = shark_dat->MultiAdsList[k].getActivityEnum();
			switch (surf_act)
			{
				case IDEAL_ADS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;

				case FLORY_HUGGINS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t FLORY-HUGGINS\n");
					break;

				case UNIQUAC_ACT:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t UNIQUAC\n");
					break;

				default:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;
			}

			fprintf(shark_dat->OutputFile, "\nList of Ligands and Reactions \n-----------------------------------------------------\n");

			//Loop through ligands
			for (int n=0; n<shark_dat->MultiAdsList[k].getNumberLigands(); n++)
			{
				fprintf(shark_dat->OutputFile, "Ligand:\t %s \n", shark_dat->MultiAdsList[k].getLigandName(n).c_str());
				fprintf(shark_dat->OutputFile, "Specific Molality (mol/kg) = \t%.6g\n", shark_dat->MultiAdsList[k].getAdsorptionObject(n).getSpecificMolality());

				//Loop through all reactions
				for (int j=0; j<shark_dat->MultiAdsList[k].getAdsorptionObject(n).getNumberRxns(); j++)
				{
					fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Equilibrium());
					bool first = true;
					for (int i=0; i<shark_dat->numvar; i++)
					{
						if (shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i) < 0.0)
						{
							if (i == 0 || first == true)
							{
								if (fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								first = false;
							}
							else
							{
								if (fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							}
						}
						if (i == shark_dat->numvar-1)
						{
							if (shark_dat->MultiAdsList[k].getAdsorptionObject(n).isAreaBasis() == true || shark_dat->MultiAdsList[k].getAdsorptionObject(n).getMolarFactor(j) == 1.0)
								fprintf(shark_dat->OutputFile, " + { %s }_%i = ",shark_dat->MultiAdsList[k].getAdsorptionObject(n).getAdsorbentName().c_str(),n+1);
							else
								fprintf(shark_dat->OutputFile, " + %g x { %s }_%i = ",shark_dat->MultiAdsList[k].getAdsorptionObject(n).getMolarFactor(j),shark_dat->MultiAdsList[k].getAdsorptionObject(n).getAdsorbentName().c_str(),n+1);
						}
					}

					first = true;
					for (int i=0; i<shark_dat->numvar; i++)
					{
						if (shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i) > 0.0)
						{
							if (i == 0 || first == true)
							{
								if (fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								first = false;
							}
							else
							{
								if (fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->MultiAdsList[k].getAdsorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							}
						}
						if (i == shark_dat->numvar-1)
							fprintf(shark_dat->OutputFile, "\n\n");
					}
				}

			}//END ligand Loop

		}//END Obj Loop
	}

	//Loop for all ChemisorptionReaction Objects (steady-state)
	if (shark_dat->ChemisorptionList.size() > 0)
	{
		fprintf(shark_dat->OutputFile, "Steady-State Chemisorption Reactions\n-----------------------------------------------------\n");
		fprintf(shark_dat->OutputFile, "Number of Different Adsorbents = %i\n\n", (int) shark_dat->ChemisorptionList.size());
		for (int n=0; n<shark_dat->ChemisorptionList.size(); n++)
		{
			fprintf(shark_dat->OutputFile, "Site Balance for Adsorbent %i \n-----------------------------------------------------\n", n+1);
			fprintf(shark_dat->OutputFile, "Specific Molality (mol/kg) = \t%.6g\t = \t", shark_dat->ChemisorptionList[n].getSpecificMolality());
			bool first = true;
			for (int i=0; i<shark_dat->numvar; i++)
			{
				if (i==0 || first==true)
				{
					if (shark_dat->ChemisorptionList[n].getDelta(i) != 0.0)
					{
						if (shark_dat->ChemisorptionList[n].getDelta(i) > 1.0)
							fprintf(shark_dat->OutputFile,"%.6g x [%s]", shark_dat->ChemisorptionList[n].getDelta(i),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						else
							fprintf(shark_dat->OutputFile,"[%s]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						first = false;
					}
				}
				else
				{
					if (shark_dat->ChemisorptionList[n].getDelta(i) != 0.0)
					{
						if (shark_dat->ChemisorptionList[n].getDelta(i) > 1.0)
							fprintf(shark_dat->OutputFile," + %.6g x [%s]", shark_dat->ChemisorptionList[n].getDelta(i),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						else
							fprintf(shark_dat->OutputFile," + [%s]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					}
				}
			}
			fprintf(shark_dat->OutputFile,"\n\n");

			fprintf(shark_dat->OutputFile, "Reactions for Adsorbent %i \n-----------------------------------------------------\n", n+1);
			fprintf(shark_dat->OutputFile, "Total Mass (kg) = \t%.6g\n", shark_dat->ChemisorptionList[n].getTotalMass());

			if (shark_dat->ChemisorptionList[n].includeSurfaceCharge() == false)
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t FALSE\n");
			else
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t TRUE\n");
			int surf_act = shark_dat->ChemisorptionList[n].getActivityEnum();
			switch (surf_act)
			{
				case IDEAL_ADS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;

				case FLORY_HUGGINS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t FLORY-HUGGINS\n");
					break;

				case UNIQUAC_ACT:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t UNIQUAC\n");
					break;

				default:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;
			}

			for (int j=0; j<shark_dat->ChemisorptionList[n].getNumberRxns(); j++)
			{
				fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->ChemisorptionList[n].getReaction(j).Get_Equilibrium());
				first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i) < 0.0)
					{
						if (i == 0 || first == true)
						{
							if (fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
						else
						{
							if (fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
					if (i == shark_dat->numvar-1)
						fprintf(shark_dat->OutputFile, " = ");
				}

				first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i) > 0.0)
					{
						if (i == 0 || first == true)
						{
							if (fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
						else
						{
							if (fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)) > 1.0)
								fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->ChemisorptionList[n].getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
					if (i == shark_dat->numvar-1)
						fprintf(shark_dat->OutputFile, "\n\n");
				}

			}//END RXN LOOP
		}

	}//End if Chemisorption

	//Print out all multisite chemisorption objects
	if (shark_dat->MultiChemList.size() > 0)
	{
		fprintf(shark_dat->OutputFile, "Steady-State Multi-ligand Chemisorption Objects\n-----------------------------------------------------\n");
		fprintf(shark_dat->OutputFile, "Number of Different Adsorbents = %i\n\n", (int) shark_dat->MultiChemList.size());

		for (int k=0; k<shark_dat->MultiChemList.size(); k++)
		{
			fprintf(shark_dat->OutputFile, "Adsorbent:\t %s \n-----------------------------------------------------\n", shark_dat->MultiChemList[k].getAdsorbentName().c_str());
			fprintf(shark_dat->OutputFile, "Specific Area (m^2/kg) = \t%.6g\n", shark_dat->MultiChemList[k].getSpecificArea());
			fprintf(shark_dat->OutputFile, "Total Mass (kg) = \t%.6g\n", shark_dat->MultiChemList[k].getTotalMass());
			if (shark_dat->MultiChemList[k].includeSurfaceCharge() == false)
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t FALSE\n");
			else
				fprintf(shark_dat->OutputFile, "Include Surface Charging =\t TRUE\n");
			int surf_act = shark_dat->MultiChemList[k].getActivityEnum();
			switch (surf_act)
			{
				case IDEAL_ADS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;

				case FLORY_HUGGINS:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t FLORY-HUGGINS\n");
					break;

				case UNIQUAC_ACT:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t UNIQUAC\n");
					break;

				default:
					fprintf(shark_dat->OutputFile, "Surface Activity =\t IDEAL\n");
					break;
			}

			fprintf(shark_dat->OutputFile, "\nList of Site Balances and Reactions \n-----------------------------------------------------\n");

			//Loop through ligands
			for (int n=0; n<shark_dat->MultiChemList[k].getNumberLigands(); n++)
			{
				fprintf(shark_dat->OutputFile, "Site Balance for Ligand:\t %s \n-----------------------------------------------------\n", shark_dat->MultiChemList[k].getLigandName(n).c_str());
				fprintf(shark_dat->OutputFile, "Specific Molality (mol/kg) = \t%.6g\t = \t", shark_dat->MultiChemList[k].getChemisorptionObject(n).getSpecificMolality());
				bool first = true;
				for (int i=0; i<shark_dat->numvar; i++)
				{
					if (i==0 || first==true)
					{
						if (shark_dat->MultiChemList[k].getChemisorptionObject(n).getDelta(i) != 0.0)
						{
							if (shark_dat->MultiChemList[k].getChemisorptionObject(n).getDelta(i) > 1.0)
								fprintf(shark_dat->OutputFile,"%.6g x [%s]", shark_dat->MultiChemList[k].getChemisorptionObject(n).getDelta(i),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile,"[%s]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							first = false;
						}
					}
					else
					{
						if (shark_dat->MultiChemList[k].getChemisorptionObject(n).getDelta(i) != 0.0)
						{
							if (shark_dat->MultiChemList[k].getChemisorptionObject(n).getDelta(i) > 1.0)
								fprintf(shark_dat->OutputFile," + %.6g x [%s]", shark_dat->MultiChemList[k].getChemisorptionObject(n).getDelta(i),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							else
								fprintf(shark_dat->OutputFile," + [%s]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
						}
					}
				}
				fprintf(shark_dat->OutputFile,"\n\n");
				fprintf(shark_dat->OutputFile, "Reactions for Ligand:\t %s \n-----------------------------------------------------\n", shark_dat->MultiChemList[k].getLigandName(n).c_str());

				//Loop through all reactions
				for (int j=0; j<shark_dat->MultiChemList[k].getChemisorptionObject(n).getNumberRxns(); j++)
				{
					fprintf(shark_dat->OutputFile, "logK = \t%.6g\t:\t",shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Equilibrium());
					first = true;
					for (int i=0; i<shark_dat->numvar; i++)
					{
						if (shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i) < 0.0)
						{
							if (i == 0 || first == true)
							{
								if (fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								first = false;
							}
							else
							{
								if (fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							}
						}
						if (i == shark_dat->numvar-1)
							fprintf(shark_dat->OutputFile, " = ");
					}

					first = true;
					for (int i=0; i<shark_dat->numvar; i++)
					{
						if (shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i) > 0.0)
						{
							if (i == 0 || first == true)
							{
								if (fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, "%g x { %s }",fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, "{ %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								first = false;
							}
							else
							{
								if (fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)) > 1.0)
									fprintf(shark_dat->OutputFile, " + %g x { %s }",fabs(shark_dat->MultiChemList[k].getChemisorptionObject(n).getReaction(j).Get_Stoichiometric(i)),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
								else
									fprintf(shark_dat->OutputFile, " + { %s }",shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
							}
						}
						if (i == shark_dat->numvar-1)
							fprintf(shark_dat->OutputFile, "\n\n");
					}

				}

			}//END ligand Loop

		}//END Obj Loop

	}//End if MultiChemi


	//Loop for all Mass Balances
	fprintf(shark_dat->OutputFile, "Mass Balances\n-----------------------------------------------------\n");
	for (int j=0; j<shark_dat->num_mbe; j++)
	{
		fprintf(shark_dat->OutputFile, "%s\n",shark_dat->MassBalanceList[j].Get_Name().c_str());
		if (shark_dat->reactor_type == BATCH)
			fprintf(shark_dat->OutputFile, "%.6g\t= ", shark_dat->MassBalanceList[j].Get_TotalConcentration());
		else
			fprintf(shark_dat->OutputFile, "%.6g\t= ", shark_dat->MassBalanceList[j].Get_InletConcentration());
		bool first = true;
		for (int i=0; i<shark_dat->numvar; i++)
		{
			if (i == 0 || first == true)
			{
				if (shark_dat->MassBalanceList[j].Get_Delta(i) != 0.0)
				{
					if (shark_dat->MassBalanceList[j].Get_Delta(i) != 1.0)
						fprintf(shark_dat->OutputFile, "%.6g x [ %s ]",shark_dat->MassBalanceList[j].Get_Delta(i),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, "[ %s ]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					first = false;
				}
			}
			else
			{
				if (shark_dat->MassBalanceList[j].Get_Delta(i) != 0.0)
				{
					if (shark_dat->MassBalanceList[j].Get_Delta(i) != 1.0)
						fprintf(shark_dat->OutputFile, " + %.6g x [ %s ]",shark_dat->MassBalanceList[j].Get_Delta(i),shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
					else
						fprintf(shark_dat->OutputFile, " + [ %s ]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
				}
			}
		}
		fprintf(shark_dat->OutputFile, "\n\n");
	}
	fprintf(shark_dat->OutputFile, "---------------END OF SOLUTION HEADER-----------------\n\n");
}

//Header for the result output from shark simulations
void print2file_shark_header(SHARK_DATA *shark_dat)
{
	if (shark_dat->File_Output == false)
		return;
	if (shark_dat->OutputFile == NULL)
	{
		mError(nullptr_error);
		return;
	}

	fprintf(shark_dat->OutputFile, "-----------------SHARK SIMULATION CONDITIONS----------------\n\n");
	switch (shark_dat->reactor_type)
	{
		case BATCH:
			fprintf(shark_dat->OutputFile, "Reactor Type = BATCH\n");
			break;

		case CSTR:
			fprintf(shark_dat->OutputFile, "Reactor Type = CSTR\n");
			break;

		case PFR:
			fprintf(shark_dat->OutputFile, "Reactor Type = PFR\n");
			break;

		default:
			fprintf(shark_dat->OutputFile, "Reactor Type = BATCH\n");
			break;
	}
	if (shark_dat->reactor_type == CSTR)
	{
		fprintf(shark_dat->OutputFile, "Flow Rate (L/hr) = \t%.6g\n", shark_dat->flow_rate);
	}
	if (shark_dat->reactor_type == PFR)
	{
		fprintf(shark_dat->OutputFile, "Flow Rate (L/hr) = \t%.6g\n", shark_dat->flow_rate);
		fprintf(shark_dat->OutputFile, "Cross Section (m^2) = \t%.6g\n", shark_dat->xsec_area);
	}
	fprintf(shark_dat->OutputFile, "Total Volume (L) = \t%.6g\n", shark_dat->volume);
	if (shark_dat->steadystate == true)
		fprintf(shark_dat->OutputFile, "Steady-State = TRUE\n");
	else
	{
		fprintf(shark_dat->OutputFile, "Steady-State = FALSE\n");
		fprintf(shark_dat->OutputFile, "\tSimulation Time = %.6g\n",shark_dat->simulationtime);
	}
	if (shark_dat->act_fun == IDEAL)
		fprintf(shark_dat->OutputFile, "Activity Function = IDEAL\n");
	else if (shark_dat->act_fun == DAVIES)
		fprintf(shark_dat->OutputFile, "Activity Function = DAVIES\n");
	else if (shark_dat->act_fun == DEBYE_HUCKEL)
		fprintf(shark_dat->OutputFile, "Activity Function = DEBYE_HUCKEL\n");
	else if (shark_dat->act_fun == SIT)
		fprintf(shark_dat->OutputFile, "Activity Function = SIT\n");
	else if (shark_dat->act_fun == PITZER)
		fprintf(shark_dat->OutputFile, "Activity Function = PITZER\n");
	else
		fprintf(shark_dat->OutputFile, "Activity Function = UNREGISTERED\n");
	if (shark_dat->Contains_pH == true)
	{
		if (shark_dat->const_pH == true)
		{
			fprintf(shark_dat->OutputFile, "Constant pH = TRUE\n");
			if (shark_dat->SpeciationCurve == false)
				fprintf(shark_dat->OutputFile, "\tpH = %.6g\n",shark_dat->pH);
			else
				fprintf(shark_dat->OutputFile, "\tpH to be varied between 1 and 14 to produce speciation curves\n");
		}
		else
			fprintf(shark_dat->OutputFile, "Constant pH = FALSE\n");
	}
	else
		fprintf(shark_dat->OutputFile, "pH Not Registered!\n");
	fprintf(shark_dat->OutputFile, "Residual Alkalinity = %.6g (M)\n",shark_dat->MasterList.alkalinity());
	fprintf(shark_dat->OutputFile, "\n-----------------END OF SIMULATION CONDITIONS----------------\n");

	fprintf(shark_dat->OutputFile, "\n----------------PJFNK SOLVER OPTIONS----------------\n\n");
	if (shark_dat->Newton_data.Bounce == true)
	{
		fprintf(shark_dat->OutputFile, "Line Search = TRUE\n");
		fprintf(shark_dat->OutputFile, "\tSearch Type = Bouncing Backtrack\n");
	}
	else if (shark_dat->Newton_data.Bounce == false && shark_dat->Newton_data.LineSearch == true)
	{
		fprintf(shark_dat->OutputFile, "Line Search = TRUE\n");
		fprintf(shark_dat->OutputFile, "\tSearch Type = Standard Backtracking\n");
	}
	else
		fprintf(shark_dat->OutputFile, "Line Search = FALSE\n");
	if (shark_dat->lin_precon == NULL)
		fprintf(shark_dat->OutputFile, "Preconditioning = FALSE\n");
	else
		fprintf(shark_dat->OutputFile, "Preconditioning = USER DEFINED\n");
	fprintf(shark_dat->OutputFile, "Linear Solver = ");
	if (shark_dat->Newton_data.linear_solver == PCG)
		fprintf(shark_dat->OutputFile, "PCG\n");
	else if (shark_dat->Newton_data.linear_solver == CGS)
		fprintf(shark_dat->OutputFile, "CGS\n");
	else if (shark_dat->Newton_data.linear_solver == BiCGSTAB)
		fprintf(shark_dat->OutputFile, "BiCGSTAB\n");
	else if (shark_dat->Newton_data.linear_solver == FOM)
		fprintf(shark_dat->OutputFile, "FOM\n");
	else if (shark_dat->Newton_data.linear_solver == GMRESLP)
		fprintf(shark_dat->OutputFile, "GMRESLP\n");
	else if (shark_dat->Newton_data.linear_solver == GMRESRP)
		fprintf(shark_dat->OutputFile, "GMRESRP\n");
	else if (shark_dat->Newton_data.linear_solver == GCR)
		fprintf(shark_dat->OutputFile, "GCR\n");
	else if (shark_dat->Newton_data.linear_solver == GMRESR)
		fprintf(shark_dat->OutputFile, "GMRESR\n");
	else if (shark_dat->Newton_data.linear_solver == KMS)
		fprintf(shark_dat->OutputFile, "KMS\n");
	else if (shark_dat->Newton_data.linear_solver == QR)
		fprintf(shark_dat->OutputFile, "QR\n");
	else
		fprintf(shark_dat->OutputFile, "UNDEFINED!\n");
	if (shark_dat->Newton_data.linear_solver == GMRESR)
	{
		fprintf(shark_dat->OutputFile, "\tgmres_tol = %.6g\n",shark_dat->Newton_data.gmresr_dat.gmres_tol);
		fprintf(shark_dat->OutputFile, "\tgmres_restart = %i\n",shark_dat->Newton_data.gmresr_dat.gmres_restart);
	}
	if (shark_dat->Newton_data.linear_solver == KMS)
	{
		fprintf(shark_dat->OutputFile, "\tinner_tol = %.6g\n",shark_dat->Newton_data.kms_dat.inner_reltol);
		fprintf(shark_dat->OutputFile, "\tmax_level = %i\n",shark_dat->Newton_data.kms_dat.max_level);
	}
	fprintf(shark_dat->OutputFile, "Maximum Non-Linear Iterations = %i\n", shark_dat->Newton_data.nl_maxit);
	fprintf(shark_dat->OutputFile, "Absolute Non-Linear Tolerance = %.6g\n", shark_dat->Newton_data.nl_tol_abs);
	fprintf(shark_dat->OutputFile, "Relative Non-Linear Tolerance = %.6g\n", shark_dat->Newton_data.nl_tol_rel);
	fprintf(shark_dat->OutputFile, "Absolute Linear Tolerance = %.6g\n", shark_dat->Newton_data.lin_tol_abs);
	fprintf(shark_dat->OutputFile, "Relative Linear Tolerance = %.6g\n", shark_dat->Newton_data.lin_tol_rel);
	fprintf(shark_dat->OutputFile, "\n-----------------END SOLVER OPTIONS-------------------\n\n");

	fprintf(shark_dat->OutputFile, "-----------------SHARK SIMULATION RESULTS-------------------\n\n");

	if (shark_dat->steadystate == false)
		fprintf(shark_dat->OutputFile, "(hr)");
	else
		fprintf(shark_dat->OutputFile, "(M)");
	fprintf(shark_dat->OutputFile, "\t(K)\t(-)");
	for (int i=0; i<shark_dat->MasterList.list_size(); i++)
	{
		switch (shark_dat->MasterList.get_species(i).MoleculePhaseID())
		{
			case AQUEOUS:
				fprintf(shark_dat->OutputFile, "\t(mol/L)");
				break;

			case LIQUID:
				fprintf(shark_dat->OutputFile, "\t(-)");
				break;

			case GAS:
				fprintf(shark_dat->OutputFile, "\t(kPa)");
				break;

			case SOLID:
				fprintf(shark_dat->OutputFile, "\t(mol/kg)");
				break;

			case ADSORBED:
				fprintf(shark_dat->OutputFile, "\t(mol/kg)");
				break;

			default:
				fprintf(shark_dat->OutputFile, "\t(-)");
				break;
		}
	}
	fprintf(shark_dat->OutputFile, "\t(-)\t(-)\t(-)\t(-)");
	fprintf(shark_dat->OutputFile, "\n");

	if (shark_dat->steadystate == false)
		fprintf(shark_dat->OutputFile, "Time");
	else
		fprintf(shark_dat->OutputFile, "Ionic-Strength");
	fprintf(shark_dat->OutputFile, "\tT\tpH");
	for (int i=0; i<shark_dat->MasterList.list_size(); i++)
		fprintf(shark_dat->OutputFile, "\t[ %s ]", shark_dat->MasterList.get_species(i).MolecularFormula().c_str());
	fprintf(shark_dat->OutputFile, "\tConverged?\tE.Norm\tNL_iter\tL_iter");
	fprintf(shark_dat->OutputFile, "\n");
}

//Printout new results
void print2file_shark_results_new(SHARK_DATA *shark_dat)
{
	if (shark_dat->File_Output == false)
		return;
	if (shark_dat->OutputFile == NULL)
	{
		mError(nullptr_error);
		return;
	}

	if (shark_dat->time < 0.0)
		fprintf(shark_dat->OutputFile, "inf\t%.6g\t%.6g", shark_dat->temperature, shark_dat->pH);
	else
	{
		if (shark_dat->steadystate == false)
			fprintf(shark_dat->OutputFile, "%.6g\t%.6g\t%.6g", shark_dat->time, shark_dat->temperature, shark_dat->pH);
		else
			fprintf(shark_dat->OutputFile, "%.6g\t%.6g\t%.6g", shark_dat->ionic_strength, shark_dat->temperature, shark_dat->pH);
	}
	for (int i=0; i<shark_dat->MasterList.list_size(); i++)
		fprintf(shark_dat->OutputFile, "\t%.6g",shark_dat->Conc_new(i,0));
	if (shark_dat->Converged == true)
		fprintf(shark_dat->OutputFile, "\tTRUE");
	else
	{
		if (shark_dat->LocalMin == true)
			fprintf(shark_dat->OutputFile, "\tLOCAL_MIN");
		else
			fprintf(shark_dat->OutputFile, "\tFALSE");
	}
	fprintf(shark_dat->OutputFile, "\t%.6g\t%i\t%i\n",shark_dat->Norm,shark_dat->Newton_data.nl_iter, shark_dat->Newton_data.l_iter);
}

//Printout old results
void print2file_shark_results_old(SHARK_DATA *shark_dat)
{
	if (shark_dat->File_Output == false)
		return;
	if (shark_dat->OutputFile == NULL)
	{
		mError(nullptr_error);
		return;
	}

	if (shark_dat->time < 0.0)
		fprintf(shark_dat->OutputFile, "inf\t%.6g\t%.6g", shark_dat->temperature, shark_dat->pH);
	else
	{
		if (shark_dat->steadystate == false)
			fprintf(shark_dat->OutputFile, "%.6g\t%.6g\t%.6g", shark_dat->time, shark_dat->temperature, shark_dat->pH);
		else
			fprintf(shark_dat->OutputFile, "%.6g\t%.6g\t%.6g", shark_dat->ionic_strength, shark_dat->temperature, shark_dat->pH);
	}
	for (int i=0; i<shark_dat->MasterList.list_size(); i++)
		fprintf(shark_dat->OutputFile, "\t%.6g",shark_dat->Conc_old(i,0));
	if (shark_dat->Converged == true)
		fprintf(shark_dat->OutputFile, "\tTRUE");
	else
		fprintf(shark_dat->OutputFile, "\tFALSE");
	fprintf(shark_dat->OutputFile, "\t%.6g\t%i\t%i\n",shark_dat->Norm,shark_dat->Newton_data.nl_iter, shark_dat->Newton_data.l_iter);
}

//Calculation of ionic strength
double calculate_ionic_strength(const Matrix<double> &x, MasterSpeciesList &MasterList)
{
	double I = 0.0;

	//Loop to calculate ionic strength
	for (int i=0; i<MasterList.list_size(); i++)
	{
		if (MasterList.get_species(i).MoleculePhaseID() == AQUEOUS || MasterList.get_species(i).MoleculePhaseID() == LIQUID)
			I = I + (pow(10.0,x(i,0)) * MasterList.charge(i) * MasterList.charge(i));
	}
	I = I / 2.0;
	if (isinf(I) || isnan(I))
		I = DBL_MAX_10_EXP;

	return I;
}

//Flory-Huggins Surface Activity Model (for AdsorptionReaction)
int FloryHuggins(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	AdsorptionReaction *dat = (AdsorptionReaction *) data;
	double logp = 0.0, invp = 0.0;
	double total = 0.0, lnact = 0.0;

	for (int i=0; i<F.rows(); i++)
		F.edit(i, 0, 1.0);

	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		total = total + pow(10.0, x(dat->getAdsorbIndex(i),0));
	}

	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		logp = 0.0;
		invp = 0.0;
		for (int j=0; j<dat->getNumberRxns(); j++)
		{
			double alpha = (dat->getAreaFactor(dat->getAdsorbIndex(i))/dat->getAreaFactor(dat->getAdsorbIndex(j))) - 1.0;
			logp = logp + ( (pow(10.0, x(dat->getAdsorbIndex(j),0))/total)/ (alpha + 1.0) );
			invp = logp;
		}
		logp = log(logp);
		invp = 1.0 / invp;
		lnact = 1.0 - logp - invp;

		F.edit(dat->getAdsorbIndex(i), 0, exp(lnact));
		if (isinf(F(dat->getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorbIndex(i),0)))
			F.edit(dat->getAdsorbIndex(i), 0, DBL_MAX);
	}

	return success;
}

//Flory-Huggins Surface Activity Model (for UnsteadyAdsorption)
int FloryHuggins_unsteady(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	UnsteadyAdsorption *dat = (UnsteadyAdsorption *) data;
	double logp = 0.0, invp = 0.0;
	double total = 0.0, lnact = 0.0;

	for (int i=0; i<F.rows(); i++)
		F.edit(i, 0, 1.0);

	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		total = total + pow(10.0, x(dat->getAdsorbIndex(i),0));
	}

	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		logp = 0.0;
		invp = 0.0;
		for (int j=0; j<dat->getNumberRxns(); j++)
		{
			double alpha = (dat->getAreaFactor(dat->getAdsorbIndex(i))/dat->getAreaFactor(dat->getAdsorbIndex(j))) - 1.0;
			logp = logp + ( (pow(10.0, x(dat->getAdsorbIndex(j),0))/total)/ (alpha + 1.0) );
			invp = logp;
		}
		logp = log(logp);
		invp = 1.0 / invp;
		lnact = 1.0 - logp - invp;

		F.edit(dat->getAdsorbIndex(i), 0, exp(lnact));
		if (isinf(F(dat->getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorbIndex(i),0)))
			F.edit(dat->getAdsorbIndex(i), 0, DBL_MAX);
	}

	return success;
}

//Flory-Huggins Surface Activity Model (for MultiligandAdsorption)
int FloryHuggins_multiligand(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;

	MultiligandAdsorption *dat = (MultiligandAdsorption *) data;
	double logp = 0.0, invp = 0.0;
	double total = 0.0, lnact = 0.0;

	for (int i=0; i<F.rows(); i++)
		F.edit(i, 0, 1.0);

	for (int l=0; l<dat->getNumberLigands(); l++)
	{
		for (int i=0; i<dat->getAdsorptionObject(l).getNumberRxns(); i++)
		{
			total = total + pow(10.0, x(dat->getAdsorptionObject(l).getAdsorbIndex(i),0));
		}
	}

	//Main loop to calculate activities
	for (int l=0; l<dat->getNumberLigands(); l++)
	{
		for (int i=0; i<dat->getAdsorptionObject(l).getNumberRxns(); i++)
		{
			logp = 0.0;
			invp = 0.0;

			for (int k=0; k<dat->getNumberLigands(); k++)
			{
				for (int j=0; j<dat->getAdsorptionObject(k).getNumberRxns(); j++)
				{
					double alpha = (dat->getAdsorptionObject(l).getAreaFactor(dat->getAdsorptionObject(l).getAdsorbIndex(i))/dat->getAdsorptionObject(k).getAreaFactor(dat->getAdsorptionObject(k).getAdsorbIndex(j))) - 1.0;
					logp = logp + ( (pow(10.0, x(dat->getAdsorptionObject(k).getAdsorbIndex(j),0))/total)/ (alpha + 1.0) );
					invp = logp;
				}
			}

			logp = log(logp);
			invp = 1.0 / invp;
			lnact = 1.0 - logp - invp;

			F.edit(dat->getAdsorptionObject(l).getAdsorbIndex(i), 0, exp(lnact));
			if (isinf(F(dat->getAdsorptionObject(l).getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorptionObject(l).getAdsorbIndex(i),0)))
				F.edit(dat->getAdsorptionObject(l).getAdsorbIndex(i), 0, DBL_MAX);
		}
	}

	return success;
}

//Flory-Huggins Surface Activity Model (for ChemisorptionReaction)
int FloryHuggins_chemi(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	ChemisorptionReaction *dat = (ChemisorptionReaction *) data;
	double logp = 0.0, invp = 0.0;
	double total = 0.0, lnact = 0.0;

	for (int i=0; i<F.rows(); i++)
		F.edit(i, 0, 1.0);

	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		total = total + pow(10.0, x(dat->getAdsorbIndex(i),0));
	}

	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		logp = 0.0;
		invp = 0.0;
		for (int j=0; j<dat->getNumberRxns(); j++)
		{
			double alpha = (dat->getAreaFactor(dat->getAdsorbIndex(i))/dat->getAreaFactor(dat->getAdsorbIndex(j))) - 1.0;
			logp = logp + ( (pow(10.0, x(dat->getAdsorbIndex(j),0))/total)/ (alpha + 1.0) );
			invp = logp;
		}
		logp = log(logp);
		invp = 1.0 / invp;
		lnact = 1.0 - logp - invp;

		F.edit(dat->getAdsorbIndex(i), 0, exp(lnact));
		if (isinf(F(dat->getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorbIndex(i),0)))
			F.edit(dat->getAdsorbIndex(i), 0, DBL_MAX);
	}

	return success;
}

//Flory-Huggins Surface Activity Model (for MultiligandChemisorption)
int FloryHuggins_multichemi(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;

	MultiligandChemisorption *dat = (MultiligandChemisorption *) data;
	double logp = 0.0, invp = 0.0;
	double total = 0.0, lnact = 0.0;

	for (int i=0; i<F.rows(); i++)
		F.edit(i, 0, 1.0);

	for (int l=0; l<dat->getNumberLigands(); l++)
	{
		for (int i=0; i<dat->getChemisorptionObject(l).getNumberRxns(); i++)
		{
			total = total + pow(10.0, x(dat->getChemisorptionObject(l).getAdsorbIndex(i),0));
		}
	}

	//Main loop to calculate activities
	for (int l=0; l<dat->getNumberLigands(); l++)
	{
		for (int i=0; i<dat->getChemisorptionObject(l).getNumberRxns(); i++)
		{
			logp = 0.0;
			invp = 0.0;

			for (int k=0; k<dat->getNumberLigands(); k++)
			{
				for (int j=0; j<dat->getChemisorptionObject(k).getNumberRxns(); j++)
				{
					double alpha = (dat->getChemisorptionObject(l).getAreaFactor(dat->getChemisorptionObject(l).getAdsorbIndex(i))/dat->getChemisorptionObject(k).getAreaFactor(dat->getChemisorptionObject(k).getAdsorbIndex(j))) - 1.0;
					logp = logp + ( (pow(10.0, x(dat->getChemisorptionObject(k).getAdsorbIndex(j),0))/total)/ (alpha + 1.0) );
					invp = logp;
				}
			}

			logp = log(logp);
			invp = 1.0 / invp;
			lnact = 1.0 - logp - invp;

			F.edit(dat->getChemisorptionObject(l).getAdsorbIndex(i), 0, exp(lnact));
			if (isinf(F(dat->getChemisorptionObject(l).getAdsorbIndex(i),0)) || isnan(F(dat->getChemisorptionObject(l).getAdsorbIndex(i),0)))
				F.edit(dat->getChemisorptionObject(l).getAdsorbIndex(i), 0, DBL_MAX);
		}
	}

	return success;
}

//UNIQUAC Surface Activity Model (for AdsorptionReaction)
int UNIQUAC(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	AdsorptionReaction *dat = (AdsorptionReaction *) data;
	double total = 0.0;
	double Temp = 0.0;
	std::vector<double> r;
	std::vector<double> s;
	std::vector<double> l;
	std::vector<double> u;
	std::vector<double> frac;
	std::vector<double> theta;
	r.resize(dat->getNumberRxns());
	s.resize(dat->getNumberRxns());
	l.resize(dat->getNumberRxns());
	u.resize(dat->getNumberRxns());
	frac.resize(dat->getNumberRxns());
	theta.resize(dat->getNumberRxns());
	Temp = -dat->getReaction(0).Get_Energy()/(Rstd*log(pow(10.0,dat->getReaction(0).Get_Equilibrium())));

	//Loop to calculate total adsorption and calculate surface area and volume constants
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		total = total + pow(10.0, x(dat->getAdsorbIndex(i),0));
	}

	//Loop to calculate fractions and fraction dependent functions
	double rx_sum = 0.0;
	double sx_sum = 0.0;
	double lx_sum = 0.0;
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		frac[i] = pow(10.0, x(dat->getAdsorbIndex(i),0)) / total;
		r[i] = dat->getVolumeFactor(dat->getAdsorbIndex(i)) / VolumeSTD;
		s[i] = dat->getAreaFactor(dat->getAdsorbIndex(i)) / AreaSTD;
		l[i] = LengthFactor(CoordSTD, r[i], s[i]);
		u[i] = -dat->getReaction(i).Get_Energy()/s[i]/CoordSTD/CoordSTD;
		rx_sum = rx_sum + (r[i]*frac[i]);
		sx_sum = sx_sum + (s[i]*frac[i]);
		lx_sum = lx_sum + (l[i]*frac[i]);
	}

	//Loop to gather the thetas
	for (int i=0; i<dat->getNumberRxns(); i++)
		theta[i] = (s[i]*frac[i])/sx_sum;

	//i Loop to fill in activities
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		double theta_tau_i = 0.0;
		double theta_tau_rat = 0.0;
		double lnact = 0.0;

		//Inner j loop
		for (int j=0; j<dat->getNumberRxns(); j++)
		{
			theta_tau_i = theta_tau_i + ( theta[j]*exp(-(sqrt(fabs(u[j]*u[i])) - u[i])/(Rstd*Temp)) );

			double theta_tau_k = 0.0;
			// k loop
			for (int k=0; k<dat->getNumberRxns(); k++)
			{
				theta_tau_k = theta_tau_k + ( theta[k]*exp(-(sqrt(fabs(u[k]*u[j])) - u[j])/(Rstd*Temp)) );

			}// End k Loop

			theta_tau_rat = theta_tau_rat + ( (theta[j]*exp(-(sqrt(fabs(u[i]*u[j])) - u[j])/(Rstd*Temp)) )/theta_tau_k);

		}// End j loop

		lnact = log(r[i]/rx_sum) + ( (CoordSTD/2.0)*s[i]*log((s[i]/r[i])*(rx_sum/sx_sum)) ) + l[i] - ( r[i]*(lx_sum/rx_sum) ) - ( s[i]*log(theta_tau_i) ) + s[i] - ( s[i]*theta_tau_rat );
		F.edit(dat->getAdsorbIndex(i), 0, exp(lnact));

		if (isinf(F(dat->getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorbIndex(i),0)))
			F.edit(dat->getAdsorbIndex(i), 0, DBL_MAX);
		if (F(dat->getAdsorbIndex(i),0) <= sqrt(DBL_MIN))
			F.edit(dat->getAdsorbIndex(i), 0, sqrt(DBL_MIN));

	}// End i Loop

	return success;
}

//UNIQUAC Surface Activity Model (for UnsteadyAdsorption)
int UNIQUAC_unsteady(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	UnsteadyAdsorption *dat = (UnsteadyAdsorption *) data;
	double total = 0.0;
	double Temp = 0.0;
	std::vector<double> r;
	std::vector<double> s;
	std::vector<double> l;
	std::vector<double> u;
	std::vector<double> frac;
	std::vector<double> theta;
	r.resize(dat->getNumberRxns());
	s.resize(dat->getNumberRxns());
	l.resize(dat->getNumberRxns());
	u.resize(dat->getNumberRxns());
	frac.resize(dat->getNumberRxns());
	theta.resize(dat->getNumberRxns());
	Temp = -dat->getReaction(0).Get_Energy()/(Rstd*log(pow(10.0,dat->getReaction(0).Get_Equilibrium())));

	//Loop to calculate total adsorption and calculate surface area and volume constants
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		total = total + pow(10.0, x(dat->getAdsorbIndex(i),0));
	}

	//Loop to calculate fractions and fraction dependent functions
	double rx_sum = 0.0;
	double sx_sum = 0.0;
	double lx_sum = 0.0;
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		frac[i] = pow(10.0, x(dat->getAdsorbIndex(i),0)) / total;
		r[i] = dat->getVolumeFactor(dat->getAdsorbIndex(i)) / VolumeSTD;
		s[i] = dat->getAreaFactor(dat->getAdsorbIndex(i)) / AreaSTD;
		l[i] = LengthFactor(CoordSTD, r[i], s[i]);
		u[i] = -dat->getReaction(i).Get_Energy()/s[i]/CoordSTD/CoordSTD;
		rx_sum = rx_sum + (r[i]*frac[i]);
		sx_sum = sx_sum + (s[i]*frac[i]);
		lx_sum = lx_sum + (l[i]*frac[i]);
	}

	//Loop to gather the thetas
	for (int i=0; i<dat->getNumberRxns(); i++)
		theta[i] = (s[i]*frac[i])/sx_sum;

	//i Loop to fill in activities
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		double theta_tau_i = 0.0;
		double theta_tau_rat = 0.0;
		double lnact = 0.0;

		//Inner j loop
		for (int j=0; j<dat->getNumberRxns(); j++)
		{
			theta_tau_i = theta_tau_i + ( theta[j]*exp(-(sqrt(fabs(u[j]*u[i])) - u[i])/(Rstd*Temp)) );

			double theta_tau_k = 0.0;
			// k loop
			for (int k=0; k<dat->getNumberRxns(); k++)
			{
				theta_tau_k = theta_tau_k + ( theta[k]*exp(-(sqrt(fabs(u[k]*u[j])) - u[j])/(Rstd*Temp)) );

			}// End k Loop

			theta_tau_rat = theta_tau_rat + ( (theta[j]*exp(-(sqrt(fabs(u[i]*u[j])) - u[j])/(Rstd*Temp)) )/theta_tau_k);

		}// End j loop

		lnact = log(r[i]/rx_sum) + ( (CoordSTD/2.0)*s[i]*log((s[i]/r[i])*(rx_sum/sx_sum)) ) + l[i] - ( r[i]*(lx_sum/rx_sum) ) - ( s[i]*log(theta_tau_i) ) + s[i] - ( s[i]*theta_tau_rat );
		F.edit(dat->getAdsorbIndex(i), 0, exp(lnact));

		if (isinf(F(dat->getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorbIndex(i),0)))
			F.edit(dat->getAdsorbIndex(i), 0, DBL_MAX);
		if (F(dat->getAdsorbIndex(i),0) <= sqrt(DBL_MIN))
			F.edit(dat->getAdsorbIndex(i), 0, sqrt(DBL_MIN));

	}// End i Loop

	return success;
}

//UNIQUAC Surface Activity Model (for MultiligandAdsorption)
int UNIQUAC_multiligand(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;

	MultiligandAdsorption *dat = (MultiligandAdsorption *) data;
	double total = 0.0;
	double Temp = 0.0;
	std::vector< std::vector<double> > r;
	std::vector< std::vector<double> > s;
	std::vector< std::vector<double> > l;
	std::vector< std::vector<double> > u;
	std::vector< std::vector<double> > frac;
	std::vector< std::vector<double> > theta;

	r.resize(dat->getNumberLigands());
	s.resize(dat->getNumberLigands());
	l.resize(dat->getNumberLigands());
	u.resize(dat->getNumberLigands());
	frac.resize(dat->getNumberLigands());
	theta.resize(dat->getNumberLigands());

	for (int k=0; k<dat->getNumberLigands(); k++)
	{
		r[k].resize(dat->getAdsorptionObject(k).getNumberRxns());
		s[k].resize(dat->getAdsorptionObject(k).getNumberRxns());
		l[k].resize(dat->getAdsorptionObject(k).getNumberRxns());
		u[k].resize(dat->getAdsorptionObject(k).getNumberRxns());
		frac[k].resize(dat->getAdsorptionObject(k).getNumberRxns());
		theta[k].resize(dat->getAdsorptionObject(k).getNumberRxns());
	}
	Temp = -dat->getAdsorptionObject(0).getReaction(0).Get_Energy()/(Rstd*log(pow(10.0,dat->getAdsorptionObject(0).getReaction(0).Get_Equilibrium())));

	//Loop to calculate total adsorption and calculate surface area and volume constants
	for (int k=0; k<dat->getNumberLigands(); k++)
	{
		for (int i=0; i<dat->getAdsorptionObject(k).getNumberRxns(); i++)
		{
			total = total + pow(10.0, x(dat->getAdsorptionObject(k).getAdsorbIndex(i),0));
		}
	}

	//Loop to calculate fractions and fraction dependent functions
	double rx_sum = 0.0;
	double sx_sum = 0.0;
	double lx_sum = 0.0;
	for (int k=0; k<dat->getNumberLigands(); k++)
	{
		for (int i=0; i<dat->getAdsorptionObject(k).getNumberRxns(); i++)
		{
			frac[k][i] = pow(10.0, x(dat->getAdsorptionObject(k).getAdsorbIndex(i),0)) / total;
			r[k][i] = dat->getAdsorptionObject(k).getVolumeFactor(dat->getAdsorptionObject(k).getAdsorbIndex(i)) / VolumeSTD;
			s[k][i] = dat->getAdsorptionObject(k).getAreaFactor(dat->getAdsorptionObject(k).getAdsorbIndex(i)) / AreaSTD;
			l[k][i] = LengthFactor(CoordSTD, r[k][i], s[k][i]);
			u[k][i] = -dat->getAdsorptionObject(k).getReaction(i).Get_Energy()/s[k][i]/CoordSTD/CoordSTD;
			rx_sum = rx_sum + (r[k][i]*frac[k][i]);
			sx_sum = sx_sum + (s[k][i]*frac[k][i]);
			lx_sum = lx_sum + (l[k][i]*frac[k][i]);
		}
	}

	//Loop to gather the thetas
	for (int k=0; k<dat->getNumberLigands(); k++)
		for (int i=0; i<dat->getAdsorptionObject(k).getNumberRxns(); i++)
			theta[k][i] = (s[k][i]*frac[k][i])/sx_sum;

	//i Loop to fill in activities
	for (int n=0; n<dat->getNumberLigands(); n++)
	{
		for (int i=0; i<dat->getAdsorptionObject(n).getNumberRxns(); i++)
		{
			double theta_tau_i = 0.0;
			double theta_tau_rat = 0.0;
			double lnact = 0.0;

			//Loop on m ligands....(note: all i's are associated with n's, while j's would be associated with m's)
			for (int m=0; m<dat->getNumberLigands(); m++)
			{
				//Inner j loop
				for (int j=0; j<dat->getAdsorptionObject(m).getNumberRxns(); j++)
				{
					theta_tau_i = theta_tau_i + ( theta[m][j]*exp(-(sqrt(fabs(u[m][j]*u[n][i])) - u[n][i])/(Rstd*Temp)) );

					double theta_tau_k = 0.0;

					//Loop on p ligands...(note: all i's are associated with n's, while j's would be associated with m's and k's with p's)
					for (int p=0; p<dat->getNumberLigands(); p++)
					{
						// k loop
						for (int k=0; k<dat->getAdsorptionObject(p).getNumberRxns(); k++)
						{
							theta_tau_k = theta_tau_k + ( theta[p][k]*exp(-(sqrt(fabs(u[p][k]*u[m][j])) - u[m][j])/(Rstd*Temp)) );

						}// End k Loop

					}//End p Loop

					theta_tau_rat = theta_tau_rat + ( (theta[m][j]*exp(-(sqrt(fabs(u[n][i]*u[m][j])) - u[m][j])/(Rstd*Temp)) )/theta_tau_k);

				}// End j loop

			}//END m loop

			lnact = log(r[n][i]/rx_sum) + ( (CoordSTD/2.0)*s[n][i]*log((s[n][i]/r[n][i])*(rx_sum/sx_sum)) ) + l[n][i] - ( r[n][i]*(lx_sum/rx_sum) ) - ( s[n][i]*log(theta_tau_i) ) + s[n][i] - ( s[n][i]*theta_tau_rat );
			F.edit(dat->getAdsorptionObject(n).getAdsorbIndex(i), 0, exp(lnact));

			if (isinf(F(dat->getAdsorptionObject(n).getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorptionObject(n).getAdsorbIndex(i),0)))
				F.edit(dat->getAdsorptionObject(n).getAdsorbIndex(i), 0, DBL_MAX);
			if (F(dat->getAdsorptionObject(n).getAdsorbIndex(i),0) <= sqrt(DBL_MIN))
				F.edit(dat->getAdsorptionObject(n).getAdsorbIndex(i), 0, sqrt(DBL_MIN));

		}// End i Loop
	}

	return success;
}

//UNIQUAC Surface Activity Model (for ChemisorptionReaction)
int UNIQUAC_chemi(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	ChemisorptionReaction *dat = (ChemisorptionReaction *) data;
	double total = 0.0;
	double Temp = 0.0;
	std::vector<double> r;
	std::vector<double> s;
	std::vector<double> l;
	std::vector<double> u;
	std::vector<double> frac;
	std::vector<double> theta;
	r.resize(dat->getNumberRxns());
	s.resize(dat->getNumberRxns());
	l.resize(dat->getNumberRxns());
	u.resize(dat->getNumberRxns());
	frac.resize(dat->getNumberRxns());
	theta.resize(dat->getNumberRxns());
	Temp = -dat->getReaction(0).Get_Energy()/(Rstd*log(pow(10.0,dat->getReaction(0).Get_Equilibrium())));

	//Loop to calculate total adsorption and calculate surface area and volume constants
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		total = total + pow(10.0, x(dat->getAdsorbIndex(i),0));
	}

	//Loop to calculate fractions and fraction dependent functions
	double rx_sum = 0.0;
	double sx_sum = 0.0;
	double lx_sum = 0.0;
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		frac[i] = pow(10.0, x(dat->getAdsorbIndex(i),0)) / total;
		r[i] = dat->getVolumeFactor(dat->getAdsorbIndex(i)) / VolumeSTD;
		s[i] = dat->getAreaFactor(dat->getAdsorbIndex(i)) / AreaSTD;
		l[i] = LengthFactor(CoordSTD, r[i], s[i]);
		u[i] = -dat->getReaction(i).Get_Energy()/s[i]/CoordSTD/CoordSTD;
		rx_sum = rx_sum + (r[i]*frac[i]);
		sx_sum = sx_sum + (s[i]*frac[i]);
		lx_sum = lx_sum + (l[i]*frac[i]);
	}

	//Loop to gather the thetas
	for (int i=0; i<dat->getNumberRxns(); i++)
		theta[i] = (s[i]*frac[i])/sx_sum;

	//i Loop to fill in activities
	for (int i=0; i<dat->getNumberRxns(); i++)
	{
		double theta_tau_i = 0.0;
		double theta_tau_rat = 0.0;
		double lnact = 0.0;

		//Inner j loop
		for (int j=0; j<dat->getNumberRxns(); j++)
		{
			theta_tau_i = theta_tau_i + ( theta[j]*exp(-(sqrt(fabs(u[j]*u[i])) - u[i])/(Rstd*Temp)) );

			double theta_tau_k = 0.0;
			// k loop
			for (int k=0; k<dat->getNumberRxns(); k++)
			{
				theta_tau_k = theta_tau_k + ( theta[k]*exp(-(sqrt(fabs(u[k]*u[j])) - u[j])/(Rstd*Temp)) );

			}// End k Loop

			theta_tau_rat = theta_tau_rat + ( (theta[j]*exp(-(sqrt(fabs(u[i]*u[j])) - u[j])/(Rstd*Temp)) )/theta_tau_k);

		}// End j loop

		lnact = log(r[i]/rx_sum) + ( (CoordSTD/2.0)*s[i]*log((s[i]/r[i])*(rx_sum/sx_sum)) ) + l[i] - ( r[i]*(lx_sum/rx_sum) ) - ( s[i]*log(theta_tau_i) ) + s[i] - ( s[i]*theta_tau_rat );
		F.edit(dat->getAdsorbIndex(i), 0, exp(lnact));

		if (isinf(F(dat->getAdsorbIndex(i),0)) || isnan(F(dat->getAdsorbIndex(i),0)))
			F.edit(dat->getAdsorbIndex(i), 0, DBL_MAX);
		if (F(dat->getAdsorbIndex(i),0) <= sqrt(DBL_MIN))
			F.edit(dat->getAdsorbIndex(i), 0, sqrt(DBL_MIN));

	}// End i Loop

	return success;
}

//UNIQUAC Surface Activity Model (for MultiligandChemisorption)
int UNIQUAC_multichemi(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;

	MultiligandChemisorption *dat = (MultiligandChemisorption *) data;
	double total = 0.0;
	double Temp = 0.0;
	std::vector< std::vector<double> > r;
	std::vector< std::vector<double> > s;
	std::vector< std::vector<double> > l;
	std::vector< std::vector<double> > u;
	std::vector< std::vector<double> > frac;
	std::vector< std::vector<double> > theta;

	r.resize(dat->getNumberLigands());
	s.resize(dat->getNumberLigands());
	l.resize(dat->getNumberLigands());
	u.resize(dat->getNumberLigands());
	frac.resize(dat->getNumberLigands());
	theta.resize(dat->getNumberLigands());

	for (int k=0; k<dat->getNumberLigands(); k++)
	{
		r[k].resize(dat->getChemisorptionObject(k).getNumberRxns());
		s[k].resize(dat->getChemisorptionObject(k).getNumberRxns());
		l[k].resize(dat->getChemisorptionObject(k).getNumberRxns());
		u[k].resize(dat->getChemisorptionObject(k).getNumberRxns());
		frac[k].resize(dat->getChemisorptionObject(k).getNumberRxns());
		theta[k].resize(dat->getChemisorptionObject(k).getNumberRxns());
	}
	Temp = -dat->getChemisorptionObject(0).getReaction(0).Get_Energy()/(Rstd*log(pow(10.0,dat->getChemisorptionObject(0).getReaction(0).Get_Equilibrium())));

	//Loop to calculate total adsorption and calculate surface area and volume constants
	for (int k=0; k<dat->getNumberLigands(); k++)
	{
		for (int i=0; i<dat->getChemisorptionObject(k).getNumberRxns(); i++)
		{
			total = total + pow(10.0, x(dat->getChemisorptionObject(k).getAdsorbIndex(i),0));
		}
	}

	//Loop to calculate fractions and fraction dependent functions
	double rx_sum = 0.0;
	double sx_sum = 0.0;
	double lx_sum = 0.0;
	for (int k=0; k<dat->getNumberLigands(); k++)
	{
		for (int i=0; i<dat->getChemisorptionObject(k).getNumberRxns(); i++)
		{
			frac[k][i] = pow(10.0, x(dat->getChemisorptionObject(k).getAdsorbIndex(i),0)) / total;
			r[k][i] = dat->getChemisorptionObject(k).getVolumeFactor(dat->getChemisorptionObject(k).getAdsorbIndex(i)) / VolumeSTD;
			s[k][i] = dat->getChemisorptionObject(k).getAreaFactor(dat->getChemisorptionObject(k).getAdsorbIndex(i)) / AreaSTD;
			l[k][i] = LengthFactor(CoordSTD, r[k][i], s[k][i]);
			u[k][i] = -dat->getChemisorptionObject(k).getReaction(i).Get_Energy()/s[k][i]/CoordSTD/CoordSTD;
			rx_sum = rx_sum + (r[k][i]*frac[k][i]);
			sx_sum = sx_sum + (s[k][i]*frac[k][i]);
			lx_sum = lx_sum + (l[k][i]*frac[k][i]);
		}
	}

	//Loop to gather the thetas
	for (int k=0; k<dat->getNumberLigands(); k++)
		for (int i=0; i<dat->getChemisorptionObject(k).getNumberRxns(); i++)
			theta[k][i] = (s[k][i]*frac[k][i])/sx_sum;

	//i Loop to fill in activities
	for (int n=0; n<dat->getNumberLigands(); n++)
	{
		for (int i=0; i<dat->getChemisorptionObject(n).getNumberRxns(); i++)
		{
			double theta_tau_i = 0.0;
			double theta_tau_rat = 0.0;
			double lnact = 0.0;

			//Loop on m ligands....(note: all i's are associated with n's, while j's would be associated with m's)
			for (int m=0; m<dat->getNumberLigands(); m++)
			{
				//Inner j loop
				for (int j=0; j<dat->getChemisorptionObject(m).getNumberRxns(); j++)
				{
					theta_tau_i = theta_tau_i + ( theta[m][j]*exp(-(sqrt(fabs(u[m][j]*u[n][i])) - u[n][i])/(Rstd*Temp)) );

					double theta_tau_k = 0.0;

					//Loop on p ligands...(note: all i's are associated with n's, while j's would be associated with m's and k's with p's)
					for (int p=0; p<dat->getNumberLigands(); p++)
					{
						// k loop
						for (int k=0; k<dat->getChemisorptionObject(p).getNumberRxns(); k++)
						{
							theta_tau_k = theta_tau_k + ( theta[p][k]*exp(-(sqrt(fabs(u[p][k]*u[m][j])) - u[m][j])/(Rstd*Temp)) );

						}// End k Loop

					}//End p Loop

					theta_tau_rat = theta_tau_rat + ( (theta[m][j]*exp(-(sqrt(fabs(u[n][i]*u[m][j])) - u[m][j])/(Rstd*Temp)) )/theta_tau_k);

				}// End j loop

			}//END m loop

			lnact = log(r[n][i]/rx_sum) + ( (CoordSTD/2.0)*s[n][i]*log((s[n][i]/r[n][i])*(rx_sum/sx_sum)) ) + l[n][i] - ( r[n][i]*(lx_sum/rx_sum) ) - ( s[n][i]*log(theta_tau_i) ) + s[n][i] - ( s[n][i]*theta_tau_rat );
			F.edit(dat->getChemisorptionObject(n).getAdsorbIndex(i), 0, exp(lnact));

			if (isinf(F(dat->getChemisorptionObject(n).getAdsorbIndex(i),0)) || isnan(F(dat->getChemisorptionObject(n).getAdsorbIndex(i),0)))
				F.edit(dat->getChemisorptionObject(n).getAdsorbIndex(i), 0, DBL_MAX);
			if (F(dat->getChemisorptionObject(n).getAdsorbIndex(i),0) <= sqrt(DBL_MIN))
				F.edit(dat->getChemisorptionObject(n).getAdsorbIndex(i), 0, sqrt(DBL_MIN));

		}// End i Loop
	}

	return success;
}

//Default Activity function is ideal solution
int ideal_solution (const Matrix<double>& x, Matrix<double> &F, const void *data)
{
	int success = 0;
	for (int i = 0; i<F.rows(); i++)
	{
		F.edit(i,0,1.0);
	}
	return success;
}

//Default non-ideal solution will use the Davies Equations
int Davies_equation (const Matrix<double>& x, Matrix<double> &F, const void *data)
{
	int success = 0;
	SHARK_DATA *dat = (SHARK_DATA *) data;
	double ionic_strength = calculate_ionic_strength(x,dat->MasterList);

	//Calculate log10(activity) from ionic_strength and variable A
	double a = 1.82E+6 * pow((dat->dielectric_const*dat->temperature), (-3.0/2.0));
	double log_gama = 0.0;
	for (int i=0; i<dat->numvar; i++)
	{
		log_gama = -a * dat->MasterList.charge(i) * dat->MasterList.charge(i) * ( (sqrt(ionic_strength)/(1.0+sqrt(ionic_strength))) - (0.2*ionic_strength));
		F.edit(i, 0, pow(10.0, log_gama));
		if (isinf(F(i,0)) || isnan(F(i,0)))
			F.edit(i, 0, DBL_MAX);
	}

	return success;
}

//Debye-Huckel Equation for activity coefficients
int DebyeHuckel_equation (const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	SHARK_DATA *dat = (SHARK_DATA *) data;
	double ionic_strength = calculate_ionic_strength(x,dat->MasterList);

	//Calculate log10(activity) from ionic_strength and variable A
	double a = 1.82E+6 * pow((dat->dielectric_const*dat->temperature), (-3.0/2.0));
	double log_gama = 0.0;
	for (int i=0; i<dat->numvar; i++)
	{
		log_gama = -a * dat->MasterList.charge(i) * dat->MasterList.charge(i) * sqrt(ionic_strength);
		F.edit(i, 0, pow(10.0, log_gama));
		if (isinf(F(i,0)) || isnan(F(i,0)))
			F.edit(i, 0, DBL_MAX);
	}

	return success;
}

/*
int Sit_equation (const Matrix<double>& x, Matrix<double> &F, const void *data)
{
	int success = 0;
	double K=0.0;
	SHARK_DATA *dat = (SHARK_DATA *) data;
	double ionic_strength = calculate_ionic_strength(x,dat->MasterList);

	//Calculate log10(activity) from ionic_strength and interaction coefficient

	double b[dat->numvar][dat->numvar], log_gama[dat->numvar];

	for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            b[i][j]=0.0;
        }

    }

	b[0][1]=0.0;
	b[0][2]=0.0;
	b[0][3]=0.038;
	b[0][4]=0.0;
	b[0][5]=0.0;
	b[0][6]=0.0;
	b[0][7]=0.0;
	b[0][8]=0.0;
	b[0][9]=0.0;
	b[0][10]=0.0;
	b[0][11]=0.0;
	b[0][12]=0.0;
	b[0][13]=0.0;
	b[0][14]=0.0;
	b[0][15]=0.0;
	b[0][16]=0.0;
	b[0][17]=0.0;
	b[0][18]=0.0;
	b[0][19]=0.0;
	b[0][20]=0.0;
	b[0][21]=0.125;
	b[0][22]=0.0;
	b[0][23]=0.0;
	b[0][24]=0.0;
	b[0][25]=0.0;
	b[3][4]=0.0;
	b[3][5]=0.0;
	b[3][6]=0.0;
	b[3][7]=0.0;
	b[3][8]=-0.08;
	b[3][9]=0.0;
	b[3][10]=0.0;
	b[3][11]=0.0;
	b[3][12]=0.0;
	b[3][13]=-0.09;
	b[3][14]=0.0;
	b[3][15]=0.0;
	b[3][16]=0.0;
	b[3][17]=-0.02;
	b[3][18]=-0.01;
	b[3][19]=0.0;
	b[3][20]=0.0;
	b[3][21]=0.0;
	b[3][22]=0.0;
	b[3][23]=0.0;
	b[3][24]=0.0;
	b[3][25]=0.0;




	for (int i=0; i<dat->numvar; i++)
    {
        for (int j=0; j<dat->numvar; j++)
        {
            if (i==j)
               break;
            b[i][j]=b[j][i];
        }
    }


	for (int i=0; i<dat->numvar; i++)
	{
        if (dat->MasterList.charge(i)==0)
            log_gama[i]=K*ionic_strength;
        else
        {
	    double log_gama0=0.0;
		log_gama0 = -dat->MasterList.charge(i)*dat->MasterList.charge(i) * 0.51 * sqrt(ionic_strength)/(1+1.5*sqrt(ionic_strength));

		//Species sum for molefractions
		double sum = 0.0;
		for (int j=0; j<dat->numvar; j++)
		{
			sum = sum + b[i][j]*pow(10,x(i,0));
		}
		log_gama[i] = log_gama0 + sum;

		F.edit(i, 0, pow(10.0, log_gama[i]));
		if (isinf(F(i,0)) || isnan(F(i,0)))
			F.edit(i, 0, DBL_MAX);
        }
	}

	return success;
}

int pitzer_equation (const Matrix<double>& x, Matrix<double> &F, const void *data)
{
	int success = 0;

	SHARK_DATA *dat = (SHARK_DATA *) data;
	double ionic_strength = calculate_ionic_strength(x,dat->MasterList);

	//Calculate log10(activity) from ionic_strength and interaction coefficient
	double log_gama[dat->numvar], Z=0, p0=0.0, p1=0.0, p2=0.0, p3=0.0, p4=0.0, p5=0.0, p6=0.0, F0=0.0, F1=0.0, F2=0.0;
	double F3=0.0, Fa=0.0, alpha=0.0, alpha1=0.0, alpha2=0.0, G1=0.0, G2=0.0, Gminus1=0.0, Gminus2=0.0;
	double minx1=0.0, minx2=0.0, minx3=0.0;
	int minum1=0, minum2=0, minum3=0;
	double beta0[dat->numvar][dat->numvar], beta1[dat->numvar][dat->numvar],beta2[dat->numvar][dat->numvar], Cphi[dat->numvar][dat->numvar],C[dat->numvar][dat->numvar],B[dat->numvar][dat->numvar],Bminus[dat->numvar][dat->numvar];
	double xMN[dat->numvar][dat->numvar],xMM[dat->numvar][dat->numvar],xNN[dat->numvar][dat->numvar],zeta[dat->numvar][dat->numvar], phi[dat->numvar][dat->numvar],phiminus[dat->numvar][dat->numvar], ezeta[dat->numvar][dat->numvar], ezetaminus[dat->numvar][dat->numvar];
	double xminus1[dat->numvar][64],xminus2[dat->numvar][64], xminus3[dat->numvar][64], jf[64],jminus[64];
	double psi[dat->numvar][dat->numvar][dat->numvar];

    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            beta0[i][j]=0.0;
            beta1[i][j]=0.0;
            beta2[i][j]=0.0;
            Cphi[i][j]=0.0;
            zeta[i][j]=0.0;

        }
    }


    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            beta0[i][j]=beta0[j][i];
            beta1[i][j]=beta1[j][i];
            beta2[i][j]=beta2[j][i];
            Cphi[i][j]=Cphi[j][i];
            zeta[i][j]=zeta[j][i];
        }
    }

    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            for(int k=0; k<dat->numvar; k++)
            {
                psi[i][j][k]=0.0;
            }

        }
    }

    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            for(int k=0; k<dat->numvar; k++)
            {
                psi[i][j][k]=psi[j][i][k];
                psi[i][j][k]=psi[i][k][j];
                psi[i][j][k]=psi[j][k][i];
                psi[i][j][k]=psi[k][i][j];
                psi[i][j][k]=psi[k][j][i];


            }
        }
    }




    jf[0]=0.0000706;
    jminus[0]=0.0127;
    jf[1]=0.0002387;
    jminus[1]=0.0207;
    jf[2]=0.0004806;
    jminus[2]=0.0275;
    jf[3]=0.0007850;
    jminus[3]=0.0333;
    jf[4]=0.0011443;
    jminus[4]=0.0385;
    jf[5]=0.0015529;
    jminus[5]=0.0432;
    jf[6]=0.0020063;
    jminus[6]=0.0475;
    jf[7]=0.0025010;
    jminus[7]=0.0514;
    jf[8]=0.0030340;
    jminus[8]=0.0551;
    jf[9]=0.0039028;
    jminus[9]=0.0586;
    jf[10]=0.0048393;
    jminus[10]=0.0649;
    jf[11]=0.0061961;
    jminus[11]=0.0706;
    jf[12]=0.0076615;
    jminus[12]=0.0758;
    jf[13]=0.0092260;
    jminus[13]=0.0806;
    jf[14]=0.010882;
    jminus[14]=0.0850;
    jf[15]=0.014441;
    jminus[15]=0.0928;
    jf[16]=0.018295;
    jminus[16]=0.0997;
    jf[17]=0.022409;
    jminus[17]=0.1059;
    jf[18]=0.026775;
    jminus[18]=0.1114;
    jf[19]=0.031313;
    jminus[19]=0.1164;
    jf[20]=0.036061;
    jminus[20]=0.1210;
    jf[21]=0.040985;
    jminus[21]=0.1252;
    jf[22]=0.046070;
    jminus[22]=0.1297;
    jf[23]=0.051306;
    jminus[23]=0.1327;
    jf[24]=0.056680;
    jminus[24]=0.1360;
    jf[25]=0.085346;
    jminus[25]=0.1499;
    jf[26]=0.11644;
    jminus[26]=0.1605;
    jf[27]=0.14941;
    jminus[27]=0.1689;
    jf[28]=0.18390;
    jminus[28]=0.1758;
    jf[29]=0.21965;
    jminus[29]=0.1815;
    jf[30]=0.25645;
    jminus[30]=0.1864;
    jf[31]=0.29416;
    jminus[31]=0.1906;
    jf[32]=0.49283;
    jminus[32]=0.2053;
    jf[33]=0.70293;
    jminus[33]=0.2142;
    jf[34]=0.92035;
    jminus[34]=0.2202;
    jf[35]=1.14288;
    jminus[35]=0.2246;
    jf[36]=1.36918;
    jminus[36]=0.2279;
    jf[37]=1.59839;
    jminus[37]=0.2304;
    jf[38]=1.82990;
    jminus[38]=0.2325;
    jf[39]=2.06328;
    jminus[39]=0.2342;
    jf[40]=2.53446;
    jminus[40]=0.2368;
    jf[41]=3.48916;
    jminus[41]=0.2402;
    jf[42]=4.45453;
    jminus[42]=0.2423;
    jf[43]=5.57865;
    jminus[43]=0.2374;
    jf[44]=6.40378;
    jminus[44]=0.2447;
    jf[45]=7.38429;
    jminus[45]=0.2455;
    jf[46]=8.36745;
    jminus[46]=0.2461;
    jf[47]=9.35270;
    jminus[47]=0.2465;
    jf[48]=11.82248;
    jminus[48]=0.2474;
    jf[49]=14.29890;
    jminus[49]=0.2479;
    jf[50]=16.77979;
    jminus[50]=0.2483;
    jf[51]=19.26387;
    jminus[51]=0.2485;
    jf[52]=21.75033;
    jminus[52]=0.2487;
    jf[53]=24.23861;
    jminus[53]=0.2489;
    jf[54]=49.17099;
    jminus[54]=0.2496;
    jf[55]=99.11907;
    jminus[55]=0.2498;
    jf[56]=149.9520;
    jminus[56]=0.2499;
    jf[57]=199.08083;
    jminus[57]=0.2499;
    jf[58]=249.07101;
    jminus[58]=0.2500;
    jf[59]=499.04682;
    jminus[59]=0.2500;
    jf[60]=999.03028;
    jminus[60]=0.2500;
    jf[61]=1499.03028;
    jminus[61]=0.2500;
    jf[62]=1999.01925;
    jminus[62]=0.2500;
    jf[63]=2499.01659;
    jminus[63]=0.2500;

    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            if (dat->MasterList.charge(i)==0 || dat->MasterList.charge(j)==0)
            {
                ezeta[i][j]=0.0;
                ezetaminus[i][j]=0.0;
                phi[i][j]=0.0;
                phiminus[i][j]=0.0;
                break;
            }

            xMN[i][j]=fabs(6*dat->MasterList.charge(i)*dat->MasterList.charge(j)*0.392*sqrt(ionic_strength));
            xMM[i][j]=6*dat->MasterList.charge(i)*dat->MasterList.charge(i)*0.392*sqrt(ionic_strength);
            xNN[i][j]=6*dat->MasterList.charge(i)*dat->MasterList.charge(i)*0.392*sqrt(ionic_strength);

                xminus1[i][0]=fabs(xMN[i][j]-0.01);
                xminus1[i][1]=fabs(xMN[i][j]-0.02);
                xminus1[i][2]=fabs(xMN[i][j]-0.03);
                xminus1[i][3]=fabs(xMN[i][j]-0.04);
                xminus1[i][4]=fabs(xMN[i][j]-0.05);
                xminus1[i][5]=fabs(xMN[i][j]-0.06);
                xminus1[i][6]=fabs(xMN[i][j]-0.07);
                xminus1[i][7]=fabs(xMN[i][j]-0.08);
                xminus1[i][8]=fabs(xMN[i][j]-0.09);
                xminus1[i][9]=fabs(xMN[i][j]-0.10);
                xminus1[i][10]=fabs(xMN[i][j]-0.12);
                xminus1[i][11]=fabs(xMN[i][j]-0.14);
                xminus1[i][12]=fabs(xMN[i][j]-0.16);
                xminus1[i][13]=fabs(xMN[i][j]-0.18);
                xminus1[i][14]=fabs(xMN[i][j]-0.20);
                xminus1[i][15]=fabs(xMN[i][j]-0.24);
                xminus1[i][16]=fabs(xMN[i][j]-0.28);
                xminus1[i][17]=fabs(xMN[i][j]-0.32);
                xminus1[i][18]=fabs(xMN[i][j]-0.36);
                xminus1[i][19]=fabs(xMN[i][j]-0.40);
                xminus1[i][20]=fabs(xMN[i][j]-0.44);
                xminus1[i][21]=fabs(xMN[i][j]-0.48);
                xminus1[i][22]=fabs(xMN[i][j]-0.52);
                xminus1[i][23]=fabs(xMN[i][j]-0.56);
                xminus1[i][24]=fabs(xMN[i][j]-0.60);
                xminus1[i][25]=fabs(xMN[i][j]-0.80);
                xminus1[i][26]=fabs(xMN[i][j]-1.00);
                xminus1[i][27]=fabs(xMN[i][j]-1.20);
                xminus1[i][28]=fabs(xMN[i][j]-1.40);
                xminus1[i][29]=fabs(xMN[i][j]-1.60);
                xminus1[i][30]=fabs(xMN[i][j]-1.80);
                xminus1[i][31]=fabs(xMN[i][j]-2.00);
                xminus1[i][32]=fabs(xMN[i][j]-3.00);
                xminus1[i][33]=fabs(xMN[i][j]-4.00);
                xminus1[i][34]=fabs(xMN[i][j]-5.00);
                xminus1[i][35]=fabs(xMN[i][j]-6.00);
                xminus1[i][36]=fabs(xMN[i][j]-7.00);
                xminus1[i][37]=fabs(xMN[i][j]-8.00);
                xminus1[i][38]=fabs(xMN[i][j]-9.00);
                xminus1[i][39]=fabs(xMN[i][j]-10.00);
                xminus1[i][40]=fabs(xMN[i][j]-12.00);
                xminus1[i][41]=fabs(xMN[i][j]-16.00);
                xminus1[i][42]=fabs(xMN[i][j]-20.00);
                xminus1[i][43]=fabs(xMN[i][j]-24.00);
                xminus1[i][44]=fabs(xMN[i][j]-28.00);
                xminus1[i][45]=fabs(xMN[i][j]-32.00);
                xminus1[i][46]=fabs(xMN[i][j]-36.00);
                xminus1[i][47]=fabs(xMN[i][j]-40.00);
                xminus1[i][48]=fabs(xMN[i][j]-50.00);
                xminus1[i][49]=fabs(xMN[i][j]-60.00);
                xminus1[i][50]=fabs(xMN[i][j]-70.00);
                xminus1[i][51]=fabs(xMN[i][j]-80.00);
                xminus1[i][52]=fabs(xMN[i][j]-90.00);
                xminus1[i][53]=fabs(xMN[i][j]-100.00);
                xminus1[i][54]=fabs(xMN[i][j]-200.00);
                xminus1[i][55]=fabs(xMN[i][j]-400.00);
                xminus1[i][56]=fabs(xMN[i][j]-600.00);
                xminus1[i][57]=fabs(xMN[i][j]-800.00);
                xminus1[i][58]=fabs(xMN[i][j]-1000.00);
                xminus1[i][59]=fabs(xMN[i][j]-2000.00);
                xminus1[i][60]=fabs(xMN[i][j]-4000.00);
                xminus1[i][61]=fabs(xMN[i][j]-6000.00);
                xminus1[i][62]=fabs(xMN[i][j]-8000.00);
                xminus1[i][63]=fabs(xMN[i][j]-10000.00);
                minx1=xminus1[i][0];
                minum1=0;
                for(int k=0; k<64; k++)
                {
                    if (xminus1[i][k]<minx1)
                    {
                        minx1=xminus1[i][k];
                        minum1=k;
                    }
                    else
                        break;
                }

                xminus2[i][0]=fabs(xMM[i][j]-0.01);
                xminus2[i][1]=fabs(xMM[i][j]-0.02);
                xminus2[i][2]=fabs(xMM[i][j]-0.05);
                xminus2[i][5]=fabs(xMM[i][j]-0.06);
                xminus2[i][6]=fabs(xMM[i][j]-0.07);
                xminus2[i][7]=fabs(xMM[i][j]-0.08);
                xminus2[i][8]=fabs(xMM[i][j]-0.09);
                xminus2[i][9]=fabs(xMM[i][j]-0.10);
                xminus2[i][10]=fabs(xMM[i][j]-0.12);
                xminus2[i][11]=fabs(xMM[i][j]-0.14);
                xminus2[i][12]=fabs(xMM[i][j]-0.16);
                xminus2[i][13]=fabs(xMM[i][j]-0.18);
                xminus2[i][14]=fabs(xMM[i][j]-0.20);
                xminus2[i][15]=fabs(xMM[i][j]-0.24);
                xminus2[i][16]=fabs(xMM[i][j]-0.28);
                xminus2[i][17]=fabs(xMM[i][j]-0.32);
                xminus2[i][18]=fabs(xMM[i][j]-0.36);
                xminus2[i][19]=fabs(xMM[i][j]-0.40);
                xminus2[i][20]=fabs(xMM[i][j]-0.44);
                xminus2[i][21]=fabs(xMM[i][j]-0.48);
                xminus2[i][22]=fabs(xMM[i][j]-0.52);
                xminus2[i][23]=fabs(xMM[i][j]-0.56);
                xminus2[i][24]=fabs(xMM[i][j]-0.60);
                xminus2[i][25]=fabs(xMM[i][j]-0.80);
                xminus2[i][26]=fabs(xMM[i][j]-1.00);
                xminus2[i][27]=fabs(xMM[i][j]-1.20);
                xminus2[i][28]=fabs(xMM[i][j]-1.40);
                xminus2[i][29]=fabs(xMM[i][j]-1.60);
                xminus2[i][30]=fabs(xMM[i][j]-1.80);
                xminus2[i][31]=fabs(xMM[i][j]-2.00);
                xminus2[i][32]=fabs(xMM[i][j]-3.00);
                xminus2[i][33]=fabs(xMM[i][j]-4.00);
                xminus2[i][34]=fabs(xMM[i][j]-5.00);
                xminus2[i][35]=fabs(xMM[i][j]-6.00);
                xminus2[i][36]=fabs(xMM[i][j]-7.00);
                xminus2[i][37]=fabs(xMM[i][j]-8.00);
                xminus2[i][38]=fabs(xMM[i][j]-9.00);
                xminus2[i][39]=fabs(xMM[i][j]-10.00);
                xminus2[i][40]=fabs(xMM[i][j]-12.00);
                xminus2[i][41]=fabs(xMM[i][j]-16.00);
                xminus2[i][42]=fabs(xMM[i][j]-20.00);
                xminus2[i][43]=fabs(xMM[i][j]-24.00);
                xminus2[i][44]=fabs(xMM[i][j]-28.00);
                xminus2[i][45]=fabs(xMM[i][j]-32.00);
                xminus2[i][46]=fabs(xMM[i][j]-36.00);
                xminus2[i][47]=fabs(xMM[i][j]-40.00);
                xminus2[i][48]=fabs(xMM[i][j]-50.00);
                xminus2[i][49]=fabs(xMM[i][j]-60.00);
                xminus2[i][50]=fabs(xMM[i][j]-70.00);
                xminus2[i][51]=fabs(xMM[i][j]-80.00);
                xminus2[i][52]=fabs(xMM[i][j]-90.00);
                xminus2[i][53]=fabs(xMM[i][j]-100.00);
                xminus2[i][54]=fabs(xMM[i][j]-200.00);
                xminus2[i][55]=fabs(xMM[i][j]-400.00);
                xminus2[i][56]=fabs(xMM[i][j]-600.00);
                xminus2[i][57]=fabs(xMM[i][j]-800.00);
                xminus2[i][58]=fabs(xMM[i][j]-1000.00);
                xminus2[i][59]=fabs(xMM[i][j]-2000.00);
                xminus2[i][60]=fabs(xMM[i][j]-4000.00);
                xminus2[i][61]=fabs(xMM[i][j]-6000.00);
                xminus2[i][62]=fabs(xMM[i][j]-8000.00);
                xminus2[i][63]=fabs(xMM[i][j]-10000.00);
                minx2=xminus2[i][0];
                minum2=0;
                for(int k=0; k<64; k++)
                {
                    if (xminus2[i][k]<minx2)
                    {
                        minx2=xminus2[i][k];
                        minum2=k;
                    }
                    else
                        break;
                }

                xminus3[i][0]=fabs(xNN[i][j]-0.01);
                xminus3[i][1]=fabs(xNN[i][j]-0.02);
                xminus3[i][2]=fabs(xNN[i][j]-0.03);
                xminus3[i][3]=fabs(xNN[i][j]-0.04);
                xminus3[i][4]=fabs(xNN[i][j]-0.05);
                xminus3[i][5]=fabs(xNN[i][j]-0.06);
                xminus3[i][6]=fabs(xNN[i][j]-0.07);
                xminus3[i][7]=fabs(xNN[i][j]-0.08);
                xminus3[i][8]=fabs(xNN[i][j]-0.09);
                xminus3[i][9]=fabs(xNN[i][j]-0.10);
                xminus3[i][10]=fabs(xNN[i][j]-0.12);
                xminus3[i][11]=fabs(xNN[i][j]-0.14);
                xminus3[i][12]=fabs(xNN[i][j]-0.16);
                xminus3[i][13]=fabs(xNN[i][j]-0.18);
                xminus3[i][14]=fabs(xNN[i][j]-0.20);
                xminus3[i][15]=fabs(xNN[i][j]-0.24);
                xminus3[i][16]=fabs(xNN[i][j]-0.28);
                xminus3[i][17]=fabs(xNN[i][j]-0.32);
                xminus3[i][18]=fabs(xNN[i][j]-0.36);
                xminus3[i][19]=fabs(xNN[i][j]-0.40);
                xminus3[i][20]=fabs(xNN[i][j]-0.44);
                xminus3[i][21]=fabs(xNN[i][j]-0.48);
                xminus3[i][22]=fabs(xNN[i][j]-0.52);
                xminus3[i][23]=fabs(xNN[i][j]-0.56);
                xminus3[i][24]=fabs(xNN[i][j]-0.60);
                xminus3[i][25]=fabs(xNN[i][j]-0.80);
                xminus3[i][26]=fabs(xNN[i][j]-1.00);
                xminus3[i][27]=fabs(xNN[i][j]-1.20);
                xminus3[i][28]=fabs(xNN[i][j]-1.40);
                xminus3[i][29]=fabs(xNN[i][j]-1.60);
                xminus3[i][30]=fabs(xNN[i][j]-1.80);
                xminus3[i][31]=fabs(xNN[i][j]-2.00);
                xminus3[i][32]=fabs(xNN[i][j]-3.00);
                xminus3[i][33]=fabs(xNN[i][j]-4.00);
                xminus3[i][34]=fabs(xNN[i][j]-5.00);
                xminus3[i][35]=fabs(xNN[i][j]-6.00);
                xminus3[i][36]=fabs(xNN[i][j]-7.00);
                xminus3[i][37]=fabs(xNN[i][j]-8.00);
                xminus3[i][38]=fabs(xNN[i][j]-9.00);
                xminus3[i][39]=fabs(xNN[i][j]-10.00);
                xminus3[i][40]=fabs(xNN[i][j]-12.00);
                xminus3[i][41]=fabs(xNN[i][j]-16.00);
                xminus3[i][42]=fabs(xNN[i][j]-20.00);
                xminus3[i][43]=fabs(xNN[i][j]-24.00);
                xminus3[i][44]=fabs(xNN[i][j]-28.00);
                xminus3[i][45]=fabs(xNN[i][j]-32.00);
                xminus3[i][46]=fabs(xNN[i][j]-36.00);
                xminus3[i][47]=fabs(xNN[i][j]-40.00);
                xminus3[i][48]=fabs(xNN[i][j]-50.00);
                xminus3[i][49]=fabs(xNN[i][j]-60.00);
                xminus3[i][50]=fabs(xNN[i][j]-70.00);
                xminus3[i][51]=fabs(xNN[i][j]-80.00);
                xminus3[i][52]=fabs(xNN[i][j]-90.00);
                xminus3[i][53]=fabs(xNN[i][j]-100.00);
                xminus3[i][54]=fabs(xNN[i][j]-200.00);
                xminus3[i][55]=fabs(xNN[i][j]-400.00);
                xminus3[i][56]=fabs(xNN[i][j]-600.00);
                xminus3[i][57]=fabs(xNN[i][j]-800.00);
                xminus3[i][58]=fabs(xNN[i][j]-1000.00);
                xminus3[i][59]=fabs(xNN[i][j]-2000.00);
                xminus3[i][60]=fabs(xNN[i][j]-4000.00);
                xminus3[i][61]=fabs(xNN[i][j]-6000.00);
                xminus3[i][62]=fabs(xNN[i][j]-8000.00);
                xminus3[i][63]=fabs(xNN[i][j]-10000.00);
                minx3=xminus3[i][0];
                minum3=0;
                for(int k=0; k<64; k++)
                {
                    if (xminus3[i][k]<minx3)
                    {
                        minx3=xminus3[i][k];
                        minum3=k;
                    }
                    else
                        break;
                }

                ezeta[i][j]=((dat->MasterList.charge(i))*(dat->MasterList.charge(j))/(4*ionic_strength))*(jf[minum1]-0.5*jf[minum2]-0.5*jf[minum3]);
                ezetaminus[i][j]=-(ezeta[i][j]/ionic_strength)+(dat->MasterList.charge(i))*(dat->MasterList.charge(j))/(8*ionic_strength*ionic_strength)*(xMN[i][j]*jminus[minum1]-0.5*xMM[i][j]*jminus[minum2]-0.5*xNN[i][j]*jminus[minum3]);
                phi[i][j]=zeta[i][j]+ezeta[i][j];
                phiminus[i][j]=ezetaminus[i][j];

        }
    }



	for(int i=0; i<dat->numvar; i++)
	{
	    Z=pow(10.0,x(i,0))*fabs(dat->MasterList.charge(i));
	}

    F0=-0.392*(sqrt(ionic_strength)/(1.0+(1.2*sqrt(ionic_strength)))+(2/1.2)*log(1+(1.2*sqrt(ionic_strength))));

    for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)<0)
            {
              for(int j=0; j<dat->numvar; j++)
              {
                 if ((dat->MasterList.charge(i)<0) || (i==j))
                    break;
                 F1=F1+pow(10.0,x(i,0))*pow(10.0,x(j,0))*Bminus[i][j];
              }
            }
            else if (dat->MasterList.charge(i)==0)
               break;
            else
             for (int k=0; k<dat->numvar; k++)
             {
                 if ((dat->MasterList.charge(i)>0) || (i==k))
                    break;
                 F1=F1+pow(10.0,x(i,0))*pow(10.0,x(k,0))*Bminus[i][k];
             }
        }

    for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)<0)
               break;
               for (int j=0; j<dat->numvar; j++)
               {
                   if ((dat->MasterList.charge(i)<0) || i==j)
                    break;
                   F2=F2+pow(10.0,x(i,0))*pow(10.0,x(j,0))*phiminus[i][j];
               }
        }

    for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)>0)
               break;
               for (int j=0; j<dat->numvar; j++)
               {
                   if ((dat->MasterList.charge(i)<0) || i==j)
                    break;
                   F3=F3+pow(10.0,x(i,0))*pow(10.0,x(j,0))*phiminus[i][j];
               }
        }

    Fa=F0+F1+F2+F3;
    F1=0.0;
    F2=0.0;
    F3=0.0;


    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            if ((fabs(dat->MasterList.charge(i))==1) || (fabs(dat->MasterList.charge(j))==1))
            {
                alpha=2;
                G1=2*(1-(1+alpha*sqrt(ionic_strength))*exp(-alpha*sqrt(ionic_strength)))/(alpha*alpha*ionic_strength);
                Gminus1=-2*(1-(1+alpha*sqrt(ionic_strength)+0.5*alpha*alpha*ionic_strength)*exp(-alpha*sqrt(ionic_strength)))/(alpha*alpha*ionic_strength);
                B[i][j]=beta0[i][j]+beta1[i][j]*G1;
                Bminus[i][j]=beta1[i][j]*Gminus1/ionic_strength;
                alpha=0.0;
                G1=0.0;
                Gminus1=0.0;
            }
            else if ((fabs(dat->MasterList.charge(i))==2) && (fabs(dat->MasterList.charge(j))==2))
            {
                alpha1=1.4;
                alpha2=12;
                G1=2*(1-(1+alpha1*sqrt(ionic_strength))*exp(-alpha1*sqrt(ionic_strength)))/(alpha1*alpha1*ionic_strength);
                Gminus1=-2*(1-(1+alpha1*sqrt(ionic_strength)+0.5*alpha1*alpha1*ionic_strength)*exp(-alpha1*sqrt(ionic_strength)))/(alpha1*alpha1*ionic_strength);
                G2=2*(1-(1+alpha2*sqrt(ionic_strength))*exp(-alpha2*sqrt(ionic_strength)))/(alpha2*alpha2*ionic_strength);
                Gminus2=-2*(1-(1+alpha2*sqrt(ionic_strength)+0.5*alpha2*alpha2*ionic_strength)*exp(-alpha2*sqrt(ionic_strength)))/(alpha2*alpha2*ionic_strength);
                B[i][j]=beta0[i][j]+beta1[i][j]*G1+beta2[i][j]*G2;
                Bminus[i][j]=beta1[i][j]*Gminus1/ionic_strength+beta2[i][j]*Gminus2/ionic_strength;
                alpha1=0.0;
                alpha2=0.0;
                G1=0.0;
                Gminus1=0.0;
                G2=0.0;
                Gminus2=0.0;

            }
            else
            {
                alpha1=2.0;
                alpha2=50.0;
                G1=2*(1-(1+alpha1*sqrt(ionic_strength))*exp(-alpha1*sqrt(ionic_strength)))/(alpha1*alpha1*ionic_strength);
                Gminus1=-2*(1-(1+alpha1*sqrt(ionic_strength)+0.5*alpha1*alpha1*ionic_strength)*exp(-alpha1*sqrt(ionic_strength)))/(alpha1*alpha1*ionic_strength);
                G2=2*(1-(1+alpha2*sqrt(ionic_strength))*exp(-alpha2*sqrt(ionic_strength)))/(alpha2*alpha2*ionic_strength);
                Gminus2=-2*(1-(1+alpha2*sqrt(ionic_strength)+0.5*alpha2*alpha2*ionic_strength)*exp(-alpha2*sqrt(ionic_strength)))/(alpha2*alpha2*ionic_strength);
                B[i][j]=beta0[i][j]+beta1[i][j]*G1+beta2[i][j]*G2;
                Bminus[i][j]=beta1[i][j]*Gminus1/ionic_strength+beta2[i][j]*Gminus2/ionic_strength;
                alpha1=0.0;
                alpha2=0.0;
                G1=0.0;
                Gminus1=0.0;
                G2=0.0;
                Gminus2=0.0;
            }

        }
    }

    for(int i=0; i<dat->numvar; i++)
    {
        for(int j=0; j<dat->numvar; j++)
        {
            if (dat->MasterList.charge(i)==0 || dat->MasterList.charge(j)==0)
                break;
            C[i][j]=Cphi[i][j]/(2*sqrt(fabs(dat->MasterList.charge(i)*dat->MasterList.charge(j))));
        }
    }

	for(int n=0; n<dat->numvar; n++)
    {
        p0=dat->MasterList.charge(n)*dat->MasterList.charge(n)*Fa;
        if (dat->MasterList.charge(n)<0)
        {
        for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)<0)
              for(int j=0; j<dat->numvar; j++)
              {
                 if ((dat->MasterList.charge(i)<0) || (i==j))
                    break;
                 p4=p4+pow(10.0,x(i,0))*pow(10.0,x(j,0))*C[i][j];
              }
            else if (dat->MasterList.charge(i)==0)
               break;
            else
             for (int k=0; k<dat->numvar; k++)
             {
                 if ((dat->MasterList.charge(i)>0) || (i==k))
                    break;
                 p4=p4+pow(10.0,x(i,0))*pow(10.0,x(k,0))*C[i][k];
             }
        }


        for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)<0)
               break;
               for (int j=0; j<dat->numvar; j++)
               {
                   if ((dat->MasterList.charge(i)<0) || i==j)
                    break;
                   p3=p3+pow(10.0,x(i,0))*pow(10.0,x(j,0))*psi[i][j][n];
               }
        }

        for (int i=0; i<dat->numvar; i++)
        {
           if (dat->MasterList.charge(i)>0)
            break;
            for (int j=0; j<dat->numvar; j++)
            {
               if (dat->MasterList.charge(i)<0)
                break;
               p5=p5+pow(10.0,x(j,0))*psi[n][i][j];
            }
           p2=p2+pow(10.0,x(i,0))*(2*phi[n][i]+p5);
        }

        for (int i=0; i<dat->numvar; i++)
        {
            p6=2*B[i][n]+Z*C[i][n];
            if (dat->MasterList.charge(i)<0)
                break;
            p1=p1+pow(10.0,x(i,0))*p6;
        }
        log_gama[n]=0.43*(p0+p1+p2+p3+fabs(dat->MasterList.charge(n))*p4);
        p1=0.0;
        p2=0.0;
        p3=0.0;
        p4=0.0;
        p5=0.0;
        p6=0.0;
        }
        else if (dat->MasterList.charge(n)>0)
        {
            for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)<0)
              for(int j=0; j<dat->numvar; j++)
              {
                 if ((dat->MasterList.charge(i)<0) || (i==j))
                    break;
                 p4=p4+pow(10.0,x(i,0))*pow(10.0,x(j,0))*C[i][j];
              }
            else if (dat->MasterList.charge(i)==0)
               break;
            else
             for (int k=0; k<dat->numvar; k++)
             {
                 if ((dat->MasterList.charge(i)>0) || (i==k))
                    break;
                 p4=p4+pow(10.0,x(i,0))*pow(10.0,x(k,0))*C[i][k];
             }
        }

        for (int i=0; i<dat->numvar; i++)
        {
            if (dat->MasterList.charge(i)>0)
               break;
               for (int j=0; j<dat->numvar; j++)
               {
                   if ((dat->MasterList.charge(i)<0) || i==j)
                    break;
                   p3=p3+pow(10.0,x(i,0))*pow(10.0,x(j,0))*psi[i][j][n];
               }
        }
        for (int i=0; i<dat->numvar; i++)
        {
           if (dat->MasterList.charge(i)<0)
            break;
            for (int j=0; j<dat->numvar; j++)
            {
               if (dat->MasterList.charge(i)>0)
                break;
               p5=p5+pow(10.0,x(j,0))*psi[n][i][j];
            }
           p2=p2+pow(10.0,x(i,0))*(2*phi[n][i]+p5);
        }

         for (int i=0; i<dat->numvar; i++)
        {
            p6=2*B[n][i]+Z*C[n][i];
            if (dat->MasterList.charge(i)<0)
                break;
            p1=p1+pow(10.0,x(i,0))*p6;
        }

        log_gama[n]=0.43*(p0+p1+p2+p3+fabs(dat->MasterList.charge(n))*p4);
        p1=0.0;
        p2=0.0;
        p3=0.0;
        p4=0.0;
        p5=0.0;
        p6=0.0;
        //std::cout << "n = " << n << "\tlog_gama = " << log_gama[n] << std::endl;
        }
        else
            log_gama[n]=0.0;
        p0=0.0;
        F.edit(n, 0, pow(10.0, log_gama[n]));
        //F.edit(n, 0, pow(10.0, 0.0));
		if (isinf(F(n,0)) || isnan(F(n,0)))
			F.edit(n, 0, DBL_MAX);
    }



	return success;
}
*/

//Determine the surface activity function choosen by user
int surf_act_choice(const std::string &input)
{
	int act = IDEAL_ADS;

	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);

	if (copy == "ideal")
		act = IDEAL_ADS;
	else if (copy == "floryhuggins" || copy == "flory-huggins")
		act = FLORY_HUGGINS;
	else if (copy == "uniquac")
		act = UNIQUAC_ACT;
	else
		act = IDEAL_ADS;

	return act;
}

//Determine the activity function choosen by user
int act_choice(const std::string &input)
{
	int act = IDEAL;

	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);

	if (copy == "davies")
		act = DAVIES;
	else if (copy == "ideal")
		act = IDEAL;
	else if (copy == "debye-huckel")
		act = DEBYE_HUCKEL;
	else if (copy == "pitzer")
		act = PITZER;
	else if (copy == "sit")
		act = SIT;
	else
		act = IDEAL;

	return act;
}

//Determine the reactor type choosen by user
int reactor_choice(const std::string &input)
{
	int reactor = BATCH;

	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);

	if (copy == "batch")
		reactor = BATCH;
	else if (copy == "cstr")
		reactor = CSTR;
	else if (copy == "pfr")
		reactor = PFR;
	else
		reactor = BATCH;

	return reactor;
}

//Determine the line search option
bool linesearch_choice(const std::string &input)
{
	bool Bounce = false;

	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);

	if (copy == "standard")
		Bounce = false;
	else if (copy == "bounce" || copy == "bouncing")
		Bounce = true;
	else
		Bounce = false;

	return Bounce;
}

//Determine what linear solver to use
int linearsolve_choice(const std::string &input)
{
	int choice = GMRESRP;

	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);

	if (copy == "gmreslp")
		choice = GMRESLP;
	else if (copy == "pcg")
		choice = PCG;
	else if (copy == "bicgstab")
		choice = BiCGSTAB;
	else if (copy == "cgs")
		choice = CGS;
	else if (copy == "fom")
		choice = FOM;
	else if (copy == "gmresrp")
		choice = GMRESRP;
	else if (copy == "gcr")
		choice = GCR;
	else if (copy == "gmresr")
		choice = GMRESR;
	else if (copy == "kms")
		choice = KMS;
	else if (copy == "gmres")
		choice = GMRESRP;
	else if (copy == "qr")
		choice = QR;
	else
		choice = GMRESRP;

	return choice;
}

//Make a conversion from x to logx
int Convert2LogConcentration(const Matrix<double> &x, Matrix<double> &logx)
{
	int success = 0;
	for (int i=0; i<logx.rows(); i++)
	{
		if (x(i,0) <= 0.0)
			logx.edit(i, 0, -DBL_MAX_10_EXP);
		else
			logx.edit(i, 0, log10(x(i,0)));
	}
	return success;
}

//Make a conversion from logx to x
int Convert2Concentration(const Matrix<double> &logx, Matrix<double> &x)
{
	int success = 0;
	for (int i=0; i<x.rows(); i++)
	{
		if (logx(i,0) >= DBL_MAX_10_EXP)
			x.edit(i, 0, DBL_MAX);
		else if (logx(i,0) <= (-DBL_MAX_10_EXP+1.0))
			x.edit(i, 0, 0.0);
		else
			x.edit(i, 0, pow(10.0, logx(i,0)));
	}
	return success;
}

//Check and read the scenario portion of the yaml object
int read_scenario(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Find all required info from Scenario Document in Header vars_fun
	try
	{
		shark_dat->numvar = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["numvar"].getInt();
		shark_dat->num_ssr = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_ssr"].getInt();
		shark_dat->num_mbe = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_mbe"].getInt();
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}
	//Find all optional info from Scenario Document in Header vars_fun
	try
	{
		shark_dat->num_usr = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_usr"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_usr = 0;
	}
	try
	{
		shark_dat->num_ssao = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_ssao"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_ssao = 0;
	}
	shark_dat->num_ssar.resize(shark_dat->num_ssao);
	shark_dat->ss_ads_names.resize(shark_dat->num_ssao);
	try
	{
		shark_dat->num_usao = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_usao"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_usao = 0;
	}
	shark_dat->num_usar.resize(shark_dat->num_usao);
	shark_dat->us_ads_names.resize(shark_dat->num_usao);
	try
	{
		shark_dat->num_multi_ssao = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_multi_ssao"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_multi_ssao = 0;
	}
	shark_dat->num_multi_ssar.resize(shark_dat->num_multi_ssao);
	shark_dat->ssmulti_names.resize(shark_dat->num_multi_ssao);
	shark_dat->MultiAdsList.resize(shark_dat->num_multi_ssao);
	try
	{
		shark_dat->num_sschem = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_sschem"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_sschem = 0;
	}
	shark_dat->num_sschem_rxns.resize(shark_dat->num_sschem);
	shark_dat->ss_chem_names.resize(shark_dat->num_sschem);
	try
	{
		shark_dat->num_multi_sschem = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_multi_sschem"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_multi_sschem = 0;
	}
	shark_dat->num_multichem_rxns.resize(shark_dat->num_multi_sschem);
	shark_dat->ssmultichem_names.resize(shark_dat->num_multi_sschem);
	shark_dat->MultiChemList.resize(shark_dat->num_multi_sschem);
	try
	{
		shark_dat->num_other = shark_dat->yaml_object.getYamlWrapper()("Scenario")("vars_fun")["num_other"].getInt();
	}
	catch (std::out_of_range)
	{
		shark_dat->num_other = 0;
	}

	//Read through adsorption objects for all adsorption reactions
	if (shark_dat->num_ssao > 0)
	{
		int object_check;
		try
		{
			object_check = (int)shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_ads_objs").getSubMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (object_check != shark_dat->num_ssao)
		{
			mError(missing_information);
			return -1;
		}

		int obj = 0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_ads_objs").getSubMap())
		{
			try
			{
				shark_dat->num_ssar[obj] = x.second.getMap().getInt("num_rxns");
				shark_dat->ss_ads_names[obj] = x.second.getMap().getString("name");
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			obj++;
		}
	}

	//Read through unsteady adsorption objects for all adsorption reactions
	if (shark_dat->num_usao > 0)
	{
		int object_check;
		try
		{
			object_check = (int)shark_dat->yaml_object.getYamlWrapper()("Scenario")("us_ads_objs").getSubMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (object_check != shark_dat->num_usao)
		{
			mError(missing_information);
			return -1;
		}

		int obj = 0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("Scenario")("us_ads_objs").getSubMap())
		{
			try
			{
				shark_dat->num_usar[obj] = x.second.getMap().getInt("num_rxns");
				shark_dat->us_ads_names[obj] = x.second.getMap().getString("name");
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			obj++;
		}
	}

	//Read through all multiligand objects for the names of all adsorbents (and number of ligands each contains)
	if (shark_dat->num_multi_ssao > 0)
	{
		int object_check;
		try
		{
			object_check = (int)shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_multiads_objs").getSubMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (object_check != shark_dat->num_multi_ssao)
		{
			mError(missing_information);
			return -1;
		}

		int obj = 0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_multiads_objs").getSubMap())
		{
			try
			{
				shark_dat->MultiAdsList[obj].setAdsorbentName( x.second.getMap().getString("name") );
				shark_dat->num_multi_ssar[obj].resize(x.second.getMap().getInt("num_ligands"));
				shark_dat->ssmulti_names[obj].resize(x.second.getMap().getInt("num_ligands"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			obj++;
		}
	}

	//Read through chemisorption objects for all adsorption reactions
	if (shark_dat->num_sschem > 0)
	{
		int object_check;
		try
		{
			object_check = (int)shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_chem_objs").getSubMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (object_check != shark_dat->num_sschem)
		{
			mError(missing_information);
			return -1;
		}

		int obj = 0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_chem_objs").getSubMap())
		{
			try
			{
				shark_dat->num_sschem_rxns[obj] = x.second.getMap().getInt("num_rxns");
				shark_dat->ss_chem_names[obj] = x.second.getMap().getString("ligand");
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			obj++;
		}
	}

	//Read through all multiligand chemisorption objects for the names of all adsorbents (and number of ligands each contains)
	if (shark_dat->num_multi_sschem > 0)
	{
		int object_check;
		try
		{
			object_check = (int)shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_multichemi_objs").getSubMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (object_check != shark_dat->num_multi_sschem)
		{
			mError(missing_information);
			return -1;
		}

		int obj = 0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("Scenario")("ss_multichemi_objs").getSubMap())
		{
			try
			{
				shark_dat->MultiChemList[obj].setAdsorbentName( x.second.getMap().getString("name") );
				shark_dat->num_multichem_rxns[obj].resize(x.second.getMap().getInt("num_ligands"));
				shark_dat->ssmultichem_names[obj].resize(x.second.getMap().getInt("num_ligands"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			obj++;
		}
	}

	//All sys_data is optional, missing info will be replaced with defaults
	try
	{
		shark_dat->const_pH = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["const_pH"].getBool();
	} catch (std::out_of_range)
	{
		shark_dat->const_pH = false;
	}

	if (shark_dat->const_pH == true)
	{
		try
		{
			shark_dat->pH = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["pH"].getDouble();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
	}

	try
	{
		shark_dat->temperature = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["temp"].getDouble();
	} catch (std::out_of_range)
	{
		shark_dat->temperature = 298.15;
	}

	try
	{
		shark_dat->dielectric_const = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["dielec"].getDouble();
	} catch (std::out_of_range)
	{
		shark_dat->dielectric_const = 78.325;
	}

	try
	{
		shark_dat->relative_permittivity = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["rel_perm"].getDouble();
	} catch (std::out_of_range)
	{
		shark_dat->relative_permittivity = WaterRelPerm;
	}

	try
	{
		shark_dat->MasterList.set_alkalinity(shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["res_alk"].getDouble());
	}
	catch (std::out_of_range)
	{
		shark_dat->MasterList.set_alkalinity(0.0);
	}

	try
	{
		std::string act = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["act_fun"].getString();
		shark_dat->act_fun = act_choice(act);
	} catch (std::out_of_range)
	{
		shark_dat->act_fun = IDEAL;
	}

	try
	{
		shark_dat->volume = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["volume"].getDouble();
	} catch (std::out_of_range)
	{
		shark_dat->volume = 1.0;
	}

	try
	{
		std::string reactor = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["reactor"].getString();
		shark_dat->reactor_type = reactor_choice(reactor);
	} catch (std::out_of_range)
	{
		shark_dat->reactor_type = BATCH;
	}

	try
	{
		shark_dat->flow_rate = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["flow_rate"].getDouble();
	} catch (std::out_of_range)
	{
		shark_dat->flow_rate = 1.0;
	}

	try
	{
		shark_dat->xsec_area = shark_dat->yaml_object.getYamlWrapper()("Scenario")("sys_data")["xsec_area"].getDouble();
	} catch (std::out_of_range)
	{
		shark_dat->xsec_area = 1.0;
	}

	//Now read in the run_time Header
	try
	{
		shark_dat->steadystate = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["steady"].getBool();
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	try
	{
		shark_dat->ZeroInitialSolids = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["zero_initial_solids"].getBool();
	}
	catch (std::out_of_range)
	{
		shark_dat->ZeroInitialSolids = false;
	}

	if (shark_dat->num_usr == 0 && shark_dat->reactor_type == BATCH)
		shark_dat->steadystate = true;

	if (shark_dat->steadystate == true)
	{
		try
		{
			shark_dat->TemperatureCurve = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["temp_curve"].getBool();
		}
		catch (std::out_of_range)
		{
			shark_dat->TemperatureCurve = false;
		}

		try
		{
			shark_dat->SpeciationCurve = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["specs_curve"].getBool();
		}
		catch (std::out_of_range)
		{
			shark_dat->SpeciationCurve = false;
		}

		if (shark_dat->SpeciationCurve == true)
		{
			shark_dat->TemperatureCurve = false;
			shark_dat->pH = 1.0;
			try
			{
				shark_dat->pH_step = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["pH_step"].getDouble();
				if (shark_dat->pH_step < 0.01)
					shark_dat->pH_step = 0.01;
				if (shark_dat->pH_step > 1.0)
					shark_dat->pH_step = 1.0;
			}
			catch (std::out_of_range)
			{
				shark_dat->pH_step = 0.5;
			}
			try
			{
				shark_dat->pH = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["pH_start"].getDouble();
			}
			catch (std::out_of_range)
			{
				shark_dat->pH = 1.0;
			}
			try
			{
				shark_dat->pH_end = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["pH_end"].getDouble();
			}
			catch (std::out_of_range)
			{
				shark_dat->pH_end = 14.0;
			}
		}

		if (shark_dat->TemperatureCurve == true)
		{
			shark_dat->SpeciationCurve = false;
			try
			{
				shark_dat->start_temp = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["start_temp"].getDouble();
			}
			catch (std::out_of_range)
			{
				shark_dat->start_temp = 277.15;
			}
			try
			{
				shark_dat->end_temp = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["end_temp"].getDouble();
			}
			catch (std::out_of_range)
			{
				shark_dat->end_temp = 323.15;
			}
			try
			{
				shark_dat->temp_step = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["temp_step"].getDouble();
			}
			catch (std::out_of_range)
			{
				shark_dat->temp_step = 10.0;
			}

			shark_dat->temperature = shark_dat->start_temp;
		}
	}
	else
	{
		try
		{
			shark_dat->dt = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["dt"].getDouble();
			shark_dat->simulationtime = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["sim_time"].getDouble();
			shark_dat->t_out = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["t_out"].getDouble();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}

		try
		{
			shark_dat->TimeAdaptivity = shark_dat->yaml_object.getYamlWrapper()("Scenario")("run_time")["time_adapt"].getBool();
		}
		catch (std::out_of_range)
		{
			shark_dat->TimeAdaptivity = false;
		}
	}


	return success;
}

//Read in the information about each multisite adsorbent fiber (must be done before calling setup)
int read_multiligand_scenario(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_multi_ssao > 0)
	{
		//Loop through ligand names and match Document name with object name
		for (int i=0; i<shark_dat->num_multi_ssao; i++)
		{
			//Set the values of the fiber constants
			try
			{
				shark_dat->MultiAdsList[i].setSpecificArea( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName())["spec_area"].getDouble() );
				shark_dat->MultiAdsList[i].setTotalVolume(shark_dat->volume);
				shark_dat->MultiAdsList[i].setTotalMass( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName())["total_mass"].getDouble() );
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->MultiAdsList[i].setSurfaceChargeBool( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName())["include_surfcharge"].getBool() );
			}
			catch (std::out_of_range)
			{
				shark_dat->MultiAdsList[i].setSurfaceChargeBool( true );
			}

			//Set the activity function
			int surf_act;
			try
			{
				surf_act = surf_act_choice( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName())["surf_activity"].getString() );
				shark_dat->MultiAdsList[i].setActivityEnum(surf_act);

				switch (surf_act)
				{
					case IDEAL_ADS:
						shark_dat->MultiAdsList[i].setActivityModelInfo(ideal_solution, NULL);
						break;

					case FLORY_HUGGINS:
						shark_dat->MultiAdsList[i].setActivityModelInfo(FloryHuggins_multiligand, &shark_dat->MultiAdsList[i]);
						break;

					case UNIQUAC_ACT:
						shark_dat->MultiAdsList[i].setActivityModelInfo(UNIQUAC_multiligand, &shark_dat->MultiAdsList[i]);
						break;

					default:
						shark_dat->MultiAdsList[i].setActivityModelInfo(ideal_solution, NULL);
						break;
				}
			} catch (std::out_of_range)
			{
				shark_dat->MultiAdsList[i].setActivityModelInfo(ideal_solution, NULL);
				surf_act = IDEAL_ADS;
			}

			//Check to make sure that this Multiligand object has the correct number of ligands
			int num_ligands;
			try
			{
				num_ligands = (int)shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName()).getHeadMap().size();
				if (num_ligands != shark_dat->num_multi_ssar[i].size())
				{
					mError(missing_information);
					return -1;
				}
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			//Now iterate over the Headers to find the ligand names and number of reactions for each ligand
			int ligand = 0;
			for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName()).getHeadMap())
			{
				try
				{
					shark_dat->num_multi_ssar[i][ligand] = shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName())(x.first)["num_rxns"].getInt();
					shark_dat->ssmulti_names[i][ligand] = shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiAdsList[i].getAdsorbentName())(x.first)["name"].getString();
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}
				ligand++;
			}
		}
	}

	return success;
}

//Read in the information about each multisite chemisorbent object (must be done before calling setup)
int read_multichemi_scenario(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_multi_sschem > 0)
	{
		//Loop through ligand names and match Document name with object name
		for (int i=0; i<shark_dat->num_multi_sschem; i++)
		{
			//Set the values of the fiber constants
			try
			{
				shark_dat->MultiChemList[i].setSpecificArea( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName())["spec_area"].getDouble() );
				shark_dat->MultiChemList[i].setTotalVolume(shark_dat->volume);
				shark_dat->MultiChemList[i].setTotalMass( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName())["total_mass"].getDouble() );
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->MultiChemList[i].setSurfaceChargeBool( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName())["include_surfcharge"].getBool() );
			}
			catch (std::out_of_range)
			{
				shark_dat->MultiChemList[i].setSurfaceChargeBool( true );
			}

			//Set the activity function
			int surf_act;
			try
			{
				surf_act = surf_act_choice( shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName())["surf_activity"].getString() );
				shark_dat->MultiChemList[i].setActivityEnum(surf_act);

				switch (surf_act)
				{
					case IDEAL_ADS:
						shark_dat->MultiChemList[i].setActivityModelInfo(ideal_solution, NULL);
						break;

					case FLORY_HUGGINS:
						shark_dat->MultiChemList[i].setActivityModelInfo(FloryHuggins_multichemi, &shark_dat->MultiChemList[i]);
						break;

					case UNIQUAC_ACT:
						shark_dat->MultiChemList[i].setActivityModelInfo(UNIQUAC_multichemi, &shark_dat->MultiChemList[i]);
						break;

					default:
						shark_dat->MultiChemList[i].setActivityModelInfo(ideal_solution, NULL);
						break;
				}
			} catch (std::out_of_range)
			{
				shark_dat->MultiChemList[i].setActivityModelInfo(ideal_solution, NULL);
				surf_act = IDEAL_ADS;
			}

			//Check to make sure that this Multiligand object has the correct number of ligands
			int num_ligands;
			try
			{
				num_ligands = (int)shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName()).getHeadMap().size();
				if (num_ligands != shark_dat->num_multichem_rxns[i].size())
				{
					mError(missing_information);
					return -1;
				}
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			//Now iterate over the Headers to find the ligand names and number of reactions for each ligand
			int ligand = 0;
			for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName()).getHeadMap())
			{
				try
				{
					shark_dat->num_multichem_rxns[i][ligand] = shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName())(x.first)["num_rxns"].getInt();
					shark_dat->ssmultichem_names[i][ligand] = shark_dat->yaml_object.getYamlWrapper()(shark_dat->MultiChemList[i].getAdsorbentName())(x.first)["name"].getString();
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}
				ligand++;
			}
		}
	}

	return success;
}

//Check for user options to set variables
int read_options(SHARK_DATA *shark_dat)
{
	int success = 0;

	try
	{
		shark_dat->Newton_data.LineSearch = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["line_search"].getBool();
	}
	catch (std::out_of_range)
	{
		shark_dat->Newton_data.LineSearch = true;
	}

	try
	{
		shark_dat->LocalMin = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["local_min"].getBool();
	}
	catch (std::out_of_range)
	{
		shark_dat->LocalMin = true;
	}

	if (shark_dat->Newton_data.LineSearch == true)
	{
		try
		{
			shark_dat->Newton_data.Bounce = linesearch_choice(shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["search_type"].getString());
		}
		catch (std::out_of_range)
		{
			shark_dat->Newton_data.Bounce =false;
		}
	}

	try
	{
		shark_dat->Newton_data.NL_Output = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["nl_print"].getBool();
	}
	catch (std::out_of_range)
	{
		if (shark_dat->steadystate == true)
			shark_dat->Newton_data.NL_Output = true;
		else
			shark_dat->Newton_data.NL_Output = false;
	}

	try
	{
		shark_dat->Newton_data.L_Output = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["l_print"].getBool();
	}
	catch (std::out_of_range)
	{
		shark_dat->Newton_data.L_Output = false;
	}

	try
	{
		shark_dat->Newton_data.linear_solver = linearsolve_choice(shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["linear_solve"].getString());
	}
	catch (std::out_of_range)
	{
		if ( shark_dat->lin_precon == NULL && shark_dat->numvar >= 100)
		{
			shark_dat->Newton_data.linear_solver = GMRESRP;
		}
		else
		{
			shark_dat->Newton_data.linear_solver = QR;
		}
	}

	if (shark_dat->Newton_data.linear_solver == GMRESRP || shark_dat->Newton_data.linear_solver == GMRESLP || shark_dat->Newton_data.linear_solver == GCR || shark_dat->Newton_data.linear_solver == GMRESR || shark_dat->Newton_data.linear_solver == KMS)
	{
		int restart;
		try
		{
			restart = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["restart"].getInt();
		}
		catch (std::out_of_range)
		{
			if (shark_dat->numvar <= 100)
				restart = shark_dat->numvar;
			else
				restart = 100;

			if (restart > shark_dat->numvar)
				restart = shark_dat->numvar;
		}

		switch (shark_dat->Newton_data.linear_solver)
		{
			case GMRESRP:
				shark_dat->Newton_data.gmresrp_dat.restart = restart;
				break;

			case GMRESLP:
				shark_dat->Newton_data.gmreslp_dat.restart = restart;
				break;

			case GMRESR:
				shark_dat->Newton_data.gmresr_dat.gcr_restart = restart;
				break;

			case GCR:
				shark_dat->Newton_data.gcr_dat.restart = restart;
				break;

			case KMS:
				shark_dat->Newton_data.kms_dat.restart = restart;
				break;

			default:
				break;
		}
	}

	if (shark_dat->Newton_data.linear_solver == GMRESR)
	{
		try
		{
			shark_dat->Newton_data.gmresr_dat.gmres_tol = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")("gmresr_options")["inner_tol"].getDouble();
		}
		catch (std::out_of_range)
		{
			shark_dat->Newton_data.gmresr_dat.gmres_tol = 0.1;
		}

		try
		{
			shark_dat->Newton_data.gmresr_dat.gmres_restart = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")("gmresr_options")["inner_restart"].getInt();
		}
		catch (std::out_of_range)
		{
			shark_dat->Newton_data.gmresr_dat.gmres_restart = -1;
		}
	}

	if (shark_dat->Newton_data.linear_solver == KMS)
	{
		try
		{
			shark_dat->Newton_data.kms_dat.max_level = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")("kms_options")["max_level"].getInt();
		}
		catch (std::out_of_range)
		{
			shark_dat->Newton_data.kms_dat.max_level = 0;
		}

		try
		{
			shark_dat->Newton_data.kms_dat.inner_reltol = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")("kms_options")["inner_tol"].getDouble();
		}
		catch (std::out_of_range)
		{
			shark_dat->Newton_data.kms_dat.inner_reltol = 0.1;
		}
	}

	try
	{
		shark_dat->Newton_data.nl_maxit = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["nl_maxit"].getInt();
	}
	catch (std::out_of_range)
	{
		if (shark_dat->steadystate == true)
			shark_dat->Newton_data.nl_maxit = 2*shark_dat->numvar;
		else
			shark_dat->Newton_data.nl_maxit = shark_dat->numvar;
	}

	try
	{
		shark_dat->Newton_data.nl_tol_abs = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["nl_abstol"].getDouble();
	}
	catch (std::out_of_range)
	{
		shark_dat->Newton_data.nl_tol_abs = 1e-6;
	}

	try
	{
		shark_dat->Newton_data.nl_tol_rel = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["nl_reltol"].getDouble();
	}
	catch (std::out_of_range)
	{
		shark_dat->Newton_data.nl_tol_rel = 1e-8;
	}

	try
	{
		shark_dat->Newton_data.lin_tol_rel = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["lin_reltol"].getDouble();
	}
	catch (std::out_of_range)
	{
		shark_dat->Newton_data.lin_tol_rel = 1e-6;
	}

	try
	{
		shark_dat->Newton_data.lin_tol_abs = shark_dat->yaml_object.getYamlWrapper()("SolverOptions")["lin_abstol"].getDouble();
	}
	catch (std::out_of_range)
	{
		shark_dat->Newton_data.lin_tol_abs = 1e-6;
	}

	return success;
}

//Check and read the species information
int read_species(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Start reading in the registered species from the reg Header
	int reg_species;
	try
	{
		reg_species = shark_dat->yaml_object.getYamlWrapper()("MasterSpecies")("reg").getDataMap().size();
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	try
	{
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("MasterSpecies")("reg").getDataMap().getMap())
		{
			int index = atoi(x.first.c_str());
			shark_dat->MasterList.set_species(index, x.second.getString());
		}
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	//Now read in the unregistered species and their corresponding information
	if (reg_species < shark_dat->numvar)
	{
		int unreg_species;
		try
		{
			unreg_species = (int) shark_dat->yaml_object.getYamlWrapper()("MasterSpecies")("unreg").getSubMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}

		if ((unreg_species+reg_species) != shark_dat->numvar)
		{
			mError(missing_information);
			return -1;
		}
		try
		{
			std::string formula, phase, name, lin_form;
			double charge, enthalpy, entropy, energy;
			bool haveHS, haveG;

			for (auto &x: shark_dat->yaml_object.getYamlWrapper()("MasterSpecies")("unreg").getSubMap())
			{
				int index = atoi(x.first.c_str());
				formula = x.second.getMap().getString("formula");
				charge = x.second.getMap().getDouble("charge");
				enthalpy = x.second.getMap().getDouble("enthalpy");
				entropy = x.second.getMap().getDouble("entropy");
				haveHS = x.second.getMap().getBool("have_HS");
				energy = x.second.getMap().getDouble("energy");
				haveG = x.second.getMap().getBool("have_G");
				phase = x.second.getMap().getString("phase");
				name = x.second.getMap().getString("name");
				lin_form = x.second.getMap().getString("lin_form");
				shark_dat->MasterList.set_species(index, charge, enthalpy, entropy, energy, haveHS, haveG, phase, name, formula, lin_form);
			}
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
	}

	//Lastly, check the newly registered objects in MasterSpecies list to ensure that no register errors occured
	bool reg_error = false;
	for (int i=0; i<shark_dat->numvar; i++)
	{
		if (shark_dat->MasterList.get_species(i).isRegistered() == false)
		{
			reg_error = true;
			break;
		}
	}
	if (reg_error == true)
	{
		mError(missing_information);
		return -1;
	}

	return success;
}

//Read and check the mass balance document
int read_massbalance(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Check the number of headers in MassBalance Doc and make sure it matches number of equations specified earlier
	int mbes;
	try
	{
		mbes = (int) shark_dat->yaml_object.getYamlWrapper()("MassBalance").getHeadMap().size();
	}
	catch (std::out_of_range)
	{
		mError(missing_information);
		return -1;
	}

	if (mbes != shark_dat->num_mbe)
	{
		mError(missing_information);
		return -1;
	}
	else
	{
		int i=0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("MassBalance").getHeadMap())
		{
			shark_dat->MassBalanceList[i].Set_Name(x.first);
			if (shark_dat->reactor_type == BATCH)
			{
				try
				{
					shark_dat->MassBalanceList[i].Set_TotalConcentration(shark_dat->yaml_object.getYamlWrapper()("MassBalance")(x.first)["total_conc"].getDouble());
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}
			}

			if (shark_dat->reactor_type == CSTR || shark_dat->reactor_type == PFR)
			{
				if (shark_dat->steadystate == false)
				{
					try
					{
						shark_dat->MassBalanceList[i].Set_InitialConcentration(shark_dat->yaml_object.getYamlWrapper()("MassBalance")(x.first)["initial_conc"].getDouble());
					}
					catch (std::out_of_range)
					{
						mError(missing_information);
						return -1;
					}
				}

				try
				{
					shark_dat->MassBalanceList[i].Set_InletConcentration(shark_dat->yaml_object.getYamlWrapper()("MassBalance")(x.first)["inlet_conc"].getDouble());
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}
			}

			int species;
			try
			{
				species = (int) shark_dat->yaml_object.getYamlWrapper()("MassBalance")(x.first)("delta").getMap().size();
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			if (species < 1)
			{
				mError(missing_information);
				return -1;
			}

			for (auto &y: shark_dat->yaml_object.getYamlWrapper()("MassBalance")(x.first)("delta").getMap())
			{
				int index = shark_dat->MasterList.get_index(y.first);
				if (index < 0 || index > (shark_dat->numvar-1))
				{
					std::cout << "\nInvalid Name in List: " << y.first << std::endl;
					mError(read_error);
					return -1;
				}
				else
				{
					try
					{
						shark_dat->MassBalanceList[i].Set_Delta(index, y.second.getDouble());
					}
					catch (std::out_of_range)
					{
						mError(read_error);
						return -1;
					}
				}
			}

			i++;
		}
	}

	return success;
}

//Check and read the equilibiurm reaction document
int read_equilrxn(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Check the number of headers to make sure it matches the expected number of reaction objects
	if (shark_dat->num_ssr > 0)
	{
		int ssr;
		try
		{
			ssr = (int) shark_dat->yaml_object.getYamlWrapper()("EquilRxn").getHeadMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (ssr != shark_dat->num_ssr)
		{
			mError(missing_information);
			return -1;
		}

		int i=0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("EquilRxn").getHeadMap())
		{
			try
			{
				shark_dat->ReactionList[i].Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()("EquilRxn")(x.first)["logK"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			int count = 0;
			double dH, dS;
			try
			{
				dH = shark_dat->yaml_object.getYamlWrapper()("EquilRxn")(x.first)["enthalpy"].getDouble();
				shark_dat->ReactionList[i].Set_Enthalpy(dH);
				count++;
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				dS = shark_dat->yaml_object.getYamlWrapper()("EquilRxn")(x.first)["entropy"].getDouble();
				shark_dat->ReactionList[i].Set_Entropy(dS);
				count++;
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}
			if (count == 2)
				shark_dat->ReactionList[i].Set_EnthalpyANDEntropy(dH, dS);

			try
			{
				shark_dat->ReactionList[i].Set_Energy(shark_dat->yaml_object.getYamlWrapper()("EquilRxn")(x.first)["energy"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			int stoich;
			try
			{
				stoich = shark_dat->yaml_object.getYamlWrapper()("EquilRxn")(x.first)("stoichiometry").getMap().size();
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			if (stoich < 2)
			{
				mError(missing_information);
				return -1;
			}

			for (auto &y: shark_dat->yaml_object.getYamlWrapper()("EquilRxn")(x.first)("stoichiometry").getMap())
			{
				int index = shark_dat->MasterList.get_index(y.first);
				if (index < 0 || index > (shark_dat->numvar-1))
				{
					mError(read_error);
					return -1;
				}
				else
				{
					try
					{
						shark_dat->ReactionList[i].Set_Stoichiometric(index, y.second.getDouble());
					}
					catch (std::out_of_range)
					{
						mError(read_error);
						return -1;
					}
				}
			}

			i++;
		}
	}

	return success;
}

//Check and read the unsteady reaction object
int read_unsteadyrxn(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Check the number of headers to make sure it matches the expected number of reaction objects
	if (shark_dat->num_usr > 0)
	{
		int usr;
		try
		{
			usr = (int) shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn").getHeadMap().size();
		}
		catch (std::out_of_range)
		{
			mError(missing_information);
			return -1;
		}
		if (usr != shark_dat->num_usr)
		{
			mError(missing_information);
			return -1;
		}

		int i=0;
		for (auto &x: shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn").getHeadMap())
		{
			int var_index;
			try
			{
				var_index = shark_dat->MasterList.get_index(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["unsteady_var"].getString());
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			if (var_index < 0 || var_index > (shark_dat->numvar-1))
			{
				mError(read_error);
				return -1;
			}
			shark_dat->UnsteadyList[i].Set_Species_Index(var_index);

			try
			{
				shark_dat->UnsteadyList[i].Set_InitialValue(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["initial_condition"].getDouble());
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["logK"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_Forward(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["forward"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_Reverse(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["reverse"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_ReverseRef(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["reverse_ref"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_ForwardRef(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["forward_ref"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_ActivationEnergy(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["activation_energy"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				shark_dat->UnsteadyList[i].Set_Affinity(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["temp_affinity"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			int count = 0;
			double dH, dS;
			try
			{
				dH = shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["enthalpy"].getDouble();
				shark_dat->UnsteadyList[i].Set_Enthalpy(dH);
				count++;
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			try
			{
				dS = shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["entropy"].getDouble();
				shark_dat->UnsteadyList[i].Set_Entropy(dS);
				count++;
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}
			if (count == 2)
				shark_dat->UnsteadyList[i].Set_EnthalpyANDEntropy(dH, dS);

			try
			{
				shark_dat->UnsteadyList[i].Set_Energy(shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)["energy"].getDouble());
			}
			catch (std::out_of_range)
			{
				//At this point, it is unknown as to whether or not this is an actual error
				//It will be checked later whether or not this causes a problem
			}

			int stoich;
			try
			{
				stoich = shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)("stoichiometry").getMap().size();
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			if (stoich < 2)
			{
				mError(missing_information);
				return -1;
			}

			for (auto &y: shark_dat->yaml_object.getYamlWrapper()("UnsteadyRxn")(x.first)("stoichiometry").getMap())
			{
				int index = shark_dat->MasterList.get_index(y.first);
				if (index < 0 || index > (shark_dat->numvar-1))
				{
					mError(read_error);
					return -1;
				}
				else
				{
					try
					{
						shark_dat->UnsteadyList[i].Set_Stoichiometric(index, y.second.getDouble());
					}
					catch (std::out_of_range)
					{
						mError(read_error);
						return -1;
					}
				}
			}

			i++;
		}
	}


	return success;
}

//Read the adsorption objects
int read_adsorbobjects(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_ssao > 0)
	{
		for (int i=0; i<shark_dat->num_ssao; i++)
		{
			//Check for existance of the necessary object and quit if necessary
			try
			{
				shark_dat->AdsorptionList[i].setAdsorbentName( shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getName() );
				shark_dat->AdsorptionList[i].setTotalVolume(shark_dat->volume);
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			// Other required pieces of information
			try
			{
				shark_dat->AdsorptionList[i].setBasis(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getString("basis"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->AdsorptionList[i].setTotalMass(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getDouble("total_mass"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->AdsorptionList[i].setSpecificArea(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getDouble("spec_area"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			// Some optional pieces of information
			try
			{
				shark_dat->AdsorptionList[i].setSpecificMolality(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getDouble("spec_mole"));
			}
			catch (std::out_of_range)
			{
				shark_dat->AdsorptionList[i].setSpecificMolality(1.0);
			}
			try
			{
				shark_dat->AdsorptionList[i].setSurfaceCharge(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getDouble("surf_charge"));
			}
			catch (std::out_of_range)
			{
				shark_dat->AdsorptionList[i].setSurfaceCharge(0.0);
			}
			try
			{
				shark_dat->AdsorptionList[i].setSurfaceChargeBool(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getBool("include_surfcharge"));
			}
			catch (std::out_of_range)
			{
				shark_dat->AdsorptionList[i].setSurfaceChargeBool(true);
			}
			int surf_act;
			try
			{
				surf_act = surf_act_choice(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getDataMap().getString("surf_activity"));
				shark_dat->AdsorptionList[i].setActivityEnum(surf_act);

				switch (surf_act)
				{
					case IDEAL_ADS:
						shark_dat->AdsorptionList[i].setActivityModelInfo(ideal_solution, NULL);
						break;

					case FLORY_HUGGINS:
						shark_dat->AdsorptionList[i].setActivityModelInfo(FloryHuggins, &shark_dat->AdsorptionList[i]);
						break;

					case UNIQUAC_ACT:
						shark_dat->AdsorptionList[i].setActivityModelInfo(UNIQUAC, &shark_dat->AdsorptionList[i]);
						break;

					default:
						shark_dat->AdsorptionList[i].setActivityModelInfo(ideal_solution, NULL);
						break;
				}
			} catch (std::out_of_range)
			{
				shark_dat->AdsorptionList[i].setActivityModelInfo(ideal_solution, NULL);
				surf_act = IDEAL_ADS;
			}

			//Read in all reaction information
			bool ContainsVolumeFactors;
			bool ContainsAreaFactors;
			std::string vol_check;
			std::string area_check;
			try
			{
				vol_check = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])("volume_factors").getName();
				ContainsVolumeFactors = true;
			}
			catch (std::out_of_range)
			{
				ContainsVolumeFactors = false;
			}
			try
			{
				area_check = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])("area_factors").getName();
				ContainsAreaFactors = true;
			}
			catch (std::out_of_range)
			{
				ContainsAreaFactors = false;
			}
			int num_head;
			try
			{
				num_head = (int)shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getHeadMap().size();
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			if (ContainsVolumeFactors == true && ContainsAreaFactors == false)
			{
				if (num_head != shark_dat->num_ssar[i]+1)
				{
					mError(missing_information);
					return -1;
				}
			}
			else if (ContainsVolumeFactors == true && ContainsAreaFactors == true)
			{
				if (num_head != shark_dat->num_ssar[i]+2)
				{
					mError(missing_information);
					return -1;
				}
			}
			else if (ContainsVolumeFactors == false && ContainsAreaFactors == true)
			{
				if (num_head != shark_dat->num_ssar[i]+1)
				{
					mError(missing_information);
					return -1;
				}
			}
			else
			{
				if (num_head != shark_dat->num_ssar[i])
				{
					mError(missing_information);
					return -1;
				}
			}

			//Loop over all headers
			int rxn = 0;
			for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i]).getHeadMap())
			{
				if (x.second.getName() != "volume_factors" && x.second.getName() != "area_factors")
				{
					if (shark_dat->AdsorptionList[i].isAreaBasis() == false)
					{
						try
						{
							shark_dat->AdsorptionList[i].setMolarFactor(rxn, shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)["mole_factor"].getDouble());
						}
						catch (std::out_of_range)
						{
							mError(missing_information);
							return -1;
						}
					}

					try
					{
						shark_dat->AdsorptionList[i].getReaction(rxn).Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)["logK"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					int count = 0;
					double dH, dS;
					try
					{
						dH = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)["enthalpy"].getDouble();
						shark_dat->AdsorptionList[i].getReaction(rxn).Set_Enthalpy(dH);
						count++;
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						dS = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)["entropy"].getDouble();
						shark_dat->AdsorptionList[i].getReaction(rxn).Set_Entropy(dS);
						count++;
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}
					if (count == 2)
						shark_dat->AdsorptionList[i].getReaction(rxn).Set_EnthalpyANDEntropy(dH, dS);

					try
					{
						shark_dat->AdsorptionList[i].getReaction(rxn).Set_Energy(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)["energy"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					int stoich;
					try
					{
						stoich = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)("stoichiometry").getMap().size();
					}
					catch (std::out_of_range)
					{
						mError(missing_information);
						return -1;
					}
					if (stoich < 2)
					{
						mError(missing_information);
						return -1;
					}

					for (auto &y: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])(x.first)("stoichiometry").getMap())
					{
						int index = shark_dat->MasterList.get_index(y.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->AdsorptionList[i].getReaction(rxn).Set_Stoichiometric(index, y.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}
					}

					rxn++;
				}
				else
				{
					//No Action
				}
			}

			// Read volume factors and check for errors
			if (surf_act != IDEAL_ADS || shark_dat->AdsorptionList[i].isAreaBasis() == true)
			{
				int num_fact = 0;
				bool HaveVol = false;

				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])("volume_factors").getDataMap().size();
					HaveVol = true;

					if (num_fact != shark_dat->num_ssar[i])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall volume_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])("volume_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->AdsorptionList[i].setVolumeFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					HaveVol = false;
					//Loop to set volumes based on mola object
					for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
					{
						if (shark_dat->MasterList.get_species(index).MolarVolume() <= 0.0)
							shark_dat->MasterList.get_species(index).setMolarVolume(7.24);
						shark_dat->AdsorptionList[i].setVolumeFactor(index, shark_dat->MasterList.get_species(index).MolarVolume()*0.602);
					}
				}

				num_fact = 0;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])("area_factors").getDataMap().size();

					if (num_fact != shark_dat->num_ssar[i])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall area_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_ads_names[i])("area_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->AdsorptionList[i].setAreaFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					if (HaveVol == true)
						shark_dat->AdsorptionList[i].calculateAreaFactors();
					else
					{
						//Loop to set area factors based on mola object
						for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
						{
							if (shark_dat->MasterList.get_species(index).MolarArea() <= 0.0)
								shark_dat->MasterList.get_species(index).setMolarArea(18.10);
							shark_dat->AdsorptionList[i].setAreaFactor(index, shark_dat->MasterList.get_species(index).MolarArea()*6020.0);
						}
					}
				}

			}

		}
	}

	return success;
}

/// Function to go through the yaml object for each Unsteady Adsorption Object
int read_unsteadyadsorbobjects(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_usao > 0)
	{
		for (int i=0; i<shark_dat->num_usao; i++)
		{
			//Check for existance of the necessary object and quit if necessary
			try
			{
				shark_dat->UnsteadyAdsList[i].setAdsorbentName( shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getName() );
				shark_dat->UnsteadyAdsList[i].setTotalVolume(shark_dat->volume);
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			// Other required pieces of information
			try
			{
				shark_dat->UnsteadyAdsList[i].setBasis(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getString("basis"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->UnsteadyAdsList[i].setTotalMass(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getDouble("total_mass"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->UnsteadyAdsList[i].setSpecificArea(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getDouble("spec_area"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			// Some optional pieces of information
			try
			{
				shark_dat->UnsteadyAdsList[i].setSpecificMolality(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getDouble("spec_mole"));
			}
			catch (std::out_of_range)
			{
				shark_dat->UnsteadyAdsList[i].setSpecificMolality(1.0);
			}
			try
			{
				shark_dat->UnsteadyAdsList[i].setSurfaceCharge(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getDouble("surf_charge"));
			}
			catch (std::out_of_range)
			{
				shark_dat->UnsteadyAdsList[i].setSurfaceCharge(0.0);
			}
			try
			{
				shark_dat->UnsteadyAdsList[i].setSurfaceChargeBool(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getBool("include_surfcharge"));
			}
			catch (std::out_of_range)
			{
				shark_dat->UnsteadyAdsList[i].setSurfaceChargeBool(true);
			}
			int surf_act;
			try
			{
				surf_act = surf_act_choice(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getDataMap().getString("surf_activity"));
				shark_dat->UnsteadyAdsList[i].setActivityEnum(surf_act);

				switch (surf_act)
				{
					case IDEAL_ADS:
						shark_dat->UnsteadyAdsList[i].setActivityModelInfo(ideal_solution, NULL);
						break;

					case FLORY_HUGGINS:
						shark_dat->UnsteadyAdsList[i].setActivityModelInfo(FloryHuggins_unsteady, &shark_dat->UnsteadyAdsList[i]);
						break;

					case UNIQUAC_ACT:
						shark_dat->UnsteadyAdsList[i].setActivityModelInfo(UNIQUAC_unsteady, &shark_dat->UnsteadyAdsList[i]);
						break;

					default:
						shark_dat->UnsteadyAdsList[i].setActivityModelInfo(ideal_solution, NULL);
						break;
				}
			} catch (std::out_of_range)
			{
				shark_dat->UnsteadyAdsList[i].setActivityModelInfo(ideal_solution, NULL);
				surf_act = IDEAL_ADS;
			}

			//Read in all reaction information
			bool ContainsVolumeFactors;
			bool ContainsAreaFactors;
			std::string vol_check;
			std::string area_check;
			try
			{
				vol_check = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])("volume_factors").getName();
				ContainsVolumeFactors = true;
			}
			catch (std::out_of_range)
			{
				ContainsVolumeFactors = false;
			}
			try
			{
				area_check = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])("area_factors").getName();
				ContainsAreaFactors = true;
			}
			catch (std::out_of_range)
			{
				ContainsAreaFactors = false;
			}
			int num_head;
			try
			{
				num_head = (int)shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getHeadMap().size();
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			if (ContainsVolumeFactors == true && ContainsAreaFactors == false)
			{
				if (num_head != shark_dat->num_usar[i]+1)
				{
					mError(missing_information);
					return -1;
				}
			}
			else if (ContainsVolumeFactors == true && ContainsAreaFactors == true)
			{
				if (num_head != shark_dat->num_usar[i]+2)
				{
					mError(missing_information);
					return -1;
				}
			}
			else if (ContainsVolumeFactors == false && ContainsAreaFactors == true)
			{
				if (num_head != shark_dat->num_usar[i]+1)
				{
					mError(missing_information);
					return -1;
				}
			}
			else
			{
				if (num_head != shark_dat->num_usar[i])
				{
					mError(missing_information);
					return -1;
				}
			}

			//Loop over all headers
			int rxn = 0;
			for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i]).getHeadMap())
			{
				if (x.second.getName() != "volume_factors" && x.second.getName() != "area_factors")
				{
					if (shark_dat->UnsteadyAdsList[i].isAreaBasis() == false)
					{
						try
						{
							shark_dat->UnsteadyAdsList[i].setMolarFactor(rxn, shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["mole_factor"].getDouble());
						}
						catch (std::out_of_range)
						{
							mError(missing_information);
							return -1;
						}
					}

					int var_index;
					try
					{
						var_index = shark_dat->MasterList.get_index(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["unsteady_var"].getString());
					}
					catch (std::out_of_range)
					{
						mError(missing_information);
						return -1;
					}
					if (var_index < 0 || var_index > (shark_dat->numvar-1))
					{
						mError(read_error);
						return -1;
					}
					shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Species_Index(var_index);

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_InitialValue(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["initial_condition"].getDouble());
					}
					catch (std::out_of_range)
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_InitialValue(0.0);
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["logK"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Forward(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["forward"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Reverse(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["reverse"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_ReverseRef(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["reverse_ref"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_ForwardRef(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["forward_ref"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_ActivationEnergy(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["activation_energy"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Affinity(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["temp_affinity"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					int count = 0;
					double dH, dS;
					try
					{
						dH = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["enthalpy"].getDouble();
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Enthalpy(dH);
						count++;
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						dS = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["entropy"].getDouble();
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Entropy(dS);
						count++;
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}
					if (count == 2)
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_EnthalpyANDEntropy(dH, dS);

					try
					{
						shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Energy(shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)["energy"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					int stoich;
					try
					{
						stoich = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)("stoichiometry").getMap().size();
					}
					catch (std::out_of_range)
					{
						mError(missing_information);
						return -1;
					}
					if (stoich < 2)
					{
						mError(missing_information);
						return -1;
					}

					for (auto &y: shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])(x.first)("stoichiometry").getMap())
					{
						int index = shark_dat->MasterList.get_index(y.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->UnsteadyAdsList[i].getReaction(rxn).Set_Stoichiometric(index, y.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}
					}

					rxn++;
				}
				else
				{
					//No Action
				}
			}

			// Read volume factors and check for errors
			if (surf_act != IDEAL_ADS || shark_dat->UnsteadyAdsList[i].isAreaBasis() == true)
			{
				int num_fact = 0;
				bool HaveVol = false;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])("volume_factors").getDataMap().size();
					HaveVol = true;

					if (num_fact != shark_dat->num_usar[i])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall volume_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])("volume_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->UnsteadyAdsList[i].setVolumeFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					HaveVol = false;
					//Loop to set volumes based on mola object
					for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
					{
						if (shark_dat->MasterList.get_species(index).MolarVolume() <= 0.0)
							shark_dat->MasterList.get_species(index).setMolarVolume(7.24);
						shark_dat->UnsteadyAdsList[i].setVolumeFactor(index, shark_dat->MasterList.get_species(index).MolarVolume()*0.602);
					}
				}

				num_fact = 0;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])("area_factors").getDataMap().size();

					if (num_fact != shark_dat->num_usar[i])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall area_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->us_ads_names[i])("area_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->UnsteadyAdsList[i].setAreaFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					if (HaveVol == true)
						shark_dat->UnsteadyAdsList[i].calculateAreaFactors();
					else
					{
						//Loop to set area factors based on mola object
						for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
						{
							if (shark_dat->MasterList.get_species(index).MolarArea() <= 0.0)
								shark_dat->MasterList.get_species(index).setMolarArea(18.10);
							shark_dat->UnsteadyAdsList[i].setAreaFactor(index, shark_dat->MasterList.get_species(index).MolarArea()*6020.0);
						}
					}
				}

			}

		}
	}

	return success;
}

//Function to read in multiligand object information from input file
int read_multiligandobjects(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_multi_ssao > 0)
	{
		//Loop over all multiligand adsorbent objects
		for (int i=0; i<shark_dat->MultiAdsList.size(); i++)
		{
			//Loop over all ligands for the give object
			for (int l=0; l<shark_dat->MultiAdsList[i].getNumberLigands(); l++)
			{
				std::string docname = shark_dat->MultiAdsList[i].getAdsorbentName() + "-" + shark_dat->ssmulti_names[i][l];
				shark_dat->MultiAdsList[i].setLigandName(l, shark_dat->ssmulti_names[i][l]);
				try
				{
					shark_dat->MultiAdsList[i].setSpecificMolality(l, shark_dat->yaml_object.getYamlWrapper()(docname)["spec_mole"].getDouble() );
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}
				try
				{
					shark_dat->MultiAdsList[i].setSurfaceCharge(l, shark_dat->yaml_object.getYamlWrapper()(docname)["ligand_charge"].getDouble() );
				}
				catch (std::out_of_range)
				{
					shark_dat->MultiAdsList[i].setSurfaceCharge(l, 0.0);
				}

				//Check to see what kind of information is given is as headers
				bool ContainsVolumeFactors;
				bool ContainsAreaFactors;
				std::string vol_check;
				std::string area_check;
				try
				{
					vol_check = shark_dat->yaml_object.getYamlWrapper()(docname)("volume_factors").getName();
					ContainsVolumeFactors = true;
				}
				catch (std::out_of_range)
				{
					ContainsVolumeFactors = false;
				}
				try
				{
					area_check = shark_dat->yaml_object.getYamlWrapper()(docname)("area_factors").getName();
					ContainsAreaFactors = true;
				}
				catch (std::out_of_range)
				{
					ContainsAreaFactors = false;
				}
				int num_header;
				try
				{
					num_header = (int)shark_dat->yaml_object.getYamlWrapper()(docname).getHeadMap().size();
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}

				//Check for missing information
				if (ContainsVolumeFactors == true && ContainsAreaFactors == false)
				{
					if (num_header != shark_dat->num_multi_ssar[i][l]+1)
					{
						mError(missing_information);
						return -1;
					}
				}
				else if (ContainsVolumeFactors == true && ContainsAreaFactors == true)
				{
					if (num_header != shark_dat->num_multi_ssar[i][l]+2)
					{
						mError(missing_information);
						return -1;
					}
				}
				else if (ContainsVolumeFactors == false && ContainsAreaFactors == true)
				{
					if (num_header != shark_dat->num_multi_ssar[i][l]+1)
					{
						mError(missing_information);
						return -1;
					}
				}
				else
				{
					if (num_header != shark_dat->num_multi_ssar[i][l])
					{
						mError(missing_information);
						return -1;
					}
				}

				//Iterate over the headers to input the reaction information
				int rxn = 0;
				for (auto &x: shark_dat->yaml_object.getYamlWrapper()(docname).getHeadMap())
				{
					if (x.second.getName() != "volume_factors" && x.second.getName() != "area_factors")
					{
						try
						{
							shark_dat->MultiAdsList[i].setMolarFactor(l, rxn, shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["mole_factor"].getDouble());
						}
						catch (std::out_of_range)
						{
							mError(missing_information);
							return -1;
						}

						try
						{
							shark_dat->MultiAdsList[i].getAdsorptionObject(l).getReaction(rxn).Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["logK"].getDouble());
						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}

						int count = 0;
						double dH, dS;
						try
						{
							dH = shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["enthalpy"].getDouble();
							shark_dat->MultiAdsList[i].getAdsorptionObject(l).getReaction(rxn).Set_Enthalpy(dH);
							count++;
						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}

						try
						{
							dS = shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["entropy"].getDouble();
							shark_dat->MultiAdsList[i].getAdsorptionObject(l).getReaction(rxn).Set_Entropy(dS);
							count++;
						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}
						if (count == 2)
							shark_dat->MultiAdsList[i].getAdsorptionObject(l).getReaction(rxn).Set_EnthalpyANDEntropy(dH, dS);

						try
						{
							shark_dat->MultiAdsList[i].getAdsorptionObject(l).getReaction(rxn).Set_Energy(shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["energy"].getDouble());

						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}

						int stoich;
						try
						{
							stoich = shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)("stoichiometry").getMap().size();
						}
						catch (std::out_of_range)
						{
							mError(missing_information);
							return -1;
						}
						if (stoich < 2)
						{
							mError(missing_information);
							return -1;
						}

						for (auto &y: shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)("stoichiometry").getMap())
						{
							int index = shark_dat->MasterList.get_index(y.first);
							if (index < 0 || index > (shark_dat->numvar-1))
							{
								mError(read_error);
								return -1;
							}
							else
							{
								try
								{
									shark_dat->MultiAdsList[i].getAdsorptionObject(l).getReaction(rxn).Set_Stoichiometric(index, y.second.getDouble());
								}
								catch (std::out_of_range)
								{
									mError(read_error);
									return -1;
								}
							}
						}

						rxn++;
					}
					else
					{
						//No Action
					}
				}

				//Try to read in any volume and area factors given
				int num_fact = 0;
				bool HaveVol = false;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(docname)("volume_factors").getDataMap().size();
					HaveVol = true;

					if (num_fact != shark_dat->num_multi_ssar[i][l])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall volume_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(docname)("volume_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->MultiAdsList[i].setVolumeFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					HaveVol = false;
					//Loop to set volumes based on mola object
					for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
					{
						if (shark_dat->MasterList.get_species(index).MolarVolume() <= 0.0)
							shark_dat->MasterList.get_species(index).setMolarVolume(7.24);
						shark_dat->MultiAdsList[i].setVolumeFactor(index, shark_dat->MasterList.get_species(index).MolarVolume()*0.602);
					}
				}

				num_fact = 0;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(docname)("area_factors").getDataMap().size();

					if (num_fact != shark_dat->num_multi_ssar[i][l])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall area_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(docname)("area_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->MultiAdsList[i].setAreaFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					if (HaveVol == true)
						shark_dat->MultiAdsList[i].calculateAreaFactors();
					else
					{
						//Loop to set area factors based on mola object
						for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
						{
							if (shark_dat->MasterList.get_species(index).MolarArea() <= 0.0)
								shark_dat->MasterList.get_species(index).setMolarArea(18.10);
							shark_dat->MultiAdsList[i].setAreaFactor(index, shark_dat->MasterList.get_species(index).MolarArea()*6020.0);
						}
					}
				}


			}//END ligand loop

		}//END object loop
	}

	return success;
}

//Read in chemisorption document
int read_chemisorbobjects(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_sschem > 0)
	{
		for (int i=0; i<shark_dat->num_sschem; i++)
		{
			//Check for existance of the necessary object and quit if necessary
			try
			{
				shark_dat->ChemisorptionList[i].setAdsorbentName( shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getName() );
				shark_dat->ChemisorptionList[i].setTotalVolume(shark_dat->volume);
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			// Other required pieces of information
			try
			{
				shark_dat->ChemisorptionList[i].setTotalMass(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getDataMap().getDouble("total_mass"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->ChemisorptionList[i].setSpecificArea(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getDataMap().getDouble("spec_area"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}
			try
			{
				shark_dat->ChemisorptionList[i].setSpecificMolality(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getDataMap().getDouble("spec_mole"));
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			//Set some optional information
			try
			{
				shark_dat->ChemisorptionList[i].setSurfaceChargeBool(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getDataMap().getBool("include_surfcharge"));
			}
			catch (std::out_of_range)
			{
				shark_dat->ChemisorptionList[i].setSurfaceChargeBool(true);
			}
			int surf_act;
			try
			{
				surf_act=surf_act_choice(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getDataMap().getString("surf_activity"));
				shark_dat->ChemisorptionList[i].setActivityEnum(surf_act);

				switch (surf_act)
				{
					case IDEAL_ADS:
						shark_dat->ChemisorptionList[i].setActivityModelInfo(ideal_solution, NULL);
						break;

					case FLORY_HUGGINS:
						shark_dat->ChemisorptionList[i].setActivityModelInfo(FloryHuggins_chemi, &shark_dat->ChemisorptionList[i]);
						break;

					case UNIQUAC_ACT:
						shark_dat->ChemisorptionList[i].setActivityModelInfo(UNIQUAC_chemi, &shark_dat->ChemisorptionList[i]);
						break;

					default:
						shark_dat->ChemisorptionList[i].setActivityModelInfo(ideal_solution, NULL);
						break;
				}
			} catch (std::out_of_range)
			{
				shark_dat->ChemisorptionList[i].setActivityModelInfo(ideal_solution, NULL);
				surf_act = IDEAL_ADS;
			}

			//Read in all reaction information
			bool ContainsVolumeFactors;
			bool ContainsAreaFactors;
			std::string vol_check;
			std::string area_check;
			try
			{
				vol_check = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])("volume_factors").getName();
				ContainsVolumeFactors = true;
			}
			catch (std::out_of_range)
			{
				ContainsVolumeFactors = false;
			}
			try
			{
				area_check = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])("area_factors").getName();
				ContainsAreaFactors = true;
			}
			catch (std::out_of_range)
			{
				ContainsAreaFactors = false;
			}
			int num_head;
			try
			{
				num_head = (int)shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getHeadMap().size();
			}
			catch (std::out_of_range)
			{
				mError(missing_information);
				return -1;
			}

			if (ContainsVolumeFactors == true && ContainsAreaFactors == false)
			{
				if (num_head != shark_dat->num_sschem_rxns[i]+1)
				{
					mError(missing_information);
					return -1;
				}
			}
			else if (ContainsVolumeFactors == true && ContainsAreaFactors == true)
			{
				if (num_head != shark_dat->num_sschem_rxns[i]+2)
				{
					mError(missing_information);
					return -1;
				}
			}
			else if (ContainsVolumeFactors == false && ContainsAreaFactors == true)
			{
				if (num_head != shark_dat->num_sschem_rxns[i]+1)
				{
					mError(missing_information);
					return -1;
				}
			}
			else
			{
				if (num_head != shark_dat->num_sschem_rxns[i])
				{
					mError(missing_information);
					return -1;
				}
			}

			//Loop over all headers
			int rxn = 0;
			for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i]).getHeadMap())
			{
				if (x.second.getName() != "volume_factors" && x.second.getName() != "area_factors")
				{

					try
					{
						shark_dat->ChemisorptionList[i].getReaction(rxn).Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])(x.first)["logK"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					int count = 0;
					double dH, dS;
					try
					{
						dH = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])(x.first)["enthalpy"].getDouble();
						shark_dat->ChemisorptionList[i].getReaction(rxn).Set_Enthalpy(dH);
						count++;
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					try
					{
						dS = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])(x.first)["entropy"].getDouble();
						shark_dat->ChemisorptionList[i].getReaction(rxn).Set_Entropy(dS);
						count++;
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}
					if (count == 2)
						shark_dat->ChemisorptionList[i].getReaction(rxn).Set_EnthalpyANDEntropy(dH, dS);

					try
					{
						shark_dat->ChemisorptionList[i].getReaction(rxn).Set_Energy(shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])(x.first)["energy"].getDouble());
					}
					catch (std::out_of_range)
					{
						//At this point, it is unknown as to whether or not this is an actual error
						//It will be checked later whether or not this causes a problem
					}

					int stoich;
					try
					{
						stoich = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])(x.first)("stoichiometry").getMap().size();
					}
					catch (std::out_of_range)
					{
						mError(missing_information);
						return -1;
					}
					if (stoich < 2)
					{
						mError(missing_information);
						return -1;
					}

					for (auto &y: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])(x.first)("stoichiometry").getMap())
					{
						int index = shark_dat->MasterList.get_index(y.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->ChemisorptionList[i].getReaction(rxn).Set_Stoichiometric(index, y.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}
					}

					rxn++;
				}
				else
				{
					//No Action
				}
			}

			//Set the Ligand Index Here
			success = shark_dat->ChemisorptionList[i].setLigandIndex();
			if (success != 0) {mError(missing_information); return -1;}

			//Set the adsorb indices here
			success = shark_dat->ChemisorptionList[i].setAdsorbIndices();
			if (success != 0) {mError(missing_information); return -1;}

			//Set the Site Balance info for the object
			success = shark_dat->ChemisorptionList[i].setDeltas();
			if (success != 0) {mError(missing_information); return -1;}

			// Read volume factors and check for errors
			int num_fact = 0;
			bool HaveVol = false;
			try
			{
				num_fact = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])("volume_factors").getDataMap().size();
				HaveVol = true;

				if (num_fact != shark_dat->num_sschem_rxns[i])
				{
					mError(missing_information);
					return -1;
				}

				//Loop overall volume_factors
				for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])("volume_factors").getDataMap().getMap())
				{
					int index = shark_dat->MasterList.get_index(x.first);
					if (index < 0 || index > (shark_dat->numvar-1))
					{
						mError(read_error);
						return -1;
					}
					else
					{
						try
						{
							shark_dat->ChemisorptionList[i].setVolumeFactor(index, x.second.getDouble());
						}
						catch (std::out_of_range)
						{
							mError(read_error);
							return -1;
						}
					}

				}
			}
			catch (std::out_of_range)
			{
				HaveVol = false;
				//Loop to set volumes based on mola object
				for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
				{
					if (shark_dat->MasterList.get_species(index).MolarVolume() <= 0.0)
						shark_dat->MasterList.get_species(index).setMolarVolume(7.24);
					shark_dat->ChemisorptionList[i].setVolumeFactor(index, shark_dat->MasterList.get_species(index).MolarVolume()*0.602);
				}
			}

			num_fact = 0;
			try
			{
				num_fact = shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])("area_factors").getDataMap().size();

				if (num_fact != shark_dat->num_sschem_rxns[i])
				{
					mError(missing_information);
					return -1;
				}

				//Loop overall area_factors
				for (auto &x: shark_dat->yaml_object.getYamlWrapper()(shark_dat->ss_chem_names[i])("area_factors").getDataMap().getMap())
				{
					int index = shark_dat->MasterList.get_index(x.first);
					if (index < 0 || index > (shark_dat->numvar-1))
					{
						mError(read_error);
						return -1;
					}
					else
					{
						try
						{
							shark_dat->ChemisorptionList[i].setAreaFactor(index, x.second.getDouble());
						}
						catch (std::out_of_range)
						{
							mError(read_error);
							return -1;
						}
					}

				}
			}
			catch (std::out_of_range)
			{
				if (HaveVol == true)
					shark_dat->ChemisorptionList[i].calculateAreaFactors();
				else
				{
					//Loop to set area factors based on mola object
					for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
					{
						if (shark_dat->MasterList.get_species(index).MolarArea() <= 0.0)
							shark_dat->MasterList.get_species(index).setMolarArea(18.10);
						shark_dat->ChemisorptionList[i].setAreaFactor(index, shark_dat->MasterList.get_species(index).MolarArea()*6020.0);
					}
				}
			}
		}
	}

	return success;
}

//Function to read in multiligand chemisorption object information from input file
int read_multichemiobjects(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->num_multi_sschem > 0)
	{
		//Loop over all multiligand adsorbent objects
		for (int i=0; i<shark_dat->MultiChemList.size(); i++)
		{
			//Loop over all ligands for the give object
			for (int l=0; l<shark_dat->MultiChemList[i].getNumberLigands(); l++)
			{
				std::string docname = shark_dat->MultiChemList[i].getAdsorbentName() + "-" + shark_dat->ssmultichem_names[i][l];
				shark_dat->MultiChemList[i].setLigandName(l, shark_dat->ssmultichem_names[i][l]);

				try
				{
					shark_dat->MultiChemList[i].setSpecificMolality(l, shark_dat->yaml_object.getYamlWrapper()(docname)["spec_mole"].getDouble() );
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}

				//Check to see what kind of information is given is as headers
				bool ContainsVolumeFactors;
				bool ContainsAreaFactors;
				std::string vol_check;
				std::string area_check;
				try
				{
					vol_check = shark_dat->yaml_object.getYamlWrapper()(docname)("volume_factors").getName();
					ContainsVolumeFactors = true;
				}
				catch (std::out_of_range)
				{
					ContainsVolumeFactors = false;
				}
				try
				{
					area_check = shark_dat->yaml_object.getYamlWrapper()(docname)("area_factors").getName();
					ContainsAreaFactors = true;
				}
				catch (std::out_of_range)
				{
					ContainsAreaFactors = false;
				}
				int num_header;
				try
				{
					num_header = (int)shark_dat->yaml_object.getYamlWrapper()(docname).getHeadMap().size();
				}
				catch (std::out_of_range)
				{
					mError(missing_information);
					return -1;
				}

				//Check for missing information
				if (ContainsVolumeFactors == true && ContainsAreaFactors == false)
				{
					if (num_header != shark_dat->num_multichem_rxns[i][l]+1)
					{
						mError(missing_information);
						return -1;
					}
				}
				else if (ContainsVolumeFactors == true && ContainsAreaFactors == true)
				{
					if (num_header != shark_dat->num_multichem_rxns[i][l]+2)
					{
						mError(missing_information);
						return -1;
					}
				}
				else if (ContainsVolumeFactors == false && ContainsAreaFactors == true)
				{
					if (num_header != shark_dat->num_multichem_rxns[i][l]+1)
					{
						mError(missing_information);
						return -1;
					}
				}
				else
				{
					if (num_header != shark_dat->num_multichem_rxns[i][l])
					{
						mError(missing_information);
						return -1;
					}
				}

				//Iterate over the headers to input the reaction information
				int rxn = 0;
				for (auto &x: shark_dat->yaml_object.getYamlWrapper()(docname).getHeadMap())
				{
					if (x.second.getName() != "volume_factors" && x.second.getName() != "area_factors")
					{

						try
						{
							shark_dat->MultiChemList[i].getChemisorptionObject(l).getReaction(rxn).Set_Equilibrium(shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["logK"].getDouble());
						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}

						int count = 0;
						double dH, dS;
						try
						{
							dH = shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["enthalpy"].getDouble();
							shark_dat->MultiChemList[i].getChemisorptionObject(l).getReaction(rxn).Set_Enthalpy(dH);
							count++;
						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}

						try
						{
							dS = shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["entropy"].getDouble();
							shark_dat->MultiChemList[i].getChemisorptionObject(l).getReaction(rxn).Set_Entropy(dS);
							count++;
						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}
						if (count == 2)
							shark_dat->MultiChemList[i].getChemisorptionObject(l).getReaction(rxn).Set_EnthalpyANDEntropy(dH, dS);

						try
						{
							shark_dat->MultiChemList[i].getChemisorptionObject(l).getReaction(rxn).Set_Energy(shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)["energy"].getDouble());

						}
						catch (std::out_of_range)
						{
							//At this point, it is unknown as to whether or not this is an actual error
							//It will be checked later whether or not this causes a problem
						}

						int stoich;
						try
						{
							stoich = shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)("stoichiometry").getMap().size();
						}
						catch (std::out_of_range)
						{
							mError(missing_information);
							return -1;
						}
						if (stoich < 2)
						{
							mError(missing_information);
							return -1;
						}

						for (auto &y: shark_dat->yaml_object.getYamlWrapper()(docname)(x.first)("stoichiometry").getMap())
						{
							int index = shark_dat->MasterList.get_index(y.first);
							if (index < 0 || index > (shark_dat->numvar-1))
							{
								mError(read_error);
								return -1;
							}
							else
							{
								try
								{
									shark_dat->MultiChemList[i].getChemisorptionObject(l).getReaction(rxn).Set_Stoichiometric(index, y.second.getDouble());
								}
								catch (std::out_of_range)
								{
									mError(read_error);
									return -1;
								}
							}
						}

						rxn++;
					}
					else
					{
						//No Action
					}
				}

				//Set the Ligand Index Here
				success = shark_dat->MultiChemList[i].getChemisorptionObject(l).setLigandIndex();
				if (success != 0) {mError(missing_information); return -1;}

				//Set the adsorb indices here
				success = shark_dat->MultiChemList[i].getChemisorptionObject(l).setAdsorbIndices();
				if (success != 0) {mError(missing_information); return -1;}

				//Set the Site Balance info for the object
				success = shark_dat->MultiChemList[i].getChemisorptionObject(l).setDeltas();
				if (success != 0) {mError(missing_information); return -1;}

				//Try to read in any volume and area factors given
				int num_fact = 0;
				bool HaveVol = false;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(docname)("volume_factors").getDataMap().size();
					HaveVol = true;

					if (num_fact != shark_dat->num_multichem_rxns[i][l])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall volume_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(docname)("volume_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->MultiChemList[i].getChemisorptionObject(l).setVolumeFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					HaveVol = false;
					//Loop to set volumes based on mola object
					for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
					{
						if (shark_dat->MasterList.get_species(index).MolarVolume() <= 0.0)
							shark_dat->MasterList.get_species(index).setMolarVolume(7.24);
						shark_dat->MultiChemList[i].getChemisorptionObject(l).setVolumeFactor(index, shark_dat->MasterList.get_species(index).MolarVolume()*0.602);
					}
				}

				num_fact = 0;
				try
				{
					num_fact = shark_dat->yaml_object.getYamlWrapper()(docname)("area_factors").getDataMap().size();

					if (num_fact != shark_dat->num_multichem_rxns[i][l])
					{
						mError(missing_information);
						return -1;
					}

					//Loop overall area_factors
					for (auto &x: shark_dat->yaml_object.getYamlWrapper()(docname)("area_factors").getDataMap().getMap())
					{
						int index = shark_dat->MasterList.get_index(x.first);
						if (index < 0 || index > (shark_dat->numvar-1))
						{
							mError(read_error);
							return -1;
						}
						else
						{
							try
							{
								shark_dat->MultiChemList[i].getChemisorptionObject(l).setAreaFactor(index, x.second.getDouble());
							}
							catch (std::out_of_range)
							{
								mError(read_error);
								return -1;
							}
						}

					}
				}
				catch (std::out_of_range)
				{
					if (HaveVol == true)
						shark_dat->MultiChemList[i].getChemisorptionObject(l).calculateAreaFactors();
					else
					{
						//Loop to set area factors based on mola object
						for (int index = 0; index<shark_dat->MasterList.list_size(); index++)
						{
							if (shark_dat->MasterList.get_species(index).MolarArea() <= 0.0)
								shark_dat->MasterList.get_species(index).setMolarArea(18.10);
							shark_dat->MultiChemList[i].getChemisorptionObject(l).setAreaFactor(index, shark_dat->MasterList.get_species(index).MolarArea()*6020.0);
						}
					}
				}


			}//END ligand loop

		}//END object loop
	}

	return success;
}

//Setup function for a SHARK object
int setup_SHARK_DATA( FILE *file, int (*residual) (const Matrix<double> &x, Matrix<double> &res, const void *data),
					 int (*activity) (const Matrix<double> &x, Matrix<double> &gama, const void *data),
					 int (*precond) (const Matrix<double> &r, Matrix<double> &p, const void *data),
					 SHARK_DATA *dat, const void *activity_data, const void *residual_data,
					 const void *precon_data, const void *other_data)
{
	int success = 0;

	//Check input args
	if ( file == NULL)
		dat->File_Output = false;
	else
	{
		dat->File_Output = true;
		dat->OutputFile = file;
	}
	if ( (*activity) == NULL)
		dat->EvalActivity = ideal_solution;
	else
		dat->EvalActivity = (*activity);
	if ( (*residual) == NULL)
		dat->Residual = shark_residual;
	else
		dat->Residual = (*residual);
	dat->lin_precon = (*precond);
	if (dat->numvar < 2)
	{
		mError(matrix_too_small);
		std::cout << "Problem size is too small for this object!\n";
		return -1;
	}
	int ssao_sum = 0;
	for (int i=0; i<dat->num_ssao; i++)
	{
		if (dat->num_ssar[i] == 0)
		{
			mError(missing_information);
			return -1;
		}
		ssao_sum += dat->num_ssar[i];
	}

	int usao_sum = 0;
	for (int i=0; i<dat->num_usao; i++)
	{
		if (dat->num_usar[i] == 0)
		{
			mError(missing_information);
			return -1;
		}
		usao_sum += dat->num_usar[i];
	}

	int multiligand_sum = 0;
	for (int i=0; i<dat->num_multi_ssao; i++)
	{
		if (dat->num_multi_ssar[i].size() == 0)
		{
			mError(missing_information);
			return -1;
		}

		for (int j=0; j<dat->num_multi_ssar[i].size(); j++)
		{
			if (dat->num_multi_ssar[i][j] == 0)
			{
				mError(missing_information);
				return -1;
			}
			multiligand_sum += dat->num_multi_ssar[i][j];
		}
	}

	int chem_sum=0;
	for (int i=0; i<dat->num_sschem; i++)
	{
		if (dat->num_sschem_rxns[i] == 0)
		{
			mError(missing_information);
			return -1;
		}
		chem_sum += dat->num_sschem_rxns[i] + 1;
	}

	int multichem_sum = 0;
	for (int i=0; i<dat->num_multi_sschem; i++)
	{
		if (dat->num_multichem_rxns[i].size() == 0)
		{
			mError(missing_information);
			return -1;
		}

		for (int j=0; j<dat->num_multichem_rxns[i].size(); j++)
		{
			if (dat->num_multichem_rxns[i][j] == 0)
			{
				mError(missing_information);
				return -1;
			}
			multichem_sum += dat->num_multichem_rxns[i][j] + 1;
		}
	}

	if (dat->numvar != (dat->num_mbe+dat->num_ssr+dat->num_usr+dat->num_other+1+ssao_sum+usao_sum+multiligand_sum+chem_sum+multichem_sum))
	{
		mError(dim_mis_match);
		std::cout << "Number of equations and variables do not match!\n";
		return -1;
	}
	if (dat->SpeciationCurve == true)
	{
		dat->steadystate = true;
		dat->const_pH = true;
	}
	if (dat->steadystate == false && dat->simulationtime <= 0.0)
	{
		dat->steadystate = true;
	}
	if (dat->steadystate == true)
	{
		dat->simulationtime = 0.0;
	}
	if (dat->dt <= 0.0 && dat->steadystate == false)
	{
		dat->dt = 0.1;
	}
	if (dat->dt <= sqrt(DBL_EPSILON) && dat->steadystate == false)
	{
		dat->dt = sqrt(DBL_EPSILON);
	}
	if (dat->pH < 1.0)
		dat->pH = 1.0;
	else if (dat->pH > 14.0)
		dat->pH = 14.0;
	if ( (activity_data) == NULL)
		dat->activity_data = dat;
	else
		dat->activity_data = activity_data;
	if ( (residual_data) == NULL)
		dat->residual_data = dat;
	else
		dat->residual_data = residual_data;
	dat->precon_data = precon_data;
	dat->other_data = other_data;
	dat->totalsteps = 0;
	dat->totalcalls = 0;
	dat->timesteps = 0;
	dat->time_old = 0.0;
	dat->time = 0.0;
	dat->t_count = 0.0;
	dat->dt_min = sqrt(DBL_EPSILON);

	//Edit the activity function if necessary
	if (dat->act_fun == IDEAL)
	{
		dat->EvalActivity = ideal_solution;
		dat->activity_data = dat;
	}
	else if (dat->act_fun == DAVIES)
	{
		dat->EvalActivity = Davies_equation;
		dat->activity_data = dat;
	}
	else if (dat->act_fun == DEBYE_HUCKEL)
	{
		dat->EvalActivity = DebyeHuckel_equation;
		dat->activity_data = dat;
	}
	//NOTE: The below lines of code will change after we define the PITZER and SIT models and what data structures they need
	else if (dat->act_fun == PITZER)
	{
		dat->EvalActivity = ideal_solution;
		dat->activity_data = dat;
	}
	else if (dat->act_fun == SIT)
	{
		dat->EvalActivity = ideal_solution;
		dat->activity_data = dat;
	}
	else
	{
		mError(invalid_boolean);
		std::cout << "Converting to ideal solution assumption.\n";
		dat->act_fun = IDEAL;
		dat->EvalActivity = ideal_solution;
	}

	//Setup PJFNK Solver options (May edit after calling setup function if desired)
	dat->Newton_data.LineSearch = true;
	dat->Newton_data.Bounce = false;
	if (dat->steadystate == true)
		dat->Newton_data.nl_maxit = 2 * dat->numvar;
	else
		dat->Newton_data.nl_maxit = dat->numvar;
	dat->Newton_data.nl_tol_abs = 1e-6;
	dat->Newton_data.nl_tol_rel = 1e-8;
	dat->Newton_data.lin_tol_abs = 1e-6;
	dat->Newton_data.lin_tol_rel = 1e-6;
	dat->Newton_data.NL_Output = dat->Console_Output;
	if (dat->steadystate == false)
		dat->Newton_data.NL_Output = false;
	if ( (*precond) == NULL && dat->numvar >= 100)
	{
		dat->Newton_data.linear_solver = GMRESRP;
	}
	else
	{
		dat->Newton_data.linear_solver = QR;
	}

	//Setup the memory working space for the problem
	dat->MasterList.set_list_size(dat->numvar);
	dat->X_old.set_size(dat->numvar, 1);
	dat->X_new.set_size(dat->numvar, 1);
	dat->Conc_old.set_size(dat->numvar, 1);
	dat->Conc_new.set_size(dat->numvar, 1);
	dat->activity_new.set_size(dat->numvar, 1);
	dat->activity_old.set_size(dat->numvar, 1);

	dat->ReactionList.resize(dat->num_ssr);
	dat->MassBalanceList.resize(dat->num_mbe);
	dat->UnsteadyList.resize(dat->num_usr);
	dat->AdsorptionList.resize(dat->num_ssao);
	dat->UnsteadyAdsList.resize(dat->num_usao);
	dat->ChemisorptionList.resize(dat->num_sschem);
	//dat->MultiAdsList.resize(dat->num_multi_ssao); //Action performed during read (this line may be redundant)
	//dat->MultiChemList.resize(dat->num_multi_sschem);	//Action performed during read (this line may be redundant)
	dat->OtherList.resize(dat->num_other);

	for (int i=0; i<dat->ReactionList.size(); i++)
	{
		dat->ReactionList[i].Initialize_Object(dat->MasterList);
	}

	for (int i=0; i<dat->MassBalanceList.size(); i++)
	{
		dat->MassBalanceList[i].Initialize_Object(dat->MasterList);
		dat->MassBalanceList[i].Set_SteadyState(dat->steadystate);
		dat->MassBalanceList[i].Set_TimeStep(dat->dt);
		dat->MassBalanceList[i].Set_Volume(dat->volume);
		dat->MassBalanceList[i].Set_FlowRate(dat->flow_rate);
		dat->MassBalanceList[i].Set_Area(dat->xsec_area);
		dat->MassBalanceList[i].Set_Type(dat->reactor_type);
		dat->MassBalanceList[i].Set_ZeroInitialSolids(dat->ZeroInitialSolids);
	}

	for (int i=0; i<dat->UnsteadyList.size(); i++)
	{
		dat->UnsteadyList[i].Initialize_Object(dat->MasterList);
	}

	for (int i=0; i<dat->AdsorptionList.size(); i++)
	{
		dat->AdsorptionList[i].Initialize_Object(dat->MasterList, dat->num_ssar[i]);
	}

	for (int i=0; i<dat->UnsteadyAdsList.size(); i++)
	{
		dat->UnsteadyAdsList[i].Initialize_Object(dat->MasterList, dat->num_usar[i]);
	}

	for (int i=0; i<dat->MultiAdsList.size(); i++)
	{
		for (int j=0; j<dat->num_multi_ssar[i].size(); j++)
		{
			dat->MultiAdsList[i].Initialize_Object(dat->MasterList, (int)dat->num_multi_ssar[i].size(), dat->num_multi_ssar[i]);
		}
	}

	for (int i=0; i<dat->ChemisorptionList.size(); i++)
	{
		dat->ChemisorptionList[i].Initialize_Object(dat->MasterList, dat->num_sschem_rxns[i]);
	}

	for (int i=0; i<dat->MultiChemList.size(); i++)
	{
		for (int j=0; j<dat->num_multichem_rxns[i].size(); j++)
		{
			dat->MultiChemList[i].Initialize_Object(dat->MasterList, (int)dat->num_multichem_rxns[i].size(), dat->num_multichem_rxns[i]);
		}
	}

	return success;
}

//Function to add a custom residual function to the solver
int shark_add_customResidual(int i, double (*other_res) (const Matrix<double> &x, SHARK_DATA *shark_dat, const void *other_data),
							 SHARK_DATA *shark_dat)
{
	int success = 0;

	if (i >= shark_dat->num_other || i < 0)
	{
		mError(out_of_bounds);
		return -1;
	}
	else if ( (*other_res) == NULL)
	{
		mError(nullptr_func);
		return -1;
	}
	else
	{
		shark_dat->OtherList[i] = other_res;
	}

	return success;
}

//Check for missing information or errors in given data
int shark_parameter_check(SHARK_DATA *shark_dat)
{
	int success = 0;

	for (int i=0; i<shark_dat->ReactionList.size(); i++)
	{
		shark_dat->ReactionList[i].checkSpeciesEnergies();
		if (shark_dat->ReactionList[i].haveEquilibrium() == false)
		{
			mError(missing_information);
			return -1;
		}
	}
	for (int i=0; i<shark_dat->UnsteadyList.size(); i++)
	{
		shark_dat->UnsteadyList[i].checkSpeciesEnergies();
		if (shark_dat->UnsteadyList[i].haveEquilibrium() == false || shark_dat->UnsteadyList[i].haveRate() == false)
		{
			mError(missing_information);
			return -1;
		}
	}
	for (int i=0; i<shark_dat->AdsorptionList.size(); i++)
	{
		for (int n=0; n<shark_dat->AdsorptionList[i].getNumberRxns(); n++)
		{
			shark_dat->AdsorptionList[i].getReaction(n).checkSpeciesEnergies();
			if (shark_dat->AdsorptionList[i].getReaction(n).haveEquilibrium() == false)
			{
				mError(missing_information);
				return -1;
			}
		}
		success = shark_dat->AdsorptionList[i].setAdsorbIndices();
		if (success != 0) {mError(missing_information); return -1;}
		success = shark_dat->AdsorptionList[i].setAqueousIndexAuto();
		if (success != 0) {mError(missing_information); return -1;}
		shark_dat->AdsorptionList[i].setTotalVolume(shark_dat->volume);
		success = shark_dat->AdsorptionList[i].checkAqueousIndices();
		if (success != 0) {mError(missing_information); return -1;}
		for (int m=0; m<shark_dat->MassBalanceList.size(); m++)
		{
			shark_dat->AdsorptionList[i].modifyDeltas(shark_dat->MassBalanceList[m]);
		}
	}

	for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
	{
		for (int n=0; n<shark_dat->UnsteadyAdsList[i].getNumberRxns(); n++)
		{
			shark_dat->UnsteadyAdsList[i].getReaction(n).checkSpeciesEnergies();
			if (shark_dat->UnsteadyAdsList[i].getReaction(n).haveEquilibrium() == false && shark_dat->UnsteadyAdsList[i].getReaction(n).haveRate() == false)
			{
				mError(missing_information);
				return -1;
			}
		}
		success = shark_dat->UnsteadyAdsList[i].setAdsorbIndices();
		if (success != 0) {mError(missing_information); return -1;}
		success = shark_dat->UnsteadyAdsList[i].setAqueousIndexAuto();
		if (success != 0) {mError(missing_information); return -1;}
		shark_dat->UnsteadyAdsList[i].setTotalVolume(shark_dat->volume);
		success = shark_dat->UnsteadyAdsList[i].checkAqueousIndices();
		if (success != 0) {mError(missing_information); return -1;}
		for (int m=0; m<shark_dat->MassBalanceList.size(); m++)
		{
			shark_dat->UnsteadyAdsList[i].modifyDeltas(shark_dat->MassBalanceList[m]);
		}
	}

	for (int i=0; i<shark_dat->MultiAdsList.size(); i++)
	{
		success = shark_dat->MultiAdsList[i].setAdsorbIndices();
		if (success != 0) {mError(missing_information); return -1;}
		success = shark_dat->MultiAdsList[i].setAqueousIndexAuto();
		if (success != 0) {mError(missing_information); return -1;}
		shark_dat->MultiAdsList[i].setTotalVolume(shark_dat->volume);
		success = shark_dat->MultiAdsList[i].checkAqueousIndices();
		if (success != 0) {mError(missing_information); return -1;}
		for (int m=0; m<shark_dat->MassBalanceList.size(); m++)
		{
			shark_dat->MultiAdsList[i].modifyDeltas(shark_dat->MassBalanceList[m]);
		}

		for (int j=0; j<shark_dat->MultiAdsList[i].getNumberLigands(); j++)
		{
			for (int n=0; n<shark_dat->MultiAdsList[i].getAdsorptionObject(j).getNumberRxns(); n++)
			{
				shark_dat->MultiAdsList[i].getAdsorptionObject(j).getReaction(n).checkSpeciesEnergies();
				if (shark_dat->MultiAdsList[i].getAdsorptionObject(j).getReaction(n).haveEquilibrium() == false)
				{
					mError(missing_information);
					return -1;
				}
			}
		}
	}

	for (int i=0; i<shark_dat->ChemisorptionList.size(); i++)
	{
		for (int n=0; n<shark_dat->ChemisorptionList[i].getNumberRxns(); n++)
		{
			shark_dat->ChemisorptionList[i].getReaction(n).checkSpeciesEnergies();
			if (shark_dat->ChemisorptionList[i].getReaction(n).haveEquilibrium() == false)
			{
				mError(missing_information);
				return -1;
			}
		}
		success = shark_dat->ChemisorptionList[i].setAdsorbIndices();
		if (success != 0) {mError(missing_information); return -1;}
		shark_dat->ChemisorptionList[i].setTotalVolume(shark_dat->volume);
		for (int m=0; m<shark_dat->MassBalanceList.size(); m++)
		{
			shark_dat->ChemisorptionList[i].modifyMBEdeltas(shark_dat->MassBalanceList[m]);
		}
	}

	for (int i=0; i<shark_dat->MultiChemList.size(); i++)
	{
		success = shark_dat->MultiChemList[i].setAdsorbIndices();
		if (success != 0) {mError(missing_information); return -1;}
		shark_dat->MultiChemList[i].setTotalVolume(shark_dat->volume);

		for (int m=0; m<shark_dat->MassBalanceList.size(); m++)
		{
			shark_dat->MultiChemList[i].modifyMBEdeltas(shark_dat->MassBalanceList[m]);
		}

		for (int j=0; j<shark_dat->MultiChemList[i].getNumberLigands(); j++)
		{
			for (int n=0; n<shark_dat->MultiChemList[i].getChemisorptionObject(j).getNumberRxns(); n++)
			{
				shark_dat->MultiChemList[i].getChemisorptionObject(j).getReaction(n).checkSpeciesEnergies();
				if (shark_dat->MultiChemList[i].getChemisorptionObject(j).getReaction(n).haveEquilibrium() == false)
				{
					mError(missing_information);
					return -1;
				}
			}
		}
	}


	return success;
}

//Calculate all reaction energy constants for the system
int shark_energy_calculations(SHARK_DATA *shark_dat)
{
	int success = 0;

	for (int i=0; i<shark_dat->ReactionList.size(); i++)
	{
		shark_dat->ReactionList[i].calculateEnergies();
	}
	for (int i=0; i<shark_dat->UnsteadyList.size(); i++)
	{
		shark_dat->UnsteadyList[i].calculateEnergies();
	}
	for (int i=0; i<shark_dat->AdsorptionList.size(); i++)
	{
		for (int n=0; n<shark_dat->AdsorptionList[i].getNumberRxns(); n++)
		{
			shark_dat->AdsorptionList[i].getReaction(n).calculateEnergies();
		}
	}
	for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
	{
		for (int n=0; n<shark_dat->UnsteadyAdsList[i].getNumberRxns(); n++)
		{
			shark_dat->UnsteadyAdsList[i].getReaction(n).calculateEnergies();
		}
	}

	for (int i=0; i<shark_dat->MultiAdsList.size(); i++)
	{
		for (int j=0; j<shark_dat->MultiAdsList[i].getNumberLigands(); j++)
		{
			for (int n=0; n<shark_dat->MultiAdsList[i].getAdsorptionObject(j).getNumberRxns(); n++)
			{
				shark_dat->MultiAdsList[i].getAdsorptionObject(j).getReaction(n).calculateEnergies();
			}
		}

	}

	for (int i=0; i<shark_dat->ChemisorptionList.size(); i++)
	{
		for (int n=0; n<shark_dat->ChemisorptionList[i].getNumberRxns(); n++)
		{
			shark_dat->ChemisorptionList[i].getReaction(n).calculateEnergies();
		}
	}

	for (int i=0; i<shark_dat->MultiChemList.size(); i++)
	{
		for (int j=0; j<shark_dat->MultiChemList[i].getNumberLigands(); j++)
		{
			for (int n=0; n<shark_dat->MultiChemList[i].getChemisorptionObject(j).getNumberRxns(); n++)
			{
				shark_dat->MultiChemList[i].getChemisorptionObject(j).getReaction(n).calculateEnergies();
			}
		}

	}

	return success;
}

//Calculate all constants needed given the system temperature
int shark_temperature_calculations(SHARK_DATA *shark_dat)
{
	int success = 0;

	for (int i=0; i<shark_dat->ReactionList.size(); i++)
	{
		shark_dat->ReactionList[i].calculateEquilibrium(shark_dat->temperature);
	}
	for (int i=0; i<shark_dat->UnsteadyList.size(); i++)
	{
		shark_dat->UnsteadyList[i].calculateRate(shark_dat->temperature);
	}
	for (int i=0; i<shark_dat->AdsorptionList.size(); i++)
	{
		shark_dat->AdsorptionList[i].calculateEquilibria(shark_dat->temperature);
	}
	for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
	{
		shark_dat->UnsteadyAdsList[i].calculateRates(shark_dat->temperature);
	}

	for (int i=0; i<shark_dat->MultiAdsList.size(); i++)
	{
		for (int j=0; j<shark_dat->MultiAdsList[i].getNumberLigands(); j++)
		{
			shark_dat->MultiAdsList[i].getAdsorptionObject(j).calculateEquilibria(shark_dat->temperature);
		}

	}
	for (int i=0; i<shark_dat->ChemisorptionList.size(); i++)
	{
		shark_dat->ChemisorptionList[i].calculateEquilibria(shark_dat->temperature);
	}

	for (int i=0; i<shark_dat->MultiChemList.size(); i++)
	{
		for (int j=0; j<shark_dat->MultiChemList[i].getNumberLigands(); j++)
		{
			shark_dat->MultiChemList[i].getChemisorptionObject(j).calculateEquilibria(shark_dat->temperature);
		}

	}

	return success;
}

//Go through SHARK data to initialize locations of pH and pOH
int shark_pH_finder(SHARK_DATA *shark_dat)
{
	int success = 0;

	bool found_pH = false, found_pOH = false;
	for (int i=0; i<shark_dat->MasterList.list_size(); i++)
	{
		if (shark_dat->MasterList.speciesName(i) == "H + (aq)")
		{
			shark_dat->Contains_pH = true;
			shark_dat->pH_index = i;
			if (found_pH == true)
			{
				mError(duplicate_variable);
				return -1;
			}
			found_pH = true;
		}
		else if (shark_dat->MasterList.speciesName(i) == "OH - (aq)")
		{
			shark_dat->Contains_pOH = true;
			shark_dat->pOH_index = i;
			if (found_pOH == true)
			{
				mError(duplicate_variable);
				return -1;
			}
			found_pOH = true;
		}
		else {/*No Action*/}

	}

	return success;
}

//Provide an initial guess to the non-linear system
int shark_guess(SHARK_DATA *shark_dat)
{
	int success = 0;

	shark_dat->Conc_new.ConstantICFill(0.0);
	shark_dat->activity_new.ConstantICFill(1.0);

	for (int i=0; i<shark_dat->MassBalanceList.size(); i++)
	{
		double delta_sum = shark_dat->MassBalanceList[i].Sum_Delta();
		double distribution;
		if (shark_dat->reactor_type == BATCH)
			distribution = shark_dat->MassBalanceList[i].Get_TotalConcentration() / delta_sum;
		else
		{
			if (shark_dat->steadystate == true)
				distribution = shark_dat->MassBalanceList[i].Get_InletConcentration() / delta_sum;
			else
				distribution = shark_dat->MassBalanceList[i].Get_InitialConcentration() / delta_sum;
		}
		for (int j=0; j<shark_dat->MasterList.list_size(); j++)
		{
			if (shark_dat->MassBalanceList[i].Get_Delta(j) > 0.0 && shark_dat->Conc_new(j,0) == 0.0)
				shark_dat->Conc_new.edit(j, 0, distribution);
			else if (shark_dat->MassBalanceList[i].Get_Delta(j) > 0.0 && shark_dat->Conc_new(j,0) != 0.0)
			{
				if (shark_dat->reactor_type == BATCH)
				{
					if (shark_dat->Conc_new(j,0) >= shark_dat->MassBalanceList[i].Get_TotalConcentration())
						shark_dat->Conc_new.edit(j, 0, distribution);
				}
				else
				{
					if (shark_dat->Conc_new(j,0) >= shark_dat->MassBalanceList[i].Get_InletConcentration())
						shark_dat->Conc_new.edit(j, 0, distribution);
				}
			}
		}
	}

	double chemi_sum = 0.0;
	for (int i=0; i<shark_dat->ChemisorptionList.size(); i++)
	{
		chemi_sum = 1.0;
		for (int j=0; j<shark_dat->MasterList.list_size(); j++)
		{
			chemi_sum = chemi_sum + shark_dat->ChemisorptionList[i].getDelta(j);
		}
		double dist = shark_dat->ChemisorptionList[i].getSpecificMolality() / chemi_sum;
		shark_dat->Conc_new.edit(shark_dat->ChemisorptionList[i].getLigandIndex(),0,dist);
		for (int j=0; j<shark_dat->MasterList.list_size(); j++)
		{
			if (shark_dat->ChemisorptionList[i].getDelta(j) > 0.0 && shark_dat->Conc_new(j,0) == 0.0)
				shark_dat->Conc_new.edit(j, 0, dist);
		}
		chemi_sum = 0.0;
	}

	for (int i=0; i<shark_dat->MultiChemList.size(); i++)
	{
		for (int l=0; l<shark_dat->MultiChemList[i].getNumberLigands(); l++)
		{
			chemi_sum = 1.0;
			for (int j=0; j<shark_dat->MasterList.list_size(); j++)
			{
				chemi_sum = chemi_sum + shark_dat->MultiChemList[i].getChemisorptionObject(l).getDelta(j);
			}
			double dist = shark_dat->MultiChemList[i].getChemisorptionObject(l).getSpecificMolality() / chemi_sum;
			shark_dat->Conc_new.edit(shark_dat->MultiChemList[i].getChemisorptionObject(l).getLigandIndex(),0,dist);
			for (int j=0; j<shark_dat->MasterList.list_size(); j++)
			{
				if (shark_dat->MultiChemList[i].getChemisorptionObject(l).getDelta(j) > 0.0 && shark_dat->Conc_new(j,0) == 0.0)
					shark_dat->Conc_new.edit(j, 0, dist);
			}
			chemi_sum = 0.0;
		}
	}

	if (shark_dat->steadystate == false)
	{
		for (int i=0; i<shark_dat->UnsteadyList.size(); i++)
		{
			double max = 0.0, current = 0.0;
			shark_dat->Conc_new.edit(shark_dat->UnsteadyList[i].Get_Species_Index(),0,shark_dat->UnsteadyList[i].Get_InitialValue());

			//Loop through MassBalanceList and MasterList to find max value for unsteady species
			for (int j=0; j<shark_dat->MassBalanceList.size(); j++)
			{
				if (shark_dat->MassBalanceList[j].Get_Delta(shark_dat->UnsteadyList[i].Get_Species_Index()) != 0.0)
				{
					if (shark_dat->reactor_type == BATCH)
					{
						if (max == 0)
							max = shark_dat->MassBalanceList[j].Get_TotalConcentration();
						current = shark_dat->MassBalanceList[j].Get_TotalConcentration();
						if (max > current)
							max = shark_dat->MassBalanceList[j].Get_TotalConcentration();
					}
					else
					{
						if (max == 0)
							max = shark_dat->MassBalanceList[j].Get_InletConcentration();
						current = shark_dat->MassBalanceList[j].Get_InletConcentration();
						if (max > current)
							max = shark_dat->MassBalanceList[j].Get_InletConcentration();
					}
				}
			}
			shark_dat->UnsteadyList[i].Set_MaximumValue(max);
		}

		//Loop for all Unsteady Adsorption Objects
		for (int n=0; n<shark_dat->UnsteadyAdsList.size(); n++)
		{
			//Loop for all reactions in the nth unsteady adsorption object
			for (int i=0; i<shark_dat->UnsteadyAdsList[n].getNumberRxns(); i++)
			{
				double max = 0.0, current = 0.0;
				shark_dat->Conc_new.edit(shark_dat->UnsteadyAdsList[n].getReaction(i).Get_Species_Index(),0,shark_dat->UnsteadyAdsList[n].getReaction(i).Get_InitialValue());

				//Loop through MassBalanceList and MasterList to find max value for unsteady species
				for (int j=0; j<shark_dat->MassBalanceList.size(); j++)
				{
					if (shark_dat->MassBalanceList[j].Get_Delta(shark_dat->UnsteadyAdsList[n].getReaction(i).Get_Species_Index()) != 0.0)
					{
						if (shark_dat->reactor_type == BATCH)
						{
							if (max == 0)
								max = shark_dat->MassBalanceList[j].Get_TotalConcentration()/shark_dat->UnsteadyAdsList[n].getBulkDensity();
							current = shark_dat->MassBalanceList[j].Get_TotalConcentration()/shark_dat->UnsteadyAdsList[n].getBulkDensity();
							if (max > current)
								max = shark_dat->MassBalanceList[j].Get_TotalConcentration()/shark_dat->UnsteadyAdsList[n].getBulkDensity();
						}
						else
						{
							if (max == 0)
								max = shark_dat->MassBalanceList[j].Get_InletConcentration()/shark_dat->UnsteadyAdsList[n].getBulkDensity();
							current = shark_dat->MassBalanceList[j].Get_InletConcentration()/shark_dat->UnsteadyAdsList[n].getBulkDensity();
							if (max > current)
								max = shark_dat->MassBalanceList[j].Get_InletConcentration()/shark_dat->UnsteadyAdsList[n].getBulkDensity();
						}
					}
				}
				shark_dat->UnsteadyAdsList[n].getReaction(i).Set_MaximumValue(max);

			}
		}
	}

	//Set pH information if needed
	if (shark_dat->Contains_pH == true)
	{
		if (shark_dat->const_pH == true)
			shark_dat->Conc_new.edit(shark_dat->pH_index, 0, pow(10.0,-shark_dat->pH));
		else
			shark_dat->Conc_new.edit(shark_dat->pH_index, 0, 10.0);
	}
	if (shark_dat->Contains_pOH == true)
	{
		if (shark_dat->const_pH == true)
			shark_dat->Conc_new.edit(shark_dat->pOH_index, 0, pow(10.0,-(14.0-shark_dat->pH)));
		else
			shark_dat->Conc_new.edit(shark_dat->pOH_index, 0, 10.0);
	}

	success = Convert2LogConcentration(shark_dat->Conc_new, shark_dat->X_new);
	if (success != 0) {mError(simulation_fail); return -1;}

	return success;
}

//Provide the initial conditions for the system
int shark_initial_conditions(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Printout the header for simulation info
	print2file_shark_header(shark_dat);

	if (shark_dat->Console_Output == true)
	{
		if (shark_dat->steadystate == false)
			std::cout << "------------Establishing Initial Conditions for SHARK Simulations--------------\n";
		else
			std::cout << "------------Establishing Initial Guess for SHARK Simulations--------------\n";
	}

	//Make initial guess for the system
	if (shark_dat->steadystate == true)
	{
		success = shark_guess(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

		shark_dat->activity_old = shark_dat->activity_new;
		shark_dat->Conc_old = shark_dat->Conc_new;
		shark_dat->X_old = shark_dat->X_new;

	}
	//Make an initial guess and solve the initial speciation
	else
	{
		shark_dat->time = 0.0;
		success = shark_guess(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

		success = pjfnk(shark_dat->Residual,shark_dat->lin_precon,shark_dat->X_new,&shark_dat->Newton_data,shark_dat->residual_data,shark_dat->precon_data);
		shark_dat->totalsteps = shark_dat->totalsteps + shark_dat->Newton_data.nl_iter + shark_dat->Newton_data.l_iter;
		shark_dat->totalcalls = shark_dat->totalcalls + shark_dat->Newton_data.fun_call;
		if (success != 0) {mError(simulation_fail); return -1;}

		shark_dat->Norm = shark_dat->Newton_data.nl_res;
		if (shark_dat->Norm <= shark_dat->Newton_data.nl_tol_abs || shark_dat->Newton_data.nl_relres <= shark_dat->Newton_data.nl_tol_rel)
			shark_dat->Converged = true;
		else
			shark_dat->Converged = false;

		if (shark_dat->ZeroInitialSolids == true)
		{
			//Loop and remove solids from solution
			for (int i=0; i<shark_dat->MasterList.list_size(); i++)
			{
				if (shark_dat->MasterList.get_species(i).MoleculePhaseID() == SOLID || shark_dat->MasterList.get_species(i).MoleculePhaseID() == ADSORBED)
					shark_dat->X_new.edit(i, 0, -DBL_MAX_10_EXP);
			}
		}

		success = Convert2Concentration(shark_dat->X_new, shark_dat->Conc_new);
		if (success != 0) {mError(simulation_fail); return -1;}

		shark_dat->activity_old = shark_dat->activity_new;
		shark_dat->Conc_old = shark_dat->Conc_new;
		shark_dat->X_old = shark_dat->X_new;

		//Printout ICs for Unsteady case
		print2file_shark_results_new(shark_dat);
	}

	if (shark_dat->Console_Output == true)
	{
		shark_dat->MasterList.DisplayConcentrations(shark_dat->Conc_new);
		std::cout << "E. Norm =\t" << shark_dat->Norm << "\nIterations =\t" << shark_dat->Newton_data.nl_iter << std::endl;
		if (shark_dat->steadystate == false)
			std::cout << "\n-------------Initial Conditions Set!-------------\n\n";
		else
			std::cout << "\n-------------Initial Guess Set!-------------\n\n";
	}

	return success;
}

//Executioner function for shark
int shark_executioner(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Loop the solver till solution is obtained or dt_min is reached
	do
	{
		//Call the preprocess function
		success = shark_preprocesses(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

		//Call the solver function
		success = shark_solver(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

	} while (shark_dat->dt > shark_dat->dt_min && shark_dat->Converged == false && shark_dat->steadystate == false);

	//Call the postprocess function
	success = shark_postprocesses(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	return success;
}

//Time step function for shark
int shark_timestep_const(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->steadystate == false)
	{
		shark_dat->time = shark_dat->time_old + shark_dat->dt;
		if (shark_dat->time > shark_dat->simulationtime)
		{
			shark_dat->time = shark_dat->simulationtime;
			shark_dat->dt = shark_dat->time - shark_dat->time_old;
		}
	}
	for (int i=0; i<shark_dat->UnsteadyList.size(); i++)
	{
		shark_dat->UnsteadyList[i].Set_TimeStep(shark_dat->dt);
	}
	for (int i=0; i<shark_dat->MassBalanceList.size(); i++)
	{
		shark_dat->MassBalanceList[i].Set_TimeStep(shark_dat->dt);
	}
	for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
	{
		for (int n=0; n<shark_dat->UnsteadyAdsList[i].getNumberRxns(); n++)
		{
			shark_dat->UnsteadyAdsList[i].getReaction(n).Set_TimeStep(shark_dat->dt);
		}
	}
	return success;
}

//Adaptive time step function for difficult problems
int shark_timestep_adapt(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->steadystate == false)
	{
		if (shark_dat->Converged == true)
		{
			shark_dat->dt = shark_dat->dt * 1.5;
			if (shark_dat->dt >= shark_dat->dt_max)
				shark_dat->dt = shark_dat->dt_max;
		}
		else
		{
			shark_dat->dt = shark_dat->dt * 0.5;
			if (shark_dat->dt <= shark_dat->dt_min)
				shark_dat->dt = shark_dat->dt_min;
		}
		shark_dat->time = shark_dat->time_old + shark_dat->dt;
		if (shark_dat->time > shark_dat->simulationtime)
		{
			shark_dat->time = shark_dat->simulationtime;
			shark_dat->dt = shark_dat->time - shark_dat->time_old;
		}
	}
	for (int i=0; i<shark_dat->UnsteadyList.size(); i++)
	{
		shark_dat->UnsteadyList[i].Set_TimeStep(shark_dat->dt);
	}
	for (int i=0; i<shark_dat->MassBalanceList.size(); i++)
	{
		shark_dat->MassBalanceList[i].Set_TimeStep(shark_dat->dt);
	}
	for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
	{
		for (int n=0; n<shark_dat->UnsteadyAdsList[i].getNumberRxns(); n++)
		{
			shark_dat->UnsteadyAdsList[i].getReaction(n).Set_TimeStep(shark_dat->dt);
		}
	}

	return success;
}

//Preprocess function for shark
int shark_preprocesses(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Call the timestep function
	if ((shark_dat->TimeAdaptivity == false && shark_dat->Converged == true) || shark_dat->time == 0.0)
	{
		success = shark_timestep_const(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	else
	{
		success = shark_timestep_adapt(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}
	}

	//Function to calculate equilibrium and rate constants as a function of temperature
	success = shark_temperature_calculations(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	return success;
}

//Solver function for shark
int shark_solver(SHARK_DATA *shark_dat)
{
	int success = 0;

	if (shark_dat->steadystate == false && shark_dat->num_usr > 0)
	{
		if (shark_dat->Console_Output == true)
			std::cout << "Explicit Approximation to Unsteady Aqueous Species...\n-----------------------------------------------\n";
		for (int i=0; i<shark_dat->num_usr; i++)
		{
			shark_dat->Conc_new.edit(shark_dat->UnsteadyList[i].Get_Species_Index(), 0, shark_dat->UnsteadyList[i].Explicit_Eval(shark_dat->X_old, shark_dat->activity_old));

			if (shark_dat->Console_Output == true)
			{
				std::cout << "[ " << shark_dat->MasterList.get_species(shark_dat->UnsteadyList[i].Get_Species_Index()).MolecularFormula() << " ] =\t" << shark_dat->Conc_new(shark_dat->UnsteadyList[i].Get_Species_Index(),0) << std::endl;
			}
		}
		if (shark_dat->Console_Output == true)
			std::cout << "\n";
		success = Convert2LogConcentration(shark_dat->Conc_new, shark_dat->X_new);
		if (success != 0) {mError(simulation_fail); return -1;}
	}
	if (shark_dat->steadystate == false && shark_dat->num_usao > 0)
	{
		if (shark_dat->Console_Output == true)
			std::cout << "Explicit Approximation to Unsteady Adsorbed Species...\n-----------------------------------------------\n";
		for (int i=0; i<shark_dat->num_usao; i++)
		{
			for (int j=0; j<shark_dat->UnsteadyAdsList[i].getNumberRxns(); j++)
			{
				shark_dat->Conc_new.edit(shark_dat->UnsteadyAdsList[i].getAdsorbIndex(j),0, shark_dat->UnsteadyAdsList[i].Explicit_Eval(shark_dat->X_old, shark_dat->activity_old, shark_dat->temperature, shark_dat->relative_permittivity, j));

				if (shark_dat->Console_Output == true)
				{
					std::cout << "[ " << shark_dat->MasterList.get_species(shark_dat->UnsteadyAdsList[i].getAdsorbIndex(j)).MolecularFormula() << " ] =\t" << shark_dat->Conc_new(shark_dat->UnsteadyAdsList[i].getAdsorbIndex(j),0) << std::endl;
				}
			}
		}
		if (shark_dat->Console_Output == true)
			std::cout << "\n";
		success = Convert2LogConcentration(shark_dat->Conc_new, shark_dat->X_new);
		if (success != 0) {mError(simulation_fail); return -1;}
	}

	if (shark_dat->steadystate == false && shark_dat->time_old == 0.0)
	{
		//Loop through all mass balances and make corrections to those who were initially zero
		for (int i=0; i<shark_dat->MassBalanceList.size(); i++)
		{
			if (shark_dat->MassBalanceList[i].Get_InitialConcentration() <= DBL_MIN && shark_dat->reactor_type != BATCH)
			{
				double delta_sum = shark_dat->MassBalanceList[i].Sum_Delta();
				double distribution;
				distribution = shark_dat->MassBalanceList[i].Get_InletConcentration() / delta_sum;
				for (int j=0; j<shark_dat->MasterList.list_size(); j++)
				{
					if (shark_dat->MassBalanceList[i].Get_Delta(j) > 0.0)
						shark_dat->Conc_new.edit(j, 0, distribution);
				}
			}
		}
		success = Convert2LogConcentration(shark_dat->Conc_new, shark_dat->X_new);
		if (success != 0) {mError(simulation_fail); return -1;}

	}

	success= pjfnk(shark_dat->Residual,shark_dat->lin_precon,shark_dat->X_new,&shark_dat->Newton_data,shark_dat->residual_data,shark_dat->precon_data);
	shark_dat->ionic_strength = calculate_ionic_strength(shark_dat->X_new, shark_dat->MasterList);
	shark_dat->totalsteps = shark_dat->totalsteps + shark_dat->Newton_data.nl_iter + shark_dat->Newton_data.l_iter;
	shark_dat->totalcalls = shark_dat->totalcalls + shark_dat->Newton_data.fun_call;
	if (success != 0) {mError(simulation_fail); return -1;}

	success = Convert2Concentration(shark_dat->X_new, shark_dat->Conc_new);
	if (success != 0) {mError(simulation_fail); return -1;}


	if (shark_dat->steadystate == true && shark_dat->reactor_type != BATCH)
	{
		//Loop through all mass balances and make corrections to those who have zero inlet
		for (int i=0; i<shark_dat->MassBalanceList.size(); i++)
		{
			if (shark_dat->MassBalanceList[i].Get_InletConcentration() <= DBL_MIN)
			{
				for (int j=0; j<shark_dat->MasterList.list_size(); j++)
				{
					if (shark_dat->MassBalanceList[i].Get_Delta(j) > 0.0)
						shark_dat->Conc_new.edit(j, 0, 0.0);
				}
			}
		}
		success = Convert2LogConcentration(shark_dat->Conc_new, shark_dat->X_new);
		if (success != 0) {mError(simulation_fail); return -1;}
	}

	shark_dat->Norm = shark_dat->Newton_data.nl_res;
	if (shark_dat->Norm <= shark_dat->Newton_data.nl_tol_abs || shark_dat->Newton_data.nl_relres <= shark_dat->Newton_data.nl_tol_rel)
		shark_dat->Converged = true;
	else
		shark_dat->Converged = false;

	if (shark_dat->Console_Output == true)
	{
		std::cout << "Evaluation Information...\n-------------------------\n";
		if (shark_dat->steadystate == false)
			std::cout << "dt =\t" << shark_dat->dt << std::endl;
		std::cout << "E. Norm =\t" << shark_dat->Norm << "\nIterations =\t" << shark_dat->Newton_data.nl_iter << std::endl;

		if (shark_dat->AdsorptionList.size() > 0)
			std::cout << "\nAdsorption info...\n-------------------------\n";
		for (int i=0; i<shark_dat->AdsorptionList.size(); i++)
		{
			std::cout << "Active Surface Fraction for " << shark_dat->AdsorptionList[i].getAdsorbentName() << " =\t" << shark_dat->AdsorptionList[i].calculateActiveFraction(shark_dat->X_new) << std::endl;
			std::cout << "Surface Charge Density (C/m^2) for " << shark_dat->AdsorptionList[i].getAdsorbentName() << " =\t" << shark_dat->AdsorptionList[i].getChargeDensity() << std::endl;
			std::cout << "logK values for " << shark_dat->AdsorptionList[i].getAdsorbentName() << "...\n";

			for (int j=0; j<shark_dat->AdsorptionList[i].getNumberRxns(); j++)
			{
				double logK = shark_dat->AdsorptionList[i].getReaction(j).Get_Equilibrium();

				if (shark_dat->AdsorptionList[i].includeSurfaceCharge() == true)
				{
					logK = logK + ((shark_dat->AdsorptionList[i].calculateEquilibriumCorrection(shark_dat->AdsorptionList[i].getChargeDensity(), shark_dat->temperature, shark_dat->AdsorptionList[i].getIonicStrength(), shark_dat->relative_permittivity, j))/log(10.0));
				}
				std::cout << "\tlogK(" << j << ") =\t" << logK << std::endl;
			}
			std::cout << "\ngama values for " << shark_dat->AdsorptionList[i].getAdsorbentName() << "...\n";
			for (int j=0; j<shark_dat->AdsorptionList[i].getNumberRxns(); j++)
			{
				double gama = shark_dat->AdsorptionList[i].getActivity(shark_dat->AdsorptionList[i].getAdsorbIndex(j));
				std::cout << "\tgama(" << j << ") =\t" << gama << std::endl;
			}
			std::cout << std::endl;
		}

		if (shark_dat->UnsteadyAdsList.size() > 0)
			std::cout << "\nUnsteady Adsorption info...\n-------------------------\n";
		for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
		{
			std::cout << "Active Surface Fraction for " << shark_dat->UnsteadyAdsList[i].getAdsorbentName() << " =\t" << shark_dat->UnsteadyAdsList[i].calculateActiveFraction(shark_dat->X_new) << std::endl;
			std::cout << "Surface Charge Density (C/m^2) for " << shark_dat->UnsteadyAdsList[i].getAdsorbentName() << " =\t" << shark_dat->UnsteadyAdsList[i].getChargeDensity() << std::endl;
			std::cout << "logK values for " << shark_dat->UnsteadyAdsList[i].getAdsorbentName() << "...\n";

			for (int j=0; j<shark_dat->UnsteadyAdsList[i].getNumberRxns(); j++)
			{
				double logK = shark_dat->UnsteadyAdsList[i].getReaction(j).Get_Equilibrium();

				if (shark_dat->UnsteadyAdsList[i].includeSurfaceCharge() == true)
				{
					logK = logK + ((shark_dat->UnsteadyAdsList[i].calculateEquilibriumCorrection(shark_dat->UnsteadyAdsList[i].getChargeDensity(), shark_dat->temperature, shark_dat->UnsteadyAdsList[i].getIonicStrength(), shark_dat->relative_permittivity, j))/log(10.0));
				}
				std::cout << "\tlogK(" << j << ") =\t" << logK << std::endl;
			}
			std::cout << "\ngama values for " << shark_dat->UnsteadyAdsList[i].getAdsorbentName() << "...\n";
			for (int j=0; j<shark_dat->UnsteadyAdsList[i].getNumberRxns(); j++)
			{
				double gama = shark_dat->UnsteadyAdsList[i].getActivity(shark_dat->UnsteadyAdsList[i].getAdsorbIndex(j));
				std::cout << "\tgama(" << j << ") =\t" << gama << std::endl;
			}
			std::cout << std::endl;
		}

		if (shark_dat->MultiAdsList.size() > 0)
			std::cout << "\nMultiligand Adsorption info...\n-------------------------\n";
		for (int i=0; i<shark_dat->MultiAdsList.size(); i++)
		{
			std::cout << "Adsorbent: " << shark_dat->MultiAdsList[i].getAdsorbentName() << "\n";
			std::cout << "Specific Surface Area (m^2/kg) = \t" << shark_dat->MultiAdsList[i].getSpecificArea() << std::endl;
			std::cout << "Surface Charge Density (C/m^2) =\t" << shark_dat->MultiAdsList[i].getChargeDensity() << std::endl;
			std::cout << "Electric Surface Potential (V) =\t" << shark_dat->MultiAdsList[i].getElectricPotential() << std::endl;
			std::cout << std::endl;
			for (int j=0; j<shark_dat->MultiAdsList[i].getNumberLigands(); j++)
			{
				std::cout << "Specific Molality for " << shark_dat->MultiAdsList[i].getLigandName(j) << " =\t" << shark_dat->MultiAdsList[i].getAdsorptionObject(j).getSpecificMolality() << std::endl;
				std::cout << "Active Surface Fraction for " << shark_dat->MultiAdsList[i].getLigandName(j) << " =\t" << shark_dat->MultiAdsList[i].getAdsorptionObject(j).calculateActiveFraction(shark_dat->X_new) << std::endl;
				std::cout << "logK values for this ligand...\n";

				for (int n=0; n<shark_dat->MultiAdsList[i].getAdsorptionObject(j).getNumberRxns(); n++)
				{
					double logK = shark_dat->MultiAdsList[i].getAdsorptionObject(j).getReaction(n).Get_Equilibrium();

					if (shark_dat->MultiAdsList[i].includeSurfaceCharge() == true)
					{
						logK = logK + ((shark_dat->MultiAdsList[i].calculateEquilibriumCorrection(shark_dat->MultiAdsList[i].getChargeDensity(), shark_dat->temperature, shark_dat->MultiAdsList[i].getIonicStrength(), shark_dat->relative_permittivity, n,j))/log(10.0));
					}
					std::cout << "\tlogK(" << n << ") =\t" << logK << std::endl;
				}
				std::cout << "\ngama values for this ligand...\n";
				for (int n=0; n<shark_dat->MultiAdsList[i].getAdsorptionObject(j).getNumberRxns(); n++)
				{
					double gama = shark_dat->MultiAdsList[i].getAdsorptionObject(j).getActivity(shark_dat->MultiAdsList[i].getAdsorptionObject(j).getAdsorbIndex(n));
					std::cout << "\tgama(" << n << ") =\t" << gama << std::endl;
				}
				std::cout << std::endl;
			}
		}

		if (shark_dat->ChemisorptionList.size() > 0)
			std::cout << "\nChemisorption info...\n-------------------------\n";
		for (int i=0; i<shark_dat->ChemisorptionList.size(); i++)
		{
			std::cout << "Surface Charge Density (C/m^2) for " << shark_dat->ChemisorptionList[i].getAdsorbentName() << " =\t" << shark_dat->ChemisorptionList[i].getChargeDensity() << std::endl;
			std::cout << "logK values for " << shark_dat->ChemisorptionList[i].getAdsorbentName() << "...\n";

			for (int j=0; j<shark_dat->ChemisorptionList[i].getNumberRxns(); j++)
			{
				double logK = shark_dat->ChemisorptionList[i].getReaction(j).Get_Equilibrium();

				if (shark_dat->ChemisorptionList[i].includeSurfaceCharge() == true)
				{
					logK = logK + ((shark_dat->ChemisorptionList[i].calculateEquilibriumCorrection(shark_dat->ChemisorptionList[i].getChargeDensity(), shark_dat->temperature, shark_dat->ChemisorptionList[i].getIonicStrength(), shark_dat->relative_permittivity, j))/log(10.0));
				}
				std::cout << "\tlogK(" << j << ") =\t" << logK << std::endl;
			}
			std::cout << "\ngama values for " << shark_dat->ChemisorptionList[i].getAdsorbentName() << "...\n";
			for (int j=0; j<shark_dat->ChemisorptionList[i].getNumberRxns(); j++)
			{
				double gama = shark_dat->ChemisorptionList[i].getActivity(shark_dat->ChemisorptionList[i].getAdsorbIndex(j));
				std::cout << "\tgama(" << j << ") =\t" << gama << std::endl;
			}
			std::cout << std::endl;
		}

		if (shark_dat->MultiChemList.size() > 0)
			std::cout << "\nMultiligand Chemisorption info...\n-------------------------\n";
		for (int i=0; i<shark_dat->MultiChemList.size(); i++)
		{
			std::cout << "Adsorbent: " << shark_dat->MultiChemList[i].getAdsorbentName() << "\n";
			std::cout << "Specific Surface Area (m^2/kg) = \t" << shark_dat->MultiChemList[i].getSpecificArea() << std::endl;
			std::cout << "Surface Charge Density (C/m^2) =\t" << shark_dat->MultiChemList[i].getChargeDensity() << std::endl;
			std::cout << "Electric Surface Potential (V) =\t" << shark_dat->MultiChemList[i].getElectricPotential() << std::endl;
			std::cout << std::endl;
			for (int j=0; j<shark_dat->MultiChemList[i].getNumberLigands(); j++)
			{
				std::cout << "Specific Molality for " << shark_dat->MultiChemList[i].getLigandName(j) << " =\t" << shark_dat->MultiChemList[i].getChemisorptionObject(j).getSpecificMolality() << std::endl;
				std::cout << "logK values for this ligand...\n";

				for (int n=0; n<shark_dat->MultiChemList[i].getChemisorptionObject(j).getNumberRxns(); n++)
				{
					double logK = shark_dat->MultiChemList[i].getChemisorptionObject(j).getReaction(n).Get_Equilibrium();

					if (shark_dat->MultiChemList[i].includeSurfaceCharge() == true)
					{
						logK = logK + ((shark_dat->MultiChemList[i].calculateEquilibriumCorrection(shark_dat->MultiChemList[i].getChargeDensity(), shark_dat->temperature, shark_dat->MultiChemList[i].getIonicStrength(), shark_dat->relative_permittivity, n,j))/log(10.0));
					}
					std::cout << "\tlogK(" << n << ") =\t" << logK << std::endl;
				}
				std::cout << "\ngama values for this ligand...\n";
				for (int n=0; n<shark_dat->MultiChemList[i].getChemisorptionObject(j).getNumberRxns(); n++)
				{
					double gama = shark_dat->MultiChemList[i].getChemisorptionObject(j).getActivity(shark_dat->MultiChemList[i].getChemisorptionObject(j).getAdsorbIndex(n));
					std::cout << "\tgama(" << n << ") =\t" << gama << std::endl;
				}
				std::cout << std::endl;
			}
		}

		if (shark_dat->Converged == false)
		{
			if (shark_dat->dt > shark_dat->dt_min)
			{
				if (shark_dat->steadystate == false)
				{
					std::cout << "Time step failed... Reducing dt...\n\n";
					return success;
				}
			}
			else
			{
				std::cout << "Time step cannot be reduced further!";
				shark_dat->LocalMin = false;
			}

			if (shark_dat->LocalMin == false)
			{
				std::cout << "\n--------------Force Quiting SHARK----------------\n\n";
				shark_dat->Newton_data.F.Display("Residual Vector");
				shark_dat->Newton_data.x.Display("Solution Vector");

				//Conversion to concentration units
				success = Convert2Concentration(shark_dat->X_new, shark_dat->Conc_new);
				if (success != 0) {mError(simulation_fail); return -1;}
				shark_dat->MasterList.DisplayConcentrations(shark_dat->Conc_new);

				std::cout << "Ionic Strength (M) = " << calculate_ionic_strength(shark_dat->X_new, shark_dat->MasterList) << std::endl;

				for (int i=0; i<shark_dat->AdsorptionList.size(); i++)
				{
					std::cout << "Active Surface Fraction " << i << " = " << shark_dat->AdsorptionList[i].calculateActiveFraction(shark_dat->X_new) << std::endl;
					std::cout << "Surface Charge Density (C/m^2) " << i << " = " << shark_dat->AdsorptionList[i].getChargeDensity() << std::endl;
					for (int n=0; n<shark_dat->AdsorptionList[i].getNumberRxns(); n++)
					{
						std::cout << "Area Factor for Species " << shark_dat->AdsorptionList[i].getAdsorbIndex(n) << " = " << shark_dat->AdsorptionList[i].getAreaFactor(shark_dat->AdsorptionList[i].getAdsorbIndex(n)) << std::endl;
					}
					for (int n=0; n<shark_dat->AdsorptionList[i].getNumberRxns(); n++)
					{
						std::cout << "logK for Reaction " << n << " = " << shark_dat->AdsorptionList[i].getReaction(n).Get_Equilibrium() << std::endl;
					}
				}
				for (int i=0; i<shark_dat->MultiAdsList.size(); i++)
				{
					std::cout << "Adsorbent: " << shark_dat->MultiAdsList[i].getAdsorbentName() << "\n";
					std::cout << "Specific Surface Area (m^2/kg) = \t" << shark_dat->MultiAdsList[i].getSpecificArea() << std::endl;
					std::cout << "Surface Charge Density (C/m^2) for =\t" << shark_dat->MultiAdsList[i].getChargeDensity() << std::endl;

					for (int j=0; j<shark_dat->MultiAdsList[i].getNumberLigands(); j++)
					{
						std::cout << "Active Surface Fraction for " << shark_dat->MultiAdsList[i].getLigandName(j) << " =\t" << shark_dat->MultiAdsList[i].getAdsorptionObject(j).calculateActiveFraction(shark_dat->X_new) << std::endl;
						std::cout << "logK values for this ligand...\n";

						for (int n=0; n<shark_dat->MultiAdsList[i].getAdsorptionObject(j).getNumberRxns(); n++)
						{
							double logK = shark_dat->MultiAdsList[i].getAdsorptionObject(j).getReaction(n).Get_Equilibrium();

							if (shark_dat->MultiAdsList[i].includeSurfaceCharge() == true)
							{
								logK = logK + ((shark_dat->MultiAdsList[i].calculateEquilibriumCorrection(shark_dat->MultiAdsList[i].getChargeDensity(), shark_dat->temperature, shark_dat->MultiAdsList[i].getIonicStrength(), shark_dat->relative_permittivity, n,j))/log(10.0));
							}
							std::cout << "\tlogK(" << n << ") =\t" << logK << std::endl;
						}
					}
				}
				std::cout << "\n";
				shark_dat->activity_new.Display("activities");

				//Form a Numerical Jacobian and print out
				if (shark_dat->numvar <= 100)
				{
					Matrix<double> J(shark_dat->numvar,shark_dat->numvar);
					NUM_JAC_DATA num_jac;
					success = NumericalJacobian(shark_dat->Residual, shark_dat->X_new, J, shark_dat->numvar, shark_dat->numvar, &num_jac, shark_dat->residual_data);

					FILE *Jacobian;

					Jacobian = fopen("output/SHARK_Jacobian.txt","w+");
					if (Jacobian == nullptr)
					{
						success = system("mkdir output");
						Jacobian = fopen("output/SHARK_Jacobian.txt", "w+");
					}

					//Loop to print out Jacobian to file
					fprintf(Jacobian,"Jacobian Matrix = \n\n");
					for (int i=0; i<shark_dat->numvar; i++)
					{
						for (int j=0; j< shark_dat->numvar; j++)
						{
							fprintf(Jacobian, "%.6g\t",J(i,j));
						}
						fprintf(Jacobian,"\n");
					}
					fprintf(Jacobian,"\n");

					fclose(Jacobian);
				}

				mError(simulation_fail);
				return -1;

			}
			else
			{
				std::cout << "\nLocal minimum was found...\n\n";
				success = 0;

			}
		}
	}

	return success;
}

//Postprocess function for shark
int shark_postprocesses(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Conversion to concentration units
	success = Convert2Concentration(shark_dat->X_new, shark_dat->Conc_new);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Check for failure
	if (shark_dat->Converged == false && shark_dat->Console_Output == true)
	{
		if (shark_dat->LocalMin == false)
			std::cout << "\nSHARK CONVERGENCE FAILURE!\n\n";
		else
		{
			if (shark_dat->Newton_data.nl_relres <= 0.1)
				std::cout << "\nLOCAL MINIMUM ACCEPTED AS SUCCESS!\n\n";
			else
			{
				std::cout << "\nWARNING! BAD LOCAL MINIMUM!\n\n";
				success = -1;
			}
		}
	}

	//Print to file
	if (shark_dat->steadystate == true)
		print2file_shark_results_new(shark_dat);
	else
	{
		shark_dat->timesteps++;
		shark_dat->t_count = shark_dat->t_count + shark_dat->dt;
		if (shark_dat->t_count >= (shark_dat->t_out+sqrt(DBL_EPSILON))
			|| shark_dat->t_count >= (shark_dat->t_out-sqrt(DBL_EPSILON))
			 || shark_dat->time == shark_dat->simulationtime)
		{
			print2file_shark_results_new(shark_dat);
			shark_dat->t_count = 0.0;
		}
	}

	return success;
}

//Reset function for shark
int shark_reset(SHARK_DATA *shark_dat)
{
	int success = 0;

	shark_dat->Conc_old = shark_dat->Conc_new;
	shark_dat->X_old = shark_dat->X_new;
	shark_dat->activity_old = shark_dat->activity_new;
	shark_dat->time_old = shark_dat->time;

	for (int i=0; i<shark_dat->UnsteadyAdsList.size(); i++)
	{
			shark_dat->UnsteadyAdsList[i].updateActivities();
	}

	return success;
}

//Default residual function for shark evaluations
int shark_residual(const Matrix<double> &x, Matrix<double> &F, const void *data)
{
	int success = 0;
	SHARK_DATA *dat = (SHARK_DATA *) data;

	//Call activity function
	success = dat->EvalActivity(x,dat->activity_new,dat->activity_data);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Form all residuals
	int index = 0;
	for (int i=0; i<dat->ReactionList.size(); i++)
	{
		F(index,0) = dat->ReactionList[i].Eval_Residual(x, dat->activity_new);
		index++;
	}
	for (int i=0; i<dat->UnsteadyList.size(); i++)
	{
		if (dat->steadystate == false)
		{
			if (dat->time == 0.0)
				F(index,0) = dat->UnsteadyList[i].Eval_IC_Residual(x);
			else
				F(index,0) = dat->UnsteadyList[i].Eval_Residual(x, dat->X_old, dat->activity_new, dat->activity_old);
		}
		else
		{
			F(index,0) = dat->UnsteadyList[i].Eval_Residual(x, dat->activity_new);
		}
		index++;
	}
	for (int i=0; i<dat->AdsorptionList.size(); i++)
	{
		success = dat->AdsorptionList[i].callSurfaceActivity(x);
		if (success != 0) {mError(simulation_fail); return -1;}

		dat->AdsorptionList[i].setIonicStrength(x);
		dat->AdsorptionList[i].setChargeDensity(x);

		for (int n=0; n<dat->AdsorptionList[i].getNumberRxns(); n++)
		{
			F(index,0) = dat->AdsorptionList[i].Eval_Residual(x, dat->activity_new, dat->temperature, dat->relative_permittivity, n);
			index++;
		}
	}
	for (int i=0; i<dat->UnsteadyAdsList.size(); i++)
	{
		success = dat->UnsteadyAdsList[i].callSurfaceActivity(x);
		if (success != 0) {mError(simulation_fail); return -1;}

		dat->UnsteadyAdsList[i].setIonicStrength(x);
		dat->UnsteadyAdsList[i].setChargeDensity(x);

		for (int n=0; n<dat->UnsteadyAdsList[i].getNumberRxns(); n++)
		{
			if (dat->steadystate == false)
			{
				if (dat->time == 0.0)
					F(index,0) = dat->UnsteadyAdsList[i].Eval_IC_Residual(x, n);
				else
					F(index,0) = dat->UnsteadyAdsList[i].Eval_Residual(x, dat->X_old, dat->activity_new, dat->activity_old, dat->temperature, dat->relative_permittivity, n);
			}
			else
			{
				F(index,0) = dat->UnsteadyAdsList[i].Eval_Residual(x, dat->activity_new, dat->temperature, dat->relative_permittivity, n);
			}
		index++;
		}
	}
	for (int i=0; i<dat->MultiAdsList.size(); i++)
	{
		success = dat->MultiAdsList[i].callSurfaceActivity(x);
		if (success != 0) {mError(simulation_fail); return -1;}

		dat->MultiAdsList[i].setIonicStrength(x);
		dat->MultiAdsList[i].setChargeDensity(x);

		for (int j=0; j<dat->MultiAdsList[i].getNumberLigands(); j++)
		{
			for (int n=0; n<dat->MultiAdsList[i].getAdsorptionObject(j).getNumberRxns(); n++)
			{
				F(index,0) = dat->MultiAdsList[i].Eval_Residual(x, dat->activity_new, dat->temperature, dat->relative_permittivity, n, j);
				index++;
			}
		}
	}
	for (int i=0; i<dat->ChemisorptionList.size(); i++)
	{
		success = dat->ChemisorptionList[i].callSurfaceActivity(x);
		if (success != 0) {mError(simulation_fail); return -1;}

		dat->ChemisorptionList[i].setIonicStrength(x);
		dat->ChemisorptionList[i].setChargeDensity(x);

		F(index,0) = dat->ChemisorptionList[i].Eval_SiteBalanceResidual(x);
		index++;

		for (int n=0; n<dat->ChemisorptionList[i].getNumberRxns(); n++)
		{
			F(index,0) = dat->ChemisorptionList[i].Eval_RxnResidual(x, dat->activity_new, dat->temperature, dat->relative_permittivity, n);
			index++;
		}
	}
	for (int i=0; i<dat->MultiChemList.size(); i++)
	{
		success = dat->MultiChemList[i].callSurfaceActivity(x);
		if (success != 0) {mError(simulation_fail); return -1;}

		dat->MultiChemList[i].setIonicStrength(x);
		dat->MultiChemList[i].setChargeDensity(x);

		for (int j=0; j<dat->MultiChemList[i].getNumberLigands(); j++)
		{
			F(index,0) = dat->MultiChemList[i].Eval_SiteBalanceResidual(x, j);
			index++;

			for (int n=0; n<dat->MultiChemList[i].getChemisorptionObject(j).getNumberRxns(); n++)
			{
				F(index,0) = dat->MultiChemList[i].Eval_RxnResidual(x, dat->activity_new, dat->temperature, dat->relative_permittivity, n, j);
				index++;
			}
		}
	}
	for (int i=0; i<dat->MassBalanceList.size(); i++)
	{
		if (dat->steadystate == false && dat->time == 0.0)
			F(index,0) = dat->MassBalanceList[i].Eval_IC_Residual(x);
		else
			F(index,0) = dat->MassBalanceList[i].Eval_Residual(x, dat->X_old);
		index++;
	}
	for (int i=0; i<dat->OtherList.size(); i++)
	{
		F(index,0) = dat->OtherList[i] (x, dat, dat->other_data);
		index++;
	}
	if (dat->const_pH == false)
	{
		F(index,0) = dat->MasterList.Eval_ChargeResidual(x);
		if (dat->Contains_pH == true)
			dat->pH = -log10(dat->activity_new(dat->pH_index,0))-x(dat->pH_index,0);
		else
			dat->pH = -1.0;
	}
	else
	{
		if (dat->Contains_pH == true)
			F(index,0) = x(dat->pH_index,0) + log10(dat->activity_new(dat->pH_index,0)) + dat->pH;
		else
		{
			F(index,0) = 0.0;
			mError(invalid_species);
			return -1;
		}
	}

	return success;
}

//Full run of shark simulations
int SHARK(SHARK_DATA *shark_dat)
{
	int success = 0;

	//Function to check for missing information
	success = shark_parameter_check(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Function to establish temperature independent constants
	success = shark_energy_calculations(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Function to calculate equilibrium and rate constants as a function of temperature
	success = shark_temperature_calculations(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Print out the header
	print2file_shark_info(shark_dat);

	//Find the indices of pH and pOH if they exist
	success = shark_pH_finder(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Setup the problem to be solved
	success = shark_initial_conditions(shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Iteratively run the simulation cases
	bool overage = false;
	do
	{
		//Start console messages
		if (shark_dat->Console_Output == true)
		{
			if (shark_dat->const_pH == true)
				std::cout << "----- Performing Simulation @ pH = " << shark_dat->pH << " @ T = " << shark_dat->temperature << " K  ---------\n\n";
			else
				std::cout << "----- Performing Simulation with Electro-Neutrality-Equation @ T = " << shark_dat->temperature << " K  ------\n\n";
		}

		//Call the executioner function
		success = shark_executioner(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

		//Call the reset function
		success = shark_reset(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

		//Conditionally change options
		if (shark_dat->SpeciationCurve == true)
		{
			shark_dat->pH = shark_dat->pH + shark_dat->pH_step;
			if (shark_dat->pH >= shark_dat->pH_end && overage == false)
			{
				overage = true;
				shark_dat->pH = shark_dat->pH_end;
			}
		}
		if (shark_dat->TemperatureCurve == true)
		{
			shark_dat->temperature = shark_dat->temperature + shark_dat->temp_step;
			if (shark_dat->temperature >= shark_dat->end_temp && overage == false)
			{
				overage = true;
				shark_dat->temperature = shark_dat->end_temp;
			}

			//Function to re-calculate equilibrium and rate constants as a function of temperature
			success = shark_temperature_calculations(shark_dat);
			if (success != 0) {mError(simulation_fail); return -1;}
		}

		//End console messages
		if (shark_dat->Console_Output == true)
		{
			if (shark_dat->steadystate == false)
				std::cout << "\n---------------- End of Simulation @ t = " << shark_dat->time << " hours  -----------------\n\n";
			else
				std::cout << "\n---------------- End of Simulation -----------------\n\n";
		}

	} while (	(shark_dat->steadystate == false && shark_dat->simulationtime > (shark_dat->time+DBL_EPSILON))
			    || (shark_dat->SpeciationCurve == true && shark_dat->pH <= shark_dat->pH_end)
				|| (shark_dat->TemperatureCurve == true && shark_dat->temperature <= shark_dat->end_temp) );

	//Call solver one last time to establish the steady-state solution
	if (shark_dat->steadystate == false)
	{
		//Console ouput
		if (shark_dat->Console_Output == true)
			std::cout << "---------------- Establishing the Steady-State Solution ------------------\n\n";

		shark_dat->steadystate = true;
		shark_dat->Newton_data.nl_maxit = 2 * shark_dat->numvar;
		shark_dat->time = -1;
		for (int i=0; i<shark_dat->MassBalanceList.size(); i++)
		{
			shark_dat->MassBalanceList[i].Set_SteadyState(shark_dat->steadystate);
		}

		//Call the executioner function
		success = shark_executioner(shark_dat);
		if (success != 0) {mError(simulation_fail); return -1;}

		shark_dat->steadystate = false;
		shark_dat->time = shark_dat->time_old;

		if (shark_dat->Console_Output == true)
			std::cout << "\n---------------- End of Steady-State Simulation -----------------\n\n";
	}

	return success;
}

//Run a SHARK simulation based on the yaml input file
int SHARK_SCENARIO(const char *yaml_input)
{
	int success = 0;

	//Declarations
	//NOTE: we may have to declare objects for the pitzer and sit models here
	double time;
	SHARK_DATA shark_dat;
	FILE *Output;

	//Initializations
	time = clock();
	Output = fopen("output/SHARK_Output.txt","w+");
	if (Output == nullptr)
	{
		success = system("mkdir output");
		Output = fopen("output/SHARK_Output.txt", "w+");
	}

	//Read the input file
	success = shark_dat.yaml_object.executeYamlRead(yaml_input);
	if (success != 0) {mError(file_dne); return -1;}

	//Read and check Scenario document
	success = read_scenario(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check all Multiligand Adsorbent Objects
	success = read_multiligand_scenario(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check all Multiligand Chemisorbent Objects
	success = read_multichemi_scenario(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//NOTE: we may have to build a custom read option for pitzer and sit models to establish their structures here
	if (shark_dat.act_fun == PITZER || shark_dat.act_fun == SIT || shark_dat.reactor_type == PFR)
	{
		std::cout << "\n--------WARNING!!!---------\n\nUser selected unavailable options! Revert to defaul settings...\n\n";
		if (shark_dat.act_fun == PITZER || shark_dat.act_fun == SIT)
			shark_dat.act_fun = IDEAL;
		if (shark_dat.reactor_type == PFR)
			shark_dat.reactor_type = BATCH;
	}

	//Call the setup function
	//NOTE: The NULL option between &shark_dat and (void *)&shark_dat needs to reflect the data structure for the PITZER model
	success = setup_SHARK_DATA(Output,shark_residual, NULL, NULL, &shark_dat, NULL, (void *)&shark_dat, NULL, NULL);
	if (success != 0) {mError(initial_error); return -1;}

	//Read and check all Multiligand Adsorbent Objects (NOTE: This is a redundant call to fix initialization issues)
	success = read_multiligand_scenario(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check all Multiligand Chemisorbent Objects (NOTE: This is a redundant call to fix initialization issues)
	success = read_multichemi_scenario(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read options
	success = read_options(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check the MasterSpeciesList
	success = read_species(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check the MassBalance
	success = read_massbalance(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check the EquilRxn
	success = read_equilrxn(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check the UnsteadyRxn
	success = read_unsteadyrxn(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check the steady adsorption objects
	success = read_adsorbobjects(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check the unsteady adsorption objects
	success = read_unsteadyadsorbobjects(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check all steady-state adsorbent-ligand paired documents
	success = read_multiligandobjects(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check all steady-state chemisorption objects
	success = read_chemisorbobjects(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Read and check all steady-state chemisorbent-ligand paired documents
	success = read_multichemiobjects(&shark_dat);
	if (success != 0) {mError(read_error); return -1;}

	//Call the SHARK routine
	success = SHARK(&shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Close files and display end messages
	fclose(Output);
	time = clock() - time;
	std::cout << "\nSimulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Time Steps: " << shark_dat.timesteps << "\n";
	std::cout << "Total Iterations: " << shark_dat.totalsteps << "\n";
	std::cout << "Total Function Calls: " << shark_dat.totalcalls << "\n";
	std::cout << "Evaluations/sec: " << shark_dat.totalcalls/(time / CLOCKS_PER_SEC) << "\n";

	return success;
}

//Test of Shark
int SHARK_TESTS()
{
	int success = 0;
	double time;

	SHARK_DATA shark_dat;
	FILE *TestOutput;
	time = clock();

	//Read problem size info here
	TestOutput = fopen("output/SHARK_Test.txt", "w+");
	if (TestOutput == nullptr)
	{
		success = system("mkdir output");
		TestOutput = fopen("output/SHARK_Test.txt", "w+");
	}
	shark_dat.numvar = 25;
	shark_dat.num_ssr = 15;
	shark_dat.num_mbe = 6;

	shark_dat.num_ssao = 0;
	shark_dat.num_usao = 0;
	shark_dat.num_sschem = 1;

	shark_dat.num_usr = 0;
	shark_dat.num_other = 0;
	shark_dat.act_fun = DAVIES;
	shark_dat.steadystate = true;
	shark_dat.simulationtime = 96.0;
	shark_dat.dt = 0.1;
	shark_dat.t_out = shark_dat.simulationtime / 1000.0;
	shark_dat.const_pH = false;
	shark_dat.SpeciationCurve = true;
	shark_dat.TimeAdaptivity = true;
	shark_dat.pH = 7.80;
	shark_dat.dielectric_const = 78.325;
	shark_dat.temperature = 293.15;
	shark_dat.volume = 1.0;

	//Done before setup
	shark_dat.num_sschem_rxns.resize(shark_dat.num_sschem);
	shark_dat.ss_chem_names.resize(shark_dat.num_sschem);
	shark_dat.num_sschem_rxns[0] = 2;

	//Temporary Variables to modify test case
	double NaHCO3 = 1.786E-1;
	double UO2 = 1.079E-6;
	double NaCl = 0.0;
	double NaOH = 0.0;
	double HCl = 0.0;
	double logK_UO2CO3 = -0.3355;
	double logK_UO2 = 4.303;
	double ads_area = 15000.0; // m^2/kg
	double ads_mol = 8.5;      // mol/kg
	double ads_mass = 1.5E-5;  // kg
	double volume = 1.0;       // L

	// ------------------ 1-L Various Carbonate Concentrations -------------------------
	NaOH = 0.0;
	HCl = 0.0;
	NaCl = 0.43;	  // 25.155 g/L
	NaHCO3 = 0.00233; // 140 ppm
	//NaHCO3 = 1e-10; // ~0 ppm
	UO2 = 3.227E-5;   // ~6 ppm
	ads_area = 15000.0; // m^2/kg
	ads_mol = 3.3;     // mol/kg
	ads_mass = 1.5E-5;  // kg
	volume = 0.75;       // L
	shark_dat.simulationtime = 96.0; //hours

	logK_UO2CO3 = -0.92;				// molar basis - slava
	logK_UO2 = -2.72;					// molar basis - slava

	shark_dat.dt = 0.001;			 //hours
	shark_dat.t_out = shark_dat.simulationtime / 1000.0;
	shark_dat.volume = volume;

	//Call the setup function
	success = setup_SHARK_DATA(TestOutput,shark_residual, NULL, NULL, &shark_dat, NULL, (void *)&shark_dat, NULL, NULL);
	if (success != 0) {mError(simulation_fail); return -1;}

	shark_dat.Newton_data.linear_solver = QR;
	shark_dat.Newton_data.LineSearch = true;
	shark_dat.Newton_data.nl_maxit = 100;
	shark_dat.Newton_data.NL_Output = true;

	//Read problem specific info here --------------------------------------------------------
	shark_dat.MasterList.set_species(0, "NaHCO3 (aq)");
	shark_dat.MasterList.set_species(1, "NaCO3 - (aq)");
	shark_dat.MasterList.set_species(2, "Na + (aq)");
	shark_dat.MasterList.set_species(3, "HNO3 (aq)");
	shark_dat.MasterList.set_species(4, "NO3 - (aq)");
	shark_dat.MasterList.set_species(5, "H2CO3 (aq)");
	shark_dat.MasterList.set_species(6, "HCO3 - (aq)");
	shark_dat.MasterList.set_species(7, "CO3 2- (aq)");
	shark_dat.MasterList.set_species(8, "UO2 2+ (aq)");
	shark_dat.MasterList.set_species(9, "UO2NO3 + (aq)");
	shark_dat.MasterList.set_species(10, "UO2(NO3)2 (aq)");
	shark_dat.MasterList.set_species(11, "UO2OH + (aq)");
	shark_dat.MasterList.set_species(12, "UO2(OH)3 - (aq)");
	shark_dat.MasterList.set_species(13, "(UO2)2(OH)2 2+ (aq)");
	shark_dat.MasterList.set_species(14, "(UO2)3(OH)5 + (aq)");
	shark_dat.MasterList.set_species(15, "UO2CO3 (aq)");
	shark_dat.MasterList.set_species(16, "UO2(CO3)2 2- (aq)");
	shark_dat.MasterList.set_species(17, "UO2(CO3)3 4- (aq)");
	shark_dat.MasterList.set_species(18, "H2O (l)");
	shark_dat.MasterList.set_species(19, "OH - (aq)");
	shark_dat.MasterList.set_species(20, "H + (aq)");
	shark_dat.MasterList.set_species(21, "Cl - (aq)");

	shark_dat.MasterList.set_species(22, 0, 0, 0, 0, false, false, "Adsorbed", "Uranyl-amidoxime", "UO2(AO)2 (ad)", "UO2");
	shark_dat.MasterList.set_species(23, -2, 0, 0, 0, false, false, "Adsorbed", "Uranyl-carbonate-amidoxime", "UO2CO3(AO)2 2- (ad)", "UO2CO3");
	shark_dat.MasterList.set_species(24, 0, 0, 0, 0, false, false, "Adsorbed", "amidoxime", "HAO (ad)", "H");

	shark_dat.ReactionList[0].Set_Equilibrium(-14.0);
	shark_dat.ReactionList[0].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(18, -1);
	shark_dat.ReactionList[0].Set_Stoichiometric(19, 1);
	shark_dat.ReactionList[0].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[0].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[1].Set_Equilibrium(-6.35);
	shark_dat.ReactionList[1].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(5, -1);
	shark_dat.ReactionList[1].Set_Stoichiometric(6, 1);
	shark_dat.ReactionList[1].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[1].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[2].Set_Equilibrium(-10.33);
	shark_dat.ReactionList[2].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(6, -1);
	shark_dat.ReactionList[2].Set_Stoichiometric(7, 1);
	shark_dat.ReactionList[2].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[2].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[3].Set_Equilibrium(-10.14);
	shark_dat.ReactionList[3].Set_Stoichiometric(0, -1);
	shark_dat.ReactionList[3].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(2, 1);
	shark_dat.ReactionList[3].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(7, 1);
	shark_dat.ReactionList[3].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[3].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[4].Set_Equilibrium(-1.02);
	shark_dat.ReactionList[4].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(1, -1);
	shark_dat.ReactionList[4].Set_Stoichiometric(2, 1);
	shark_dat.ReactionList[4].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(7, 1);
	shark_dat.ReactionList[4].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[5].Set_Equilibrium(1.4);
	shark_dat.ReactionList[5].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(3, -1);
	shark_dat.ReactionList[5].Set_Stoichiometric(4, 1);
	shark_dat.ReactionList[5].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[5].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[6].Set_Equilibrium(-0.3);
	shark_dat.ReactionList[6].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(4, -1);
	shark_dat.ReactionList[6].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[6].Set_Stoichiometric(9, 1);
	shark_dat.ReactionList[6].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[7].Set_Equilibrium(-12.15);
	shark_dat.ReactionList[7].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(4, -2);
	shark_dat.ReactionList[7].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[7].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(10, 1);
	shark_dat.ReactionList[7].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[8].Set_Equilibrium(-6.2);
	shark_dat.ReactionList[8].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[8].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(11, 1);
	shark_dat.ReactionList[8].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(18, -1);
	shark_dat.ReactionList[8].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[8].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[9].Set_Equilibrium(-20.2);
	shark_dat.ReactionList[9].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[9].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(12, 1);
	shark_dat.ReactionList[9].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(18, -3);
	shark_dat.ReactionList[9].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(20, 3);
	shark_dat.ReactionList[9].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[10].Set_Equilibrium(-5.87);
	shark_dat.ReactionList[10].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(8, -2);
	shark_dat.ReactionList[10].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(13, 1);
	shark_dat.ReactionList[10].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(18, -2);
	shark_dat.ReactionList[10].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(20, 2);
	shark_dat.ReactionList[10].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[11].Set_Equilibrium(-16.5);
	shark_dat.ReactionList[11].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(8, -3);
	shark_dat.ReactionList[11].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(14, 1);
	shark_dat.ReactionList[11].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(18, -5);
	shark_dat.ReactionList[11].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(20, 5);
	shark_dat.ReactionList[11].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[12].Set_Equilibrium(8.4);
	shark_dat.ReactionList[12].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(7, -1);
	shark_dat.ReactionList[12].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[12].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(15, 1);
	shark_dat.ReactionList[12].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[13].Set_Equilibrium(15.7);
	shark_dat.ReactionList[13].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(7, -2);
	shark_dat.ReactionList[13].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[13].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(16, 1);
	shark_dat.ReactionList[13].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[14].Set_Equilibrium(21.6);
	shark_dat.ReactionList[14].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(7, -3);
	shark_dat.ReactionList[14].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[14].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(17, 1);
	shark_dat.ReactionList[14].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(23, 0);

	shark_dat.MassBalanceList[0].Set_Name("Water Balance");
	shark_dat.MassBalanceList[0].Set_TotalConcentration(1.0);
	shark_dat.MassBalanceList[0].Set_Delta(0, 0);
	shark_dat.MassBalanceList[0].Set_Delta(1, 0);
	shark_dat.MassBalanceList[0].Set_Delta(2, 0);
	shark_dat.MassBalanceList[0].Set_Delta(3, 0);
	shark_dat.MassBalanceList[0].Set_Delta(4, 0);
	shark_dat.MassBalanceList[0].Set_Delta(5, 0);
	shark_dat.MassBalanceList[0].Set_Delta(6, 0);
	shark_dat.MassBalanceList[0].Set_Delta(7, 0);
	shark_dat.MassBalanceList[0].Set_Delta(8, 0);
	shark_dat.MassBalanceList[0].Set_Delta(9, 0);
	shark_dat.MassBalanceList[0].Set_Delta(10, 0);
	shark_dat.MassBalanceList[0].Set_Delta(11, 0);
	shark_dat.MassBalanceList[0].Set_Delta(12, 0);
	shark_dat.MassBalanceList[0].Set_Delta(13, 0);
	shark_dat.MassBalanceList[0].Set_Delta(14, 0);
	shark_dat.MassBalanceList[0].Set_Delta(15, 0);
	shark_dat.MassBalanceList[0].Set_Delta(16, 0);
	shark_dat.MassBalanceList[0].Set_Delta(17, 0);
	shark_dat.MassBalanceList[0].Set_Delta(18, 1);
	shark_dat.MassBalanceList[0].Set_Delta(19, 0);
	shark_dat.MassBalanceList[0].Set_Delta(20, 0);
	shark_dat.MassBalanceList[0].Set_Delta(21, 0);
	shark_dat.MassBalanceList[0].Set_Delta(22, 0);
	shark_dat.MassBalanceList[0].Set_Delta(23, 0);

	shark_dat.MassBalanceList[1].Set_Name("Carbonate Balance");
	shark_dat.MassBalanceList[1].Set_TotalConcentration(NaHCO3);
	shark_dat.MassBalanceList[1].Set_Delta(0, 1);
	shark_dat.MassBalanceList[1].Set_Delta(1, 1);
	shark_dat.MassBalanceList[1].Set_Delta(2, 0);
	shark_dat.MassBalanceList[1].Set_Delta(3, 0);
	shark_dat.MassBalanceList[1].Set_Delta(4, 0);
	shark_dat.MassBalanceList[1].Set_Delta(5, 1);
	shark_dat.MassBalanceList[1].Set_Delta(6, 1);
	shark_dat.MassBalanceList[1].Set_Delta(7, 1);
	shark_dat.MassBalanceList[1].Set_Delta(8, 0);
	shark_dat.MassBalanceList[1].Set_Delta(9, 0);
	shark_dat.MassBalanceList[1].Set_Delta(10, 0);
	shark_dat.MassBalanceList[1].Set_Delta(11, 0);
	shark_dat.MassBalanceList[1].Set_Delta(12, 0);
	shark_dat.MassBalanceList[1].Set_Delta(13, 0);
	shark_dat.MassBalanceList[1].Set_Delta(14, 0);
	shark_dat.MassBalanceList[1].Set_Delta(15, 1);
	shark_dat.MassBalanceList[1].Set_Delta(16, 2);
	shark_dat.MassBalanceList[1].Set_Delta(17, 3);
	shark_dat.MassBalanceList[1].Set_Delta(18, 0);
	shark_dat.MassBalanceList[1].Set_Delta(19, 0);
	shark_dat.MassBalanceList[1].Set_Delta(20, 0);
	shark_dat.MassBalanceList[1].Set_Delta(21, 0);
	shark_dat.MassBalanceList[1].Set_Delta(22, 0);
	shark_dat.MassBalanceList[1].Set_Delta(23, 1);

	shark_dat.MassBalanceList[2].Set_Name("Nitrate Balance");
	shark_dat.MassBalanceList[2].Set_TotalConcentration(2.0 * UO2);
	shark_dat.MassBalanceList[2].Set_Delta(0, 0);
	shark_dat.MassBalanceList[2].Set_Delta(1, 0);
	shark_dat.MassBalanceList[2].Set_Delta(2, 0);
	shark_dat.MassBalanceList[2].Set_Delta(3, 1);
	shark_dat.MassBalanceList[2].Set_Delta(4, 1);
	shark_dat.MassBalanceList[2].Set_Delta(5, 0);
	shark_dat.MassBalanceList[2].Set_Delta(6, 0);
	shark_dat.MassBalanceList[2].Set_Delta(7, 0);
	shark_dat.MassBalanceList[2].Set_Delta(8, 0);
	shark_dat.MassBalanceList[2].Set_Delta(9, 1);
	shark_dat.MassBalanceList[2].Set_Delta(10, 2);
	shark_dat.MassBalanceList[2].Set_Delta(11, 0);
	shark_dat.MassBalanceList[2].Set_Delta(12, 0);
	shark_dat.MassBalanceList[2].Set_Delta(13, 0);
	shark_dat.MassBalanceList[2].Set_Delta(14, 0);
	shark_dat.MassBalanceList[2].Set_Delta(15, 0);
	shark_dat.MassBalanceList[2].Set_Delta(16, 0);
	shark_dat.MassBalanceList[2].Set_Delta(17, 0);
	shark_dat.MassBalanceList[2].Set_Delta(18, 0);
	shark_dat.MassBalanceList[2].Set_Delta(19, 0);
	shark_dat.MassBalanceList[2].Set_Delta(20, 0);
	shark_dat.MassBalanceList[2].Set_Delta(21, 0);
	shark_dat.MassBalanceList[2].Set_Delta(22, 0);
	shark_dat.MassBalanceList[2].Set_Delta(23, 0);

	shark_dat.MassBalanceList[3].Set_Name("Sodium Balance");
	shark_dat.MassBalanceList[3].Set_TotalConcentration(NaHCO3+NaCl+NaOH);
	shark_dat.MassBalanceList[3].Set_Delta(0, 1);
	shark_dat.MassBalanceList[3].Set_Delta(1, 1);
	shark_dat.MassBalanceList[3].Set_Delta(2, 1);
	shark_dat.MassBalanceList[3].Set_Delta(3, 0);
	shark_dat.MassBalanceList[3].Set_Delta(4, 0);
	shark_dat.MassBalanceList[3].Set_Delta(5, 0);
	shark_dat.MassBalanceList[3].Set_Delta(6, 0);
	shark_dat.MassBalanceList[3].Set_Delta(7, 0);
	shark_dat.MassBalanceList[3].Set_Delta(8, 0);
	shark_dat.MassBalanceList[3].Set_Delta(9, 0);
	shark_dat.MassBalanceList[3].Set_Delta(10, 0);
	shark_dat.MassBalanceList[3].Set_Delta(11, 0);
	shark_dat.MassBalanceList[3].Set_Delta(12, 0);
	shark_dat.MassBalanceList[3].Set_Delta(13, 0);
	shark_dat.MassBalanceList[3].Set_Delta(14, 0);
	shark_dat.MassBalanceList[3].Set_Delta(15, 0);
	shark_dat.MassBalanceList[3].Set_Delta(16, 0);
	shark_dat.MassBalanceList[3].Set_Delta(17, 0);
	shark_dat.MassBalanceList[3].Set_Delta(18, 0);
	shark_dat.MassBalanceList[3].Set_Delta(19, 0);
	shark_dat.MassBalanceList[3].Set_Delta(20, 0);
	shark_dat.MassBalanceList[3].Set_Delta(21, 0);
	shark_dat.MassBalanceList[3].Set_Delta(22, 0);
	shark_dat.MassBalanceList[3].Set_Delta(23, 0);

	shark_dat.MassBalanceList[4].Set_Name("Uranium Balance");
	shark_dat.MassBalanceList[4].Set_TotalConcentration(UO2);
	shark_dat.MassBalanceList[4].Set_Delta(0, 0);
	shark_dat.MassBalanceList[4].Set_Delta(1, 0);
	shark_dat.MassBalanceList[4].Set_Delta(2, 0);
	shark_dat.MassBalanceList[4].Set_Delta(3, 0);
	shark_dat.MassBalanceList[4].Set_Delta(4, 0);
	shark_dat.MassBalanceList[4].Set_Delta(5, 0);
	shark_dat.MassBalanceList[4].Set_Delta(6, 0);
	shark_dat.MassBalanceList[4].Set_Delta(7, 0);
	shark_dat.MassBalanceList[4].Set_Delta(8, 1);
	shark_dat.MassBalanceList[4].Set_Delta(9, 1);
	shark_dat.MassBalanceList[4].Set_Delta(10, 1);
	shark_dat.MassBalanceList[4].Set_Delta(11, 1);
	shark_dat.MassBalanceList[4].Set_Delta(12, 1);
	shark_dat.MassBalanceList[4].Set_Delta(13, 2);
	shark_dat.MassBalanceList[4].Set_Delta(14, 3);
	shark_dat.MassBalanceList[4].Set_Delta(15, 1);
	shark_dat.MassBalanceList[4].Set_Delta(16, 1);
	shark_dat.MassBalanceList[4].Set_Delta(17, 1);
	shark_dat.MassBalanceList[4].Set_Delta(18, 0);
	shark_dat.MassBalanceList[4].Set_Delta(19, 0);
	shark_dat.MassBalanceList[4].Set_Delta(20, 0);
	shark_dat.MassBalanceList[4].Set_Delta(21, 0);
	shark_dat.MassBalanceList[4].Set_Delta(22, 1);
	shark_dat.MassBalanceList[4].Set_Delta(23, 1);

	shark_dat.MassBalanceList[5].Set_Name("Chlorine Balance");
	shark_dat.MassBalanceList[5].Set_TotalConcentration(NaCl+HCl);
	shark_dat.MassBalanceList[5].Set_Delta(0, 0);
	shark_dat.MassBalanceList[5].Set_Delta(1, 0);
	shark_dat.MassBalanceList[5].Set_Delta(2, 0);
	shark_dat.MassBalanceList[5].Set_Delta(3, 0);
	shark_dat.MassBalanceList[5].Set_Delta(4, 0);
	shark_dat.MassBalanceList[5].Set_Delta(5, 0);
	shark_dat.MassBalanceList[5].Set_Delta(6, 0);
	shark_dat.MassBalanceList[5].Set_Delta(7, 0);
	shark_dat.MassBalanceList[5].Set_Delta(8, 0);
	shark_dat.MassBalanceList[5].Set_Delta(9, 0);
	shark_dat.MassBalanceList[5].Set_Delta(10, 0);
	shark_dat.MassBalanceList[5].Set_Delta(11, 0);
	shark_dat.MassBalanceList[5].Set_Delta(12, 0);
	shark_dat.MassBalanceList[5].Set_Delta(13, 0);
	shark_dat.MassBalanceList[5].Set_Delta(14, 0);
	shark_dat.MassBalanceList[5].Set_Delta(15, 0);
	shark_dat.MassBalanceList[5].Set_Delta(16, 0);
	shark_dat.MassBalanceList[5].Set_Delta(17, 0);
	shark_dat.MassBalanceList[5].Set_Delta(18, 0);
	shark_dat.MassBalanceList[5].Set_Delta(19, 0);
	shark_dat.MassBalanceList[5].Set_Delta(20, 0);
	shark_dat.MassBalanceList[5].Set_Delta(21, 1);
	shark_dat.MassBalanceList[5].Set_Delta(22, 0);
	shark_dat.MassBalanceList[5].Set_Delta(23, 0);


	shark_dat.ChemisorptionList[0].setSpecificArea(ads_area);
	shark_dat.ChemisorptionList[0].setSpecificMolality(ads_mol);
	shark_dat.ChemisorptionList[0].setTotalMass(ads_mass);
	shark_dat.ChemisorptionList[0].setAdsorbentName("HAO (ad)");
	shark_dat.ChemisorptionList[0].setActivityModelInfo(UNIQUAC_chemi, &shark_dat.ChemisorptionList[0]);
	shark_dat.ChemisorptionList[0].setActivityEnum(UNIQUAC_ACT);

	shark_dat.ChemisorptionList[0].getReaction(0).Set_Equilibrium(logK_UO2);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(0, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(1, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(2, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(3, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(4, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(5, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(6, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(7, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(8, -1);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(9, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(10, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(11, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(12, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(13, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(14, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(15, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(16, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(17, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(18, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(19, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(20, 2);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(21, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(22, 1);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(23, 0);
	shark_dat.ChemisorptionList[0].getReaction(0).Set_Stoichiometric(24, -2);


	shark_dat.ChemisorptionList[0].getReaction(1).Set_Equilibrium(logK_UO2CO3);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(0, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(1, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(2, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(3, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(4, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(5, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(6, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(7, -1);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(8, -1);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(9, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(10, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(11, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(12, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(13, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(14, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(15, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(16, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(17, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(18, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(19, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(20, 2);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(21, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(22, 0);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(23, 1);
	shark_dat.ChemisorptionList[0].getReaction(1).Set_Stoichiometric(24, -2);

	shark_dat.ChemisorptionList[0].setVolumeFactor(22, 5);
	shark_dat.ChemisorptionList[0].setVolumeFactor(23, 8);
	shark_dat.ChemisorptionList[0].setVolumeFactor(24, 1);

	//shark_dat.ChemisorptionList[0].setDelta(22, 2);			//Done in the input file and linked to stoicheometry
	//shark_dat.ChemisorptionList[0].setDelta(23, 2);

	shark_dat.ChemisorptionList[0].setLigandIndex();			//NOTE: Required here because this is done during the reading step
	shark_dat.ChemisorptionList[0].setAdsorbIndices();
	shark_dat.ChemisorptionList[0].setDeltas();
	shark_dat.ChemisorptionList[0].calculateAreaFactors();    //NOTE: Required here because this is done during the reading step

	/*
	shark_dat.ChemisorptionList[0].Display_Info();

	//Function to check for missing information
	success = shark_parameter_check(&shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Function to establish temperature independent constants
	success = shark_energy_calculations(&shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Function to calculate equilibrium and rate constants as a function of temperature
	success = shark_temperature_calculations(&shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Print out the header
	print2file_shark_info(&shark_dat);
	*/

	// END problem specific info here --------------------------------------------------------

	//Call the SHARK routine
	success = SHARK(&shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Close files and display end messages
	fclose(TestOutput);
	time = clock() - time;
	std::cout << "\nSimulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Time Steps: " << shark_dat.timesteps << "\n";
	std::cout << "Total Iterations: " << shark_dat.totalsteps << "\n";
	std::cout << "Total Function Calls: " << shark_dat.totalcalls << "\n";
	std::cout << "Evaluations/sec: " << shark_dat.totalcalls/(time / CLOCKS_PER_SEC) << "\n";

	return success;
}


//Test of Shark (old version)
int SHARK_TESTS_OLD()
{
	int success = 0;
	double time;

	SHARK_DATA shark_dat;
	FILE *TestOutput;
	time = clock();

	//Read problem size info here
	TestOutput = fopen("output/SHARK_Test.txt", "w+");
	if (TestOutput == nullptr)
	{
		success = system("mkdir output");
		TestOutput = fopen("output/SHARK_Test.txt", "w+");
	}
	shark_dat.numvar = 24;
	shark_dat.num_ssr = 15;
	shark_dat.num_mbe = 6;

	shark_dat.num_ssao = 0;
	shark_dat.num_usao = 1;

	shark_dat.num_usr = 0;
	shark_dat.num_other = 0;
	shark_dat.act_fun = DAVIES;
	shark_dat.steadystate = true;
	shark_dat.simulationtime = 96.0;
	shark_dat.dt = 0.1;
	shark_dat.t_out = shark_dat.simulationtime / 1000.0;
	shark_dat.const_pH = false;
	shark_dat.SpeciationCurve = true;
	shark_dat.TimeAdaptivity = false;
	//shark_dat.reactor_type = CSTR;
	shark_dat.pH = 7.80;
	shark_dat.dielectric_const = 78.325;
	shark_dat.temperature = 293.15;
	shark_dat.volume = 1.0;							//SHOULD BE REQUIRED IF WE HAVE ADSORPTION OBJECTS!!!

	shark_dat.num_usar.resize(shark_dat.num_usao);	//Required to set this up PRIOR to calling the setup function for shark
	shark_dat.us_ads_names.resize(shark_dat.num_usao);
	shark_dat.num_usar[0] = 2;						//Required to set this up PRIOR to calling the setup function for shark


	//Temporary Variables to modify test case
	double NaHCO3 = 1.786E-1;
	double UO2 = 1.079E-6;
	double NaCl = 0.0;
	double NaOH = 0.0;
	double HCl = 0.0;
	double logK_UO2CO3 = -0.3355;
	double logK_UO2 = 4.303;
	double ads_area = 15000.0; // m^2/kg
	double ads_mol = 8.5;      // mol/kg
	double ads_mass = 1.5E-5;  // kg
	double volume = 1.0;       // L

	// ------------------ 1-L Various Carbonate Concentrations -------------------------
	NaOH = 0.0;
	HCl = 0.0;
	NaCl = 0.43;	  // 25.155 g/L
	NaHCO3 = 0.00233; // 140 ppm
	//NaHCO3 = 1e-10; // ~0 ppm
	UO2 = 3.227E-5;   // ~6 ppm
	ads_area = 15000.0; // m^2/kg
	ads_mol = 3.3;     // mol/kg
	ads_mass = 1.5E-5;  // kg
	volume = 0.75;       // L
	shark_dat.simulationtime = 96.0; //hours

	logK_UO2CO3 = -0.92;				// molar basis - slava
	logK_UO2 = -2.72;					// molar basis - slava

	shark_dat.dt = 0.001;			 //hours
	shark_dat.t_out = shark_dat.simulationtime / 1000.0;
	shark_dat.volume = volume;

	//Call the setup function
	success = setup_SHARK_DATA(TestOutput,shark_residual, NULL, NULL, &shark_dat, NULL, (void *)&shark_dat, NULL, NULL);
	if (success != 0) {mError(simulation_fail); return -1;}

	shark_dat.Newton_data.linear_solver = QR;
	shark_dat.Newton_data.LineSearch = true;
	shark_dat.Newton_data.nl_maxit = 100;
	shark_dat.Newton_data.NL_Output = true;

	//Read problem specific info here --------------------------------------------------------
	shark_dat.MasterList.set_species(0, "NaHCO3 (aq)");
	shark_dat.MasterList.set_species(1, "NaCO3 - (aq)");
	shark_dat.MasterList.set_species(2, "Na + (aq)");
	shark_dat.MasterList.set_species(3, "HNO3 (aq)");
	shark_dat.MasterList.set_species(4, "NO3 - (aq)");
	shark_dat.MasterList.set_species(5, "H2CO3 (aq)");
	shark_dat.MasterList.set_species(6, "HCO3 - (aq)");
	shark_dat.MasterList.set_species(7, "CO3 2- (aq)");
	shark_dat.MasterList.set_species(8, "UO2 2+ (aq)");
	shark_dat.MasterList.set_species(9, "UO2NO3 + (aq)");
	shark_dat.MasterList.set_species(10, "UO2(NO3)2 (aq)");
	shark_dat.MasterList.set_species(11, "UO2OH + (aq)");
	shark_dat.MasterList.set_species(12, "UO2(OH)3 - (aq)");
	shark_dat.MasterList.set_species(13, "(UO2)2(OH)2 2+ (aq)");
	shark_dat.MasterList.set_species(14, "(UO2)3(OH)5 + (aq)");
	shark_dat.MasterList.set_species(15, "UO2CO3 (aq)");
	shark_dat.MasterList.set_species(16, "UO2(CO3)2 2- (aq)");
	shark_dat.MasterList.set_species(17, "UO2(CO3)3 4- (aq)");
	shark_dat.MasterList.set_species(18, "H2O (l)");
	shark_dat.MasterList.set_species(19, "OH - (aq)");
	shark_dat.MasterList.set_species(20, "H + (aq)");
	shark_dat.MasterList.set_species(21, "Cl - (aq)");

	shark_dat.MasterList.set_species(22, 0, 0, 0, 0, false, false, "Adsorbed", "Uranyl-amidoxime", "UO2(AO)2 (ad)", "UO2");
	shark_dat.MasterList.set_species(23, -2, 0, 0, 0, false, false, "Adsorbed", "Uranyl-carbonate-amidoxime", "UO2CO3(AO)2 2- (ad)", "UO2CO3");

	shark_dat.ReactionList[0].Set_Equilibrium(-14.0);
	shark_dat.ReactionList[0].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(18, -1);
	shark_dat.ReactionList[0].Set_Stoichiometric(19, 1);
	shark_dat.ReactionList[0].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[0].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[0].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[1].Set_Equilibrium(-6.35);
	shark_dat.ReactionList[1].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(5, -1);
	shark_dat.ReactionList[1].Set_Stoichiometric(6, 1);
	shark_dat.ReactionList[1].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[1].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[1].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[2].Set_Equilibrium(-10.33);
	shark_dat.ReactionList[2].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(6, -1);
	shark_dat.ReactionList[2].Set_Stoichiometric(7, 1);
	shark_dat.ReactionList[2].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[2].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[2].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[3].Set_Equilibrium(-10.14);
	shark_dat.ReactionList[3].Set_Stoichiometric(0, -1);
	shark_dat.ReactionList[3].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(2, 1);
	shark_dat.ReactionList[3].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(7, 1);
	shark_dat.ReactionList[3].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[3].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[3].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[4].Set_Equilibrium(-1.02);
	shark_dat.ReactionList[4].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(1, -1);
	shark_dat.ReactionList[4].Set_Stoichiometric(2, 1);
	shark_dat.ReactionList[4].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(7, 1);
	shark_dat.ReactionList[4].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[4].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[5].Set_Equilibrium(1.4);
	shark_dat.ReactionList[5].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(3, -1);
	shark_dat.ReactionList[5].Set_Stoichiometric(4, 1);
	shark_dat.ReactionList[5].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(8, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[5].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[5].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[6].Set_Equilibrium(-0.3);
	shark_dat.ReactionList[6].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(4, -1);
	shark_dat.ReactionList[6].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[6].Set_Stoichiometric(9, 1);
	shark_dat.ReactionList[6].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[6].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[7].Set_Equilibrium(-12.15);
	shark_dat.ReactionList[7].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(4, -2);
	shark_dat.ReactionList[7].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[7].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(10, 1);
	shark_dat.ReactionList[7].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[7].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[8].Set_Equilibrium(-6.2);
	shark_dat.ReactionList[8].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[8].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(11, 1);
	shark_dat.ReactionList[8].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(18, -1);
	shark_dat.ReactionList[8].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(20, 1);
	shark_dat.ReactionList[8].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[8].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[9].Set_Equilibrium(-20.2);
	shark_dat.ReactionList[9].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[9].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(12, 1);
	shark_dat.ReactionList[9].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(18, -3);
	shark_dat.ReactionList[9].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(20, 3);
	shark_dat.ReactionList[9].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[9].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[10].Set_Equilibrium(-5.87);
	shark_dat.ReactionList[10].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(8, -2);
	shark_dat.ReactionList[10].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(13, 1);
	shark_dat.ReactionList[10].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(18, -2);
	shark_dat.ReactionList[10].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(20, 2);
	shark_dat.ReactionList[10].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[10].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[11].Set_Equilibrium(-16.5);
	shark_dat.ReactionList[11].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(7, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(8, -3);
	shark_dat.ReactionList[11].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(14, 1);
	shark_dat.ReactionList[11].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(18, -5);
	shark_dat.ReactionList[11].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(20, 5);
	shark_dat.ReactionList[11].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[11].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[12].Set_Equilibrium(8.4);
	shark_dat.ReactionList[12].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(7, -1);
	shark_dat.ReactionList[12].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[12].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(15, 1);
	shark_dat.ReactionList[12].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[12].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[13].Set_Equilibrium(15.7);
	shark_dat.ReactionList[13].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(7, -2);
	shark_dat.ReactionList[13].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[13].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(16, 1);
	shark_dat.ReactionList[13].Set_Stoichiometric(17, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[13].Set_Stoichiometric(23, 0);

	shark_dat.ReactionList[14].Set_Equilibrium(21.6);
	shark_dat.ReactionList[14].Set_Stoichiometric(0, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(1, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(2, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(3, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(4, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(5, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(6, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(7, -3);
	shark_dat.ReactionList[14].Set_Stoichiometric(8, -1);
	shark_dat.ReactionList[14].Set_Stoichiometric(9, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(10, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(11, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(12, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(13, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(14, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(15, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(16, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(17, 1);
	shark_dat.ReactionList[14].Set_Stoichiometric(18, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(19, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(20, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(21, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(22, 0);
	shark_dat.ReactionList[14].Set_Stoichiometric(23, 0);

	shark_dat.MassBalanceList[0].Set_Name("Water Balance");
	shark_dat.MassBalanceList[0].Set_TotalConcentration(1.0);
	shark_dat.MassBalanceList[0].Set_Delta(0, 0);
	shark_dat.MassBalanceList[0].Set_Delta(1, 0);
	shark_dat.MassBalanceList[0].Set_Delta(2, 0);
	shark_dat.MassBalanceList[0].Set_Delta(3, 0);
	shark_dat.MassBalanceList[0].Set_Delta(4, 0);
	shark_dat.MassBalanceList[0].Set_Delta(5, 0);
	shark_dat.MassBalanceList[0].Set_Delta(6, 0);
	shark_dat.MassBalanceList[0].Set_Delta(7, 0);
	shark_dat.MassBalanceList[0].Set_Delta(8, 0);
	shark_dat.MassBalanceList[0].Set_Delta(9, 0);
	shark_dat.MassBalanceList[0].Set_Delta(10, 0);
	shark_dat.MassBalanceList[0].Set_Delta(11, 0);
	shark_dat.MassBalanceList[0].Set_Delta(12, 0);
	shark_dat.MassBalanceList[0].Set_Delta(13, 0);
	shark_dat.MassBalanceList[0].Set_Delta(14, 0);
	shark_dat.MassBalanceList[0].Set_Delta(15, 0);
	shark_dat.MassBalanceList[0].Set_Delta(16, 0);
	shark_dat.MassBalanceList[0].Set_Delta(17, 0);
	shark_dat.MassBalanceList[0].Set_Delta(18, 1);
	shark_dat.MassBalanceList[0].Set_Delta(19, 0);
	shark_dat.MassBalanceList[0].Set_Delta(20, 0);
	shark_dat.MassBalanceList[0].Set_Delta(21, 0);
	shark_dat.MassBalanceList[0].Set_Delta(22, 0);
	shark_dat.MassBalanceList[0].Set_Delta(23, 0);

	shark_dat.MassBalanceList[1].Set_Name("Carbonate Balance");
	shark_dat.MassBalanceList[1].Set_TotalConcentration(NaHCO3);
	shark_dat.MassBalanceList[1].Set_Delta(0, 1);
	shark_dat.MassBalanceList[1].Set_Delta(1, 1);
	shark_dat.MassBalanceList[1].Set_Delta(2, 0);
	shark_dat.MassBalanceList[1].Set_Delta(3, 0);
	shark_dat.MassBalanceList[1].Set_Delta(4, 0);
	shark_dat.MassBalanceList[1].Set_Delta(5, 1);
	shark_dat.MassBalanceList[1].Set_Delta(6, 1);
	shark_dat.MassBalanceList[1].Set_Delta(7, 1);
	shark_dat.MassBalanceList[1].Set_Delta(8, 0);
	shark_dat.MassBalanceList[1].Set_Delta(9, 0);
	shark_dat.MassBalanceList[1].Set_Delta(10, 0);
	shark_dat.MassBalanceList[1].Set_Delta(11, 0);
	shark_dat.MassBalanceList[1].Set_Delta(12, 0);
	shark_dat.MassBalanceList[1].Set_Delta(13, 0);
	shark_dat.MassBalanceList[1].Set_Delta(14, 0);
	shark_dat.MassBalanceList[1].Set_Delta(15, 1);
	shark_dat.MassBalanceList[1].Set_Delta(16, 2);
	shark_dat.MassBalanceList[1].Set_Delta(17, 3);
	shark_dat.MassBalanceList[1].Set_Delta(18, 0);
	shark_dat.MassBalanceList[1].Set_Delta(19, 0);
	shark_dat.MassBalanceList[1].Set_Delta(20, 0);
	shark_dat.MassBalanceList[1].Set_Delta(21, 0);
	shark_dat.MassBalanceList[1].Set_Delta(22, 0);
	shark_dat.MassBalanceList[1].Set_Delta(23, 1);

	shark_dat.MassBalanceList[2].Set_Name("Nitrate Balance");
	shark_dat.MassBalanceList[2].Set_TotalConcentration(2.0 * UO2);
	shark_dat.MassBalanceList[2].Set_Delta(0, 0);
	shark_dat.MassBalanceList[2].Set_Delta(1, 0);
	shark_dat.MassBalanceList[2].Set_Delta(2, 0);
	shark_dat.MassBalanceList[2].Set_Delta(3, 1);
	shark_dat.MassBalanceList[2].Set_Delta(4, 1);
	shark_dat.MassBalanceList[2].Set_Delta(5, 0);
	shark_dat.MassBalanceList[2].Set_Delta(6, 0);
	shark_dat.MassBalanceList[2].Set_Delta(7, 0);
	shark_dat.MassBalanceList[2].Set_Delta(8, 0);
	shark_dat.MassBalanceList[2].Set_Delta(9, 1);
	shark_dat.MassBalanceList[2].Set_Delta(10, 2);
	shark_dat.MassBalanceList[2].Set_Delta(11, 0);
	shark_dat.MassBalanceList[2].Set_Delta(12, 0);
	shark_dat.MassBalanceList[2].Set_Delta(13, 0);
	shark_dat.MassBalanceList[2].Set_Delta(14, 0);
	shark_dat.MassBalanceList[2].Set_Delta(15, 0);
	shark_dat.MassBalanceList[2].Set_Delta(16, 0);
	shark_dat.MassBalanceList[2].Set_Delta(17, 0);
	shark_dat.MassBalanceList[2].Set_Delta(18, 0);
	shark_dat.MassBalanceList[2].Set_Delta(19, 0);
	shark_dat.MassBalanceList[2].Set_Delta(20, 0);
	shark_dat.MassBalanceList[2].Set_Delta(21, 0);
	shark_dat.MassBalanceList[2].Set_Delta(22, 0);
	shark_dat.MassBalanceList[2].Set_Delta(23, 0);

	shark_dat.MassBalanceList[3].Set_Name("Sodium Balance");
	shark_dat.MassBalanceList[3].Set_TotalConcentration(NaHCO3+NaCl+NaOH);
	shark_dat.MassBalanceList[3].Set_Delta(0, 1);
	shark_dat.MassBalanceList[3].Set_Delta(1, 1);
	shark_dat.MassBalanceList[3].Set_Delta(2, 1);
	shark_dat.MassBalanceList[3].Set_Delta(3, 0);
	shark_dat.MassBalanceList[3].Set_Delta(4, 0);
	shark_dat.MassBalanceList[3].Set_Delta(5, 0);
	shark_dat.MassBalanceList[3].Set_Delta(6, 0);
	shark_dat.MassBalanceList[3].Set_Delta(7, 0);
	shark_dat.MassBalanceList[3].Set_Delta(8, 0);
	shark_dat.MassBalanceList[3].Set_Delta(9, 0);
	shark_dat.MassBalanceList[3].Set_Delta(10, 0);
	shark_dat.MassBalanceList[3].Set_Delta(11, 0);
	shark_dat.MassBalanceList[3].Set_Delta(12, 0);
	shark_dat.MassBalanceList[3].Set_Delta(13, 0);
	shark_dat.MassBalanceList[3].Set_Delta(14, 0);
	shark_dat.MassBalanceList[3].Set_Delta(15, 0);
	shark_dat.MassBalanceList[3].Set_Delta(16, 0);
	shark_dat.MassBalanceList[3].Set_Delta(17, 0);
	shark_dat.MassBalanceList[3].Set_Delta(18, 0);
	shark_dat.MassBalanceList[3].Set_Delta(19, 0);
	shark_dat.MassBalanceList[3].Set_Delta(20, 0);
	shark_dat.MassBalanceList[3].Set_Delta(21, 0);
	shark_dat.MassBalanceList[3].Set_Delta(22, 0);
	shark_dat.MassBalanceList[3].Set_Delta(23, 0);

	shark_dat.MassBalanceList[4].Set_Name("Uranium Balance");
	shark_dat.MassBalanceList[4].Set_TotalConcentration(UO2);
	shark_dat.MassBalanceList[4].Set_Delta(0, 0);
	shark_dat.MassBalanceList[4].Set_Delta(1, 0);
	shark_dat.MassBalanceList[4].Set_Delta(2, 0);
	shark_dat.MassBalanceList[4].Set_Delta(3, 0);
	shark_dat.MassBalanceList[4].Set_Delta(4, 0);
	shark_dat.MassBalanceList[4].Set_Delta(5, 0);
	shark_dat.MassBalanceList[4].Set_Delta(6, 0);
	shark_dat.MassBalanceList[4].Set_Delta(7, 0);
	shark_dat.MassBalanceList[4].Set_Delta(8, 1);
	shark_dat.MassBalanceList[4].Set_Delta(9, 1);
	shark_dat.MassBalanceList[4].Set_Delta(10, 1);
	shark_dat.MassBalanceList[4].Set_Delta(11, 1);
	shark_dat.MassBalanceList[4].Set_Delta(12, 1);
	shark_dat.MassBalanceList[4].Set_Delta(13, 2);
	shark_dat.MassBalanceList[4].Set_Delta(14, 3);
	shark_dat.MassBalanceList[4].Set_Delta(15, 1);
	shark_dat.MassBalanceList[4].Set_Delta(16, 1);
	shark_dat.MassBalanceList[4].Set_Delta(17, 1);
	shark_dat.MassBalanceList[4].Set_Delta(18, 0);
	shark_dat.MassBalanceList[4].Set_Delta(19, 0);
	shark_dat.MassBalanceList[4].Set_Delta(20, 0);
	shark_dat.MassBalanceList[4].Set_Delta(21, 0);
	shark_dat.MassBalanceList[4].Set_Delta(22, 1);
	shark_dat.MassBalanceList[4].Set_Delta(23, 1);

	shark_dat.MassBalanceList[5].Set_Name("Chlorine Balance");
	shark_dat.MassBalanceList[5].Set_TotalConcentration(NaCl+HCl);
	shark_dat.MassBalanceList[5].Set_Delta(0, 0);
	shark_dat.MassBalanceList[5].Set_Delta(1, 0);
	shark_dat.MassBalanceList[5].Set_Delta(2, 0);
	shark_dat.MassBalanceList[5].Set_Delta(3, 0);
	shark_dat.MassBalanceList[5].Set_Delta(4, 0);
	shark_dat.MassBalanceList[5].Set_Delta(5, 0);
	shark_dat.MassBalanceList[5].Set_Delta(6, 0);
	shark_dat.MassBalanceList[5].Set_Delta(7, 0);
	shark_dat.MassBalanceList[5].Set_Delta(8, 0);
	shark_dat.MassBalanceList[5].Set_Delta(9, 0);
	shark_dat.MassBalanceList[5].Set_Delta(10, 0);
	shark_dat.MassBalanceList[5].Set_Delta(11, 0);
	shark_dat.MassBalanceList[5].Set_Delta(12, 0);
	shark_dat.MassBalanceList[5].Set_Delta(13, 0);
	shark_dat.MassBalanceList[5].Set_Delta(14, 0);
	shark_dat.MassBalanceList[5].Set_Delta(15, 0);
	shark_dat.MassBalanceList[5].Set_Delta(16, 0);
	shark_dat.MassBalanceList[5].Set_Delta(17, 0);
	shark_dat.MassBalanceList[5].Set_Delta(18, 0);
	shark_dat.MassBalanceList[5].Set_Delta(19, 0);
	shark_dat.MassBalanceList[5].Set_Delta(20, 0);
	shark_dat.MassBalanceList[5].Set_Delta(21, 1);
	shark_dat.MassBalanceList[5].Set_Delta(22, 0);
	shark_dat.MassBalanceList[5].Set_Delta(23, 0);

/*
	shark_dat.AdsorptionList[0].setSpecificArea(ads_area);
	shark_dat.AdsorptionList[0].setSpecificMolality(ads_mol);
	shark_dat.AdsorptionList[0].setTotalMass(ads_mass);
	shark_dat.AdsorptionList[0].setSurfaceCharge(0.0);
	shark_dat.AdsorptionList[0].setAdsorbentName("A(OH)2");
	shark_dat.AdsorptionList[0].setActivityModelInfo(FloryHuggins, &shark_dat.AdsorptionList[0]);
	shark_dat.AdsorptionList[0].setBasis("area");
	//shark_dat.AdsorptionList[0].setBasis("molar");

	shark_dat.AdsorptionList[0].setMolarFactor(0, 1.0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Equilibrium(logK_UO2);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(0, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(1, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(2, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(3, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(4, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(5, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(6, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(7, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(8, -1);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(9, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(10, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(11, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(12, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(13, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(14, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(15, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(16, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(17, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(18, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(19, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(20, 2);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(21, 0);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(22, 1);
	shark_dat.AdsorptionList[0].getReaction(0).Set_Stoichiometric(23, 0);

	shark_dat.AdsorptionList[0].setMolarFactor(1, 1.0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Equilibrium(logK_UO2CO3);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(0, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(1, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(2, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(3, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(4, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(5, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(6, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(7, -1);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(8, -1);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(9, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(10, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(11, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(12, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(13, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(14, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(15, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(16, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(17, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(18, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(19, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(20, 2);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(21, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(22, 0);
	shark_dat.AdsorptionList[0].getReaction(1).Set_Stoichiometric(23, 1);


	shark_dat.AdsorptionList[0].setVolumeFactor(0, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(1, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(2, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(3, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(4, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(5, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(6, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(7, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(8, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(9, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(10, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(11, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(12, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(13, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(14, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(15, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(16, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(17, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(18, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(19, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(20, 0);
	shark_dat.AdsorptionList[0].setVolumeFactor(21, 0);
	 */
	// Only need to give volume factors for adsorbed species
	//shark_dat.AdsorptionList[0].setVolumeFactor(22, 5.0); //MADE UP NUMBER - given to the adsorbed species
	//shark_dat.AdsorptionList[0].setVolumeFactor(23, 8.0); //MADE UP NUMBER - given to the adsorbed species

	shark_dat.UnsteadyAdsList[0].setSpecificArea(ads_area);
	shark_dat.UnsteadyAdsList[0].setSpecificMolality(ads_mol);
	shark_dat.UnsteadyAdsList[0].setTotalMass(ads_mass);
	shark_dat.UnsteadyAdsList[0].setSurfaceCharge(0.0);
	shark_dat.UnsteadyAdsList[0].setAdsorbentName("HAO");
	shark_dat.UnsteadyAdsList[0].setActivityModelInfo(UNIQUAC_unsteady, &shark_dat.UnsteadyAdsList[0]);
	shark_dat.UnsteadyAdsList[0].setActivityEnum(UNIQUAC_ACT);
	shark_dat.UnsteadyAdsList[0].setBasis("molar");

	shark_dat.UnsteadyAdsList[0].setMolarFactor(0, 2.0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Equilibrium(logK_UO2);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Forward(1.0E6);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_InitialValue(0.0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(0, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(1, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(2, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(3, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(4, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(5, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(6, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(7, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(8, -1);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(9, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(10, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(11, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(12, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(13, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(14, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(15, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(16, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(17, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(18, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(19, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(20, 2);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(21, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(22, 1);
	shark_dat.UnsteadyAdsList[0].getReaction(0).Set_Stoichiometric(23, 0);


	shark_dat.UnsteadyAdsList[0].setMolarFactor(1, 2.0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Equilibrium(logK_UO2CO3);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Forward(1.0E3);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_InitialValue(0.0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(0, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(1, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(2, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(3, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(4, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(5, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(6, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(7, -1);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(8, -1);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(9, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(10, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(11, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(12, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(13, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(14, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(15, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(16, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(17, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(18, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(19, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(20, 2);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(21, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(22, 0);
	shark_dat.UnsteadyAdsList[0].getReaction(1).Set_Stoichiometric(23, 1);

	shark_dat.UnsteadyAdsList[0].setVolumeFactor(22, 5);
	shark_dat.UnsteadyAdsList[0].setVolumeFactor(23, 8);
	shark_dat.UnsteadyAdsList[0].calculateAreaFactors();    //NOTE: Required here because this is done during the reading step


	// END problem specific info here --------------------------------------------------------

	//Call the SHARK routine
	success = SHARK(&shark_dat);
	if (success != 0) {mError(simulation_fail); return -1;}

	//Close files and display end messages
	fclose(TestOutput);
	time = clock() - time;
	std::cout << "\nSimulation Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	std::cout << "Total Time Steps: " << shark_dat.timesteps << "\n";
	std::cout << "Total Iterations: " << shark_dat.totalsteps << "\n";
	std::cout << "Total Function Calls: " << shark_dat.totalcalls << "\n";
	std::cout << "Evaluations/sec: " << shark_dat.totalcalls/(time / CLOCKS_PER_SEC) << "\n";

	return success;
}
