//----------------------------------------
//  Created by Austin Ladshaw on 05/27/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 *		SHARK = Speciation-object Hierarchy for Aqueous Reaction Kinetics
 */

#include "mola.h"
#include "macaw.h"
#include "lark.h"
#include "yaml_wrapper.h"

#ifndef SHARK_HPP_
#define SHARK_HPP_

#ifndef Rstd
#define Rstd 8.3144621						// J/K/mol (or) L*kPa/K/mol (Standard Units)
#endif

//Master Species List Object
class MasterSpeciesList : Molecule
{
public:
	MasterSpeciesList();										//Default constructor
	~MasterSpeciesList();										//Default destructor
	MasterSpeciesList(const MasterSpeciesList &msl);			//Copy Constructor
	
	MasterSpeciesList& operator=(const MasterSpeciesList &msl);	//Equal operator
	
	void set_list_size(int i);									//Initialize the size of the list
	void set_species(int i, std::string formula);				//Set the ith species in the list
	void set_species(int i, int charge,
					 double enthalpy,
					 double entropy,
					 double energy,
					 bool HS,
					 bool G,
					 std::string Phase,
					 std::string Name,
					 std::string Formula,
					 std::string lin_formula);					//Set the ith species in the list to a custom object
	
	void DisplayInfo(int i);									//Function to display information of ith object
	void DisplayAll();											//Function to display all information of list
	void DisplayConcentrations(Matrix<double> &C);						//Function to display the concentrations of species in list
	
	void set_alkalinity(double alk);							//Set the alkalinity of the solution
	
	int list_size();											//Returns size of list
	Molecule& get_species(int i);								//Returns a reference to the ith species in master list
	int get_index(std::string name);							//Returns an integer representing location of the species in the object
	
	double charge(int i);										//Fetch and return charge of ith species in list
	double alkalinity();										//Fetch the value of alkalinity of the solution (mol/L)
	
	std::string speciesName(int i);								//Function to return the name of the ith species
	
	double Eval_ChargeResidual(const Matrix<double> &x);				//Calculate charge balance residual
	
protected:
	int size;													//Size of the list
	std::vector<Molecule> species;								//List of Molecule Objects
	double residual_alkalinity;									//Conc of strong base - conc of strong acid in solution (mol/L)
	
private:
	
	
};

//Reaction Object
class Reaction
{
public:
	Reaction();											//Default constructor
	~Reaction();										//Default destructor
	
	void Initialize_List (MasterSpeciesList &List);		//Initialize the Reaction object from the Species
	void Display_Info ();								//Display the species information in the list
	void Set_Stoichiometric (int i, double v);			//Set the ith stoichiometric value
	void Set_Equilibrium (double v);					//Set the equilibrium constant (logK)
	void Set_Enthalpy (double H);						//Set the enthalpy of the reaction (J/mol)
	void Set_Entropy (double S);						//Set the entropy of the reaction (J/K/mol)
	void Set_EnthalpyANDEntropy (double H, double S);	//Set both the enthalpy and entropy (J/mol) & (J/K/mol)
	void Set_Energy (double G);							//Set the Gibb's free energy of reaction (J/mol)
	
	void checkSpeciesEnergies();						//Function to check MasterList Reference for species energy info
	void calculateEnergies();							//Function to calculate the energy of the reaction
	void calculateEquilibrium(double T);				//Function to calculate the equilibrium constant
	bool haveEquilibrium();								//Function to return true if equilibrium constant is given or calculated
	
	double Get_Stoichiometric (int i);					//Fetch the ith stoichiometric value
	double Get_Equilibrium ();							//Fetch the equilibrium constant (logK)
	double Get_Enthalpy();								//Fetch the enthalpy of the reaction
	double Get_Entropy();								//Fetch the entropy of the reaction
	double Get_Energy();								//Fetch the energy of the reaction
	
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama);	//Evaluate a residual for the reaction given variable x and activity gama
	
protected:
	
	MasterSpeciesList *List;							//Pointer to a master species object
	std::vector<double> Stoichiometric;					//Vector of stoichiometric constants corresponding to species list
	double Equilibrium;									//Equilibrium constant for the reaction (logK)
	double enthalpy;									//Reaction enthalpy (J/mol)
	double entropy;										//Reaction entropy (J/K/mol)
	double energy;										//Gibb's Free energy of reaction (J/mol)
	bool CanCalcHS;										//True if all molecular info is avaiable to calculate dH and dS
	bool CanCalcG;										//True if all molecular info is available to calculate dG
	bool HaveHS;										//True if dH and dS is given, or can be calculated
	bool HaveG;											//True if dG is given, or can be calculated
	bool HaveEquil;										//True as long as Equilibrium is given, or can be calculated
	
private:
	
};

//Mass Balance Object
class MassBalance
{
public:
	MassBalance();
	~MassBalance();
	
	void Initialize_List (MasterSpeciesList &List);		//Initialize the Reaction object from the Species
	void Display_Info ();								//Display the species information in the list
	void Set_Delta (int i, double v);					//Set the ith delta value
	void Set_TotalConcentration (double v);				//Set the total concentration (mol/L)
	void Set_Name (std::string name);
	
	double Get_Delta (int i);							//Fetch the ith delta value
	double Sum_Delta ();								//Sums up the delta values and returns the total
	double Get_TotalConcentration ();					//Fetch the total concentration (mol/L)
	std::string Get_Name ();							//Return name of mass balance species
	
	double Eval_Residual(const Matrix<double> &x);				//Evaluate a residual for the reaction given variable x
	
protected:
	MasterSpeciesList *List;							//Pointer to a master species object
	std::vector<double> Delta;							//Vector of deltas used in mass balaces
	double TotalConcentration;							//Total concentration of specific object (mol/L)
	
private:
	std::string Name;									//Name designation used in mass balance
	
};

//Unsteady Reaction Object 
class UnsteadyReaction : Reaction
{
public:
	UnsteadyReaction();
	~UnsteadyReaction();
	
	void Initialize_List (MasterSpeciesList &List);		//Initialize the Reaction object from the Species
	void Display_Info ();								//Display the species information in the list
	
	void Set_Species_Index(int i);						//Set the Unsteady species index by number
	void Set_Species_Index(std::string formula);		//Set the Unsteady species index by formula
	void Set_Stoichiometric (int i, double v);			//Set the ith stoichiometric value
	void Set_Equilibrium (double v);					//Set the equilibrium constant (logK)
	void Set_Enthalpy (double H);						//Set the enthalpy of the reaction (J/mol)
	void Set_Entropy (double S);						//Set the entropy of the reaction (J/K/mol)
	void Set_EnthalpyANDEntropy (double H, double S);	//Set both the enthalpy and entropy (J/mol) & (J/K/mol)
	void Set_Energy (double G);							//Set the Gibb's free energy of reaction (J/mol)
	void Set_InitialValue (double ic);					//Set the initial value of the unsteady variable
	void Set_MaximumValue (double max);					//Set the maximum value of the unsteady variable
	void Set_Forward(double forward);					//Set the forward rate
	void Set_Reverse(double reverse);					//Set the reverse rate
	void Set_ForwardRef(double Fref);					//Set the forward reference rate
	void Set_ReverseRef(double Rref);					//Set the reverse reference rate
	void Set_ActivationEnergy(double E);				//Set the activation energy for the reaction
	void Set_Affinity(double b);						//Set the temperature affinity parameter
	void Set_TimeStep (double dt);						//Set the time step
	
	void checkSpeciesEnergies();						//Function to check MasterList Reference for species energy info
	void calculateEnergies();							//Function to calculate the energy of the reaction
	void calculateEquilibrium(double T);				//Function to calculate the equilibrium constant
	void calculateRate(double T);						//Function to calculate the rate constant
	bool haveEquilibrium();								//Function to return true if equilibrium constant is given or calculated
	bool haveRate();									//Function to return true if you have the forward or reverse rate calculated
	
	int Get_Species_Index();							//Fetch the index of the Unsteady species
	double Get_Stoichiometric (int i);					//Fetch the ith stoichiometric value
	double Get_Equilibrium ();							//Fetch the equilibrium constant (logK)
	double Get_Enthalpy();								//Fetch the enthalpy of the reaction
	double Get_Entropy();								//Fetch the entropy of the reaction
	double Get_Energy();								//Fetch the energy of the reaction
	double Get_InitialValue();							//Fetch the initial value of the variable
	double Get_MaximumValue();							//Fetch the maximum value of the variable
	double Get_Forward();								//Fetch the forward rate
	double Get_Reverse();								//Fetch the reverse rate
	double Get_ForwardRef();							//Fetch the forward reference rate
	double Get_ReverseRef();							//Fetch the reverse reference rate
	double Get_ActivationEnergy();						//Fetch the activation energy for the reaction
	double Get_Affinity();								//Fetch the temperature affinity for the reaction
	double Get_TimeStep();								//Fetch the time step
	
	double Eval_ReactionRate(const Matrix<double> &x, const Matrix<double> &gama);	//Calculate reation rate from concentrations and activities
	double Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old,
						 const Matrix<double> &gama_new, const Matrix<double> &gama_old); //Calculate the unsteady residual
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama);		//Calculate the steady-state residual
	double Eval_IC_Residual(const Matrix<double> &x);						//Calculate the unsteady residual for initial conditions
	double Explicit_Eval(const Matrix<double> &x, const Matrix<double> &gama);		//Return an approximate explicit solution
	
protected:
	double initial_value;								//Initial value given at t=0 (in mol/L)
	double max_value;									//Maximum value plausible (in mol/L)
	double forward_rate;								//Forward reaction rate constant (in (mol/L)^n/time)
	double reverse_rate;								//Reverse reaction rate constant (in (mol/L)^n/time)
	double forward_ref_rate;							//Forward reference rate constant (in (mol/L)^n/time)
	double reverse_ref_rate;							//Reverse reference rate constant (in (mol/L)^n/time)
	double activation_energy;							//Activation or barrier energy for the reaction (J/mol)
	double temperature_affinity;						//Temperature affinity parameter (dimensionless)
	double time_step;									//Time step size for current step
	bool HaveForward;									//True if can calculate, or was given the forward rate
	bool HaveReverse;									//True if can calculate, or was given the reverse rate
	bool HaveForRef;									//True if given the forward reference rate
	bool HaveRevRef;									//True if given the reverse reference rate
	int species_index;									//Index in MasterList of Unsteady Species
	
private:
	
};

//Reaction Mechanism Object
class Mechanism
{
public:
	
protected:
	// --------------------- NOTE: May have to radically change -----------------------------------------------
	MasterSpeciesList *List;							//Reference to the Master List of Species
	std::vector<UnsteadyReaction> reactions;			//Vector of individual reactions making up the mechanism
	std::vector<double> weight;							//The weight contributed by each Unsteady Reaction object
	int species_index;									//Index of the unsteady species of interest
	
private:
	
	
};

//Precipitation Reaction Object
class Precipitation : Reaction
{
public:
	
protected:
	
private:
	
	
};

//Unsteady Precipitation Reaction Object
class UnsteadyPrecipitation : Precipitation
{
public:
	
protected:
	
private:
	
	
};

typedef enum {IDEAL, DAVIES, DEBYE_HUCKEL, DAVIES_LADSHAW} valid_act;

typedef struct SHARK_DATA
{
	MasterSpeciesList MasterList;					//Master List of species object
	std::vector<Reaction> ReactionList;				//Equilibrium reaction objects
	std::vector<MassBalance> MassBalanceList;		//Mass balance objects
	std::vector<UnsteadyReaction> UnsteadyList;		//Unsteady Reaction objects
	
	std::vector<
		double (*) (const Matrix<double> &x,
					SHARK_DATA *shark_dat,
					const void *data) > OtherList;	//Array of Other Residual functions to be defined by user
													//NOTE: These must be setup individually by the user
	
	int numvar;										//Total number of functions and species
	int num_ssr;									//Number of steady-state reactions
	int num_mbe;									//Number of mass balance equations
	int num_usr;									//Number of unsteady-state reactions
	int num_other = 0;								//Number of other functions to be used (default is always 0)
	int act_fun = IDEAL;							//Number denoting the activity function to use (default is ideal)
	int totalsteps = 0;								//Number of iterations and function calls
	int timesteps = 0;								//Number of time steps taken to complete simulation
	int pH_index = -1;								//Contains the index of the pH variable
	int pOH_index = -1;								//Contains the index of the pOH variable
	
	double simulationtime = 0.0;					//Time to simulate unsteady reactions for (default = 0.0 for steady problem)
	double dt = 0.1;								//Time step size
	double dt_min = sqrt(DBL_EPSILON);				//Minimum allowable step size
	double t_out = 0.0;								//Time increment by which file output is made (default = always print)
	double t_count = 0.0;							//Running count of time increments
	double time = 0.0;								//Current value of time
	double time_old = 0.0;							//Previous value of time (start from t = 0.0)
	double pH = 7.0;								//Value of pH if needed (default = 7)
	double Norm = 0.0;								//Current value of euclidean norm in solution
	
	double dielectric_const = 78.325;				//Dielectric constant used in many activity models (default: water = 78.325 (1/K))
	double temperature = 298.15;					//Solution temperature (default = 25 oC or 298 K)
	
	bool steadystate = true;						//True = solve steady problem; False = solve transient problem
	bool TimeAdaptivity = false;					//True = solve using variable time step
	bool const_pH = false;							//True = set pH to a constant; False = solve for pH
	bool SpeciationCurve = false;					//True = runs a series of constant pH steady-state problems to produce curves
	bool Console_Output = true;						//True = display output to console
	bool File_Output = false;						//True = write output to a file
	bool Contains_pH = false;						//True = system contains pH as a variable
	bool Contains_pOH = false;						//True = system contains pOH as a variable
	bool Converged = false;							//True = system converged within tolerance 
	
	Matrix<double> X_old;									//Solution vector for old time step
	Matrix<double> X_new;									//Solution vector for current time step
	Matrix<double> Conc_old;								//Concentration vector for old time step
	Matrix<double> Conc_new;								//Concentration vector for current time step
	Matrix<double> activity_new;							//Activity matrix for current time step
	Matrix<double> activity_old;							//Activity matrix from prior time step
	
	//Function pointers for activity and residual evaluations
	int (*EvalActivity) (const Matrix<double>& x, Matrix<double> &F, const void *data);
	int (*Residual) (const Matrix<double>& x, Matrix<double> &F, const void *data);
	int (*lin_precon) (const Matrix<double> &r, Matrix<double> &p, const void *data);
	
	PJFNK_DATA Newton_data;							//Data structure for the Newton-Krylov solver
	const void *activity_data;						//User defined data structure for an activity model
	const void *residual_data;						//User defined data structure for the residual function
	const void *precon_data;						//User defined data structure for preconditioning
	const void *other_data;							//User define data structure used for user defined residuals
	FILE *OutputFile;								//Output File pointer
	
	yaml_cpp_class yaml_object;						//yaml object to read and access digitized yaml documents
	
} SHARK_DATA;

void print2file_shark_info(SHARK_DATA *shark_dat);

void print2file_shark_header(SHARK_DATA *shark_dat);

void print2file_shark_results_new(SHARK_DATA *shark_dat);

void print2file_shark_results_old(SHARK_DATA *shark_dat);

int ideal_solution (const Matrix<double>& x, Matrix<double> &F, const void *data);

int Davies_equation (const Matrix<double>& x, Matrix<double> &F, const void *data);

int DebyeHuckel_equation (const Matrix<double> &x, Matrix<double> &F, const void *data);

int DaviesLadshaw_equation (const Matrix<double>& x, Matrix<double> &F, const void *data);

int act_choice(const std::string &input);

bool linesearch_choice(const std::string &input);

int linearsolve_choice(const std::string &input);

int Convert2LogConcentration(const Matrix<double> &x, Matrix<double> &logx);

int Convert2Concentration(const Matrix<double> &logx, Matrix<double> &x);

int read_scenario(SHARK_DATA *shark_dat);

int read_options(SHARK_DATA *shark_dat);

int read_species(SHARK_DATA *shark_dat);

int read_massbalance(SHARK_DATA *shark_dat);

int read_equilrxn(SHARK_DATA *shark_dat);

int read_unsteadyrxn(SHARK_DATA *shark_dat);

int setup_SHARK_DATA( FILE *file, int (*residual) (const Matrix<double> &x, Matrix<double> &res, const void *data),
					  int (*activity) (const Matrix<double> &x, Matrix<double> &gama, const void *data),
					  int (*precond) (const Matrix<double> &r, Matrix<double> &p, const void *data),
					  SHARK_DATA *dat, const void *activity_data, const void *residual_data,
					  const void *precon_data, const void *other_data);

int shark_add_customResidual(int i, double (*other_res) (const Matrix<double> &x, SHARK_DATA *shark_dat, const void *other_data),
							 SHARK_DATA *shark_dat);

int shark_parameter_check(SHARK_DATA *shark_dat);

int shark_energy_calculations(SHARK_DATA *shark_dat);

int shark_temperature_calculations(SHARK_DATA *shark_dat);

int shark_pH_finder(SHARK_DATA *shark_dat);

int shark_guess(SHARK_DATA *shark_dat);

int shark_initial_conditions(SHARK_DATA *shark_dat);

int shark_executioner(SHARK_DATA *shark_dat);

int shark_timestep_const(SHARK_DATA *shark_dat);

int shark_timestep_adapt(SHARK_DATA *shark_dat);

int shark_preprocesses(SHARK_DATA *shark_dat);

int shark_solver(SHARK_DATA *shark_dat);

int shark_postprocesses(SHARK_DATA *shark_dat);

int shark_reset(SHARK_DATA *shark_dat);

int shark_residual(const Matrix<double> &x, Matrix<double> &F, const void *data);

int SHARK(SHARK_DATA *shark_dat);

int SHARK_SCENARIO(const char *yaml_input);

int SHARK_TESTS();

#endif
