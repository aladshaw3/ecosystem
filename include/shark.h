/*!
 *  \file shark.h shark.cpp
 *	\brief Speciation-object Hierarchy for Aqueous Reaction Kinetics
 *	\details This file contains structures and functions associated with solving speciation and kinetic
 *			problems in aqueous systems. The primary aim for the development of these algorithms was to
 *			solve speciation and adsorption problems for the recovery of uranium resources from seawater.
 *			Seawater is an extradorinarily complex medium in which to work, which is why these algorithms
 *			are being constructed in a piece-wise, object-oriented fashion. This allows us to displace
 *			much of the complexity of the problem by breaking it down into smaller, more managable pieces.
 *
 *			Each piece of SHARK contributes to a residual function when solving the overall speciation,
 *			reaction, kinetic chemical problem. These residuals are then fed into the PJFNK solver function
 *			in lark.h. The variables of the system are the log(C) concentration values of each species in 
 *			the system. We solve for log(C) concentrations, rather than just C, because the PJFNK method
 *			is an unbounded solution algorithm. So to prevent the algorithm from producing negative values
 *			for concentration, we reformulate all residuals in terms of the log(C) values. In this way,
 *			regardless of the value found for log(C), the concentration C will always be greater than 0.
 *
 *			Currenty, SHARK supports standard aqueous speciation problems with simple kinetic models based
 *			on an unsteady form of the standard reaction stoichiometry. As more methods and algorithms are
 *			completed, the SHARK simulations will be capable of doing much, much more.
 *
 *	\warning Much of this is still underconstruction and many methods or interfaces may change. Use
 *			with caution.
 *  \author Austin Ladshaw
 *	\date 05/27/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved.
 */

#include "mola.h"
#include "macaw.h"
#include "lark.h"
#include "yaml_wrapper.h"

#ifndef SHARK_HPP_
#define SHARK_HPP_

#ifndef Rstd
#define Rstd 8.3144621						///< Gas Law Constant in J/K/mol (or) L*kPa/K/mol (Standard Units)
#endif

/// Master Species List Object
/** C++ style object that holds data and function associated with solving multi-species problems. This
	object contains a vector of Molecule objects from mola.h and uses those objects to help setup speciation
	problems that need to be solved. One of the primary functions in this object is the contribution of
	electroneutrality (Eval_ChargeResidual). However, we only need this constraint if the pH of our aqueous
	system is unknown. */
class MasterSpeciesList
{
public:
	MasterSpeciesList();										///< Default constructor
	~MasterSpeciesList();										///< Default destructor
	MasterSpeciesList(const MasterSpeciesList &msl);			///< Copy Constructor
	
	MasterSpeciesList& operator=(const MasterSpeciesList &msl);	///< Equals operator
	
	void set_list_size(int i);									///< Function to initialize the size of the list
	
	/// Function to register the ith species in the list based on a registered molecular formula (see mola.h)
	void set_species(int i, std::string formula);
	
	/// Function to register the ith species in the list based on custom molecule information (see mola.h)
	void set_species(int i, int charge,
					 double enthalpy,
					 double entropy,
					 double energy,
					 bool HS,
					 bool G,
					 std::string Phase,
					 std::string Name,
					 std::string Formula,
					 std::string lin_formula);
	
	void DisplayInfo(int i);									///< Function to display information of ith object
	
	void DisplayAll();											///< Function to display all information of list
	
	/// Function to display the concentrations of species in list
	/** This function will print to the console the species list in order with each species associated
		concentration from the matrix C. The concentrations and species list MUST be in the same order
		and the units of C are assumed to be mol/L.
	 
		\param C matrix of concentrations of species in the list in mol/L*/
	void DisplayConcentrations(Matrix<double> &C);
	
	/// Set the alkalinity of the solution (Default = 0 M)
	/** This function is used to set the value of residual alkalinity used in the electroneutrality calculations.
		Typically, this value will be 0 M (mol/L) if all species in the system are present as variables. However,
		occasionally, one may want to set the alkalinity of the solution to a constant in order to restrict the
		pH of the solution. 
	 
		\param alk Residual alkalinity in M (mol/L)*/
	void set_alkalinity(double alk);
	
	int list_size();											///< Returns size of list
	
	/// Returns a reference to the ith species in master list
	/** This function will return a Molecule object for the ith species in the list of molecules.
		Once returned, the user then can operate on that molecule using the functions define in mola.h. */
	Molecule& get_species(int i);
	
	int get_index(std::string name);							///< Returns an integer representing location of the named species in the list
	
	double charge(int i);										///< Fetch and return charge of ith species in list
	
	double alkalinity();										///< Fetch the value of alkalinity of the solution (mol/L)
	
	std::string speciesName(int i);								///< Function to return the name of the ith species
	
	/// Calculate charge balance residual for the electroneutrality constraint
	/** This function returns the value of the residual for the electroneutrality equation in the system.
		Electroneutrality is based on the concentrations and charges of each species in the system so the
		charges of each molecule must be appropriately set. Concentrations of those species are fed into 
		this function via the argument x, but come in as the log(C) values (i.e., x = log(C)).
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double Eval_ChargeResidual(const Matrix<double> &x);
	
protected:
	int size;													///< Size of the list
	std::vector<Molecule> species;								///< List of Molecule Objects
	double residual_alkalinity;									///< Conc of strong base - conc of strong acid in solution (mol/L)
	
private:
	
};

/// Reaction Object
/** C++ style object that holds data and functions associated with standard chemical reactions... \n
 
	i.e., aA + bB <=> cC + dD \n
 
	These reactions are assumed steady state and are characterized by stoichiometry coefficients
	and equilibrium/stability constants. Types of reactions that these are valid for would be
	acid/base reactions, metal-ligand complexation reactions, oxidation-reduction reactions, 
	Henry's Law phase changes, and more. Reactions that this may not be suitable for include
	mechanisms, adsorption, and precipitation. Those types of reactions would be better handled
	by more specific objects that inherit from this object. \n
 
	If all species in the reaction are registered and known species in mola.h AND have known
	formation energies, then the equilibrium constants for that particular reaction will be 
	calculated based on the species involved in the reaction. However, if using some custom
	molecule objects, then the reaction equilibrium may not be able to be automatically formed
	by the routine. In this case, you would need to also supply the equilibrium constant for
	the particular reaction.
	*/
class Reaction
{
public:
	Reaction();											///< Default constructor
	~Reaction();										///< Default destructor
	
	void Initialize_List (MasterSpeciesList &List);		///< Function to initialize the Reaction object from the MasterSpeciesList
	void Display_Info ();								///< Display the reaction information
	
	/// Set the ith stoichiometric value
	/** This function will set the stoichiometric constant of the ith species in the master list to 
		the given value of v. All values of v are set to zero unless overriden by this function. 
	 
		\param i index of the species in the MasterSpeciesList
		\param v value of the stoichiometric constant for that species in the reaction*/
	void Set_Stoichiometric (int i, double v);
	void Set_Equilibrium (double logK);					///< Set the equilibrium constant in log(K) units
	void Set_Enthalpy (double H);						///< Set the enthalpy of the reaction (J/mol)
	void Set_Entropy (double S);						///< Set the entropy of the reaction (J/K/mol)
	void Set_EnthalpyANDEntropy (double H, double S);	///< Set both the enthalpy and entropy (J/mol) & (J/K/mol)
	void Set_Energy (double G);							///< Set the Gibb's free energy of reaction (J/mol)
	
	/// Function to check MasterList Reference for species energy info
	/** This function will go through the stoichiometry of this reaction and check the molecules
		in the MasterSpeciesList that correspond to the species present in this reaction for the
		existance of their formation energies. Based on the states of those energies, it will 
		note internally whether or not it can determine the equilibrium constants based soley
		on individual species information. If it cannot, then the user must provide either the
		reaction energies to form the equilibrium constant or the equilibrium constant itself. */
	void checkSpeciesEnergies();
	
	///< Function to calculate and set the energy of the reaction
	/** If the energies of the reaction can be determined from the individual species in the 
		reaction, then this function uses that information. Otherwise, it sets the energies
		equal to the constants given to the object by the user. */
	void calculateEnergies();
	
	void calculateEquilibrium(double T);				///< Function to calculate the equilibrium constant based on temperature in K
	
	bool haveEquilibrium();								///< Function to return true if equilibrium constant is given or can be calculated
	
	double Get_Stoichiometric (int i);					///< Fetch the ith stoichiometric value
	double Get_Equilibrium ();							///< Fetch the equilibrium constant (logK)
	double Get_Enthalpy();								///< Fetch the enthalpy of the reaction (J/mol)
	double Get_Entropy();								///< Fetch the entropy of the reaction (J/K/mol)
	double Get_Energy();								///< Fetch the energy of the reaction (J/mol)
	
	///< Evaluate a residual for the reaction given variable x=log(C) and activity coefficients gama
	/** This function will calculate the reaction residual from this object's stoichiometry, equilibrium constant,
		log(C) concentrations, and activity coefficients.
	 
		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step*/
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama);
	
protected:
	
	MasterSpeciesList *List;							///< Pointer to a master species object
	std::vector<double> Stoichiometric;					///< Vector of stoichiometric constants corresponding to species list
	double Equilibrium;									///< Equilibrium constant for the reaction (logK)
	double enthalpy;									///< Reaction enthalpy (J/mol)
	double entropy;										///< Reaction entropy (J/K/mol)
	double energy;										///< Gibb's Free energy of reaction (J/mol)
	bool CanCalcHS;										///< True if all molecular info is avaiable to calculate dH and dS
	bool CanCalcG;										///< True if all molecular info is available to calculate dG
	bool HaveHS;										///< True if dH and dS is given, or can be calculated
	bool HaveG;											///< True if dG is given, or can be calculated
	bool HaveEquil;										///< True as long as Equilibrium is given, or can be calculated
	
private:
	
};

/// Mass Balance Object
/** C++ style object that holds data and functions associated with mass balances of primary species
	in a system. The mass balances involve a total concentration (in mol/L) and a vector of weighted
	contributions to that total concentration from each species in the MasterSpeciesList. This object
	only considers mass balances in a batch type of system (i.e., not input or output of mass). However,
	one could inherit from this object to create mass balances for flow systems as well. */
class MassBalance
{
public:
	MassBalance();										///< Default Constructor
	~MassBalance();										///< Default Destructor
	
	void Initialize_List (MasterSpeciesList &List);		///< Function to initialize the MassBalance object from the MasterSpeciesList
	void Display_Info ();								///< Display the mass balance information
	
	/// Function to set the ith weight (delta) for the mass balance
	/** This function sets the weight (i.e., delta value) of the ith species in the list
		to the value of v. That value represents the weighting of that species in the 
		determination of the total mass for the primary species set. 
	 
		\param i index of the species in the MasterSpeciesList
		\param v value of the weigth (or delta) applied to the mass balance*/
	void Set_Delta (int i, double v);
	
	void Set_TotalConcentration (double v);				///< Set the total concentration of the mass balance to v (mol/L)
	
	void Set_Name (std::string name);					///< Set the name of the mass balance (i.e., Uranium, Carbonate, etc.)
	
	double Get_Delta (int i);							///< Fetch the ith weight (i.e., delta) value
	double Sum_Delta ();								///< Sums up the delta values and returns the total (should never be zero)
	double Get_TotalConcentration ();					///< Fetch the total concentration (mol/L)
	std::string Get_Name ();							///< Return name of mass balance object
	
	/// Evaluate the residual for the mass balance object given the log(C) concentrations
	/** This function calculates and provides the residual for this mass balance object based on the total
		concentration in the system and the weighted contributions from each species. Concentrations are 
		given as the log(C) values.
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double Eval_Residual(const Matrix<double> &x);
	
protected:
	MasterSpeciesList *List;							///< Pointer to a master species object
	std::vector<double> Delta;							///< Vector of weights (i.e., deltas) used in the mass balance
	double TotalConcentration;							///< Total concentration of specific object (mol/L)
	
private:
	std::string Name;									///< Name designation used in mass balance
	
};

/// Unsteady Reaction Object (inherits from Reaction)
/** C++ style object that holds data and functions associated with unsteady chemical reactions... \n
	
	i.e., aA + bB <-- reverse  :  forward --> cC + dD \n
 
	This is essentially the same as the steady reaction, but we now have a forward and reverse
	reaction rate to deal with. It should be noted that this is a very simple kinetic reaction
	model based on splitting an overall equilibrium reaction into an overall forward and reverse
	reaction model. Therefore, it is not expected that this representation of the reaction will
	provide high accuracy results for reaction kinetics, but should at least provide an overall
	idea of the process occuring. */
class UnsteadyReaction : Reaction
{
public:
	UnsteadyReaction();									///< Default Constructor
	~UnsteadyReaction();								///< Default Destructor
	
	void Initialize_List (MasterSpeciesList &List);		///< Function to initialize the UnsteadyReaction object from the MasterSpeciesList
	void Display_Info ();								///< Display the unsteady reaction information
	
	/// Set the Unsteady species index by number
	/** This function will set the unsteady species index by the index i given. That
		given index must correspond to the index of the species in the MasterSpeciesList
		that is being considered as the unsteady species. 
	 
		\param i index of the unsteady species in the MasterSpeciesList*/
	void Set_Species_Index(int i);
	
	/// Set the Unsteady species index by formula
	/** This function will check the MasterSpeciesList for the molecule object that has
		the given formula, then set the unsteady species index based on the index of that
		species in the master list. 
		
		\param formula molecular formula of the unsteady species (see mola.h for standard formatting) */
	void Set_Species_Index(std::string formula);
	
	void Set_Stoichiometric (int i, double v);			///< Set the ith stoichiometric value (see Reaction object)
	void Set_Equilibrium (double v);					///< Set the equilibrium constant (logK) (see Reaction object)
	void Set_Enthalpy (double H);						///< Set the enthalpy of the reaction (J/mol) (see Reaction object)
	void Set_Entropy (double S);						///< Set the entropy of the reaction (J/K/mol) (see Reaction object)
	void Set_EnthalpyANDEntropy (double H, double S);	///< Set both the enthalpy and entropy (J/mol) & (J/K/mol) (see Reaction object)
	void Set_Energy (double G);							///< Set the Gibb's free energy of reaction (J/mol) (see Reaction object)
	
	/// Set the initial value of the unsteady variable
	/** This function sets the initial concentration value for the unsteady species to the given value ic
		(mol/L). Only unsteady species need to be given an initial value. All other species initial values
		for the overall system is setup based on a speciation calculation performed while holding the unsteady
		variables constant at their respective initial values. 
	 
		\param ic initial concentration value for the unsteady object (mol/L)*/
	void Set_InitialValue (double ic);
	
	/// Set the maximum value of the unsteady variable to a given value max (mol/L)
	/** This function will be called internally to help bound the unsteady variable to reasonable
		maximum values. That maximum is usually based on the mass balances for the current non-linear
		iteration. 
	 
		\param max maximum allowable value for the unsteady variable (mol/L)*/
	void Set_MaximumValue (double max);
	
	void Set_Forward(double forward);					///< Set the forward rate for the reaction (mol/L/hr)
	void Set_Reverse(double reverse);					///< Set the reverse rate for the reaction (mol/L/hr)
	
	/// Set the forward reference rate (mol/L/hr)
	/** Unlike just setting the forward rate, this function sets a reference forward rate of the reaction
		that can be used to correct the overall forward rate based on system temperature and Arrhenius Rate
		Equation constants. 
	 
		\param Fref forward reference rate constant (mol/L/hr)*/
	void Set_ForwardRef(double Fref);
	
	/// Set the reverse reference rate (mol/L/hr)
	/** Unlike just setting the reverse rate, this function sets a reference reverse rate of the reaction
		that can be used to correct the overall reverse rate based on system temperature and Arrhenius Rate
		Equation constants.
	 
		\param Rref reverse reference rate constant (mol/L/hr)*/
	void Set_ReverseRef(double Rref);
	
	/// Set the activation energy for the reaction (J/mol)
	/** This function will set the activation energy for the reaction to the given value of E. Note that
		we will only set one value for activation energy, even though there are rates for forward and
		reverse reactions. This is because we use the ratio of the rates and the equilibrium constant to
		establish the other rate. Therefore, we only need either the forward or reverse rate and the 
		equilibrium constant to set all the rates. 
	 
		\param E activation energy for the forward or reverse rate, depending on which was given*/
	void Set_ActivationEnergy(double E);
	
	/// Set the temperature affinity parameter for the reaction
	/** This function will set the temperature affinity for the reaction to the given value of b. Note that
		we will only set one value for temperature affinity, even though there are rates for forward and
		reverse reactions. This is because we use the ratio of the rates and the equilibrium constant to
		establish the other rate. Therefore, we only need either the forward or reverse rate and the
		equilibrium constant to set all the rates.
	 
		\param b temperature affinity for the forward or reverse rate, depending on which was given*/
	void Set_Affinity(double b);
	
	void Set_TimeStep (double dt);						///< Set the time step for the current simulation
	
	void checkSpeciesEnergies();						///< Function to check MasterSpeciesList for species energy info (see Reaction object)
	void calculateEnergies();							///< Function to calculate the energy of the reaction (see Reaction object)
	void calculateEquilibrium(double T);				///< Function to calculate the equilibrium constant (see Reaction object)
	
	/// Function to calculate the rate constant based on given temperature
	/** This function will calculate and set either the forward or reverse rate for the unsteady
		reaction based on what information was given. If the forward rate information was given,
		then it sets the reverse rate and visa versa. If nothing was set correctly, an error will
		occur.
	 
		\param T temperature of the system in Kelvin*/
	void calculateRate(double T);
	bool haveEquilibrium();								///< True if equilibrium constant is given or can be calculated (see Reaction object)
	bool haveRate();									///< Function to return true if you have the forward or reverse rate calculated
	
	int Get_Species_Index();							///< Fetch the index of the Unsteady species
	double Get_Stoichiometric (int i);					///< Fetch the ith stoichiometric value
	double Get_Equilibrium ();							///< Fetch the equilibrium constant (logK)
	double Get_Enthalpy();								///< Fetch the enthalpy of the reaction
	double Get_Entropy();								///< Fetch the entropy of the reaction
	double Get_Energy();								///< Fetch the energy of the reaction
	double Get_InitialValue();							///< Fetch the initial value of the variable
	double Get_MaximumValue();							///< Fetch the maximum value of the variable
	double Get_Forward();								///< Fetch the forward rate
	double Get_Reverse();								///< Fetch the reverse rate
	double Get_ForwardRef();							///< Fetch the forward reference rate
	double Get_ReverseRef();							///< Fetch the reverse reference rate
	double Get_ActivationEnergy();						///< Fetch the activation energy for the reaction
	double Get_Affinity();								///< Fetch the temperature affinity for the reaction
	double Get_TimeStep();								///< Fetch the time step
	
	/// Calculate reation rate (dC/dt) from concentrations and activities
	/** This function calculates the right hand side of the unsteady reaction equation based on the available
		rates, the current values of the non-linear variables (x=log(C)), and the activity coefficients (gama).
  
		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step*/
	double Eval_ReactionRate(const Matrix<double> &x, const Matrix<double> &gama);
	
	/// Calculate the unsteady residual for the reaction using and implicit time discretization
	/** This function uses the current time step and states of the non-linear variables and activities
		to form the residual contribution of the unsteady reaction. The time dependent functions are 
		discretized using an implicit finite difference for best stability. 
	 
		\param x_new matrix of the log(C) concentration values at the current non-linear step
		\param gama_new matrix of activity coefficients for each species at the current non-linear step
		\param x_old matrix of the log(C) concentration values at the previous non-linear step
		\param gama_old matrix of activity coefficients for each species at the previous non-linear step*/
	double Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old,
						 const Matrix<double> &gama_new, const Matrix<double> &gama_old);
	
	/// Calculate the steady-state residual for this reaction (see Reaction object)
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama);
	
	/// Calculate the unsteady residual for initial conditions
	/** Setting the intial conditions for all variables in the system requires a speciation calculation.
		However, we want the unsteady variables to be set to their respective initial conditions. Using this
		residual function imposes an equality constraint on those non-linear, unsteady variables allowing the
		rest of the speciation problem to be solved via PJFNK iterations. 
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double Eval_IC_Residual(const Matrix<double> &x);
	
	/// Return an approximate explicit solution to our unsteady variable (mol/L)
	/** This function will approximate the concentration of the unsteady variables based on an explicit time
		discretization. The purpose of this function is to try to provide the PJFNK method with a good initial
		guess for the values of the non-linear, unsteady variables. If we do not provide a good initial guess
		to these variables, then the PJFNK method may not converge to the correct solution, because the unsteady
		problem is the most difficult to solve. 
	 
		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step*/
	double Explicit_Eval(const Matrix<double> &x, const Matrix<double> &gama);
	
protected:
	double initial_value;								///< Initial value given at t=0 (in mol/L)
	double max_value;									///< Maximum value plausible (in mol/L)
	double forward_rate;								///< Forward reaction rate constant (in (mol/L)^n/hr)
	double reverse_rate;								///< Reverse reaction rate constant (in (mol/L)^n/hr)
	double forward_ref_rate;							///< Forward reference rate constant (in (mol/L)^n/hr)
	double reverse_ref_rate;							///< Reverse reference rate constant (in (mol/L)^n/hr)
	double activation_energy;							///< Activation or barrier energy for the reaction (J/mol)
	double temperature_affinity;						///< Temperature affinity parameter (dimensionless)
	double time_step;									///< Time step size for current step
	bool HaveForward;									///< True if can calculate, or was given the forward rate
	bool HaveReverse;									///< True if can calculate, or was given the reverse rate
	bool HaveForRef;									///< True if given the forward reference rate
	bool HaveRevRef;									///< True if given the reverse reference rate
	int species_index;									///< Index in MasterList of Unsteady Species
	
private:
	
};

/// \cond

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

/// \endcond

/// Enumeration for the list of valid activity models for non-ideal solutions
/** \note The SIT and PITZER models are not currently supported. */
typedef enum {IDEAL, DAVIES, DEBYE_HUCKEL, DAVIES_LADSHAW, SIT, PITZER} valid_act;

/// Data structure for SHARK simulations
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
