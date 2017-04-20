/*!
 *  \file shark.h shark.cpp
 *	\brief Speciation-object Hierarchy for Adsorption Reactions and Kinetics
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
#include "dogfish.h"

#ifndef SHARK_HPP_
#define SHARK_HPP_

#ifndef Rstd
#define Rstd 8.3144621						///< Gas Law Constant in J/K/mol (or) L*kPa/K/mol (Standard Units)
#endif

#ifndef	Na
#define Na 6.0221413E+23					///< Avagadro's Number - Units: molecules/mol
#endif

#ifndef	kB
#define kB 1.3806488E-23					///< Boltzmann's Constant - Units: J/K or C*V/K
#endif

#ifndef e
#define e 1.6021766208E-19					///< Elementary Electric Charge - Units: C
#endif

#ifndef Faraday
#define Faraday 96485.33289					///< Faraday's Constant - C/mol
#endif

#ifndef VolumeSTD
#define VolumeSTD 15.17						///< Standard Segment Volume - cm^3/mol
#endif

#ifndef AreaSTD
#define AreaSTD 2.5E5						///< Standard Segment Area - m^2/mol
#endif

#ifndef CoordSTD
#define CoordSTD 10							///< Standard Coordination Number
#endif

#ifndef LengthFactor
#define LengthFactor(z,r,s) (((z/2.0)*(r-s)) - (r-1.0))	///< Calculation of the Length Factor Parameter in UNIQUAC
#endif

#ifndef VacuumPermittivity
#define VacuumPermittivity 8.8541878176E-12	///< Vacuum Permittivity Constant - F/m or C/V/m
#endif

#ifndef WaterRelPerm
#define WaterRelPerm 80.1					///< Approximate Relative Permittivity for water - Unitless
#endif

#ifndef AbsPerm
#define AbsPerm(Rel) (Rel*VacuumPermittivity)	///< Calculation of Absolute Permittivity of a medium - F/m or C/V/m
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

	void Initialize_Object(MasterSpeciesList &List);	///< Function to initialize the Reaction object from the MasterSpeciesList
	void Display_Info();								///< Display the reaction information

	/// Set the ith stoichiometric value
	/** This function will set the stoichiometric constant of the ith species in the master list to
		the given value of v. All values of v are set to zero unless overriden by this function.

		\param i index of the species in the MasterSpeciesList
		\param v value of the stoichiometric constant for that species in the reaction*/
	void Set_Stoichiometric(int i, double v);
	void Set_Equilibrium(double logK);					///< Set the equilibrium constant in log(K) units
	void Set_Enthalpy(double H);						///< Set the enthalpy of the reaction (J/mol)
	void Set_Entropy(double S);							///< Set the entropy of the reaction (J/K/mol)
	void Set_EnthalpyANDEntropy (double H, double S);	///< Set both the enthalpy and entropy (J/mol) & (J/K/mol)
	void Set_Energy(double G);							///< Set the Gibb's free energy of reaction (J/mol)

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

	double Get_Stoichiometric(int i);					///< Fetch the ith stoichiometric value
	double Get_Equilibrium();							///< Fetch the equilibrium constant (logK)
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

/// Enumeration for the list of valid activity models for non-ideal solutions
/** \note The SIT and PITZER models are not currently supported. */
typedef enum {BATCH, CSTR, PFR} valid_mb;

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

	void Initialize_Object(MasterSpeciesList &List);	///< Function to initialize the MassBalance object from the MasterSpeciesList
	void Display_Info();								///< Display the mass balance information

	/// Function to set the ith weight (delta) for the mass balance
	/** This function sets the weight (i.e., delta value) of the ith species in the list
		to the value of v. That value represents the weighting of that species in the
		determination of the total mass for the primary species set.

		\param i index of the species in the MasterSpeciesList
		\param v value of the weigth (or delta) applied to the mass balance*/
	void Set_Delta(int i, double v);

	void Set_TotalConcentration(double v);				///< Set the total concentration of the mass balance to v (mol/L)
	void Set_Type(int type);							///< Set the Mass Balance type to BATCH, CSTR, or PFR
	void Set_Volume(double v);							///< Set the volume of the reactor
	void Set_FlowRate(double v);						///< Set the flow rate for the CSTR or PFR
	void Set_Area(double v);							///< Set the cross sectional area for the PFR
	void Set_TimeStep(double v);						///< Set the time step for the CSTR or PFR
	void Set_InitialConcentration(double v);			///< Set the initial concentration for the mass balance
	void Set_InletConcentration(double v);				///< Set the inlet concentration for the CSTR or PFR
	void Set_SteadyState(bool ss);						///< Set the boolean for Steady-State simulation
	void Set_ZeroInitialSolids(bool solids);			///< Set the boolean for initial solids in solution

	void Set_Name(std::string name);					///< Set the name of the mass balance (i.e., Uranium, Carbonate, etc.)

	double Get_Delta(int i);							///< Fetch the ith weight (i.e., delta) value
	double Sum_Delta();									///< Sums up the delta values and returns the total (should never be zero)
	double Get_TotalConcentration();					///< Fetch the total concentration (mol/L)
	int Get_Type();										///< Fetch the reactor type
	double Get_Volume();								///< Fetch the reactor volume
	double Get_FlowRate();								///< Fetch the reactor flow rate
	double Get_Area();									///< Fetch the reactor cross section area
	double Get_TimeStep();								///< Fetch the time step
	double Get_InitialConcentration();					///< Fetch the initial concentration
	double Get_InletConcentration();					///< Fetch the inlet concentration
	bool isSteadyState();								///< Fetch the steady-state condition
	bool isZeroInitialSolids();							///< Fetch the initial solids condition
	std::string Get_Name();								///< Return name of mass balance object

	/// Evaluate the residual for the mass balance object given the log(C) concentrations
	/** This function calculates and provides the residual for this mass balance object based on the total
		concentration in the system and the weighted contributions from each species. Concentrations are
		given as the log(C) values.

		\param x_new matrix of the log(C) concentration values at the current non-linear step
	    \param x_old matrix of the old log(C) concentration values for transient simulations*/
	double Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old);
	
	/// Evaluate the initial residual for the unsteady mass balance object given the log(C) concentrations
	/** This function calculates and provides the initial residual for this mass balance object based on the initial
		concentration in the system and the weighted contributions from each species. Concentrations are
		given as the log(C) values.
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double Eval_IC_Residual(const Matrix<double> &x);

protected:
	MasterSpeciesList *List;							///< Pointer to a master species object
	std::vector<double> Delta;							///< Vector of weights (i.e., deltas) used in the mass balance
	double TotalConcentration;							///< Total concentration of specific object (mol/L)
	
	int Type;											///< Type of mass balance object (default = BATCH)
	double volume;										///< Volume of the reactor (L)
	double flow_rate;									///< Volumetric flow rate in reactor (L/hr)
	double xsec_area;									///< Cross sectional area in PFR configuration (m^2)
	double dt;											///< Time step for non-batch case (hrs)
	double InitialConcentration;						///< Concentration initially in the domain (mol/L)
	double InletConcentration;							///< Concentration in the inlet of the domain (mol/L)
	bool SteadyState;									///< True if running steady-state simulation
	bool ZeroInitialSolids;								///< True if zero solids present for initial condition

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

	void Initialize_Object(MasterSpeciesList &List);	///< Function to initialize the UnsteadyReaction object from the MasterSpeciesList
	void Display_Info();								///< Display the unsteady reaction information

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

	void Set_Stoichiometric(int i, double v);			///< Set the ith stoichiometric value (see Reaction object)
	void Set_Equilibrium(double v);						///< Set the equilibrium constant (logK) (see Reaction object)
	void Set_Enthalpy(double H);						///< Set the enthalpy of the reaction (J/mol) (see Reaction object)
	void Set_Entropy(double S);							///< Set the entropy of the reaction (J/K/mol) (see Reaction object)
	void Set_EnthalpyANDEntropy(double H, double S);	///< Set both the enthalpy and entropy (J/mol) & (J/K/mol) (see Reaction object)
	void Set_Energy(double G);							///< Set the Gibb's free energy of reaction (J/mol) (see Reaction object)

	/// Set the initial value of the unsteady variable
	/** This function sets the initial concentration value for the unsteady species to the given value ic
		(mol/L). Only unsteady species need to be given an initial value. All other species initial values
		for the overall system is setup based on a speciation calculation performed while holding the unsteady
		variables constant at their respective initial values.

		\param ic initial concentration value for the unsteady object (mol/L)*/
	void Set_InitialValue(double ic);

	/// Set the maximum value of the unsteady variable to a given value max (mol/L)
	/** This function will be called internally to help bound the unsteady variable to reasonable
		maximum values. That maximum is usually based on the mass balances for the current non-linear
		iteration.

		\param max maximum allowable value for the unsteady variable (mol/L)*/
	void Set_MaximumValue(double max);

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

	void Set_TimeStep(double dt);						///< Set the time step for the current simulation

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
	
	bool haveForwardRef();								///< Function to return true if you have the forward reference rate
	bool haveReverseRef();								///< Function to return true if you have the reverse reference rate
	bool haveForward();									///< Function to return true if you have the forward rate
	bool haveReverse();									///< Function to return true if you have the reverse rate

	int Get_Species_Index();							///< Fetch the index of the Unsteady species
	double Get_Stoichiometric(int i);					///< Fetch the ith stoichiometric value
	double Get_Equilibrium();							///< Fetch the equilibrium constant (logK)
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

/// Adsorption Reaction Object
/** C++ Object to handle data and functions associated with forumlating adsorption equilibrium reactions
	in a aqueous mixture. Each unique surface in a system will require an instance of this structure.

*/
class AdsorptionReaction
{
public:
	AdsorptionReaction();								///< Default Constructor
	~AdsorptionReaction();								///< Default Destructor

	void Initialize_Object(MasterSpeciesList &List, int n); ///< Function to call the initialization of objects sequentially
	void Display_Info();								///< Display the adsorption reaction information (PLACE HOLDER)

	/// Modify the Deltas in the MassBalance Object
	/** This function will take a mass balance object as an argument and modify the deltas in that object to
		correct for how adsorption affects that particular mass balance. Since adsorption can effect multiple
		mass balances, this function must be called for each mass balance in the system.

		\param mbo reference to the MassBalance Object the adsorption is acting on*/
	void modifyDeltas(MassBalance &mbo);

	/// Find and set the adsorbed species indices for each reaction object
	/** This function searches through the Reaction objects in AdsorptionReaction to find the solid species
		and their indices to set that information in the adsorb_index structure. That information will be used
		later to approximate maximum capacities and equilibrium parameters for use in a modified extended Langmuir
		type expression. Function will return 0 if successful and -1 on a failure.*/
	int setAdsorbIndices();

	int checkAqueousIndices();							///< Function to check and report errors in the aqueous species indices

	/// Function to set the surface activity model and data pointer
	/** This function will setup the surface activity model based on the given pointer arguments. If no arguments
		are given, or are given as NULL, then the activity model will default to ideal solution assumption.*/
	void setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
							  const void *act_data);

	void setAqueousIndex(int rxn_i, int species_i);		///< Set the primary aqueous species index for the ith reaction

	/// Automatically sets the primary aqueous species index based on reactions
	/** This function will go through all species and all reactions in the adsorption object and automatically set the
		primary aqueous species index based on the stoicheometry of the reaction. It will also check and make sure that
		the primary aqueous index species appears opposite of the adsorbed species in the reactions. Note: This function
		assumes that the adsorbed indices have already been set. */
	int setAqueousIndexAuto();
	void setActivityEnum(int act);						///< Set the surface activity enum value
	void setMolarFactor(int rxn_i, double m);			///< Set the molar factor for the ith reaction (mol/mol)
	void setVolumeFactor(int i, double v);				///< Set the ith volume factor for the species list (cm^3/mol)
	void setAreaFactor(int i, double a);				///< Set the ith area factor for the species list (m^2/mol)
	void setSpecificArea(double a);						///< Set the specific area for the adsorbent (m^2/kg)
	void setSpecificMolality(double a);					///< Set the specific molality for the adsorbent (mol/kg)
	void setSurfaceCharge(double c);					///< Set the surface charge of the uncomplexed ligands
	void setTotalMass(double m);						///< Set the total mass of the adsorbent (kg)
	void setTotalVolume(double v);						///< Set the total volume of the system (L)
	void setAreaBasisBool(bool opt);					///< Set the basis boolean directly
	void setSurfaceChargeBool(bool opt);				///< Set the boolean for inclusion of surface charging 
	void setBasis(std::string option);					///< Set the basis of the adsorption problem from the given string arg
	void setAdsorbentName(std::string name);			///< Set the name of the adsorbent to the given string
	
	void setChargeDensityValue(double a);				///< Set the value of the charge density parameter to a (C/m^2)
	void setIonicStrengthValue(double a);				///< Set the value of the ionic strength parameter to a (mol/L)
	void setActivities(Matrix<double> &x);		///< Set the values of activities in the activity matrix

	void calculateAreaFactors();						///< Calculates the area factors used from the van der Waals volumes
	void calculateEquilibria(double T);					///< Calculates all equilibrium parameters as a function of temperature
	void setChargeDensity(const Matrix<double> &x);		///< Calculates and sets the current value of charge density
	void setIonicStrength(const Matrix<double> &x);		///< Calculates and sets the current value of ionic strength
	int callSurfaceActivity(const Matrix<double> &x);	///< Calls the activity model and returns an int flag for success or failure
	double calculateActiveFraction(const Matrix<double> &x);	///< Calculates the fraction of the surface that is active and available
	
	/// Function to calculate the surface charge density based on concentrations
	/** This function is used to calculate the surface charge density of the adsorbed species based on
		the charges and concentrations of the adsorbed species. The calculation is used to correct the
		adsorption equilibria constant based on a localized surface charge balance. This requires that 
		you know the molality of the uncomplexed ligand species on the surface, as well as the specific
		surface area for the adsorbent. 
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double calculateSurfaceChargeDensity( const Matrix<double> &x);

	/// Calculates the theoretical maximum capacity for adsorption  in reaction i
	/** This function is used to calculate the current maximum capacity of a species for a given
		adsorption reaction using the concentrations and activities of other species in the system.
		You must pass the index of the reaction of interest. The index of the species of interest
		is determined from the adsorb_index object. Note: This is only true if the stoicheometry for
		the adsorbed species is 1.

		\param i index of the reaction of interest for the adsorption object*/
	double calculateLangmuirMaxCapacity(int i);

	/// Calculates the equivalent Langmuir isotherm equilibrium parameter
	/** This function will take in the current aqueous activities and calculate an effective
		Langmuir adsorption parameter for use in determining the adsorption in the system. It
		uses the system temperature as well to calculate equilibrium. Note: This is only true
		if the stoicheometry for the adsorbed species is 1.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param i index of the reaction of interest for the adsorption object*/
	double calculateLangmuirEquParam(const Matrix<double> &x, const Matrix<double> &gama, int i);

	/// Calculates the equivalent Langmuir adsorption by forming the Langmuir-like parameters
	/** This function will use the calculateLangmuirMaxCapacity and calculateLangmuirEquParam functions to
		approximate the adsorption of the ith reaction given the concentration of aqueous species, activities,
		and temperature. Note: This is only true if the stoicheometry for the adsorbed species is 1.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param i index of the reaction of interest for the adsorption object*/
	double calculateLangmuirAdsorption(const Matrix<double> &x, const Matrix<double> &gama, int i);

	/// Function calculates the Psi (electric surface potential) given a set of arguments
	/** This function will calculate the electric surface potential of the adsorbent under the current
		conditions of charge density, temperature, ionic strength, and relative permittivity.

		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless) */
	double calculatePsi(double sigma, double T, double I, double rel_epsilon);

	/// Function to calculate the net exchange of charges of the aqeous species involved in a given reaction
	/** This function will look at all aqueous species involved in the ith adsorption reaction and sum up
		their stoicheometries and charges to see what the net change in charge is caused by the adsorption
		of charged species in solution. It is then used to adjust or correct the equilibrium constant for
		the given adsorption reaction.

		\param i index of the reaction of interest for the adsorption object*/
	double calculateAqueousChargeExchange(int i);

	/// Function to calculate the correction term for the equilibrium parameter
	/** This function calculates the correction term that gets applied to the equilibrium parameter to
		correct for surface charge and charge accumulation/depletion effects. It will call the psi approximation
		and charge exchange functions, therefore it needs to have those functions arguments passed to it as well.

		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int i);

	/// Calculates the residual for the ith reaction in the system
	/** This function will provide a system residual for the ith reaction object involved in the Adsorption
		Reaction. The residual is fed into the SHARK solver to find the solution to solid and aqueous phase
		concentrations simultaneously. This function will also adjust the equilibrium parameter for the reaction
		

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i);

	Reaction& getReaction(int i);				///< Return reference to the ith reaction object in the adsorption object
	double getMolarFactor(int i);				///< Get the ith reaction's molar factor for adsorption (mol/mol)
	double getVolumeFactor(int i);				///< Get the ith volume factor (species not involved return zeros) (cm^3/mol)
	double getAreaFactor(int i);				///< Get the ith area factor (species not involved return zeros) (m^2/mol)
	double getActivity(int i);					///< Get the ith activity factor for the surface species
	double getSpecificArea();					///< Get the specific area of the adsorbent (m^2/kg) or (mol/kg)
	double getSpecificMolality();				///< Get the specific molality of the adsorbent (mol/kg)
	double getSurfaceCharge();					///< Get the surface charge of the adsorbent
	double getBulkDensity();					///< Calculate and return bulk density of adsorbent in system (kg/L)
	double getTotalMass();						///< Get the total mass of adsorbent in the system (kg)
	double getTotalVolume();					///< Get the total volume of the system (L)
	double getChargeDensity();					///< Get the value of the surface charge density (C/m^2)
	double getIonicStrength();					///< Get the value of the ionic strength of solution (mol/L)
	int getNumberRxns();						///< Get the number of reactions involved in the adsorption object
	int getAdsorbIndex(int i);					///< Get the index of the adsorbed species in the ith reaction
	int getAqueousIndex(int i);					///< Get the index of the primary aqueous species in the ith reaction
	int getActivityEnum();						///< Return the enum representing the choosen activity function
	bool isAreaBasis();							///< Returns true if we are in the Area Basis, False if in Molar Basis
	bool includeSurfaceCharge();				///< Returns true if we are considering surface charging during adsorption
	std::string getAdsorbentName();				///< Returns the name of the adsorbent as a string 

protected:
	MasterSpeciesList *List;					///< Pointer to the MasterSpeciesList object

	/// Pointer to a surface activity model
	/** This is a function pointer for a surface activity model. The function must accept the log of the
		surface concentrations as an argument (logq) and provide the activities for each species (activity).
		The pointer data is used to pass any additional arguments needed.

		\param logq matrix of the log (base 10) of surface concentrations of all species
		\param activity matrix of activity coefficients for all surface species (must be overriden)
		\param data pointer to a data structure needed to calculate activities*/
	int (*surface_activity) (const Matrix<double>& logq, Matrix<double> &activity, const void *data);

	const void *activity_data;					///< Pointer to the data structure needed for surface activities.
	int act_fun;								///< Enumeration of the activity function being used for the surface phase
	std::vector<double> area_factors;			///< List of the van der Waals areas associated with surface species (m^2/mol)
	std::vector<double> volume_factors;			///< List of the van der Waals volumes of each surface species (cm^3/mol)
	std::vector<int> adsorb_index;				///< List of the indices for the adsorbed species in the reactions
	std::vector<int> aqueous_index;				///< List of the indices for the primary aqueous species in the reactions
	std::vector<double> molar_factor;			///< List of the number of ligands needed to form one mole of adsorption in each reaction
	Matrix<double> activities;					///< List of the activities calculated by the activity model
	double specific_area;						///< Specific surface area of the adsorbent (m^2/kg)
	double specific_molality;					///< Specific molality of the adsorbent - moles of ligand per kg sorbent (mol/kg)
	double surface_charge;						///< Charge of the uncomplexed surface ligand species
	double total_mass;							///< Total mass of the adsorbent in the system (kg)
	double total_volume;						///< Total volume of the system (L)
	double ionic_strength;						///< Ionic Strength of the system used to adjust equilibria constants (mol/L)
	double charge_density;						///< Surface charge density of the adsorbent used to adjust equilbria (C/m^2)
	int num_rxns;								///< Number of reactions involved in the adsorption equilibria
	bool AreaBasis;								///< True = Adsorption on an area basis, False = Adsorption on a ligand basis
	bool IncludeSurfCharge;						///< True = Includes surface charging corrections, False = Does not consider surface charge
	std::string adsorbent_name;					///< Name of the adsorbent for this object

private:
	std::vector<Reaction> ads_rxn;				///< List of reactions involved with adsorption

};

/// Unsteady Adsorption Reaction Object
/** C++ Object to handle data and functions associated with forumlating unsteady adsorption reactions
	in a aqueous mixture. Each unique surface in a system will require an instance of this structure.

*/
class UnsteadyAdsorption : AdsorptionReaction
{
public:
	UnsteadyAdsorption();								///< Default Constructor
	~UnsteadyAdsorption();								///< Default Destructor

	void Initialize_Object(MasterSpeciesList &List, int n); ///< Function to call the initialization of objects sequentially
	void Display_Info();								///< Display the adsorption reaction information (PLACE HOLDER)

	/// Modify the Deltas in the MassBalance Object
	/** This function will take a mass balance object as an argument and modify the deltas in that object to
		correct for how adsorption affects that particular mass balance. Since adsorption can effect multiple
		mass balances, this function must be called for each mass balance in the system.

		\param mbo reference to the MassBalance Object the adsorption is acting on*/
	void modifyDeltas(MassBalance &mbo);

	/// Find and set the adsorbed species indices for each reaction object
	/** This function searches through the Reaction objects in UnsteadyAdsorption to find the solid species
		and their indices to set that information in the adsorb_index structure. That information will be used
		later to approximate maximum capacities and equilibrium parameters for use in a modified extended Langmuir
		type expression. Function will return 0 if successful and -1 on a failure.*/
	int setAdsorbIndices();

	int checkAqueousIndices();							///< Function to check and report errors in the aqueous species indices

	/// Function to set the surface activity model and data pointer
	/** This function will setup the surface activity model based on the given pointer arguments. If no arguments
		are given, or are given as NULL, then the activity model will default to ideal solution assumption.*/
	void setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
							  const void *act_data);

	void setAqueousIndex(int rxn_i, int species_i);		///< Set the primary aqueous species index for the ith reaction

	/// Automatically sets the primary aqueous species index based on reactions
	/** This function will go through all species and all reactions in the adsorption object and automatically set the
		primary aqueous species index based on the stoicheometry of the reaction. It will also check and make sure that
		the primary aqueous index species appears opposite of the adsorbed species in the reactions. Note: This function
		assumes that the adsorbed indices have already been set. */
	int setAqueousIndexAuto();
	void setActivityEnum(int act);						///< Set the surface activity enum value
	void setMolarFactor(int rxn_i, double m);			///< Set the molar factor for the ith reaction (mol/mol)
	void setVolumeFactor(int i, double v);				///< Set the ith volume factor for the species list (cm^3/mol)
	void setAreaFactor(int i, double a);				///< Set the ith area factor for the species list (m^2/mol)
	void setSpecificArea(double a);						///< Set the specific area for the adsorbent (m^2/kg)
	void setSpecificMolality(double a);					///< Set the specific molality for the adsorbent (mol/kg)
	void setSurfaceCharge(double c);					///< Set the surface charge of the uncomplexed ligands
	void setTotalMass(double m);						///< Set the total mass of the adsorbent (kg)
	void setTotalVolume(double v);						///< Set the total volume of the system (L)
	void setAreaBasisBool(bool opt);					///< Set the basis boolean directly
	void setSurfaceChargeBool(bool opt);				///< Set the boolean for inclusion of surface charging
	void setBasis(std::string option);					///< Set the basis of the adsorption problem from the given string arg
	void setAdsorbentName(std::string name);			///< Set the name of the adsorbent to the given string

	void updateActivities();							///< Set the old activities as the new activities before doing next time step
	void calculateAreaFactors();						///< Calculates the area factors used from the van der Waals volumes
	void calculateEquilibria(double T);					///< Calculates all equilibrium parameters as a function of temperature
	void calculateRates(double T);						///< Calculates all reaction rate parameters as a function of temperature
	void setChargeDensity(const Matrix<double> &x);		///< Calculates and sets the current value of charge density
	void setIonicStrength(const Matrix<double> &x);		///< Calculates and sets the current value of ionic strength
	int callSurfaceActivity(const Matrix<double> &x);	///< Calls the activity model and returns an int flag for success or failure
	double calculateActiveFraction(const Matrix<double> &x);	///< Calculates the fraction of the surface that is active and available

	/// Function to calculate the surface charge density based on concentrations
	/** This function is used to calculate the surface charge density of the adsorbed species based on
		the charges and concentrations of the adsorbed species. The calculation is used to correct the
		adsorption equilibria constant based on a localized surface charge balance. This requires that
		you know the molality of the uncomplexed ligand species on the surface, as well as the specific
		surface area for the adsorbent.

		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double calculateSurfaceChargeDensity( const Matrix<double> &x);

	/// Function calculates the Psi (electric surface potential) given a set of arguments
	/** This function will calculate the electric surface potential of the adsorbent under the current
		conditions of charge density, temperature, ionic strength, and relative permittivity.

		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless) */
	double calculatePsi(double sigma, double T, double I, double rel_epsilon);

	/// Function to calculate the net exchange of charges of the aqeous species involved in a given reaction
	/** This function will look at all aqueous species involved in the ith adsorption reaction and sum up
		their stoicheometries and charges to see what the net change in charge is caused by the adsorption
		of charged species in solution. It is then used to adjust or correct the equilibrium constant for
		the given adsorption reaction.

		\param i index of the reaction of interest for the adsorption object*/
	double calculateAqueousChargeExchange(int i);

	/// Function to calculate the correction term for the equilibrium parameter
	/** This function calculates the correction term that gets applied to the equilibrium parameter to
		correct for surface charge and charge accumulation/depletion effects. It will call the psi approximation
		and charge exchange functions, therefore it needs to have those functions arguments passed to it as well.

		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int i);

	/// Calculates the residual for the ith reaction in the system
	/** This function will provide a system residual for the ith reaction object involved in the Adsorption
		Reaction. The residual is fed into the SHARK solver to find the solution to solid and aqueous phase
		concentrations simultaneously. This function will also adjust the equilibrium parameter for the reaction


		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i);

	/// Calculates the unsteady residual for the ith reaction in the system
	/** This function will provide a system residual for the ith reaction object involved in the Unsteady Adsorption
		Reaction. The residual is fed into the SHARK solver to find the solution to solid and aqueous phase
		concentrations simultaneously. This function will also adjust the equilibrium parameter for the reaction


		\param x_new matrix of the current log(C) concentration values at the current non-linear step
		\param gama_new matrix of current activity coefficients for each species at the current non-linear step
		\param x_old matrix of the old log(C) concentration values at the current non-linear step
		\param gama_old matrix of old activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double Eval_Residual(const Matrix<double> &x_new, const Matrix<double> &x_old, const Matrix<double> &gama_new, const Matrix<double> &gama_old, double T, double rel_perm, int i);

	/// Function to calculate the explicit or implicit rate of reaction
	/** This function will calculate the rate/extent of the unsteady adsorption reaction given the log(C) concentrations
	 * and aqueous activities, as well as temperature and permittivity. The temperature and permittivity are used to make
	 * surface charge corrections to the equilibria and rate constants.
	 *
	 	\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param i index of the reaction of interest for the adsorption object
	 */
	double Eval_ReactionRate(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i);

	/// Calculate the unsteady residual for initial conditions
	/** Setting the intial conditions for all variables in the system requires a speciation calculation.
		However, we want the unsteady variables to be set to their respective initial conditions. Using this
		residual function imposes an equality constraint on those non-linear, unsteady variables allowing the
		rest of the speciation problem to be solved via PJFNK iterations.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param i index of the reaction of interest for the adsorption object*/
	double Eval_IC_Residual(const Matrix<double> &x, int i);

	/// Return an approximate explicit solution to our unsteady adsorption variable (mol/kg)
	/** This function will approximate the concentration of the unsteady variables based on an explicit time
		discretization. The purpose of this function is to try to provide the PJFNK method with a good initial
		guess for the values of the non-linear, unsteady variables. If we do not provide a good initial guess
		to these variables, then the PJFNK method may not converge to the correct solution, because the unsteady
		problem is the most difficult to solve.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double Explicit_Eval(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i);

	UnsteadyReaction& getReaction(int i);		///< Return reference to the ith reaction object in the adsorption object
	double getMolarFactor(int i);				///< Get the ith reaction's molar factor for adsorption (mol/mol)
	double getVolumeFactor(int i);				///< Get the ith volume factor (species not involved return zeros) (cm^3/mol)
	double getAreaFactor(int i);				///< Get the ith area factor (species not involved return zeros) (m^2/mol)
	double getActivity(int i);					///< Get the ith activity factor for the surface species
	double getOldActivity(int i);				///< Get the ith old activity factor for the surface species
	double getSpecificArea();					///< Get the specific area of the adsorbent (m^2/kg) or (mol/kg)
	double getSpecificMolality();				///< Get the specific molality of the adsorbent (mol/kg)
	double getSurfaceCharge();					///< Get the surface charge of the adsorbent
	double getBulkDensity();					///< Calculate and return bulk density of adsorbent in system (kg/L)
	double getTotalMass();						///< Get the total mass of adsorbent in the system (kg)
	double getTotalVolume();					///< Get the total volume of the system (L)
	double getChargeDensity();					///< Get the value of the surface charge density (C/m^2)
	double getIonicStrength();					///< Get the value of the ionic strength of solution (mol/L)
	int getNumberRxns();						///< Get the number of reactions involved in the adsorption object
	int getAdsorbIndex(int i);					///< Get the index of the adsorbed species in the ith reaction
	int getAqueousIndex(int i);					///< Get the index of the primary aqueous species in the ith reaction
	int getActivityEnum();						///< Return the enum representing the choosen activity function
	bool isAreaBasis();							///< Returns true if we are in the Area Basis, False if in Molar Basis
	bool includeSurfaceCharge();				///< Returns true if we are considering surface charging during adsorption
	std::string getAdsorbentName();				///< Returns the name of the adsorbent as a string

protected:
	Matrix<double> activities_old;				///< List of the old activities calculated by the activity model

private:
	std::vector<UnsteadyReaction> ads_rxn;				///< List of reactions involved with adsorption

};

/// Multi-ligand Adsorption Reaction Object
/** C++ Object to handle data and functions associated with forumlating multi-ligand adsorption reactions
	in a aqueous mixture. Each unique surface in a system will require an instance of this structure. This
	object is made from a vector of AdsorptionReaction objects, but differentiate between different ligands
	that exist on the surface.
 
 */
class MultiligandAdsorption
{
public:
	MultiligandAdsorption();								///< Default Constructor
	~MultiligandAdsorption();								///< Default Destructor
	
	/// Function to call the initialization of objects sequentially
	/** Function will initialize each ligand adsorption object. 
	 
		\param List reference to MasterSpeciesList object
		\param l number of ligands on the surface
		\param n number of reactions for each ligand (ligands must be correctly indexed)*/
	void Initialize_Object(MasterSpeciesList &List, int l, std::vector<int> n);
	
	/// Modify the Deltas in the MassBalance Object
	/** This function will take a mass balance object as an argument and modify the deltas in that object to
		correct for how adsorption affects that particular mass balance. Since adsorption can effect multiple
		mass balances, this function must be called for each mass balance in the system.
	 
		\param mbo reference to the MassBalance Object the adsorption is acting on*/
	void modifyDeltas(MassBalance &mbo);
	
	/// Find and set the adsorbed species indices for each reaction object in each ligand object
	/** This function searches through the Reaction objects in AdsorptionReaction to find the solid species
		and their indices to set that information in the adsorb_index structure. That information will be used
		later to approximate maximum capacities and equilibrium parameters for use in a modified extended Langmuir
		type expression. Function will return 0 if successful and -1 on a failure.*/
	int setAdsorbIndices();
	
	/// Function to check and report errors in the aqueous species indices
	int checkAqueousIndices();
	
	/// Function to set the surface activity model and data pointer
	/** This function will setup the surface activity model based on the given pointer arguments. If no arguments
		are given, or are given as NULL, then the activity model will default to ideal solution assumption.*/
	void setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
							  const void *act_data);
	
	/// Automatically sets the primary aqueous species index based on reactions for each ligand
	/** This function will go through all species and all reactions in each adsorption object and automatically set the
		primary aqueous species index based on the stoicheometry of the reaction. It will also check and make sure that
		the primary aqueous index species appears opposite of the adsorbed species in the reactions. Note: This function
		assumes that the adsorbed indices have already been set. */
	int setAqueousIndexAuto();
	
	void setActivityEnum(int act);					///< Set the activity enum to the value of act
	void setMolarFactor(int ligand, int rxn, double m);///< Set the molar factor for the rxn reaction of the ligand ligand to a value of m
	void setVolumeFactor(int i, double v);			///< Set all ith volume factors for the species list (cm^3/mol)
	void setAreaFactor(int i, double a);			///< Set all ith area factors for the species list (m^2/mol)
	void setSpecificMolality(int ligand, double a);	///< Set the specific molality for the ligand (mol/kg)
	void setSurfaceCharge(int ligand, double c);	//< Set the surface charge of the uncomplexed ligand
	void setAdsorbentName(std::string name);		///< Set the name of the adsorbent material or particle
	void setLigandName(int i, std::string name);	///< Set the name of the ith ligand
	void setSpecificArea(double area);				///< Set the specific area of the adsorbent
	void setTotalMass(double mass);					///< Set the mass of the adsorbent
	void setTotalVolume(double volume);				///< Set the total volume of the system
	void setSurfaceChargeBool(bool opt);			///< Set the surface charge boolean
	void setElectricPotential(double a);			///< Set the surface electric potential
	
	void calculateAreaFactors();						///< Calculates the area factors used from the van der Waals volumes
	void calculateEquilibria(double T);					///< Calculates all equilibrium parameters as a function of temperature
	void setChargeDensity(const Matrix<double> &x);		///< Calculates and sets the current value of charge density
	void setIonicStrength(const Matrix<double> &x);		///< Calculates and sets the current value of ionic strength
	int callSurfaceActivity(const Matrix<double> &x);	///< Calls the activity model and returns an int flag for success or failure
	
	/// Function calculates the Psi (electric surface potential) given a set of arguments
	/** This function will calculate the electric surface potential of the adsorbent under the current
		conditions of charge density, temperature, ionic strength, and relative permittivity.
	 
		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless) */
	void calculateElecticPotential(double sigma, double T, double I, double rel_epsilon);
	
	/// Function to calculate the correction term for the equilibrium parameter
	/** This function calculates the correction term that gets applied to the equilibrium parameter to
		correct for surface charge and charge accumulation/depletion effects. It will call the psi approximation
		and charge exchange functions, therefore it needs to have those functions arguments passed to it as well.
	 
		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless)
		\param rxn index of the reaction of interest for the adsorption object
		\param ligand index of the ligand of interest for the adsorption object*/
	double calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int rxn, int ligand);
	
	/// Calculates the residual for the ith reaction and lth ligand in the system
	/** This function will provide a system residual for the ith reaction object involved in the lth ligand's Adsorption
		Reaction object. The residual is fed into the SHARK solver to find the solution to solid and aqueous phase
		concentrations simultaneously. This function will also adjust the equilibrium parameter for the reaction
		
	 
		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param rxn index of the reaction of interest for the adsorption object
		\param ligand index of the ligand of interest for the adsorption object*/
	double Eval_Residual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int rxn, int ligand);
	
	AdsorptionReaction& getAdsorptionObject(int i);	///< Return reference to the adsortpion object corresponding to ligand i
	int getNumberLigands();							///< Get the number of ligands involved with the surface
	int getActivityEnum();							///< Get the value of the activity enum set by user
	double getActivity(int i);						///< Get the ith activity coefficient from the matrix object
	double getSpecificArea();						///< Get the specific area of the adsorbent (m^2/kg) or (mol/kg)
	double getBulkDensity();						///< Calculate and return bulk density of adsorbent in system (kg/L)
	double getTotalMass();							///< Get the total mass of adsorbent in the system (kg)
	double getTotalVolume();						///< Get the total volume of the system (L)
	double getChargeDensity();						///< Get the value of the surface charge density (C/m^2)
	double getIonicStrength();						///< Get the value of the ionic strength of solution (mol/L)
	double getElectricPotential();					///< Get the value of the electric surface potential (V)
	bool includeSurfaceCharge();					///< Returns true if we are considering surface charging during adsorption
	std::string getLigandName(int i);				///< Get the name of the ligand object indexed by i
	std::string getAdsorbentName();					///< Get the name of the adsorbent
	
protected:
	MasterSpeciesList *List;						///< Pointer to the MasterSpeciesList object
	int num_ligands;								///< Number of different ligands to consider
	std::string adsorbent_name;						///< Name of the adsorbent
	
	/// Pointer to a surface activity model
	/** This is a function pointer for a surface activity model. The function must accept the log of the
		surface concentrations as an argument (logq) and provide the activities for each species (activity).
		The pointer data is used to pass any additional arguments needed.
	 
		\param logq matrix of the log (base 10) of surface concentrations of all species
		\param activity matrix of activity coefficients for all surface species (must be overriden)
		\param data pointer to a data structure needed to calculate activities*/
	int (*surface_activity) (const Matrix<double>& logq, Matrix<double> &activity, const void *data);
	
	const void *activity_data;					///< Pointer to the data structure needed for surface activities.
	int act_fun;								///< Enumeration to represent the choosen surface activity function
	Matrix<double> activities;					///< List of the activities calculated by the activity model
	double specific_area;						///< Specific surface area of the adsorbent (m^2/kg)
	double total_mass;							///< Total mass of the adsorbent in the system (kg)
	double total_volume;						///< Total volume of the system (L)
	double ionic_strength;						///< Ionic Strength of the system used to adjust equilibria constants (mol/L)
	double charge_density;						///< Surface charge density of the adsorbent used to adjust equilbria (C/m^2)
	double electric_potential;					///< Electric surface potential of the adsorbent used to adjust equilibria (V)
	bool IncludeSurfCharge;						///< True = Includes surface charging corrections, False = Does not consider surface charge
	
private:
	std::vector<AdsorptionReaction> ligand_obj;		///< List of the ligands and reactions they have on the surface
	
};

/// Chemisorption Reaction Object
/** C++ Object to handle data and functions associated with forumlating adsorption equilibrium reactions
	in a aqueous mixture based on chemisorption mechanisms. Each unique surface in a system will require 
	an instance of this structure. This is very similar to AdsorptionReaction, however, this will include
	a site balance residual that will allow us to consider protonation and deprotonation of the ligands.
 
 */
class ChemisorptionReaction : AdsorptionReaction
{
public:
	ChemisorptionReaction();								///< Default Constructor
	~ChemisorptionReaction();								///< Default Destructor
	
	void Initialize_Object(MasterSpeciesList &List, int n); ///< Function to call the initialization of objects sequentially
	void Display_Info();									///< Display the adsorption reaction information (PLACE HOLDER)
	
	/// Modify the Deltas in the MassBalance Object
	/** This function will take a mass balance object as an argument and modify the deltas in that object to
		correct for how adsorption affects that particular mass balance. Since adsorption can effect multiple
		mass balances, this function must be called for each mass balance in the system.
	 
		\param mbo reference to the MassBalance Object the adsorption is acting on*/
	void modifyMBEdeltas(MassBalance &mbo);
	
	/// Find and set the adsorbed species indices for each reaction object
	/** This function searches through the Reaction objects in ChemisorptionReaction to find the adsorbed species
		and their indices to set that information in the adsorb_index structure. Function will return 0 if successful 
		and -1 on a failure.*/
	int setAdsorbIndices();
	
	/// Find and set the ligand species index
	/** This function searches through the Reaction objects in ChemisorptionReaction to find the ligand species
		and its index to set that information in the ligand_index structure. Function will return 0 if successful
		and -1 on a failure.*/
	int setLigandIndex();
	
	/// Function to set the surface activity model and data pointer
	/** This function will setup the surface activity model based on the given pointer arguments. If no arguments
		are given, or are given as NULL, then the activity model will default to ideal solution assumption.*/
	void setActivityModelInfo( int (*act) (const Matrix<double>& logq, Matrix<double> &activity, const void *data),
							  const void *act_data);
	
	void setActivityEnum(int act);						///< Set the surface activity enum value
	void setDelta(int i, double v);						///< Set the ith delta factor for the site balance
	void setVolumeFactor(int i, double v);				///< Set the ith volume factor for the species list (cm^3/mol)
	void setAreaFactor(int i, double a);				///< Set the ith area factor for the species list (m^2/mol)
	void setSpecificArea(double a);						///< Set the specific area for the adsorbent (m^2/kg)
	void setSpecificMolality(double a);					///< Set the specific molality for the adsorbent (mol/kg)
	void setTotalMass(double m);						///< Set the total mass of the adsorbent (kg)
	void setTotalVolume(double v);						///< Set the total volume of the system (L)
	void setSurfaceChargeBool(bool opt);				///< Set the boolean for inclusion of surface charging
	void setAdsorbentName(std::string name);			///< Set the name of the adsorbent to the given string
	
	void setChargeDensityValue(double a);				///< Set the value of the charge density parameter to a (C/m^2)
	void setIonicStrengthValue(double a);				///< Set the value of the ionic strength parameter to a (mol/L)
	void setActivities(Matrix<double> &x);				///< Set the values of activities in the activity matrix
	
	void calculateAreaFactors();						///< Calculates the area factors used from the van der Waals volumes
	void calculateEquilibria(double T);					///< Calculates all equilibrium parameters as a function of temperature
	void setChargeDensity(const Matrix<double> &x);		///< Calculates and sets the current value of charge density
	void setIonicStrength(const Matrix<double> &x);		///< Calculates and sets the current value of ionic strength
	int callSurfaceActivity(const Matrix<double> &x);	///< Calls the activity model and returns an int flag for success or failure
	
	/// Function to calculate the surface charge density based on concentrations
	/** This function is used to calculate the surface charge density of the adsorbed species based on
		the charges and concentrations of the adsorbed species. The calculation is used to correct the
		adsorption equilibria constant based on a localized surface charge balance. This requires that
		you know the molality of the uncomplexed ligand species on the surface, as well as the specific
		surface area for the adsorbent.
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double calculateSurfaceChargeDensity( const Matrix<double> &x);
	
	/// Function calculates the Psi (electric surface potential) given a set of arguments
	/** This function will calculate the electric surface potential of the adsorbent under the current
		conditions of charge density, temperature, ionic strength, and relative permittivity.
	 
		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless) */
	double calculateElecticPotential(double sigma, double T, double I, double rel_epsilon);
	
	/// Function to calculate the net exchange of charges of the aqeous species involved in a given reaction
	/** This function will look at all aqueous species involved in the ith adsorption reaction and sum up
		their stoicheometries and charges to see what the net change in charge is caused by the adsorption
		of charged species in solution. It is then used to adjust or correct the equilibrium constant for
		the given adsorption reaction.
	 
		\param i index of the reaction of interest for the adsorption object*/
	double calculateAqueousChargeExchange(int i);
	
	/// Function to calculate the correction term for the equilibrium parameter
	/** This function calculates the correction term that gets applied to the equilibrium parameter to
		correct for surface charge and charge accumulation/depletion effects. It will call the psi approximation
		and charge exchange functions, therefore it needs to have those functions arguments passed to it as well.
	 
		\param sigma charge density of the surface (C/m^2)
		\param T temperature of the system in question (K)
		\param I ionic strength of the medium the surface is in (mol/L)
		\param rel_epsilon relative permittivity of the medium (Unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double calculateEquilibriumCorrection(double sigma, double T, double I, double rel_epsilon, int i);
	
	/// Calculates the residual for the ith reaction in the system
	/** This function will provide a system residual for the ith reaction object involved in the Adsorption
		Reaction. The residual is fed into the SHARK solver to find the solution to solid and aqueous phase
		concentrations simultaneously. This function will also adjust the equilibrium parameter for the reaction
		
	 
		\param x matrix of the log(C) concentration values at the current non-linear step
		\param gama matrix of activity coefficients for each species at the current non-linear step
		\param T temperature of the system in question (K)
		\param rel_perm relative permittivity of the media (unitless)
		\param i index of the reaction of interest for the adsorption object*/
	double Eval_RxnResidual(const Matrix<double> &x, const Matrix<double> &gama, double T, double rel_perm, int i);
	
	/// Calculates the residual for the overall site balance
	/** This function will provide a system residual for the site/ligand balance for the Chemisorption
		Reaction object. The residual is fed into the SHARK solver to find the solution to solid and aqueous phase
		concentrations simultaneously.
		
	 
		\param x matrix of the log(C) concentration values at the current non-linear step*/
	double Eval_SiteBalanceResidual(const Matrix<double> &x);
	
	Reaction& getReaction(int i);				///< Return reference to the ith reaction object in the adsorption object
	double getDelta(int i);						///< Get the ith delta factor for the site balance
	double getVolumeFactor(int i);				///< Get the ith volume factor (species not involved return zeros) (cm^3/mol)
	double getAreaFactor(int i);				///< Get the ith area factor (species not involved return zeros) (m^2/mol)
	double getActivity(int i);					///< Get the ith activity factor for the surface species
	double getSpecificArea();					///< Get the specific area of the adsorbent (m^2/kg) or (mol/kg)
	double getSpecificMolality();				///< Get the specific molality of the adsorbent (mol/kg)
	double getBulkDensity();					///< Calculate and return bulk density of adsorbent in system (kg/L)
	double getTotalMass();						///< Get the total mass of adsorbent in the system (kg)
	double getTotalVolume();					///< Get the total volume of the system (L)
	double getChargeDensity();					///< Get the value of the surface charge density (C/m^2)
	double getIonicStrength();					///< Get the value of the ionic strength of solution (mol/L)
	int getNumberRxns();						///< Get the number of reactions involved in the adsorption object
	int getAdsorbIndex(int i);					///< Get the index of the adsorbed species in the ith reaction
	int getLigandIndex();						///< Get the index of the ligand species
	int getActivityEnum();						///< Return the enum representing the choosen activity function
	bool includeSurfaceCharge();				///< Returns true if we are considering surface charging during adsorption
	std::string getAdsorbentName();				///< Returns the name of the adsorbent as a string
	
protected:
	int ligand_index;						///< Index of the ligand for all reactions
	std::vector<double> Delta;				///< Vector of weights (i.e., deltas) used in the site balance
	
private:
	std::vector<Reaction> ads_rxn;				///< List of reactions involved with adsorption
	
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
typedef enum {IDEAL, DAVIES, DEBYE_HUCKEL, SIT, PITZER} valid_act;

/// Enumeration for the list of valid surface activity models for non-ideal adsorption
/** \note We had to create an IDEAL_ADS option to replace the IDEAL enum already in use for non-ideal solution
	or aqueous phases. (ADS => adsorption) */
typedef enum {IDEAL_ADS, FLORY_HUGGINS, UNIQUAC_ACT} valid_surf_act;

/// Data structure for SHARK simulations
/** C-style object holding data and function pointers associated with solving aqueous speciation and reaction
	kinetics. This object couples all other objects available in shark.h in order to provide residual calculations
	for each individual function that makes up the overall system model. Those residuals are brought together inside
	the residual function and fed into the lark.h PJFNK solver routine. That solver then attempts to find a solution
	to all non-linear variables simultaneously. Any function or data pointers in this structure can be overriden
	to change how you interface with and solve the problem. Users may also provide a set of custom residual functions
	through the "OtherList" vector object. Those residual function must all have the same format. */
typedef struct SHARK_DATA
{
	MasterSpeciesList MasterList;							///< Master List of species object
	std::vector<Reaction> ReactionList;						///< Equilibrium reaction objects
	std::vector<MassBalance> MassBalanceList;				///< Mass balance objects
	std::vector<UnsteadyReaction> UnsteadyList;				///< Unsteady Reaction objects
	std::vector<AdsorptionReaction> AdsorptionList;			///< Equilibrium Adsorption Reaction Objects
	std::vector<UnsteadyAdsorption> UnsteadyAdsList;		///< Unsteady Adsorption Reaction Objects
	std::vector<MultiligandAdsorption> MultiAdsList;		///< Multiligand Adsorptioin Objects
	std::vector<ChemisorptionReaction> ChemisorptionList;	///< Chemisorption Reaction objects

	/// Array of Other Residual functions to be defined by user
	/** This list of function pointers can be declared and set up by the user in order to add to or change
		the behavior of the SHARK system. Each one must be declared setup individually by the user. They will
		be called by the shark_residual function when needed. Alternatively, the user is free to provide their
		own shark_residual function for the overall system.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param shark_dat pointer to the SHARK_DATA data structure
		\param data pointer to a user defined data structure that is used to evaluate this residual*/
	std::vector<
		double (*) (const Matrix<double> &x,
					SHARK_DATA *shark_dat,
					const void *data) > OtherList;

	int numvar;										///< Total number of functions and species
	int num_ssr;									///< Number of steady-state reactions
	int num_mbe;									///< Number of mass balance equations
	int num_usr = 0;								///< Number of unsteady-state reactions
	int num_ssao = 0;								///< Number of steady-state adsorption objects
	int num_usao = 0;								///< Number of unsteady adsorption objects
	int num_multi_ssao = 0;							///< Number of multiligand steady-state adsorption objects
	int num_sschem = 0;								///< Number of steady-state chemisorption objects
	std::vector<int> num_ssar;						///< List of the numbers of reactions in each adsorption object
	std::vector<int> num_usar;						///< List of the numbers of reactions in each unsteady adsorption object
	std::vector<int> num_sschem_rxns;				///< List of the numbers of reactions in each steady-state chemisorption object
	std::vector< std::vector<int> > num_multi_ssar; ///< List of all multiligand objects -> List of ligands and rxns of that ligand
	std::vector<std::string> ss_ads_names;			///< List of the steady-state adsorbent object names
	std::vector<std::string> us_ads_names;			///< List of the unsteady adsorption object names
	std::vector<std::string> ss_chem_names;			///< List of the steady-state chemisorption object names
	std::vector< std::vector<std::string> > ssmulti_names;	///< List of the names of the ligands in each multiligand object
	int num_other = 0;								///< Number of other functions to be used (default is always 0)
	int act_fun = IDEAL;							///< Flag denoting the activity function to use (default is IDEAL)
	int reactor_type = BATCH;						///< Flag denoting the type of reactor considered for the system (default is BATCH)
	int totalsteps = 0;								///< Total number of iterations
	int totalcalls = 0;								///< Total number of residual function calls
	int timesteps = 0;								///< Number of time steps taken to complete simulation
	int pH_index = -1;								///< Contains the index of the pH variable (set internally)
	int pOH_index = -1;								///< Contains the index of the pOH variable (set internally)

	double simulationtime = 0.0;					///< Time to simulate unsteady reactions for (default = 0.0 hrs)
	double dt = 0.1;								///< Time step size (hrs)
	double dt_min = sqrt(DBL_EPSILON);				///< Minimum allowable step size
	double dt_max = 744.0;							///< Maximum allowable step size (~1 month in time)
	double t_out = 0.0;								///< Time increment by which file output is made (default = print all time steps)
	double t_count = 0.0;							///< Running count of time increments
	double time = 0.0;								///< Current value of time (starts from t = 0.0 hrs)
	double time_old = 0.0;							///< Previous value of time (start from t = 0.0 hrs)
	double pH = 7.0;								///< Value of pH if needed (default = 7)
	double pH_step = 0.5;							///< Value by which to increment pH when doing a speciation curve (default = 0.5)
	double start_temp = 277.15;						///< Value of the starting temperature used for Temperature Curves (default = 277.15 K)
	double end_temp = 323.15;						///< Value of the ending temperature used for Temperature Curves (default = 323.15 K)
	double temp_step = 10.0;						///< Size of the step changes to use for Temperature Curves (default = 10.0 K);
	double volume = 1.0;							///< Volume of the domain in liters (default = 1 L)
	double flow_rate = 1.0;							///< Flow rate in the reactor in L/hr (default = 1 L/hr)
	double xsec_area = 1.0;							///< Cross sectional area of the reactor in m^2 (default = 1 m^2)
	double Norm = 0.0;								///< Current value of euclidean norm in solution

	double dielectric_const = 78.325;				///< Dielectric constant used in many activity models (default: water = 78.325 (1/K))
	double relative_permittivity = 80.1;			///< Relative permittivity of the medium (default: water = 80.1 (-))
	double temperature = 298.15;					///< Solution temperature (default = 25 oC or 298.15 K)
	double ionic_strength = 0.0;					///< Solution ionic strength in Molar (calculated internally)

	bool steadystate = true;						///< True = solve steady problem; False = solve transient problem
	bool ZeroInitialSolids = false;					///< True = no solids or adsorption initially in the reactor
	bool TimeAdaptivity = false;					///< True = solve using variable time step
	bool const_pH = false;							///< True = set pH to a constant; False = solve for pH
	bool SpeciationCurve = false;					///< True = runs a series of constant pH steady-state problems to produce curves
	bool TemperatureCurve = false;					///< True = runs a series of constant temperature steady-state problmes to produce curves
	bool Console_Output = true;						///< True = display output to console
	bool File_Output = false;						///< True = write output to a file
	bool Contains_pH = false;						///< True = system contains pH as a variable (set internally)
	bool Contains_pOH = false;						///< True = system contains pOH as a variable (set internally)
	bool Converged = false;							///< True = system converged within tolerance
	bool LocalMin = true;							///< True = allow the system to settle for a local minimum if tolerance not reached

	Matrix<double> X_old;							///< Solution vector for old time step - log(C)
	Matrix<double> X_new;							///< Solution vector for current time step - log(C)
	Matrix<double> Conc_old;						///< Concentration vector for old time step - 10^x
	Matrix<double> Conc_new;						///< Concentration vector for current time step - 10^x
	Matrix<double> activity_new;					///< Activity matrix for current time step
	Matrix<double> activity_old;					///< Activity matrix from prior time step

	//Function pointers for activity and residual evaluations

	/// Function pointer to evaluate activity coefficients
	/** This function pointer is called within the shark_residual function to calculate and modify the activity_new
		matrix entries. When using the SHARK default options, this function pointer will be automatically set to a
		cooresponding activity function for the list of valid functions from the valid_act enum. User may override
		this function pointer if they desire. Must be overriden after calling the setup function.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param F matrix of activity coefficients that are to be altered by this function
		\param data pointer to a data structure needed to evaluate the activity model*/
	int (*EvalActivity) (const Matrix<double>& x, Matrix<double> &F, const void *data);

	/// Function pointer to evaluate all residuals in the system
	/** This function will be fed into the PJFNK solver (see lark.h) to solve the non-linear system of equations.
		By default, this pointer will be the shark_residual function (see below). However, the user may override
		the function and provide their own residuals for the PJFNK solver to operate on.

		\param x matrix of the log(C) concentration values at the current non-linear step
		\param F matrix of residuals that are to be altered from the functions in the system
		\param data pointer to a data structure needed to evaluate the activity model*/
	int (*Residual) (const Matrix<double>& x, Matrix<double> &F, const void *data);

	/// Function pointer to form a linear preconditioning operation for the Jacobian
	/** This function will be fed into the linear solver used for each non-linear step in PJFNK (see lark.h). By
		default, we cannot provide any linear preconditioner, because we do not know the form or sparcity of the
		Jacobian before hand. It will be the user's responsibility to form their own preconditioner until we can
		figure out a generic way to precondition the system. */
	int (*lin_precon) (const Matrix<double> &r, Matrix<double> &p, const void *data);

	PJFNK_DATA Newton_data;							///< Data structure for the Newton-Krylov solver (see lark.h)
	const void *activity_data;						///< User defined data structure for an activity model
	const void *residual_data;						///< User defined data structure for the residual function
	const void *precon_data;						///< User defined data structure for preconditioning
	const void *other_data;							///< User define data structure used for user defined residuals
	FILE *OutputFile;								///< Output File pointer

	yaml_cpp_class yaml_object;						///< yaml object to read and access digitized yaml documents (see yaml_wrapper.h)

} SHARK_DATA;

/// Function to print out simulation conditions and options to the output file
void print2file_shark_info(SHARK_DATA *shark_dat);

/// Function to print out the head of species and time stamps to the output file
void print2file_shark_header(SHARK_DATA *shark_dat);

/// Function to print out the simulation results for the current time step
void print2file_shark_results_new(SHARK_DATA *shark_dat);

/// Function to print out the simulation results for the previous time step
void print2file_shark_results_old(SHARK_DATA *shark_dat);

/// Function to calculate the ionic strength of the solution
/** This function calculates the ionic strength of a system given the concentrations of the species present
	in solution, as well as any other relavent information from SHARK_DATA such as charge. 
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param MasterList reference to the MasterSpeciesList object holding species information*/
double calculate_ionic_strength(const Matrix<double> &x, MasterSpeciesList &MasterList);

///Surface Activity function for simple non-ideal adsorption (for adsorption reaction object)
/** This is a simple surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's only parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. NOTE: Only for AdsorptionReaction!

	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the AdsorptionReaction object holding parameter information*/
int FloryHuggins(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for simple non-ideal adsorption (for unsteady adsorption object)
/** This is a simple surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's only parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. NOTE: Only for UnsteadyAdsorption!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the UnsteadyAdsorption object holding parameter information*/
int FloryHuggins_unsteady(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for simple non-ideal adsorption (for multiligand adsorption object)
/** This is a simple surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's only parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. NOTE: Only for MultiligandAdsorption!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the MultiligandAdsorption object holding parameter information*/
int FloryHuggins_multiligand(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for simple non-ideal adsorption (for chemisorption reaction object)
/** This is a simple surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's only parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. NOTE: Only for ChemisorptionReaction!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the AdsorptionReaction object holding parameter information*/
int FloryHuggins_chemi(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for the UNIQUAC model for non-ideal adsorption (for adsorption reaction object)
/** This is a more complex surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's primary parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. However, future development may require some additional parameters, which
	will be accessed later. NOTE: Only for AdsorptionReaction!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the AdsorptionReaction object holding parameter information*/
int UNIQUAC(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for the UNIQUAC model for non-ideal adsorption (for unsteady adsorption object)
/** This is a more complex surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's primary parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. However, future development may require some additional parameters, which
	will be accessed later. NOTE: Only for UnsteadyAdsorption!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the UnsteadyAdsorption object holding parameter information*/
int UNIQUAC_unsteady(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for the UNIQUAC model for non-ideal adsorption (for multiligand adsorption object)
/** This is a more complex surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's primary parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. However, future development may require some additional parameters, which
	will be accessed later. NOTE: Only for MultiligandAdsorption!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the MultiligandAdsorption object holding parameter information*/
int UNIQUAC_multiligand(const Matrix<double> &x, Matrix<double> &F, const void *data);

///Surface Activity function for the UNIQUAC model for non-ideal adsorption (for chemisorption reaction object)
/** This is a more complex surface activity model to be used with the Adsorption objects to evaluate the non-ideal
	behavoir of the surface phase. The model's primary parameters are the shape factors in adsorption and the
	relative concentrations of each surface species. Therefore, we will pass the Adsorption Object itself
	as the const void *data structure. However, future development may require some additional parameters, which
	will be accessed later. NOTE: Only for ChemisorptionReaction!
 
	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to the AdsorptionReaction object holding parameter information*/
int UNIQUAC_chemi(const Matrix<double> &x, Matrix<double> &F, const void *data);

/// Activity function for Ideal Solution
/** This is one of the default activity models available. It assumes the system behaves ideally and sets the
	activity coefficients to 1 for all species.

	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to a data structure needed to evaluate the activity model*/
int ideal_solution (const Matrix<double>& x, Matrix<double> &F, const void *data);

/// Activity function for Davies Equation
/** This is one of the default activity models available. It uses the Davies semi-empirical model to calculate
	average activities of each species in solution. This model is typically valid for systems involving high
	ionic strengths upto 0.5 M (mol/L).

	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to a data structure needed to evaluate the activity model*/
int Davies_equation (const Matrix<double>& x, Matrix<double> &F, const void *data);

/// Activity function for Debye-Huckel Equation
/** This is one of the default activity models available. It uses the Debye-Huckel limiting model to calculate
	average activities of each species in solution. This model is typically valid for systems involving low
	ionic strengths and is only good for solutions between 0 and 0.01 M.

	\param x matrix of the log(C) concentration values at the current non-linear step
	\param F matrix of activity coefficients that are to be altered by this function
	\param data pointer to a data structure needed to evaluate the activity model*/
int DebyeHuckel_equation (const Matrix<double> &x, Matrix<double> &F, const void *data);

/// First test of SIT Model
//int Sit_equation (const Matrix<double>& x, Matrix<double> &F, const void *data);

/// Function takes a given string and returns a flag denoting which surface activity model was choosen
/** This function returns an integer flag that will be one of the valid surface activity model flags from the
	valid_surf_act enum. If the input string was not recognized, then it defaults to returning the IDEAL_ADS flag.
 
	\param input string for the name of the surface activity model*/
int surf_act_choice(const std::string &input);

/// Function takes a given string and returns a flag denoting which activity model was choosen
/** This function returns an integer flag that will be one of the valid activity model flags from the
	valid_act enum. If the input string was not recognized, then it defaults to returning the IDEAL flag.

	\param input string for the name of the activity model*/
int act_choice(const std::string &input);

/// Function takes a give string and returns a flag denoting which type of reactor was choosen for the system
/** This function returns an integer flag that will be one of the valid reactor type flags from the
	valid_mb enum. If the input string was not recognized, then it defaults to returning the BATCH flag.
 
	\param input string for the name of the activity model*/
int reactor_choice(const std::string &input);

/// Function returns a bool to determine the form of line search requested
/** This function returns true if the user requests a bouncing line search algorithm and false if the
	user wants a standard line search. If the input string is unrecognized, then it returns false.

	\param input string for the line search method option*/
bool linesearch_choice(const std::string &input);

/// Function returns the linear solver flag for the PJFNK method
/** This function takes in a string argument and returns the integer flag for the appropriate linear
	solver in PJFNK. If the input string was unrecognized, then it returns the GMRESRP flag.

	\param input string for the linear solver method option*/
int linearsolve_choice(const std::string &input);

/// Function to convert the given values of variables (x) to the log of those variables (logx)
/** This function returns an integer flag to denote success of failure. It takes a constant matrix argument x
	and replaces the elements of the matrix logx with the base 10 log of those x values. This is used mainly
	to convert a set of concentrations (x) to their respective log(C) values (logx).

	\param x matrix of values to take the base 10 log of
	\param logx matrix whose entries are to be changed to base 10 log(x)*/
int Convert2LogConcentration(const Matrix<double> &x, Matrix<double> &logx);

/// Function to convert the given log values of variables (logx) to the values of those variables (x)
/** This function returns an integer flag to denote success of failure. It takes a constant matrix argument logx
	and replaces the elements of the matrix x with 10^logx. This is used mainly to convert a set of log(C) values
	(logx) to their respective concentration values (x).

	\param logx matrix of values to apply as the power of 10 (i.e., 10^logx)
	\param x matrix whose entries are to be changed to the result of 10^logx*/
int Convert2Concentration(const Matrix<double> &logx, Matrix<double> &x);

/// Function to go through the yaml object for the scenario document
/** This function checks the yaml object for the expected keys and values of the scenario document
	to setup the shark simulation for the input given in the input file. */
int read_scenario(SHARK_DATA *shark_dat);

/// Function to go through the yaml object to setup memory space for multiligand objects
/** This function checks the yaml object for the expected keys and values of the multiligand scenario documents
	to setup the shark simulation for the input given in the input file.
	*/
int read_multiligand_scenario(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for the solver options document
/** This function checks the yaml object for the expected keys and values of the solver options document
	to setup the shark simulation for the input given in the input file. */
int read_options(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for the master species document
/** This function checks the yaml object for the expected keys and values of the master species document
	to setup the shark simulation for the input given in the input file. */
int read_species(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for the mass balance document
/** This function checks the yaml object for the expected keys and values of the mass balance document
	to setup the shark simulation for the input given in the input file. */
int read_massbalance(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for the equilibrium reaction document
/** This function checks the yaml object for the expected keys and values of the equilibrium reaction document
	to setup the shark simulation for the input given in the input file. */
int read_equilrxn(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for the unsteady reaction document
/** This function checks the yaml object for the expected keys and values of the unsteady reaction document
	to setup the shark simulation for the input given in the input file. */
int read_unsteadyrxn(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for each Adsorption Object
/** This function checks the yaml object for the expected keys and values of the adsorption object documents
	to setup the shark simulation for the input given in the input file.
	
	\note Each adsorption object will have its own document header by the name of that object*/
int read_adsorbobjects(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for each Unsteady Adsorption Object
/** This function checks the yaml object for the expected keys and values of the unsteady adsorption object documents
	to setup the shark simulation for the input given in the input file.
	
	\note Each unsteady adsorption object will have its own document header by the name of that object*/
int read_unsteadyadsorbobjects(SHARK_DATA *shark_dat);

/// Function to go through the yaml object for each MultiligandAdsorption Object
/** This function checks the yaml object for the expected keys and values of the multiligand object documents
	to setup the shark simulation for the input given in the input file.
	
	\note Each ligand object will have its own document header by the name of that object*/
int read_multiligandobjects(SHARK_DATA *shark_dat);

/// Function to setup the memory and pointers for the SHARK_DATA structure for the current simulation
/** This function will be called after reading the scenario file and is used to setup the memory and other pointers
	for the user requested simulation. This function must be called before running a simulation or trying to read in
	the remander of the yaml formatted input file. Options may be overriden manually after calling this function.

	\param file pointer for the output file where shark results will be stored
	\param residual pointer to the residual function that will be fed into the PJFNK solver
	\param activity pointer to the activity function that will determine the activity coefficients
	\param precond pointer to the linear preconditioning operation to be applied to the Jacobian
	\param dat pointer to the SHARK_DATA data structure
	\param activity_data optional pointer for data needed in activity functions
	\param residual_data optional pointer for data needed in residual functions
	\param precon_data optional pointer for data needed in preconditioning functions
	\param other_data optional pointer for data needed in the evaluation of user defined residual functions*/
int setup_SHARK_DATA( FILE *file, int (*residual) (const Matrix<double> &x, Matrix<double> &res, const void *data),
					  int (*activity) (const Matrix<double> &x, Matrix<double> &gama, const void *data),
					  int (*precond) (const Matrix<double> &r, Matrix<double> &p, const void *data),
					  SHARK_DATA *dat, const void *activity_data, const void *residual_data,
					  const void *precon_data, const void *other_data);

/// Function to add user defined custom residual functions to the OtherList vector object in SHARK_DATA
/** This function will need to be used if the user wants to include custom residuals into the system via the OtherList
	object in SHARK_DATA. For each i residual you want to add, you must call this function passing your residual function and
	the SHARK_DATA structure pointer. The order that those functions are executed in are determined by the integer i.

	\param i index that the other_res function will appear at in the OtherList object
	\param other_res function pointer for the user's custom residual function
	\param shark_dat pointer to the SHARK_DATA data structure*/
int shark_add_customResidual(int i, double (*other_res) (const Matrix<double> &x, SHARK_DATA *shark_dat, const void *other_data),
							 SHARK_DATA *shark_dat);

/// Function to check the Reaction and UnsteadyReaction objects for missing info
/** This function checks the Reaction and UnsteadyReaction objects for missing information. If information is missing, this
	function will return an error that will cause the program to force quit. */
int shark_parameter_check(SHARK_DATA *shark_dat);

/// Function to calculate all Reaction and UnsteadyReaction energies
/** This function will call the calculate energy functions for Reaction and UnsteadyReaction objects.*/
int shark_energy_calculations(SHARK_DATA *shark_dat);

/// Function to calculate all Reaction and UnsteadyReaction parameters as a function of temperature
/** This function will call all temperature dependent functions in Reaction and UnsteadyReaction to calculate equilibirium
	and reaction rate parameters as a function of system temperature. */
int shark_temperature_calculations(SHARK_DATA *shark_dat);

/// Function will search MasterSpeciesList for existance of H + (aq) and OH - (aq) molecules
/** This function searches all molecules in the MasterSpeciesList object for the H + (aq) and OH - (aq) molecules. If they
	are found, then it sets the pH_index and pOH_index of the SHARK_DATA structure and indicates that the system contains
	these variables. */
int shark_pH_finder(SHARK_DATA *shark_dat);

/// Function provides a rough initial guess for the values of all non-linear variables
/** This function constructs an rough initial guess for the values of all non-linear variables in the system. The guess
	is based primarily off of trying to statisfy all mass balance constraints, initial conditions, and pH constraints if
	any apply. */
int shark_guess(SHARK_DATA *shark_dat);

/// Function to establish the initial conditions of the shark simulation
/** This function will establish the initial conditions for a transient problem by solving the speciation of the system
	while holding the transient/unsteady variables constant at their respective initial values. However, if the system
	we are trying to solve is steady, then this function just calls the shark_guess function. */
int shark_initial_conditions(SHARK_DATA *shark_dat);

/// Function to execute a shark simulation at a single time step or pH value
/** This function calls the preprocess, solver, and postprocess functions in order. If a particular solve did not converge,
	then it will retry the solver routine until it runs out of tries or attains convergence. */
int shark_executioner(SHARK_DATA *shark_dat);

/// Function to set up all time steps in the simulation to a specified constant
/** This function will set all time steps for the current simulation to a constant that is specified in the input file.
	The time step will not be changed unless the simulation fails, then it will be reduced in order to try to get the
	system to converge. */
int shark_timestep_const(SHARK_DATA *shark_dat);

/// Function to set up all time steps in the simulation based on success or failure to converge
/** This function will set all time steps for the current simulation based on some factor multiple of the prior time
	step used and whether or not the previous solution step was successful. If the previous step converged, then the
	new time step will be 1.5x the old time step. If it failed, then the simulation will be retried with a new time
	step of 0.5x the old time step. */
int shark_timestep_adapt(SHARK_DATA *shark_dat);

/// Function to call other functions for calculation of parameters and setting of time steps
/** This function will call the shark_temperature_calculations function and the appropriate time step function. If the
	user requests a constant time step, it will call the shark_timestep_const function. Otherwise, it calls the
	shark_timestep_adapt function. */
int shark_preprocesses(SHARK_DATA *shark_dat);

/// Function to call the PJFNK solver routine given the current SHARK_DATA information
/** This function will perform the necessary steps before and after calling the PJFNK solver routine. Based on the
	simulation flags, the solver function will perform an intial guess for unsteady variables, call the PJFNK method,
	and the printout a console message about the performance. If a terminal failure occurs during the solver, it will
	print out the current state of residuals, variables, and the Jacobian matrix to the console. Analyzing this information
	could provide clues as to why failure occured. */
int shark_solver(SHARK_DATA *shark_dat);

/// Function to convert PJFNK solutions to concentration values and print to the output file
/** This function will convert the non-linear variables to their respective concentration values, then print the solve
	information out to the output file. */
int shark_postprocesses(SHARK_DATA *shark_dat);

/// Function to reset the values of all stateful information in SHARK_DATA
/** This function will reset all stateful matrix data in the SHARK_DATA structure in preparation of the next time
	step simulation. */
int shark_reset(SHARK_DATA *shark_dat);

/// Default residual function for shark evaluations
/** This function calls each individual object's residual function to formulate the overall residual function used
	in the PJFNK solver routine. It will also call the activity function. The order in which these function calls
	occurs is as follows: (i) activities, (ii) Reaction, (iii) UnsteadyReaction, (iv) MassBalance, (v) OtherList,
	and (vi) MasterSpeciesList. If a constant pH is specified, then the MasterSpeciesList residual call is replaced
	with a constraint on the H + (aq) variable (if one exists). */
int shark_residual(const Matrix<double> &x, Matrix<double> &F, const void *data);

/// Function to call all above functions to perform a shark simulation
/** This function is called after reading in all inputs, setting all constants, and calling the setup function. It
	will call all the necessary functions and subroutines iteratively until the desired simulation is complete. */
int SHARK(SHARK_DATA *shark_dat);

/// Function to perform a shark simulation based on the conditions in a yaml formatted input file
/** This is the primary function used to run shark simulations from the UI. It requires that ths user provide one
	input file that is formatted with yaml keys, symbols, and spacing so that it can be recognized by the parser.
	This style of input file is much easier to use and understand than the input files used for SCOPSOWL or SKUA.
	Below shows an example of a typical input file. Note that the # symbol is used in the input file to comment
	out lines of text that the parser does not need to read. \n

	Example Yaml Input for SHARK
	----------------------------
	\#This will serve as a test input file for shark to demo how to structure the document \n
	\#In practice, this section should be listed first, but it doesn't really matter \n
	\#DO NOT USE TABS IN THESE INPUT FILES \n
	\#--- Starts a document ... Ends a document \n
	\#All keys must be proceeded by a : \n
	\#All lists/header must be preceeded by a - \n
	\#Spacing of the keys will indicate which list/header they belong to \n
	Scenario: \n
	--- \n
	- vars_fun: \n
	  numvar: 25 \n
	  num_ssr: 15 \n
	  num_mbe: 7 \n
	  num_usr: 2 \n
      num_other: 0      \#Not required or used in current version \n

	- sys_data: \n
	  act_fun: davies \n
	  const_pH: false \n
	  pH: 7              \#Only required if we are specifying a const_pH \n
	  temp: 298.15       \#Units must be in Kelvin \n
	  dielec: 78.325     \#Units must be in (1/Kelvin) \n
	  rel_perm: 80.1     \#Unitless number \n
	  res_alk: 0         \#Units must be in mol/L (Residual Alkalinity) \n
	  volume: 1.0		 \#Units must be in L \n

	- run_time: \n
      steady: false         \#NOTE: All time must be represented in hours \n
	  specs_curve: false    \#Only needed if steady = true, and will default to false \n
	  dt: 0.001             \#Only required if steady = false \n
	  time_adapt: true      \#Only needed if steady = false, and will default to false \n
	  sim_time: 96.0        \#Only required if steady = false \n
	  t_out: 0.01           \#Only required if steady = false \n
	...

	\#The following header is entirely optional, but is used to set solver options \n
	SolverOptions: \n
	--- \n
	line_search: true      \#Default = true, and is recommended to be true \n
	search_type: standard \n
	linear_solve: gmresrp      \#Note: FOM will be fastest for small problems \n
	restart: 25                \#Note: restart only used if using GMRES or GCR type solvers \n
	nl_maxit: 50 \n
	nl_abstol: 1e-5 \n
	nl_reltol: 1e-8 \n
	lin_reltol: 1e-10 	    \#Min Tol = 1e-15 \n
	lin_abstol: 1e-10		\#Min Tol = 1e-15 \n
	nl_print: true \n
	l_print: true \n
	...

	\#After the Scenario read, shark will call the setup_function, then read info below \n
	MasterSpecies: \n
	--- \n
	\#Header names are specific \n
	\#Keys are chosen by user, but must span numbers 0 through numvar-1 \n
	\#Keys will denote the ordering of the variables \n
	\#Note: Currently, the number of reg molecules is very limited \n
	- reg: \n
	  0: Cl - (aq) \n
	  1: NaHCO3 (aq) \n
	  2: NaCO3 - (aq) \n
	  3: Na + (aq) \n
	  4: HNO3 (aq) \n
	  5: NO3 - (aq) \n
	  6: H2CO3 (aq) \n
      7: HCO3 - (aq) \n
	  8: CO3 2- (aq) \n
	  9: UO2 2+ (aq) \n
	  10: UO2NO3 + (aq) \n
	  11: UO2(NO3)2 (aq) \n
	  12: UO2OH + (aq) \n
	  13: UO2(OH)3 - (aq) \n
	  14: (UO2)2(OH)2 2+ (aq) \n
	  15: (UO2)3(OH)5 + (aq) \n
	  16: UO2CO3 (aq) \n
	  17: UO2(CO3)2 2- (aq) \n
	  18: UO2(CO3)3 4- (aq) \n
	  19: H2O (l) \n
	  20: OH - (aq) \n
	  21: H + (aq) \n

	\#Keys for the sub-headers must follow same rules as keys from above \n
	- unreg: \n
	  - 22: \n
	    formula: A(OH)2 (aq) \n
	    charge: 0 \n
	    enthalpy: 0 \n
	    entropy: 0 \n
	    have_HS: false \n
	    energy: 0 \n
	    have_G: false \n
	    phase: Aqueous \n
	    name: Amidoxime \n
	   lin_form: none \n

	  - 23: \n
	    formula: UO2AO2 (aq) \n
	    charge: 0 \n
	    enthalpy: 0 \n
	    entropy: 0 \n
	    have_HS: false \n
	    energy: 0 \n
	    have_G: false \n
	    phase: Aqueous \n
	    name: Uranyl-amidoximate \n
	    lin_form: none \n

	  - 24: \n
	    formula: UO2CO3AO2 2- (aq) \n
	    charge: -2 \n
	    enthalpy: 0 \n
	    entropy: 0 \n
	    have_HS: false \n
	    energy: 0 \n
	    have_G: false \n
	    phase: Aqueous \n
	    name: Uranyl-carbonate-amidoximate \n
	    lin_form: none \n
	... \n

	\#NOTE: Total concentrations must be given in mol/L \n
	MassBalance: \n
	--- \n
	\#Header names under MassBalance are choosen by the user \n
	\#All other keys will be checked \n
	- water: \n
	  total_conc: 1 \n
	  - delta: \n
	    "H2O (l)": 1 \n

	- carbonate: \n
	  total_conc: 0.0004175 \n
	  - delta: \n
	   "NaHCO3 (aq)": 1 \n
	   "NaCO3 - (aq)": 1 \n
	   "H2CO3 (aq)": 1 \n
	   "HCO3 - (aq)": 1 \n
	   "CO3 2- (aq)": 1 \n
	   "UO2CO3 (aq)": 1 \n
	   "UO2(CO3)2 2- (aq)": 2 \n
	   "UO2(CO3)3 4- (aq)": 3 \n
	   "UO2CO3AO2 2- (aq)": 1 \n

	\#Other mass balances skipped for demo purposes... \n
	... \n

	\#Document for equilibrium or steady reactions \n
	EquilRxn: \n
	--- \n
	\#Headers under EquilRxn separate out each reaction object \n
	\#Keys for these headers only factor into the order of the equations \n
	\#Stoichiometry follows the convention that products are pos(+) and reactants are neg(-) \n
	\#Note: logK is only required if any species in stoichiometry is unregistered \n
	\#Example: below represents - {H2O (l)} --> {H + (aq)} + {OH - (aq)} \n
	\#Note: a valid reaction statement requires at least 1 stoichiometry args \n
	\#Note: You can also provide reaction energies: enthalpy, entropy, and energy \n

	- rxn00: \n
	  logK: -14 \n
	  - stoichiometry: \n
	   "H2O (l)": -1 \n
	   "OH - (aq)": 1 \n
	   "H + (aq)": 1 \n

	- rxn01: \n
	  logK: -6.35 \n
	   - stoichiometry: \n
	     "H2CO3 (aq)": -1 \n
	     "HCO3 - (aq)": 1 \n
	     "H + (aq)": 1 \n

	\#Other reactions skipped for demo purposes... \n
	... \n

	\#Document for unsteady reactions \n
	UnsteadyRxn: \n
	--- \n
	\#Same basic standards for this doc as the EquilRxn \n
	\#Main difference is the inclusion of rate information \n
	\#You are required to give at least 1 rate \n
	\#You are also required to denote which variable is unsteady \n
	\#You must give the initial concentration for the variable in mol/L \n
	\#Rate units are in (L/mol)^n/hr \n
	\#Note: we also have keys for forward_ref, reverse_ref, \n
	\#activation_energy, and temp_affinity. \n
	\#These are optional if forward and/or reverse are given \n
	\#Note: You can also provide reaction energies: enthalpy, entropy, and energy \n

	- rxn00: \n
	  unsteady_var: UO2AO2 (aq) \n
	  initial_condition: 0 \n
	  logK: -1.35 \n
	  forward: 4.5e+6 \n
	  reverse: 1.00742e+8 \n
	  - stoichiometry: \n
	    "UO2 2+ (aq)": -1 \n
	    "A(OH)2 (aq)": -1 \n
	    "UO2AO2 (aq)": 1 \n
	    "H + (aq)": 2 \n

	- rxn01: \n
	  unsteady_var: UO2CO3AO2 2- (aq) \n
	  initial_condition: 0 \n
	  logK: 3.45 \n
	  forward: 2.55e+15 \n
	  reverse: 9.04774e+11 \n
	  - stoichiometry: \n
	    "UO2 2+ (aq)": -1 \n
	    "CO3 2- (aq)": -1 \n
	    "A(OH)2 (aq)": -1 \n
	    "UO2CO3AO2 2- (aq)": 1 \n
	    "H + (aq)": 2 \n

	... \n

	\note It may be advantageous to look at some other shark input file examples. More input files are provided in the
	input_files/SHARK directory of the ecosystem project folder. Please refer to your own source file location for more
	input file examples for SHARK.
 */
int SHARK_SCENARIO(const char *yaml_input);

/// Function to perform a series of shark calculation tests
/** This function sets up and solves a test problem for shark. It is callable from the UI. */
int SHARK_TESTS();

/// Function to perform a series of shark calculation tests (older version)
/** This function sets up and solves a test problem for shark. It is NOT callable from the UI. */
int SHARK_TESTS_OLD();

#endif
