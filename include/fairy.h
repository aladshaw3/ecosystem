/*!
 *  \file fairy.h fairy.cpp
 *	\brief Fission-products from Atomic Incident and their Respective Yields
 *	\details This file contains a C++ object for determining fission products and
 *			their yields from some nuclear event based on: (i) type of fission,
 *			either neutron-induced or spontaneous, (ii) energy level of neutron
 *			source or bomb yield, (iii) extent of fission, and (iv) initial mass
 *			and composition of fuel or bomb materials. Data for fission products
 *			comes from ENDF-6 data libraries that were read with a python script
 *			and output into a yaml format (see 'scripts/fission-product-yields').
 *
 *  \author Austin Ladshaw
 *	\date 12/07/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for research in the area
 *             of radioactive particle decay and transport. Copyright (c) 2018, 
 *				all rights reserved.
 */

#ifndef FAIRY_HPP_
#define FAIRY_HPP_

#include "ibis.h"

/// Enumeration for Fission Product Yield Type
typedef enum {neutron, spontaneous, explosion} fiss_type;

/// Function to determine the fission type based on given string
fiss_type fisstype_choice(std::string &choice);

/// FissionProducts class object to create decay chains from fission yields
/** This object inherits from DecayChain and will formulate decay chains based on starting weapon or
	fuel materials, as well as the type of fission taking place and the extent of fission for a given
	mass of the starting material. Fission Product Yields (FPYs) are based on the Evaluated Nuclear 
	Data Format (ENDF) libraries, which contain FPYs for a number of fuel and weapon materials at
	various energy levels for two different types of fission: (i) neutron-induced and (ii) spontaneous.
	*/
class FissionProducts : DecayChain
{
public:
	FissionProducts();								///< Default Constructor
	~FissionProducts();								///< Default Destructor
	
	void DisplayInfo();								///< Display the FAIRY information, initial materials and fractions
	
	void loadNuclides(yaml_cpp_class &data);		///< Function to load the nuclide library into the pointer
	void unloadNuclides();							///< Delete the pointer to nuclide library to free space
	
	void loadFissionProductYields();				///< Function to read in the library and put into the Yaml object
	void unloadFissionProductYields();				///< Function to delete the items in the Yaml object
	
	void setFissionType(fiss_type opt);				///< Set the fission type for the Fission Products
	void setTotalMass(double mass);					///< Set the total mass of fissible materials (kg)
	void setFissionExtent(double per);				///< Set the extent of fission parameter (%)
	void setEnergyLevel(double el);					///< Set the energy level for neutron source (eV)
	void setThreshold(double val);					///< Set the threshold value for half-life (in sec)
	
	/// Automatically sets the threshold value given time (in sec) after detonation and percent conversion rate (default 99%)
	void timeSinceDetonation(double time, double per);
	
	void setTimeUnits(time_units units);			///< Set the time units for the simulation
	void setEndTime(double end);					///< Set the end time for the simulation (in choosen units)
	void setTimeSteps(int steps);					///< Set the number of time steps to simulate
	void setVerifySoln(bool opt);					///< Set boolean to verify eigen solution
	void setPrintSparsity(bool opt);				///< Set boolean to print sparsity on output
	void setPrintChain(bool opt);					///< Set boolean to print chain info
	void setPrintResults(bool opt);					///< Set boolean to print results of simulation
	void setConsoleOut(bool opt);					///< Set boolean to print output to console
	
	int run_simulation(std::string file_name);		///< Run a decay simulation from Fission Products
	
	int print_yields();								///< Print yield data for weapon or fuel to output file 
	
	/// NOTE: For each of the below functions, you must first call the loadNuclides function
	void addIsotopeMaterial(std::string iso, double percent);				///< Add an isotope for the fissible material (checks string value)
	void addIsotopeMaterial(std::string sym, int iso_num, double percent);	///< Add an isotope for the fissible material (by symbol and mass num)
	void addIsotopeMaterial(int atom_num, int iso_num, double percent);		///< Add an isotope for the fissible material (by atom num and mass num)
	
	void checkFractions();							///< Check fractions of materials in fuel/weapone and correct if necessary
	int evaluateYields();							///< Read yield data and set isotope fractionation based on yields
	
protected:
	fiss_type type;								///< Type of fission products to be produced
	std::vector<Isotope> InitialMat;			///< Initial materials/isotopes to undergoe fission (Conc. in moles)
	std::vector<double> MatFrac;				///< Material fractionation of the initial material (%)
	std::map<int, double> Yields;				///< Map of fission yields of material based on atomic mass
	double total_mass;							///< Total mass of the weapon or fuel rod (kg)
	double fiss_extent;							///< Percentage of the starting material that undergoes fission (%)
	double energy_level;						///< Energy level of neutron source (eV)
	yaml_cpp_class fpy_data;					///< Yaml object to read and store the FPY library files
	
private:
	
};

/// Test function for FAIRY
int FAIRY_TESTS();

#endif /* FAIRY_HPP_ */
