/*!
 *  \file crane.h
 *	\brief Cloud Rise After Nuclear Explosion
 *  \author Austin Ladshaw
 *	\date 05/30/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of radioactive particle aggregation and transport. Copyright (c) 2018,
 *             all rights reserved.
 */

#include "crane.h"

/*
 *								Start: Crane Class Definitions
 *	-------------------------------------------------------------------------------------
 */

//Default constructor
Crane::Crane()
{
	eps = 18.0/29.0;
	grav = 9.8;
	k = 0.0;
	k2 = 0.075;
	k3 = 0.175;
	k6 = 1.0;
	mu = 0.12;
	k_temp = 1.0;
	apparent_temp = 298.0;
	apparent_amb_temp = 298.0;
	q_x = 1.0;
	q_xe = 1.0;
	beta_prime = 1.0;
	xe = 0.0;
	char_vel = 0.0;
	vert_rad = 0.0;
	horz_rad = 0.0;
	virtual_mass = 0.0;
	gas_const = 287.0;
	latent_heat = 2.5e6;
	sigma_turbulence = 0.0;
	mean_spec_heat = 1005.0;
	actual_spec_heat = 1005.0;
	spec_heat_water = 2108.0;
	spec_heat_conds = 1480.0;
	spec_heat_entrain = 1005.0;
	spec_heat_entrain_integral = 0.0;
	cloud_volume = 0.0;
	equil_temp = 1000.0;
	total_mass_fallout_rate = 0.0;
	surf_area = 0.0;
	shear_ratio = 0.0;
	shear_vel = 0.0;
	part_density = 2600.0;
	adjusted_height = 0.0;
	bomb_yield = 0.0;
	force_factor = 0.0;
	det_alt = 0.0;
	burst_height = 0.0;
	ground_alt = 0.0;
	energy_frac = 0.4;
	eccentricity = 0.75;
	air_density = 1.225;
	air_viscosity = 1.81e-5;
	slip_factor = 1.0;
	davies_num = 1.0;
	vapor_pressure = 0.0;
	sat_vapor_pressure = 1.0;
	initial_soil_mass = 0.0;
	initial_water_mass = 0.0;
	initial_air_mass = 0.0;
	current_amb_temp = 298.0;
	current_atm_press = 101325.0;
	includeShearVel = false;
	isSaturated = false;
	ConsoleOut = false;
	FileOut = false;
	isTight = true;
	saturation_time = 0.0;
	min_dia = 0.0001;
	max_dia = 10000.0;
	mean_dia = 0.407;
	std_dia = 4.0;
	num_bins = 22;
	cloud_mass = 0.0;
	entrained_mass = 0.0;
	cloud_rise = 0.0;
	cloud_alt = 0.0;
	x_water_vapor = 0.0;
	w_water_conds = 0.0;
	s_soil = 0.0;
	temperature = 298;
	energy = 0.0;
	current_time = 0.0;
	
	create_default_atmosphere();
}

//Default destructor
Crane::~Crane()
{
	delete_atmosphere();
	delete_particles();
}

// Below are some display functions used for testing different functions
void Crane::display_part_hist()
{
	std::cout << "Normalized Particle Distribution by Size\n";
	std::cout << "----------------------------------------\n";
	std::cout << "Size (um)\tNormalDist\n";
	//Iterate through map
	for (std::map<double,double>::iterator it=this->part_hist.begin(); it!=this->part_hist.end(); ++it)
	{
		std::cout << it->first << "\t" << it->second << std::endl;
	}
}

void Crane::display_part_conc()
{
	std::cout << "Particle Concentration Distribution by Size\n";
	std::cout << "----------------------------------------\n";
	std::cout << "Size (um)\tConcDist\n";
	//Iterate through map
	for (std::map<double,double>::iterator it=this->part_conc.begin(); it!=this->part_conc.end(); ++it)
	{
		std::cout << it->first << "\t" << it->second << std::endl;
	}
}

///< Below are the set functions for various parameters
void Crane::set_eps(double val)
{
	this->eps = val;
}

void Crane::set_grav(double val)
{
	this->grav = val;
}

void Crane::set_k(double val)
{
	this->k = val;
}

void Crane::set_k2(double val)
{
	this->k2 = val;
}

void Crane::set_k3(double val)
{
	this->k3 = val;
}

void Crane::set_k6(double val)
{
	this->k6 = val;
}

void Crane::set_mu(double val)
{
	this->mu = val;
}

void Crane::set_k_temp(double val)
{
	this->k_temp = val;
}

void Crane::set_apparent_temp(double val)
{
	this->apparent_temp = val;
}

void Crane::set_apparent_amb_temp(double val)
{
	this->apparent_amb_temp = val;
}

void Crane::set_q_x(double val)
{
	this->q_x = val;
}

void Crane::set_q_xe(double val)
{
	this->q_xe = val;
}

void Crane::set_beta_prime(double val)
{
	this->beta_prime = val;
}

void Crane::set_xe(double val)
{
	this->xe = val;
}

void Crane::set_char_vel(double val)
{
	this->char_vel = val;
}

void Crane::set_vert_rad(double val)
{
	this->vert_rad = val;
}

void Crane::set_horz_rad(double val)
{
	this->horz_rad = val;
}

void Crane::set_virtual_mass(double val)
{
	this->virtual_mass = val;
}

void Crane::set_gas_const(double val)
{
	this->gas_const = val;
}

void Crane::set_latent_heat(double val)
{
	this->latent_heat = val;
}

void Crane::set_sigma_turbulence(double val)
{
	this->sigma_turbulence = val;
}

void Crane::set_mean_spec_heat(double val)
{
	this->mean_spec_heat = val;
}

void Crane::set_actual_spec_heat(double val)
{
	this->actual_spec_heat = val;
}

void Crane::set_spec_heat_water(double val)
{
	this->spec_heat_water = val;
}

void Crane::set_spec_heat_entrain(double val)
{
	this->spec_heat_entrain = val;
}

void Crane::set_spec_heat_entrain_integral(double val)
{
	this->spec_heat_entrain_integral = val;
}

void Crane::set_spec_heat_conds(double val)
{
	this->spec_heat_conds = val;
}

void Crane::set_cloud_volume(double val)
{
	this->cloud_volume = val;
}

void Crane::set_equil_temp(double val)
{
	this->equil_temp = val;
}

void Crane::set_total_mass_fallout_rate(double val)
{
	this->total_mass_fallout_rate = val;
}

void Crane::set_surf_area(double val)
{
	this->surf_area = val;
}

void Crane::set_shear_ratio(double val)
{
	this->shear_ratio = val;
}

void Crane::set_shear_vel(double val)
{
	this->shear_vel = val;
}

void Crane::set_part_density(double val)
{
	this->part_density = val;
}

void Crane::set_adjusted_height(double val)
{
	this->adjusted_height = val;
}

void Crane::set_bomb_yield(double val)
{
	this->bomb_yield = val;
}

void Crane::set_force_factor(double val)
{
	this->force_factor = val;
}

void Crane::set_det_alt(double val)
{
	this->det_alt = val;
}

void Crane::set_burst_height(double val)
{
	this->burst_height = val;
}

void Crane::set_ground_alt(double val)
{
	this->ground_alt = val;
}

void Crane::set_energy_frac(double val)
{
	this->energy_frac = val;
}

void Crane::set_eccentricity(double val)
{
	this->eccentricity = val;
}

void Crane::set_air_density(double val)
{
	this->air_density = val;
}

void Crane::set_air_viscosity(double val)
{
	this->air_viscosity = val;
}

void Crane::set_slip_factor(double val)
{
	this->slip_factor = val;
}

void Crane::set_davies_num(double val)
{
	this->davies_num = val;
}

void Crane::set_min_dia(double val)
{
	this->min_dia = val;
}

void Crane::set_max_dia(double val)
{
	this->max_dia = val;
}

void Crane::set_mean_dia(double val)
{
	this->mean_dia = val;
}

void Crane::set_std_dia(double val)
{
	this->std_dia = val;
}

void Crane::set_num_bins(int val)
{
	this->num_bins = val;
}

void Crane::set_cloud_mass(double val)
{
	this->cloud_mass = val;
}

void Crane::set_entrained_mass(double val)
{
	this->entrained_mass = val;
}

void Crane::set_cloud_rise(double val)
{
	this->cloud_rise = val;
}

void Crane::set_cloud_alt(double val)
{
	this->cloud_alt = val;
}

void Crane::set_x_water_vapor(double val)
{
	this->x_water_vapor = val;
}

void Crane::set_w_water_conds(double val)
{
	this->w_water_conds = val;
}

void Crane::set_s_soil(double val)
{
	this->s_soil = val;
}

void Crane::set_temperature(double val)
{
	this->temperature = val;
}

void Crane::set_energy(double val)
{
	this->energy = val;
}

void Crane::set_current_time(double val)
{
	this->current_time = val;
}

void Crane::set_vapor_pressure(double val)
{
	this->vapor_pressure = val;
}

void Crane::set_sat_vapor_pressure(double val)
{
	this->sat_vapor_pressure = val;
}

void Crane::set_initial_soil_mass(double val)
{
	this->initial_soil_mass = val;
}

void Crane::set_initial_water_mass(double val)
{
	this->initial_water_mass = val;
}

void Crane::set_initial_air_mass(double val)
{
	this->initial_air_mass = val;
}

void Crane::set_current_amb_temp(double val)
{
	this->current_amb_temp = val;
}

void Crane::set_current_atm_press(double val)
{
	this->current_atm_press = val;
}

void Crane::set_includeShearVel(bool val)
{
	this->includeShearVel = val;
}

void Crane::set_isSaturated(bool val)
{
	this->isSaturated = val;
}

void Crane::set_ConsoleOut(bool val)
{
	this->ConsoleOut = val;
}

void Crane::set_FileOut(bool val)
{
	this->FileOut = val;
}

void Crane::set_saturation_time(double val)
{
	this->saturation_time = val;
}

void Crane::set_isTight(bool val)
{
	this->isTight = val;
}

///< Below are the get functions for various parameters

double Crane::get_eps()
{
	return this->eps;
}

double Crane::get_grav()
{
	return this->grav;
}

double Crane::get_k()
{
	return this->k;
}

double Crane::get_k2()
{
	return this->k2;
}

double Crane::get_k3()
{
	return this->k3;
}

double Crane::get_k6()
{
	return this->k6;
}

double Crane::get_mu()
{
	return this->mu;
}

double Crane::get_k_temp()
{
	return this->k_temp;
}

double Crane::get_apparent_temp()
{
	return this->apparent_temp;
}

double Crane::get_apparent_amb_temp()
{
	return this->apparent_amb_temp;
}

double Crane::get_q_x()
{
	return this->q_x;
}

double Crane::get_q_xe()
{
	return this->q_xe;
}

double Crane::get_beta_prime()
{
	return this->beta_prime;
}

double Crane::get_xe()
{
	return this->xe;
}

double Crane::get_char_vel()
{
	return this->char_vel;
}

double Crane::get_vert_rad()
{
	return this->vert_rad;
}

double Crane::get_horz_rad()
{
	return this->horz_rad;
}

double Crane::get_virtual_mass()
{
	return this->virtual_mass;
}

double Crane::get_gas_const()
{
	return this->gas_const;
}

double Crane::get_latent_heat()
{
	return this->latent_heat;
}

double Crane::get_sigma_turbulence()
{
	return this->sigma_turbulence;
}

double Crane::get_mean_spec_heat()
{
	return this->mean_spec_heat;
}

double Crane::get_actual_spec_heat()
{
	return this->actual_spec_heat;
}

double Crane::get_spec_heat_water()
{
	return this->spec_heat_water;
}

double Crane::get_spec_heat_entrain()
{
	return this->spec_heat_entrain;
}

double Crane::get_spec_heat_entrain_integral()
{
	return this->spec_heat_entrain_integral;
}

double Crane::get_spec_heat_conds()
{
	return this->spec_heat_conds;
}

double Crane::get_cloud_volume()
{
	return this->cloud_volume;
}

double Crane::get_equil_temp()
{
	return this->equil_temp;
}

double Crane::get_total_mass_fallout_rate()
{
	return this->total_mass_fallout_rate;
}

double Crane::get_surf_area()
{
	return this->surf_area;
}

double Crane::get_shear_ratio()
{
	return this->shear_ratio;
}

double Crane::get_shear_vel()
{
	return this->shear_vel;
}

double Crane::get_part_density()
{
	return this->part_density;
}

double Crane::get_adjusted_height()
{
	return this->adjusted_height;
}

double Crane::get_bomb_yield()
{
	return this->bomb_yield;
}

double Crane::get_force_factor()
{
	return this->force_factor;
}

double Crane::get_det_alt()
{
	return this->det_alt;
}

double Crane::get_burst_height()
{
	return this->burst_height;
}

double Crane::get_ground_alt()
{
	return this->ground_alt;
}

double Crane::get_energy_frac()
{
	return this->energy_frac;
}

double Crane::get_eccentricity()
{
	return this->eccentricity;
}

double Crane::get_air_density()
{
	return this->air_density;
}

double Crane::get_air_viscosity()
{
	return this->air_viscosity;
}

double Crane::get_slip_factor()
{
	return this->slip_factor;
}

double Crane::get_davies_num()
{
	return this->davies_num;
}

double Crane::get_min_dia()
{
	return this->min_dia;
}

double Crane::get_max_dia()
{
	return this->max_dia;
}

double Crane::get_mean_dia()
{
	return this->mean_dia;
}

double Crane::get_std_dia()
{
	return this->std_dia;
}

int Crane::get_num_bins()
{
	return this->num_bins;
}

double Crane::get_cloud_mass()
{
	return this->cloud_mass;
}

double Crane::get_entrained_mass()
{
	return this->entrained_mass;
}

double Crane::get_cloud_rise()
{
	return this->cloud_rise;
}

double Crane::get_cloud_alt()
{
	return this->cloud_alt;
}

double Crane::get_x_water_vapor()
{
	return this->x_water_vapor;
}

double Crane::get_w_water_conds()
{
	return this->w_water_conds;
}

double Crane::get_s_soil()
{
	return this->s_soil;
}

double Crane::get_temperature()
{
	return this->temperature;
}

double Crane::get_energy()
{
	return this->energy;
}

double Crane::get_current_time()
{
	return this->current_time;
}

double Crane::get_vapor_pressure()
{
	return this->vapor_pressure;
}

double Crane::get_sat_vapor_pressure()
{
	return this->sat_vapor_pressure;
}

double Crane::get_initial_soil_mass()
{
	return this->initial_soil_mass;
}

double Crane::get_initial_water_mass()
{
	return this->initial_water_mass;
}

double Crane::get_initial_air_mass()
{
	return this->initial_air_mass;
}

double Crane::get_current_amb_temp()
{
	return this->current_amb_temp;
}

double Crane::get_current_atm_press()
{
	return this->current_atm_press;
}

bool Crane::get_includeShearVel()
{
	return this->includeShearVel;
}

bool Crane::get_isSaturated()
{
	return this->isSaturated;
}

double Crane::get_part_size(int i)
{
	if (i < 0 || i >= this->part_size.size())
		return 0.0;
	
	return this->part_size[i];
}

double Crane::get_settling_rate(double Dj)
{
	return this->settling_rate[Dj];
}

bool Crane::get_ConsoleOut()
{
	return this->ConsoleOut;
}

bool Crane::get_FileOut()
{
	return this->FileOut;
}

Matrix<double> & Crane::get_part_conc_var()
{
	return this->part_conc_var;
}

double Crane::get_saturation_time()
{
	return this->saturation_time;
}

bool Crane::get_isTight()
{
	return this->isTight;
}

// Below are listed all the compute function for various parameters
void Crane::compute_beta_prime(double x, double s, double w)
{
	double val = (1.0 + x) / (1.0 + x + s + w);
	this->set_beta_prime(val);
}

void Crane::compute_q_x(double x)
{
	double val = (1.0 + (this->get_eps()/x)) / (1.0 + x);
	this->set_q_x(val);
}

void Crane::compute_q_xe(double xe)
{
	double val = (1.0 + (this->get_eps()/xe)) / (1.0 + xe);
	this->set_q_xe(val);
}

void Crane::compute_apparent_temp(double T, double x)
{
	this->compute_q_x(x);
	this->set_apparent_temp(T*this->get_q_x());
}

void Crane::compute_apparent_amb_temp(double Te, double xe)
{
	this->compute_q_xe(xe);
	this->set_apparent_amb_temp(Te*this->get_q_xe());
}

void Crane::compute_char_vel(double u, double E)
{
	this->set_char_vel(fmax(fabs(u), 2.0*E));
}

void Crane::compute_air_viscosity(double T)
{
	this->set_air_viscosity( (145.8*pow(10.0, -8.0)*pow(T, (3.0/2.0))) / (110.4 + T) );
}

void Crane::compute_vapor_pressure(double P, double x)
{
	this->set_vapor_pressure( (P*x) / (this->get_eps() + x) );
}

void Crane::compute_sat_vapor_pressure(double T)
{
	//this->set_sat_vapor_pressure( 611.0*pow(T/273.0, -5.13)*exp(25.0*(T-273.0)/T) );
	
	//NOTE: Changed DEFLIC's model because it failed to produce good results at high temperatures
	
	T = T - 273.15;
	double Pws = 618.8*exp(17.27*T/(T+237.3));
	this->set_sat_vapor_pressure(Pws);
}

void Crane::compute_xe(double Te, double P, double HR)
{
	double val = (109.98*HR/29.0/P)*pow(Te/273.0, -5.13)*exp(25.0*(Te-273.0)/Te);
	this->set_xe(val);
}

void Crane::compute_air_density(double P, double x, double T)
{
	//Assume P and Pws come in as Pa --> convert to mBar
	//this->compute_sat_vapor_pressure(T);
	//P = P*0.01;
	//double Pws = this->get_sat_vapor_pressure()*0.01;
	//double val = ( P - (Pws*HR*(1.0-this->eps)/100.0) ) / (2.8679*T);
	//this->set_air_density(val);
	
	//NOTE: Changed DEFLIC's model because it failed to produce good results at high humidity
	
	this->set_air_density( (P/this->get_gas_const()/T)*(1.0+x)/(1.0+x*1.609) );
}

void Crane::compute_spec_heat_entrain(double T)
{
	if (T <= 2300.0)
	{
		this->set_spec_heat_entrain(946.6+(0.1971*T));
	}
	else
	{
		this->set_spec_heat_entrain(-3587.5+(2.125*T));
	}
}

void Crane::compute_spec_heat_water(double T)
{
	this->set_spec_heat_water(1697.66+(1.144174*T));
}

void Crane::compute_spec_heat_conds(double T)
{
	if (T <= 848.0)
	{
		this->set_spec_heat_conds(781.6+(0.5612*T)-(1.881e7/T/T));
	}
	else
	{
		this->set_spec_heat_conds(1003.8+(0.1351*T));
	}
}

void Crane::compute_actual_spec_heat(double T, double x)
{
	this->compute_spec_heat_entrain(T);
	this->compute_spec_heat_water(T);
	this->set_actual_spec_heat( (this->get_spec_heat_entrain() + (x*this->get_spec_heat_water())) / (1.0+x) );
}

void Crane::compute_k_temp(double T)
{
	if (T > this->equil_temp)
	{
		this->set_k_temp(0.0);
	}
	else
	{
		this->set_k_temp(1.0);
	}
}

void Crane::compute_mean_spec_heat(double T, double x, double s, double w)
{
	this->compute_k_temp(T);
	this->compute_beta_prime(x, s, w);
	this->compute_actual_spec_heat(T, x);
	this->compute_spec_heat_conds(T);
	this->set_mean_spec_heat( this->get_beta_prime()*this->get_actual_spec_heat() + (1.0-this->get_beta_prime())*this->get_k_temp()*this->get_spec_heat_conds() );
}

void Crane::compute_cloud_volume(double m, double x, double s, double w, double T, double P)
{
	m = m * 1000.0; //Convert from Mg to kg
	this->compute_beta_prime(x, s, w);
	this->compute_apparent_temp(T, x);
	this->set_cloud_volume( m*this->get_beta_prime()*this->get_gas_const()*this->get_apparent_temp()/P );
}

void Crane::compute_vert_rad(double z)
{
	this->set_vert_rad( this->get_mu()*(z - this->get_adjusted_height()) );
}

void Crane::compute_horz_rad(double m, double x, double s, double w, double T, double P, double z)
{
	this->compute_cloud_volume(m, x, s, w, T, P);
	this->compute_vert_rad(z);
	this->set_horz_rad( sqrt(3.0*this->get_cloud_volume()/4.0/M_PI/this->get_vert_rad()) );
}

void Crane::compute_sigma_turbulence(double E, double z)
{
	this->compute_vert_rad(z);
	this->set_sigma_turbulence( this->get_k3()*pow(2.0*E, 3.0/2.0)/this->get_vert_rad() );
}

void Crane::compute_surf_area(double m, double x, double s, double w, double T, double P, double z)
{
	this->compute_horz_rad(m, x, s, w, T, P, z);
	this->set_surf_area( 4.0*M_PI*pow(this->get_horz_rad(), 2.0) );
}

void Crane::compute_shear_vel(double z, Matrix<double> v)
{
	this->compute_vert_rad(z);
	this->set_shear_vel( 2.0*this->get_vert_rad()*v.norm() );
}

void Crane::compute_shear_ratio(double m, double x, double s, double w, double T, double P, double z, double u, double E, Matrix<double> v)
{
	this->compute_surf_area(m, x, s, w, T, P, z);
	this->compute_cloud_volume(m, x, s, w, T, P);
	this->compute_char_vel(u, E);
	this->compute_shear_vel(z, v);
	
	if (this->get_includeShearVel() == true)
	{
		this->set_shear_ratio( this->get_mu()*((this->get_surf_area()*this->get_char_vel()/this->get_cloud_volume())+(this->get_k6()*1.5*this->get_shear_vel()/this->get_horz_rad())) );
	}
	else
	{
		this->set_shear_ratio( this->get_mu()*(this->get_surf_area()*this->get_char_vel()/this->get_cloud_volume()) );
	}
}

void Crane::compute_slip_factor(double Dj, double T, double P)
{
	//NOTE: Dj comes in as um --> convert to m
	double dj = Dj/1.0E+6;
	this->compute_air_viscosity(T);
	this->set_slip_factor( 1.0 + (54.088*this->get_air_viscosity()*pow(T, 0.5)/dj/P) );
}

void Crane::compute_davies_num(double Dj, double P, double x, double T)
{
	//NOTE: Dj comes in as um --> convert to m
	double dj = Dj/1.0E+6;
	this->compute_air_density(P, x, T);
	this->compute_air_viscosity(T);
	this->set_davies_num( (4.0*this->get_air_density()*(this->get_part_density()-this->get_air_density())*this->get_grav()*pow(dj,3.0)) / (3.0*pow(this->get_air_viscosity(),2.0)) );
}

void Crane::compute_settling_rate(double Dj, double P, double x, double T)
{
	//NOTE: Dj comes in as um --> convert to m
	double dj = Dj/1.0E+6;
	this->compute_slip_factor(Dj, T, P);
	this->compute_davies_num(Dj, P, x, T);
	
	//If statements for flow conditions
	if (this->get_davies_num() <= 0.3261)
	{
		this->settling_rate[Dj] = this->get_air_viscosity()*this->get_davies_num()*this->get_slip_factor()/24.0/this->get_air_density()/dj;
	}
	else if (this->get_davies_num() <= 84.175)
	{
		double Y = log(this->get_davies_num());
		double exp = -3.18657 + (0.992696*Y) - (0.153193E-2*pow(Y,2.0)) - (0.987059E-3*pow(Y,3.0)) - (0.578878E-3*pow(Y,4.0)) + (0.855176E-4*pow(Y,5.0)) - (0.327815E-5*pow(Y,6.0));
		this->settling_rate[Dj] = this->get_air_viscosity()*exp*this->get_slip_factor()/this->get_air_density()/dj;
	}
	else if (this->get_davies_num() < 140.0)
	{
		double poly = 4.1667E-2 - (2.3363E-4*this->get_davies_num()) + (2.0154E-6*this->get_davies_num()*this->get_davies_num()) - (6.9105E-9*this->get_davies_num()*this->get_davies_num()*this->get_davies_num());
		this->settling_rate[Dj] = this->get_air_viscosity()*poly*this->get_slip_factor()*this->get_davies_num()/this->get_air_density()/dj;
	}
	else if (this->get_davies_num() < 4.5E+7)
	{
		double X = log10(this->get_davies_num());
		double poly = -1.29536 + (0.986*X) - (0.046677*X*X) + (1.1235E-3*X*X*X);
		this->settling_rate[Dj] = this->get_air_viscosity()*pow(10.0,poly)/this->get_air_density()/dj;
	}
	else
	{
		double X = log10(this->get_davies_num());
		double poly = -1.29536 + (0.986*X) - (0.046677*X*X) + (1.1235E-3*X*X*X);
		this->settling_rate[Dj] = this->get_air_viscosity()*pow(10.0,poly)/this->get_air_density()/dj;
	}
}

void Crane::compute_total_mass_fallout_rate(double m, double x, double s, double w, double T, double P, double z, const Matrix<double> &n)
{
	this->compute_horz_rad(m, x, s, w, T, P, z);
	this->settling_rate.clear();
	double sum = 0.0;
	int i = 0;
	
	//Iterate through part_hist map for summation term (NOTE: first N variables in matrix will correspond to particle concentrations)
	for (std::map<double,double>::iterator it=this->part_hist.begin(); it!=this->part_hist.end(); ++it)
	{
		//NOTE: Dj comes in as um --> convert to m
		double dj = it->first/1.0E+6;
		this->compute_settling_rate(it->first, P, x, T);
		sum += this->settling_rate[it->first]*M_PI*dj*dj*dj*n(i,0)*1.0e9/6.0;
		i++;
	}
	sum = sum*M_PI*this->get_horz_rad()*this->get_horz_rad()*this->get_part_density();
	this->set_total_mass_fallout_rate(sum);
}

// Below are listed compute functions specific for initial conditions
void Crane::compute_k(double W)
{
	this->set_k(595.0*pow(W, -0.0527));
}

void Crane::compute_k2(double W)
{
	this->set_k2( fmax(0.004, fmin(0.1, 0.1*pow(W, -(1.0/3.0)))) );
}

void Crane::compute_mu(double W)
{
	this->set_mu( fmax(fmax(0.12, 0.1*pow(W, 0.1)), 0.01*pow(W, (1.0/3.0)))	);
}

void Crane::compute_force_factor(double W)
{
	this->set_force_factor(0.44*pow(W, 0.014));
}

void Crane::compute_equil_temp(double W)
{
	this->set_equil_temp(200.0*log10(W)+1000.0);
}

void Crane::create_part_hist(double min, double max, int size, double avg, double std)
{
	if (max <= min || max <= 0.0 || min <= 0.0 || size < 1 || avg <= 0.0 || std <= 0.0)
	{
		mError(distribution_impossible);
		return;
	}
	
	this->part_hist.clear();
	this->settling_rate.clear();
	this->part_size.clear();
	this->set_min_dia(min);
	this->set_max_dia(max);
	this->set_num_bins(size);
	this->set_mean_dia(avg);
	this->set_std_dia(std);
	this->part_size.resize(size);
	this->part_conc_var.set_size(size, 1);
	
	double distance = log10(this->get_max_dia()) - log10(this->get_min_dia());
	double logstep = distance / ((double)this->get_num_bins());
	double current_log = log10(this->get_min_dia());
	double sum = 0.0;
	
	//Loop to create initial map
	for (int i=0; i<this->get_num_bins(); i++)
	{
		double next_log = current_log + logstep;
		double Dj = sqrt(pow(10.0,current_log)*pow(10.0,next_log));
		double Nj = ( 1.0 / sqrt(2.0*M_PI) / Dj / log(this->get_std_dia()) ) * exp( -0.5*pow( log(Dj/this->get_mean_dia())/log(this->get_std_dia()) ,2.0) );
		this->part_hist[Dj] = Nj*(pow(10.0,next_log)-pow(10.0,current_log));
		this->settling_rate[Dj] = 0.0;
		this->part_size[i] = Dj;
		sum += Nj*(pow(10.0,next_log)-pow(10.0,current_log));
		current_log += logstep;
	}
	
	//Iterate through map to normalize
	for (std::map<double,double>::iterator it=this->part_hist.begin(); it!=this->part_hist.end(); ++it)
	{
		it->second = it->second/sum;
	}
}

void Crane::compute_det_alt(double gz, double hb)
{
	this->set_det_alt(gz+hb);
}

void Crane::compute_initial_cloud_alt(double W, double gz, double hb)
{
	this->compute_det_alt(gz, hb);
	this->set_cloud_alt(this->get_det_alt() + 108.0*pow(W,0.349));
}

void Crane::compute_initial_current_time(double W, double gz, double hb)
{
	this->compute_det_alt(gz, hb);
	double scaled = this->get_det_alt() / pow(W,1.0/3.0);
	double t2m;
	
	if (scaled <= 180.0)
	{
		t2m = 0.037*pow(1.216,scaled/180.0)*pow(W, 0.49 - (0.07*scaled/180.0));
	}
	else
	{
		t2m = 0.045*pow(W, 0.42);
	}
	
	this->set_current_time(56.0*t2m*pow(W, -0.3));
}

void Crane::compute_initial_temperature(double W, double gz, double hb)
{
	this->compute_initial_current_time(W, gz, hb);
	double scaled = this->get_det_alt() / pow(W,1.0/3.0);
	double t2m;
	double K, n;
	
	if (scaled <= 180.0)
	{
		t2m = 0.037*pow(1.216,scaled/180.0)*pow(W, 0.49 - (0.07*scaled/180.0));
		K = 5980.0*pow(1.145, scaled/180.0)*pow(W, -0.0395 + (0.0264*scaled/180.0));
		n = -0.4473*pow(W, 0.0436);
	}
	else
	{
		t2m = 0.045*pow(W, 0.42);
		K = 6847.0*pow(W, -0.0131);
		n = -0.4473*pow(W, 0.0436);
	}
	
	this->set_temperature( K*pow(this->get_current_time()/t2m, n) + 1500.0 );
}

void Crane::compute_initial_soil_mass(double W, double gz, double hb)
{
	this->compute_det_alt(gz, hb);
	double scaled;
	//Check for underground detonation
	if (hb < 0.0)
	{
		scaled = fabs(hb)/pow(W, 1.0/3.4);
		double Rad = 112.5 + 0.755*scaled - 9.6e-6*scaled*scaled*scaled - 9.11e-12*scaled*scaled*scaled*scaled*scaled;
		double D = 32.7 + 0.851*scaled - 2.52e-5*scaled*scaled*scaled - 1.78e-10*scaled*scaled*scaled*scaled*scaled;
		
		this->set_initial_soil_mass( 2.182*pow(W, 3.0/3.4)*Rad*Rad*D );
	}
	else
	{
		scaled = this->get_det_alt() / pow(W, 1.0/3.4);
		
		if (scaled <= 180.0)
		{
			this->set_initial_soil_mass( 0.07741*pow(W, 3.0/3.4)*pow(180.0-scaled,2.0)*(360.0+scaled) );
		}
		else
		{
			this->set_initial_soil_mass( 90.7 );
		}
	}
}

void Crane::compute_initial_part_hist(double W, double gz, double hb, int size)
{
	this->compute_det_alt(gz, hb);
	double scaled = this->get_det_alt() / pow(W, 1.0/3.4);
	
	if (scaled < 180.0)
	{
		this->create_part_hist(0.0001, 100, size, 0.407, 4.0);
	}
	else
	{
		this->create_part_hist(0.0001, 100, size, 0.15, 2.0);
	}
}

void Crane::compute_initial_air_mass(double W, double gz, double hb)
{
	this->compute_initial_cloud_alt(W, gz, hb);
	this->compute_force_factor(W);
	this->compute_initial_soil_mass(W, gz, hb);
	this->compute_initial_temperature(W, gz, hb);
	this->compute_equil_temp(W);
	double Tei = this->return_amb_temp(this->get_cloud_alt());
	double P = this->return_atm_press(this->get_cloud_alt());
	double HR = this->return_rel_humid(this->get_cloud_alt());
	this->compute_xe(Tei, P, HR);
	double cs_int = this->return_spec_heat_conds_integral(this->get_equil_temp(), Tei);
	this->compute_spec_heat_entrain_integral(this->get_temperature(), Tei);
	double cpw_int = this->return_spec_heat_water_integral(this->get_temperature(), Tei);
	
	this->set_initial_air_mass( this->get_energy_frac()*(4.18e12*this->get_force_factor()*W - this->get_initial_soil_mass()*cs_int)/(this->get_spec_heat_entrain_integral() + (this->get_xe()*cpw_int)) );
}

void Crane::compute_initial_water_mass(double W, double gz, double hb)
{
	this->compute_initial_air_mass(W, gz, hb);
	double Tei = this->return_amb_temp(this->get_cloud_alt());
	double cs_int = this->return_spec_heat_conds_integral(this->get_equil_temp(), Tei);
	double cpw_int = this->return_spec_heat_water_integral(this->get_temperature(), Tei);
	
	this->set_initial_water_mass( ( (1.0-this->get_energy_frac())*(4.18e12*this->get_force_factor()*W - this->get_initial_soil_mass()*cs_int)/(cpw_int + this->get_latent_heat()) ) + (this->get_xe()*this->get_initial_air_mass()) );
}

void Crane::compute_initial_entrained_mass(double W, double gz, double hb)
{
	this->compute_initial_water_mass(W, gz, hb);
	this->set_entrained_mass(this->get_initial_air_mass()+this->get_initial_water_mass());
}

void Crane::compute_initial_cloud_mass(double W, double gz, double hb)
{
	this->compute_initial_entrained_mass(W, gz, hb);
	this->set_cloud_mass( (this->get_entrained_mass()+this->get_initial_soil_mass())/1000.0 );
}

void Crane::compute_initial_s_soil(double W, double gz, double hb)
{
	this->compute_initial_air_mass(W, gz, hb);
	this->set_s_soil(this->get_initial_soil_mass()/this->get_initial_air_mass());
}

void Crane::compute_initial_x_water_vapor(double W, double gz, double hb)
{
	this->compute_initial_water_mass(W, gz, hb);
	this->set_x_water_vapor(this->get_initial_water_mass()/this->get_initial_air_mass());
}

void Crane::compute_initial_cloud_volume(double W, double gz, double hb)
{
	this->compute_initial_entrained_mass(W, gz, hb);
	this->compute_initial_x_water_vapor(W, gz, hb);
	this->compute_apparent_temp(this->get_temperature(), this->get_x_water_vapor());
	this->set_cloud_volume( (this->get_initial_air_mass()+this->get_initial_water_mass()) * this->get_gas_const() * this->get_apparent_temp() / this->return_atm_press(this->get_cloud_alt()) );
}

void Crane::compute_initial_horz_rad(double W, double gz, double hb)
{
	this->compute_initial_cloud_volume(W, gz, hb);
	double top = 3.0*this->get_cloud_volume();
	double bot = 4.0*M_PI*sqrt(1.0 - (this->get_eccentricity()*this->get_eccentricity()));
	this->set_horz_rad( pow(top/bot, 1.0/3.0) );
}

void Crane::compute_initial_vert_rad(double W, double gz, double hb)
{
	this->compute_initial_horz_rad(W, gz, hb);
	this->set_vert_rad( sqrt(this->get_horz_rad()*this->get_horz_rad()*(1.0 - (this->get_eccentricity()*this->get_eccentricity()))) );
}

void Crane::compute_initial_cloud_rise(double W, double gz, double hb)
{
	this->compute_initial_horz_rad(W, gz, hb);
	this->set_cloud_rise( 1.2*sqrt(this->get_grav()*this->get_horz_rad()) );
}

void Crane::compute_initial_energy(double W, double gz, double hb)
{
	this->compute_initial_cloud_rise(W, gz, hb);
	this->set_energy(this->get_cloud_rise()*this->get_cloud_rise()/2.0);
}

void Crane::compute_initial_part_conc(double W, double gz, double hb, int size)
{
	this->compute_initial_cloud_volume(W, gz, hb);
	this->compute_initial_part_hist(W, gz, hb, size);
	double conc = this->get_initial_soil_mass() / this->get_cloud_volume();
	
	//Iterate through map
	int i=0;
	for (std::map<double,double>::iterator it=this->part_hist.begin(); it!=this->part_hist.end(); ++it)
	{
		//NOTE: Dj comes in as um --> convert to m
		double dj = it->first/1.0E+6;
		double mass = this->get_part_density() * M_PI * dj*dj*dj / 6.0;
		this->part_conc[it->first] = it->second * conc / mass / 1.0e9;
		this->part_conc_var(i,0) = this->part_conc[it->first];
		i++;
	}
}

void Crane::compute_initial_virtual_mass(double W, double gz, double hb)
{
	this->compute_initial_cloud_mass(W, gz, hb);
	this->compute_initial_s_soil(W, gz, hb);
	this->compute_initial_x_water_vapor(W, gz, hb);
	this->compute_apparent_temp(this->get_temperature(), this->get_x_water_vapor());
	double Tei = this->return_amb_temp(this->get_cloud_alt());
	this->compute_apparent_amb_temp(Tei, this->get_xe());
	this->compute_beta_prime(this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds());
	
	this->set_virtual_mass( this->get_cloud_mass()*this->get_beta_prime()*this->get_apparent_temp()/2.0/this->get_apparent_amb_temp() );
}

void Crane::compute_adjusted_height(double W, double gz, double hb)
{
	this->compute_initial_vert_rad(W, gz, hb);
	this->compute_mu(W);
	this->compute_initial_cloud_alt(W, gz, hb);
	this->set_adjusted_height( this->get_cloud_alt() - (this->get_vert_rad()/this->get_mu()) );
}

void Crane::delete_particles()
{
	part_hist.clear();
	settling_rate.clear();
	part_conc.clear();
	part_size.clear();
}

// Below are listed return functions specific for temperature integral related values

void Crane::compute_spec_heat_entrain_integral(double T, double Te)
{
	//Integration by parts for a piecewise polynomial
	double upper = 0.0, lower = 0.0;
	if (T > 2300.0)
		upper = (-3587.5*T + (2.125/2.0)*T*T) - (-3587.5*2300.0 + (2.125/2.0)*2300.0*2300.0);
	else
		upper = (946.6*T + (0.1971/2.0)*T*T) - (946.6*2300.0 + (0.1971/2.0)*2300.0*2300.0);
	
	if (Te > 2300.0)
		lower = (-3587.5*2300.0 + (2.125/2.0)*2300.0*2300.0) - (-3587.5*Te + (2.125/2.0)*Te*Te);
	else
		lower = (946.6*2300.0 + (0.1971/2.0)*2300.0*2300.0) - (946.6*Te + (0.1971/2.0)*Te*Te);
	
	this->set_spec_heat_entrain_integral(upper + lower);
}

double Crane::return_spec_heat_water_integral(double T, double Te)
{
	return (1697.66*T + (1.144174/2.0)*T*T) - (1697.66*Te + (1.144174/2.0)*Te*Te);
}

double Crane::return_spec_heat_conds_integral(double T, double Te)
{
	//Integration by parts for a piecewise polynomial
	double upper = 0.0, lower = 0.0;
	if (T > 848.0)
		upper = (1003.8*T + (0.1351/2.0)*T*T) - (1003.8*848.0 + (0.1351/2.0)*848.0*848.0);
	else
		upper = (781.6*T + (0.5612/2.0)*T*T + 1.881e+7/T) - (781.6*848.0 + (0.5612/2.0)*848.0*848.0 + 1.881e+7/848.0);
	
	if (Te > 848.0)
		lower = (1003.8*848.0 + (0.1351/2.0)*848.0*848.0) - (1003.8*Te + (0.1351/2.0)*Te*Te);
	else
		lower = (781.6*848.0 + (0.5612/2.0)*848.0*848.0 + 1.881e+7/848.0) - (781.6*Te + (0.5612/2.0)*Te*Te + 1.881e+7/Te);
	
	return upper+lower;
}

// Below are listed functions specific for air profile related operations

void Crane::add_amb_temp(double z, double Te)
{
	this->amb_temp[z] = Te;
}

void Crane::add_atm_press(double z, double P)
{
	this->atm_press[z] = P;
}

void Crane::add_rel_humid(double z, double HR)
{
	this->rel_humid[z] = HR;
}

void Crane::add_wind_vel(double z, double vx, double vy)
{
	Matrix<double> v(2,1);
	v(0,0) = vx;
	v(1,0) = vy;
	this->wind_vel[z] = v;
}

void Crane::create_default_atmosphere()
{
	this->delete_atmosphere();
	
	this->add_amb_temp(-1000, 21.5+273.15);
	this->add_atm_press(-1000, 11.39*10000.0);
	this->add_rel_humid(-1000, 0.86);
	this->add_wind_vel(-1000, 8.19, 0.0);
	
	this->add_amb_temp(0, 15.0+273.15);
	this->add_atm_press(0, 10.13*10000.0);
	this->add_rel_humid(0, 0.91);
	this->add_wind_vel(0, 9.04, 0.0);
	
	this->add_amb_temp(1000, 8.5+273.15);
	this->add_atm_press(1000, 8.988*10000.0);
	this->add_rel_humid(1000, 1.29);
	this->add_wind_vel(1000, 9.58, 0.0);
	
	this->add_amb_temp(2000, 2.0+273.15);
	this->add_atm_press(2000, 7.95*10000.0);
	this->add_rel_humid(2000, 3.21);
	this->add_wind_vel(2000, 9.98, 0.0);
	
	this->add_amb_temp(3000, -4.49+273.15);
	this->add_atm_press(3000, 7.012*10000.0);
	this->add_rel_humid(3000, 4.80);
	this->add_wind_vel(3000, 14.06, 0.0);
	
	this->add_amb_temp(4000, -10.98+273.15);
	this->add_atm_press(4000, 6.166*10000.0);
	this->add_rel_humid(4000, 5.39);
	this->add_wind_vel(4000, 18.14, 0.0);
	
	this->add_amb_temp(5000, -17.47+273.15);
	this->add_atm_press(5000, 5.405*10000.0);
	this->add_rel_humid(5000, 13.83);
	this->add_wind_vel(5000, 20.68, 0.0);
	
	this->add_amb_temp(6000, -23.96+273.15);
	this->add_atm_press(6000, 4.722*10000.0);
	this->add_rel_humid(6000, 27.61);
	this->add_wind_vel(6000, 23.21, 0.0);
	
	this->add_amb_temp(7000, -30.45+273.15);
	this->add_atm_press(7000, 4.111*10000.0);
	this->add_rel_humid(7000, 41.67);
	this->add_wind_vel(7000, 26.37, 0.0);
	
	this->add_amb_temp(8000, -36.94+273.15);
	this->add_atm_press(8000, 3.565*10000.0);
	this->add_rel_humid(8000, 49.80);
	this->add_wind_vel(8000, 29.52, 0.0);
	
	this->add_amb_temp(9000, -43.42+273.15);
	this->add_atm_press(9000, 3.080*10000.0);
	this->add_rel_humid(9000, 54.80);
	this->add_wind_vel(9000, 32.37, 0.0);
	
	this->add_amb_temp(10000, -49.90+273.15);
	this->add_atm_press(10000, 2.650*10000.0);
	this->add_rel_humid(10000, 79.90);
	this->add_wind_vel(10000, 35.21, 0.0);
	
	this->add_amb_temp(15000, -56.5+273.15);
	this->add_atm_press(15000, 1.211*10000.0);
	this->add_rel_humid(15000, 83.27);
	this->add_wind_vel(15000, 28.64, 0.0);
	
	this->add_amb_temp(20000, -56.5+273.15);
	this->add_atm_press(20000, 0.5529*10000.0);
	this->add_rel_humid(20000, 46.74);
	this->add_wind_vel(20000, 12.26, 0.0);
	
	this->add_amb_temp(25000, -51.6+273.15);
	this->add_atm_press(25000, 0.2549*10000.0);
	this->add_rel_humid(25000, 16.33);
	this->add_wind_vel(25000, 7.51, 0.0);
	
	this->add_amb_temp(30000, -46.64+273.15);
	this->add_atm_press(30000, 0.1197*10000.0);
	this->add_rel_humid(30000, 3.02);
	this->add_wind_vel(30000, 20.49, 0.0);
	
	this->add_amb_temp(40000, -22.80+273.15);
	this->add_atm_press(40000, 0.0287*10000.0);
	this->add_rel_humid(40000, 0.15);
	this->add_wind_vel(40000, 61.79, 0.0);
	
	this->add_amb_temp(50000, -2.5+273.15);
	this->add_atm_press(50000, 0.007978*10000.0);
	this->add_rel_humid(50000, 0.003);
	this->add_wind_vel(50000, 88.44, 0.0);
	
	this->add_amb_temp(60000, -26.13+273.15);
	this->add_atm_press(60000, 0.002196*10000.0);
	this->add_rel_humid(60000, 0.009);
	this->add_wind_vel(60000, 74.65, 0.0);
	
	this->add_amb_temp(70000, -53.57+273.15);
	this->add_atm_press(70000, 0.00052*10000.0);
	this->add_rel_humid(70000, 0.65);
	this->add_wind_vel(70000, 32.79, 0.0);
	
	this->add_amb_temp(80000, -74.51+273.15);
	this->add_atm_press(80000, 0.00011*10000.0);
	this->add_rel_humid(80000, 1.28);
	this->add_wind_vel(80000, 56.37, 0.0);
}

void Crane::delete_atmosphere()
{
	this->amb_temp.clear();
	this->atm_press.clear();
	this->rel_humid.clear();
	this->wind_vel.clear();
}

double Crane::return_amb_temp(double z)
{
	double Te = 0.0;
	
	//Setup the iterators
	std::map<double,double>::iterator it=this->amb_temp.begin();
	std::map<double,double>::reverse_iterator rit=this->amb_temp.rbegin();
	
	//Special Case 1: z less than lowest value in map
	if (z <= it->first)
		return it->second;
	
	//Special Case 2: z greater than highest value in map
	if (z >= rit->first)
		return rit->second;
	
	//Iterate through map
	double old_z = it->first;
	double old_Te = it->second;
	for (it=this->amb_temp.begin(); it!=this->amb_temp.end(); ++it)
	{
		if (it->first > z)
		{
			double slope = (it->second - old_Te) / (it->first - old_z);
			double inter = it->second - (slope*it->first);
			Te = slope*z + inter;
			break;
		}
		else
		{
			old_z = it->first;
			old_Te = it->second;
		}
	}
	
	return Te;
}

double Crane::return_atm_press(double z)
{
	double P = 0.0;
	
	//Setup the iterators
	std::map<double,double>::iterator it=this->atm_press.begin();
	std::map<double,double>::reverse_iterator rit=this->atm_press.rbegin();
	
	//Special Case 1: z less than lowest value in map
	if (z <= it->first)
		return it->second;
	
	//Special Case 2: z greater than highest value in map
	if (z >= rit->first)
		return rit->second;
	
	//Iterate through map
	double old_z = it->first;
	double old_Te = it->second;
	for (it=this->atm_press.begin(); it!=this->atm_press.end(); ++it)
	{
		if (it->first > z)
		{
			double slope = (it->second - old_Te) / (it->first - old_z);
			double inter = it->second - (slope*it->first);
			P = slope*z + inter;
			break;
		}
		else
		{
			old_z = it->first;
			old_Te = it->second;
		}
	}
	
	return P;
}

double Crane::return_rel_humid(double z)
{
	double HR = 0.0;
	
	//Setup the iterators
	std::map<double,double>::iterator it=this->rel_humid.begin();
	std::map<double,double>::reverse_iterator rit=this->rel_humid.rbegin();
	
	//Special Case 1: z less than lowest value in map
	if (z <= it->first)
		return it->second;
	
	//Special Case 2: z greater than highest value in map
	if (z >= rit->first)
		return rit->second;
	
	//Iterate through map
	double old_z = it->first;
	double old_Te = it->second;
	for (it=this->rel_humid.begin(); it!=this->rel_humid.end(); ++it)
	{
		if (it->first > z)
		{
			double slope = (it->second - old_Te) / (it->first - old_z);
			double inter = it->second - (slope*it->first);
			HR = slope*z + inter;
			break;
		}
		else
		{
			old_z = it->first;
			old_Te = it->second;
		}
	}
	
	return HR;
}

Matrix<double> Crane::return_wind_vel(double z)
{
	Matrix<double> HR;
	
	//Setup the iterators
	std::map<double,Matrix<double>>::iterator it=this->wind_vel.begin();
	std::map<double,Matrix<double>>::reverse_iterator rit=this->wind_vel.rbegin();
	
	//Special Case 1: z less than lowest value in map
	if (z <= it->first)
		return it->second;
	
	//Special Case 2: z greater than highest value in map
	if (z >= rit->first)
		return rit->second;
	
	//Iterate through map
	double old_z = it->first;
	Matrix<double> old_Te = it->second;
	for (it=this->wind_vel.begin(); it!=this->wind_vel.end(); ++it)
	{
		if (it->first > z)
		{
			Matrix<double> slope = (it->second - old_Te) / (it->first - old_z);
			Matrix<double> inter = it->second - (slope*it->first);
			HR = slope*z + inter;
			break;
		}
		else
		{
			old_z = it->first;
			old_Te = it->second;
		}
	}
	
	return HR;
}

double Crane::return_wind_speed(double z)
{
	Matrix<double> v = this->return_wind_vel(z);
	return v.norm();
}

void Crane::compute_current_amb_temp(double z)
{
	this->set_current_amb_temp( this->return_amb_temp(z) );
}

void Crane::compute_current_atm_press(double z)
{
	this->set_current_atm_press( this->return_atm_press(z) );
}

// Below are listed functions to feed to DOVE as residuals

double rate_cloud_rise(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	double m = u(dove.getVariableIndex("m (Mg)"),0);
	double U = u(dove.getVariableIndex("u (m/s)"),0);
	double dm_dt = dove.coupledTimeDerivative("m (Mg)",u);
	double mass_rat = m / (m + dat->get_virtual_mass());
	double T = u(dove.getVariableIndex("T (K)"),0);
	double x = u(dove.getVariableIndex("x (kg/kg)"),0);
	double s = u(dove.getVariableIndex("s (kg/kg)"),0);
	double w = u(dove.getVariableIndex("w (kg/kg)"),0);
	double z = u(dove.getVariableIndex("z (m)"),0);
	
	if (dat->get_isTight() == true)
	{
		dat->compute_apparent_temp(T, x); //NOTE: be aware of potential nan or inf residuals
		dat->compute_beta_prime(x, s, w); //NOTE: be aware of potential nan or inf residuals
		dat->compute_vert_rad(z);
	}
	
	double p1 = (((dat->get_apparent_temp()/dat->get_apparent_amb_temp())*dat->get_beta_prime()) - 1.0) * (dat->get_grav()/(1.0-dat->get_mu()));
	double p2 = (2.0*dat->get_k2()*dat->get_char_vel()*dat->get_apparent_temp()*dat->get_beta_prime()*(1.0-dat->get_mu()))/(dat->get_vert_rad()*dat->get_apparent_amb_temp());
	double p3 = (1.0/m)*dm_dt;
	
	res = (p1 - ((p2+p3)*U))*mass_rat;
	
	return res;
}

double rate_cloud_alt(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	return u(dove.getVariableIndex("u (m/s)"),0);
}

double rate_x_water_vapor(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	// Check for saturation
	if (dat->get_isSaturated() == true)
	{
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double T = u(dove.getVariableIndex("T (K)"),0);
		double U = u(dove.getVariableIndex("u (m/s)"),0);
		double dT_dt = dove.coupledTimeDerivative("T (K)",u);
		
		double p1 = (1.0 + (x/dat->get_eps()))*(dat->get_latent_heat()*dat->get_eps()/dat->get_gas_const()/T/T)*dT_dt;
		double p2 = (1.0 + (x/dat->get_eps()))*dat->get_grav()*U/dat->get_gas_const()/dat->get_apparent_amb_temp();
		
		res = x*(p1+p2);
	}
	else
	{
		double m = u(dove.getVariableIndex("m (Mg)"),0);
		double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double s = u(dove.getVariableIndex("s (kg/kg)"),0);
		
		res = -((1.0+x+s)/(1.0+dat->get_xe()))*(x - dat->get_xe())*dment_dt/m;
	}
	
	return res;
}

double rate_temperature(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	// Check for saturation
	if (dat->get_isSaturated() == true)
	{
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double T = u(dove.getVariableIndex("T (K)"),0);
		double m = u(dove.getVariableIndex("m (Mg)"),0);
		double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
		double E = u(dove.getVariableIndex("E (J/kg)"),0);
		double s = u(dove.getVariableIndex("s (kg/kg)"),0);
		double w = u(dove.getVariableIndex("w (kg/kg)"),0);
		double U = u(dove.getVariableIndex("u (m/s)"),0);
		double z = u(dove.getVariableIndex("z (m)"),0);
		
		if (dat->get_isTight() == true)
		{
			dat->compute_apparent_temp(T, x); //NOTE: be aware of potential nan or inf residuals
			dat->compute_beta_prime(x, s, w); //NOTE: be aware of potential nan or inf residuals
			dat->compute_sigma_turbulence(E, z);
			dat->compute_actual_spec_heat(T, x);
		}
		
		double p1 = dat->get_beta_prime() / (1.0 + (dat->get_latent_heat()*dat->get_latent_heat()*x*dat->get_eps()/dat->get_actual_spec_heat()/dat->get_gas_const()/T/T));
		double p2 = ( T - dat->get_current_amb_temp() + (dat->get_latent_heat()*(x-dat->get_xe())/dat->get_actual_spec_heat()) ) * dment_dt/m/dat->get_beta_prime();
		double p3 = (dat->get_apparent_temp()*dat->get_grav()*U/dat->get_apparent_amb_temp()/dat->get_actual_spec_heat())*(1.0+(dat->get_latent_heat()*x/dat->get_gas_const()/T));
		double p4 = dat->get_sigma_turbulence()/dat->get_actual_spec_heat();
		
		res = -p1*(p2+p3-p4);
		
	}
	else
	{
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double T = u(dove.getVariableIndex("T (K)"),0);
		double m = u(dove.getVariableIndex("m (Mg)"),0);
		double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
		double E = u(dove.getVariableIndex("E (J/kg)"),0);
		double s = u(dove.getVariableIndex("s (kg/kg)"),0);
		double w = u(dove.getVariableIndex("w (kg/kg)"),0);
		double U = u(dove.getVariableIndex("u (m/s)"),0);
		double z = u(dove.getVariableIndex("z (m)"),0);
		
		if (dat->get_isTight() == true)
		{
			dat->compute_apparent_temp(T, x); //NOTE: be aware of potential nan or inf residuals
			dat->compute_beta_prime(x, s, w); //NOTE: be aware of potential nan or inf residuals
			dat->compute_sigma_turbulence(E, z);
			dat->compute_mean_spec_heat(T, x, s, w);
		}
		
		double p1 = dat->get_beta_prime()/dat->get_mean_spec_heat();
		double p2 = dat->get_apparent_temp()*dat->get_grav()*U/dat->get_apparent_amb_temp();
		double p3 = dat->get_spec_heat_entrain_integral()*dment_dt/dat->get_beta_prime()/m;
		
		res = -p1*(p2 + p3 - dat->get_sigma_turbulence());
	}
	
	return res;
}

double rate_w_water_conds(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	// Check for saturation
	if (dat->get_isSaturated() == true)
	{
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double T = u(dove.getVariableIndex("T (K)"),0);
		double m = u(dove.getVariableIndex("m (Mg)"),0);
		double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
		double dx_dt = dove.coupledTimeDerivative("x (kg/kg)",u);
		double s = u(dove.getVariableIndex("s (kg/kg)"),0);
		double w = u(dove.getVariableIndex("w (kg/kg)"),0);
		double z = u(dove.getVariableIndex("z (m)"),0);
		
		if (dat->get_isTight() == true)
		{
			dat->compute_beta_prime(x, s, w); //NOTE: be aware of potential nan or inf residuals
			dat->compute_total_mass_fallout_rate(m, x, s, w, T, dat->get_current_atm_press(), z, dat->get_part_conc_var());
		}
		
		double p1 = (1.0/dat->get_beta_prime())*((1.0+x)/(1.0+dat->get_xe()))*(w+x-dat->get_xe())*dment_dt/m;
		double p2 = ((1.0+x+s+w)/m)*(w/(s+w))*dat->get_total_mass_fallout_rate()/1000.0;
		
		res = -p1 - dx_dt - p2;
	}
	else
	{
		res = 0.0;
	}
	
	return res;
}

double rate_energy(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	double x = u(dove.getVariableIndex("x (kg/kg)"),0);
	double T = u(dove.getVariableIndex("T (K)"),0);
	double m = u(dove.getVariableIndex("m (Mg)"),0);
	double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
	double E = u(dove.getVariableIndex("E (J/kg)"),0);
	double s = u(dove.getVariableIndex("s (kg/kg)"),0);
	double w = u(dove.getVariableIndex("w (kg/kg)"),0);
	double U = u(dove.getVariableIndex("u (m/s)"),0);
	double z = u(dove.getVariableIndex("z (m)"),0);
	
	if (dat->get_isTight() == true)
	{
		dat->compute_apparent_temp(T, x);
		dat->compute_beta_prime(x, s, w);
		dat->compute_vert_rad(z);
		dat->compute_sigma_turbulence(E, z);
	}
	
	double p1 = 2.0*dat->get_k2()*dat->get_apparent_temp()*dat->get_beta_prime()*U*U*dat->get_char_vel()/dat->get_apparent_amb_temp()/dat->get_vert_rad();
	double p2 = U*U*dment_dt/2.0/m;
	double p3 = E*dment_dt/m;
	
	res = p1+p2-p3-dat->get_sigma_turbulence();
	
	return res;
}

double rate_cloud_mass(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	double x = u(dove.getVariableIndex("x (kg/kg)"),0);
	double T = u(dove.getVariableIndex("T (K)"),0);
	double m = u(dove.getVariableIndex("m (Mg)"),0);
	double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
	double s = u(dove.getVariableIndex("s (kg/kg)"),0);
	double w = u(dove.getVariableIndex("w (kg/kg)"),0);
	double z = u(dove.getVariableIndex("z (m)"),0);
	
	if (dat->get_isTight() == true)
		dat->compute_total_mass_fallout_rate(m, x, s, w, T, dat->get_current_atm_press(), z, dat->get_part_conc_var());
	
	res = dment_dt - dat->get_total_mass_fallout_rate()/1000.0;
	
	return res;
}

double rate_s_soil(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	double x = u(dove.getVariableIndex("x (kg/kg)"),0);
	double T = u(dove.getVariableIndex("T (K)"),0);
	double m = u(dove.getVariableIndex("m (Mg)"),0);
	double dment_dt = rate_entrained_mass(0, u, dove.getCurrentTime(), data, dove);
	double s = u(dove.getVariableIndex("s (kg/kg)"),0);
	double w = u(dove.getVariableIndex("w (kg/kg)"),0);
	double z = u(dove.getVariableIndex("z (m)"),0);
	
	if (dat->get_isTight() == true)
	{
		dat->compute_total_mass_fallout_rate(m, x, s, w, T, dat->get_current_atm_press(), z, dat->get_part_conc_var());
		dat->compute_beta_prime(x, s, w);
	}
	
	double p1 = ((1.0+x)/(1.0+dat->get_xe()))*s*dment_dt/m/dat->get_beta_prime();
	double p2 = ((1.0+x+s+w)/m)*(s/(s+w))*dat->get_total_mass_fallout_rate()/1000.0;
	
	res = -p1-p2;
	
	return res;
}

double rate_entrained_mass(int i, const Matrix<double> &u, double t, const void *data, const Dove &dove)
{
	double res = 0.0;
	
	Crane *dat = (Crane *) data;
	
	// Check for saturation
	if (dat->get_isSaturated() == true)
	{
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double T = u(dove.getVariableIndex("T (K)"),0);
		double m = u(dove.getVariableIndex("m (Mg)"),0);
		double E = u(dove.getVariableIndex("E (J/kg)"),0);
		double s = u(dove.getVariableIndex("s (kg/kg)"),0);
		double w = u(dove.getVariableIndex("w (kg/kg)"),0);
		double U = u(dove.getVariableIndex("u (m/s)"),0);
		double z = u(dove.getVariableIndex("z (m)"),0);
		
		if (dat->get_isTight() == true)
		{
			dat->compute_apparent_temp(T, x);
			dat->compute_beta_prime(x, s, w);
			dat->compute_sigma_turbulence(E, z);
			dat->compute_actual_spec_heat(T, x);
		}
		
		double p1 = T - dat->get_current_amb_temp() + (dat->get_latent_heat()*(x-dat->get_xe())/dat->get_actual_spec_heat());
		p1 = p1 / (1.0+(dat->get_latent_heat()*dat->get_latent_heat()*x*dat->get_eps()/dat->get_actual_spec_heat()/dat->get_gas_const()/T/T));
		p1 = dat->get_beta_prime()*m/(1.0 - (dat->get_beta_prime()*p1/dat->get_apparent_temp()));
		
		double p2 = ((dat->get_grav()*U*dat->get_apparent_temp()/dat->get_actual_spec_heat()/dat->get_apparent_amb_temp())*(1.0+(dat->get_latent_heat()*x/dat->get_gas_const()/T))) - (dat->get_sigma_turbulence()/dat->get_actual_spec_heat());
		p2 = p2 / (1.0+(dat->get_latent_heat()*dat->get_latent_heat()*x*dat->get_eps()/dat->get_actual_spec_heat()/dat->get_gas_const()/T/T));
		p2 = dat->get_beta_prime()*p2/dat->get_apparent_temp();
		
		double p3 = dat->get_grav()*U/dat->get_gas_const()/dat->get_apparent_amb_temp();
		
		res = p1*(dat->get_shear_ratio()+p2-p3);
	}
	else
	{
		double x = u(dove.getVariableIndex("x (kg/kg)"),0);
		double T = u(dove.getVariableIndex("T (K)"),0);
		double m = u(dove.getVariableIndex("m (Mg)"),0);
		double E = u(dove.getVariableIndex("E (J/kg)"),0);
		double s = u(dove.getVariableIndex("s (kg/kg)"),0);
		double w = u(dove.getVariableIndex("w (kg/kg)"),0);
		double U = u(dove.getVariableIndex("u (m/s)"),0);
		double z = u(dove.getVariableIndex("z (m)"),0);
		
		if (dat->get_isTight() == true)
		{
			dat->compute_apparent_temp(T, x);
			dat->compute_beta_prime(x, s, w);
			dat->compute_sigma_turbulence(E, z);
			dat->compute_mean_spec_heat(T, x, s, w);
		}
		
		double p1 = dat->get_beta_prime()*dat->get_spec_heat_entrain_integral()/dat->get_apparent_temp()/dat->get_mean_spec_heat();
		p1 = dat->get_beta_prime()*m/(1.0 - p1);
		
		double p2 = (dat->get_apparent_temp()*dat->get_grav()*U/dat->get_apparent_amb_temp()) - dat->get_sigma_turbulence();
		p2 = dat->get_beta_prime()*p2/dat->get_apparent_temp()/dat->get_mean_spec_heat();
		
		double p3 = dat->get_grav()*U/dat->get_gas_const()/dat->get_apparent_amb_temp();
		
		res = p1*(dat->get_shear_ratio()+p2-p3);
	}
	
	return res;
}

// Below are list functions associated with actions completed outside of the solver in DOVE

void Crane::establish_initial_conditions(double W, double gz, double hb, int bins, bool includeShear, bool isTight, Dove &dove)
{
	this->compute_initial_current_time(W, gz, hb);
	this->compute_initial_temperature(W, gz, hb);
	this->compute_initial_cloud_alt(W, gz, hb);
	this->compute_initial_cloud_rise(W, gz, hb);
	this->compute_initial_cloud_mass(W, gz, hb);
	this->compute_initial_x_water_vapor(W, gz, hb);
	this->compute_initial_s_soil(W, gz, hb);
	this->compute_initial_energy(W, gz, hb);
	this->compute_initial_part_conc(W, gz, hb, bins);
	this->compute_initial_virtual_mass(W, gz, hb);
	this->compute_initial_entrained_mass(W, gz, hb);
	this->compute_initial_cloud_volume(W, gz, hb);
	this->compute_initial_horz_rad(W, gz, hb);
	this->compute_initial_vert_rad(W, gz, hb);
	this->compute_adjusted_height(W, gz, hb);
	this->compute_equil_temp(W);
	this->compute_k(W);
	this->compute_k2(W);
	this->compute_mu(W);
	this->set_includeShearVel(includeShear);
	this->set_isTight(isTight);
	this->set_isSaturated(false);
	
	///*** Override ICs **///
	this->set_energy( 1993.76 );
	this->set_cloud_alt( 1645.96 );
	this->set_cloud_mass( 3974.0 );
	this->set_s_soil(1.35e-5);
	this->set_temperature(2803.7);
	this->set_cloud_rise(65.01);
	this->set_w_water_conds(0);
	this->set_x_water_vapor(0.002626);
	
	// Setup data
	dove.set_userdata(this);
	dove.set_numfunc(8);
	
	// Name variables

	dove.set_variableName(0, "z (m)");
	dove.set_variableName(1, "w (kg/kg)");
	dove.set_variableName(2, "x (kg/kg)");
	dove.set_variableName(3, "s (kg/kg)");
	dove.set_variableName(4, "u (m/s)");
	dove.set_variableName(5, "m (Mg)");
	dove.set_variableName(6, "E (J/kg)");
	dove.set_variableName(7, "T (K)");
	
	// Register rate functions
	
	dove.registerFunction("z (m)", rate_cloud_alt);
    dove.registerCoeff("z (m)", default_coeff);
	dove.registerFunction("w (kg/kg)", rate_w_water_conds);
    dove.registerCoeff("w (kg/kg)", default_coeff);
	dove.registerFunction("x (kg/kg)", rate_x_water_vapor);
    dove.registerCoeff("x (kg/kg)", default_coeff);
	dove.registerFunction("s (kg/kg)", rate_s_soil);
    dove.registerCoeff("s (kg/kg)", default_coeff);
	dove.registerFunction("u (m/s)", rate_cloud_rise);
    dove.registerCoeff("u (m/s)", default_coeff);
	dove.registerFunction("m (Mg)", rate_cloud_mass);
    dove.registerCoeff("m (Mg)", default_coeff);
	dove.registerFunction("E (J/kg)", rate_energy);
    dove.registerCoeff("E (J/kg)", default_coeff);
	dove.registerFunction("T (K)", rate_temperature);
    dove.registerCoeff("T (K)", default_coeff);
	
	// Set initial conditions
	
	dove.set_starttime(this->get_current_time());
	dove.set_initialcondition("z (m)", this->get_cloud_alt());
	dove.set_initialcondition("w (kg/kg)", this->get_w_water_conds());
	dove.set_initialcondition("x (kg/kg)", this->get_x_water_vapor());
	dove.set_initialcondition("s (kg/kg)", this->get_s_soil());
	dove.set_initialcondition("u (m/s)", this->get_cloud_rise());
	dove.set_initialcondition("m (Mg)", this->get_cloud_mass());
	dove.set_initialcondition("E (J/kg)", this->get_energy());
	dove.set_initialcondition("T (K)", this->get_temperature());
	dove.set_timestep(1.0/this->get_cloud_rise());
}

void Crane::establish_dove_options(Dove &dove, FILE *file, bool fileout, bool consoleout, integrate_subtype inttype, timestep_type timetype,
								   precond_type type, double tol, double dtmin, double dtmax, double dtmin_conv, double t_out, double endtime)
{
	dove.set_outputfile(file);
	dove.set_fileoutput(fileout);
	dove.set_headeroutput(fileout);
	this->set_FileOut(fileout);
	dove.set_output(consoleout);
	this->set_ConsoleOut(consoleout);
	dove.set_integrationtype(inttype);
	dove.set_timestepper(timetype);
	dove.set_preconditioner(type);
	dove.set_tolerance(tol);
	dove.set_timestepmin(dtmin);
	dove.set_timestepmax(dtmax);
	dove.set_timestepmin_converged(dtmin_conv);
	dove.set_t_out(t_out);
	dove.set_endtime(endtime);
}

void Crane::establish_pjfnk_options(Dove &dove, krylov_method lin_method, linesearch_type linesearch, bool linear, bool precon, bool nl_out,
		bool l_out, int max_nlit, int max_lit, int restart, int recursive, double nl_abstol, double nl_reltol, double l_abstol, double l_reltol)
{
	dove.set_LinearMethod(lin_method);
	dove.set_LineSearchMethod(linesearch);
	dove.set_LinearStatus(linear);
	dove.set_Preconditioning(precon);
	dove.set_NonlinearOutput(nl_out);
	dove.set_LinearOutput(l_out);
	dove.set_MaxNonLinearIterations(max_nlit);
	dove.set_MaxLinearIterations(max_lit);
	dove.set_RestartLimit(restart);
	dove.set_RecursionLevel(recursive);
	dove.set_NonlinearAbsTol(nl_abstol);
	dove.set_NonlinearRelTol(nl_reltol);
	dove.set_LinearAbsTol(l_abstol);
	dove.set_LinearRelTol(l_reltol);
}

void Crane::estimate_parameters(Dove &dove)
{
	this->compute_beta_prime(this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds());
	this->compute_q_x(this->get_x_water_vapor());
	
	double Te, P, HR;
	Matrix<double> v;
	Te = this->return_amb_temp(this->get_cloud_alt());
	P = this->return_atm_press(this->get_cloud_alt());
	HR = this->return_rel_humid(this->get_cloud_alt());
	v = this->return_wind_vel(this->get_cloud_alt());
	
	this->set_current_amb_temp(Te);
	this->set_current_atm_press(P);
	this->compute_xe(Te, P, HR);
	
	this->compute_q_xe(this->get_xe());
	this->compute_apparent_temp(this->get_temperature(), this->get_x_water_vapor());
	this->compute_apparent_amb_temp(Te, this->get_xe());
	this->compute_char_vel(this->get_cloud_rise(), this->get_energy());
	this->compute_air_viscosity(this->get_temperature());
	this->compute_vapor_pressure(P, this->get_x_water_vapor());
	this->compute_sat_vapor_pressure(this->get_temperature());
	
	if (this->get_vapor_pressure() > this->get_sat_vapor_pressure())
	{
		this->set_isSaturated(true);
		if (this->get_saturation_time() <= 0.0)
			this->set_saturation_time(this->get_current_time());
	}
	else
		this->set_isSaturated(false);
	
	this->compute_air_density(P, this->get_x_water_vapor(), this->get_temperature());
	this->compute_spec_heat_entrain(this->get_temperature());
	this->compute_spec_heat_water(this->get_temperature());
	this->compute_spec_heat_conds(this->get_temperature());
	this->compute_actual_spec_heat(this->get_temperature(), this->get_x_water_vapor());
	this->compute_k_temp(this->get_temperature());
	this->compute_mean_spec_heat(this->get_temperature(), this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds());
	this->compute_cloud_volume(this->get_cloud_mass(), this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds(), this->get_temperature(), P);
	this->compute_vert_rad(this->get_cloud_alt());
	this->compute_horz_rad(this->get_cloud_mass(), this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds(), this->get_temperature(), P, this->get_cloud_alt());
	this->compute_sigma_turbulence(this->get_energy(), this->get_cloud_alt());
	this->compute_surf_area(this->get_cloud_mass(), this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds(), this->get_temperature(), P, this->get_cloud_alt());
	this->compute_shear_vel(this->get_cloud_alt(), v);
	this->compute_shear_ratio(this->get_cloud_mass(), this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds(), this->get_temperature(), P, this->get_cloud_alt(), this->get_cloud_rise(), this->get_energy(), v);
	
	this->compute_total_mass_fallout_rate(this->get_cloud_mass(), this->get_x_water_vapor(), this->get_s_soil(), this->get_w_water_conds(), this->get_temperature(), P, this->get_cloud_alt(), this->get_part_conc_var());
	
	for (int i=0; i<this->part_conc_var.rows(); i++)
	{
		double rate = M_PI*this->get_horz_rad()*this->get_horz_rad()*this->get_settling_rate(this->get_part_size(i))/this->get_cloud_volume();
		
		rate = exp(-rate*dove.getTimeStep());
		
		this->part_conc_var(i,0) = this->part_conc_var(i,0)*rate;
	}
}

void Crane::store_variables(Dove &dove)
{
	this->set_current_time(dove.getCurrentTime());
	this->set_cloud_mass( dove.getNewU("m (Mg)", dove.getNewU()) );
	this->set_cloud_rise( dove.getNewU("u (m/s)", dove.getNewU()) );
	this->set_cloud_alt( dove.getNewU("z (m)", dove.getNewU()) );
	this->set_x_water_vapor( dove.getNewU("x (kg/kg)", dove.getNewU()) );
	if (this->get_isSaturated() == false)
		dove.getNewU()(dove.getVariableIndex("w (kg/kg)"),0) = 0.0;
	this->set_w_water_conds( dove.getNewU("w (kg/kg)", dove.getNewU()) );
	this->set_s_soil(dove.getNewU("s (kg/kg)", dove.getNewU()));
	this->set_temperature(dove.getNewU("T (K)", dove.getNewU()));
	this->set_energy( dove.getNewU("E (J/kg)", dove.getNewU()) );
}

int Crane::run_crane_simulation(Dove &dove)
{
	int success = 0;
	
	dove.reset_all();
	
	if (this->get_FileOut() == true)
	{
		dove.print_header();
		dove.print_result();
	}
	if (this->get_ConsoleOut() == true)
	{
		std::cout << "Dove Scheme: ";
		switch (dove.getIntegrationMethod())
		{
			case BE:
				std::cout << "Backwards-Euler method.";
				break;
				
			case FE:
				std::cout << "Forwards-Euler method.";
				break;
				
			case CN:
				std::cout << "Crank-Nicholson method. ";
				break;
				
			case BDF2:
				std::cout << "Backwards-Differentiation 2nd Order method.";
				break;
				
			case RK4:
				std::cout << "Runge-Kutta 4th Order method.";
				break;
				
			case RKF:
				std::cout << "Runge-Kutta-Fehlberg method.";
				break;
				
			default:
				std::cout << "Backwards-Euler method.";
				break;
				
		}
		std::cout << "\n------------------------------------------------------";
	}
	
	//Do-while loop
	this->estimate_parameters(dove);
	do
	{
		dove.update_timestep();
		if (this->get_ConsoleOut() == true)
		{
			std::cout << "\nSolving (" << dove.getNumFunc() << ") equation(s) at time (" << dove.getCurrentTime() << ") with time step (" << dove.getTimeStep() << "). Please wait...\n";
		}
		success = dove.solve_timestep();
		if (success != 0)
		{
			mError(simulation_fail);
			return -1;
		}
		this->store_variables(dove);
		if (this->get_FileOut() == true)
			dove.print_newresult();
		dove.update_states();
		this->estimate_parameters(dove);
	} while (dove.getEndTime() > (dove.getCurrentTime()+dove.getMinTimeStep()));
	if (this->get_ConsoleOut() == true)
		std::cout << "------------------------------------------------------\n\n";
	
	return success;
}

/*
 *	-------------------------------------------------------------------------------------
 *								End: Crane Class Definitions
 */

//Test function
int CRANE_TESTS()
{
	int success = 0;
	double time;
	
	Crane test;
	Dove dove;
	time = clock();
	
	FILE *file;
	file = fopen("output/CRANE_Tests.txt", "w+");
	if (file == nullptr)
	{
		system("mkdir output");
		file = fopen("output/CRANE_Tests.txt", "w+");
	}
	
	double W = 12.0; //kT
	double hb = 500.0*0.3048;// 500 ft
	double gz = 1155; //m (Nevada Test Site)
	int bins = 10;
	bool includeShear = true;
	bool isTight = true;
	
	std::cout << "\nTesting of the CRANE for Plumbbob Boltzman Bomb at Nevada Test Site\n";
	std::cout <<   "-------------------------------------------------------------------\n\n";
	std::cout << "Bomb Yield (kT) =\t" << W << std::endl;
	std::cout << "Burst Height (m) =\t" << hb << std::endl;
	std::cout << "Ground Altitude (m) =\t" << gz << std::endl;
	std::cout << "\n";
	
	test.establish_initial_conditions(W, gz, hb, bins, includeShear, isTight, dove);
	
	bool fileout = true;
	bool consoleout = true;
	double tol = 1e-8;
	double dtmin = 1e-8;
	double dtmax = 1.0;
	double dtmin_conv = 0.001;
	double t_out = 1.0;
	double endtime = 5.0;
	
	test.establish_dove_options(dove, file, fileout, consoleout, BE, CONSTANT, SGS, tol, dtmin, dtmax, dtmin_conv, t_out, endtime);
	
	bool isLinear = false;
	bool isPrecon = false;
	bool nl_out = true;
	bool l_out = false;
	int max_nlit = 50;
	int max_lit = 200;
	int restart = 20;
	int recursive = 2;
	double nl_abstol = 1e-8;
	double nl_reltol = 1e-8;
	double l_abstol = 1e-6;
	double l_reltol = 1e-6;
	
	test.establish_pjfnk_options(dove, QR, BT, isLinear, isPrecon, nl_out, l_out, max_nlit, max_lit, restart, recursive, nl_abstol, nl_reltol, l_abstol, l_reltol);
	
	//test.run_crane_simulation(dove);
	for (int i=0; i<8; i++)
	{
		std::cout << dove.getVariableName(i) << " =\t " << dove.getNewU(i, dove.getNewU()) << "\t: rate =\t " << dove.Eval_Func(i, dove.getNewU(), dove.getCurrentTime()) << std::endl;
	}
	
	test.run_crane_simulation(dove);
	
    dove.createJacobian();
    Matrix<double> J = dove.getNumJacobian();
    J.Display("J");
	
	std::cout << "Saturation Time (s) =\t" << test.get_saturation_time() << std::endl;
	
	test.get_part_conc_var().Display("n");
	
	time = clock() - time;
	std::cout << "\nTest Runtime: " << (time / CLOCKS_PER_SEC) << " seconds\n";
	
	if (file!= nullptr)
		fclose(file);
	
	return success;
}
