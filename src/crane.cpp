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
	energy_frac = 0.0;
	eccentricity = 0.75;
	air_density = 1.225;
	air_viscosity = 1.81e-5;
	slip_factor = 1.0;
	davies_num = 1.0;
	vapor_pressure = 0.0;
	sat_vapor_pressure = 1.0;
	includeShearVel = false;
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
}

//Default destructor
Crane::~Crane()
{
	amb_temp.clear();
	atm_press.clear();
	rel_humid.clear();
	wind_vel.clear();
	part_hist.clear();
	settling_rate.clear();
	part_conc.clear();
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

void Crane::set_includeShearVel(bool val)
{
	this->includeShearVel = val;
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

bool Crane::get_includeShearVel()
{
	return this->includeShearVel;
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
	this->set_sat_vapor_pressure( 611.0*pow(T/273.0, -5.13)*exp(25.0*(T-273.0)/T) );
}

void Crane::compute_xe(double Te, double P, double HR)
{
	double val = (109.98*HR/29.0/P)*pow(Te/273.0, -5.13)*exp(25.0*(Te-273.0)/Te);
	this->set_xe(val);
}

void Crane::compute_air_density(double P, double Pws, double HR, double T)
{
	//Assume P and Pws come in as Pa --> convert to mBar
	P = P*0.01;
	Pws = Pws*0.01;
	double val = ( P - (Pws*HR*(1.0-this->eps)/100.0) ) / (2.8679*T);
	this->set_air_density(val);
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

// Below are listed return functions specific for temperature integral related values

// Below are listed return functions specific for air profile related values

// Below are listed return functions specific for particle histograms

/*
 *	-------------------------------------------------------------------------------------
 *								End: Crane Class Definitions
 */

//Test function
int CRANE_TESTS()
{
	int success = 0;
	
	std::cout << "here is the test\n";
	
	return success;
}
