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
	energy_frac = 0.5;
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
		sum += this->settling_rate[it->first]*M_PI*dj*dj*dj*n(i,0)/6.0;
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
	this->set_min_dia(min);
	this->set_max_dia(max);
	this->set_num_bins(size);
	this->set_mean_dia(avg);
	this->set_std_dia(std);
	
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
	this->amb_temp.clear();
	this->atm_press.clear();
	this->rel_humid.clear();
	this->wind_vel.clear();
	
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

/*
 *	-------------------------------------------------------------------------------------
 *								End: Crane Class Definitions
 */

//Test function
int CRANE_TESTS()
{
	int success = 0;
	
	Crane test;
	
	std::cout << "Test of particle histogram\n\n";
	
	test.create_part_hist(0.0001, 100, 10, 0.15, 2);
	test.display_part_hist();
	
	std::cout << "\nTest of air viscosity and air density\n\n";
	std::cout << "Temp(oC)\tVis.(cP)\tDens.(kg/m^3)\tPws(Pa)\tPv(Pa)\tCpw(J/kg/K)\tCpa(J/kg/K)\tCs(J/kg/K)\n";
	for (int i=0; i<=100 ; i++)
	{
		double temp = 273 + ((double)i*10);
		test.compute_air_viscosity(temp);
		
		test.compute_xe(temp, 101325, 50.0);
		double x = test.get_xe();
		
		//NOTE: Pressure must be given in Pa (Pressure (Pa), x (kg/kg), temp (K))
		test.compute_air_density(101325, x, temp);
		
		test.compute_sat_vapor_pressure(temp);
		
		test.compute_vapor_pressure(101325, x);
		
		test.compute_spec_heat_water(temp);
		
		test.compute_spec_heat_entrain(temp);
		
		test.compute_spec_heat_conds(temp);
		
		std::cout << temp-273.0 << "\t" << test.get_air_viscosity()*1000.0 << "\t" << test.get_air_density() << "\t" << test.get_sat_vapor_pressure() << "\t" << test.get_vapor_pressure() << "\t" << test.get_spec_heat_water() << "\t" << test.get_spec_heat_entrain() << "\t" << test.get_spec_heat_conds() << std::endl;
	}
	std::cout << std::endl;
	
	test.create_default_atmosphere();
	
	
	std::cout << "\nTest of air profiler\n\n";
	std::cout << "Alt(m)\tT(K)\tP(Pa)\tHR(%)\tv(m/s)\t\n";
	for (int i=0; i<=80; i++)
	{
		double alt = (double)i*1000.0;
		std::cout << alt << "\t" << test.return_amb_temp(alt) << "\t" << test.return_atm_press(alt) << "\t" << test.return_rel_humid(alt) << "\t" << test.return_wind_speed(alt) << std::endl;
	}
	
	std::cout << "\nTest of integrals\n\n";
	
	double T = 900.0;
	double Te = 200.0;
	test.compute_spec_heat_entrain_integral(T, Te);
	std::cout << "Integral of c_pa(T) from T = " << T << " to Te = " << Te << " ==>\t" << test.get_spec_heat_entrain_integral() << std::endl;
	std::cout << "Integral of c_s(T) from T = " << T << " to Te = " << Te << " ==>\t" << test.return_spec_heat_conds_integral(T,Te) << std::endl;
	
	return success;
}
