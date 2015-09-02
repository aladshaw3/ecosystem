//----------------------------------------
//  Created by Austin Ladshaw on 4/28/14
//  Copyright (c) 2014
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef ERROR_HPP_
#define ERROR_HPP_

#include <iostream>				  //Line to allow for read/write to the console using cpp functions

#ifndef mError
#define mError(i) \
{error(i);         \
std::cout << "Source: " << __FILE__ << "\nLine: " << __LINE__ << std::endl;}
#endif

//List of names for error type
typedef enum
{
generic_error,
file_dne,
indexing_error,
magpie_reverse_error,
simulation_fail,
invalid_components,
invalid_boolean,
invalid_molefraction,
invalid_gas_sum,
invalid_solid_sum,
scenario_fail,
out_of_bounds,
non_square_matrix,
dim_mis_match,
empty_matrix,
opt_no_support,
invalid_fraction,
ortho_check_fail,
unstable_matrix,
no_diffusion,
negative_mass,
negative_time,
matvec_mis_match,
arg_matrix_same,
singular_matrix,
matrix_too_small,
invalid_size,
nullptr_func,
invalid_norm,
vector_out_of_bounds,
zero_vector,
tensor_out_of_bounds,
non_real_edge,
nullptr_error,
invalid_atom,
invalid_proton,
invalid_neutron,
invalid_electron,
invalid_valence,
string_parse_error,
unregistered_name,
rxn_rate_error,
invalid_species,
duplicate_variable,
missing_information,
invalid_type,
key_not_found,
anchor_alias_dne,
initial_error,
not_a_token,
read_error,
invalid_console_input
} error_type;


void error(int flag);

#endif /* ERROR_HPP_ */
