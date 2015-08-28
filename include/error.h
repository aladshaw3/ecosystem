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

//List of macro names for error type
#define generic_error 0
#define file_dne 1
#define indexing_error 2
#define magpie_reverse_error 3
#define simulation_fail 4
#define invalid_components 5
#define invalid_boolean 6
#define invalid_molefraction 7
#define invalid_gas_sum 8
#define invalid_solid_sum 9
#define scenario_fail 10
#define out_of_bounds 11
#define non_square_matrix 12
#define dim_mis_match 13
#define empty_matrix 14
#define opt_no_support 15
#define invalid_fraction 16
#define ortho_check_fail 17
#define unstable_matrix 18
#define no_diffusion 19
#define negative_mass 20
#define negative_time 21
#define matvec_mis_match 22
#define arg_matrix_same 23
#define singular_matrix 24
#define matrix_too_small 25
#define invalid_size 26
#define nullptr_func 27
#define invalid_norm 28
#define vector_out_of_bounds 29
#define zero_vector 30
#define tensor_out_of_bounds 31
#define non_real_edge 32
#define nullptr_error 33
#define invalid_atom 34
#define invalid_proton 35
#define invalid_neutron 36
#define invalid_electron 37
#define invalid_valence 38
#define string_parse_error 39
#define unregistered_name 40
#define rxn_rate_error 41
#define invalid_species 42
#define duplicate_variable 43
#define missing_information 44
#define invalid_type 45
#define key_not_found 46
#define anchor_alias_dne 47
#define initial_error 48
#define not_a_token 49
#define read_error 50
#define invalid_console_input 51


void error(int flag);

#endif /* ERROR_HPP_ */
