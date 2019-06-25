# Python script writted by Austin Ladshaw for the purpose of developing a
# simple nuclide data file containing half-lifes and decay modes of nuclides
#
# Austin Ladshaw does not claim copyright to any of the nuclide_data.py code
# or methods. This script was simply developed to take advantage of existing
# works in python for the express purposes of developing a more convenient
# nuclide library for ECOSYSTEM.
#
# nuclide_data.py and related files were developed by GitHub user jhykes
# https://github.com/jhykes/nuclide-data
#
# Nuclide data comes from 2 sources:
# ----------------------------------
#	1. Jagdish K. Tuli, Nuclear Wallet Cards, National Nuclear Data Center,
#		April 2005. http://www.nndc.bnl.gov/wallet
#	2. J. S. Coursey, D. J. Schwab, J. J. Tsai, and R. A. Dragoset, Atomic Weights
#		and Isotopic Compositions with Relative Atomic Masses, NIST Physical Measurement
#		Laboratory, updated July 2010. http://www.nist.gov/pml/data/comp.cfm
#



# Library of nuclides from n-ENDF-B-VII.1.endf.list
import nuclide_data as data

#Object nuclides is a map of nuclide data
A = 0				#mass number
Z = 0				#atomic number
hl = 0				#half-life in seconds
hl_s = ' '			#half-life string
decay_mode = ' '	#decay modes
branch_frac = 0		#branching fractions
ex_mass = 0       	#Mass excess energy in MeV
excite_energy = []	#List of excitation energies for the nuclide
n_modes = 0			#number of decay modes
symbol = ' '		#Element symbol
units = 'seconds'	#units for half-life
hl_in_units = 0		#half-life in the specific units
stable = False		#stability of the nuclide
AW = 0				#atomic weight of the nuclide
iso = False			#True if nuclide is an excited state
iso_name = ' '		#Name of isotope for specific-isotope decay
daughter = ' '		#Name of the daughter isotope produced from decay mode (check to make sure daughter is real)
particle_em = ' '	#Name of the particle emitted (if any) to produce the daughter
part_num = 0		#Number of those particles emitted

react_sum = 0.0		#add all protons and neutrons for the isotope
prod_sum = 0.0		#add all protons and neutrons for each daughter and emitted particle (modulated by branch fraction)
combo_sum = 0.0		#Sum of combo decay modes
d_A = 0				#mass number of the daughter
ep_A = 0			#mass number of the emitted particle

modes = {}			#create a map object where keys are decay modes and values are branch fractions
q_val = {}			#create a map object where keys are decay modes and values are Q-values
daughters = {}		#create a map object where keys are decay modes and values are daughter particles
emissions = {}		#create a map object where keys are decay modes and values are emitted particles
number_parts = {}	#create a map object where keys are decay modes and values are number of emitted particles
daughter_A = {}
emission_A = {}

#Open yaml file to write to
file = open('../../database/NuclideLibrary.yml', 'w')

#Creat a list of nuclide keys to sort
key_list = []
for n in data.nuclides:
	key_list.append(n)

key_list.sort()

#Iterate through the sorted list
problems = 0
for n in key_list:
	Z = n[0] #atomic num
	A = n[1] #mass num
	react_sum = float(A)
	prod_sum = 0
	combo_sum = 0.0
	hl = data.nuclides[n][0]['half-life']
	stable = data.nuclides[n][0]['stable']
	ex_mass = str(data.nuclides[n][0]['mass excess']).split('+')[0]
        excite_energy = data.isomers(Z,A)
	try:
		AW = data.weight(Z,A,0)
	except:
		AW = A
	iso = data.nuclides[n][0]['isomeric']

	#Manual correction for blanks
	if not AW: AW = A
	
	#Specify unit basis for half-life
	if hl <= 60: units = 'seconds'
	if hl > 60: units = 'minutes'
	if hl > 3600: units = 'hours'
	if hl > 86400: units = 'days'
	if hl > 31557600: units = 'years'
	
	#Perform unit conversion
	if units == 'seconds': hl = hl
	if units == 'minutes': hl = hl/60
	if units == 'hours': hl = hl/3600
	if units == 'days': hl = hl/86400
	if units == 'years': hl = hl/31557600
	
	#Make sure a value is given
	if hl == 0:
		if stable == False:
			hl = 1E-20
		else:
			hl = inf
	
	hl_s = data.nuclides[n][0]['half-life string']
	symbol = data.nuclides[n][0]['symbol']
	
	#Correct some symbols from database
	if Z == 0: symbol = 'n'
	if Z == 112: symbol = 'Cn'
	if Z == 113: symbol = 'Nh'
	if Z == 114: symbol = 'Fl'
	if Z == 115: symbol = 'Mc'
	if Z == 116: symbol = 'Lv'
	if Z == 117: symbol = 'Ts'
	if Z == 118: symbol = 'Og'
	
	#write out basic info
	#	file.write(str(symbol) + '-' + str(A) + ':\n')
	#file.write('---\n')
	#file.write('symbol: ' + str(symbol) + '\n')
	#file.write('atom_num: ' + str(Z) + '\n')
	#file.write('mass_num: ' + str(A) + '\n')
	#file.write('atom_weight: ' + str(AW) + '\n')
	#file.write('isomeric: ' + str(iso) + '\n')
	#file.write('half_life: ' + str(hl) + '\n')
	#file.write('hl_units: ' + str(units) + '\n')
	
	#Loop through the decay modes
	i = 0
	combo_sum = 0.0
	for m in data.nuclides[n][0]['decay modes']:
		n_modes = len(data.nuclides[n][0]['decay modes'])
		decay_mode = m
		
		#Rename and group some decay modes
		if m == None: decay_mode = 'stable'
		if stable == True: decay_mode = 'stable'
		if m == 'A' or m == 'A<': decay_mode = 'alpha'
		if m == 'SF' or m == 'EF' or m == 'BF': decay_mode = 'spontaneous-fission'
		if m == 'EC': decay_mode = 'beta+'
		if m == 'B-': decay_mode = 'beta-'
		if m == 'IT': decay_mode = 'isomeric-transition'
		if m == 'N': decay_mode = 'neutron-emission'
		if m == 'BN' or m == 'BNA': decay_mode = 'beta-/neutron-emission'
		if m == 'EP': decay_mode = 'beta+/proton-emission'
		if m == 'P' or m == '2A': decay_mode = 'proton-emission'
		if m == 'EA': decay_mode = 'beta+/alpha'
		if m == '2EC': decay_mode = 'beta+/beta+'
		if m == '2B-': decay_mode = 'beta-/beta-'
		if m == 'B2N': decay_mode = 'beta-/neutron-emission/neutron-emission'
		if m == 'BA' or m == 'B3A': decay_mode = 'beta-/alpha'
		if m == '2P': decay_mode = 'proton-emission/proton-emission'
		if m == '2N' or m == '2N?': decay_mode = 'neutron-emission/neutron-emission'
		if m == 'B3N': decay_mode = 'beta-/neutron-emission/neutron-emission/neutron-emission'
		if m == 'B4N': decay_mode = 'beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission'
		if m == 'E2P': decay_mode = 'beta+/proton-emission/proton-emission'
		if m == 'E3P': decay_mode = 'beta+/proton-emission/proton-emission/proton-emission'
		
		#Determine daughters and particle_emissions
		if decay_mode == 'stable':
			daughter = 'None'
			particle_em = 'None'
			part_num = 0
			d_A = 0
			ep_A = 0
		
		if decay_mode == 'alpha':
			#daughter
			Zn = Z - 2
			An = A - 4
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'He-4'
				part_num = 1
				d_A = An
				ep_A = 4
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0
					
		if decay_mode == 'spontaneous-fission':
			#daughter
			Zn = int(Z/2)
			An = int(A/2)
			Zm = Z - Zn
			Am = A - An
			try:
				sym1 = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym1 = 'n'
				if Zn == 112: sym1 = 'Cn'
				if Zn == 113: sym1 = 'Nh'
				if Zn == 114: sym1 = 'Fl'
				if Zn == 115: sym1 = 'Mc'
				if Zn == 116: sym1 = 'Lv'
				if Zn == 117: sym1 = 'Ts'
				if Zn == 118: sym1 = 'Og'
				sym2 = data.nuclides[(Zm,Am)][0]['symbol']
				#Correct some symbols from database
				if Zm == 0: sym2 = 'n'
				if Zm == 112: sym2 = 'Cn'
				if Zm == 113: sym2 = 'Nh'
				if Zm == 114: sym2 = 'Fl'
				if Zm == 115: sym2 = 'Mc'
				if Zm == 116: sym2 = 'Lv'
				if Zm == 117: sym2 = 'Ts'
				if Zm == 118: sym2 = 'Og'
				daughter = str(sym1) + '-' + str(An)
				particle_em = str(sym2) + '-' + str(Am)
				part_num = 1
				d_A = An
				ep_A = Am
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta+':
			#daughter
			Zn = Z - 1
			An = A
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'None'
				part_num = 0
				d_A = An
				ep_A = 0
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-':
			#daughter
			Zn = Z + 1
			An = A
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'None'
				part_num = 0
				d_A = An
				ep_A = 0
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'isomeric-transition':
			#daughter
			Zn = Z
			An = A
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'None'
				part_num = 0
				d_A = An
				ep_A = 0
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'neutron-emission':
			#daughter
			Zn = Z
			An = A - 1
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'n-1'
				part_num = 1
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-/neutron-emission':
			#daughter
			Zn = Z + 1
			An = A - 1
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'n-1'
				part_num = 1
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta+/proton-emission':
			#daughter
			Zn = Z - 2
			An = A - 1
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'H-1'
				part_num = 1
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'proton-emission':
			#daughter
			Zn = Z - 1
			An = A - 1
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'H-1'
				part_num = 1
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta+/alpha':
			#daughter
			Zn = Z - 3
			An = A - 4
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'He-4'
				part_num = 1
				d_A = An
				ep_A = 4
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta+/beta+':
			#daughter
			Zn = Z - 2
			An = A
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'None'
				part_num = 0
				d_A = An
				ep_A = 0
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-/beta-':
			#daughter
			Zn = Z + 2
			An = A
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'None'
				part_num = 0
				d_A = An
				ep_A = 0
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-/neutron-emission/neutron-emission':
			#daughter
			Zn = Z + 1
			An = A - 2
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'n-1'
				part_num = 2
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-/alpha':
			#daughter
			Zn = Z + 1 - 2
			An = A - 4
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'He-4'
				part_num = 1
				d_A = An
				ep_A = 4
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'proton-emission/proton-emission':
			#daughter
			Zn = Z - 2
			An = A - 2
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'H-1'
				part_num = 2
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'neutron-emission/neutron-emission':
			#daughter
			Zn = Z
			An = A - 2
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'n-1'
				part_num = 2
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-/neutron-emission/neutron-emission/neutron-emission':
			#daughter
			Zn = Z + 1
			An = A - 3
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'n-1'
				part_num = 3
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission':
			#daughter
			Zn = Z + 1
			An = A - 4
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'n-1'
				part_num = 4
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta+/proton-emission/proton-emission':
			#daughter
			Zn = Z - 1 - 2
			An = A - 2
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'H-1'
				part_num = 2
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0

		if decay_mode == 'beta+/proton-emission/proton-emission/proton-emission':
			#daughter
			Zn = Z - 1 - 3
			An = A - 3
			try:
				sym = data.nuclides[(Zn,An)][0]['symbol']
				#Correct some symbols from database
				if Zn == 0: sym = 'n'
				if Zn == 112: sym = 'Cn'
				if Zn == 113: sym = 'Nh'
				if Zn == 114: sym = 'Fl'
				if Zn == 115: sym = 'Mc'
				if Zn == 116: sym = 'Lv'
				if Zn == 117: sym = 'Ts'
				if Zn == 118: sym = 'Og'
				daughter = str(sym) + '-' + str(An)
				particle_em = 'H-1'
				part_num = 3
				d_A = An
				ep_A = 1
			except:
				decay_mode = 'undefined'
				daughter = 'None'
				particle_em = 'None'
				part_num = 0
				d_A = 0
				ep_A = 0
		
		#End if statements
	
		#Make additional corrections to data
		branch_frac = data.nuclides[n][0]['decay modes'][m]['branch fraction']
		if decay_mode == 'stable' or decay_mode == 'undefined': branch_frac = 0
		if decay_mode == 'stable': stable = True
		if branch_frac == None: branch_frac = 0.01
		if decay_mode == 'isomeric-transition': branch_frac = 0

		#Combo sum
		if decay_mode == 'beta-/neutron-emission' or decay_mode == 'beta+/proton-emission' or decay_mode == 'beta+/alpha' or decay_mode == 'beta-/neutron-emission/neutron-emission' or decay_mode == 'beta-/alpha' or decay_mode == 'beta-/neutron-emission/neutron-emission/neutron-emission' or decay_mode == 'beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission' or decay_mode == 'beta+/proton-emission/proton-emission' or decay_mode == 'beta+/proton-emission/proton-emission/proton-emission':
			
			combo_sum = combo_sum + branch_frac
		
		#write out stability condition (only on first iteration)
		#if i == 0: file.write('stable: ' + str(stable) + '\n')
		
		#Check for specific atom emission
		if decay_mode != 'stable' and decay_mode != 'alpha' and decay_mode != 'spontaneous-fission' and decay_mode != 'beta+' and decay_mode != 'beta-' and decay_mode != 'isomeric-transition' and decay_mode != 'neutron-emission' and decay_mode != 'beta-/neutron-emission' and decay_mode != 'beta+/proton-emission' and decay_mode != 'proton-emission' and decay_mode != 'beta+/alpha' and decay_mode != 'beta+/beta+' and decay_mode != 'beta-/beta-' and decay_mode != 'beta-/neutron-emission/neutron-emission' and decay_mode != 'beta-/alpha' and decay_mode != 'proton-emission/proton-emission' and decay_mode != 'neutron-emission/neutron-emission' and decay_mode != 'beta-/neutron-emission/neutron-emission/neutron-emission' and decay_mode != 'beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission' and decay_mode != 'beta+/proton-emission/proton-emission' and decay_mode != 'beta+/proton-emission/proton-emission/proton-emission' and decay_mode != 'undefined':
			
			decay_mode = 'specific-isotope'
			
			#Set the isotope name for the specific isotope emission during decay
			#NOTE: particle emitted is the specific isotope, while daughter is the result of the emission
			try:
				tempA = data.Nuclide(m).A
				tempEle = data.Nuclide(m).element
				iso_name = str(tempEle) + '-' + str(tempA)
				#daughter
				Zn = Z - data.Nuclide(m).Z
				An = A - tempA
				try:
					sym = data.nuclides[(Zn,An)][0]['symbol']
					#Correct some symbols from database
					if Zn == 0: sym = 'n'
					if Zn == 112: sym = 'Cn'
					if Zn == 113: sym = 'Nh'
					if Zn == 114: sym = 'Fl'
					if Zn == 115: sym = 'Mc'
					if Zn == 116: sym = 'Lv'
					if Zn == 117: sym = 'Ts'
					if Zn == 118: sym = 'Og'
					daughter = str(sym) + '-' + str(An)
					particle_em = iso_name
					part_num = 1
					d_A = An
					ep_A = tempA
				except:
					decay_mode = 'undefined'
					daughter = 'None'
					particle_em = 'None'
					part_num = 0
					d_A = 0
					ep_A = 0
			except:
				tempA = int(data.weight(m))
				tempEle = m
				iso_name = str(tempEle) + '-' + str(tempA)
				try:
					Zn = Z - data.Nuclide(iso_name).Z
					An = A - tempA
					try:
						sym = data.nuclides[(Zn,An)][0]['symbol']
						#Correct some symbols from database
						if Zn == 0: sym = 'n'
						if Zn == 112: sym = 'Cn'
						if Zn == 113: sym = 'Nh'
						if Zn == 114: sym = 'Fl'
						if Zn == 115: sym = 'Mc'
						if Zn == 116: sym = 'Lv'
						if Zn == 117: sym = 'Ts'
						if Zn == 118: sym = 'Og'
						daughter = str(sym) + '-' + str(An)
						particle_em = iso_name
						part_num = 1
						d_A = An
						ep_A = tempA
					except:
						decay_mode = 'undefined'
						daughter = 'None'
						particle_em = 'None'
						part_num = 0
						d_A = 0
						ep_A = 0
				except:
					decay_mode = 'undefined'
					daughter = 'None'
					particle_em = 'None'
					part_num = 0
					d_A = 0
					ep_A = 0
	
		#End if statement
		
		#Set the modes and branch fractions
		modes[decay_mode] = branch_frac
		q_val[decay_mode] = data.nuclides[n][0]['decay modes'][m]['Q-value']
		daughters[decay_mode] = daughter
		emissions[decay_mode] = particle_em
		number_parts[decay_mode] = part_num
		daughter_A[decay_mode] = d_A
		emission_A[decay_mode] = ep_A
		
		#Check the summation of particles
		if daughter != 'None':
			prod_sum = float(prod_sum + branch_frac*d_A + branch_frac*part_num*ep_A)
		#End if
		
		i = i + 1

	#END decay modes loop

	#Check summations outside of the loop, then may have to go back and fix
	SumError = False
	if daughter != 'None':
		if abs(react_sum - prod_sum) > 1e-10:
			SumError = True
	#End if
	
	i = 0
	#Start new decay loop for editing
	prod_sum = 0.0
	branch_sum = 0.0;
	all_modes_bad = True
	for m in modes:
		
		if SumError == True:
			if m == 'beta-' or m == 'beta+':
				#manual corrections
				if Z == 3 and A == 11:
					modes[m] = 0.95873
				elif Z == 51 and A == 104:
					modes[m] = 0.92
				elif Z == 52 and A == 108:
					modes[m] = 0.486
				elif Z == 26 and A == 45:
					modes[m] = 0.129
				elif Z == 52 and A == 109:
					modes[m] = 0.8670
				elif Z == 53 and A == 110:
					modes[m] = 0.709
				elif Z == 54 and A == 113:
					modes[m] = 0.92983
				elif Z == 54 and A == 115:
					modes[m] = 0.996597
				elif Z == 55 and A == 114:
					modes[m] = 0.9109
				elif Z == 56 and A == 114:
					modes[m] = 0.790966
				elif Z == 61 and A == 128:
					modes[m] = 1.0
				elif Z == 80 and A == 179:
					modes[m] = 0.4485
				elif Z == 80 and A == 181:
					modes[m] = 0.72989991
				elif Z == 80 and A == 183:
					modes[m] = 0.8829974
				#automatic corrections (works for 90% of errors with combination decay)
				else:
					modes[m] = 1.0 - combo_sum
			if m == 'proton-emission' or m == 'neutron-emission':
				#manual corrections
				if Z == 39 and A == 77:
					modes[m] = 0
				elif Z == 47 and A == 93:
					modes[m] = 0
			if m == 'proton-emission/proton-emission':
				#manual corrections
				if Z == 26 and A == 45:
					modes[m] = 0.57
			if m == 'alpha':
				#manual corrections
				if Z == 52 and A == 109:
					modes[m] = 0.03895
				elif Z == 61 and A == 128:
					modes[m] = 0
			if m == 'beta+/proton-emission':
				if Z == 54 and A == 110:
					modes[m] = 0.36
				elif Z == 61 and A == 128:
					modes[m] = 0.0
	
		branch_sum = branch_sum + modes[m]
		
		#Check the summation of particles
		if m != 'None':
			prod_sum = float(prod_sum + modes[m]*daughter_A[m] + modes[m]*number_parts[m]*emission_A[m])

		if all_modes_bad == True:
			if m != 'undefined' and m!= 'stable' and m != 'isomeric-transition': all_modes_bad = False
			else: all_modes_bad = True
	
		i = i + 1
	#End new decay loop
	
	#redefine stability
	if all_modes_bad == True: stable = 'True'
	
	#Print out basic info
	file.write(str(symbol) + '-' + str(A) + ':\n')
	file.write('---\n')
	file.write('symbol: ' + str(symbol) + '\n')
	file.write('atom_num: ' + str(Z) + '\n')
	file.write('mass_num: ' + str(A) + '\n')
	file.write('atom_weight: ' + str(AW) + '\n')
	file.write('isomeric: ' + str(iso) + '\n')
        file.write('mass-excess: ' + str(ex_mass) + '\n')
	file.write('half_life: ' + str(hl) + '\n')
	file.write('hl_units: ' + str(units) + '\n')
        file.write('\n- excitations:\n')
        i = 0
        for en in excite_energy:
		file.write('  level0'+str(i)+': ' + str(en) + '\n')
		i += 1
	
	i = 0
	#Start new decay loop for final editing
	prod_sum = 0.0
	for m in modes:
		
		#Fix modes
		if branch_sum > 0.0:
			modes[m] = modes[m]/branch_sum
		else: modes[m] = 1.0/len(modes)

		if Z == 75 and A == 193: modes[m] = 1.0
	
		#Check the summation of particles
		if m != 'None':
			prod_sum = float(prod_sum + modes[m]*daughter_A[m] + modes[m]*number_parts[m]*emission_A[m])
		
		#write out remaining decay info
		#write out stability condition (only on first iteration)
		if i == 0: file.write('\nstable: ' + str(stable) + '\n')
		if i == 0: file.write('\n- decay_modes:\n')
		file.write('  - mode' + str(i) + ':\n')
		file.write('    type: ' + str(m) + '\n')
		file.write('    branch_frac: ' + str(modes[m]) + '\n')
		file.write('    Q-value: ' + str(q_val[m]) + '\n')
		file.write('    daughter: ' + str(daughters[m]) + '\n')
		file.write('    part_emitted: ' + str(emissions[m]) + '\n')
		file.write('    num_parts: ' + str(number_parts[m]) + '\n')
		file.write('\n')
		i = i + 1
	#End new decay loop

	#Check summations outside of the loop, then may have to go back and fix
	SumError = False
	if daughter != 'None':
		if abs(react_sum - prod_sum) > 1e-10:
			SumError = True
			problems = problems + 1
			print str(symbol) + '-' + str(A) + ': ' + str(react_sum) + '-' + str(prod_sum)
			print str(combo_sum)
			print modes
	#End if

	modes = {}
	daughters = {}
	emissions = {}
	number_parts = {}
	daughter_A = {}
	emission_A = {}
	file.write('...\n\n')
#END key_list of nuclides loop

print '\nYaml Library Construction Complete!\n'
print 'Number problems = ' + str(problems)

#Close the yaml file
file.close()
