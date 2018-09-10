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

#data.nuclides.keys()

#data.nuclides[92,238]

#data.nuclides[92,238][0].keys()

#data.nuclides[92,238][0]['decay modes']

#data.nuclides[92,238][0]['decay modes'].keys() #what are valid keys for decay modes?

# Pick the modes that correspond to keys
#data.nuclides[92,238][0]['decay modes']['A']['branch fraction']

#data.nuclides[92,238][0]['decay modes']['key for the mode']['branch fraction']


#Object nuclides is a map of nuclide data
A = 0				#atomic number
Z = 0				#mass number
hl = 0				#half-life in seconds
hl_s = ' '			#half-life string
decay_mode = ' '	#decay modes
branch_frac = 0		#branching fractions
n_modes = 0			#number of decay modes
symbol = ' '		#Element symbol

#Creat a list of nuclide keys to sort
key_list = []
for n in data.nuclides:
	key_list.append(n)

key_list.sort()

#Iterate through the sorted list
for n in key_list:
	A = n[0]
	Z = n[1]
	hl = data.nuclides[n][0]['half-life']
	hl_s = data.nuclides[n][0]['half-life string']
	symbol = data.nuclides[n][0]['symbol']
	
	#Correct some symbols from database
	if A == 0: symbol = 'n'
	if A == 112: symbol = 'Cn'
	if A == 113: symbol = 'Nh'
	if A == 114: symbol = 'Fl'
	if A == 115: symbol = 'Mc'
	if A == 116: symbol = 'Lv'
	if A == 117: symbol = 'Ts'
	if A == 118: symbol = 'Og'
	
	print 'Symbol = ' + str(symbol) + ' A = ' + str(A) + ' Z = ' + str(Z) + ' half-life = ' + str(hl_s)
	
	#Loop through the decay modes
	for m in data.nuclides[n][0]['decay modes']:
		n_modes = len(data.nuclides[n][0]['decay modes'])
		decay_mode = m
		branch_frac = data.nuclides[n][0]['decay modes'][m]['branch fraction']
		print '\t Num modes: ' + str(n_modes) + ' Mode: ' + str(m) + ' branch_frac: ' + str(branch_frac)

