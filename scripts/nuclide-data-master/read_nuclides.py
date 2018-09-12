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
n_modes = 0			#number of decay modes
symbol = ' '		#Element symbol
units = 'seconds'	#units for half-life
hl_in_units = 0		#half-life in the specific units
stable = False		#stability of the nuclide
AW = 0				#atomic weight of the nuclide
iso = False			#True if nuclide is an excited state
iso_name = ' '		#Name of isotope for specific-isotope decay

#Open yaml file to write to
file = open('../../database/NuclideLibrary.yml', 'w')

#Creat a list of nuclide keys to sort
key_list = []
for n in data.nuclides:
	key_list.append(n)

key_list.sort()

#Iterate through the sorted list
for n in key_list:
	Z = n[0] #atomic num
	A = n[1] #mass num
	hl = data.nuclides[n][0]['half-life']
	stable = data.nuclides[n][0]['stable']
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
	file.write(str(symbol) + '-' + str(A) + ':\n')
	file.write('---\n')
	file.write('symbol: ' + str(symbol) + '\n')
	file.write('atom_num: ' + str(Z) + '\n')
	file.write('mass_num: ' + str(A) + '\n')
	file.write('atom_weight: ' + str(AW) + '\n')
	file.write('isomeric: ' + str(iso) + '\n')
	file.write('half_life: ' + str(hl) + '\n')
	file.write('hl_units: ' + str(units) + '\n')
	
	#Loop through the decay modes
	i = 0
	for m in data.nuclides[n][0]['decay modes']:
		n_modes = len(data.nuclides[n][0]['decay modes'])
		decay_mode = m
		
		#Rename and group some decay modes
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
		
		#Make additional corrections to data
		if m == None: decay_mode = 'stable'
		branch_frac = data.nuclides[n][0]['decay modes'][m]['branch fraction']
		if decay_mode == 'stable': branch_frac = 0
		if decay_mode == 'stable': stable = True
		if branch_frac == None: branch_frac = 0
		
		#write out stability condition (only on first iteration)
		if i == 0: file.write('stable: ' + str(stable) + '\n')
		
		#Check for specific atom emission
		if decay_mode != 'stable' and decay_mode != 'alpha' and decay_mode != 'spontaneous-fission' and decay_mode != 'beta+' and decay_mode != 'beta-' and decay_mode != 'isomeric-transition' and decay_mode != 'neutron-emission' and decay_mode != 'beta-/neutron-emission' and decay_mode != 'beta+/proton-emission' and decay_mode != 'proton-emission' and decay_mode != 'beta+/alpha' and decay_mode != 'beta+/beta+' and decay_mode != 'beta-/beta-' and decay_mode != 'beta-/neutron-emission/neutron-emission' and decay_mode != 'beta-/alpha' and decay_mode != 'proton-emission/proton-emission' and decay_mode != 'neutron-emission/neutron-emission' and decay_mode != 'beta-/neutron-emission/neutron-emission/neutron-emission' and decay_mode != 'beta-/neutron-emission/neutron-emission/neutron-emission/neutron-emission' and decay_mode != 'beta+/proton-emission/proton-emission' and decay_mode != 'beta+/proton-emission/proton-emission/proton-emission':
			
			decay_mode = 'specific-isotope'
			
			#Set the isotope name for the specific isotope emission during decay
			try:
				tempA = data.Nuclide(m).A
				tempEle = data.Nuclide(m).element
				iso_name = str(tempEle) + '-' + str(tempA)
			except:
				tempA = int(data.weight(m))
				tempEle = m
				iso_name = str(tempEle) + '-' + str(tempA)
		
		#End if statement
		
		#write out remaining decay info
		if i == 0: file.write('\n- decay_modes:\n')
		file.write('  - mode' + str(i) + ':\n')
		file.write('    type: ' + str(decay_mode) + '\n')
		file.write('    branch_frac: ' + str(branch_frac) + '\n')
		if decay_mode == 'specific-isotope': file.write('    isotope: ' + str(iso_name) + '\n')
		file.write('\n')
		
		i = i + 1

	#END decay modes loop

	file.write('...\n\n')
#END key_list of nuclides loop

print '\nYaml Library Construction Complete!\n'

#Close the yaml file
file.close()
