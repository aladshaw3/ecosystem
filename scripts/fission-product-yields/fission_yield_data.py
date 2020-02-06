#
#Python script to read ENDF-6 formatted Fission Yield Data.
#
#Data from National Nuclear Data Center
# https://www.nndc.bnl.gov/endf/b7.1/index.html
#
#Download Library --> zipFile for Neutron Induced Fission Product Yields
#
#Format of data file explained clearly in JENDL FP Decay Data File 2011 and
#	Fission Yields Data File 2011 (https://www.oecd-nea.org/science/wpec/sg37/Meetings/2013_May/13_JENDL.pdf)
#
#	See Page 51 of above reference for explaination of formatting.

import os.path
import warnings
import string
import numpy as np

# ----------- Class Object to Hold ENDF Data for single energy level --------------
class FPY_Energy:
	Energy = 0		#Energy of neutron source (eV)
	Yield = {}		#Map of ZAIDs --> Yield Values
	DY = {}			#Map of ZAIDs --> Yield Uncertainties

# --------------------- End Class ------------------------------

# ----------- Class Object to Hold ENDF Data for single isotope --------------
class FPY_Iso:
	ZAID = ' '		#ZAID of the fissioning isotope
	NIFP = True		#True = Neutron Induced, False = Spontaneous
	Levels = {}		#List of all FPY_Energy levels and associated data

# --------------------- End Class ------------------------------

# --------------------- Utility Functions ------------------------------
def parse_num(s_flt):
	ZA = {}
	if int(s_flt[10]) < 5:
		ZA[1] = int(s_flt[1]+s_flt[3])
		ZA[2] = int(s_flt[4]+s_flt[5]+s_flt[6])
	else:
		ZA[1] = int(s_flt[1]+s_flt[3]+s_flt[4])
		ZA[2] = int(s_flt[5]+s_flt[6]+s_flt[7])
	return ZA

def sci_notation(s_flt):
	num = ''
	for char in s_flt:
		if char == '-' or char == '+':
			num+='E'
		num+=char
	return num
# --------------------- End Utility Func ------------------------------

# ----------- Function to read ENDF-6 format for FPY --------------
def read_fpy(filepath):
	#Define all z code
	z = {}
	z[1] = 'H'
	z[2] = 'He'
	z[3] = 'Li'
	z[4] = 'Be'
	z[5] = 'B'
	z[6] = 'C'
	z[7] = 'N'
	z[8] = 'O'
	z[9] = 'F'
	z[10] = 'Ne'
	z[11] = 'Na'
	z[12] = 'Mg'
	z[13] = 'Al'
	z[14] = 'Si'
	z[15] = 'P'
	z[16] = 'S'
	z[17] = 'Cl'
	z[18] = 'Ar'
	z[19] = 'K'
	z[20] = 'Ca'
	z[21] = 'Sc'
	z[22] = 'Ti'
	z[23] = 'V'
	z[24] = 'Cr'
	z[25] = 'Mn'
	z[26] = 'Fe'
	z[27] = 'Co'
	z[28] = 'Ni'
	z[29] = 'Cu'
	z[30] = 'Zn'
	z[31] = 'Ga'
	z[32] = 'Ge'
	z[33] = 'As'
	z[34] = 'Se'
	z[35] = 'Br'
	z[36] = 'Kr'
	z[37] = 'Rb'
	z[38] = 'Sr'
	z[39] = 'Y'
	z[40] = 'Zr'
	z[41] = 'Nb'
	z[42] = 'Mo'
	z[43] = 'Tc'
	z[44] = 'Ru'
	z[45] = 'Rh'
	z[46] = 'Pd'
	z[47] = 'Ag'
	z[48] = 'Cd'
	z[49] = 'In'
	z[50] = 'Sn'
	z[51] = 'Sb'
	z[52] = 'Te'
	z[53] = 'I'
	z[54] = 'Xe'
	z[55] = 'Cs'
	z[56] = 'Ba'
	z[57] = 'La'
	z[58] = 'Ce'
	z[59] = 'Pr'
	z[60] = 'Nd'
	z[61] = 'Pm'
	z[62] = 'Sm'
	z[63] = 'Eu'
	z[64] = 'Gd'
	z[65] = 'Tb'
	z[66] = 'Dy'
	z[67] = 'Ho'
	z[68] = 'Er'
	z[69] = 'Tm'
	z[70] = 'Yb'
	z[71] = 'Lu'
	z[72] = 'Hf'
	z[73] = 'Ta'
	z[74] = 'W'
	z[75] = 'Re'
	z[76] = 'Os'
	z[77] = 'Ir'
	z[78] = 'Pt'
	z[79] = 'Au'
	z[80] = 'Hg'
	z[81] = 'Tl'
	z[82] = 'Pb'
	z[83] = 'Bi'
	z[84] = 'Po'
	z[85] = 'At'
	z[86] = 'Rn'
	z[87] = 'Fr'
	z[88] = 'Ra'
	z[89] = 'Ac'
	z[90] = 'Th'
	z[91] = 'Pa'
	z[92] = 'U'
	z[93] = 'Np'
	z[94] = 'Pu'
	z[95] = 'Am'
	z[96] = 'Cm'
	z[97] = 'Bk'
	z[98] = 'Cf'
	z[99] = 'Es'
	z[100] = 'Fm'

	#Initialize
	FPs = FPY_Iso()
	FPs.Levels = {}
	start = False
	c = {}
	NN = 0		#Total number of items to read for given energy level (=4*NFP)
	NFP = 0		#Number of states to read (i.e., state = 4 items)
	LE = -1		#Number of energy levels (if =0, then only 1 energy)
	cur = 0		#Current energy level
	iNN = 0		#index for NN
	NewLevel = True   #Changes to false after reading first level, then changes back to true
	odd = True	#Boolean to determine if we are odd or even line
	id = ''		#id value carried over between lines
	includeState = False		#True to differentiate between states, False to combine states
	name = ''	#FP name

	if filepath[0] == "n":
		FPs.NIFP = True
	else:
		FPs.NIFP = False

	#Loop over lines in file
	for line in open(filepath):
		#Skip over all lines before initial data
		c[1] = line[0:11]
		c[2] = line[11:22]
		c[3] = line[22:33]
		c[4] = line[33:44]
		c[5] = line[44:55]
		c[6] = line[55:66]
		c[7] = line[66:70]
		c[8] = line[70:72]
		c[9] = line[72:75]
		c[10] = line[75:80]

		#check for reading line
		if start == True and int(c[10]) == 1 and int(c[9]) == 454:
			#print "First Data Line"
			ZA = parse_num(c[1])
			FPs.ZAID = z[ZA[1]] + '-' + str(ZA[2])
			LE = int(c[3])
		elif start == True and int(c[10]) > 1 and int(c[9]) == 454 and NewLevel == True:
			#print "Second Data Line: Repeats for each energy level"
			FPs.Levels[cur] = FPY_Energy()
			num = sci_notation(c[1])
			FPs.Levels[cur].Energy = float(num)
			FPs.Levels[cur].Yield = {}
			FPs.Levels[cur].DY = {}
			NN = int(c[5])
			NFP = int(c[6])
			NewLevel = False
			cur = cur+1
		elif start == True and int(c[10]) > 1 and int(c[9]) == 454 and NewLevel == False:
			#Here, every 4 data entries is something, but they cross over into other lines...
			if odd == True:
				#First full set
				if c[1].strip():
					ZA = parse_num(c[1])
					name = z[ZA[1]] + '-' + str(ZA[2])
					State = sci_notation(c[2])
					if float(State) > 0.0 and includeState == True:
						name+='m'

					if name in FPs.Levels[cur-1].Yield:
						FPs.Levels[cur-1].Yield[name] += float(sci_notation(c[3]))
						FPs.Levels[cur-1].DY[name] += float(sci_notation(c[4]))
					else:
						FPs.Levels[cur-1].Yield[name] = float(sci_notation(c[3]))
						FPs.Levels[cur-1].DY[name] = float(sci_notation(c[4]))
					id = ''

				#First half of second set
				if c[5].strip():
					ZA = parse_num(c[5])
					name = z[ZA[1]] + '-' + str(ZA[2])
					State = sci_notation(c[6])
					if float(State) > 0.0 and includeState == True:
						name+='m'

					if name in FPs.Levels[cur-1].Yield:
						FPs.Levels[cur-1].Yield[name] += 0.0
						FPs.Levels[cur-1].DY[name] += 0.0
					else:
						FPs.Levels[cur-1].Yield[name] = 0.0
						FPs.Levels[cur-1].DY[name] = 0.0
					id = name

				odd = False
			else:
				#Second half of second set
				if c[1].strip():
					if id in FPs.Levels[cur-1].Yield:
						FPs.Levels[cur-1].Yield[id] += float(sci_notation(c[1]))
						FPs.Levels[cur-1].DY[id] += float(sci_notation(c[2]))
					else:
						FPs.Levels[cur-1].Yield[id] = float(sci_notation(c[1]))
						FPs.Levels[cur-1].DY[id] = float(sci_notation(c[2]))

				#Last full set
				if c[3].strip():
					ZA = parse_num(c[3])
					name = z[ZA[1]] + '-' + str(ZA[2])
					State = sci_notation(c[4])
					if float(State) > 0.0 and includeState == True:
						name+='m'

					if name in FPs.Levels[cur-1].Yield:
						FPs.Levels[cur-1].Yield[name] += float(sci_notation(c[5]))
						FPs.Levels[cur-1].DY[name] += float(sci_notation(c[6]))
					else:
						FPs.Levels[cur-1].Yield[name] = float(sci_notation(c[5]))
						FPs.Levels[cur-1].DY[name] = float(sci_notation(c[6]))
					id = ''

				odd = True
			#End odd check

			iNN+=6 #Every line has 6 data entries
			if iNN >= NN:
				NewLevel = True
				odd = True
				iNN = 0
		#End if

		#Check for first starting point
		if int(c[10]) == 99999 and start == False:
			start = True
		#End c[10] check

	#End line loop

	#Done reading

	return FPs
# --------------------- End read_fpy ------------------------------

# ----------- Function to  Write output to a file --------------
def write_fpy(nfile, iso):
	#Loop through all isotopes and yields and make corrections
	yield_sum = 0.0
	cor = 1.0

	nfile.write(iso.ZAID + ':\n')
	nfile.write('---\n')
	print 'Expanding Isotope data... ' + iso.ZAID
	for levels in iso.Levels:
		yield_sum = 0.0
		unc_sum = 0.0
		nfile.write('- EL-' + str(levels) + ':\n')
		nfile.write('  Energy: ' + str(iso.Levels[levels].Energy) + '  # (eV)\n')
		nfile.write('\n')

		#Make Corrections
		for keys in iso.Levels[levels].Yield:
			yield_sum+= iso.Levels[levels].Yield[keys]
		cor = 2.0/yield_sum
		yield_sum = 0.0
		for keys in iso.Levels[levels].Yield:
			iso.Levels[levels].Yield[keys] = iso.Levels[levels].Yield[keys]*cor
			yield_sum+= iso.Levels[levels].Yield[keys]
		#End Corrections
		if abs(yield_sum-2.0) > 1E-10:
			print 'Yield Sum Error\n'
			print yield_sum
		else:
			for keys in iso.Levels[levels].Yield:
				if iso.Levels[levels].Yield[keys] != 0.0:
					nfile.write('  - ' + keys + ':\n')
					nfile.write('    Yield: ' + str(iso.Levels[levels].Yield[keys]) + '\n')
					nfile.write('    DY: ' + str(iso.Levels[levels].DY[keys]) + '\n')
					nfile.write('\n')
			#End keys loop

		nfile.write('\n')
	#End level loop

	nfile.write('...\n\n')

	return
# --------------------- End write_fpy ------------------------------

#Open yaml file to write to
nfile = open('../../database/NeutronFissionProductYields.yml', 'w')
sfile = open('../../database/SpontaneousFissionProductYields.yml', 'w')

# Create pathing variables relative to this directory
basepath = os.path.dirname(__file__)
nfy_path = os.path.join(basepath, "nfy")
sfy_path = os.path.join(basepath, "sfy")

#Preserve below for testing
#path = "nfy/nfy-092_U_235.endf"
#print 'Reading Neutron file... ' + path
#iso = read_fpy(path)
#write_fpy(nfile,iso)
#path = "nfy/nfy-092_U_238.endf"
#print 'Reading Neutron file... ' + path
#iso = read_fpy(path)
#write_fpy(nfile,iso)
#nfile.close()
#sfile.close()


# Loop through all files in path
i = 0
for filename in os.listdir(nfy_path):
	path = os.path.join(nfy_path, filename)
	print 'Reading Neutron file... ' + path
	iso = read_fpy(path)
	write_fpy(nfile,iso)
	i+=1
#End filename Loop

#Close the yaml file
print '\nNeutron FPY Library Complete!\n'
nfile.close()

# Loop through all files in path
i = 0
for filename in os.listdir(sfy_path):
	path = os.path.join(sfy_path, filename)
	print 'Reading Spontaneous file... ' + path
	iso = read_fpy(path)
	write_fpy(sfile,iso)
	i+=1
#End filename Loop

#Close the yaml file
print '\nSpontaneous FPY Library Complete!\n'
sfile.close()
