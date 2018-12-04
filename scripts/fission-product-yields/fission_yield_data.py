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

# ----------- Function to read ENDF-6 format for FPY --------------
def read_fpy(filepath):
	start = False
	c = {}
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
		#print c[1] + ' ' + c[2] + ' ' + c[3] + ' ' + c[4] + ' ' + c[5] + ' ' + c[6] + ' ' + c[7] + ' ' + c[8] + ' ' + c[9] + ' ' + c[10]
		
		#check for reading line
		if start == True and int(c[10]) == 1:
			print "First Data Line"
			print c[1] + ' ' + c[2] + ' ' + c[3] + ' ' + c[4] + ' ' + c[5] + ' ' + c[6] + ' ' + c[7] + ' ' + c[8] + ' ' + c[9] + ' ' + c[10]
		elif start == True and int(c[10]) == 2:
			print "Second Data Line"
			print c[1] + ' ' + c[2] + ' ' + c[3] + ' ' + c[4] + ' ' + c[5] + ' ' + c[6] + ' ' + c[7] + ' ' + c[8] + ' ' + c[9] + ' ' + c[10]
		elif start == True and int(c[10]) > 2:
			#Here, every 4 data entries is something, but they cross over into other lines...
			print c[1] + ' ' + c[2] + ' ' + c[3] + ' ' + c[4] + ' ' + c[5] + ' ' + c[6] + ' ' + c[7] + ' ' + c[8] + ' ' + c[9] + ' ' + c[10]
		
		#Check for first starting point
		if int(c[10]) == 99999 and start == False:
			start = True
		#End c10 check

	#End line loop
	
	return
# --------------------- End read_fpy ------------------------------

# Create pathing variables relative to this directory
basepath = os.path.dirname(__file__)
nfy_path = os.path.join(basepath, "nfy")
sfy_path = os.path.join(basepath, "sfy")

read_fpy("nfy/nfy-092_U_235.endf")

# Loop through all files in path
#for filename in os.listdir(nfy_path):
#	path = os.path.join(nfy_path, filename)
#	read_fpy(path)
#End filename Loop

# NIST data -------------------------------------------------------------

#def split_line(line):
#    return map(str.strip, line.split('='))
#
#def parse_one_chunk(chunk):
#    d = {}
#    for line in chunk:
#        k,v = split_line(line)
#        if v.find('.') >= 0:
#            if v.endswith('#'):
#                v = v[:-1]
#            d[k] = unc.ufloat_fromstr(v)
#        else:
#            try:
#                d[k] = int(v)
#            except ValueError:
#                d[k] = v
#
#    return d

# NIST data file
#data_file = os.path.join(basepath, "nist-nuclide-data.txt")


# ENDF-6 MAT data -------------------------------------------------------
#  mats is dictionary with
#    key : (Z, A, metastable), Z, A are int, metastable is bool
#    value : MAT nuclide id, integer, from ENDF-6 list

#mat_file = os.path.join(basepath, "n-ENDF-B-VII.1.endf.list")
#mats = {}
#for line in open(mat_file):
#    # Skip comment line
#    if line.startswith('#'): continue
#
#    # Grab Z, A, and MAT
#    Z = int(line[6:9])
#    A = int(line[13:16])
#    mat = int(line[72:76])
#
#    # Is it metastable?
#   metastable = (line[16] == 'M')
#
#   key = (Z, A, metastable)
#
#   mats[key] = int(mat)

# ---------------------------------------------------------------------------- #
# means intended for public access of data
