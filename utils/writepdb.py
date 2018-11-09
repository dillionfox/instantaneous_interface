
#--- WRITE OUTPUT
def write_pdb(coor,fr):
	"""
	Simple function to write interface coordinates as a PDB

	"""
	outfile = open(str(fr)+"_II.pdb","w")
	for i in range(len(coor)):
		if (coor[i][0]!=0 and coor[i][1]!=0 and coor[i][2]!=0):
			t1 = "ATOM"					# ATOM
			t2 = 1						# INDEX
			t3 = "C"					# ATOM NAME
			t4 = ""						# ALTERNATE LOCATION INDICATOR
			t5 = "AAA"					# RESIDUE NAME
			t6 = "X"					# CHAIN
			t7 = 0						# RESIDUE NUMBER
			t8 = ""						# INSERTION CODE
			t9 = float(coor[i][0])				# X
			t10 = float(coor[i][1])				# Y
			t11 = float(coor[i][2])				# Z
			t12 = 0.0					# OCCUPANCY
			t13 = 0.0					# TEMPERATURE FACTOR
			t14 = ""					# ELEMENT SYMBOL
			t15 = ""					# CHARGE ON ATOM
			outfile.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15))
	outfile.close()
	return 0
