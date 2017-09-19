from __future__ import division
import numpy as np
from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory
import time
import MDAnalysis
import CW_interface

#--- WRITE OUTPUT
def write_pdb(coor,fr):
	"""
	This function writes the coordinates of the Willard-Chandler instantaneous interface as a pdb,
	and populates the Beta field with the long-range electrostatic potential
	"""
	global name_modifier
	global verbose

	if verbose >= 1:
		print 'writing pdb...'
	outfile = open(str(fr)+"_II.pdb","w")
	count_zeros = 0
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

#--- EXTRACT COORDINATES
def extract_traj_info(PSF,DCD,selection_key):
	"""
	This function uses MDAnalysis to extract coordinates from the trajectory
	"""
	global verbose

	if verbose >= 1:
        	print 'loading coordinates...'
        # load some variables into global namespace
        global n_heavy_atoms
        global pbc
        global pbc_fr
	global box_shift
	global first_frame
	global last_frame
	global water_names

        uni = MDAnalysis.Universe(PSF,DCD)
        nframes = len(uni.trajectory)					# number of frames
	box_shift = np.zeros((nframes,3))
	pbc_fr = np.zeros((nframes,3))

        protein = uni.select_atoms(selection_key)			# identify atoms to build interface around
        heavy_atoms = protein.select_atoms('not name H*')		# Only need to consider heavy atoms
        protein_indices = heavy_atoms.indices 
        n_heavy_atoms = len(heavy_atoms.atoms)                          # number of heavy protein atoms
        positions = np.zeros((nframes,n_heavy_atoms,3))

        water = uni.select_atoms("resname TIP3")			# identify atoms to build interface around
        water_indices = water.indices 
	water_names = water.names
        n_water = len(water.atoms)					# number of heavy protein atoms
        water_pos = np.zeros((nframes,n_water,3))

	frames = range(nframes)
        for fr in frames:                                       # save coordinates for each frame
                uni.trajectory[fr]
        	pbc = uni.dimensions[0:3]				# retrieve periodic bounds
        	pbc_fr[fr] = pbc
        	sel = uni.select_atoms('all')
		box_shift[fr] = -sel.atoms.center_of_mass()+pbc/2.0	# first center at origin, then put vertx of quadrant 7 at origin
        	sel.atoms.translate(box_shift[fr])			# center selection

        	protein = uni.select_atoms(selection_key)		# identify atoms to build interface around
        	heavy_atoms = protein.select_atoms('not name H*')	# Only need to consider heavy atoms
		positions[fr] = heavy_atoms.positions

        	water = uni.select_atoms("resname TIP3")		# identify atoms to build interface around
		water_pos[fr] = water.positions

        pbc = uni.dimensions[0:3]					# retrieve periodic bounds
        return [nframes,positions,water_pos]

def run_emaps(fr):
	"""
	This function takes in a set of II points for each frame and calculates the LREP for each frame using the
	pertinent set of coordinates.
	"""
	global verbose
	global box_shift
	global positions
	global water
	global nframes

	#--- EXTRACT INFO FOR FRAME, FR
	if verbose >= 3:
		print 'working on frame', fr+1, ' of', nframes
	pos=positions[fr] 
	water_coor = water[fr]
	
	#--- COMPUTE RHO
	coarse_grain_start = time.time()
	rho = CW_interface.compute_coarse_grain_density(pos,dL,pbc,n_heavy_atoms) # defines global variables: n_grid_pts, grid_spacing
	coarse_grain_stop = time.time()
	if verbose >= 2:
		print 'elapsed time to compute coarse grain density:', coarse_grain_stop-coarse_grain_start
	
	#--- MARCHING CUBES
	marching_cubes_start = time.time()
	interface_coors = CW_interface.marching_cubes(rho) # defines global variables: cube_coor
	marching_cubes_stop = time.time()
	if verbose >= 2:
		print 'elapsed time to run marching cubes:', marching_cubes_stop-marching_cubes_start
	
	write_pdb(interface_coors,fr)
	return 0

def instant_interface(grofile,trajfile,**kwargs):
	"""
	This is the MAIN function. It reads in all data and decides which functions should be called to
	compute either averaged or individual electrostatic maps.
	"""
	#--- READ IN PARAMETERS FROM YAML FILE
	global selection_key
	global water_resname
	global dL
	global name_modifier
	global verbose
	global nthreads
	global first_frame
	global last_frame

	selection_key = "resname GOL or resname TRG or resname TAS or resname GLB"
	dL = 1 
	av_LREP = 'n'
	water_resname = "resname TIP3"
	name_modifier = ""
	verbose = 3
	nthreads = 1
	writepdb = 'n'

	#--- LOAD VARIABLES INTO GLOBAL NAMESPACE
	global box_shift
	global positions
	global water
	global nframes

	#--- READ DATA
	[nframes, positions, water] = extract_traj_info(grofile,trajfile,selection_key)
	Parallel(n_jobs=nthreads)(delayed(run_emaps,has_shareable_memory)(fr) for fr in range(nframes))

traj = "/home/dillion/research/HRF/constrained_atoms/TAS/TAS.dcd"
psf =  "/home/dillion/research/HRF/constrained_atoms/TAS/beta.psf"

instant_interface(psf,traj)
 
