from __future__ import division
import numpy as np
from utils import grid_interp as MC
from utils import writepdb as pdb

"""
Compute the Willard-Chandler instantaneous interface 

Reference:
Willard AP, Chandler D. Instantaneous liquid interfaces. J Phys Chem B. 2010;114(5):1954-8.

This code implements the algorithm described in the above paper to compute the
instantaneous liquid interface from atomic coordinates.

See example at the end of the code to get an idea of how to use this code.

"""

class CW:
	def __init__(self,pos,dL,pbc,n_heavy_atoms):
		self.pos = pos
		self.dL = dL
		self.pbc = pbc
		self.n_heavy_atoms = n_heavy_atoms
		self.n_grid_pts = None
		self.grid_spacing = None
		self.frame = 0

	@staticmethod
	def erf(x):
		"""
		This is a straight-forward implementation of a Numerical Recipe code for computing the error function
	
		"""
	
		sign = 1 if x >= 0 else -1
		a = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429, 0.3275911] # constants
		x = abs(x) ; t = 1.0/(1.0 + a[5]*x)
		return sign*(1.0 - (((((a[4]*t + a[3])*t) + a[2])*t + a[1])*t + a[0])*t*np.exp(-x*x))

	def phi(self, x, sig, cutoff):
		"""
		Equation 2 from Chandler paper

		"""
	
		phic = np.exp(-cutoff*cutoff/(2.0*sig*sig))
		C = 1.0 / ( (2*np.pi)**(0.5) * sig * self.erf(cutoff / (2.0**(0.5) * sig)) - 2.0*cutoff*phic )
		if np.abs(x) <= cutoff:
			return C * ( np.exp(-x*x/(2.0*sig*sig)) - phic )
		else: 
			return 0.0
	
	@staticmethod
	def gaussian_convolution(voxel_i,N,ngrid,pos,grid_spacing_i,dl,phi_bar):
		"""
		Convolve the density field with Gaussians

		"""
	
		nx = voxel_i-N			# N = int(cutoff/dL)-xg-1, xg: (0,2*int(cutoff/dL))
		if nx<0: nx+=ngrid		# wrap around periodic bounds
		elif nx>=ngrid: nx-=ngrid
		rx = np.abs(N*grid_spacing_i+(pos-voxel_i*grid_spacing_i))
		nrx = int(rx/dl)
		phix = phi_bar[nrx] + (phi_bar[nrx+1] - phi_bar[nrx]) * (rx - nrx*dl) / dl
		return phix,nx
	
	def compute_coarse_grain_density(self):
		"""
		This function takes the positions of the protein/lipids/whatever and computes the coarse grain
		density field

		"""
		nconf = 1							# number of conformations to consider
		rho_pro = 0.05							# bulk density of protein
		cutoff = 7 							# cutoff for Gaussian (for convolution)
		dl = 0.1							# grid spacing for coarse grain density
		npoints = int(cutoff/dl)+1					# number of density points
		sigp = 2.4 							# width of Gaussian smoothing: protein
		phi_bar_p = [self.phi(i*dl, sigp, cutoff) for i in range(npoints*2)]	# coarse grain density
		Ninc = int(cutoff/self.dL)						# number of density points
	
		# define voxels
		self.n_grid_pts = [int(p/self.dL) for p in self.pbc]				# number of bins that naturally fit along each axis
		self.grid_spacing =  [self.pbc[i]/self.n_grid_pts[i] for i in range(len(self.pbc))]	# distance between grid points for each direction
		rho = np.zeros((self.n_grid_pts[0],self.n_grid_pts[1],self.n_grid_pts[2]))	# dummy arrays store coarse grain density
	
		for i in range(self.n_heavy_atoms):
			pos_i = self.pos[i]
			voxel_i = [int(pos_i[dim]/self.grid_spacing[dim]) for dim in range(3)]	# convert xyz to voxel
			for xg in range(2*Ninc): # 
				phix,nx = self.gaussian_convolution(voxel_i[0],Ninc-xg-1,self.n_grid_pts[0],pos_i[0],self.grid_spacing[0],dl,phi_bar_p)	
				for yg in range(2*Ninc):
					phiy,ny = self.gaussian_convolution(voxel_i[1],Ninc-yg-1,self.n_grid_pts[1],pos_i[1],self.grid_spacing[1],dl,phi_bar_p)	
					for zg in range(2*Ninc):
						phiz,nz = self.gaussian_convolution(voxel_i[2],Ninc-zg-1,self.n_grid_pts[2],pos_i[2],self.grid_spacing[2],dl,phi_bar_p)	
						rho[int(nx)][int(ny)][int(nz)] += phix*phiy*phiz/(rho_pro*nconf)	# Equation 3
		return rho
	
	def marching_cubes(self,rho): 
		"""
		Main implementation of the Marching Cubes algorithm. This is used to smooth out the coarse grain
		density field. This code was largely translated from the original IBM code. Some of it is so
		compact that it may be difficult to interpret. I recommend tracking down the original IBM code
		if you are confused.

		"""
		# constants
		II = 1 ; epsilon = 0.000001 ; rhoc = 0.1 ; trip = np.zeros((5,3,3)) ; ii_coor = np.zeros((II,3))

		# Function to mirror points in periodic grid
		def mirror(ind,c):
			if ind[c] >= self.n_grid_pts[c_]: return ind[c]-self.n_grid_pts[c]
			else: return ind[c]
	
		# iterate through each grid point
		for i in range(self.n_grid_pts[0]):
			for j in range(self.n_grid_pts[1]):
				for k in range(self.n_grid_pts[2]):
					# gridv contains the rho values at the 8 neighboring voxels
					i1,j1,k1 = [mirror([i+1,j+1,k+1],c_) for c_ in range(3)]
					ind2 = [[i,j,k],[i,j1,k],[i1,j1,k],[i1,j,k],[i,j,k1],[i,j1,k1],[i1,j1,k1],[i1,j,k1]]
					gridv = [rho[q1][q2][q3] for q1,q2,q3 in ind2]
					
					# find if the cube is inside bubble, and whether it is near a heavy atom
					cubefactor = sum([1 if gridv[v] <= rhoc else 0 for v in range(8)])

					# if next to heavy atom
					if cubefactor >=4: cube_coor = np.einsum('i,i->i',[i+0.5,j+0.5,k+0.5],self.grid_spacing)
	
					# define grid 
					ind = [[0,0,0],[0,1,0],[1,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,1],[1,0,1]]
					gridp = [np.einsum('i,i->i',[i+q1,j+q2,k+q3],self.grid_spacing) for q1,q2,q3 in ind]
					[ntri,trip] = MC.MC_table(gridv, gridp, rhoc, trip)
	
					for nt in range(ntri):
						for qi in range(3):
							vertexflag = False
	                                                for ii in range(II):
								vertexdist = np.dot(trip[nt][qi]-ii_coor[ii],trip[nt][qi]-ii_coor[ii])
								if vertexdist < epsilon: vertexflag = True ; break
	                                                if vertexflag == True: break
	                                                else:
								ii_coor = np.vstack((ii_coor, np.zeros(3)))
								ii_coor[II] = trip[nt][qi]
								II+=1
		return ii_coor

	def run(self):
		rho = self.compute_coarse_grain_density()
		interface_coors = self.marching_cubes(rho)
		pdb.write_pdb(interface_coors,self.frame)
		return interface_coors

if __name__ == "__main__":
	"""
	Example of how you can use MDAnalysis to extract necessary information to use this class
	Note: MDAnalysis is just one way to extract the required information.

	"""
	import MDAnalysis
	DCD = "/home/dillion/data/HRF/constrained_atoms/GOL/GOL_10000-15000_100.dcd"
	PSF =  "/home/dillion/data/HRF/constrained_atoms/GOL/beta.psf"

        uni = MDAnalysis.Universe(PSF,DCD)							# instantiate MDAanalysis universe
        uni.trajectory[0] ; pbc = uni.dimensions[0:3] ; dL = 1					# go to first frame, get pbc, set constant
	sel_all = uni.select_atoms('all')							# select all atoms
        sel_all.atoms.translate(-sel_all.atoms.center_of_mass()+pbc/2.0)			# center all atoms in system. This is important!!
        heavy_atoms = uni.select_atoms("resname GOL").select_atoms('not name H*')		# Only need to consider heavy atoms

	ii = CW(heavy_atoms.positions,dL,pbc,len(heavy_atoms.atoms))				# instantiate CW class
	interface_coors = ii.run()								# compute the instantaneous interface
	print interface_coors									# print interface coordinates
