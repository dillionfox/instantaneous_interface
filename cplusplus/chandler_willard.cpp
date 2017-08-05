#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <ctime>
#include <malloc.h>
#include </home/dillion/pkg/surf-master/depends/xdrfile-1.1.4/include/xdrfile_xtc.h> // xdr include file

using namespace std;
#include "calphi.h"
#include "array_3D.h"
#include "readpro.h"
#include "Marching_Cube.h"

//******************************************************************
// ::Version 2::
// Dillion Fox, 2017
//
// Original code written by R. Remsing, 2014
//
// This code will calculate the instantaneous interface
// for the protein & water atoms following Willard and Chandler
// (Willard, Chandler, J. Phys. Chem. B 2010, 114, 1954).
//
// R. Remsing, 2014
//******************************************************************

// To run, make sure the .gro and .xtc files listed below are correct
//
// Then compile with the compile script: compile_xtc.sh protein_proIIcal
// and run: ./protein_proIIcal.out > IIprocal.dat
//
// NOTE: location of libraries, etc. in compile script may need to be modified.



int main(int argc, char **argv) {
  // USER INPUT
  int npro = 14;		// Number of protein atoms (total)
  int nheavy = 7;		// Number of protein heavy atoms
  int natoms = 16645;           // number of atoms
  int lastfr = 6;		// Last frame of the trajectory to analyze. Here, we keep the protein rigid,
  int NII = 2500;		// Number of points on the instantaneous interface


  /* XTC variables */
  XDRFILE *xdin;                // the xdr file holder
  XDRFILE *xdout;
  int step;                     // the step counter for the simulation
  float time;                   // simulation time
  matrix box;                   // box coordinates in a 3x3 matrix
  rvec *xtc_coor;               // atom coordinates in a 2-D matrix
  float prec;                   // precision of the xtc file

  /* Program variables */
  int frame = 0;
  int cutfr = 0;               // Number of frames to skip at the beginning of the trajectory
    
  int countfr = 0;
  int skipfr = 1;              // Skip every <skipfr> frames during the analysis
  double nconf = 2.0;          // Total number of configurations to analyze.

  // look up table
  double l;
  double sigw = 0.24;          // width of the Gaussian for smoothing the water density
  double sigp = 0.24;          // width of the Gaussian for smoothing the protein heavy atom density
  double cutoff = 0.7;         // ///// Cutoff for the Gaussian
  double dl = 0.01;            // ///// grid spacing
  int npoints = int(cutoff/dl) + 1;
  double LUT_phiw[npoints];
  double LUT_phip[npoints];
  for (int i=0; i<npoints; i++)
  {
    l = i*dl;
    LUT_phiw[i] = calphi(l, sigw, cutoff); // Precalculate the smoothing functions
    LUT_phip[i] = calphi(l, sigp, cutoff);
  }

  // find the index of the protein heavy atoms
  FILE* infile;
  infile = fopen (argv[2], "r");                              // Gro file to read in protein heavy atoms, etc.
  int* index_pro = readpro(infile, nheavy, npro);
  int iheavy;
  
  // set up grids
  double dL = 0.1;
  double dgrid[3];
  int ngrids[3];
  int Ncubes;
  int cubeid;
  double* cubedata;

  // calculate density
  double*** rhow;
  double*** rhop;
  double*** rho;
  int index[3];
  int Ninc = int(cutoff/dL);
  double phix, phiy, phiz;
  double rx, ry, rz;
  int nx, ny, nz;
  int nrx, nry, nrz;

  // search interface
  double rho_water = 33.4; // #/nm^3            // Bulk water density
  double rho_pro = 50.0;                        // Protein heavy atom "bulk density," determined through trial and error by examining smoothed densities
  double rhoc = 0.1;                            // Cutoff density in units of the bulk density, usually 0.5, but here set to 0.1. Small changes in this will not qualitatively change results.
  int II = 0;
  rvec* ii_coor;
  int i1, j1, k1;
  int cubefactor;
  int vertexflag;
  int cubecount;
  double vertexdist;
  double avg_nii = 0.0;
  double avg_area = 0.0;
  double avg_vol = 0.0;
  double area, volume;
  double a2, b2, c2;
  double cube_coor[3];
  double rwater, rcube;
  double* heavy_ncube = (double*) calloc(nheavy, sizeof(double));
  double* heavy_nwater = (double*) calloc(nheavy, sizeof(double));

  /* Marching Cube variables */
  int ntri;
  CUBEINFO grid;
  TRIANGLE tri[5];
  
  /* Memory allocation to read coordinates */
  read_xtc_natoms(argv[1], &natoms);
  xtc_coor = (rvec *) malloc(natoms*sizeof(rvec));
  if(xtc_coor==0){
    cout<<"Insufficient memory to load .xtc file.\n";
    return 0;
  }

  /* Open the xtc file and loop through each frame. */
  xdin=xdrfile_open(argv[1],"r");
  xdout=xdrfile_open("protII.xtc","w");

  printf ("# Time\tNii\tavg_nii \t area\tavg_area\tvolume\tavg_vol\n");

  while( ! read_xtc(xdin, natoms, &step, &time, box, xtc_coor, &prec) ) 
  {
    frame ++;
    // using fixed grid in order to calculate avg. (WARNING: this will leave a defect on the pbc boundries, only applicable when boundary density is not relavent)
    if (frame == 1)
    {
      for (int i=0; i<3; i++)
      {
        ngrids[i] = int(box[i][i]/dL);
        dgrid[i] = box[i][i]/ngrids[i];
      }
      Ncubes = ngrids[0] * ngrids[1] * ngrids[2];
      cubedata = (double *) calloc(Ncubes, sizeof(double));
      rhow = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
      rhop = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
      rho = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
    }
      
    if (frame > cutfr && frame <= lastfr && (frame-1)%skipfr == 0 )
    {
      countfr ++;
      for (int i=0; i<nheavy; i++)  // compute coarse grained density from protein heavy atoms
      {                             // the indices of the heavy atoms is prepared in advance (above)
        iheavy = index_pro[i] - 1;
        for (int dim=0; dim<3; dim++)
          index[dim] = int( xtc_coor[iheavy][dim]/dgrid[dim] );
	cout << index[0] << " " << index[1] << " " << index[2] << endl;
	cout << xtc_coor[iheavy][0] << " " << xtc_coor[iheavy][1] << " " << xtc_coor[iheavy][2] << endl;
        for (int xg=0; xg<2*Ninc; xg++)
        {
          nx = index[0] - Ninc + 1 + xg;
          if (nx < 0) 
            nx += ngrids[0];
          else if (nx >= ngrids[0]) 
            nx -= ngrids[0];

          rx = abs( (Ninc-1-xg)*dgrid[0] + (xtc_coor[iheavy][0] - index[0]*dgrid[0]) );
          nrx = int(rx/dl);
          phix = LUT_phip[nrx] + (LUT_phip[nrx+1] - LUT_phip[nrx]) * (rx - nrx*dl) / dl;
          for (int yg=0; yg<2*Ninc; yg++)
          {
            ny = index[1] - Ninc + 1 + yg;
            if (ny < 0) 
              ny += ngrids[1];
            else if (ny >= ngrids[1]) 
              ny -= ngrids[1];

            ry = abs( (Ninc-1-yg)*dgrid[1] + (xtc_coor[iheavy][1] - index[1]*dgrid[1]) );
            nry = int(ry/dl);
            phiy = LUT_phip[nry] + (LUT_phip[nry+1] - LUT_phip[nry]) * (ry - nry*dl) / dl;
            for (int zg=0; zg<2*Ninc; zg++)
            {
              nz = index[2] - Ninc + 1 + zg;
              if (nz < 0) 
                nz += ngrids[2];
              else if (nz >= ngrids[2]) 
                nz -= ngrids[2];

              rz = abs( (Ninc-1-zg)*dgrid[2] + (xtc_coor[iheavy][2] - index[2]*dgrid[2]) );
              nrz = int(rz/dl);
              phiz = LUT_phip[nrz] + (LUT_phip[nrz+1] - LUT_phip[nrz]) * (rz - nrz*dl) / dl;
              rhop[nx][ny][nz] += phix*phiy*phiz;
            }
          }
        }
      }

      // normalize and combine the density distribution 
      for (int i=0; i<ngrids[0]; i++)
      {
        for (int j=0; j<ngrids[1]; j++)
        {
          for (int k=0; k<ngrids[2]; k++)
          {
            rho[i][j][k] = rhop[i][j][k]/rho_pro; //rhow[i][j][k]/rho_water;
          }
        }
      }
    } //} ~~~~~ ENDS IF LOOP ~~~~~

    for (int i=0; i<ngrids[0]; i++)
    {
      for (int j=0; j<ngrids[1]; j++)
      {
        for (int k=0; k<ngrids[2]; k++)
        {
          rho[i][j][k] = rho[i][j][k]/(nconf);
        }
      }
    }
  
    // search for instantaneous interface using the marching cube algorithm
    ii_coor = (rvec *) calloc(NII, sizeof(rvec));
    II = 0;
    area = 0.0;
    cubeid = 0;
    cubecount = 0;
  
    for (int i=0; i<ngrids[0]; i++)
    {
      for (int j=0; j<ngrids[1]; j++)
      {
        for (int k=0; k<ngrids[2]; k++)
  	{
          i1 = i + 1;
          j1 = j + 1;
          k1 = k + 1;
          if (i1 >= ngrids[0]) i1 -= ngrids[0];
          if (j1 >= ngrids[1]) j1 -= ngrids[1];
          if (k1 >= ngrids[2]) k1 -= ngrids[2];
          
          grid.v[0] = rho[i][j][k];
          grid.v[1] = rho[i][j1][k];
          grid.v[2] = rho[i1][j1][k];
          grid.v[3] = rho[i1][j][k];
          grid.v[4] = rho[i][j][k1];
          grid.v[5] = rho[i][j1][k1];
          grid.v[6] = rho[i1][j1][k1];
          grid.v[7] = rho[i1][j][k1];
          
          // find if the cube is inside bubble, and whether it is near a heavy atom
          cubefactor = 0;
          for (int v=0; v<8; v++)
          {
            if (grid.v[v] <= rhoc)
            {
              cubefactor ++;
            }
          }
          if (cubefactor >=4)
          {
            cubecount ++;
            cubedata[cubeid] ++;
            cube_coor[0] = (i+0.5) * dgrid[0];
            cube_coor[1] = (j+0.5) * dgrid[1];
            cube_coor[2] = (k+0.5) * dgrid[2];
            
          for (int l=0; l<nheavy; l++)
            {
              iheavy = index_pro[l] - 1;
              rcube = 0.0;
              rcube += (cube_coor[0] - xtc_coor[iheavy][0]) * (cube_coor[0] - xtc_coor[iheavy][0]);
              rcube += (cube_coor[1] - xtc_coor[iheavy][1]) * (cube_coor[1] - xtc_coor[iheavy][1]);
              rcube += (cube_coor[2] - xtc_coor[iheavy][2]) * (cube_coor[2] - xtc_coor[iheavy][2]);
              if (rcube <= 0.36)
                heavy_ncube[l] += 1.0;
            }
          
          }
          
          grid.p[0].x = i*dgrid[0];
          grid.p[0].y = j*dgrid[1];
          grid.p[0].z = k*dgrid[2];
          
          grid.p[1].x = i*dgrid[0];
          grid.p[1].y = (j+1)*dgrid[1];
          grid.p[1].z = k*dgrid[2];
          
          grid.p[2].x = (i+1)*dgrid[0];
          grid.p[2].y = (j+1)*dgrid[1];
          grid.p[2].z = k*dgrid[2];
          
          grid.p[3].x = (i+1)*dgrid[0];
          grid.p[3].y = j*dgrid[1];
          grid.p[3].z = k*dgrid[2];
          
          grid.p[4].x = i*dgrid[0];
          grid.p[4].y = j*dgrid[1];
          grid.p[4].z = (k+1)*dgrid[2];
          
          grid.p[5].x = i*dgrid[0];
          grid.p[5].y = (j+1)*dgrid[1];
          grid.p[5].z = (k+1)*dgrid[2];
          
          grid.p[6].x = (i+1)*dgrid[0];
          grid.p[6].y = (j+1)*dgrid[1];
          grid.p[6].z = (k+1)*dgrid[2];
          
          grid.p[7].x = (i+1)*dgrid[0];
          grid.p[7].y = j*dgrid[1];
          grid.p[7].z = (k+1)*dgrid[2];
          
          ntri = marching_cube(grid, rhoc, tri);
          for (int t=0; t<ntri; t++)
          {
            a2 = (tri[t].p[0].x-tri[t].p[1].x)*(tri[t].p[0].x-tri[t].p[1].x) + (tri[t].p[0].y-tri[t].p[1].y)*(tri[t].p[0].y-tri[t].p[1].y) + (tri[t].p[0].z-tri[t].p[1].z)*(tri[t].p[0].z-tri[t].p[1].z);
            b2 = (tri[t].p[1].x-tri[t].p[2].x)*(tri[t].p[1].x-tri[t].p[2].x) + (tri[t].p[1].y-tri[t].p[2].y)*(tri[t].p[1].y-tri[t].p[2].y) + (tri[t].p[1].z-tri[t].p[2].z)*(tri[t].p[1].z-tri[t].p[2].z);
            c2 = (tri[t].p[2].x-tri[t].p[0].x)*(tri[t].p[2].x-tri[t].p[0].x) + (tri[t].p[2].y-tri[t].p[0].y)*(tri[t].p[2].y-tri[t].p[0].y) + (tri[t].p[2].z-tri[t].p[0].z)*(tri[t].p[2].z-tri[t].p[0].z);
            area += 1.0/4.0 * sqrt( 2.0*(a2*b2 + b2*c2 + c2*a2) - a2*a2 - b2*b2 - c2*c2 );
          }
          
          for (int l=0; l<ntri; l++)
          {
            for (int m=0; m<3; m++)
            {
              vertexflag = 0;
              for (int n=0; n<II; n++)
              {
                vertexdist = 0.0;
                vertexdist += (tri[l].p[m].x - ii_coor[n][0]) * (tri[l].p[m].x - ii_coor[n][0]);
                vertexdist += (tri[l].p[m].y - ii_coor[n][1]) * (tri[l].p[m].y - ii_coor[n][1]);
                vertexdist += (tri[l].p[m].z - ii_coor[n][2]) * (tri[l].p[m].z - ii_coor[n][2]);
                if (vertexdist < 1e-6)
                {
                  vertexflag = 1;
                  break;
                }
              }
            
              if (vertexflag == 0)
              {
                ii_coor[II][0] = tri[l].p[m].x;
                ii_coor[II][1] = tri[l].p[m].y;
                ii_coor[II][2] = tri[l].p[m].z;
                II ++;
              }
            }
          } 
          cubeid ++;
        }// loop over k 
      }// loop over j
    }//loop over i

    volume = cubecount*dgrid[0]*dgrid[1]*dgrid[2];           // If II is a bubble, this will compute its volume
    avg_nii = ( (countfr-1) * avg_nii + II ) / countfr;      // Avg number of points on the II
    avg_area = ( (countfr-1) * avg_area + area ) / countfr;  // Area of II (average)
    avg_vol = ( (countfr-1) * avg_vol + volume ) / countfr;  // Volume of bubble contained by II (average)
    printf ("%8.1f%10d%10.0f%10.3f%10.3f%10.3f%10.3f\n", time, II, avg_nii, area, avg_area, volume, avg_vol);
    
    // print statements to help trouble shooting
    // save II coordinates in xtc file
    if (II > NII)
      cout << "Need to assign more space to store II ! at time " << II << " " << NII << " " << time << endl;
    else
      write_xtc(xdout, NII, step, time, box, ii_coor, prec);
    // generate a gro file for the first frame counted
    cout << "~~~~~~ countfr = " << countfr << ", time = " << time << " ~~~~~~" << endl;

    FILE* outfile;
    string filename;
    string head = "protII_time";
    string ext = ".gro";
    stringstream ss;
    ss << time;
    filename = head + ss.str() + ext;
    outfile = fopen (filename.c_str(), "w");
    fprintf (outfile, "Instantaneous interfaces by Marching Cube\n");
    fprintf (outfile, "%-10d\n", NII);
    
    for (int i=0; i<NII; i++)
      fprintf (outfile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i+1, "TRI", "II", i+1, ii_coor[i][0], ii_coor[i][1], ii_coor[i][2]);
    
    fprintf (outfile, "%10.5f%10.5f%10.5f\n", box[0][0], box[1][1], box[2][2]);
    fclose (outfile);
    
    free (ii_coor);
    ii_coor = NULL;
  }
  deallocate3D(rhow, ngrids[0], ngrids[1], ngrids[2]);
  deallocate3D(rhop, ngrids[0], ngrids[1], ngrids[2]);
  deallocate3D(rho, ngrids[0], ngrids[1], ngrids[2]);

  xdrfile_close(xdout);
  xdrfile_close(xdin);
  return 0;
}

