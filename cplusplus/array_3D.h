#ifndef ARRAY_3D_H
#define ARRAY_3D_H

double*** allocate3D(int x, int y, int z)
{
  double*** the_array = new double**[x];
  for ( int i=0; i<x; i++ )
  {
    the_array[i] = new double*[y];
    for ( int j=0; j<y; j++ )
    {
      the_array[i][j] = new double[z];
      for ( int k=0; k<z; k++ )
      {
        the_array[i][j][k] = 0.0;
      }
    }
  }
  return the_array;
};

void deallocate3D(double*** the_array, int x, int y, int z)
{
  for ( int i=0; i<x; i++ )
  {
    for ( int j=0; j<y; j++ )
    {
      delete [] the_array[i][j];
    }
    delete [] the_array[i];
  }
  delete [] the_array;
};

#endif
