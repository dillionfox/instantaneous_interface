#ifndef READPRO_H
#define READPRO_H

int* readpro(FILE* infile, int nheavy, int npro)
{
  int* index = new int[nheavy];
  char str[100];
  char atmname[5];
  char hydrogen[] = "H";  
 
  int hi = 0;
  fgets(str, 100, infile);
  fscanf (infile, "%*d\n");

  for ( int i=0; i<npro; i++)
  {
    fscanf (infile, "%6s %6d", atmname,&index[hi]);
    //cout << "index: " << index[hi] << ", atom name: " << atmname << endl;
    if ( strncmp(atmname, hydrogen, 1) != 0 )
    {
      //fscanf (infile, "%5d%*8.3f%*8.3f%*8.3f%*8.4f%*8.4f%*8.4f", &index[hi]);
      hi ++;
     
      cout << "atom name: " << atmname << ", index: " << index[hi-1] << endl;
    }
    //else
    //  fscanf (infile, "%*d %*f %*f %*f %*f %*f %*f\n");
  }
  //cout << "index: " << index << ", atom name: " << atmname << endl;
  return index;
};

#endif
