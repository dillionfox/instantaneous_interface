#ifndef CALPHI_H
#define CALPHI_H

double calphi(double x, double sig, double cutoff)
{
  //double sig = 0.24; 
  //double cutoff = 0.7; 
  double phic, C;
  double phix;

  phic = exp(-cutoff*cutoff/(2.0*sig*sig));
  C = 1.0 / ( pow(2.0*M_PI,0.5) * sig * erf(cutoff / (pow(2.0,0.5) * sig)) - 2.0*cutoff*phic );
  if (abs(x) <= cutoff)
  {
    phix = C * ( exp(-x*x/(2.0*sig*sig)) - phic );
  }
  else 
  {
    phix = 0.0;
  }
  return phix;
};

#endif



