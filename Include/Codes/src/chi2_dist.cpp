#include <chi2_dist.h>

Double_t chi2dist(Double_t *xx, Double_t *p)
{
  Double_t value = 0.;
  Double_t x = xx[0];

  Double_t deg = floor(p[1]);

  if(x > 0)
  {
    value = p[0]*(pow(x,(deg/2.)-1.)*exp(-x/2.))/(pow(2,deg/2.)*tgamma(deg/2.));
  }
  else
  {
    value = 0.;
  }

  return value;
}