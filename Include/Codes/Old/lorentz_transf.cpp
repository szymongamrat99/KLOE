#include "lorentz_transf.h"

void lorentz_transf(Float_t *boost_vec, Float_t *vec_init, Float_t *vec_end)
{
  TLorentzVector before(vec_init);
  TVector3 boost(boost_vec);

  before.Boost(boost);

  for(Int_t i = 0; i < 4; i++)
    vec_end[i] = before(i);
}

void lorentz_transf(Double_t *boost_vec, Double_t *vec_init, Double_t *vec_end)
{
  TLorentzVector before(vec_init);
  TVector3 boost(boost_vec);

  before.Boost(boost);

  for(Int_t i = 0; i < 4; i++)
    vec_end[i] = before(i);
}