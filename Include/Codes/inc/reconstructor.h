#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <TMath.h>

struct Solution{
  /*
    1,2
    x,y,z,t
  */
  Double_t sol[2][4];
  // indicates error computing solution
  bool error[2];
};

class Reconstructor{
  
 public:
  void SetClu(Int_t i, Double_t x, Double_t y, Double_t z, Double_t t, Double_t E);
  Solution MySolve(Int_t * selected);
  Solution KleusbergSolve(Int_t * selected);
  Solution LeastSquaresSolve(Double_t * x0);
  Solution MinuitSolve(Double_t *x0);
  
  Double_t ResidualErr(Int_t i, Double_t const * x);
  Double_t ResidualErrTot(Double_t const * x);
  Double_t CombEnergy(Int_t * selected);
  Double_t TotalEnergy()const;
  Double_t GetInvMasses(Double_t * sol, Int_t * comb, Double_t * Minvgg)const;
  Double_t GetInvMassDiscrepancy(Double_t * sol, Int_t * comb)const;
  void GetKmomentum(const Double_t * sol, Double_t * p)const;
  
 private:
  /*
    x,y,z,t
   */
  Double_t _clu[6][4];
  Double_t _ene[6];
};

#endif
