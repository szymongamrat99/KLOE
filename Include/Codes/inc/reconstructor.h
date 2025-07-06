#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <kloe_class.h>

namespace KLOE
{
  struct Solution
  {
    /*
      1,2
      x,y,z,t
    */
    float sol[2][4];
    // indicates error computing solution
    bool error[2];
  };

  class Reconstructor
  {

  public:
    void SetClu(int i, float x, float y, float z, float t, float E);
    Solution MySolve(int *selected);
    Solution KleusbergSolve(int *selected);
    Solution LeastSquaresSolve(float *x0);
    Solution MinuitSolve(float *x0);

    float ResidualErr(int i, float const *x);
    float ResidualErrTot(float const *x);
    float CombEnergy(int *selected);
    float TotalEnergy() const;
    float GetInvMasses(float *sol, int *comb, float *Minvgg) const;
    float GetInvMassDiscrepancy(float *sol, int *comb) const;
    void GetKmomentum(const float *sol, float *p) const;

  private:
    /*
      x,y,z,t
     */
    float _clu[6][4];
    float _ene[6];
  };

}

#endif
