#ifndef CONSTRAINTSTRILATERATION_H
#define CONSTRAINTSTRILATERATION_H

#include <KinFit.h>
#include <reconstructor.h>
#include <plane_intersection.h>

/*List of variables (order is the following)*/
/*
  For 4 gamma decay, without charged decay:
  4 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 4 gamma decay, with charged decay:
  2 x 3 (CurvOld, PhivOld, CotvOld)
  4 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kch momentum)
  Charged kaon 4-mom determined during fit (using the boost method as well)
  Neutral kaon 4-vec determined during fit
  Neutral kaon dependent of Kch direction and energy

--------------------------------------------------------------------------------

  For 6 gamma decay, without charged decay:
  6 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 6 gamma decay, with charged decay:
  2 x 3 (CurvOld, PhivOld, CotvOld)
  6 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kch momentum)
  Charged kaon 4-mom determined during fit (using the boost method as well)
  Neutral kaon 4-vec determined during fit
  Neutral kaon dependent of Kch direction and energy

--------------------------------------------------------------------------------

*/

namespace KLOE
{
  /**
   * @class ConstraintsTrilateration
   * @brief Auxiliary class with the constraints for trilateration kinematic fit - do not include charged part of the event
   */
  class ConstraintsTrilateration : public KinFit, public ChargedVtxRec<Float_t, Int_t>
  {
  private:
    // Path conservation in LAB
    Double_t NeutralPathConsvLAB(Double_t *x, Double_t *p);

    // Fictitious overriders
    Double_t FourMomConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t PhotonPathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t PxConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t PyConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t PzConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t EnergyConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t Photon1PathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t Photon2PathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t Photon3PathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t Photon4PathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    Double_t MinvConsvChKaon(Double_t *x, Double_t *p) { return 0; };
    Double_t MinvConsvNeuKaon(Double_t *x, Double_t *p) { return 0; };
    Double_t MinvConsvOmega(Double_t *x, Double_t *p) { return 0; };

    Int_t
        _chosenComponent,
        _selected[4] = {1, 2, 3, 4};
    Bool_t
        _cond_detector = kFALSE;

  protected:
    void SetParameters(Double_t *p);
    void ResetParameters();
    void IntermediateReconstruction();

  public:
    /* Specific physical Constraints for trilateration kinematic fit */

    neutralParticle fphoton[4];
    chargedParticle fpionCh[2];
    kaonNeutral fKchrec, fKchboost, fKnerec, fKnerecCMPhi, fKnereclor;
    phiMeson fphi;
    std::vector<Float_t> fip;

    void IntermediateReconstruction(Double_t *p);

    // Energy conservation in CM
    Double_t EnergyConsvCM(Double_t *x, Double_t *p);

    // Invariant mass of Kaon conservation
    Double_t MinvConsv(Double_t *x, Double_t *p);

    // Gamma path of flight from IP Conservation
    Double_t NeutralXPathConsvLAB(Double_t *x, Double_t *p)
    {
      _chosenComponent = 0;
      return NeutralPathConsvLAB(x, p);
    };
    Double_t NeutralYPathConsvLAB(Double_t *x, Double_t *p)
    {
      _chosenComponent = 1;
      return NeutralPathConsvLAB(x, p);
    };
    Double_t NeutralZPathConsvLAB(Double_t *x, Double_t *p)
    {
      _chosenComponent = 2;
      return NeutralPathConsvLAB(x, p);
    };
  };

}
#endif
