#ifndef KINFIT_H
#define KINFIT_H

#include <TF1.h>
#include <TMath.h>

#include <const.h>

#include <charged_mom.h>
#include <neutral_mom.h>
#include <kloe_class.h>
#include <neu_triangle.h>

namespace KLOE
{
  class KinFit : public virtual pm00
  {
  public:
    /*Specific virtual functions for constraints*/

    // Four momentum conservation in LAB
    virtual Double_t FourMomConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t EnergyConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t PxConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t PyConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t PzConsvLAB(Double_t *, Double_t *) = 0;

    // Energy conservation in CM (trilateration only)
    virtual Double_t EnergyConsvCM(Double_t *, Double_t *) = 0;

    // Conservation of invariant mass
    virtual Double_t MinvConsv(Double_t *, Double_t *) = 0;
    virtual Double_t MinvConsvNeuKaon(Double_t *, Double_t *) = 0;
    virtual Double_t MinvConsvChKaon(Double_t *, Double_t *) = 0;
    virtual Double_t MinvConsvOmega(Double_t *, Double_t *) = 0;

    // Gamma path of flight from IP Conservation
    virtual Double_t PhotonPathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t Photon1PathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t Photon2PathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t Photon3PathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t Photon4PathConsvLAB(Double_t *, Double_t *) = 0;

    // Neutral kaon path of flight from IP Conservation
    virtual Double_t NeutralPathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t NeutralXPathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t NeutralYPathConsvLAB(Double_t *, Double_t *) = 0;
    virtual Double_t NeutralZPathConsvLAB(Double_t *, Double_t *) = 0;

    virtual void SetParameters(Double_t *) = 0;
    virtual void ResetParameters() = 0;
    virtual void IntermediateReconstruction() = 0;

    virtual ~KinFit() = default;
  };
}

#endif // !KINFIT_H
