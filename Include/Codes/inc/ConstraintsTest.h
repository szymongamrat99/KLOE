#ifndef CONSTRAINTSTEST_H
#define CONSTRAINTSTEST_H

#include <KinFit.h>

#include <ErrorLogs.h>

/*List of variables (order is the following)*/
/*
  For 4 gamma decay, without charged decay:
  4 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 4 gamma decay, with charged decay:
  2 x 3 (Curv_pi, Phiv_pi, Cotv_pi)
  4 x 5 (Xcl, Ycl, Zcl, Tcl, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, sqrt(S))

  Total: 2*3 + 4*5 + 1*4 = 6 + 20 + 4 = 30 parameters

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
   * @class ConstraintsTest
   * @brief Auxiliary class with the constraints for \omega\pi^{0} fitting
   */
  class ConstraintsTest : public KinFit, public ChargedVtxRec<Float_t, Int_t>
  {
  private:
    /**
     * @brief Four momentum conservation general method
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t FourMomConsvLAB(Double_t *x, Double_t *p) override;

    /**
     * @brief Photon path conservation general method
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t PhotonPathConsvLAB(Double_t *x, Double_t *p) override;

    /**
     * @brief Invariant mass conservation for \pi^{0} meson
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t MinvConsv(Double_t *x, Double_t *p) override;

    // Fictitious overriders
    /** Fictitious overrider of virtual method - do not use*/
    Double_t EnergyConsvCM(Double_t *x, Double_t *p) override { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralPathConsvLAB(Double_t *x, Double_t *p) override { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralXPathConsvLAB(Double_t *x, Double_t *p) override { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralYPathConsvLAB(Double_t *x, Double_t *p) override { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralZPathConsvLAB(Double_t *x, Double_t *p) override { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t MinvConsvOmega(Double_t *x, Double_t *p) override { return 0; };

    Int_t
        _chosen4MomComponent, /*!< Component of a 4-momentum to choose from FourMomConsvLAB*/
        _chosenPhoton;        /*!< Index of a photon to choose from PhotonPathConsvLAB*/

    std::string
        _minvMode;

  protected:
    void SetParameters(Double_t *p) override;
    void ResetParameters() override;
    void IntermediateReconstruction() override;

  public:
    /* Specific physical Constraints for Test-pi0 hypothesis */

    /**
     * @brief A method used to pair the photons in an event. Needed to get omega and neutral pions' parameters.
     */
    void PhotonPairing();

    // Four momentum conservation in LAB
    /**
     * @brief Conservation of total momentum's x-component
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t PxConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosen4MomComponent = 0;
      return FourMomConsvLAB(x, p);
    };
    /**
     * @brief Conservation of total momentum's y-component
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t PyConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosen4MomComponent = 1;
      return FourMomConsvLAB(x, p);
    };
    /**
     * @brief Conservation of total momentum's z-component
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t PzConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosen4MomComponent = 2;
      return FourMomConsvLAB(x, p);
    };
    /**
     * @brief Conservation of total energy
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t EnergyConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosen4MomComponent = 3;
      return FourMomConsvLAB(x, p);
    };

    // Gamma path of flight from IP Conservation
    /**
     * @brief Conservation of photon 1 path from IP
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t Photon1PathConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosenPhoton = 0;
      return PhotonPathConsvLAB(x, p);
    };
    /**
     * @brief Conservation of photon 2 path from IP
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t Photon2PathConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosenPhoton = 1;
      return PhotonPathConsvLAB(x, p);
    };
    /**
     * @brief Conservation of photon 3 path from IP
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t Photon3PathConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosenPhoton = 2;
      return PhotonPathConsvLAB(x, p);
    };
    /**
     * @brief Conservation of photon 4 path from IP
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t Photon4PathConsvLAB(Double_t *x, Double_t *p) override
    {
      _chosenPhoton = 3;
      return PhotonPathConsvLAB(x, p);
    };

    Double_t MinvConsvNeuKaon(Double_t *x, Double_t *p) override
    {
      _minvMode = "neutral";
      return MinvConsv(x, p);
    };

    Double_t MinvConsvChKaon(Double_t *x, Double_t *p) override
    {
      _minvMode = "charged";
      return MinvConsv(x, p);
    };
  };
}
#endif
