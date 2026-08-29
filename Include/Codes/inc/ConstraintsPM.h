#ifndef CONSTRAINTSPM_H
#define CONSTRAINTSPM_H

#include <KinFit.h>

#include <ErrorLogs.h>

/*List of variables (order is the following)*/
/*
  For 4 gamma decay, without charged decay:
  4 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, std::sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 4 gamma decay, with charged decay:
  2 x 3 (Curv_pi, Phiv_pi, Cotv_pi)
  4 x 5 (Xcl, Ycl, Zcl, Tcl, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, std::sqrt(S))

  Total: 2*3 + 4*5 + 1*4 = 6 + 20 + 4 = 30 parameters

  IP determined during fit (using direction of Kch momentum)
  Charged kaon 4-mom determined during fit (using the boost method as well)
  Neutral kaon 4-vec determined during fit
  Neutral kaon dependent of Kch direction and energy

--------------------------------------------------------------------------------

  For 6 gamma decay, without charged decay:
  6 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, std::sqrt(S))

  IP determined during fit (using direction of Kne momentum from clusters)
  Neutral kaon 4-vec determined during fit

--------------------------------------------------------------------------------

  For 6 gamma decay, with charged decay:
  2 x 3 (CurvOld, PhivOld, CotvOld)
  6 x 5 (Xcl, Ycl, Zcl, TclOld, EneCl)
  1 x 4 (Px_Phi, Py_Phi, Pz_Phi, std::sqrt(S))

  IP determined during fit (using direction of Kch momentum)
  Charged kaon 4-mom determined during fit (using the boost method as well)
  Neutral kaon 4-vec determined during fit
  Neutral kaon dependent of Kch direction and energy

--------------------------------------------------------------------------------

*/

namespace KLOE
{
  /**
   * @class ConstraintsPM
   * @brief Auxiliary class with the constraints for \omega\pi^{0} fitting
   */
  class ConstraintsPM : public ChargedVtxRec<Double_t, Int_t>
  {
  private:
    /**
     * @brief Energy conservation in the center-of-mass frame
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t EnergyConsvCM(Double_t *x, Double_t *p);

    /**
     * @brief Invariant mass conservation for \pi^{0} meson
     * @param x pointer to the variables table - artificial variable to use TF1
     * @param p pointer to the parameters table - physical variables to fit
     * @returns double precision value of a constraint
     */
    Double_t MinvConsv(Double_t *x, Double_t *p);

    // Fictitious overriders
    /** Fictitious overrider of virtual method - do not use*/
    Double_t PhotonPathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralPathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralXPathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralYPathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t NeutralZPathConsvLAB(Double_t *x, Double_t *p) { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t MinvConsvOmega(Double_t *x, Double_t *p) { return 0; };
    /** Fictitious overrider of virtual method - do not use*/
    Double_t FourMomConsvLAB(Double_t *x, Double_t *p) { return 0; };

    Int_t
        _chosen4MomComponent, /*!< Component of a 4-momentum to choose from FourMomConsvLAB*/
        _chosenPhoton;        /*!< Index of a photon to choose from PhotonPathConsvLAB*/

    std::string
        _minvMode;

    ErrorHandling::ErrorLogs &_logger;

  protected:
    void SetParameters(Double_t *p);
    void ResetParameters();
    void IntermediateReconstruction();


  public:
    /* Specific physical Constraints for Signal-pi0 hypothesis */

    ConstraintsPM(ErrorHandling::ErrorLogs &logger) : ChargedVtxRec<Double_t, Int_t>(logger), _logger(logger) {};

    std::array<neutralParticle, 4> fphoton;
    std::array<chargedParticle, 2> fpionCh;
    kaonNeutral fKchrec, fKchboost, fKchrecCMPhi;
    phiMeson fphi;
    std::vector<Double_t> fip;

    void IntermediateReconstruction(Double_t *p);

    /**
     * @brief A method used to pair the photons in an event. Needed to get omega and neutral pions' parameters.
     */
    void PhotonPairing();

    Double_t EnergyConsvCMChKaon(Double_t *x, Double_t *p)
    {
      return EnergyConsvCM(x, p);
    };

    Double_t MinvConsvChKaon(Double_t *x, Double_t *p)
    {
      return MinvConsv(x, p);
    };
  };
}
#endif
