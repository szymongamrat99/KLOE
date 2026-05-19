#ifndef CHARGED_MOM_H
#define CHARGED_MOM_H

#include <const.h>
#include <TMath.h>

#include <fort_common.h>
#include <kloe_class.h>
#include <MomentumSmearing.h>

namespace KLOE
{
  // Template for function definition
  template <typename F = Double_t, typename T = Int_t>
  class ChargedVtxRec : public virtual pm00
  {
  private:
    Int_t
        &_nv,
        &_ntv,
        def = 0,
        &_mode;

    T
        *_iv;

    F
        *_CurV = nullptr,
        *_PhiV = nullptr,
        *_CotV = nullptr,
        *_PxTv = nullptr,
        *_PyTv = nullptr,
        *_PzTv = nullptr,
        *_xv = nullptr,
        *_yv = nullptr,
        *_zv = nullptr,
        *_IP = nullptr,
        _dist = 0;

    std::array<const std::string, 6>
        logMessages =
            {
                "Error in the KaonMomFromBoost - discriminant is negative.",
                "Error in the findKchRec - no vertex with opposite tracks found.",
                "Error in the findKSLRec - no vertex with opposite tracks found.",
                "Error in the findKSLRec - null pointer encountered.",
                "Error in the IPBoostCorr - denominator is zero.",
                "Error in the findKClosest - no vertex with opposite tracks found."};

    ErrorHandling::ErrorLogs &_logger;

    void charged_mom(Int_t &i, F *mom_vec, ErrorHandling::ErrorLogs &logger);
    static Double_t calculateTrackEnergy(std::vector<F> &trk, T mode = 1);

    ErrorHandling::ErrorCodes setTrackMom(F *mom_vec1Tmp, F *mom_vec2Tmp, Int_t j1, Int_t j2)
    {
      if (_PxTv != nullptr && _PyTv != nullptr && _PzTv != nullptr)
      {
        mom_vec1Tmp[0] = _PxTv[j1];
        mom_vec1Tmp[1] = _PyTv[j1];
        mom_vec1Tmp[2] = _PzTv[j1];
        mom_vec2Tmp[0] = _PxTv[j2];
        mom_vec2Tmp[1] = _PyTv[j2];
        mom_vec2Tmp[2] = _PzTv[j2];
      }
      else if (_CurV != nullptr && _CotV != nullptr && _PhiV != nullptr)
      {
        ChargedVtxRec::charged_mom(j1, mom_vec1Tmp, _logger);
        ChargedVtxRec::charged_mom(j2, mom_vec2Tmp, _logger);
      }
      else
      {
        ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NULL_POINTER;
        return err;
      }
      return ErrorHandling::ErrorCodes::NO_ERROR;
    };

    void makeMomentumSmearing(T mcflag, Bool_t smearingFlag, TMatrixT<Double_t> covMatrix, TVectorT<Double_t> &momVecMC, TVectorT<Double_t> &momVecSmeared, F *mom_vec1Tmp, F *mom_vec2Tmp)
    {
      KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, covMatrix);

      if (mcflag == 1 && smearingFlag == 1)
      {
        momVecMC[0] = mom_vec1Tmp[0];
        momVecMC[1] = mom_vec1Tmp[1];
        momVecMC[2] = mom_vec1Tmp[2];
        momVecMC[3] = mom_vec2Tmp[0];
        momVecMC[4] = mom_vec2Tmp[1];
        momVecMC[5] = mom_vec2Tmp[2];

        CovMatrixCalcObj.SetMCVector(momVecMC);
        CovMatrixCalcObj.SmearMomentum();
        CovMatrixCalcObj.GetSmearedMomentum(momVecSmeared);
      }
      else
      {
        momVecSmeared[0] = mom_vec1Tmp[0];
        momVecSmeared[1] = mom_vec1Tmp[1];
        momVecSmeared[2] = mom_vec1Tmp[2];
        momVecSmeared[3] = mom_vec2Tmp[0];
        momVecSmeared[4] = mom_vec2Tmp[1];
        momVecSmeared[5] = mom_vec2Tmp[2];
      }
    }

    void bestTrackCombination(std::vector<F> &trk1, std::vector<F> &trk2, std::vector<F> &Kchrec, ErrorHandling::ErrorLogs &logger, T mode = 1) const
    {
      Double_t
          diff[2] = {std::numeric_limits<Double_t>::max(), std::numeric_limits<Double_t>::max()},
          energyPion[2] = {0, 0},
          energyLepton[2] = {0, 0};

      std::vector<F> trk[2] = {trk1, trk2};

      std::vector<F> KchrecTmp[2] = {std::vector<F>(6, 0), std::vector<F>(6, 0)};

      energyPion[0] = calculateTrackEnergy(trk[0]);
      energyPion[1] = calculateTrackEnergy(trk[1]);

      if (mode != 1)
      {
        energyLepton[0] = calculateTrackEnergy(trk[0], mode);
        energyLepton[1] = calculateTrackEnergy(trk[1], mode);

        KchrecTmp[0][0] = trk[0][0] + trk[1][0];
        KchrecTmp[0][1] = trk[0][1] + trk[1][1];
        KchrecTmp[0][2] = trk[0][2] + trk[1][2];
        KchrecTmp[0][3] = energyPion[0] + energyLepton[1];
        KchrecTmp[0][4] = std::sqrt(std::pow(KchrecTmp[0][0], 2) + std::pow(KchrecTmp[0][1], 2) + std::pow(KchrecTmp[0][2], 2));
        KchrecTmp[0][5] = std::sqrt(std::pow(KchrecTmp[0][3], 2) - std::pow(KchrecTmp[0][4], 2));

        KchrecTmp[1][0] = trk[0][0] + trk[1][0];
        KchrecTmp[1][1] = trk[0][1] + trk[1][1];
        KchrecTmp[1][2] = trk[0][2] + trk[1][2];
        KchrecTmp[1][3] = energyLepton[0] + energyPion[1];
        KchrecTmp[1][4] = std::sqrt(std::pow(KchrecTmp[1][0], 2) + std::pow(KchrecTmp[1][1], 2) + std::pow(KchrecTmp[1][2], 2));
        KchrecTmp[1][5] = std::sqrt(std::pow(KchrecTmp[1][3], 2) - std::pow(KchrecTmp[1][4], 2));

        diff[0] = std::abs(KchrecTmp[0][5] - PhysicsConstants::mK0);
        diff[1] = std::abs(KchrecTmp[1][5] - PhysicsConstants::mK0);

        if (diff[0] < diff[1])
        {
          trk1[3] = energyPion[0];
          trk2[3] = energyLepton[1];

          Kchrec[0] = KchrecTmp[0][0];
          Kchrec[1] = KchrecTmp[0][1];
          Kchrec[2] = KchrecTmp[0][2];
          Kchrec[3] = KchrecTmp[0][3];
          Kchrec[4] = KchrecTmp[0][4];
          Kchrec[5] = KchrecTmp[0][5];
        }
        else
        {
          trk1[3] = energyLepton[0];
          trk2[3] = energyPion[1];

          Kchrec[0] = KchrecTmp[1][0];
          Kchrec[1] = KchrecTmp[1][1];
          Kchrec[2] = KchrecTmp[1][2];
          Kchrec[3] = KchrecTmp[1][3];
          Kchrec[4] = KchrecTmp[1][4];
          Kchrec[5] = KchrecTmp[1][5];
        }
      }
      else
      {
        trk1[3] = energyPion[0];
        trk2[3] = energyPion[1];

        Kchrec[0] = trk1[0] + trk2[0];
        Kchrec[1] = trk1[1] + trk2[1];
        Kchrec[2] = trk1[2] + trk2[2];
        Kchrec[3] = energyPion[0] + energyPion[1];
        Kchrec[4] = std::sqrt(std::pow(Kchrec[0], 2) + std::pow(Kchrec[1], 2) + std::pow(Kchrec[2], 2));
        Kchrec[5] = std::sqrt(std::pow(Kchrec[3], 2) - std::pow(Kchrec[4], 2));
      }
    };


  public:
    ChargedVtxRec(Int_t &nv, Int_t &ntv, T *ivOld, F *IP, F *CurV, F *PhiV, F *CotV, F *xvOld, F *yvOld, F *zvOld, Int_t &mode, ErrorHandling::ErrorLogs &logger);
    ChargedVtxRec(Int_t &nv, Int_t &ntv, T *ivOld, F *IP, F *CurV, F *PxTv, F *PyTv, F *PzTv, F *xvOld, F *yvOld, F *zvOld, Int_t &mode, ErrorHandling::ErrorLogs &logger);
    ChargedVtxRec(ErrorHandling::ErrorLogs &logger);

    static void charged_mom(F CurvOld, F PhivOld, F CotvOld, F *mom_vec, Int_t mode, ErrorHandling::ErrorLogs &logger);

    ErrorHandling::ErrorCodes findKchRec(T mcflag, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKchRec(T mcflag, Bool_t smearingFlag, TMatrixT<Double_t> &covMatrix, std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger, T mode = 1);
    ErrorHandling::ErrorCodes findKSLRec(Int_t kaonFlag, Int_t KSvtx, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKSLRec(Int_t kaonFlag, Int_t KSvtx, std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKClosestRec(F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKClosestRec(std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger);

    Int_t KaonMomFromBoost(F *pKaon, F *pboost, F *pKaonBoost) const;
    Int_t KaonMomFromBoost(std::vector<F> &pKaon, F *pboost, std::vector<F> &pKaonBoost) const;

    Int_t IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, F *int_point) const;
    Int_t IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, std::vector<F> &int_point) const;

    void setMode(Int_t mode) { _mode = mode; }

    void calculateTrackTOF(std::vector<T> &vtaken, std::vector<T> &Asstrk, std::vector<T> &Asscl, std::vector<F> &Assleng, std::vector<F> &Enecl, std::vector<F> &Xcl, std::vector<F> &Ycl, std::vector<F> &Zcl, std::vector<F> &Tcl, std::vector<F> trk[2], std::vector<F> trkCluster[2], std::vector<F> trkDT[2]) const;
  };

}
#endif