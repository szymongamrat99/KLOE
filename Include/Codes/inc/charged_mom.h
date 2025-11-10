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
  template <typename F = Float_t, typename T = Int_t> 
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

    void charged_mom(Int_t &i, F *mom_vec);

  public:
    ChargedVtxRec(Int_t &nv, Int_t &ntv, T *ivOld, F *IP, F *CurV, F *PhiV, F *CotV, F *xvOld, F *yvOld, F *zvOld, Int_t &mode);
    ChargedVtxRec(Int_t &nv, Int_t &ntv, T *ivOld, F *IP, F *CurV, F *PxTv, F *PyTv, F *PzTv, F *xvOld, F *yvOld, F *zvOld, Int_t &mode);
    ChargedVtxRec();

    static void charged_mom(F CurvOld, F PhivOld, F CotvOld, F *mom_vec, Int_t mode);

    ErrorHandling::ErrorCodes findKchRec(T mcflag, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKchRec(T mcflag, Bool_t smearingFlag, TMatrixT<Double_t> covMatrix, std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKSLRec(Int_t kaonFlag, Int_t KSvtx, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKSLRec(Int_t kaonFlag, Int_t KSvtx, std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKClosestRec(F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger);
    ErrorHandling::ErrorCodes findKClosestRec(std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger);

    Int_t KaonMomFromBoost(F *pKaon, F *pboost, F *pKaonBoost) const;
    Int_t KaonMomFromBoost(std::vector<F> &pKaon, F *pboost, std::vector<F> &pKaonBoost) const;

    Int_t IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, F *int_point) const;
    Int_t IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, std::vector<F> &int_point) const;

  };

}
#endif