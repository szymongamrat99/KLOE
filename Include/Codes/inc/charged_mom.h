#ifndef CHARGED_MOM_H
#define CHARGED_MOM_H

#include <const.h>
#include <TMath.h>

#include <fort_common.h>
#include <kloe_class.h>

namespace KLOE
{
  // Template for function definition
  template <typename F = Float_t> 
  class ChargedVtxRec : public pm00
  {
  private:
    Int_t
        &_nv,
        &_ntv,
        def = 0,
        &_mode;

    UChar_t
          *_iv;

    F
        *_CurV,
        *_PhiV,
        *_CotV,
        *_xv,
        *_yv,
        *_zv,
        *_IP,
        _dist;

    void charged_mom(Int_t &i, F *mom_vec);

  public:
    ChargedVtxRec(Int_t &nv, Int_t &ntv, UChar_t *iv, F *IP, F *CurV, F *PhiV, F *CotV, F *xv, F *yv, F *zv, Int_t &mode);
    ChargedVtxRec();

    void charged_mom(F Curv, F Phiv, F Cotv, F *mom_vec, Int_t mode);

    void findKchRec(F *KchRec, F *trk1, F *trk2, Int_t *vtaken, Int_t &errFlag);
    void findKSLRec(Int_t kaonFlag, Int_t KSvtx, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, Int_t &errFlag);

    void findKClosestRec(F *KchRec, F *trk1, F *trk2, Int_t *vtaken, Int_t &errFlag);

    Int_t KaonMomFromBoost(F *pKaon, F *pboost, F *pKaonBoost);

    Int_t IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, F *int_point);
  };

}
#endif