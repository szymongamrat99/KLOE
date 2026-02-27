#pragma once

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <ConfigManager.h>
#include <charged_mom.h>
#include <neutral_mom.h>

#include <KinFitter.h>
#include <kloe_class.h>

namespace KLOE
{
  class OmegaKinFit : public KinFitter
  {
  private:
    TVectorD
        _X,
        _C,
        _X_final,
        _L,
        _CORR,
        _X_init,
        _X_min,
        _C_min,
        _L_min,
        _C_aux,
        _L_aux,
        _X_init_min,
        _X_init_aux;

    TMatrixD
        _V,
        _D,
        _D_T,
        _V_final,
        _V_aux,
        _V_min,
        _Aux,
        _V_invert,
        _V_init;

    Double_t
        _min_value_def,
        _CHISQRMIN,
        _Chi2TriKinFit;

    Bool_t
        _isConverged;

    std::array<std::vector<Float_t>, 2>
        _trkFit,
        _trackParameters,
        _trackParametersErr,
        _Pi0OmegaFit;

    std::array<std::vector<Float_t>, 4>
        _cluster,
        _photonFit;

    std::vector<Float_t>
        _omegaVtx,
        _omegaVtxErr,
        _bhabha_mom,
        _bhabha_mom_err,
        _bhabha_vtx,
        _bhabhaVtxErr,
        _Param,
        _Errors,
        _OmegaFit,
        _ipFit,
        _PhiMomFit;

    Int_t
        _offset;

    ConfigManager
        &_config = ConfigManager::getInstance();

  public:
    OmegaKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, ErrorHandling::ErrorLogs &logger);
    ~OmegaKinFit();

    void SetParameters(const std::vector<Float_t> trackParameters[2], const std::vector<Float_t> trackParametersErr[2], const std::vector<Float_t> cluster[4], const std::vector<Float_t> bhabha_mom, const std::vector<Float_t> bhabha_mom_err, const std::vector<Float_t> omegaVtx, const std::vector<Float_t> omegaVtxErr, const std::vector<Float_t> bhabha_vtx, const std::vector<Float_t> bhabhaVtxErr)
    {
      for (Int_t i = 0; i < 2; i++)
      {
        _trackParameters[i] = trackParameters[i];
        _trackParametersErr[i] = trackParametersErr[i];
      }

      for (Int_t i = 0; i < 4; i++)
        _cluster[i].assign(cluster[i].begin(), cluster[i].end());


      _bhabha_mom = bhabha_mom;
      _bhabha_mom_err = bhabha_mom_err;

      _omegaVtx = omegaVtx;
      _omegaVtxErr = omegaVtxErr;

      _bhabha_vtx = bhabha_vtx;
      _bhabhaVtxErr = bhabhaVtxErr;
    }

    void GetResults(std::vector<Float_t> &Param, std::vector<Float_t> &Errors, std::vector<Float_t> &ParamFit, std::vector<Float_t> &ErrorsFit, std::vector<Float_t> trkFit[2], std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> &OmegaFit, std::vector<Float_t> Pi0OmegaFit[2], std::vector<Float_t> &PhiMomFit, Float_t &Chi2OmegaKinFit, std::vector<Float_t> &pulls)
    {
      Param = _Param;
      Errors = _Errors;
      ParamFit.resize(_X_min.GetNrows());
      ErrorsFit.resize(_X_min.GetNrows());

      for (Int_t i = 0; i < _X_min.GetNrows(); i++)
      {
        ParamFit[i] = _X_min[i];
        ErrorsFit[i] = sqrt(_V_min[i][i]);
      }

      for (Int_t i = 0; i < 2; i++)
        trkFit[i] = _trkFit[i];

      ipFit = _ipFit;

      for (Int_t i = 0; i < 4; i++)
        photonFit[i] = _photonFit[i];

      Pi0OmegaFit[0] = _Pi0OmegaFit[0];
      Pi0OmegaFit[1] = _Pi0OmegaFit[1];

      OmegaFit = _OmegaFit;

      PhiMomFit = _PhiMomFit;

      Chi2OmegaKinFit = _CHISQRMIN;

      for (Int_t i = 0; i < _X_min.GetNrows(); i++)
      {
        pulls.push_back((_X_init_min[i] - _X_min[i]) / sqrt(_V_init[i][i] - _V_min[i][i]));
      }
    }

    ErrorHandling::ErrorCodes Reconstruct();
  };

} // namespace KLOE
