#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>
#include <chrono>

#include "../inc/trilaterationKinFit.h"

namespace KLOE
{
  TrilaterationReconstructionKinFit::TrilaterationReconstructionKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, Int_t jmin, Int_t jmax, ErrorHandling::ErrorLogs &logger) : KinFitter("Trilateration", N_free, N_const, M, 0, loopcount, chiSqrStep, logger), _jmin(jmin), _jmax(jmax)
  {
    _chargedVtxRec = new ChargedVtxRec<Float_t, Int_t>(logger);

    _V.ResizeTo(N_free + N_const, N_free + N_const),
        _D.ResizeTo(M, N_free + N_const),
        _D_T.ResizeTo(N_free + N_const, M),
        _V_final.ResizeTo(N_free + N_const, N_free + N_const),
        _V_aux.ResizeTo(N_free + N_const, N_free + N_const),
        _V_min.ResizeTo(N_free + N_const, N_free + N_const),
        _Aux.ResizeTo(M, M),
        _V_invert.ResizeTo(N_free, N_free),
        _V_init.ResizeTo(N_free + N_const, N_free + N_const);

    _X.ResizeTo(N_free + N_const);
    _C.ResizeTo(M);
    _X_final.ResizeTo(N_free + N_const);
    _L.ResizeTo(M);
    _CORR.ResizeTo(N_free + N_const);
    _X_init.ResizeTo(N_free + N_const);
    _X_min.ResizeTo(N_free + N_const);
    _C_min.ResizeTo(M);
    _L_min.ResizeTo(M);
    _C_aux.ResizeTo(M);
    _L_aux.ResizeTo(M);
    _X_init_min.ResizeTo(N_free + N_const);
    _X_init_aux.ResizeTo(N_free + N_const);

    _Param.resize(N_free + N_const);
    _Errors.resize(N_free + N_const);

    _recMode = pm00::StringToTrilaterationCode(_config.getProperty<std::string>("flags.trilaterationReconstructionMode"));

    if (_recMode == TrilaterationCode::TWO_PI0)
    {
      _selected.resize(4);
      for (Int_t i = 1; i <= 4; i++)
        _selected[i - 1] = i;
      _ind_gam.resize(4);

      for (Int_t i = 0; i < 4; i++)
        _gamma_mom_final[i].resize(8);

      _neu_vtx_min.resize(4);
      _iptri_kinfit.resize(3);
      _fourKnetri_kinfit.resize(9);
      _g4takentri_kinfit.resize(4);

      KinFitter::ConstraintSet({"EnergyConsvCM",
                                "MinvConsv",
                                "NeutralXPathConsvLAB",
                                "NeutralYPathConsvLAB",
                                "NeutralZPathConsvLAB"});
    }
    else if (_recMode == TrilaterationCode::THREE_PI0)
    {
      _selected.resize(6);
      for (Int_t i = 1; i <= 6; i++)
        _selected[i - 1] = i;
      _ind_gam.resize(6);
    }
    else
    {
      std::cout << "Wrong trilateration reconstruction mode!" << std::endl;
      exit(EXIT_FAILURE);
    }

    gErrorIgnoreLevel = 6001;
  }

  TrilaterationReconstructionKinFit::~TrilaterationReconstructionKinFit()
  {
  }

  ErrorHandling::ErrorCodes TrilaterationReconstructionKinFit::Reconstruct()
  {
    Float_t
        CHISQRTMP = 999.,
        FUNVALTMP = 999999.,
        Tcorr = 0;

    Float_t value_tmp[2] = {999999., 999999.};

    Bool_t
        clusterEnergy = true,
        cond_clus[4] = {false, false, false, false},
        cond_time_clus[2] = {false, false};

    Float_t
        gamma_mom_tmp[4][8] = {},
        fourKnetri_tmp[2][10] = {},
        kaon_vel_tmp[2] = {},
        y_axis[3] = {},
        ip_tmp[2][3] = {};

    Float_t distance[4] = {0.};

    Float_t neu_vtx[2][4] = {};
    Float_t dist_tmp[2] = {0.};

    _CHISQRMIN = 999999.;
    _isConverged = 0;

    Int_t chosenSolution = -1;

    Int_t neucluwrong = 0;
    std::vector<Int_t> neucluwrongInd;

    for (Int_t i = 0; i < _NeuClusters.size(); i++)
    {
      if (_cluster[4][_NeuClusters[i] - 1] < MIN_CLU_ENE)
      {
        neucluwrong++;
        neucluwrongInd.push_back(i);
      }
    }

    if (_NeuClusters.size() - neucluwrong < NCLMIN)
      return ErrorHandling::ErrorCodes::LESS_THAN_FOUR_CLUSTERS_WITH_GOOD_ENERGY;

    for (Int_t j1 = 0; j1 < _NeuClusters.size() - 3; j1++)
      for (Int_t j2 = j1 + 1; j2 < _NeuClusters.size() - 2; j2++)
        for (Int_t j3 = j2 + 1; j3 < _NeuClusters.size() - 1; j3++)
          for (Int_t j4 = j3 + 1; j4 < _NeuClusters.size(); j4++)
          {
            _ind_gam[0] = j1;
            _ind_gam[1] = j2;
            _ind_gam[2] = j3;
            _ind_gam[3] = j4;

            Bool_t hasBadCluster = false;
            for (const auto &idx : neucluwrongInd)
            {
              if (idx == j1 || idx == j2 || idx == j3 || idx == j4)
              {
                hasBadCluster = true;
                break;
              }
            }

            if (hasBadCluster)
              continue;

            for (Int_t k = 0; k < 4; k++)
            {
              _R.SetClu(k, _cluster[0][_NeuClusters[_ind_gam[k]] - 1],
                        _cluster[1][_NeuClusters[_ind_gam[k]] - 1],
                        _cluster[2][_NeuClusters[_ind_gam[k]] - 1],
                        _cluster[3][_NeuClusters[_ind_gam[k]] - 1],
                        _cluster[4][_NeuClusters[_ind_gam[k]] - 1]);
            }

            _S = _R.MySolve(_selected.data());

            for (Int_t k = 0; k < 4; k++)
            {
              cond_clus[k] =
                  _cluster[0][_NeuClusters[_ind_gam[k]] - 1] != 0 &&
                  _cluster[1][_NeuClusters[_ind_gam[k]] - 1] != 0 &&
                  _cluster[2][_NeuClusters[_ind_gam[k]] - 1] != 0 &&
                  _cluster[3][_NeuClusters[_ind_gam[k]] - 1] != 0;
            }

            Bool_t cond_tot = cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && clusterEnergy;

            if ((_S.error[0] && _S.error[1]) || !cond_tot)
              continue;

            for (Int_t k = 0; k < 4; k++)
            {
              _Param[k * 5] = _cluster[0][_NeuClusters[_ind_gam[k]] - 1];
              _Param[k * 5 + 1] = _cluster[1][_NeuClusters[_ind_gam[k]] - 1];
              _Param[k * 5 + 2] = _cluster[2][_NeuClusters[_ind_gam[k]] - 1];
              _Param[k * 5 + 3] = _cluster[3][_NeuClusters[_ind_gam[k]] - 1];
              _Param[k * 5 + 4] = _cluster[4][_NeuClusters[_ind_gam[k]] - 1];

              _Errors[k * 5] = clu_x_error(_Param[k * 5], _Param[k * 5 + 1], _Param[k * 5 + 2], _Param[k * 5 + 4]);
              _Errors[k * 5 + 1] = clu_y_error(_Param[k * 5], _Param[k * 5 + 1], _Param[k * 5 + 2], _Param[k * 5 + 4]);
              _Errors[k * 5 + 2] = clu_z_error(_Param[k * 5], _Param[k * 5 + 1], _Param[k * 5 + 2], _Param[k * 5 + 4]);
              // cm
              _Errors[k * 5 + 3] = clu_time_error(_Param[k * 5 + 4]); // ns
              _Errors[k * 5 + 4] = clu_ene_error(_Param[k * 5 + 4]);  // MeV

              _Param[20 + k] = _bhabha_mom[k];
              _Errors[20 + k] = _bhabha_mom_err[k];

              if (k < 3)
              {
                _Param[24 + k] = _bhabha_vtx[k];
                _Errors[24 + k] = 0.;
              }
            }

            for (Int_t k1 = _jmin; k1 <= _jmax; k1++)
            {
              KinFitter::ParameterInitialization(_Param.data(), _Errors.data());

              Tcorr = k1 * T0;

              CHISQRTMP = KinFitter::FitFunction(Tcorr);


              KinFitter::GetResults(_X_min, _V_min, _X_init_min, _V_init, _ipFitTri, _photonFitTri, _KnerecFitTri, _PhiFitTri);

              Bool_t hasBetterChi2 = (CHISQRTMP < _CHISQRMIN);

              Bool_t condTime = _KnerecFitTri[9] < _X_min(3) &&
                                _KnerecFitTri[9] < _X_min(8) &&
                                _KnerecFitTri[9] < _X_min(13) &&
                                _KnerecFitTri[9] < _X_min(18);

              if ((hasBetterChi2 && condTime) || (hasBetterChi2 && _isConverged == 0))
              {
                _isConverged = 1;
                _FUNVALMIN = FUNVALTMP;
                _CHISQRMIN = CHISQRTMP;

                _Chi2TriKinFit = _CHISQRMIN;

                _g4takentri_kinfit = _ind_gam;

                _neu_vtx_min[0] = _KnerecFitTri[6];
                _neu_vtx_min[1] = _KnerecFitTri[7];
                _neu_vtx_min[2] = _KnerecFitTri[8];
                _neu_vtx_min[3] = _KnerecFitTri[9];

                _gamma_mom_final[0] = _photonFitTri[0];
                _gamma_mom_final[1] = _photonFitTri[1];
                _gamma_mom_final[2] = _photonFitTri[2];
                _gamma_mom_final[3] = _photonFitTri[3];

                _fourKnetri_kinfit = _KnerecFitTri;

                _iptri_kinfit = _ipFitTri;

                _bunchnum = k1;
              }
            }
          }

    if (_isConverged == 1)
    {
      return ErrorHandling::ErrorCodes::NO_ERROR;
    }
    else
    {
      _neu_vtx_min[0] = 0.;
      _neu_vtx_min[1] = 999;
      _neu_vtx_min[2] = 999;
      _neu_vtx_min[3] = 999;

      for (Int_t l = 0; l < 4; l++)
      {
        _gamma_mom_final[l][0] = 999.;
        _gamma_mom_final[l][1] = 999.;
        _gamma_mom_final[l][2] = 999.;
        _gamma_mom_final[l][3] = 999.;
        _gamma_mom_final[l][4] = 999.;
        _gamma_mom_final[l][5] = 999.;
        _gamma_mom_final[l][6] = 999.;
        _gamma_mom_final[l][7] = 999.;
      }

      _fourKnetri_kinfit[0] = 999.;
      _fourKnetri_kinfit[1] = 999.;
      _fourKnetri_kinfit[2] = 999.;
      _fourKnetri_kinfit[3] = 999.;
      _fourKnetri_kinfit[4] = 999.;
      _fourKnetri_kinfit[5] = 999.;
      _fourKnetri_kinfit[6] = 999.;
      _fourKnetri_kinfit[7] = 999.;
      _fourKnetri_kinfit[8] = 999.;
      _fourKnetri_kinfit[9] = 999.;

      _iptri_kinfit[0] = 999.;
      _iptri_kinfit[1] = 999.;
      _iptri_kinfit[2] = 999.;

      _bunchnum = 999;

      _Chi2TriKinFit = 999.;

      return ErrorHandling::ErrorCodes::TRILATERATION_KIN_FIT;
    }
  }
}