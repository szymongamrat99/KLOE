
#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>
#include <boost/optional.hpp>

#include "../inc/trilateration.hpp"

ErrorHandling::ErrorCodes TrilaterationKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Float_t chiSqrStep, Int_t jmin, Int_t jmax, Int_t nclu, std::vector<Float_t> cluster[5], std::vector<Int_t> Asscl, std::vector<Float_t> bhabha_mom, std::vector<Float_t> bhabha_mom_err, std::vector<Float_t> bhabha_vtx, Int_t &bunchnum, std::vector<Float_t> &iptri_kinfit, std::vector<Int_t> &g4takentri_kinfit, std::vector<Float_t> gamma_mom_final[4], std::vector<Float_t> &fourKnetri_kinfit, std::vector<Float_t> &neu_vtx_min, Float_t &Chi2TriKinFit, ErrorHandling::ErrorLogs &logger)
{
  gErrorIgnoreLevel = 6001;

  std::unique_ptr<KLOE::ChargedVtxRec<>> eventAnalysis = std::make_unique<KLOE::ChargedVtxRec<>>(logger);

  std::vector<std::string> ConstSet = {
      "EnergyConsvCM",
      "MinvConsv",
      "NeutralXPathConsvLAB",
      "NeutralYPathConsvLAB",
      "NeutralZPathConsvLAB"};

  std::unique_ptr<KLOE::KinFitter> kinematicFitObj = std::make_unique<KLOE::KinFitter>("Trilateration", N_free, N_const, M, 0, loopcount, chiSqrStep, logger);
  kinematicFitObj->ConstraintSet(ConstSet);

  Double_t
      Param[N_free + N_const],
      Errors[N_free + N_const];

  Double_t det;
  Float_t CHISQR, CHISQRTMP, FUNVAL, FUNVALTMP, FUNVALMIN, Tcorr;
  Int_t fail;
  Bool_t cond_time_clus[2];

  Int_t selected[4] = {1, 2, 3, 4};

  Bool_t clusterEnergy, solError, cond_clus[4];
  Int_t ind_gam[4], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged, n_bunch;
  Float_t CHISQRMIN, min_value_def;
  Float_t gamma_mom_min[4][4], kaon_mom_min[4], kaon_vel[3], kaon_vel_tot, kaon_path_tot, bhabha_vtx_min[3], ip_min[3], time_diff[2][4], time_diff_fin, gamma_path[2][4], neu_vtx[2][4], dist_tmp[2], value[2];

  TMatrixD
      *V = new TMatrixD(N_free + N_const, N_free + N_const),
      *D = new TMatrixD(M, N_free + N_const),
      *D_T = new TMatrixD(N_free + N_const, M),
      *V_final = new TMatrixD(N_free + N_const, N_free + N_const),
      *V_aux = new TMatrixD(N_free + N_const, N_free + N_const),
      *V_min = new TMatrixD(N_free + N_const, N_free + N_const),
      *Aux = new TMatrixD(M, M),
      *V_invert = new TMatrixD(N_free, N_free),
      *V_init = new TMatrixD(N_free + N_const, N_free + N_const);

  TVectorD
      *X = new TVectorD(N_free + N_const),
      *C = new TVectorD(M),
      *X_final = new TVectorD(N_free + N_const),
      *L = new TVectorD(M),
      *CORR = new TVectorD(N_free + N_const),
      *X_init = new TVectorD(N_free + N_const),
      *X_min = new TVectorD(N_free + N_const),
      *C_min = new TVectorD(M),
      *L_min = new TVectorD(M),
      *C_aux = new TVectorD(M),
      *L_aux = new TVectorD(M),
      *X_init_min = new TVectorD(N_free + N_const),
      *X_init_aux = new TVectorD(N_free + N_const);

  min_value_def = 999999.;
  FUNVALMIN = 999999.;
  CHISQRMIN = 999999.;

  isConverged = 0;

  for (Int_t j1 = 0; j1 < Asscl.size() - 3; j1++)
    for (Int_t j2 = j1 + 1; j2 < Asscl.size() - 2; j2++)
      for (Int_t j3 = j2 + 1; j3 < Asscl.size() - 1; j3++)
        for (Int_t j4 = j3 + 1; j4 < Asscl.size(); j4++)
        {

          ind_gam[0] = j1;
          ind_gam[1] = j2;
          ind_gam[2] = j3;
          ind_gam[3] = j4;

          Reconstructor *R = new Reconstructor();
          Solution *S = new Solution();

          for (Int_t k = 0; k < 4; k++)
          {
            R->SetClu(k, cluster[0][Asscl[ind_gam[k]] - 1],
                      cluster[1][Asscl[ind_gam[k]] - 1],
                      cluster[2][Asscl[ind_gam[k]] - 1],
                      cluster[3][Asscl[ind_gam[k]] - 1],
                      cluster[4][Asscl[ind_gam[k]] - 1]);

            R->SetClu(4, 0., 0., 0., 0., 0.);
            R->SetClu(5, 0., 0., 0., 0., 0.);
          }

          *S = R->MySolve(selected);

          clusterEnergy = cluster[4][Asscl[ind_gam[0]] - 1] > KLOE::MIN_CLU_ENE &&
                          cluster[4][Asscl[ind_gam[1]] - 1] > KLOE::MIN_CLU_ENE &&
                          cluster[4][Asscl[ind_gam[2]] - 1] > KLOE::MIN_CLU_ENE &&
                          cluster[4][Asscl[ind_gam[3]] - 1] > KLOE::MIN_CLU_ENE;

          for (Int_t k = 0; k < 4; k++)
          {
            cond_clus[k] =
                cluster[3][Asscl[ind_gam[k]] - 1] > 0 &&
                cluster[0][Asscl[ind_gam[k]] - 1] != 0 &&
                cluster[1][Asscl[ind_gam[k]] - 1] != 0 &&
                cluster[2][Asscl[ind_gam[k]] - 1] != 0;

            if (k < 2)
              cond_time_clus[k] = S->sol[k][3] < cluster[3][Asscl[ind_gam[0]] - 1] &&
                                  S->sol[k][3] < cluster[3][Asscl[ind_gam[1]] - 1] &&
                                  S->sol[k][3] < cluster[3][Asscl[ind_gam[2]] - 1] &&
                                  S->sol[k][3] < cluster[3][Asscl[ind_gam[3]] - 1];
          }

          Bool_t cond_tot = cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && clusterEnergy;

          if (cond_tot && (cond_time_clus[0] || cond_time_clus[1]))
          {
            for (Int_t k = 0; k < 4; k++)
            {
              Param[k * 5] = cluster[0][Asscl[ind_gam[k]] - 1];
              Param[k * 5 + 1] = cluster[1][Asscl[ind_gam[k]] - 1];
              Param[k * 5 + 2] = cluster[2][Asscl[ind_gam[k]] - 1];
              Param[k * 5 + 3] = cluster[3][Asscl[ind_gam[k]] - 1];
              Param[k * 5 + 4] = cluster[4][Asscl[ind_gam[k]] - 1];

              Errors[k * 5] = clu_x_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]);
              Errors[k * 5 + 1] = clu_y_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]);
              Errors[k * 5 + 2] = clu_z_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]);
              // cm
              Errors[k * 5 + 3] = clu_time_error(Param[k * 5 + 4]); // ns
              Errors[k * 5 + 4] = clu_ene_error(Param[k * 5 + 4]);  // MeV

              Param[20 + k] = bhabha_mom[k];
              Errors[20 + k] = bhabha_mom_err[k];

              if (k < 3)
              {
                Param[24 + k] = bhabha_vtx[k];
                Errors[24 + k] = 0.;
              }
            }

            Float_t
                gamma_mom_tmp[4][4],
                fourKnetri_tmp[2][4],
                kaon_vel_tmp[2],
                y_axis[3],
                ip_tmp[2][3];

            for (Int_t k1 = jmin; k1 <= jmax; k1++)
            {
              kinematicFitObj->ParameterInitialization(Param, Errors);

              Tcorr = k1 * KLOE::T0;

              CHISQRTMP = kinematicFitObj->FitFunction(Tcorr);

              kinematicFitObj->GetResults(*X, *V, *X_init, *V_init);

              Reconstructor R;
              Solution S;

              for (Int_t k = 0; k < 4; k++)
              {
                R.SetClu(k, (*X)[k * 5],
                         (*X)[k * 5 + 1],
                         (*X)[k * 5 + 2],
                         (*X)[k * 5 + 3],
                         (*X)[k * 5 + 4]);

                R.SetClu(4, 0., 0., 0., 0., 0.);
                R.SetClu(5, 0., 0., 0., 0., 0.);
              }

              S = R.MySolve(selected);

              Float_t distance[4] = {0.};

              for (Int_t k = 0; k < 2; k++)
              {
                if (!S.error[k])
                {
                  for (Int_t l = 0; l < 4; l++)
                    neu_vtx[k][l] = S.sol[k][l];
                }
                else
                {
                  for (Int_t l = 0; l < 4; l++)
                    neu_vtx[k][l] = 999.;
                }

                for (Int_t l = 0; l < 4; l++)
                {
                  neutral_mom((*X)[l * 5], (*X)[l * 5 + 1], (*X)[l * 5 + 2], (*X)[l * 5 + 4], neu_vtx[k], gamma_mom_tmp[l]);

                  gamma_mom_tmp[l][4] = (*X)[l * 5];
                  gamma_mom_tmp[l][5] = (*X)[l * 5 + 1];
                  gamma_mom_tmp[l][6] = (*X)[l * 5 + 2];
                  gamma_mom_tmp[l][7] = (*X)[l * 5 + 3];
                }

                fourKnetri_tmp[k][0] = gamma_mom_tmp[0][0] + gamma_mom_tmp[1][0] + gamma_mom_tmp[2][0] + gamma_mom_tmp[3][0];
                fourKnetri_tmp[k][1] = gamma_mom_tmp[0][1] + gamma_mom_tmp[1][1] + gamma_mom_tmp[2][1] + gamma_mom_tmp[3][1];
                fourKnetri_tmp[k][2] = gamma_mom_tmp[0][2] + gamma_mom_tmp[1][2] + gamma_mom_tmp[2][2] + gamma_mom_tmp[3][2];
                fourKnetri_tmp[k][3] = gamma_mom_tmp[0][3] + gamma_mom_tmp[1][3] + gamma_mom_tmp[2][3] + gamma_mom_tmp[3][3];

                fourKnetri_tmp[k][4] = sqrt(pow(fourKnetri_tmp[k][0], 2) + pow(fourKnetri_tmp[k][1], 2) + pow(fourKnetri_tmp[k][2], 2));
                fourKnetri_tmp[k][5] = sqrt(pow(fourKnetri_tmp[k][3], 2) - pow(fourKnetri_tmp[k][4], 2));

                kaon_vel_tmp[k] = PhysicsConstants::cVel * fourKnetri_tmp[k][4] / fourKnetri_tmp[k][3];

                y_axis[0] = 0.;
                y_axis[1] = (*X)[21];
                y_axis[2] = 0.;

                eventAnalysis->IPBoostCorr(bhabha_vtx.data(), y_axis, neu_vtx[k], fourKnetri_tmp[k], ip_tmp[k]);

                ip_tmp[k][0] = bhabha_vtx[0];
                ip_tmp[k][1] = bhabha_vtx[1];
                if (abs(ip_tmp[k][2] - bhabha_vtx[2]) > 2)
                  ip_tmp[k][2] = bhabha_vtx[2];

                dist_tmp[k] = sqrt(pow(neu_vtx[k][0] - ip_tmp[k][0], 2) +
                                   pow(neu_vtx[k][1] - ip_tmp[k][1], 2) +
                                   pow(neu_vtx[k][2] - ip_tmp[k][2], 2));

                value[k] = sqrt(pow(neu_vtx[k][3] - (dist_tmp[k] / kaon_vel_tmp[k]), 2) + pow(fourKnetri_tmp[k][5] - PhysicsConstants::mK0, 2));

                if (TMath::IsNaN(value[k]))
                  value[k] = 999999.;
              }

              cond_time_clus[0] = S.sol[0][3] < (*X)(3) &&
                                  S.sol[0][3] < (*X)(8) &&
                                  S.sol[0][3] < (*X)(13) &&
                                  S.sol[0][3] < (*X)(18);

              cond_time_clus[1] = S.sol[1][3] < (*X)(3) &&
                                  S.sol[1][3] < (*X)(8) &&
                                  S.sol[1][3] < (*X)(13) &&
                                  S.sol[1][3] < (*X)(18);

              if (abs(CHISQRTMP) < abs(CHISQRMIN))
              {
                if (cond_time_clus[0] && value[0] < value[1])
                {
                  isConverged = 1;
                  FUNVALMIN = FUNVALTMP;
                  CHISQRMIN = CHISQRTMP;

                  Chi2TriKinFit = CHISQRMIN;

                  kinematicFitObj->GetResults((*X_min), (*V_min), (*X_init_min), (*V_init));

                  g4takentri_kinfit[0] = ind_gam[0];
                  g4takentri_kinfit[1] = ind_gam[1];
                  g4takentri_kinfit[2] = ind_gam[2];
                  g4takentri_kinfit[3] = ind_gam[3];

                  neu_vtx_min[0] = neu_vtx[0][0];
                  neu_vtx_min[1] = neu_vtx[0][1];
                  neu_vtx_min[2] = neu_vtx[0][2];
                  neu_vtx_min[3] = neu_vtx[0][3];

                  for (Int_t l = 0; l < 4; l++)
                  {
                    distance[l] = sqrt(pow((*X)[l * 5] - neu_vtx[0][0], 2) +
                                       pow((*X)[l * 5 + 1] - neu_vtx[0][1], 2) +
                                       pow((*X)[l * 5 + 2] - neu_vtx[0][2], 2));

                    gamma_mom_final[l][0] = (*X)[l * 5 + 4] * (((*X)[l * 5] - neu_vtx[0][0]) / distance[l]);
                    gamma_mom_final[l][1] = (*X)[l * 5 + 4] * (((*X)[l * 5 + 1] - neu_vtx[0][1]) / distance[l]);
                    gamma_mom_final[l][2] = (*X)[l * 5 + 4] * (((*X)[l * 5 + 2] - neu_vtx[0][2]) / distance[l]);
                    gamma_mom_final[l][3] = (*X)[l * 5 + 4];
                    gamma_mom_final[l][4] = (*X)[l * 5];
                    gamma_mom_final[l][5] = (*X)[l * 5 + 1];
                    gamma_mom_final[l][6] = (*X)[l * 5 + 2];
                    gamma_mom_final[l][7] = (*X)[l * 5 + 3];
                  }

                  fourKnetri_kinfit[0] = gamma_mom_final[0][0] + gamma_mom_final[1][0] + gamma_mom_final[2][0] + gamma_mom_final[3][0];
                  fourKnetri_kinfit[1] = gamma_mom_final[0][1] + gamma_mom_final[1][1] + gamma_mom_final[2][1] + gamma_mom_final[3][1];
                  fourKnetri_kinfit[2] = gamma_mom_final[0][2] + gamma_mom_final[1][2] + gamma_mom_final[2][2] + gamma_mom_final[3][2];
                  fourKnetri_kinfit[3] = gamma_mom_final[0][3] + gamma_mom_final[1][3] + gamma_mom_final[2][3] + gamma_mom_final[3][3];
                  fourKnetri_kinfit[4] = sqrt(pow(fourKnetri_kinfit[0], 2) + pow(fourKnetri_kinfit[1], 2) + pow(fourKnetri_kinfit[2], 2));
                  fourKnetri_kinfit[5] = sqrt(pow(fourKnetri_kinfit[3], 2) - pow(fourKnetri_kinfit[4], 2));
                  fourKnetri_kinfit[6] = neu_vtx_min[0];
                  fourKnetri_kinfit[7] = neu_vtx_min[1];
                  fourKnetri_kinfit[8] = neu_vtx_min[2];
                  fourKnetri_kinfit[9] = neu_vtx_min[3];

                  iptri_kinfit[0] = ip_tmp[0][0];
                  iptri_kinfit[1] = ip_tmp[0][1];
                  iptri_kinfit[2] = ip_tmp[0][2];

                  bunchnum = k1;
                }
                else if (cond_time_clus[1] && value[1] < value[0])
                {
                  isConverged = 1;
                  FUNVALMIN = FUNVALTMP;
                  CHISQRMIN = CHISQRTMP;

                  Chi2TriKinFit = CHISQRMIN;

                  kinematicFitObj->GetResults((*X_min), (*V_min), (*X_init_min), (*V_init));

                  g4takentri_kinfit[0] = ind_gam[0];
                  g4takentri_kinfit[1] = ind_gam[1];
                  g4takentri_kinfit[2] = ind_gam[2];
                  g4takentri_kinfit[3] = ind_gam[3];

                  neu_vtx_min[0] = neu_vtx[1][0];
                  neu_vtx_min[1] = neu_vtx[1][1];
                  neu_vtx_min[2] = neu_vtx[1][2];
                  neu_vtx_min[3] = neu_vtx[1][3];

                  for (Int_t l = 0; l < 4; l++)
                  {
                    distance[l] = sqrt(pow((*X)[l * 5] - neu_vtx[1][0], 2) +
                                       pow((*X)[l * 5 + 1] - neu_vtx[1][1], 2) +
                                       pow((*X)[l * 5 + 2] - neu_vtx[1][2], 2));

                    gamma_mom_final[l][0] = (*X)[l * 5 + 4] * (((*X)[l * 5] - neu_vtx[1][0]) / distance[l]);
                    gamma_mom_final[l][1] = (*X)[l * 5 + 4] * (((*X)[l * 5 + 1] - neu_vtx[1][1]) / distance[l]);
                    gamma_mom_final[l][2] = (*X)[l * 5 + 4] * (((*X)[l * 5 + 2] - neu_vtx[1][2]) / distance[l]);
                    gamma_mom_final[l][3] = (*X)[l * 5 + 4];
                    gamma_mom_final[l][4] = (*X)[l * 5];
                    gamma_mom_final[l][5] = (*X)[l * 5 + 1];
                    gamma_mom_final[l][6] = (*X)[l * 5 + 2];
                    gamma_mom_final[l][7] = (*X)[l * 5 + 3];
                  }

                  fourKnetri_kinfit[0] = gamma_mom_final[0][0] + gamma_mom_final[1][0] + gamma_mom_final[2][0] + gamma_mom_final[3][0];
                  fourKnetri_kinfit[1] = gamma_mom_final[0][1] + gamma_mom_final[1][1] + gamma_mom_final[2][1] + gamma_mom_final[3][1];
                  fourKnetri_kinfit[2] = gamma_mom_final[0][2] + gamma_mom_final[1][2] + gamma_mom_final[2][2] + gamma_mom_final[3][2];
                  fourKnetri_kinfit[3] = gamma_mom_final[0][3] + gamma_mom_final[1][3] + gamma_mom_final[2][3] + gamma_mom_final[3][3];
                  fourKnetri_kinfit[4] = sqrt(pow(fourKnetri_kinfit[0], 2) + pow(fourKnetri_kinfit[1], 2) + pow(fourKnetri_kinfit[2], 2));
                  fourKnetri_kinfit[5] = sqrt(pow(fourKnetri_kinfit[3], 2) - pow(fourKnetri_kinfit[4], 2));
                  fourKnetri_kinfit[6] = neu_vtx_min[0];
                  fourKnetri_kinfit[7] = neu_vtx_min[1];
                  fourKnetri_kinfit[8] = neu_vtx_min[2];
                  fourKnetri_kinfit[9] = neu_vtx_min[3];

                  iptri_kinfit[0] = ip_tmp[1][0];
                  iptri_kinfit[1] = ip_tmp[1][1];
                  iptri_kinfit[2] = ip_tmp[1][2];

                  bunchnum = k1;
                }
                else if (isConverged == 0)
                {
                  neu_vtx_min[0] = 0.;
                  neu_vtx_min[1] = 999;
                  neu_vtx_min[2] = 999;
                  neu_vtx_min[3] = 999;

                  for (Int_t l = 0; l < 4; l++)
                  {
                    gamma_mom_final[l][0] = 999.;
                    gamma_mom_final[l][1] = 999.;
                    gamma_mom_final[l][2] = 999.;
                    gamma_mom_final[l][3] = 999.;
                    gamma_mom_final[l][4] = 999.;
                    gamma_mom_final[l][5] = 999.;
                    gamma_mom_final[l][6] = 999.;
                    gamma_mom_final[l][7] = 999.;
                  }

                  fourKnetri_kinfit[0] = 999.;
                  fourKnetri_kinfit[1] = 999.;
                  fourKnetri_kinfit[2] = 999.;
                  fourKnetri_kinfit[3] = 999.;
                  fourKnetri_kinfit[4] = 999.;
                  fourKnetri_kinfit[5] = 999.;
                  fourKnetri_kinfit[6] = 999.;
                  fourKnetri_kinfit[7] = 999.;
                  fourKnetri_kinfit[8] = 999.;
                  fourKnetri_kinfit[9] = 999.;

                  iptri_kinfit[0] = 999.;
                  iptri_kinfit[1] = 999.;
                  iptri_kinfit[2] = 999.;

                  bunchnum = 999;
                }
              }
            }
          }

          delete S;
          delete R;
        }

  delete V;
  delete D;
  delete D_T;
  delete V_final;
  delete V_aux;
  delete V_min;
  delete Aux;
  delete V_invert;
  delete V_init;
  delete X;
  delete C;
  delete X_final;
  delete L;
  delete CORR;
  delete X_init;
  delete X_min;
  delete C_min;
  delete L_min;
  delete C_aux;
  delete L_aux;
  delete X_init_min;
  delete X_init_aux;

  if (isConverged == 0)
    return ErrorHandling::ErrorCodes::TRILATERATION_KIN_FIT;

  return ErrorHandling::ErrorCodes::NO_ERROR;
}