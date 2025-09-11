#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>

#include "../inc/signalKinFit.h"

namespace KLOE
{
	SignalKinFit::SignalKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, Int_t jmin, Int_t jmax, ErrorHandling::ErrorLogs &logger) : KinFitter("Trilateration", N_free, N_const, M, 0, loopcount, chiSqrStep, logger)
	{
		_V.ResizeTo(N_free + N_const, N_free + N_const);
		_D.ResizeTo(M, N_free + N_const);
		_D_T.ResizeTo(N_free + N_const, M);
		_V_final.ResizeTo(N_free + N_const, N_free + N_const);
		_V_aux.ResizeTo(N_free + N_const, N_free + N_const);
		_V_min.ResizeTo(N_free + N_const, N_free + N_const);
		_Aux.ResizeTo(M, M);
		_V_invert.ResizeTo(N_free, N_free);
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

		_selected.resize(4);
		for (Int_t i = 1; i <= 4; i++)
			_selected[i - 1] = i;
		_ind_gam = std::make_unique<Int_t[]>(4);

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

		gErrorIgnoreLevel = 6001;
	}

	SignalKinFit::~SignalKinFit()
	{
	}

	ErrorHandling::ErrorCodes SignalKinFit::Reconstruct()
	{
		volatile Float_t
			CHISQRTMP = 999.,
			FUNVALTMP = 999999.,
			Tcorr;

		_CHISQRMIN = 999999.;
		_isConverged = 0;

		for (Int_t j1 = 0; j1 < _NeuClusters.size() - 3; j1++)
			for (Int_t j2 = j1 + 1; j2 < _NeuClusters.size() - 2; j2++)
				for (Int_t j3 = j2 + 1; j3 < _NeuClusters.size() - 1; j3++)
					for (Int_t j4 = j3 + 1; j4 < _NeuClusters.size(); j4++)
					{
						_ind_gam[0] = j1;
						_ind_gam[1] = j2;
						_ind_gam[2] = j3;
						_ind_gam[3] = j4;

						volatile Bool_t
							clusterEnergy = false,
							cond_clus[4] = {false, false, false, false},
							cond_time_clus[2] = {false, false};

						for (Int_t k = 0; k < 4; k++)
						{
							_R.SetClu(k, _cluster[0][_NeuClusters[_ind_gam[k]] - 1],
									  _cluster[1][_NeuClusters[_ind_gam[k]] - 1],
									  _cluster[2][_NeuClusters[_ind_gam[k]] - 1],
									  _cluster[3][_NeuClusters[_ind_gam[k]] - 1],
									  _cluster[4][_NeuClusters[_ind_gam[k]] - 1]);

							_R.SetClu(4, 0., 0., 0., 0., 0.);
							_R.SetClu(5, 0., 0., 0., 0., 0.);
						}

						_S = _R.MySolve(_selected.data());

						clusterEnergy = _cluster[4][_NeuClusters[_ind_gam[0]] - 1] > MIN_CLU_ENE &&
										_cluster[4][_NeuClusters[_ind_gam[1]] - 1] > MIN_CLU_ENE &&
										_cluster[4][_NeuClusters[_ind_gam[2]] - 1] > MIN_CLU_ENE &&
										_cluster[4][_NeuClusters[_ind_gam[3]] - 1] > MIN_CLU_ENE;

						for (Int_t k = 0; k < 4; k++)
						{
							cond_clus[k] =
								_cluster[3][_NeuClusters[_ind_gam[k]] - 1] > 0 &&
								_cluster[0][_NeuClusters[_ind_gam[k]] - 1] != 0 &&
								_cluster[1][_NeuClusters[_ind_gam[k]] - 1] != 0 &&
								_cluster[2][_NeuClusters[_ind_gam[k]] - 1] != 0;

							if (k < 2)
								cond_time_clus[k] = _S.sol[k][3] < _cluster[3][_NeuClusters[_ind_gam[0]] - 1] &&
													_S.sol[k][3] < _cluster[3][_NeuClusters[_ind_gam[1]] - 1] &&
													_S.sol[k][3] < _cluster[3][_NeuClusters[_ind_gam[2]] - 1] &&
													_S.sol[k][3] < _cluster[3][_NeuClusters[_ind_gam[3]] - 1];
						}

						Bool_t cond_tot = cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && clusterEnergy;

						if (cond_tot && (cond_time_clus[0] || cond_time_clus[1]))
						{
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
								_Errors[k * 5 + 4] = clu_ene_error(_Param[k * 5 + 4]);	// MeV

								_Param[20 + k] = _bhabha_mom[k];
								_Errors[20 + k] = _bhabha_mom_err[k];

								if (k < 3)
								{
									_Param[24 + k] = _bhabha_vtx[k];
									_Errors[24 + k] = 0.;
								}
							}

							Float_t
								gamma_mom_tmp[4][4],
								fourKnetri_tmp[2][4],
								kaon_vel_tmp[2],
								y_axis[3],
								ip_tmp[2][3];

							KinFitter::ParameterInitialization(_Param.data(), _Errors.data());

							CHISQRTMP = KinFitter::FitFunction(Tcorr);

							KinFitter::GetResults(_X, _V, _X_init, _V_init, _C, _L);

							for (Int_t k = 0; k < 4; k++)
							{
								_R.SetClu(k, _X[k * 5],
										  _X[k * 5 + 1],
										  _X[k * 5 + 2],
										  _X[k * 5 + 3],
										  _X[k * 5 + 4]);

								_R.SetClu(4, 0., 0., 0., 0., 0.);
								_R.SetClu(5, 0., 0., 0., 0., 0.);
							}

							_S = _R.MySolve(_selected.data());

							Float_t distance[4] = {0.};

							Float_t neu_vtx[2][4];
							Float_t value[2] = {999999., 999999.};
							Float_t dist_tmp[2] = {0.};

							for (Int_t k = 0; k < 2; k++)
							{
								if (!_S.error[k])
								{
									for (Int_t l = 0; l < 4; l++)
										neu_vtx[k][l] = _S.sol[k][l];
								}
								else
								{
									for (Int_t l = 0; l < 4; l++)
										neu_vtx[k][l] = 999.;
								}

								for (Int_t l = 0; l < 4; l++)
								{
									neutral_mom(_X[l * 5], _X[l * 5 + 1], _X[l * 5 + 2], _X[l * 5 + 4], neu_vtx[k], gamma_mom_tmp[l]);

									gamma_mom_tmp[l][4] = _X[l * 5];
									gamma_mom_tmp[l][5] = _X[l * 5 + 1];
									gamma_mom_tmp[l][6] = _X[l * 5 + 2];
									gamma_mom_tmp[l][7] = _X[l * 5 + 3];
								}

								fourKnetri_tmp[k][0] = gamma_mom_tmp[0][0] + gamma_mom_tmp[1][0] + gamma_mom_tmp[2][0] + gamma_mom_tmp[3][0];
								fourKnetri_tmp[k][1] = gamma_mom_tmp[0][1] + gamma_mom_tmp[1][1] + gamma_mom_tmp[2][1] + gamma_mom_tmp[3][1];
								fourKnetri_tmp[k][2] = gamma_mom_tmp[0][2] + gamma_mom_tmp[1][2] + gamma_mom_tmp[2][2] + gamma_mom_tmp[3][2];
								fourKnetri_tmp[k][3] = gamma_mom_tmp[0][3] + gamma_mom_tmp[1][3] + gamma_mom_tmp[2][3] + gamma_mom_tmp[3][3];

								fourKnetri_tmp[k][4] = sqrt(pow(fourKnetri_tmp[k][0], 2) + pow(fourKnetri_tmp[k][1], 2) + pow(fourKnetri_tmp[k][2], 2));
								fourKnetri_tmp[k][5] = sqrt(pow(fourKnetri_tmp[k][3], 2) - pow(fourKnetri_tmp[k][4], 2));

								kaon_vel_tmp[k] = cVel * fourKnetri_tmp[k][4] / fourKnetri_tmp[k][3];

								y_axis[0] = 0.;
								y_axis[1] = _X[21];
								y_axis[2] = 0.;

								ChargedVtxRec::IPBoostCorr(_bhabha_vtx.data(), y_axis, neu_vtx[k], fourKnetri_tmp[k], ip_tmp[k]);

								ip_tmp[k][0] = _bhabha_vtx[0];
								ip_tmp[k][1] = _bhabha_vtx[1];
								if (abs(ip_tmp[k][2] - _bhabha_vtx[2]) > 2)
									ip_tmp[k][2] = _bhabha_vtx[2];

								dist_tmp[k] = sqrt(pow(neu_vtx[k][0] - ip_tmp[k][0], 2) +
												   pow(neu_vtx[k][1] - ip_tmp[k][1], 2) +
												   pow(neu_vtx[k][2] - ip_tmp[k][2], 2));

								value[k] = sqrt(pow(neu_vtx[k][3] - (dist_tmp[k] / kaon_vel_tmp[k]), 2) + pow(fourKnetri_tmp[k][5] - mK0, 2));

								if (TMath::IsNaN(value[k]))
									value[k] = 999999.;
							}

							cond_time_clus[0] = _S.sol[0][3] < _X(3) &&
												_S.sol[0][3] < _X(8) &&
												_S.sol[0][3] < _X(13) &&
												_S.sol[0][3] < _X(18);

							cond_time_clus[1] = _S.sol[1][3] < _X(3) &&
												_S.sol[1][3] < _X(8) &&
												_S.sol[1][3] < _X(13) &&
												_S.sol[1][3] < _X(18);

							if (abs(CHISQRTMP) < abs(_CHISQRMIN))
							{
								if (cond_time_clus[0] && value[0] < value[1])
								{
									_isConverged = 1;
									_FUNVALMIN = FUNVALTMP;
									_CHISQRMIN = CHISQRTMP;

									_Chi2TriKinFit = _CHISQRMIN;
									KinFitter::GetResults(_X_min, _V_min, _X_init_min, _V_init, _C_min, _L_min);

									_g4takentri_kinfit[0] = _ind_gam[0];
									_g4takentri_kinfit[1] = _ind_gam[1];
									_g4takentri_kinfit[2] = _ind_gam[2];
									_g4takentri_kinfit[3] = _ind_gam[3];

									_neu_vtx_min[0] = neu_vtx[0][0];
									_neu_vtx_min[1] = neu_vtx[0][1];
									_neu_vtx_min[2] = neu_vtx[0][2];
									_neu_vtx_min[3] = neu_vtx[0][3];

									for (Int_t l = 0; l < 4; l++)
									{
										distance[l] = sqrt(pow(_X[l * 5] - neu_vtx[0][0], 2) +
														   pow(_X[l * 5 + 1] - neu_vtx[0][1], 2) +
														   pow(_X[l * 5 + 2] - neu_vtx[0][2], 2));

										_gamma_mom_final[l][0] = _X[l * 5 + 4] * ((_X[l * 5] - neu_vtx[0][0]) / distance[l]);
										_gamma_mom_final[l][1] = _X[l * 5 + 4] * ((_X[l * 5 + 1] - neu_vtx[0][1]) / distance[l]);
										_gamma_mom_final[l][2] = _X[l * 5 + 4] * ((_X[l * 5 + 2] - neu_vtx[0][2]) / distance[l]);
										_gamma_mom_final[l][3] = _X[l * 5 + 4];
										_gamma_mom_final[l][4] = _X[l * 5];
										_gamma_mom_final[l][5] = _X[l * 5 + 1];
										_gamma_mom_final[l][6] = _X[l * 5 + 2];
										_gamma_mom_final[l][7] = _X[l * 5 + 3];
									}

									_fourKnetri_kinfit[0] = _gamma_mom_final[0][0] + _gamma_mom_final[1][0] + _gamma_mom_final[2][0] + _gamma_mom_final[3][0];
									_fourKnetri_kinfit[1] = _gamma_mom_final[0][1] + _gamma_mom_final[1][1] + _gamma_mom_final[2][1] + _gamma_mom_final[3][1];
									_fourKnetri_kinfit[2] = _gamma_mom_final[0][2] + _gamma_mom_final[1][2] + _gamma_mom_final[2][2] + _gamma_mom_final[3][2];
									_fourKnetri_kinfit[3] = _gamma_mom_final[0][3] + _gamma_mom_final[1][3] + _gamma_mom_final[2][3] + _gamma_mom_final[3][3];
									_fourKnetri_kinfit[4] = sqrt(pow(_fourKnetri_kinfit[0], 2) + pow(_fourKnetri_kinfit[1], 2) + pow(_fourKnetri_kinfit[2], 2));
									_fourKnetri_kinfit[5] = sqrt(pow(_fourKnetri_kinfit[3], 2) - pow(_fourKnetri_kinfit[4], 2));
									_fourKnetri_kinfit[6] = _neu_vtx_min[0];
									_fourKnetri_kinfit[7] = _neu_vtx_min[1];
									_fourKnetri_kinfit[8] = _neu_vtx_min[2];
									_fourKnetri_kinfit[9] = _neu_vtx_min[3];

									_iptri_kinfit[0] = ip_tmp[0][0];
									_iptri_kinfit[1] = ip_tmp[0][1];
									_iptri_kinfit[2] = ip_tmp[0][2];
								}
								else if (cond_time_clus[1] && value[1] < value[0])
								{
									_isConverged = 1;
									_FUNVALMIN = FUNVALTMP;
									_CHISQRMIN = CHISQRTMP;

									_Chi2TriKinFit = _CHISQRMIN;

									KinFitter::GetResults(_X_min, _V_min, _X_init_min, _V_init, _C_min, _L_min);

									_g4takentri_kinfit[0] = _ind_gam[0];
									_g4takentri_kinfit[1] = _ind_gam[1];
									_g4takentri_kinfit[2] = _ind_gam[2];
									_g4takentri_kinfit[3] = _ind_gam[3];

									_neu_vtx_min[0] = neu_vtx[1][0];
									_neu_vtx_min[1] = neu_vtx[1][1];
									_neu_vtx_min[2] = neu_vtx[1][2];
									_neu_vtx_min[3] = neu_vtx[1][3];

									for (Int_t l = 0; l < 4; l++)
									{
										distance[l] = sqrt(pow(_X[l * 5] - neu_vtx[1][0], 2) +
														   pow(_X[l * 5 + 1] - neu_vtx[1][1], 2) +
														   pow(_X[l * 5 + 2] - neu_vtx[1][2], 2));

										_gamma_mom_final[l][0] = _X[l * 5 + 4] * ((_X[l * 5] - neu_vtx[1][0]) / distance[l]);
										_gamma_mom_final[l][1] = _X[l * 5 + 4] * ((_X[l * 5 + 1] - neu_vtx[1][1]) / distance[l]);
										_gamma_mom_final[l][2] = _X[l * 5 + 4] * ((_X[l * 5 + 2] - neu_vtx[1][2]) / distance[l]);
										_gamma_mom_final[l][3] = _X[l * 5 + 4];
										_gamma_mom_final[l][4] = _X[l * 5];
										_gamma_mom_final[l][5] = _X[l * 5 + 1];
										_gamma_mom_final[l][6] = _X[l * 5 + 2];
										_gamma_mom_final[l][7] = _X[l * 5 + 3];
									}

									_fourKnetri_kinfit[0] = _gamma_mom_final[0][0] + _gamma_mom_final[1][0] + _gamma_mom_final[2][0] + _gamma_mom_final[3][0];
									_fourKnetri_kinfit[1] = _gamma_mom_final[0][1] + _gamma_mom_final[1][1] + _gamma_mom_final[2][1] + _gamma_mom_final[3][1];
									_fourKnetri_kinfit[2] = _gamma_mom_final[0][2] + _gamma_mom_final[1][2] + _gamma_mom_final[2][2] + _gamma_mom_final[3][2];
									_fourKnetri_kinfit[3] = _gamma_mom_final[0][3] + _gamma_mom_final[1][3] + _gamma_mom_final[2][3] + _gamma_mom_final[3][3];
									_fourKnetri_kinfit[4] = sqrt(pow(_fourKnetri_kinfit[0], 2) + pow(_fourKnetri_kinfit[1], 2) + pow(_fourKnetri_kinfit[2], 2));
									_fourKnetri_kinfit[5] = sqrt(pow(_fourKnetri_kinfit[3], 2) - pow(_fourKnetri_kinfit[4], 2));
									_fourKnetri_kinfit[6] = _neu_vtx_min[0];
									_fourKnetri_kinfit[7] = _neu_vtx_min[1];
									_fourKnetri_kinfit[8] = _neu_vtx_min[2];
									_fourKnetri_kinfit[9] = _neu_vtx_min[3];

									_iptri_kinfit[0] = ip_tmp[1][0];
									_iptri_kinfit[1] = ip_tmp[1][1];
									_iptri_kinfit[2] = ip_tmp[1][2];
								}
								else if (_isConverged == 0)
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
								}
							}
						}
					}

		if (_isConverged)
			return ErrorHandling::ErrorCodes::NO_ERROR;
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

			_Chi2TriKinFit = 999.;

			return ErrorHandling::ErrorCodes::TRILATERATION_KIN_FIT;
		}
	}
}