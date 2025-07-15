#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

#include "../../../Include/Codes/reconstructor.h"
#include "../../../Include/const.h"
#include "../../../Include/Codes/uncertainties.h"
#include "../../../Include/Codes/charged_mom.h"
#include "../../../Include/Codes/neutral_mom.h"
#include "../../../Include/Codes/lorentz_transf.h"
#include "../../../Include/Codes/plane_intersection.h"
#include "../../../Include/Codes/closest_approach.h"
#include "../../../Include/Codes/kloe_class.h"
#include "../../../Include/Codes/chi2_dist.h"
#include "../inc/trilateration.hpp"

Float_t T0;

const Float_t Trf = 2.715; // ns - time of a bunch

void tri_neurec_kinfit(int first_file, int last_file) //	arguments are: 1. Number of points
																											//								 2. First analyzed file
																											//                 3. Last analyzed file
{
	Reconstructor R;
	Solution S;

	const Int_t N = 27, M = 9;
	Int_t selected[4] = {1, 2, 3, 4};

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	TString name = "neuvtx_tri_kin_fit_" + std::to_string(first_file) + "_" + std::to_string(last_file) + ".root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_tri_kin_fit", "Neu vtx rec with trilateration kin fit");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_mom_err[4], bhabha_vtx[3], bhabha_vtx_err[3];

	chain->SetBranchAddress("Bpx", &bhabha_mom[0]);
	chain->SetBranchAddress("Bpy", &bhabha_mom[1]);
	chain->SetBranchAddress("Bpz", &bhabha_mom[2]);
	chain->SetBranchAddress("Broots", &bhabha_mom[3]);

	chain->SetBranchAddress("Bwidpx", &bhabha_mom_err[0]);
	chain->SetBranchAddress("Bwidpy", &bhabha_mom_err[1]);
	chain->SetBranchAddress("Bwidpz", &bhabha_mom_err[2]);
	chain->SetBranchAddress("Brootserr", &bhabha_mom_err[3]);

	chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
	chain->SetBranchAddress("By", &bhabha_vtx[1]);
	chain->SetBranchAddress("Bz", &bhabha_vtx[2]);
	chain->SetBranchAddress("Bsx", &bhabha_vtx_err[0]);
	chain->SetBranchAddress("Bsy", &bhabha_vtx_err[1]);
	chain->SetBranchAddress("Bsz", &bhabha_vtx_err[2]);

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth;
	Float_t cluster[5][500], Kchboost[9], Knerec[9], Knemc[9], ipmc[3], ip[3], Dtmc;

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster[0]);
	chain->SetBranchAddress("Ycl", cluster[1]);
	chain->SetBranchAddress("Zcl", cluster[2]);
	chain->SetBranchAddress("TclOld", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);

	chain->SetBranchAddress("Knemc", Knemc);
	chain->SetBranchAddress("Knerec", Knerec);

	chain->SetBranchAddress("Kchboost", Kchboost);
	chain->SetBranchAddress("Dtmc", &Dtmc);

	chain->SetBranchAddress("ipmc", ipmc);
	chain->SetBranchAddress("ip", ip);

	//! Parameters for fitting

	Double_t P[N], DP[N];

	Double_t CHISQR;
	Float_t CHISQRMIN;

	//!

	Int_t nentries = (Int_t)chain->GetEntries();

	Bool_t clusterEnergy, solError, cond_clus[4];
	Int_t ind_gam[4], sort_index[2][3], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged;
	Float_t neu_vtx[2][4], inv_mPi0[2][3][2], gamma_mom[2][4][4], mass_pair[2][3];

	Float_t P1[3 * N + M + 3], min_value, min_value_def, neu_vtx_min[2][4], clusters_min[4][5], neu_vtx_min_final[4],
			mom_kaon[2][4], ene_kaon[2], v_kaon[2], length_kaon[2], diff_kaon[2], length, length_mc, length_rec;

	Float_t length_ch, time_ch, velocity_ch, gamma_mom_final[4][8], fourKnetri_kinfit[10], constraints[M], iptri_kinfit[3], y_axis[3] = {0., 1., 0.}, ip_tri[2][3], bhabha_mom_fit[4];
	Int_t g4takentri_kinfit[4];

	TBranch *b_gamma1tri = tree->Branch("fourgamma1tri_kinfit", gamma_mom_final[0], "fourgamma1tri_kinfit[8]/F");
	TBranch *b_gamma2tri = tree->Branch("fourgamma2tri_kinfit", gamma_mom_final[1], "fourgamma2tri_kinfit[8]/F");
	TBranch *b_gamma3tri = tree->Branch("fourgamma3tri_kinfit", gamma_mom_final[2], "fourgamma3tri_kinfit[8]/F");
	TBranch *b_gamma4tri = tree->Branch("fourgamma4tri_kinfit", gamma_mom_final[3], "fourgamma4tri_kinfit[8]/F");
	TBranch *b_iptri = tree->Branch("iptri_kinfit", iptri_kinfit, "iptri_kinfit[3]/F");
	TBranch *b_Knetri = tree->Branch("fourKnetri_kinfit", fourKnetri_kinfit, "fourKnetri_kinfit[10]/F");
	TBranch *b_done = tree->Branch("done4_kinfit", &found_best, "done4_kinfit/I");
	TBranch *b_fourg4taken = tree->Branch("g4takentri_kinfit", g4takentri_kinfit, "g4takentri_kinfit[4]/I");

	TBranch *b_chisqr = tree->Branch("chi2min", &CHISQRMIN, "chi2min/F");

	TH1 *chi2_hist = new TH1F("chi2", "", 100, 0, 100);
	TH2 *chi2_corr = new TH2F("chi2_corr", "", N, 1, N, N, 1, N);

	for (Int_t i = 0; i < nentries/10.; i++)
	{
		chain->GetEntry(i);

		min_value_def = 99999.;
		CHISQRMIN = 99999.;
		found_best = 0;

		if (nclu >= 4 && (mctruth == 1 || mctruth == 2))
		{

			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							CHISQR = 0.;

							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;

							clusterEnergy = (cluster[4][ind_gam[0]] > MIN_CLU_ENE && cluster[4][ind_gam[1]] > MIN_CLU_ENE && cluster[4][ind_gam[2]] > MIN_CLU_ENE && cluster[4][ind_gam[3]] > MIN_CLU_ENE);

							cond_clus[0] = cluster[3][ind_gam[0]] > 0 && cluster[0][ind_gam[0]] != 0 && cluster[1][ind_gam[0]] != 0 && cluster[2][ind_gam[0]] != 0;
							cond_clus[1] = cluster[3][ind_gam[1]] > 0 && cluster[0][ind_gam[1]] != 0 && cluster[1][ind_gam[1]] != 0 && cluster[2][ind_gam[1]] != 0;
							cond_clus[2] = cluster[3][ind_gam[2]] > 0 && cluster[0][ind_gam[2]] != 0 && cluster[1][ind_gam[2]] != 0 && cluster[2][ind_gam[2]] != 0;
							cond_clus[3] = cluster[3][ind_gam[3]] > 0 && cluster[0][ind_gam[3]] != 0 && cluster[1][ind_gam[3]] != 0 && cluster[2][ind_gam[3]] != 0;

							if (cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && clusterEnergy)
							{
								// Setting clusters for a solution
								for (Int_t l = 0; l < 4; l++)
									R.SetClu(l, cluster[0][ind_gam[l]],
													 cluster[1][ind_gam[l]],
													 cluster[2][ind_gam[l]],
													 cluster[3][ind_gam[l]] - T0,
													 cluster[4][ind_gam[l]]);

								R.SetClu(4, 0., 0., 0., 0., 0.);
								R.SetClu(5, 0., 0., 0., 0., 0.);

								S = R.MySolve(selected);
								////////////////////////////////////

								solError = S.error[0] == false || S.error[1] == false;

								if (solError)
								{
									// !ITERATING OVER SOLUTIONS OF THE RECONSTRUCTION

									ROOT::Math::Minimizer *minimum;

									minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

									// set tolerance , etc...
									minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
									minimum->SetTolerance(0.01);
									minimum->SetPrintLevel(0);

									// create function wrapper for minimizer
									// a IMultiGenFunction type
									ROOT::Math::Functor f(&trilateration_chi_square, 3 * N + M);
									minimum->SetFunction(f);
									////////////////////////////////////////////////////////////

									for (Int_t l = 0; l < 3 * N + M; l++)
									{
										minimum->SetVariable(l, std::to_string(l), 0., 0.001);
										if (l >= N)
											minimum->FixVariable(l);
									}

									// Building the parameters for a fit
									for (Int_t k = 0; k < 4; k++)
									{
										P[k * 5] = cluster[0][ind_gam[k]];		 // X_k
										P[k * 5 + 1] = cluster[1][ind_gam[k]]; // Y_k
										P[k * 5 + 2] = cluster[2][ind_gam[k]]; // Z_k
										P[k * 5 + 3] = cluster[3][ind_gam[k]]; // T_k
										P[k * 5 + 4] = cluster[4][ind_gam[k]]; // E_k

										DP[k * 5] = 1.2;																						  // cm

										if(sqrt(pow(P[k * 5],2) + pow(P[k * 5 + 1],2)) > 200)
										{
											DP[k * 5 + 1] = 1.2;																				// cm
											DP[k * 5 + 2] = 1.2 / sqrt(cluster[4][ind_gam[k]] / 1000.); // cm
										}
										else
										{
											DP[k * 5 + 1] = 1.2 / sqrt(cluster[4][ind_gam[k]] / 1000.);	// cm
											DP[k * 5 + 2] = 1.2; // cm
										}

										DP[k * 5 + 3] = clu_time_error(cluster[4][ind_gam[k]]);
										DP[k * 5 + 4] = clu_ene_error(cluster[4][ind_gam[k]]);
									}

									P[20] = bhabha_mom[0]; // Px_phi
									P[21] = bhabha_mom[1]; // Py_phi
									P[22] = bhabha_mom[2]; // Pz_phi
									P[23] = bhabha_mom[3]; // Pz_phi

									DP[20] = bhabha_mom_err[0];
									DP[21] = bhabha_mom_err[1];
									DP[22] = bhabha_mom_err[2];
									DP[23] = bhabha_mom_err[3];

									P[24] = bhabha_vtx[0]; // x_phi
									P[25] = bhabha_vtx[1]; // y_phi
									P[26] = bhabha_vtx[2]; // z_phi

									DP[24] = bhabha_vtx_err[0];
									DP[25] = bhabha_vtx_err[1];
									DP[26] = bhabha_vtx_err[2];

									for (Int_t k = 0; k < N; k++)
									{
										minimum->SetVariableValue(k, P[k]);

										minimum->SetVariableLimits(k, P[k] - 30 * DP[k], P[k] + 30 * DP[k]);

										minimum->SetVariableValue(k + N, P[k]);
										minimum->SetVariableValue(k + 2 * N, DP[k]);
									}

									for (Int_t k = 0; k < M; k++)
									{
										minimum->SetVariableValue(k + 3 * N, 1.);
									}

									minimum->SetVariableStepSize(24, 0.);
									minimum->SetVariableStepSize(25, 0.);
									minimum->SetVariableStepSize(26, 0.);

									// do the minimization
									isConverged = minimum->Minimize();

									if (isConverged == 0)
									{
										min_value = minimum->MinValue();
										for (Int_t m = 0; m < N; m++)
										{
											P1[m] = minimum->X()[m];
											CHISQR += pow((P1[m] - P[m]) / DP[m], 2);
											for (Int_t n = 0; n < N; n++)
												chi2_corr->SetBinContent(m + 1, n + 1, minimum->Correlation(m, n));
										}
									}
									else
									{
										for (Int_t m = 0; m < N; m++)
										{
											CHISQR = 999999;
											min_value = 999999;
											P1[m] = -999.;
										}
									}

									if (min_value < min_value_def)
									{
										min_value_def = min_value;
										CHISQRMIN = CHISQR;
										found_best = 1;

										// Setting clusters for a solution
										for (Int_t l = 0; l < 4; l++)
										{
											clusters_min[l][0] = P1[l * 5];
											clusters_min[l][1] = P1[l * 5 + 1];
											clusters_min[l][2] = P1[l * 5 + 2];
											clusters_min[l][3] = P1[l * 5 + 3];
											clusters_min[l][4] = P1[l * 5 + 4];

											R.SetClu(l, clusters_min[l][0],
															 clusters_min[l][1],
															 clusters_min[l][2],
															 clusters_min[l][3],
															 clusters_min[l][4]);
										}

										R.SetClu(4, 0., 0., 0., 0., 0.);
										R.SetClu(5, 0., 0., 0., 0., 0.);

										S = R.MySolve(selected);

										neu_vtx_min[0][0] = S.sol[0][0];
										neu_vtx_min[0][1] = S.sol[0][1];
										neu_vtx_min[0][2] = S.sol[0][2];
										neu_vtx_min[0][3] = S.sol[0][3];

										neu_vtx_min[1][0] = S.sol[1][0];
										neu_vtx_min[1][1] = S.sol[1][1];
										neu_vtx_min[1][2] = S.sol[1][2];
										neu_vtx_min[1][3] = S.sol[1][3];

										for (Int_t k = 0; k < 2; k++)
										{
											for (Int_t l = 0; l < 4; l++)
											{
												neutral_mom(clusters_min[l][0], clusters_min[l][1], clusters_min[l][2], clusters_min[l][4], neu_vtx_min[k], gamma_mom[k][l]);
											}

											mom_kaon[k][0] = gamma_mom[k][0][0] + gamma_mom[k][1][0] + gamma_mom[k][2][0] + gamma_mom[k][3][0];
											mom_kaon[k][1] = gamma_mom[k][0][1] + gamma_mom[k][1][1] + gamma_mom[k][2][1] + gamma_mom[k][3][1];
											mom_kaon[k][2] = gamma_mom[k][0][2] + gamma_mom[k][1][2] + gamma_mom[k][2][2] + gamma_mom[k][3][2];
											mom_kaon[k][3] = gamma_mom[k][0][3] + gamma_mom[k][1][3] + gamma_mom[k][2][3] + gamma_mom[k][3][3];

											v_kaon[k] = cVel * sqrt(pow(mom_kaon[k][0], 2) + pow(mom_kaon[k][1], 2) + pow(mom_kaon[k][2], 2)) / mom_kaon[k][3];

											bhabha_mom_fit[0] = P1[20];
											bhabha_mom_fit[1] = P1[21];
											bhabha_mom_fit[2] = P1[22];
											bhabha_mom_fit[3] = P1[23];

											plane_intersection(bhabha_vtx, y_axis, neu_vtx_min[k], mom_kaon[k], ip_tri[k]); //! Plane rec

											ip_tri[k][0] = bhabha_vtx[0];
											ip_tri[k][1] = bhabha_vtx[1];

											if (abs(ip_tri[k][2] - bhabha_vtx[2]) > 2)
												ip_tri[k][2] = bhabha_vtx[2];

											length_kaon[k] = sqrt(pow(neu_vtx_min[k][0] - ip_tri[k][0], 2) + pow(neu_vtx_min[k][1] - ip_tri[k][1], 2) + pow(neu_vtx_min[k][2] - ip_tri[k][2], 2));

											diff_kaon[k] = v_kaon[k] * neu_vtx_min[k][3] - length_kaon[k];
										}

										if (abs(diff_kaon[0]) < abs(diff_kaon[1]))
										{
											gamma_mom_final[0][0] = gamma_mom[0][0][0];
											gamma_mom_final[0][1] = gamma_mom[0][0][1];
											gamma_mom_final[0][2] = gamma_mom[0][0][2];
											gamma_mom_final[0][3] = gamma_mom[0][0][3];
											gamma_mom_final[0][4] = clusters_min[0][0];
											gamma_mom_final[0][5] = clusters_min[0][1];
											gamma_mom_final[0][6] = clusters_min[0][2];
											gamma_mom_final[0][7] = clusters_min[0][3];

											gamma_mom_final[1][0] = gamma_mom[0][1][0];
											gamma_mom_final[1][1] = gamma_mom[0][1][1];
											gamma_mom_final[1][2] = gamma_mom[0][1][2];
											gamma_mom_final[1][3] = gamma_mom[0][1][3];
											gamma_mom_final[1][4] = clusters_min[1][0];
											gamma_mom_final[1][5] = clusters_min[1][1];
											gamma_mom_final[1][6] = clusters_min[1][2];
											gamma_mom_final[1][7] = clusters_min[1][3];

											gamma_mom_final[2][0] = gamma_mom[0][2][0];
											gamma_mom_final[2][1] = gamma_mom[0][2][1];
											gamma_mom_final[2][2] = gamma_mom[0][2][2];
											gamma_mom_final[2][3] = gamma_mom[0][2][3];
											gamma_mom_final[2][4] = clusters_min[2][0];
											gamma_mom_final[2][5] = clusters_min[2][1];
											gamma_mom_final[2][6] = clusters_min[2][2];
											gamma_mom_final[2][7] = clusters_min[2][3];

											gamma_mom_final[3][0] = gamma_mom[0][3][0];
											gamma_mom_final[3][1] = gamma_mom[0][3][1];
											gamma_mom_final[3][2] = gamma_mom[0][3][2];
											gamma_mom_final[3][3] = gamma_mom[0][3][3];
											gamma_mom_final[3][4] = clusters_min[3][0];
											gamma_mom_final[3][5] = clusters_min[3][1];
											gamma_mom_final[3][6] = clusters_min[3][2];
											gamma_mom_final[3][7] = clusters_min[3][3];

											fourKnetri_kinfit[0] = gamma_mom_final[0][0] + gamma_mom_final[1][0] +
																						 gamma_mom_final[2][0] + gamma_mom_final[3][0];
											fourKnetri_kinfit[1] = gamma_mom_final[0][1] + gamma_mom_final[1][1] +
																						 gamma_mom_final[2][1] + gamma_mom_final[3][1];
											fourKnetri_kinfit[2] = gamma_mom_final[0][2] + gamma_mom_final[1][2] +
																						 gamma_mom_final[2][2] + gamma_mom_final[3][2];
											fourKnetri_kinfit[3] = gamma_mom_final[0][3] + gamma_mom_final[1][3] +
																						 gamma_mom_final[2][3] + gamma_mom_final[3][3];
											fourKnetri_kinfit[4] = sqrt(pow(fourKnetri_kinfit[0], 2) +
																									pow(fourKnetri_kinfit[1], 2) +
																									pow(fourKnetri_kinfit[2], 2));
											fourKnetri_kinfit[5] = sqrt(pow(fourKnetri_kinfit[3], 2) - pow(fourKnetri_kinfit[4], 2));
											fourKnetri_kinfit[6] = neu_vtx_min[0][0];
											fourKnetri_kinfit[7] = neu_vtx_min[0][1];
											fourKnetri_kinfit[8] = neu_vtx_min[0][2];
											fourKnetri_kinfit[9] = neu_vtx_min[0][3];

											g4takentri_kinfit[0] = ind_gam[0];
											g4takentri_kinfit[1] = ind_gam[1];
											g4takentri_kinfit[2] = ind_gam[2];
											g4takentri_kinfit[3] = ind_gam[3];

											plane_intersection(bhabha_vtx, y_axis, neu_vtx_min[0], fourKnetri_kinfit, iptri_kinfit); //! Plane rec

											iptri_kinfit[0] = bhabha_vtx[0];
											iptri_kinfit[1] = bhabha_vtx[1];

											if (abs(iptri_kinfit[2] - bhabha_vtx[2]) > 2)
												iptri_kinfit[2] = bhabha_vtx[2];
										}
										else if (abs(diff_kaon[0]) > abs(diff_kaon[1]))
										{
											gamma_mom_final[0][0] = gamma_mom[1][0][0];
											gamma_mom_final[0][1] = gamma_mom[1][0][1];
											gamma_mom_final[0][2] = gamma_mom[1][0][2];
											gamma_mom_final[0][3] = gamma_mom[1][0][3];
											gamma_mom_final[0][4] = clusters_min[0][0];
											gamma_mom_final[0][5] = clusters_min[0][1];
											gamma_mom_final[0][6] = clusters_min[0][2];
											gamma_mom_final[0][7] = clusters_min[0][3];

											gamma_mom_final[1][0] = gamma_mom[1][1][0];
											gamma_mom_final[1][1] = gamma_mom[1][1][1];
											gamma_mom_final[1][2] = gamma_mom[1][1][2];
											gamma_mom_final[1][3] = gamma_mom[1][1][3];
											gamma_mom_final[1][4] = clusters_min[1][0];
											gamma_mom_final[1][5] = clusters_min[1][1];
											gamma_mom_final[1][6] = clusters_min[1][2];
											gamma_mom_final[1][7] = clusters_min[1][3];

											gamma_mom_final[2][0] = gamma_mom[1][2][0];
											gamma_mom_final[2][1] = gamma_mom[1][2][1];
											gamma_mom_final[2][2] = gamma_mom[1][2][2];
											gamma_mom_final[2][3] = gamma_mom[1][2][3];
											gamma_mom_final[2][4] = clusters_min[2][0];
											gamma_mom_final[2][5] = clusters_min[2][1];
											gamma_mom_final[2][6] = clusters_min[2][2];
											gamma_mom_final[2][7] = clusters_min[2][3];

											gamma_mom_final[3][0] = gamma_mom[1][3][0];
											gamma_mom_final[3][1] = gamma_mom[1][3][1];
											gamma_mom_final[3][2] = gamma_mom[1][3][2];
											gamma_mom_final[3][3] = gamma_mom[1][3][3];
											gamma_mom_final[3][4] = clusters_min[3][0];
											gamma_mom_final[3][5] = clusters_min[3][1];
											gamma_mom_final[3][6] = clusters_min[3][2];
											gamma_mom_final[3][7] = clusters_min[3][3];

											fourKnetri_kinfit[0] = gamma_mom_final[0][0] + gamma_mom_final[1][0] +
																						 gamma_mom_final[2][0] + gamma_mom_final[3][0];
											fourKnetri_kinfit[1] = gamma_mom_final[0][1] + gamma_mom_final[1][1] +
																						 gamma_mom_final[2][1] + gamma_mom_final[3][1];
											fourKnetri_kinfit[2] = gamma_mom_final[0][2] + gamma_mom_final[1][2] +
																						 gamma_mom_final[2][2] + gamma_mom_final[3][2];
											fourKnetri_kinfit[3] = gamma_mom_final[0][3] + gamma_mom_final[1][3] +
																						 gamma_mom_final[2][3] + gamma_mom_final[3][3];
											fourKnetri_kinfit[4] = sqrt(pow(fourKnetri_kinfit[2], 2));
											fourKnetri_kinfit[5] = sqrt(pow(fourKnetri_kinfit[3], 2) - pow(fourKnetri_kinfit[4], 2));
											fourKnetri_kinfit[6] = neu_vtx_min[1][0];
											fourKnetri_kinfit[7] = neu_vtx_min[1][1];
											fourKnetri_kinfit[8] = neu_vtx_min[1][2];
											fourKnetri_kinfit[9] = neu_vtx_min[1][3];

											g4takentri_kinfit[0] = ind_gam[0];
											g4takentri_kinfit[1] = ind_gam[1];
											g4takentri_kinfit[2] = ind_gam[2];
											g4takentri_kinfit[3] = ind_gam[3];

											plane_intersection(bhabha_vtx, y_axis, neu_vtx_min[1], fourKnetri_kinfit, iptri_kinfit); //! Plane rec

											iptri_kinfit[0] = bhabha_vtx[0];
											iptri_kinfit[1] = bhabha_vtx[1];

											if (abs(iptri_kinfit[2] - bhabha_vtx[2]) > 2)
												iptri_kinfit[2] = bhabha_vtx[2];
										}
									}

									delete minimum;
								}
								///////////////////////////////////////////////////////////////////
							}
						}

			if (found_best == 1)
			{
				chi2_hist->Fill(CHISQRMIN);
				std::cout << CHISQRMIN / 9. << std::endl;
			}
		}

		if (found_best == 0)
		{
			gamma_mom_final[0][0] = -999.;
			gamma_mom_final[0][1] = -999.;
			gamma_mom_final[0][2] = -999.;
			gamma_mom_final[0][3] = -999.;
			gamma_mom_final[0][4] = -999.;
			gamma_mom_final[0][5] = -999.;
			gamma_mom_final[0][6] = -999.;
			gamma_mom_final[0][7] = -999.;

			gamma_mom_final[1][0] = -999.;
			gamma_mom_final[1][1] = -999.;
			gamma_mom_final[1][2] = -999.;
			gamma_mom_final[1][3] = -999.;
			gamma_mom_final[1][4] = -999.;
			gamma_mom_final[1][5] = -999.;
			gamma_mom_final[1][6] = -999.;
			gamma_mom_final[1][7] = -999.;

			gamma_mom_final[2][0] = -999.;
			gamma_mom_final[2][1] = -999.;
			gamma_mom_final[2][2] = -999.;
			gamma_mom_final[2][3] = -999.;
			gamma_mom_final[2][4] = -999.;
			gamma_mom_final[2][5] = -999.;
			gamma_mom_final[2][6] = -999.;
			gamma_mom_final[2][7] = -999.;

			gamma_mom_final[3][0] = -999.;
			gamma_mom_final[3][1] = -999.;
			gamma_mom_final[3][2] = -999.;
			gamma_mom_final[3][3] = -999.;
			gamma_mom_final[3][4] = -999.;
			gamma_mom_final[3][5] = -999.;
			gamma_mom_final[3][6] = -999.;
			gamma_mom_final[3][7] = -999.;

			fourKnetri_kinfit[0] = -999.;
			fourKnetri_kinfit[1] = -999.;
			fourKnetri_kinfit[2] = -999.;
			fourKnetri_kinfit[3] = -999.;
			fourKnetri_kinfit[4] = -999.;
			fourKnetri_kinfit[5] = -999.;
			fourKnetri_kinfit[6] = -999.;
			fourKnetri_kinfit[7] = -999.;
			fourKnetri_kinfit[8] = -999.;
			fourKnetri_kinfit[9] = -999.;

			g4takentri_kinfit[0] = -999;
			g4takentri_kinfit[1] = -999;
			g4takentri_kinfit[2] = -999;
			g4takentri_kinfit[3] = -999;

			iptri_kinfit[0] = -999.;
			iptri_kinfit[1] = -999.;
			iptri_kinfit[2] = -999.;
		}

		tree->Fill();
	}

	TCanvas *c1 = new TCanvas("c1", "", 750, 750);

	name = "#chi^{2}(" + std::to_string(M - 4) + ")";

	chi2_hist->GetXaxis()->SetTitle(name);
	chi2_hist->GetYaxis()->SetTitle("Counts");
	chi2_hist->GetYaxis()->SetRangeUser(0, 1.2 * chi2_hist->GetMaximum());
	chi2_hist->Draw();

	c1->Print("chi2_test.png");

	TCanvas *c2 = new TCanvas("c2", "", 750, 750);
	c1->SetRightMargin(0.15);
	gStyle->SetOptStat(0);

	c2->SetLogz(0);
	chi2_corr->GetXaxis()->SetTitle("");
	chi2_corr->GetYaxis()->SetTitle("");
	chi2_corr->Draw("COLZ");

	c2->Print("chi2_corr.png");

	tree->Print();

	file->Write();
	file->Close();
	delete file;
}
