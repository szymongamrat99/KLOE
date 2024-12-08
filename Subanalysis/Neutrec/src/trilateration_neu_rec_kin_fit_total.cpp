#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TError.h"

#include "../../../Include/Codes/reconstructor.h"
#include "../../../Include/const.h"
#include "../../../Include/Codes/uncertainties.h"
#include "../../../Include/Codes/charged_mom.h"
#include "../../../Include/Codes/neutral_mom.h"
#include "../../../Include/Codes/lorentz_transf.h"
#include "../../../Include/Codes/plane_intersection.h"
#include "../../../Include/Codes/closest_approach.h"
#include "../../../Include/Codes/constraints_tri.h"
#include "../../../Include/Codes/chi2_dist.h"

#include "../inc/trilateration.hpp"

using namespace std;

const Int_t N_free = 24, N_const = 4, M = 5, jmin = 0, jmax = 0;
Int_t j_ch, k_ch;
const Float_t Trf = 2.715; // ns - time of a bunch (correction)

const Int_t loopcount = 10;

TF1 *constraints[M];

Double_t det;
Float_t CHISQR, CHISQRTMP, FUNVAL, FUNVALTMP, FUNVALMIN, T0;
Int_t fail;

Int_t selected[4] = {1, 2, 3, 4};

void tri_neurec_kinfit_corr(Short_t ind_data_mc, Int_t first_file, Int_t last_file, Short_t loopcount, Short_t jmin, Short_t jmax)
{

	gErrorIgnoreLevel = 6001;

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	TString name = "";

	const Int_t range = (jmax - jmin) + 1;

	name = "neuvtx_tri_kin_fit_" + std::to_string(first_file) + "_" + std::to_string(last_file) + "_" + loopcount + "_" + M + "_" + range + "_" + ind_data_mc + "_" + "new" + ".root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_tri_kin_fit", "Neu vtx rec with trilateration kin fit");

	// TFile *file_corr = new TFile("bunch_corr.root");
	// TTree *tree_corr = (TTree *)file_corr->Get("h_bunch_corr");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_mom_err[4], bhabha_vtx[3];

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

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth, mcflag;
	Float_t cluster[5][500], Kchboost[9], Knerec[9], Knemc[9], ipmc[3], ip[3], Dtmc, bunch_corr;

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster[0]);
	chain->SetBranchAddress("Ycl", cluster[1]);
	chain->SetBranchAddress("Zcl", cluster[2]);
	chain->SetBranchAddress("Tcl", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);
	chain->SetBranchAddress("mcflag", &mcflag);

	// tree_corr->SetBranchAddress("bunchcorr", &bunch_corr);
	// chain->AddFriend(tree_corr);

	Int_t nentries = (Int_t)chain->GetEntries();

	Bool_t clusterEnergy, solError, cond_clus[4];
	Int_t ind_gam[4], sort_index[range], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged, n_bunch;
	Float_t CHISQRMIN, min_value_def, value_first[range];
	Double_t reinitialize[(N_free + N_const) * (N_free + N_const)] = {0.};
	Float_t gamma_mom_min[4][4], neu_vtx_min[4], kaon_mom_min[4], kaon_vel[3], kaon_vel_tot, kaon_path_tot, bhabha_vtx_min[3], ip_min[3], time_diff[2][4], time_diff_fin, gamma_path[2][4], neu_vtx[2][4];

	TMatrixD V(N_free + N_const, N_free + N_const), D(M, N_free + N_const), D_T(N_free + N_const, M), V_final(N_free + N_const, N_free + N_const), V_aux(N_free + N_const, N_free + N_const), V_min(N_free + N_const, N_free + N_const), Aux(M, M), V_invert(N_free, N_free), V_init(N_free + N_const, N_free + N_const);
	TVectorD X(N_free + N_const), C(M), X_final(N_free + N_const), L(M), CORR(N_free + N_const), X_init(N_free + N_const), X_min(N_free + N_const), C_min(M), L_min(M), C_aux(M), L_aux(M), X_init_min(N_free + N_const), X_init_aux(N_free + N_const);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Float_t length_ch, time_ch, velocity_ch, gamma_mom_final[4][8], fourKnetri_kinfit[10], iptri_kinfit[3], y_axis[3], ip_tri[2][3], bhabha_mom_fit[4], gamma_mom_tmp[4][4], fourKnetri_tmp[2][4], value[range], kaon_vel_tmp[2], dist_tmp[2], ip_tmp[2][3], gamma_len[4];
	Int_t g4takentri_kinfit[4], bunchnum;

	TBranch *b_gamma1tri = tree->Branch("fourgamma1tri_kinfit", gamma_mom_final[0], "fourgamma1tri_kinfit[8]/F");
	TBranch *b_gamma2tri = tree->Branch("fourgamma2tri_kinfit", gamma_mom_final[1], "fourgamma2tri_kinfit[8]/F");
	TBranch *b_gamma3tri = tree->Branch("fourgamma3tri_kinfit", gamma_mom_final[2], "fourgamma3tri_kinfit[8]/F");
	TBranch *b_gamma4tri = tree->Branch("fourgamma4tri_kinfit", gamma_mom_final[3], "fourgamma4tri_kinfit[8]/F");
	TBranch *b_iptri = tree->Branch("iptri_kinfit", iptri_kinfit, "iptri_kinfit[3]/F");
	TBranch *b_Knetri = tree->Branch("fourKnetri_kinfit", fourKnetri_kinfit, "fourKnetri_kinfit[10]/F");
	TBranch *b_done = tree->Branch("done4_kinfit", &isConverged, "done4_kinfit/I");
	TBranch *b_fourg4taken = tree->Branch("g4takentri_kinfit", g4takentri_kinfit, "g4takentri_kinfit[4]/I");
	TBranch *b_bunchnum = tree->Branch("bunchnum", &bunchnum, "bunchnum/I");

	TBranch *b_Lmin = tree->Branch("lag_mult", "TVectorD", &L_min);
	TBranch *b_Cmin = tree->Branch("const_min", "TVectorD", &C_min);
	TBranch *b_Xmin = tree->Branch("min_vars", "TVectorD", &X_min);
	TBranch *b_Vmin = tree->Branch("min_cov", "TMatrixD", &V_min);
	TBranch *b_Xinit = tree->Branch("init_vars", "TVectorD", &X_init_min);

	TBranch *b_chisqr = tree->Branch("chi2min", &CHISQRMIN, "chi2min/F");

	constraints[0] = new TF1("Ene consv", &ene_consv, 0, 1, N_free + N_const);
	constraints[1] = new TF1("Minv consv", &minv_consv, 0, 1, N_free + N_const);
	constraints[2] = new TF1("x consv", &x_consv, 0, 1, N_free + N_const);
	constraints[3] = new TF1("y consv", &y_consv, 0, 1, N_free + N_const);
	constraints[4] = new TF1("z consv", &z_consv, 0, 1, N_free + N_const);
	constraints[5] = new TF1("gamma1 consv", &gamma1_consv, 0, 1, N_free + N_const);
	constraints[6] = new TF1("gamma2 consv", &gamma2_consv, 0, 1, N_free + N_const);
	constraints[7] = new TF1("gamma3 consv", &gamma3_consv, 0, 1, N_free + N_const);
	constraints[8] = new TF1("gamma4 consv", &gamma4_consv, 0, 1, N_free + N_const);

	TH1 *chi2 = new TH1F("chi2", "", 100, -10.0, 30.0);

	Bool_t data_flag;
	Bool_t cond_time_clus[2];

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		min_value_def = 999999.;
		FUNVALMIN = 999999.;
		CHISQRMIN = 999999.;

		isConverged = 0;

		if (ind_data_mc == 0)
			data_flag = (mctruth == 6);//(mctruth == 1 || mctruth == 2));
		else if (ind_data_mc == 1)
			data_flag = (mcflag == 1 && mctruth != 0);
		else if (ind_data_mc == 2)
			data_flag = (mcflag == 1 && mctruth != 0);

		if (nclu >= 4 && data_flag)
		{

			std::cout << int(mctruth) << std::endl;
			std::cout << 100 * i / (Float_t)nentries << "% done" << std::endl;

			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							V.SetMatrixArray(reinitialize);

							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;

							Reconstructor R;
							Solution S;

							for (Int_t k = 0; k < 4; k++)
							{
								R.SetClu(k, cluster[0][ind_gam[k]],
												 cluster[1][ind_gam[k]],
												 cluster[2][ind_gam[k]],
												 cluster[3][ind_gam[k]],
												 cluster[4][ind_gam[k]]);

								R.SetClu(4, 0., 0., 0., 0., 0.);
								R.SetClu(5, 0., 0., 0., 0., 0.);
							}

							S = R.MySolve(selected);

							cond_time_clus[0] = S.sol[0][3] < cluster[3][ind_gam[0]] && S.sol[0][3] < cluster[3][ind_gam[1]] && S.sol[0][3] < cluster[3][ind_gam[2]] && S.sol[0][3] < cluster[3][ind_gam[3]];

							cond_time_clus[1] = S.sol[1][3] < cluster[3][ind_gam[0]] && S.sol[1][3] < cluster[3][ind_gam[1]] && S.sol[1][3] < cluster[3][ind_gam[2]] && S.sol[1][3] < cluster[3][ind_gam[3]];

							clusterEnergy = (cluster[4][ind_gam[0]] > MIN_CLU_ENE && cluster[4][ind_gam[1]] > MIN_CLU_ENE && cluster[4][ind_gam[2]] > MIN_CLU_ENE && cluster[4][ind_gam[3]] > MIN_CLU_ENE);

							cond_clus[0] = cluster[3][ind_gam[0]] > 0 && cluster[0][ind_gam[0]] != 0 && cluster[1][ind_gam[0]] != 0 && cluster[2][ind_gam[0]] != 0;
							cond_clus[1] = cluster[3][ind_gam[1]] > 0 && cluster[0][ind_gam[1]] != 0 && cluster[1][ind_gam[1]] != 0 && cluster[2][ind_gam[1]] != 0;
							cond_clus[2] = cluster[3][ind_gam[2]] > 0 && cluster[0][ind_gam[2]] != 0 && cluster[1][ind_gam[2]] != 0 && cluster[2][ind_gam[2]] != 0;
							cond_clus[3] = cluster[3][ind_gam[3]] > 0 && cluster[0][ind_gam[3]] != 0 && cluster[1][ind_gam[3]] != 0 && cluster[2][ind_gam[3]] != 0;

							if (clusterEnergy && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && (cond_time_clus[0] || cond_time_clus[1]))
							{
								for (Int_t k = 0; k < 4; k++)
								{
									X_init(k * 5) = cluster[0][ind_gam[k]];
									X_init(k * 5 + 1) = cluster[1][ind_gam[k]];
									X_init(k * 5 + 2) = cluster[2][ind_gam[k]];
									X_init(k * 5 + 3) = cluster[3][ind_gam[k]];
									X_init(k * 5 + 4) = cluster[4][ind_gam[k]];

									V(k * 5, k * 5) = pow(clu_x_error(X_init(k * 5), X_init(k * 5 + 1), X_init(k * 5 + 2), X_init(k * 5 + 4)), 2);
									V(k * 5 + 1, k * 5 + 1) = pow(clu_y_error(X_init(k * 5), X_init(k * 5 + 1), X_init(k * 5 + 2), X_init(k * 5 + 4)), 2);
									V(k * 5 + 2, k * 5 + 2) = pow(clu_z_error(X_init(k * 5), X_init(k * 5 + 1), X_init(k * 5 + 2), X_init(k * 5 + 4)), 2); // cm
									V(k * 5 + 3, k * 5 + 3) = pow(clu_time_error(X_init(k * 5 + 4)), 2);																									 // ns
									V(k * 5 + 4, k * 5 + 4) = pow(clu_ene_error(X_init(k * 5 + 4)), 2);																										 // MeV
								}

								X_init(20) = bhabha_mom[0];
								X_init(21) = bhabha_mom[1];
								X_init(22) = bhabha_mom[2];
								X_init(23) = bhabha_mom[3];

								V(20, 20) = pow(bhabha_mom_err[0], 2);
								V(21, 21) = pow(bhabha_mom_err[1], 2);
								V(22, 22) = pow(bhabha_mom_err[2], 2);
								V(23, 23) = pow(bhabha_mom_err[3], 2);

								X_init(24) = bhabha_vtx[0];
								X_init(25) = bhabha_vtx[1];
								X_init(26) = bhabha_vtx[2];

								V(24, 24) = 0.;
								V(25, 25) = 0.;
								V(26, 26) = 0.;

								for (Int_t k1 = jmin; k1 <= jmax; k1++)
								{
									fail = 0;
									CHISQR = 999999.;
									CHISQRTMP = 999999.;
									FUNVALTMP = 999999.;

									X_init(27) = k1;

									X = X_init;

									V_invert = V.GetSub(0, N_free - 1, 0, N_free - 1).Invert();

									for (Int_t j = 0; j < loopcount; j++)
									{

										for (Int_t l = 0; l < M; l++)
										{
											C(l) = constraints[l]->EvalPar(0, X.GetMatrixArray());
											for (Int_t m = 0; m < N_free + N_const; m++)
											{
												constraints[l]->SetParameters(X.GetMatrixArray());
												if (m < N_free)
													D(l, m) = constraints[l]->GradientPar(m, 0, 0.01);
												// else if (m == 24)
												//	D(l, m) = constraints[l]->GradientPar(m, 0, 1.0);
												else
													D(l, m) = 0;
											}
										}

										D_T.Transpose(D);

										Aux = (D * V * D_T);

										Aux.Invert(&det);

										if (!TMath::IsNaN(det) && det != 0)
										{

											L = (Aux * C);

											CORR = V * D_T * L;

											if (1)
											{
												X_final = X - CORR;
												V_final = V - V * D_T * Aux * D * V;

												for (Int_t l = 0; l < M; l++)
													C(l) = constraints[l]->EvalPar(0, X_final.GetMatrixArray());

												FUNVAL = Dot((X_final - X_init).GetSub(0, N_free - 1), V_invert * (X_final - X_init).GetSub(0, N_free - 1)) + Dot(L, C);

												CHISQR = Dot((X_final - X_init).GetSub(0, N_free - 1), V_invert * (X_final - X_init).GetSub(0, N_free - 1));
											}
											else
											{
												break;
											}
										}
										else
										{
											fail = 1;
										}

										if (fail == 0)
										{
											X = X_final;
											X_init_aux = X_init;
											V_aux = V_final;
											L_aux = L;
											C_aux = C;
											FUNVALTMP = FUNVAL;
											CHISQRTMP = CHISQR;
										}
									}

									T0 = X(27) * Trf;

									Reconstructor R;
									Solution S;

									for (Int_t k = 0; k < 4; k++)
									{
										R.SetClu(k, X[k * 5],
														 X[k * 5 + 1],
														 X[k * 5 + 2],
														 X[k * 5 + 3] + T0,
														 X[k * 5 + 4]);

										R.SetClu(4, 0., 0., 0., 0., 0.);
										R.SetClu(5, 0., 0., 0., 0., 0.);
									}

									S = R.MySolve(selected);

									Float_t distance[4] = {0.};

									for (Int_t k = 0; k < 2; k++)
									{
										if (!S.error[k])
										{
											neu_vtx[k][0] = S.sol[k][0];
											neu_vtx[k][1] = S.sol[k][1];
											neu_vtx[k][2] = S.sol[k][2];
											neu_vtx[k][3] = S.sol[k][3];
										}
										else
										{
											neu_vtx[k][0] = 999.;
											neu_vtx[k][1] = 999.;
											neu_vtx[k][2] = 999.;
											neu_vtx[k][3] = 999.;
										}

										for (Int_t l = 0; l < 4; l++)
										{
											distance[l] = sqrt(pow(X[l * 5] - neu_vtx[k][0], 2) +
																				 pow(X[l * 5 + 1] - neu_vtx[k][1], 2) +
																				 pow(X[l * 5 + 2] - neu_vtx[k][2], 2));

											gamma_mom_tmp[l][0] = X[l * 5 + 4] * ((X[l * 5] - neu_vtx[k][0]) / distance[l]);
											gamma_mom_tmp[l][1] = X[l * 5 + 4] * ((X[l * 5 + 1] - neu_vtx[k][1]) / distance[l]);
											gamma_mom_tmp[l][2] = X[l * 5 + 4] * ((X[l * 5 + 2] - neu_vtx[k][2]) / distance[l]);
											gamma_mom_tmp[l][3] = X[l * 5 + 4];

											gamma_mom_tmp[l][4] = X[l * 5];
											gamma_mom_tmp[l][5] = X[l * 5 + 1];
											gamma_mom_tmp[l][6] = X[l * 5 + 2];
											gamma_mom_tmp[l][7] = X[l * 5 + 3] + T0;
										}

										fourKnetri_tmp[k][0] = gamma_mom_tmp[0][0] + gamma_mom_tmp[1][0] + gamma_mom_tmp[2][0] + gamma_mom_tmp[3][0];
										fourKnetri_tmp[k][1] = gamma_mom_tmp[0][1] + gamma_mom_tmp[1][1] + gamma_mom_tmp[2][1] + gamma_mom_tmp[3][1];
										fourKnetri_tmp[k][2] = gamma_mom_tmp[0][2] + gamma_mom_tmp[1][2] + gamma_mom_tmp[2][2] + gamma_mom_tmp[3][2];
										fourKnetri_tmp[k][3] = gamma_mom_tmp[0][3] + gamma_mom_tmp[1][3] + gamma_mom_tmp[2][3] + gamma_mom_tmp[3][3];

										fourKnetri_tmp[k][4] = sqrt(pow(fourKnetri_tmp[k][0], 2) + pow(fourKnetri_tmp[k][1], 2) + pow(fourKnetri_tmp[k][2], 2));
										fourKnetri_tmp[k][5] = sqrt(pow(fourKnetri_tmp[k][3], 2) - pow(fourKnetri_tmp[k][4], 2));

										kaon_vel_tmp[k] = cVel * fourKnetri_tmp[k][4] / fourKnetri_tmp[k][3];

										y_axis[0] = 0.;
										y_axis[1] = X[21];
										y_axis[2] = 0.;

										plane_intersection(bhabha_vtx, y_axis, neu_vtx[k], fourKnetri_tmp[k], ip_tmp[k]);

										ip_tmp[k][0] = bhabha_vtx[0];
										ip_tmp[k][1] = bhabha_vtx[1];
										if (abs(ip_tmp[k][2] - bhabha_vtx[2]) > 2)
											ip_tmp[k][2] = bhabha_vtx[2];

										dist_tmp[k] = sqrt(pow(neu_vtx[k][0] - ip_tmp[k][0], 2) +
																			 pow(neu_vtx[k][1] - ip_tmp[k][1], 2) +
																			 pow(neu_vtx[k][2] - ip_tmp[k][2], 2));

										value[k] = sqrt(pow(neu_vtx[k][3] - (dist_tmp[k] / kaon_vel_tmp[k]), 2) + pow(fourKnetri_tmp[k][5] - mK0, 2));

										if (TMath::IsNaN(value[k]))
											value[k] = 999999.;
									}

									cond_time_clus[0] = S.sol[0][3] < X(3) + T0 && S.sol[0][3] < X(8) + T0 && S.sol[0][3] < X(13) + T0 && S.sol[0][3] < X(18) + T0;

									cond_time_clus[1] = S.sol[1][3] < X(3) + T0 && S.sol[1][3] < X(8) + T0 && S.sol[1][3] < X(13) + T0 && S.sol[1][3] < X(18) + T0;

									if (abs(CHISQRTMP) < abs(CHISQRMIN))
									{
										if (cond_time_clus[0] && value[0] < value[1])
										{
											isConverged = 1;
											FUNVALMIN = FUNVALTMP;
											CHISQRMIN = CHISQRTMP;

											X_min = X;
											X_init_min = X_init_aux;
											V_min = V_aux;
											V_init = V;
											C_min = C_aux;
											L_min = L_aux;

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
												distance[l] = sqrt(pow(X[l * 5] - neu_vtx[0][0], 2) +
																					 pow(X[l * 5 + 1] - neu_vtx[0][1], 2) +
																					 pow(X[l * 5 + 2] - neu_vtx[0][2], 2));

												gamma_mom_final[l][0] = X[l * 5 + 4] * ((X[l * 5] - neu_vtx[0][0]) / distance[l]);
												gamma_mom_final[l][1] = X[l * 5 + 4] * ((X[l * 5 + 1] - neu_vtx[0][1]) / distance[l]);
												gamma_mom_final[l][2] = X[l * 5 + 4] * ((X[l * 5 + 2] - neu_vtx[0][2]) / distance[l]);
												gamma_mom_final[l][3] = X[l * 5 + 4];
												gamma_mom_final[l][4] = X[l * 5];
												gamma_mom_final[l][5] = X[l * 5 + 1];
												gamma_mom_final[l][6] = X[l * 5 + 2];
												gamma_mom_final[l][7] = X[l * 5 + 3] + T0;
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

											bunchnum = X_min(27);
										}
										else if (cond_time_clus[1] && value[1] < value[0])
										{
											isConverged = 1;
											FUNVALMIN = FUNVALTMP;
											CHISQRMIN = CHISQRTMP;

											X_min = X;
											X_init_min = X_init_aux;
											V_min = V_aux;
											V_init = V;
											C_min = C_aux;
											L_min = L_aux;

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
												distance[l] = sqrt(pow(X[l * 5] - neu_vtx[1][0], 2) +
																					 pow(X[l * 5 + 1] - neu_vtx[1][1], 2) +
																					 pow(X[l * 5 + 2] - neu_vtx[1][2], 2));

												gamma_mom_final[l][0] = X[l * 5 + 4] * ((X[l * 5] - neu_vtx[1][0]) / distance[l]);
												gamma_mom_final[l][1] = X[l * 5 + 4] * ((X[l * 5 + 1] - neu_vtx[1][1]) / distance[l]);
												gamma_mom_final[l][2] = X[l * 5 + 4] * ((X[l * 5 + 2] - neu_vtx[1][2]) / distance[l]);
												gamma_mom_final[l][3] = X[l * 5 + 4];
												gamma_mom_final[l][4] = X[l * 5];
												gamma_mom_final[l][5] = X[l * 5 + 1];
												gamma_mom_final[l][6] = X[l * 5 + 2];
												gamma_mom_final[l][7] = X[l * 5 + 3] + T0;
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

											bunchnum = X_min(27);
										}
										else
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
						}
			chi2->Fill(CHISQRMIN);
		}
		else
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

		tree->Fill();
	}

	gStyle->SetOptStat("iMr");
	TF1 *func = new TF1("chi2dist", chi2dist, 0, 100, 2);

	func->SetParameters(3000.0, (Double_t)M);
	func->SetParNames("Norm", "Degrees of Freedom");

	func->SetParLimits(0, 0.001, 5000.0);
	func->SetParLimits(1, 1.0, 10.0);

	TCanvas *c1 = new TCanvas("c1", "", 750, 750);
	// chi2->Fit(func);
	//c1->SetLogy(1);

	chi2->GetYaxis()->SetRangeUser(0, 1.2 * chi2->GetMaximum());
	chi2->Draw();
	c1->Print("test_chi2.png");

	tree->Print();

	file->Write();
	file->Close();
	delete file;
}
