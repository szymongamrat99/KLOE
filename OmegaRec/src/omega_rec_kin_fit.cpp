#include <string.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <Math/Functor.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>

#include <reconstructor.h>
#include "const.h"
#include "uncertainties.h"
#include "charged_mom.h"
#include "neutral_mom.h"
#include "lorentz_transf.h"
#include "plane_intersection.h"
#include "closest_approach.h"
#include "constraints_omega.h"
#include "chi2_dist.h"
#include <pi0_photon_pair.h>

#include "../inc/omegarec.hpp"

using namespace std;

const Int_t N_free = 27, N_const = 4, M = 5, jmin = 0, jmax = 0;
Int_t j_ch, k_ch;

const Int_t loopcount = 10;

TF1 *constraints[M];

Double_t det;
Float_t CHISQR, CHISQRTMP, FUNVAL, FUNVALTMP, FUNVALMIN, Tcorr;
Int_t fail;

Int_t selected[4] = {1, 2, 3, 4};

int omegarec(Int_t first_file, Int_t last_file, Short_t loopcount, Short_t jmin, Short_t jmax, Controls::DataType data_type)
{

	gErrorIgnoreLevel = 6001;

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	TString name = "";

	const Int_t range = (jmax - jmin) + 1;

	name = omegarec_dir + root_files_dir + omega_rec_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree(omegarec_tree, "Omega reconstruction with kin fit");

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

	BaseKinematics baseKin;

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster[0]);
	chain->SetBranchAddress("Ycl", cluster[1]);
	chain->SetBranchAddress("Zcl", cluster[2]);
	chain->SetBranchAddress("Tcl", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);
	chain->SetBranchAddress("mcflag", &mcflag);
	chain->SetBranchAddress("ncll", baseKin.ncll);

	// Charged tracks momenta

	chain->SetBranchAddress("trk1", baseKin.trk[0]);
	chain->SetBranchAddress("trk2", baseKin.trk[1]);
	chain->SetBranchAddress("Kchrec", baseKin.Kchrec);

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

	Float_t length_ch, time_ch, velocity_ch, gamma_mom_final[4][8], Omegarec_kinfit[10], Omegapi0_kinfit[10], pi0_kinfit[10], iptri_kinfit[3], y_axis[3], ip_tri[2][3], bhabha_mom_fit[4], gamma_mom_tmp[4][8], Omegarec_tmp[2][10], Omegapi0_tmp[2][10], pi0_tmp[2][10], value[range], kaon_vel_tmp[2], dist_tmp[2], ip_tmp[2][3], gamma_len[4];
	Int_t g4takentri_kinfit[4], bunchnum;

	Int_t
			clusterind[4],
			clusterindpi0[2][2];

	Float_t
			gamma_mom[4][4],
			gamma_mom_pi0[2][2][4],
			pi0_mom[2][4],
			omega_mom[4];

	TBranch *b_gamma1tri = tree->Branch("gamma1tri_kinfit", gamma_mom_final[0], "gamma1tri_kinfit[8]/F");
	TBranch *b_gamma2tri = tree->Branch("gamma2tri_kinfit", gamma_mom_final[1], "gamma2tri_kinfit[8]/F");
	TBranch *b_gamma3tri = tree->Branch("gamma3tri_kinfit", gamma_mom_final[2], "gamma1tri_kinfit[8]/F");
	TBranch *b_gamma4tri = tree->Branch("gamma4tri_kinfit", gamma_mom_final[3], "gamma2tri_kinfit[8]/F");

	TBranch *b_pi0omega = tree->Branch("omegapi0tri_kinfit", Omegapi0_kinfit, "omegapi0tri_kinfit[10]/F");
	TBranch *b_pi0 = tree->Branch("pi0_kinfit", pi0_kinfit, "pi0_kinfit[10]/F");

	TBranch *b_omegatri = tree->Branch("omega_kinfit", Omegarec_kinfit, "omega_kinfit[10]/F");

	TBranch *b_iptri = tree->Branch("iptri_kinfit", iptri_kinfit, "iptri_kinfit[3]/F");
	TBranch *b_done = tree->Branch("done4_kinfit", &isConverged, "done4_kinfit/I");
	TBranch *b_fourg4taken = tree->Branch("g4takentri_kinfit", g4takentri_kinfit, "g4takentri_kinfit[4]/I");
	TBranch *b_bunchnum = tree->Branch("bunchnum", &bunchnum, "bunchnum/I");

	TBranch *b_Lmin = tree->Branch("lag_mult", "TVectorD", &L_min);
	TBranch *b_Cmin = tree->Branch("const_min", "TVectorD", &C_min);
	TBranch *b_Xmin = tree->Branch("min_vars", "TVectorD", &X_min);
	TBranch *b_Vmin = tree->Branch("min_cov", "TMatrixD", &V_min);
	TBranch *b_Xinit = tree->Branch("init_vars", "TVectorD", &X_init_min);

	TBranch *b_chisqr = tree->Branch("chi2min", &CHISQRMIN, "chi2min/F");

	//constraints[0] = new TF1("Ene consv", &OmegaConstraints::ene_consv, 0, 1, N_free + N_const);
	//constraints[1] = new TF1("Px consv", &OmegaConstraints::px_consv, 0, 1, N_free + N_const);
	//constraints[2] = new TF1("Py consv", &OmegaConstraints::py_consv, 0, 1, N_free + N_const);
	//constraints[3] = new TF1("Pz consv", &OmegaConstraints::pz_consv, 0, 1, N_free + N_const);
	//constraints[4] = new TF1("Minv consv", &OmegaConstraints::minv_omega_consv, 0, 1, N_free + N_const);
	constraints[0] = new TF1("Minv pi01 consv", &OmegaConstraints::minv_pi01_consv, 0, 1, N_free + N_const);
	constraints[1] = new TF1("Minv pi02 consv", &OmegaConstraints::minv_pi02_consv, 0, 1, N_free + N_const);
	constraints[2] = new TF1("x consv", &OmegaConstraints::x_consv, 0, 1, N_free + N_const);
	constraints[3] = new TF1("y consv", &OmegaConstraints::y_consv, 0, 1, N_free + N_const);
	constraints[4] = new TF1("z consv", &OmegaConstraints::z_consv, 0, 1, N_free + N_const);

	TH1 *chi2 = new TH1F("chi2", "", 100, -10.0, 30.0);

	TH1 *pi01 = new TH1F("pi01", "", 100, 600.0, 1000.0);
	TH1 *pi02 = new TH1F("pi02", "", 100, 0.0, 200.0);

	Bool_t data_flag;
	Bool_t cond_time_clus[2];

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		min_value_def = 999999.;
		FUNVALMIN = 999999.;
		CHISQRMIN = 999999.;

		isConverged = 0;

		dataFlagSetter(data_type, data_flag, int(mcflag), int(mctruth));

		if (nclu >= 4 && data_flag)
		{
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
								R.SetClu(k, cluster[0][baseKin.ncll[ind_gam[k]] - 1],
												 cluster[1][baseKin.ncll[ind_gam[k]] - 1],
												 cluster[2][baseKin.ncll[ind_gam[k]] - 1],
												 cluster[3][baseKin.ncll[ind_gam[k]] - 1],
												 cluster[4][baseKin.ncll[ind_gam[k]] - 1]);

								R.SetClu(4, 0., 0., 0., 0., 0.);
								R.SetClu(5, 0., 0., 0., 0., 0.);
							}

							S = R.MySolve(selected);

							cond_time_clus[0] = S.sol[0][3] < cluster[3][baseKin.ncll[ind_gam[0]] - 1] && S.sol[0][3] < cluster[3][baseKin.ncll[ind_gam[1]] - 1] && S.sol[0][3] < cluster[3][baseKin.ncll[ind_gam[2]] - 1] && S.sol[0][3] < cluster[3][baseKin.ncll[ind_gam[3]] - 1];

							cond_time_clus[1] = S.sol[1][3] < cluster[3][baseKin.ncll[ind_gam[0]] - 1] && S.sol[1][3] < cluster[3][baseKin.ncll[ind_gam[1]] - 1] && S.sol[1][3] < cluster[3][baseKin.ncll[ind_gam[2]] - 1] && S.sol[1][3] < cluster[3][baseKin.ncll[ind_gam[3]] - 1];

							clusterEnergy = (cluster[4][baseKin.ncll[ind_gam[0]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[1]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[2]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[3]] - 1] > MIN_CLU_ENE);

							cond_clus[0] = cluster[3][baseKin.ncll[ind_gam[0]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[0]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[0]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[0]] - 1] != 0;
							cond_clus[1] = cluster[3][baseKin.ncll[ind_gam[1]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[1]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[1]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[1]] - 1] != 0;
							cond_clus[2] = cluster[3][baseKin.ncll[ind_gam[2]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[2]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[2]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[2]] - 1] != 0;
							cond_clus[3] = cluster[3][baseKin.ncll[ind_gam[3]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[3]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[3]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[3]] - 1] != 0;

							if (clusterEnergy && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && (cond_time_clus[0] || cond_time_clus[1]))
							{
								for (Int_t k = 0; k < 4; k++)
								{
									X_init(k * 5) = cluster[0][baseKin.ncll[ind_gam[k]] - 1];
									X_init(k * 5 + 1) = cluster[1][baseKin.ncll[ind_gam[k]] - 1];
									X_init(k * 5 + 2) = cluster[2][baseKin.ncll[ind_gam[k]] - 1];
									X_init(k * 5 + 3) = cluster[3][baseKin.ncll[ind_gam[k]] - 1];
									X_init(k * 5 + 4) = cluster[4][baseKin.ncll[ind_gam[k]] - 1];

									V(k * 5, k * 5) = pow(clu_x_error(X_init(k * 5), X_init(k * 5 + 1), X_init(k * 5 + 2), X_init(k * 5 + 4)), 2);
									V(k * 5 + 1, k * 5 + 1) = pow(clu_y_error(X_init(k * 5), X_init(k * 5 + 1), X_init(k * 5 + 2), X_init(k * 5 + 4)), 2);
									V(k * 5 + 2, k * 5 + 2) = pow(clu_z_error(X_init(k * 5), X_init(k * 5 + 1), X_init(k * 5 + 2), X_init(k * 5 + 4)), 2); // cm
									V(k * 5 + 3, k * 5 + 3) = pow(clu_time_error(X_init(k * 5 + 4)), 2);																									 // ns
									V(k * 5 + 4, k * 5 + 4) = pow(clu_ene_error(X_init(k * 5 + 4)), 2);																										 // MeV
								}

								X_init(20) = baseKin.Kchrec[6];
								X_init(21) = baseKin.Kchrec[7];
								X_init(22) = baseKin.Kchrec[8];

								V(20, 20) = pow(2.0, 2);
								V(21, 21) = pow(1.97, 2);
								V(22, 22) = pow(4.72, 2);

								X_init(23) = bhabha_mom[0];
								X_init(24) = bhabha_mom[1];
								X_init(25) = bhabha_mom[2];
								X_init(26) = bhabha_mom[3];

								V(23, 23) = pow(bhabha_mom_err[0], 2);
								V(24, 24) = pow(bhabha_mom_err[1], 2);
								V(25, 25) = pow(bhabha_mom_err[2], 2);
								V(26, 26) = pow(bhabha_mom_err[3], 2);

								X_init(27) = bhabha_vtx[0];
								X_init(28) = bhabha_vtx[1];
								X_init(29) = bhabha_vtx[2];

								V(27, 27) = 0.;
								V(28, 28) = 0.;
								V(29, 29) = 0.;

								for (Int_t k1 = jmin; k1 <= jmax; k1++)
								{
									fail = 0;
									CHISQR = 999999.;
									CHISQRTMP = 999999.;
									FUNVALTMP = 999999.;

									X_init(30) = k1;

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

									Tcorr = X(30) * T0;

									Reconstructor R;
									Solution S;

									for (Int_t k = 0; k < 4; k++)
									{
										R.SetClu(k, X[k * 5],
														 X[k * 5 + 1],
														 X[k * 5 + 2],
														 X[k * 5 + 3] + Tcorr,
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

											gamma_mom[l][0] = X[l * 5 + 4] * ((X[l * 5] - neu_vtx[k][0]) / distance[l]);
											gamma_mom[l][1] = X[l * 5 + 4] * ((X[l * 5 + 1] - neu_vtx[k][1]) / distance[l]);
											gamma_mom[l][2] = X[l * 5 + 4] * ((X[l * 5 + 2] - neu_vtx[k][2]) / distance[l]);
											gamma_mom[l][3] = X[l * 5 + 4];
										}

										Pi0PhotonPair(ind_gam, gamma_mom, clusterindpi0, gamma_mom_pi0, pi0_mom, true, baseKin.trk, omega_mom);

										for (Int_t l = 0; l < 4; l++)
										{
											gamma_mom_tmp[l][0] = gamma_mom[l][0];
											gamma_mom_tmp[l][1] = gamma_mom[l][1];
											gamma_mom_tmp[l][2] = gamma_mom[l][2];
											gamma_mom_tmp[l][3] = gamma_mom[l][3];

											gamma_mom_tmp[l][4] = X[l * 5];
											gamma_mom_tmp[l][5] = X[l * 5 + 1];
											gamma_mom_tmp[l][6] = X[l * 5 + 2];
											gamma_mom_tmp[l][7] = X[l * 5 + 3];
										}

										Omegapi0_tmp[k][0] = pi0_mom[0][0];
										Omegapi0_tmp[k][1] = pi0_mom[0][1];
										Omegapi0_tmp[k][2] = pi0_mom[0][2];
										Omegapi0_tmp[k][3] = pi0_mom[0][3];

										Omegapi0_tmp[k][4] = sqrt(pow(Omegapi0_tmp[k][0], 2) + pow(Omegapi0_tmp[k][1], 2) + pow(Omegapi0_tmp[k][2], 2));
										Omegapi0_tmp[k][5] = sqrt(pow(Omegapi0_tmp[k][3], 2) - pow(Omegapi0_tmp[k][4], 2));

										Omegapi0_tmp[k][6] = neu_vtx[k][0];
										Omegapi0_tmp[k][7] = neu_vtx[k][1];
										Omegapi0_tmp[k][8] = neu_vtx[k][2];
										Omegapi0_tmp[k][9] = neu_vtx[k][3];

										pi0_tmp[k][0] = pi0_mom[1][0];
										pi0_tmp[k][1] = pi0_mom[1][1];
										pi0_tmp[k][2] = pi0_mom[1][2];
										pi0_tmp[k][3] = pi0_mom[1][3];

										pi0_tmp[k][4] = sqrt(pow(pi0_tmp[k][0], 2) + pow(pi0_tmp[k][1], 2) + pow(pi0_tmp[k][2], 2));
										pi0_tmp[k][5] = sqrt(pow(pi0_tmp[k][3], 2) - pow(pi0_tmp[k][4], 2));

										pi0_tmp[k][6] = neu_vtx[k][0];
										pi0_tmp[k][7] = neu_vtx[k][1];
										pi0_tmp[k][8] = neu_vtx[k][2];
										pi0_tmp[k][9] = neu_vtx[k][3];

										Omegarec_tmp[k][0] = omega_mom[0];
										Omegarec_tmp[k][1] = omega_mom[1];
										Omegarec_tmp[k][2] = omega_mom[2];
										Omegarec_tmp[k][3] = omega_mom[3];

										Omegarec_tmp[k][4] = sqrt(pow(Omegarec_tmp[k][0], 2) + pow(Omegarec_tmp[k][1], 2) + pow(Omegarec_tmp[k][2], 2));
										Omegarec_tmp[k][5] = sqrt(pow(Omegarec_tmp[k][3], 2) - pow(Omegarec_tmp[k][4], 2));

										Omegarec_tmp[k][6] = neu_vtx[k][0];
										Omegarec_tmp[k][7] = neu_vtx[k][1];
										Omegarec_tmp[k][8] = neu_vtx[k][2];
										Omegarec_tmp[k][9] = neu_vtx[k][3];

										ip_tmp[k][0] = neu_vtx[k][0];
										ip_tmp[k][1] = neu_vtx[k][1];
										ip_tmp[k][2] = neu_vtx[k][2];

										value[k] = sqrt(pow(Omegarec_tmp[k][5] - mOmega, 2) +
																		pow(Omegapi0_tmp[k][5] - mPi0, 2) +
																		pow(pi0_tmp[k][5] - mPi0, 2));

										if (TMath::IsNaN(value[k]))
											value[k] = 999999.;
									}

									cond_time_clus[0] = S.sol[0][3] < X(3) + Tcorr && S.sol[0][3] < X(8) + Tcorr && S.sol[0][3] < X(13) + Tcorr && S.sol[0][3] < X(18) + Tcorr;

									cond_time_clus[1] = S.sol[1][3] < X(3) + Tcorr && S.sol[1][3] < X(8) + Tcorr && S.sol[1][3] < X(13) + Tcorr && S.sol[1][3] < X(18) + Tcorr;

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
												gamma_mom_final[l][7] = X[l * 5 + 3] + Tcorr;
											}

											Omegapi0_kinfit[0] = Omegapi0_tmp[0][0];
											Omegapi0_kinfit[1] = Omegapi0_tmp[0][1];
											Omegapi0_kinfit[2] = Omegapi0_tmp[0][2];
											Omegapi0_kinfit[3] = Omegapi0_tmp[0][3];
											Omegapi0_kinfit[4] = Omegapi0_tmp[0][4];
											Omegapi0_kinfit[5] = Omegapi0_tmp[0][5];
											Omegapi0_kinfit[6] = neu_vtx_min[0];
											Omegapi0_kinfit[7] = neu_vtx_min[1];
											Omegapi0_kinfit[8] = neu_vtx_min[2];
											Omegapi0_kinfit[9] = neu_vtx_min[3];

											pi0_kinfit[0] = pi0_tmp[0][0];
											pi0_kinfit[1] = pi0_tmp[0][1];
											pi0_kinfit[2] = pi0_tmp[0][2];
											pi0_kinfit[3] = pi0_tmp[0][3];
											pi0_kinfit[4] = pi0_tmp[0][4];
											pi0_kinfit[5] = pi0_tmp[0][5];
											pi0_kinfit[6] = neu_vtx_min[0];
											pi0_kinfit[7] = neu_vtx_min[1];
											pi0_kinfit[8] = neu_vtx_min[2];
											pi0_kinfit[9] = neu_vtx_min[3];

											Omegarec_kinfit[0] = Omegarec_tmp[0][0];
											Omegarec_kinfit[1] = Omegarec_tmp[0][1];
											Omegarec_kinfit[2] = Omegarec_tmp[0][2];
											Omegarec_kinfit[3] = Omegarec_tmp[0][3];
											Omegarec_kinfit[4] = Omegarec_tmp[0][4];
											Omegarec_kinfit[5] = Omegarec_tmp[0][5];
											Omegarec_kinfit[6] = neu_vtx_min[0];
											Omegarec_kinfit[7] = neu_vtx_min[1];
											Omegarec_kinfit[8] = neu_vtx_min[2];
											Omegarec_kinfit[9] = neu_vtx_min[3];

											iptri_kinfit[0] = ip_tmp[0][0];
											iptri_kinfit[1] = ip_tmp[0][1];
											iptri_kinfit[2] = ip_tmp[0][2];

											bunchnum = X_min(33);
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
												gamma_mom_final[l][7] = X[l * 5 + 3] + Tcorr;
											}

											Omegapi0_kinfit[0] = Omegapi0_tmp[1][0];
											Omegapi0_kinfit[1] = Omegapi0_tmp[1][1];
											Omegapi0_kinfit[2] = Omegapi0_tmp[1][2];
											Omegapi0_kinfit[3] = Omegapi0_tmp[1][3];
											Omegapi0_kinfit[4] = Omegapi0_tmp[1][4];
											Omegapi0_kinfit[5] = Omegapi0_tmp[1][5];
											Omegapi0_kinfit[6] = neu_vtx_min[0];
											Omegapi0_kinfit[7] = neu_vtx_min[1];
											Omegapi0_kinfit[8] = neu_vtx_min[2];
											Omegapi0_kinfit[9] = neu_vtx_min[3];

											pi0_kinfit[0] = pi0_tmp[1][0];
											pi0_kinfit[1] = pi0_tmp[1][1];
											pi0_kinfit[2] = pi0_tmp[1][2];
											pi0_kinfit[3] = pi0_tmp[1][3];
											pi0_kinfit[4] = pi0_tmp[1][4];
											pi0_kinfit[5] = pi0_tmp[1][5];
											pi0_kinfit[6] = neu_vtx_min[0];
											pi0_kinfit[7] = neu_vtx_min[1];
											pi0_kinfit[8] = neu_vtx_min[2];
											pi0_kinfit[9] = neu_vtx_min[3];

											Omegarec_kinfit[0] = Omegarec_tmp[1][0];
											Omegarec_kinfit[1] = Omegarec_tmp[1][1];
											Omegarec_kinfit[2] = Omegarec_tmp[1][2];
											Omegarec_kinfit[3] = Omegarec_tmp[1][3];
											Omegarec_kinfit[4] = Omegarec_tmp[1][4];
											Omegarec_kinfit[5] = Omegarec_tmp[1][5];
											Omegarec_kinfit[6] = neu_vtx_min[0];
											Omegarec_kinfit[7] = neu_vtx_min[1];
											Omegarec_kinfit[8] = neu_vtx_min[2];
											Omegarec_kinfit[9] = neu_vtx_min[3];

											iptri_kinfit[0] = ip_tmp[1][0];
											iptri_kinfit[1] = ip_tmp[1][1];
											iptri_kinfit[2] = ip_tmp[1][2];

											bunchnum = X_min(33);
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

											Omegarec_kinfit[0] = 999.;
											Omegarec_kinfit[1] = 999.;
											Omegarec_kinfit[2] = 999.;
											Omegarec_kinfit[3] = 999.;
											Omegarec_kinfit[4] = 999.;
											Omegarec_kinfit[5] = 999.;
											Omegarec_kinfit[6] = 999.;
											Omegarec_kinfit[7] = 999.;
											Omegarec_kinfit[8] = 999.;
											Omegarec_kinfit[9] = 999.;

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
			pi01->Fill(Omegarec_kinfit[5]);
			pi02->Fill(pi0_kinfit[5]);
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

			Omegapi0_kinfit[0] = 999.;
			Omegapi0_kinfit[1] = 999.;
			Omegapi0_kinfit[2] = 999.;
			Omegapi0_kinfit[3] = 999.;
			Omegapi0_kinfit[4] = 999.;
			Omegapi0_kinfit[5] = 999.;
			Omegapi0_kinfit[6] = 999.;
			Omegapi0_kinfit[7] = 999.;
			Omegapi0_kinfit[8] = 999.;
			Omegapi0_kinfit[9] = 999.;

			pi0_kinfit[0] = 999.;
			pi0_kinfit[1] = 999.;
			pi0_kinfit[2] = 999.;
			pi0_kinfit[3] = 999.;
			pi0_kinfit[4] = 999.;
			pi0_kinfit[5] = 999.;
			pi0_kinfit[6] = 999.;
			pi0_kinfit[7] = 999.;
			pi0_kinfit[8] = 999.;
			pi0_kinfit[9] = 999.;

			Omegarec_kinfit[0] = 999.;
			Omegarec_kinfit[1] = 999.;
			Omegarec_kinfit[2] = 999.;
			Omegarec_kinfit[3] = 999.;
			Omegarec_kinfit[4] = 999.;
			Omegarec_kinfit[5] = 999.;
			Omegarec_kinfit[6] = 999.;
			Omegarec_kinfit[7] = 999.;
			Omegarec_kinfit[8] = 999.;
			Omegarec_kinfit[9] = 999.;

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

	chi2->GetYaxis()->SetRangeUser(0, 1.2 * chi2->GetMaximum());
	chi2->Draw();
	c1->Print("test_chi2.png");

	TCanvas *c2 = new TCanvas("c2", "", 750, 750);

	pi01->GetYaxis()->SetRangeUser(0, 1.2 * pi01->GetMaximum());
	pi01->Draw();
	c2->Print("pi01.png");

	TCanvas *c3 = new TCanvas("c3", "", 750, 750);

	pi02->GetYaxis()->SetRangeUser(0, 1.2 * pi02->GetMaximum());
	pi02->Draw();
	c3->Print("pi02.png");

	tree->Print();

	file->Write();
	file->Close();
	delete file;

	return 0;
}
