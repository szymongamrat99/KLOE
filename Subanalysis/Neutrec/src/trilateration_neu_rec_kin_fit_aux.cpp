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

#include <boost/progress.hpp>

#include "reconstructor.h"
#include "const.h"
#include "uncertainties.h"
#include "charged_mom.h"
#include "neutral_mom.h"
#include "plane_intersection.h"
#include "closest_approach.h"
#include "chi2_dist.h"
#include <KinFitter.h>

#include "../inc/trilateration.hpp"

using namespace std;

Int_t TrilaterationNeurecKinfit(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	Bool_t
			good_clus = (Bool_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["goodClus"];

	const Short_t
			loopcount = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["loopCount"],
			jmin = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
			jmax = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["bunchMax"],
			M = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["numOfConstraints"],
			N_const = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["fixedVars"],
			N_free = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["freeVars"],
			range = Int_t(jmax - jmin) + 1;

	const Double_t
			chiSqrStep = (Double_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["chiSqrStep"];

	Double_t det;
	Float_t CHISQR, CHISQRTMP, FUNVAL, FUNVALTMP, FUNVALMIN, Tcorr;
	Int_t fail;

	Int_t selected[4] = {1, 2, 3, 4};

	TF1 *constraints[M];

	gErrorIgnoreLevel = 6001;

	// Creation of filename for the analysis step
	std::string
			datestamp = Obj.getCurrentDate(),
			name = "";

	name = Paths::neutrec_dir + Paths::root_files_dir + Filenames::neu_trilateration_kin_fit_filename + datestamp + "_" + std::to_string(N_free) + "_" + std::to_string(N_const) + "_" + std::to_string(M) + "_" + std::to_string(loopcount) + "_" + int(dataType) + Paths::ext_root;

	Utils::properties["variables"]["tree"]["filename"]["trilaterationKinFit"] = name;

	TFile *file = new TFile(name.c_str(), "recreate");
	TTree *tree = new TTree(Filenames::neutrec_tri_tree, "Trilateration reconstruction with kin fit");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_mom_err[4], bhabha_vtx[3];

	chain.SetBranchAddress("Bpx", &bhabha_mom[0]);
	chain.SetBranchAddress("Bpy", &bhabha_mom[1]);
	chain.SetBranchAddress("Bpz", &bhabha_mom[2]);
	chain.SetBranchAddress("Broots", &bhabha_mom[3]);

	chain.SetBranchAddress("Bwidpx", &bhabha_mom_err[0]);
	chain.SetBranchAddress("Bwidpy", &bhabha_mom_err[1]);
	chain.SetBranchAddress("Bwidpz", &bhabha_mom_err[2]);
	chain.SetBranchAddress("Brootserr", &bhabha_mom_err[3]);

	chain.SetBranchAddress("Bx", &bhabha_vtx[0]);
	chain.SetBranchAddress("By", &bhabha_vtx[1]);
	chain.SetBranchAddress("Bz", &bhabha_vtx[2]);

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth, mcflag;
	Float_t cluster[5][500], Kchboost[9], Knerec[9], KnemcOld[9], ipmcOld[3], ip[3], Dtmc, bunch_corr;

	KLOE::BaseKinematics baseKin;

	chain.SetBranchAddress("nclu", &nclu);
	chain.SetBranchAddress("Xcl", cluster[0]);
	chain.SetBranchAddress("Ycl", cluster[1]);
	chain.SetBranchAddress("Zcl", cluster[2]);
	chain.SetBranchAddress("TclOld", cluster[3]);
	chain.SetBranchAddress("Enecl", cluster[4]);

	chain.SetBranchAddress("mctruth", &mctruth);
	chain.SetBranchAddress("mcflag", &mcflag);
	chain.SetBranchAddress("ncll", baseKin.ncll.data());

	Int_t nentries = (Int_t)chain.GetEntries();

	Bool_t clusterEnergy, solError, cond_clus[4];
	Int_t ind_gam[4], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged, n_bunch;
	Float_t CHISQRMIN, min_value_def;
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

	TH1 *chi2 = new TH1F("chi2", "", 100, -10.0, 10.0);

	Bool_t data_flag;
	Bool_t cond_time_clus[2];

	// Progress bar
	boost::progress_display show_progress(nentries);
	// --------------------------------------------------------------------

	Int_t mode = 1;

	std::vector<std::string> ConstSet = {
			"EnergyConsvCM",
			"MinvConsv",
			"NeutralXPathConsvLAB",
			"NeutralYPathConsvLAB",
			"NeutralZPathConsvLAB"};

	KLOE::KinFitter *kinematicFitObj = new KLOE::KinFitter("Trilateration", N_free, N_const, M, 0, loopcount, chiSqrStep, logger);
	kinematicFitObj->ConstraintSet(ConstSet);

	Double_t
			Param[N_free + N_const],
			Errors[N_free + N_const];

	for (Int_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		min_value_def = 999999.;
		FUNVALMIN = 999999.;
		CHISQRMIN = 999999.;

		isConverged = 0;

		int
				mctruth_int = int(mctruth),
				mcflag_int = int(mcflag);

		Obj.dataFlagSetter(dataType, data_flag, mcflag_int, mctruth_int);

		if (nclu >= 4 && data_flag)
		{
			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
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

							clusterEnergy = cluster[4][baseKin.ncll[ind_gam[0]] - 1] > KLOE::MIN_CLU_ENE &&
															cluster[4][baseKin.ncll[ind_gam[1]] - 1] > KLOE::MIN_CLU_ENE &&
															cluster[4][baseKin.ncll[ind_gam[2]] - 1] > KLOE::MIN_CLU_ENE &&
															cluster[4][baseKin.ncll[ind_gam[3]] - 1] > KLOE::MIN_CLU_ENE;

							for (Int_t k = 0; k < 4; k++)
							{
								cond_clus[k] =
										cluster[3][baseKin.ncll[ind_gam[k]] - 1] > 0 &&
										cluster[0][baseKin.ncll[ind_gam[k]] - 1] != 0 &&
										cluster[1][baseKin.ncll[ind_gam[k]] - 1] != 0 &&
										cluster[2][baseKin.ncll[ind_gam[k]] - 1] != 0;

								if (k < 2)
									cond_time_clus[k] = S.sol[k][3] < cluster[3][baseKin.ncll[ind_gam[0]] - 1] &&
																			S.sol[k][3] < cluster[3][baseKin.ncll[ind_gam[1]] - 1] &&
																			S.sol[k][3] < cluster[3][baseKin.ncll[ind_gam[2]] - 1] &&
																			S.sol[k][3] < cluster[3][baseKin.ncll[ind_gam[3]] - 1];
							}

							Bool_t cond_tot = cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && clusterEnergy;

							if (cond_tot && (cond_time_clus[0] || cond_time_clus[1]))
							{
								for (Int_t k = 0; k < 4; k++)
								{
									Param[k * 5] = cluster[0][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 1] = cluster[1][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 2] = cluster[2][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 3] = cluster[3][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 4] = cluster[4][baseKin.ncll[ind_gam[k]] - 1];

									Errors[k * 5] = clu_x_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]);
									Errors[k * 5 + 1] = clu_y_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]);
									Errors[k * 5 + 2] = clu_z_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]);
									// cm
									Errors[k * 5 + 3] = clu_time_error(Param[k * 5 + 4]); // ns
									Errors[k * 5 + 4] = clu_ene_error(Param[k * 5 + 4]);	 // MeV

									Param[20 + k] = bhabha_mom[k];
									Errors[20 + k] = bhabha_mom_err[k];

									if (k < 3)
									{
										Param[24 + k] = bhabha_vtx[k];
										Errors[24 + k] = 0.;
									}
								}

								for (Int_t k1 = jmin; k1 <= jmax; k1++)
								{
									kinematicFitObj->ParameterInitialization(Param, Errors);

									Tcorr = k1 * KLOE::T0;

									CHISQRTMP = kinematicFitObj->FitFunction(Tcorr);

									kinematicFitObj->GetResults(X, V, X_init, V_init);

									Reconstructor R;
									Solution S;

									for (Int_t k = 0; k < 4; k++)
									{
										R.SetClu(k, X[k * 5],
														 X[k * 5 + 1],
														 X[k * 5 + 2],
														 X[k * 5 + 3],
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
											neutral_mom(X[l * 5], X[l * 5 + 1], X[l * 5 + 2], X[l * 5 + 4], neu_vtx[k], gamma_mom_tmp[l]);

											gamma_mom_tmp[l][4] = X[l * 5];
											gamma_mom_tmp[l][5] = X[l * 5 + 1];
											gamma_mom_tmp[l][6] = X[l * 5 + 2];
											gamma_mom_tmp[l][7] = X[l * 5 + 3];
										}

										fourKnetri_tmp[k][0] = gamma_mom_tmp[0][0] + gamma_mom_tmp[1][0] + gamma_mom_tmp[2][0] + gamma_mom_tmp[3][0];
										fourKnetri_tmp[k][1] = gamma_mom_tmp[0][1] + gamma_mom_tmp[1][1] + gamma_mom_tmp[2][1] + gamma_mom_tmp[3][1];
										fourKnetri_tmp[k][2] = gamma_mom_tmp[0][2] + gamma_mom_tmp[1][2] + gamma_mom_tmp[2][2] + gamma_mom_tmp[3][2];
										fourKnetri_tmp[k][3] = gamma_mom_tmp[0][3] + gamma_mom_tmp[1][3] + gamma_mom_tmp[2][3] + gamma_mom_tmp[3][3];

										fourKnetri_tmp[k][4] = sqrt(pow(fourKnetri_tmp[k][0], 2) + pow(fourKnetri_tmp[k][1], 2) + pow(fourKnetri_tmp[k][2], 2));
										fourKnetri_tmp[k][5] = sqrt(pow(fourKnetri_tmp[k][3], 2) - pow(fourKnetri_tmp[k][4], 2));

										kaon_vel_tmp[k] = PhysicsConstants::cVel * fourKnetri_tmp[k][4] / fourKnetri_tmp[k][3];

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

										value[k] = sqrt(pow(neu_vtx[k][3] - (dist_tmp[k] / kaon_vel_tmp[k]), 2) + pow(fourKnetri_tmp[k][5] - PhysicsConstants::mK0, 2));

										if (TMath::IsNaN(value[k]))
											value[k] = 999999.;
									}

									cond_time_clus[0] = S.sol[0][3] < X(3) &&
																			S.sol[0][3] < X(8) &&
																			S.sol[0][3] < X(13) &&
																			S.sol[0][3] < X(18);

									cond_time_clus[1] = S.sol[1][3] < X(3) &&
																			S.sol[1][3] < X(8) &&
																			S.sol[1][3] < X(13) &&
																			S.sol[1][3] < X(18);

									if (abs(CHISQRTMP) < abs(CHISQRMIN))
									{										
										if (cond_time_clus[0] && value[0] < value[1])
										{
											isConverged = 1;
											FUNVALMIN = FUNVALTMP;
											CHISQRMIN = CHISQRTMP;

											kinematicFitObj->GetResults(X_min, V_min, X_init_min, V_init);

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
												gamma_mom_final[l][7] = X[l * 5 + 3];
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

											kinematicFitObj->GetResults(X_min, V_min, X_init_min, V_init);

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
												gamma_mom_final[l][7] = X[l * 5 + 3];
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

			if(isConverged == 1)
			{
				chi2->Fill(CHISQRMIN);
			}
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

		++show_progress;
	}

	gStyle->SetOptStat("iMr");
	TF1 *func = new TF1("chi2dist", chi2dist, 0, 100, 2);

	func->SetParameters(3000.0, (Double_t)M);
	func->SetParNames("Norm", "Degrees of Freedom");

	func->SetParLimits(0, 0.001, 5000.0);
	func->SetParLimits(1, 1.0, 10.0);

	TCanvas *c1 = new TCanvas("c1", "", 750, 750);
	// chi2->Fit(func);
	// c1->SetLogy(1);

	chi2->GetYaxis()->SetRangeUser(0, 1.2 * chi2->GetMaximum());
	chi2->Draw();
	c1->Print("test_chi2.png");

	tree->Print();

	file->Write();
	file->Close();
	delete file;

	return 0;
}
