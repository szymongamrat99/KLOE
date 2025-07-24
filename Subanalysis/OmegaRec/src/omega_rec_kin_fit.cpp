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
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>

#include <reconstructor.h>
#include "uncertainties.h"
#include "charged_mom.h"
#include "neutral_mom.h"
#include "plane_intersection.h"
#include "closest_approach.h"
#include "chi2_dist.h"
#include <pi0_photon_pair.h>
#include <KinFitter.h>

#include "../inc/omegarec.hpp"

using namespace std;

Int_t j_ch, k_ch;

Double_t det;
Float_t CHISQR, CHISQRTMP, FUNVAL, FUNVALTMP, FUNVALMIN, Tcorr;
Int_t fail;

Int_t selected[4] = {1, 2, 3, 4};

int omegarec_kin_fit(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	const int
			loopcount = properties["variables"]["KinFit"]["Omega"]["loopCount"],
			M = properties["variables"]["KinFit"]["Omega"]["numOfConstraints"],
			N_free = properties["variables"]["KinFit"]["Omega"]["freeVars"],
			N_const = properties["variables"]["KinFit"]["Omega"]["fixedVars"],
			chiSqrStep = properties["variables"]["KinFit"]["Omega"]["chiSqrStep"];

	TF1
			*constraints[M];

	gErrorIgnoreLevel = 6001;

	// Creation of filename for the analysis step
	std::string
			datestamp = Obj.getCurrentDate(),
			name = "";

	name = omegarec_dir + root_files_dir + omega_rec_kin_fit_filename + datestamp + "_" + std::to_string(N_free) + "_" + std::to_string(N_const) + "_" + std::to_string(M) + "_" + std::to_string(loopcount) + "_" + int(dataType) + ext_root;

	properties["variables"]["tree"]["filename"]["omegarecKinFit"] = name;

	TFile *file = new TFile(name.c_str(), "recreate");
	TTree *tree = new TTree(omegarec_kin_fit_tree, "Omega reconstruction with kin fit");

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

	BaseKinematics baseKin;

	chain.SetBranchAddress("nclu", &nclu);
	chain.SetBranchAddress("Xcl", cluster[0]);
	chain.SetBranchAddress("Ycl", cluster[1]);
	chain.SetBranchAddress("Zcl", cluster[2]);
	chain.SetBranchAddress("TclOld", cluster[3]);
	chain.SetBranchAddress("Enecl", cluster[4]);

	chain.SetBranchAddress("mctruth", &mctruth);
	chain.SetBranchAddress("mcflag", &mcflag);
	chain.SetBranchAddress("ncll", baseKin.ncll);

	// Charged vars
	chain.SetBranchAddress("nv", &baseKin.nv);
	chain.SetBranchAddress("ntv", &baseKin.ntv);
	chain.SetBranchAddress("ivOld", baseKin.ivOld);
	chain.SetBranchAddress("CurvOld", baseKin.CurvOld);
	chain.SetBranchAddress("PhivOld", baseKin.PhivOld);
	chain.SetBranchAddress("CotvOld", baseKin.CotvOld);
	chain.SetBranchAddress("xvOld", baseKin.xvOld);
	chain.SetBranchAddress("yvOld", baseKin.yvOld);
	chain.SetBranchAddress("zvOld", baseKin.zvOld);

	Int_t nentries = (Int_t)chain.GetEntries();

	Bool_t clusterEnergy, solError, cond_clus[4];
	Int_t ind_gam[4], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged, n_bunch;
	Float_t CHISQRMIN, min_value_def;
	std::vector<Double_t>
			reinitialize((N_free + N_const) * (N_free + N_const));
	Float_t gamma_mom_min[4][4], neu_vtx_min[4], kaon_mom_min[4], kaon_vel[3], kaon_vel_tot, kaon_path_tot, bhabha_vtx_min[3], ip_min[3], time_diff[2][4], time_diff_fin, gamma_path[2][4], neu_vtx[2][4];

	TMatrixD V(N_free + N_const, N_free + N_const), D(M, N_free + N_const), D_T(N_free + N_const, M), V_final(N_free + N_const, N_free + N_const), V_aux(N_free + N_const, N_free + N_const), V_min(N_free + N_const, N_free + N_const), Aux(M, M), V_invert(N_free, N_free), V_init(N_free + N_const, N_free + N_const);
	TVectorD X(N_free + N_const), C(M), X_final(N_free + N_const), L(M), CORR(N_free + N_const), X_init(N_free + N_const), X_min(N_free + N_const), C_min(M), L_min(M), C_aux(M), L_aux(M), X_init_min(N_free + N_const), X_init_aux(N_free + N_const);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Float_t length_ch, time_ch, velocity_ch, gamma_mom_final[4][8], Omegarec_kinfit[10], Omegapi0_kinfit[10], pi0_kinfit[10], iptri_kinfit[3], y_axis[3], ip_tri[2][3], bhabha_mom_fit[4], gamma_mom_tmp[4][8], Omegarec_tmp[2][10], Omegapi0_tmp[2][10], pi0_tmp[2][10], value[2], kaon_vel_tmp[2], dist_tmp[2], ip_tmp[2][3], gamma_len[4];
	Int_t g4takentri_kinfit[4], bunchnum;

	Int_t
			clusterind[4],
			clusterindpi0[2][2];

	Float_t
			gamma_mom[4][4],
			gamma_mom_pi0[2][2][4],
			pi0_mom[2][4],
			omega_mom[4];

	Int_t
			doneOmega,
			g4takenomega[4];

	Float_t
			gammaomega[4][8],
			Omegapi0[6],
			Pi0[6],
			PichFourMom[2][4],
			Omegarec[6],
			lengthPhotonMin[4],
			lengthKch,
			rho_00,
			rho_pm_IP,
			rho_00_IP,
			rho,
			anglePi0KaonCM,
			anglePichKaonCM,
			anglePi0OmegaPhiCM,
			anglePhiOmega,
			neu_vtx_avg[3];

	TBranch *b_gamma1omega = tree->Branch("gamma1omega", gammaomega[0], "gamma1omega[8]/F");
	TBranch *b_gamma2omega = tree->Branch("gamma2omega", gammaomega[1], "gamma2omega[8]/F");
	TBranch *b_gamma3omega = tree->Branch("gamma3omega", gammaomega[2], "gamma1omega[8]/F");
	TBranch *b_gamma4omega = tree->Branch("gamma4omega", gammaomega[3], "gamma2omega[8]/F");

	TBranch *b_pich1 = tree->Branch("pich1", PichFourMom[0], "pich1[4]/F");
	TBranch *b_pich2 = tree->Branch("pich2", PichFourMom[1], "pich2[4]/F");

	TBranch *b_pi0omega = tree->Branch("omegapi0", Omegapi0, "omegapi0[6]/F");
	TBranch *b_pi0 = tree->Branch("pi0", Pi0, "pi0[6]/F");

	TBranch *b_omega = tree->Branch("omega", Omegarec, "omega[6]/F");

	TBranch *b_neuvtxavg = tree->Branch("NeuVtxAvg", neu_vtx_avg, "NeuVtxAvg[3]/F");

	TBranch *b_lengthphoton = tree->Branch("lengthphoton", lengthPhotonMin, "lengthphoton[4]/F");

	TBranch *b_lengthkch = tree->Branch("lengthKch", &lengthKch, "lengthKch/F");

	TBranch *b_done = tree->Branch("doneomega", &doneOmega, "doneomega/I");
	TBranch *b_g4takenomega = tree->Branch("g4takenomega", g4takenomega, "g4takenomega[4]/I");

	TBranch
			*b1 = tree->Branch("rho_00", &rho_00, "rho_00/F"),
			*b2 = tree->Branch("rho_00_IP", &rho_00_IP, "rho_00_IP/F"),
			*b3 = tree->Branch("rho_pm_IP", &rho_pm_IP, "rho_pm_IP/F"),
			*b4 = tree->Branch("rho", &rho, "rho/F"),
			*b5 = tree->Branch("anglePi0KaonCM", &anglePi0KaonCM, "anglePi0KaonCM/F"),
			*b6 = tree->Branch("anglePichKaonCM", &anglePichKaonCM, "anglePichKaonCM/F"),
			*b7 = tree->Branch("anglePi0OmegaPhiCM", &anglePi0OmegaPhiCM, "anglePi0OmegaPhiCM/F"),
			*b8 = tree->Branch("anglePhiOmega", &anglePhiOmega, "anglePhiOmega/F");

	TBranch *b_Lmin = tree->Branch("lag_mult", "TVectorD", &L_min);
	TBranch *b_Cmin = tree->Branch("const_min", "TVectorD", &C_min);
	TBranch *b_Xmin = tree->Branch("min_vars", "TVectorD", &X_min);
	TBranch *b_Vmin = tree->Branch("min_cov", "TMatrixD", &V_min);
	TBranch *b_Xinit = tree->Branch("init_vars", "TVectorD", &X_init_min);

	TBranch *b_chisqr = tree->Branch("chi2min", &CHISQRMIN, "chi2min/F");

	TH1 *chi2 = new TH1F("chi2", "", 100, -10.0, 30.0);

	TH1 *pi01 = new TH1F("pi01", "", 100, 600.0, 1000.0);
	TH1 *pi02 = new TH1F("pi02", "", 100, 0.0, 200.0);

	Bool_t data_flag;
	Bool_t cond_time_clus[2];

	// Progress bar
	boost::progress_display show_progress(nentries);
	// --------------------------------------------------------------------

	Int_t mode = 1;

	// Initialization of Charged part of decay reconstruction class
	// Constructor is below, in the loop
	boost::optional<KLOE::ChargedVtxRec<Float_t, UChar_t>> eventAnalysis;
	// -------------------------------------------------------------

	std::vector<std::string> ConstSet = {
			"EnergyConsvLAB",
			"PxConsvLAB",
			"PyConsvLAB",
			"PzConsvLAB",
			"Photon1PathLAB",
			"Photon2PathLAB",
			"Photon3PathLAB",
			"Photon4PathLAB"};

	KLOE::KinFitter kinematicFitObj("Omega", N_free, N_const, 8, 0, 10, 0.001, logger);
	kinematicFitObj.ConstraintSet(ConstSet);

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

		if (data_flag)
		{
			// Reconstruction of the charged part of the decay - vtx closest to the IP - to be included

			// Construction of the charged rec class object
			eventAnalysis.emplace(baseKin.nv, baseKin.ntv, baseKin.ivOld, bhabha_vtx, baseKin.CurvOld, baseKin.PhivOld, baseKin.CotvOld, baseKin.xvOld, baseKin.yvOld, baseKin.zvOld, mode);

			// // Finding the vtx closest to the interaction point
			// eventAnalysis->findKClosestRec(baseKin.Kchrec, baseKin.trk[0], baseKin.trk[1], baseKin.vtaken, baseKin.errFlag);
			// // ----------------------------------------------------------------------------------------

			// Looping over all clusters in the event
			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;

							clusterEnergy = cluster[4][baseKin.ncll[ind_gam[0]] - 1] > MIN_CLU_ENE &&
															cluster[4][baseKin.ncll[ind_gam[1]] - 1] > MIN_CLU_ENE &&
															cluster[4][baseKin.ncll[ind_gam[2]] - 1] > MIN_CLU_ENE &&
															cluster[4][baseKin.ncll[ind_gam[3]] - 1] > MIN_CLU_ENE;

							for (Int_t k = 0; k < 4; k++)
							{
								cond_clus[k] =
										cluster[3][baseKin.ncll[ind_gam[k]] - 1] > 0 &&
										cluster[0][baseKin.ncll[ind_gam[k]] - 1] != 0 &&
										cluster[1][baseKin.ncll[ind_gam[k]] - 1] != 0 &&
										cluster[2][baseKin.ncll[ind_gam[k]] - 1] != 0;
							}

							Bool_t cond_tot = cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3] && clusterEnergy;

							if (cond_tot)
							{
								for (Int_t k = 0; k < 4; k++)
								{
									Param[k * 5] = cluster[0][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 1] = cluster[1][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 2] = cluster[2][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 3] = cluster[3][baseKin.ncll[ind_gam[k]] - 1];
									Param[k * 5 + 4] = cluster[4][baseKin.ncll[ind_gam[k]] - 1];

									Errors[k * 5, k * 5] = pow(clu_x_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]), 2);
									Errors[k * 5 + 1, k * 5 + 1] = pow(clu_y_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]), 2);
									Errors[k * 5 + 2, k * 5 + 2] = pow(clu_z_error(Param[k * 5], Param[k * 5 + 1], Param[k * 5 + 2], Param[k * 5 + 4]), 2); // cm
									Errors[k * 5 + 3, k * 5 + 3] = pow(clu_time_error(Param[k * 5 + 4]), 2);																								// ns
									Errors[k * 5 + 4, k * 5 + 4] = pow(clu_ene_error(Param[k * 5 + 4]), 2);																									// MeV
								}

								Param[20] = baseKin.Kchrec[6];
								Param[21] = baseKin.Kchrec[7];
								Param[22] = baseKin.Kchrec[8];

								Errors[20, 20] = pow(Double_t(properties["variables"]["Resolutions"]["vtxCharged"][0]), 2);
								Errors[21, 21] = pow(Double_t(properties["variables"]["Resolutions"]["vtxCharged"][1]), 2);
								Errors[22, 22] = pow(Double_t(properties["variables"]["Resolutions"]["vtxCharged"][2]), 2);

								Param[23] = baseKin.Kchrec[0];
								Param[24] = baseKin.Kchrec[1];
								Param[25] = baseKin.Kchrec[2];
								Param[26] = baseKin.Kchrec[3];

								Errors[23, 23] = pow(Double_t(properties["variables"]["Resolutions"]["momCharged"][0]), 2);
								Errors[24, 24] = pow(Double_t(properties["variables"]["Resolutions"]["momCharged"][1]), 2);
								Errors[25, 25] = pow(Double_t(properties["variables"]["Resolutions"]["momCharged"][2]), 2);
								Errors[26, 26] = pow(Double_t(properties["variables"]["Resolutions"]["momCharged"][3]), 2);

								Param[27] = bhabha_mom[0];
								Param[28] = bhabha_mom[1];
								Param[29] = bhabha_mom[2];
								Param[30] = bhabha_mom[3];

								Errors[27, 27] = pow(bhabha_mom_err[0], 2);
								Errors[28, 28] = pow(bhabha_mom_err[1], 2);
								Errors[29, 29] = pow(bhabha_mom_err[2], 2);
								Errors[30, 30] = pow(bhabha_mom_err[3], 2);

								kinematicFitObj.ParameterInitialization(Param, Errors);

								CHISQRTMP = kinematicFitObj.FitFunction();

								if (abs(CHISQRTMP) < abs(CHISQRMIN))
								{
									isConverged = 1;
									FUNVALMIN = FUNVALTMP;
									CHISQRMIN = CHISQRTMP;

									kinematicFitObj.GetResults(X_min, V_min, X_init_min, V_min, C_min, L_min);

									g4takenomega[0] = ind_gam[0];
									g4takenomega[1] = ind_gam[1];
									g4takenomega[2] = ind_gam[2];
									g4takenomega[3] = ind_gam[3];

									neu_vtx_min[0] = baseKin.Kchrec[6];
									neu_vtx_min[1] = baseKin.Kchrec[7];
									neu_vtx_min[2] = baseKin.Kchrec[8];
									neu_vtx_min[3] = 0.;

									for (Int_t l = 0; l < 4; l++)
									{
										neutral_mom(cluster[0][baseKin.ncll[g4takenomega[l]] - 1], cluster[1][baseKin.ncll[g4takenomega[l]] - 1], cluster[2][baseKin.ncll[g4takenomega[l]] - 1], cluster[4][baseKin.ncll[g4takenomega[l]] - 1], neu_vtx_min, gammaomega[l]);

										gammaomega[l][4] = cluster[0][baseKin.ncll[g4takenomega[l]] - 1];
										gammaomega[l][5] = cluster[1][baseKin.ncll[g4takenomega[l]] - 1];
										gammaomega[l][6] = cluster[2][baseKin.ncll[g4takenomega[l]] - 1];
										gammaomega[l][7] = cluster[3][baseKin.ncll[g4takenomega[l]] - 1];
									}

									Int_t
											g4takenPi0[2][2];

									Float_t
											PhotonMomPi0[2][2][4],
											Pi0Mom[2][4],
											OmegaMom[4],
											M_omega_tmp[2],
											M_omega_diff[2],
											Pi0NonOmega[4];

									Pi0PhotonPair(g4takenomega, gammaomega, g4takenPi0, PhotonMomPi0, Pi0Mom, false, PichFourMom, OmegaMom);

									PichFourMom[0][0] = baseKin.trk[0][0];
									PichFourMom[0][1] = baseKin.trk[0][1];
									PichFourMom[0][2] = baseKin.trk[0][2];
									PichFourMom[0][3] = baseKin.trk[0][3];

									PichFourMom[1][0] = baseKin.trk[1][0];
									PichFourMom[1][1] = baseKin.trk[1][1];
									PichFourMom[1][2] = baseKin.trk[1][2];
									PichFourMom[1][3] = baseKin.trk[1][3];

									Float_t
											lengthPhoton[4] = {0.},
											lengthPhotonAvg = 0.,
											totEnergy = 0.;

									for (Int_t k = 0; k < 4; k++)
									{
										lengthPhoton[k] = cluster[3][baseKin.ncll[g4takenomega[k]] - 1] -
																			(sqrt(pow(cluster[0][baseKin.ncll[g4takenomega[k]] - 1] - bhabha_vtx[0], 2) +
																						pow(cluster[1][baseKin.ncll[g4takenomega[k]] - 1] - bhabha_vtx[1], 2) +
																						pow(cluster[2][baseKin.ncll[g4takenomega[k]] - 1] - baseKin.Kchrec[8], 2)) /
																			 cVel);

										totEnergy += cluster[4][baseKin.ncll[g4takenomega[k]] - 1];

										lengthPhotonAvg += pow(cluster[4][baseKin.ncll[g4takenomega[k]] - 1] * lengthPhoton[k], 2);
									}

									lengthPhotonAvg = sqrt(lengthPhotonAvg) / totEnergy;

									lengthPhotonMin[0] = lengthPhoton[0];
									lengthPhotonMin[1] = lengthPhoton[1];
									lengthPhotonMin[2] = lengthPhoton[2];
									lengthPhotonMin[3] = lengthPhoton[3];

									Float_t
											photonVel[3] = {0.},
											kaonMom[3] = {0.},
											photonLength = 0.,
											IP_vtx[3] = {bhabha_vtx[0], bhabha_vtx[1], baseKin.Kchrec[8]};

									kaonMom[0] = bhabha_mom[0] - baseKin.Kchrec[0];
									kaonMom[1] = bhabha_mom[1] - baseKin.Kchrec[1];
									kaonMom[2] = bhabha_mom[2] - baseKin.Kchrec[2];
									kaonMom[3] = bhabha_mom[3] - baseKin.Kchrec[3];

									Float_t
											kaonMomTot = sqrt(pow(kaonMom[0], 2) + pow(kaonMom[1], 2) + pow(kaonMom[2], 2)),
											kaonVelTot = cVel * (kaonMomTot / kaonMom[3]);

									for (Int_t k = 0; k < 3; k++)
										neu_vtx_avg[k] = lengthPhotonAvg * kaonVelTot * (kaonMom[k] / kaonMomTot) + IP_vtx[k];

									// Calculation of pi+pi-pi0 invariant mass - best option chosen by quadrature.
									for (Int_t j = 0; j < 2; j++)
									{
										M_omega_tmp[j] = sqrt(pow(PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[j][3], 2) -
																					pow(PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[j][0], 2) -
																					pow(PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[j][1], 2) -
																					pow(PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[j][2], 2));
										M_omega_diff[j] = M_omega_tmp[j] - mOmega;
									}

									if (std::isnan(M_omega_diff[0]) || std::isnan(M_omega_diff[1]))
										Omegarec[5] = 999999.;
									else if (std::isinf(M_omega_diff[0]) || std::isinf(M_omega_diff[1]))
										Omegarec[5] = 999999.;
									else
									{
										doneOmega = 1;
										if (abs(M_omega_diff[0]) < abs(M_omega_diff[1]))
										{
											Omegarec[0] = PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[0][0];
											Omegarec[1] = PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[0][1];
											Omegarec[2] = PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[0][2];
											Omegarec[3] = PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[0][3];
											Omegarec[4] = sqrt(pow(Omegarec[0], 2) + pow(Omegarec[1], 2) + pow(Omegarec[2], 2));
											Omegarec[5] = M_omega_tmp[0];

											Omegapi0[0] = Pi0Mom[0][0];
											Omegapi0[1] = Pi0Mom[0][1];
											Omegapi0[2] = Pi0Mom[0][2];
											Omegapi0[3] = Pi0Mom[0][3];
											Omegapi0[4] = sqrt(pow(Pi0Mom[0][0], 2) + pow(Pi0Mom[0][1], 2) + pow(Pi0Mom[0][2], 2));
											Omegapi0[5] = sqrt(pow(Omegapi0[3], 2) - pow(Omegapi0[4], 2));

											Pi0[0] = Pi0Mom[1][0];
											Pi0[1] = Pi0Mom[1][1];
											Pi0[2] = Pi0Mom[1][2];
											Pi0[3] = Pi0Mom[1][3];
											Pi0[4] = sqrt(pow(Pi0Mom[1][0], 2) + pow(Pi0Mom[1][1], 2) + pow(Pi0Mom[1][2], 2));
											Pi0[5] = sqrt(pow(Pi0[3], 2) - pow(Pi0[4], 2));
										}
										else
										{
											Omegarec[0] = PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[1][0];
											Omegarec[1] = PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[1][1];
											Omegarec[2] = PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[1][2];
											Omegarec[3] = PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[1][3];
											Omegarec[4] = sqrt(pow(Omegarec[0], 2) + pow(Omegarec[1], 2) + pow(Omegarec[2], 2));
											Omegarec[5] = M_omega_tmp[1];

											Omegapi0[0] = Pi0Mom[1][0];
											Omegapi0[1] = Pi0Mom[1][1];
											Omegapi0[2] = Pi0Mom[1][2];
											Omegapi0[3] = Pi0Mom[1][3];
											Omegapi0[4] = sqrt(pow(Pi0Mom[1][0], 2) + pow(Pi0Mom[1][1], 2) + pow(Pi0Mom[1][2], 2));
											Omegapi0[5] = sqrt(pow(Omegapi0[3], 2) - pow(Omegapi0[4], 2));

											Pi0[0] = Pi0Mom[0][0];
											Pi0[1] = Pi0Mom[0][1];
											Pi0[2] = Pi0Mom[0][2];
											Pi0[3] = Pi0Mom[0][3];
											Pi0[4] = sqrt(pow(Pi0Mom[0][0], 2) + pow(Pi0Mom[0][1], 2) + pow(Pi0Mom[0][2], 2));
											Pi0[5] = sqrt(pow(Pi0[3], 2) - pow(Pi0[4], 2));
										}

										Float_t
												Kne[4] = {Pi0[0] + Omegapi0[0], Pi0[1] + Omegapi0[1], Pi0[2] + Omegapi0[2], Pi0[3] + Omegapi0[3]},
												boost_vec_Kchboost[3] = {-(bhabha_mom[0] - Kne[0]) / (bhabha_mom[3] - Kne[3]),
																								 -(bhabha_mom[1] - Kne[1]) / (bhabha_mom[3] - Kne[3]),
																								 -(bhabha_mom[2] - Kne[2]) / (bhabha_mom[3] - Kne[3])},
												PichFourMomKaonCM[2][4];

										Obj.lorentz_transf(boost_vec_Kchboost, PichFourMom[0], PichFourMomKaonCM[0]);
										Obj.lorentz_transf(boost_vec_Kchboost, PichFourMom[1], PichFourMomKaonCM[1]);

										TVector3
												pich1(PichFourMomKaonCM[0][0], PichFourMomKaonCM[0][1], PichFourMomKaonCM[0][2]),
												pich2(PichFourMomKaonCM[1][0], PichFourMomKaonCM[1][1], PichFourMomKaonCM[1][2]);

										anglePichKaonCM = pich1.Angle(pich2) * 180. / M_PI;

										// Lorentz transformation of Pi0 to Kaon CM frame
										Float_t
												boost_vec_Kne[3] = {-(Kne[0]) / (Kne[3]),
																						-(Kne[1]) / (Kne[3]),
																						-(Kne[2]) / (Kne[3])},
												Pi0KaonCM[2][4];

										Obj.lorentz_transf(boost_vec_Kne, Pi0Mom[0], Pi0KaonCM[0]);
										Obj.lorentz_transf(boost_vec_Kne, Pi0Mom[1], Pi0KaonCM[1]);

										TVector3
												pi01(Pi0KaonCM[0][0], Pi0KaonCM[0][1], Pi0KaonCM[0][2]),
												pi02(Pi0KaonCM[1][0], Pi0KaonCM[1][1], Pi0KaonCM[1][2]);

										anglePi0KaonCM = pi01.Angle(pi02) * 180. / M_PI;

										// Lorentz transformation of Pi0 to Kaon CM frame
										Float_t
												boost_vec_phi[3] = {-(bhabha_mom[0]) / (bhabha_mom[3]),
																						-(bhabha_mom[1]) / (bhabha_mom[3]),
																						-(bhabha_mom[2]) / (bhabha_mom[3])},
												Pi0NonOmegaCM[4],
												OmegaMomCM[4];

										Obj.lorentz_transf(boost_vec_phi, Pi0NonOmega, Pi0NonOmegaCM);
										Obj.lorentz_transf(boost_vec_phi, OmegaMom, OmegaMomCM);

										TVector3
												pi0CM(Pi0NonOmegaCM[0], Pi0NonOmegaCM[1], Pi0NonOmegaCM[2]),
												omegaCM(OmegaMomCM[0], OmegaMomCM[1], OmegaMomCM[2]);

										anglePi0OmegaPhiCM = pi0CM.Angle(omegaCM) * 180. / M_PI;

										// Angle between phi and omega

										TVector3
												phi(bhabha_mom[0], bhabha_mom[1], bhabha_mom[2]),
												omega(OmegaMom[0], OmegaMom[1], OmegaMom[2]);

										anglePhiOmega = phi.Angle(omega) * 180. / M_PI;

										rho_00 = 0.;
										rho_00_IP = sqrt(pow(neu_vtx_avg[0] - IP_vtx[0], 2) + pow(neu_vtx_avg[1] - IP_vtx[1], 2));
										rho_pm_IP = sqrt(pow(baseKin.Kchrec[6] - IP_vtx[0], 2) + pow(baseKin.Kchrec[7] - IP_vtx[1], 2));
										rho = sqrt(pow(rho_00_IP, 2) + pow(rho_pm_IP, 2));
									}
								}
							}
							chi2->Fill(CHISQRMIN);
							pi01->Fill(Omegarec[5]);
							pi02->Fill(Pi0[5]);
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

		// tree->Fill();

		++show_progress; // Progress of the loading bar
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

	// file->Write();
	// file->Close();
	// delete file;

	std::ofstream outfile(propName);
	outfile << properties.dump(4);
	outfile.close();

	return 0;
}
