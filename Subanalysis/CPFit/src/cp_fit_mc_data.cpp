#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRatioPlot.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TBufferJSON.h>

#include "../inc/cpfit.hpp"

int cp_fit_mc_data(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj, ConfigWatcher &cfgWatcher)
{
	// =============================================================================
	BaseKinematics
			baseKin;
	NeutRec4
			neutVars;
	Int_t
			file_num;
	TFile
			*file_corr,
			*file_mctruth,
			*file_omega,
			*file;
	TTree
			*tree,
			*tree_mctruth,
			*tree_omega;
	// =============================================================================

	properties = cfgWatcher.getConfig();

	const TString cpfit_res_dir = cpfit_dir + result_dir;

	Double_t *eff_vals;

	if (check_corr == true)
	{
		file_corr = new TFile("../Efficiency_analysis/correction_factor.root");
		TGraphAsymmErrors *eff_signal = (TGraphAsymmErrors *)file_corr->Get("correction_factor");
		eff_vals = eff_signal->GetY();

		delete eff_signal;
	}

	// ===========================================================================

	chain.SetBranchAddress("mcflag", &baseKin.mcflag);

	chain.SetBranchAddress("Dtmc", &baseKin.Dtmc);

	chain.SetBranchAddress("Chi2", &baseKin.Chi2);
	chain.SetBranchAddress("minv4gam", &baseKin.minv4gam);
	chain.SetBranchAddress("Kchboost", baseKin.Kchboost);
	chain.SetBranchAddress("Kchrec", baseKin.Kchrec);
	chain.SetBranchAddress("Qmiss", &baseKin.Qmiss);
	chain.SetBranchAddress("ip", interfcommon_.ip);

	Float_t
			PichFourMom[2][4];

	chain.SetBranchAddress("trk1", PichFourMom[0]);
	chain.SetBranchAddress("trk2", PichFourMom[1]);

	chain.SetBranchAddress("Xcl", baseKin.cluster[0]);
	chain.SetBranchAddress("Ycl", baseKin.cluster[1]);
	chain.SetBranchAddress("Zcl", baseKin.cluster[2]);
	chain.SetBranchAddress("Tcl", baseKin.cluster[3]);
	chain.SetBranchAddress("Enecl", baseKin.cluster[4]);

	chain.SetBranchAddress("ip", baseKin.ip);

	chain.SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
	chain.SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
	chain.SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
	chain.SetBranchAddress("Broots", &baseKin.phi_mom[3]);

	chain.SetBranchAddress("Bx", &baseKin.bhabha_vtx[0]);
	chain.SetBranchAddress("By", &baseKin.bhabha_vtx[1]);
	chain.SetBranchAddress("Bz", &baseKin.bhabha_vtx[2]);

	// ===========================================================================

	TString
			filename_trilateration = std::string(properties["variables"]["tree"]["filename"]["trilateration"]),
			treename_trilateration = std::string(properties["variables"]["tree"]["treename"]["trilateration"]),

			filename_trilateration_kin_fit = std::string(properties["variables"]["tree"]["filename"]["trilaterationKinFit"]),
			treename_trilateration_kin_fit = std::string(properties["variables"]["tree"]["treename"]["trilaterationKinFit"]),

			filename_triangle = std::string(properties["variables"]["tree"]["filename"]["trianglefinal"]),
			treename_triangle = std::string(properties["variables"]["tree"]["treename"]["trianglefinal"]),

			filename_omega = std::string(properties["variables"]["tree"]["filename"]["omegarec"]),
			treename_omega = std::string(properties["variables"]["tree"]["treename"]["omegarec"]),

			filename_mctruth = std::string(properties["variables"]["tree"]["filename"]["mctruth"]),
			treename_mctruth = std::string(properties["variables"]["tree"]["treename"]["mctruth"]);

	std::vector<TString>
			file_name,
			tree_name;

	file_name.push_back(filename_trilateration);
	file_name.push_back(filename_trilateration_kin_fit);
	file_name.push_back(filename_triangle);

	tree_name.push_back(treename_trilateration);
	tree_name.push_back(treename_trilateration_kin_fit);
	tree_name.push_back(treename_triangle);

	Controls::Menu menu(3, file_name);

	menu.ShowOpt();
	menu.EndMenu();

	try
	{
		std::cin >> file_num;

		if (!std::cin)
			throw(ErrorHandling::ErrorCodes::DATA_TYPE);
		else if (file_num < 1 || file_num > file_name.size())
			throw(ErrorHandling::ErrorCodes::MENU_RANGE);

		file = new TFile(file_name[file_num]);

		if (file->IsZombie())
			throw(ErrorHandling::ErrorCodes::FILE_NOT_EXIST);

		tree = (TTree *)file->Get(tree_name[file_num]);

		if (tree->IsZombie())
			throw(ErrorHandling::ErrorCodes::TREE_NOT_EXIST);

		file_omega = new TFile(filename_omega);

		if (file_omega->IsZombie())
			throw(ErrorHandling::ErrorCodes::FILE_NOT_EXIST);

		tree_omega = (TTree *)file_omega->Get(treename_omega);

		if (tree_omega->IsZombie())
			throw(ErrorHandling::ErrorCodes::TREE_NOT_EXIST);

		file_mctruth = new TFile(filename_mctruth);

		if (file_mctruth->IsZombie())
			throw(ErrorHandling::ErrorCodes::FILE_NOT_EXIST);

		tree_mctruth = (TTree *)file_mctruth->Get(treename_mctruth);

		if (tree_mctruth->IsZombie())
			throw(ErrorHandling::ErrorCodes::TREE_NOT_EXIST);
	}
	catch (ErrorHandling::ErrorCodes err)
	{
		logger.getErrLog(err);
		return int(err);
	}

	if (file_num == 0)
	{
	}
	else if (file_num == 1)
	{
		tree->SetBranchAddress("done4_kinfit", &neutVars.done);
		tree->SetBranchAddress("fourKnetri_kinfit", neutVars.Knerec);
		tree->SetBranchAddress("g4takentri_kinfit", neutVars.gtaken);
		tree->SetBranchAddress("chi2min", &neutVars.chi2min);
	}
	else if (file_num == 2)
	{
		tree->SetBranchAddress("done_triangle", &neutVars.done);
		tree->SetBranchAddress("fourKnetriangle", neutVars.Knerec);
		tree->SetBranchAddress("g4taken_triangle", neutVars.gtaken);
		tree->SetBranchAddress("chi2min", &neutVars.chi2min);
	}

	Int_t
			doneOmega,
			g4takenomega[4];

	Float_t
			gammaomega[4][8],
			Omegapi0[6],
			Pi0[6],
			PichFourMom1[2][4],
			Omegarec[7],
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
			neu_vtx_avg[3],
			chi2min;

	tree_omega->SetBranchAddress("gamma1omega", gammaomega[0]);
	tree_omega->SetBranchAddress("gamma2omega", gammaomega[1]);
	tree_omega->SetBranchAddress("gamma3omega", gammaomega[2]);
	tree_omega->SetBranchAddress("gamma4omega", gammaomega[3]);

	tree_omega->SetBranchAddress("pich1", PichFourMom1[0]);
	tree_omega->SetBranchAddress("pich2", PichFourMom1[1]);

	tree_omega->SetBranchAddress("omegapi0", Omegapi0);
	tree_omega->SetBranchAddress("pi0", Pi0);

	tree_omega->SetBranchAddress("omega", Omegarec);

	tree_omega->SetBranchAddress("doneomega", &doneOmega);

	tree_omega->SetBranchAddress("g4takenomega", g4takenomega);

	tree_omega->SetBranchAddress("lengthphoton", lengthPhotonMin);

	tree_omega->SetBranchAddress("lengthKch", &lengthKch);

	tree_omega->SetBranchAddress("NeuVtxAvg", neu_vtx_avg);

	tree_omega->SetBranchAddress("rho_00", &rho_00);
	tree_omega->SetBranchAddress("rho_00_IP", &rho_00_IP);
	tree_omega->SetBranchAddress("rho_pm_IP", &rho_pm_IP);
	tree_omega->SetBranchAddress("rho", &rho);
	tree_omega->SetBranchAddress("anglePi0KaonCM", &anglePi0KaonCM);
	tree_omega->SetBranchAddress("anglePichKaonCM", &anglePichKaonCM);
	tree_omega->SetBranchAddress("anglePi0OmegaPhiCM", &anglePi0OmegaPhiCM);

	tree_mctruth->SetBranchAddress("mctruth", &baseKin.mctruth_int);

	chain.AddFriend(tree_mctruth);
	chain.AddFriend(tree_omega);
	chain.AddFriend(tree);

	UInt_t nentries = chain.GetEntries();

	const Double_t
			x_min = properties["variables"]["CPFit"]["histoResults"]["rangeX"][0],
			x_max = properties["variables"]["CPFit"]["histoResults"]["rangeX"][1],
			res_deltaT = properties["variables"]["Resolutions"]["deltaT"];
	const UInt_t
			nbins = 1 + ((x_max - x_min) / res_deltaT);

	Double_t split[3] = {-30.0, 0.0, 30.0};

	KLOE::interference event(mode, check_corr, nbins, x_min, x_max, split);

	Double_t velocity_kch, velocity_kne, tch_LAB, tch_CM, tne_LAB, tne_CM, trcv_sum, TRCV[4];
	Float_t Kch_LAB[4], Kne_LAB[4], Kch_CM[4], Kne_CM[4], Kch_CMCM[4], Kne_CMCM[4], Kne_boost[3], Kch_boost[3], Phi_boost[3], Kchmom_LAB[4], Knemom_LAB[4], Kchmom_CM[4], Knemom_CM[4];

	TH1 *sig_total = new TH1D("sig_tot", ";#Delta t;Counts", nbins, x_min, x_max);
	TH1 *sig_pass = new TH1D("sig_pass", ";#Delta t;Counts", nbins, x_min, x_max);

	std::vector<Float_t> no_cuts_sig[2];

	// Initialization of cut formulas
	Int_t cuts_number = 6;

	std::vector<TFormula> formula_vector(cuts_number);
	std::vector<Double_t *> mean_sigma_vector(cuts_number);
	std::vector<Bool_t> values_vector(cuts_number);
	std::vector<Double_t> x_vector(cuts_number);
	std::vector<TString> formula_names =
			{
					"DriftChamberRegenerationLeft",
					"BeamPipeRegenerationLeft",
					"DriftChamberRegenerationRight",
					"BeamPipeRegenerationRight",
					"RhoCutOmega",
					"ChiSqrCut"};

	Bool_t scanFlag = (Bool_t)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["flag"];
	Int_t numberOfPoints;
	Double_t
			cutLimits[2] = {0.0},
			cutStep = 0.0;

	if (scanFlag)
	{
		numberOfPoints = (Int_t)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["numberOfPoints"];
		cutLimits[0] = (Double_t)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["cutLimits"][0];
		cutLimits[1] = (Double_t)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["cutLimits"][1];
	}
	else
	{
		numberOfPoints = 1;
		cutLimits[0] = 0.0;
		cutLimits[1] = 0.0;
	}

	cutStep = abs(cutLimits[1] - cutLimits[0]) / (Double_t)numberOfPoints;

	// Vectors needed for the initialization of graphs
	TString graphMode = "FitResultErr";
	std::vector<Double_t>
			cutLimit(numberOfPoints);
	std::vector<std::vector<Double_t>>
			errValue(2),
			realValue(2),
			imaginaryValue(2);

	TString xTitle = (std::string)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["cutTitle"],
					yTitle = "|#sigma(Re(#varepsilon'/#varepsilon))/Re(#varepsilon'/#varepsilon)|",
					yRightTitle = "|#sigma(Im(#varepsilon'/#varepsilon))/Im(#varepsilon'/#varepsilon)|",
					yTitleReal = "Re(#varepsilon'/#varepsilon)",
					yRightTitleReal = "#sigma(Re(#varepsilon'/#varepsilon))",
					yTitleImaginary = "Im(#varepsilon'/#varepsilon)",
					yRightTitleImaginary = "#sigma(Im(#varepsilon'/#varepsilon))",
					cutName = (std::string)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["cutName"];

	// Values of cuts to be used in the analysis
	Double_t
			rho_cut = 0.0;

	std::vector<Double_t> par, parErr;

	for (Int_t scanIter = 1; scanIter <= numberOfPoints; scanIter++)
	{
		Bool_t
				simona_cuts = properties["variables"]["CPFit"]["cuts"]["simonaCuts"]["flag"];

		Double_t
				sigmas = 1.5,
				sigma = (Double_t)properties["variables"]["RegenRejection"]["sigma"],
				ChHigher[2] = {
						properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["mean"][1],
						sigma * (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["width"][1]},
				ChLower[2] = {properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["mean"][0], sigma * (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["width"][0]}, TrHigher[2] = {properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["mean"][1], sigma * (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["width"][1]}, TrLower[2] = {properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["mean"][0], sigma * (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["width"][0]}, OmegaCut1[2], OmegaCut2[2];

		if (properties["variables"]["CPFit"]["cuts"]["omegaRho"].is_null())
			rho_cut = 0.;
		else if (scanFlag && ToLower(cutName) == "rho")
		{
			rho_cut = scanIter * cutStep;
			cutLimit[scanIter - 1] = rho_cut;

			std::cout << rho_cut << std::endl;
		}
		else
			rho_cut = properties["variables"]["CPFit"]["cuts"]["omegaRho"];
		// ---------------------------------------------------------------
		if (scanFlag && ToLower(cutName) == "regen")
		{
			TrHigher[1] = scanIter * cutStep;
			ChHigher[1] = scanIter * cutStep;

			cutLimit[scanIter - 1] = scanIter * cutStep;
		}
		else
		{
			TrHigher[1] = 1.5;
			ChHigher[1] = 1.5;
		}
		// ---------------------------------------------------------------

		Double_t a, b, lineWidth;

		if (!simona_cuts)
		{
			OmegaCut1[0] = 0.0;
			OmegaCut1[1] = rho_cut;

			OmegaCut2[0] = 0.0;
			OmegaCut2[1] = 40.0;
		}
		else
		{
			Double_t
					meanInvMass = properties["variables"]["OmegaRec"]["invMass"]["mean"]["value"],
					// meanInvMassErr = properties["variables"]["OmegaRec"]["invMass"]["mean"]["error"],
					stdInvMass = properties["variables"]["OmegaRec"]["invMass"]["stdDev"]["value"],
					// stdInvMassErr = properties["variables"]["OmegaRec"]["invMass"]["stdDev"]["error"],
					InvMass[2] = {meanInvMass, 3 * stdInvMass},
					meanKinEne = properties["variables"]["OmegaRec"]["kinEne"]["mean"]["value"],
					// meanKinEneErr = properties["variables"]["OmegaRec"]["kinEne"]["mean"]["error"],
					stdKinEne = properties["variables"]["OmegaRec"]["kinEne"]["stdDev"]["value"],
					// stdKinEneErr = properties["variables"]["OmegaRec"]["kinEne"]["stdDev"]["error"],
					KinEne[2] = {meanKinEne, 3 * stdKinEne};

			a = properties["variables"]["OmegaRec"]["combined"]["line"]["slope"],
			b = properties["variables"]["OmegaRec"]["combined"]["line"]["inter"],
			lineWidth = properties["variables"]["OmegaRec"]["combined"]["stdDev"]["value"];

			OmegaCut1[0] = meanInvMass;
			OmegaCut1[1] = sigmas * stdInvMass;

			OmegaCut2[0] = meanKinEne;
			OmegaCut2[1] = sigmas * stdKinEne;
		}

		mean_sigma_vector[0] = ChHigher;
		mean_sigma_vector[1] = ChLower;
		mean_sigma_vector[2] = TrHigher;
		mean_sigma_vector[3] = TrLower;
		mean_sigma_vector[4] = OmegaCut1;
		mean_sigma_vector[5] = OmegaCut2;

		Obj.ConditionInitializer(formula_vector, formula_names);

		Obj.SetAllConditionParameters(formula_vector, mean_sigma_vector);
		// ---------------------------------------------------------

		Bool_t
				cut_regen_neg,
				cut_regen_pos,
				cut_line_down,
				cut_line_up,
				cut_KinEnePi0,
				cut_MinvOmega,
				cut_total;

		const Int_t numberOfMomenta = 2;

		TVectorT<Double_t>
				momVecMC(numberOfMomenta * 3),
				momVecSmeared(numberOfMomenta * 3);

		std::vector<double> elems = properties["momSmearing"]["covarianceMatrix"]["fElements"].get<std::vector<Double_t>>();

		Int_t nRows = properties["momSmearing"]["covarianceMatrix"]["fNrows"],
					nCols = properties["momSmearing"]["covarianceMatrix"]["fNcols"];

		TMatrixT<Double_t>
				covMatrix(nRows, nCols, elems.data());

		KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, covMatrix);
		KLOE::ChargedVtxRec<Float_t> BoostMethodObj;
		

		for (UInt_t i = 0; i < nentries; i++)
		{
			chain.GetEntry(i);

			Float_t KchrecSmeared[4], KchboostSmeared[4], energyPion[2];

			if (neutVars.done == 1 && doneOmega == 1)
			{
				if (baseKin.mcflag == 1)
				{
					momVecMC[0] = PichFourMom[0][0];
					momVecMC[1] = PichFourMom[0][1];
					momVecMC[2] = PichFourMom[0][2];
					momVecMC[3] = PichFourMom[1][0];
					momVecMC[4] = PichFourMom[1][1];
					momVecMC[5] = PichFourMom[1][2];

					CovMatrixCalcObj.SetMCVector(momVecMC);
					CovMatrixCalcObj.SmearMomentum();
					CovMatrixCalcObj.GetSmearedMomentum(momVecSmeared);

					KchrecSmeared[0] = momVecSmeared[0] + momVecSmeared[3];
					KchrecSmeared[1] = momVecSmeared[1] + momVecSmeared[4];
					KchrecSmeared[2] = momVecSmeared[2] + momVecSmeared[5];

					energyPion[0] = sqrt(pow(momVecSmeared[0], 2) + pow(momVecSmeared[1], 2) + pow(momVecSmeared[2], 2) + pow(mPiCh, 2));
					energyPion[1] = sqrt(pow(momVecSmeared[3], 2) + pow(momVecSmeared[4], 2) + pow(momVecSmeared[5], 2) + pow(mPiCh, 2));

					KchrecSmeared[3] = energyPion[0] + energyPion[1];

					BoostMethodObj.KaonMomFromBoost(KchrecSmeared, baseKin.phi_mom, KchboostSmeared);

					Float_t X_line[3] = {baseKin.Kchboost[6], baseKin.Kchboost[7], baseKin.Kchboost[8]},
									p[3] = {KchboostSmeared[0], KchboostSmeared[1], KchboostSmeared[2]},
									xB[3] = {baseKin.bhabha_vtx[0], baseKin.bhabha_vtx[1], baseKin.bhabha_vtx[2]},
									plane_perp[3] = {0., baseKin.phi_mom[1], 0.};

					BoostMethodObj.IPBoostCorr(X_line, p, xB, plane_perp, baseKin.ip);

					baseKin.ip[0] = baseKin.bhabha_vtx[0];
					baseKin.ip[1] = baseKin.bhabha_vtx[1];
					if (abs(baseKin.ip[2] - baseKin.bhabha_vtx[2]) > 2.0)
						baseKin.ip[2] = baseKin.bhabha_vtx[2];
				}
				else
				{
					KchboostSmeared[0] = baseKin.Kchboost[0];
					KchboostSmeared[1] = baseKin.Kchboost[1];
					KchboostSmeared[2] = baseKin.Kchboost[2];
					KchboostSmeared[3] = baseKin.Kchboost[3];
				}

				velocity_kch = cVel * sqrt(pow(KchboostSmeared[0], 2) + pow(KchboostSmeared[1], 2) + pow(KchboostSmeared[2], 2)) / KchboostSmeared[3];

				velocity_kne = cVel * sqrt(pow(neutVars.Knerec[0], 2) + pow(neutVars.Knerec[1], 2) + pow(neutVars.Knerec[2], 2)) / neutVars.Knerec[3];

				tch_LAB = sqrt(pow(baseKin.Kchboost[6] - baseKin.ip[0], 2) + pow(baseKin.Kchboost[7] - baseKin.ip[1], 2) + pow(baseKin.Kchboost[8] - baseKin.ip[2], 2)) / velocity_kch;
				tne_LAB = sqrt(pow(neutVars.Knerec[6] - baseKin.ip[0], 2) + pow(neutVars.Knerec[7] - baseKin.ip[1], 2) + pow(neutVars.Knerec[8] - baseKin.ip[2], 2)) / velocity_kne;

				Kch_LAB[0] = baseKin.Kchboost[6] - baseKin.ip[0];
				Kch_LAB[1] = baseKin.Kchboost[7] - baseKin.ip[1];
				Kch_LAB[2] = baseKin.Kchboost[8] - baseKin.ip[2];
				Kch_LAB[3] = tch_LAB * cVel;

				Kchmom_LAB[0] = KchboostSmeared[0];
				Kchmom_LAB[1] = KchboostSmeared[1];
				Kchmom_LAB[2] = KchboostSmeared[2];
				Kchmom_LAB[3] = KchboostSmeared[3];

				Kne_LAB[0] = neutVars.Knerec[6] - baseKin.ip[0];
				Kne_LAB[1] = neutVars.Knerec[7] - baseKin.ip[1];
				Kne_LAB[2] = neutVars.Knerec[8] - baseKin.ip[2];
				Kne_LAB[3] = tne_LAB * cVel;

				Knemom_LAB[0] = neutVars.Knerec[0];
				Knemom_LAB[1] = neutVars.Knerec[1];
				Knemom_LAB[2] = neutVars.Knerec[2];
				Knemom_LAB[3] = neutVars.Knerec[3];

				Phi_boost[0] = -baseKin.phi_mom[0] / baseKin.phi_mom[3];
				Phi_boost[1] = -baseKin.phi_mom[1] / baseKin.phi_mom[3];
				Phi_boost[2] = -baseKin.phi_mom[2] / baseKin.phi_mom[3];

				Obj.lorentz_transf(Phi_boost, Kch_LAB, Kch_CM);
				Obj.lorentz_transf(Phi_boost, Kne_LAB, Kne_CM);
				Obj.lorentz_transf(Phi_boost, Kchmom_LAB, Kchmom_CM);
				Obj.lorentz_transf(Phi_boost, Knemom_LAB, Knemom_CM);

				Kch_boost[0] = -Kchmom_CM[0] / Kchmom_CM[3];
				Kch_boost[1] = -Kchmom_CM[1] / Kchmom_CM[3];
				Kch_boost[2] = -Kchmom_CM[2] / Kchmom_CM[3];

				Kne_boost[0] = Kchmom_CM[0] / Kchmom_CM[3];
				Kne_boost[1] = Kchmom_CM[1] / Kchmom_CM[3];
				Kne_boost[2] = Kchmom_CM[2] / Kchmom_CM[3];

				Obj.lorentz_transf(Kch_boost, Kch_CM, Kch_CMCM);
				Obj.lorentz_transf(Kch_boost, Kne_CM, Kne_CMCM);

				baseKin.Dtboostlor = (Kch_CMCM[3] - Kne_CMCM[3]) / (cVel * tau_S_nonCPT);

				for (Int_t i = 0; i < 4; i++)
				{
					TRCV[i] = baseKin.cluster[3][neutVars.gtaken[i]] - (sqrt(pow(baseKin.cluster[0][neutVars.gtaken[i]] - neutVars.Knerec[6], 2) + pow(baseKin.cluster[1][neutVars.gtaken[i]] - neutVars.Knerec[7], 2) + pow(baseKin.cluster[2][neutVars.gtaken[i]] - neutVars.Knerec[8], 2)) / cVel) - tne_LAB;
				}

				trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

				Double_t radius[2] = {0., 0.},
								 radius_ch[2] = {0., 0.};

				Double_t sphere_bound = 10, bp_bound = 4.4;

				Bool_t cuts[4] = {false, false, false, false};

				Double_t normFactor = sigmas * lineWidth / sqrt(a * a + 1); // Norm of normal vector
				cut_line_up = (a * Omegarec[6] + b + (1 - pow(a, 2)) * normFactor) > Omegarec[5];
				cut_line_down = (a * Omegarec[6] + b - (1 - pow(a, 2)) * normFactor) < Omegarec[5];

				for (Int_t i = 0; i < 3; i++)
				{
					radius[0] += pow(neutVars.Knerec[6 + i] - interfcommon_.ip[i], 2);
					radius_ch[0] += pow(baseKin.Kchboost[6 + i] - interfcommon_.ip[i], 2);

					if (i < 2)
					{
						radius[1] += pow(neutVars.Knerec[6 + i] - interfcommon_.ip[i], 2);
						radius_ch[1] += pow(baseKin.Kchboost[6 + i] - interfcommon_.ip[i], 2);
					}
				}

				radius[0] = sqrt(radius[0]);
				radius[1] = sqrt(radius[1]);

				radius_ch[0] = sqrt(radius_ch[0]);
				radius_ch[1] = sqrt(radius_ch[1]);

				Bool_t
						cond[4],
						cond_tot;

				std::vector<Double_t>
						stdDevOmegaVtx(6);

				// Limitation for the Omega fiducial volume - based on Rayleigh distribution
				std::string decayType[2] = {"neutral", "charged"};
				for (Int_t i = 0; i < 2; i++)
					for (Int_t j = 0; j < 3; j++)
					{
						stdDevOmegaVtx[i * 3 + j] = properties["variables"]["OmegaRec"]["fiducialVolume"][decayType[i]]["stdDev"][j];
					}

				Double_t
						stdDevRho00IP = sqrt((4 - M_PI) / 2) * (stdDevOmegaVtx[0] + stdDevOmegaVtx[1]) / 2.,
						stdDevRhopmIP = sqrt((4 - M_PI) / 2) * (stdDevOmegaVtx[3] + stdDevOmegaVtx[4]) / 2.,
						meanRho00IP = sqrt(M_PI / 2.) * (stdDevOmegaVtx[0] + stdDevOmegaVtx[1]) / 2.,
						meanRhopmIP = sqrt(M_PI / 2.) * (stdDevOmegaVtx[3] + stdDevOmegaVtx[4]) / 2.;

				cond[0] = rho_00_IP - meanRho00IP < stdDevRho00IP;
				cond[1] = rho_pm_IP - meanRhopmIP < stdDevRhopmIP;
				cond[2] = abs(neu_vtx_avg[2] - baseKin.Kchboost[8]) < stdDevOmegaVtx[2];
				cond[3] = abs(baseKin.Kchboost[8] - baseKin.bhabha_vtx[2]) < stdDevOmegaVtx[5];

				cond_tot = cond[0] && cond[1] && cond[2] && cond[3];
				// ----------------------------------------------------------------------------------------------------------------

				x_vector[0] = radius_ch[0];
				x_vector[1] = radius_ch[1];
				x_vector[2] = radius[0];
				x_vector[3] = radius[1];

				// Check of Simona's cuts
				if (simona_cuts)
				{
					x_vector[4] = Omegarec[5];
					x_vector[5] = Omegarec[6];
					// ---------------------------------------------------------------
				}
				else
				{
					x_vector[4] = rho;
					x_vector[5] = 0.0;
				}

				// cuts[0] = abs(radius[0] - meanRadiusTriangleHigher) > errorRadiusTriangleHigher * sigma; // && radius[1] > 8;
				// cuts[1] = 1;																																						 // abs(radius[1] - meanRadiusTriangleLower) > errorRadiusTriangleLower * sigma;     // && radius[1] <= 8;
				// cuts[2] = abs(radius_ch[0] - meanRadiusChHigher) > errorRadiusChHigher * sigma;					 // && radius_ch[1] > 8;
				// cuts[3] = 1;																																						 // abs(radius_ch[1] - meanRadiusChLower) > errorRadiusChLower * sigma;  // && radius_ch[1] <= 8;

				cuts[0] = Obj.GetSingleConditionValue(formula_vector[0], x_vector[0]);
				cuts[2] = Obj.GetSingleConditionValue(formula_vector[2], x_vector[2]);

				cuts[1] = 1;
				cuts[3] = 1;

				if (cond_tot)
				{
					cuts[4] = Obj.GetSingleConditionValue(formula_vector[4], x_vector[4]);

					cuts[5] = Obj.GetSingleConditionValue(formula_vector[5], x_vector[5]);
				}
				else
				{
					cuts[4] = 1;
					cuts[5] = 1;
				}

				if (baseKin.mcflag == 1)
				{
					if (baseKin.mctruth_int == 1 || baseKin.mctruth_int == 2)
					{
						no_cuts_sig[0].push_back(baseKin.Dtmc);
						no_cuts_sig[1].push_back(baseKin.Dtboostlor);
					}

					if (cuts[0] && cuts[1] && cuts[2] && cuts[3] && (cuts[4] || cuts[5] || cut_line_up || cut_line_down))
					{
						if (baseKin.mctruth_int == 1)
						{
							event.time_diff_gen.push_back(baseKin.Dtmc);
							event.time_diff[0].push_back(baseKin.Dtboostlor);
						}

						if (baseKin.mctruth_int == 3)
						{
							event.time_diff[1].push_back(baseKin.Dtboostlor);
						}

						if (baseKin.mctruth_int == 4)
						{
							event.time_diff[2].push_back(baseKin.Dtboostlor);
						}

						if (baseKin.mctruth_int == 5)
						{
							event.time_diff[3].push_back(baseKin.Dtboostlor);
						}

						if (baseKin.mctruth_int == 6)
						{
							event.time_diff[4].push_back(baseKin.Dtboostlor);
						}

						if (baseKin.mctruth_int == 7)
						{
							event.time_diff[5].push_back(baseKin.Dtboostlor);
						}
					}
				}

				if (baseKin.mcflag == 0 && cuts[0] && cuts[1] && cuts[2] && cuts[3] && (cuts[4] || cuts[5] || cut_line_up || cut_line_down))
				{
					event.time_diff_data.push_back(baseKin.Dtboostlor);
				}
			}
		}

		ROOT::Math::Minimizer *minimum =
				ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

		// set tolerance , etc...
		minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
		minimum->SetTolerance(0.1);
		minimum->SetPrintLevel(1);
		minimum->SetStrategy(2);

		const UInt_t num_of_vars = 11;

		ROOT::Math::Functor minimized_function(&event, &KLOE::interference::interf_chi2, num_of_vars);

		minimum->SetFunction(minimized_function);

		const Double_t init_vars[num_of_vars] = {
				Re,
				Im_nonCPT,
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Signal"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["FarLeft"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["CloseLeft"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["CloseRight"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["FarRight"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Omegapi0"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Threepi0"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Semileptonic"],
				properties["variables"]["CPFit"]["initParams"]["Norm"]["Other"]},
									 step[num_of_vars] = {properties["variables"]["CPFit"]["step"]["Re"], properties["variables"]["CPFit"]["step"]["Im"], properties["variables"]["CPFit"]["step"]["Norm"]["Signal"], properties["variables"]["CPFit"]["step"]["Norm"]["Regeneration"]["FarLeft"], properties["variables"]["CPFit"]["step"]["Norm"]["Regeneration"]["CloseLeft"], properties["variables"]["CPFit"]["step"]["Norm"]["Regeneration"]["CloseRight"], properties["variables"]["CPFit"]["step"]["Norm"]["Regeneration"]["FarRight"], properties["variables"]["CPFit"]["step"]["Norm"]["Omegapi0"], properties["variables"]["CPFit"]["step"]["Norm"]["Threepi0"], properties["variables"]["CPFit"]["step"]["Norm"]["Semileptonic"], properties["variables"]["CPFit"]["step"]["Norm"]["Other"]};

		Double_t
				limit_upper = properties["variables"]["CPFit"]["limitPer"]["upper"],
				limit_lower = properties["variables"]["CPFit"]["limitPer"]["lower"];

		minimum->SetVariable(0, "Real part", init_vars[0], step[0]);
		minimum->SetVariable(1, "Imaginary part", init_vars[1], step[1]);
		minimum->SetLimitedVariable(2, "Norm signal", init_vars[2], step[2], init_vars[2] - limit_lower * init_vars[2], init_vars[2] + limit_upper * init_vars[2]);

		if (x_min > -30.0 && x_max < 30.0)
		{
			minimum->SetFixedVariable(3, "Norm left DC wall", init_vars[3]);
			minimum->SetLowerLimitedVariable(4, "Norm left beam pipe", init_vars[4], step[4], 0.0);
			minimum->SetLowerLimitedVariable(5, "Norm right beam pipe", init_vars[5], step[5], 0.0);
			minimum->SetFixedVariable(6, "Norm right DC wall", init_vars[6]);
		}
		else
		{
			minimum->SetLowerLimitedVariable(3, "Norm left DC wall", init_vars[3], step[3], 0.0);
			minimum->SetLowerLimitedVariable(4, "Norm left beam pipe", init_vars[4], step[4], 0.0);
			minimum->SetLowerLimitedVariable(5, "Norm right beam pipe", init_vars[5], step[5], 0.0);
			minimum->SetLowerLimitedVariable(6, "Norm right DC wall", init_vars[6], step[6], 0.0);
		}

		minimum->SetLimitedVariable(7, "Norm omega", init_vars[7], step[7], init_vars[7] - limit_lower * init_vars[7], init_vars[7] + limit_upper * init_vars[7]);
		minimum->SetLimitedVariable(8, "Norm three", init_vars[8], step[8], init_vars[8] - limit_lower * init_vars[8], init_vars[8] + limit_upper * init_vars[8]);
		minimum->SetLimitedVariable(9, "Norm semi", init_vars[9], step[9], init_vars[9] - limit_lower * init_vars[9], init_vars[9] + limit_upper * init_vars[9]);
		minimum->SetLimitedVariable(10, "Norm other bcg", init_vars[10], step[10], init_vars[10] - limit_lower * init_vars[10], init_vars[10] + limit_upper * init_vars[10]);

		minimum->Minimize();

		if (scanFlag)
		{
			errValue[0].push_back(abs(minimum->Errors()[0] / minimum->X()[0]));
			errValue[1].push_back(abs(minimum->Errors()[1] / minimum->X()[1]));

			realValue[0].push_back(abs(minimum->X()[0]));
			realValue[1].push_back(abs(minimum->Errors()[0]));

			imaginaryValue[0].push_back(abs(minimum->X()[1]));
			imaginaryValue[1].push_back(abs(minimum->Errors()[1]));

			event.time_diff[0].clear();
			event.time_diff[1].clear();
			event.time_diff[2].clear();
			event.time_diff[3].clear();
			event.time_diff[4].clear();
			event.time_diff[5].clear();
			event.time_diff_data.clear();
			event.time_diff_gen.clear();
			no_cuts_sig[0].clear();
			no_cuts_sig[1].clear();

			event.time_diff[0].shrink_to_fit();
			event.time_diff[1].shrink_to_fit();
			event.time_diff[2].shrink_to_fit();
			event.time_diff[3].shrink_to_fit();
			event.time_diff[4].shrink_to_fit();
			event.time_diff[5].shrink_to_fit();
			event.time_diff_data.shrink_to_fit();
			event.time_diff_gen.shrink_to_fit();
			no_cuts_sig[0].shrink_to_fit();
			no_cuts_sig[1].shrink_to_fit();
		}
		else
		{
			par.push_back(minimum->X()[0]);
			par.push_back(minimum->X()[1]);
			par.push_back(minimum->X()[2]);
			par.push_back(minimum->X()[3]);
			par.push_back(minimum->X()[4]);
			par.push_back(minimum->X()[5]);
			par.push_back(minimum->X()[6]);
			par.push_back(minimum->X()[7]);
			par.push_back(minimum->X()[8]);
			par.push_back(minimum->X()[9]);
			par.push_back(minimum->X()[10]);

			parErr.push_back(minimum->Errors()[0]);
			parErr.push_back(minimum->Errors()[1]);
			parErr.push_back(minimum->Errors()[2]);
			parErr.push_back(minimum->Errors()[3]);
			parErr.push_back(minimum->Errors()[4]);
			parErr.push_back(minimum->Errors()[5]);
			parErr.push_back(minimum->Errors()[6]);
			parErr.push_back(minimum->Errors()[7]);
			parErr.push_back(minimum->Errors()[8]);
			parErr.push_back(minimum->Errors()[9]);
			parErr.push_back(minimum->Errors()[10]);
		}
	}

	if (!scanFlag)
	{
		Double_t sum_of_events = 0., fractions[6] = {0.};

		for (Int_t i = 0; i < channNum; i++)
			sum_of_events += event.time_diff[i].size();

		for (Int_t i = 0; i < channNum; i++)
			fractions[i] = 100 * event.time_diff[i].size() / sum_of_events;

		std::ofstream myfile_num;
		myfile_num.open(cpfit_res_dir + "num_of_events.csv");
		myfile_num << "Channel,Number of events,Fraction\n";
		myfile_num << "Signal," << event.time_diff[0].size() << "," << fractions[0] << "%,\n";
		myfile_num << "Regeneration," << event.time_diff[1].size() << "," << fractions[1] << "%,\n";
		myfile_num << "Omega," << event.time_diff[2].size() << "," << fractions[2] << "%,\n";
		myfile_num << "Three," << event.time_diff[3].size() << "," << fractions[3] << "%,\n";
		myfile_num << "Semi," << event.time_diff[4].size() << "," << fractions[4] << "%,\n";
		myfile_num << "Other bcg," << event.time_diff[5].size() << "," << fractions[5] << "%,\n";
		myfile_num.close();

		for (UInt_t i = 0; i < no_cuts_sig[1].size(); i++)
		{
			sig_total->Fill(no_cuts_sig[1][i]);
		}

		for (UInt_t i = 0; i < channNum; i++)
		{
			for (UInt_t j = 0; j < event.time_diff[i].size(); j++)
			{
				if (i == 0)
				{
					sig_pass->Fill(event.time_diff[i][j]);

					event.frac[i]->Fill(event.time_diff[i][j], event.interf_function(event.time_diff_gen[j], 0, par.data()));
				}
				else if (i == 1)
				{
					if (event.time_diff[i][j] < event.left_x_split)
						event.frac[i]->Fill(event.time_diff[i][j], par[3]);
					else if (event.time_diff[i][j] > event.left_x_split && event.time_diff[i][j] < event.center_x_split)
						event.frac[i]->Fill(event.time_diff[i][j], par[4]);
					else if (event.time_diff[i][j] > event.center_x_split && event.time_diff[i][j] < event.right_x_split)
						event.frac[i]->Fill(event.time_diff[i][j], par[5]);
					else if (event.time_diff[i][j] > event.right_x_split)
						event.frac[i]->Fill(event.time_diff[i][j], par[6]);
				}
				else
				{
					event.frac[i]->Fill(event.time_diff[i][j]);
				}
			}
		}

		for (UInt_t j = 0; j < event.time_diff_data.size(); j++)
		{
			event.data->Fill(event.time_diff_data[j]);
		}

		event.frac[0]->Scale(par[2] * event.frac[0]->GetEntries() / event.frac[0]->Integral(0, nbins + 1));

		if (check_corr == true)
		{
			for (Int_t i = 0; i < nbins; i++)
			{
				event.frac[0]->SetBinContent(i + 1, event.frac[0]->GetBinContent(i + 1) * event.corr_vals[i]);
			}
		}

		event.frac[2]->Scale(par[7] * event.frac[2]->GetEntries() / event.frac[2]->Integral(0, nbins + 1));
		event.frac[3]->Scale(par[8] * event.frac[3]->GetEntries() / event.frac[3]->Integral(0, nbins + 1));
		event.frac[4]->Scale(par[9] * event.frac[4]->GetEntries() / event.frac[4]->Integral(0, nbins + 1));
		event.frac[5]->Scale(par[10] * event.frac[5]->GetEntries() / event.frac[5]->Integral(0, nbins + 1));

		for (UInt_t i = 0; i < channNum; i++)
		{
			event.mc_sum->Add(event.frac[i]);

			event.frac[i]->SetLineWidth(3);
			event.frac[i]->SetLineColor(channColor[i]);
		}

		event.mc_sum->SetLineWidth(3);
		event.mc_sum->SetLineColor(mcSumColor);

		event.data->SetLineWidth(3);
		event.data->SetLineColor(dataColor);

		TCanvas *c1 = new TCanvas("c1", "", 790, 1200);

		c1->SetBottomMargin(0.5);
		c1->Draw();

		TPad *padup_c1 = new TPad("pad_up_c1", "", 0.0, 0.3, 1.0, 1.0);
		TPad *paddown_c1 = new TPad("pad_down_c1", "", 0.0, 0.0, 1.0, 0.3);
		paddown_c1->SetBottomMargin(0.3);
		paddown_c1->SetRightMargin(0.10);

		gStyle->SetOptStat(0);

		TGraphAsymmErrors *sig_eff = new TGraphAsymmErrors();

		Double_t xMinRangeDisplay = x_min, xMaxRangeDisplay = x_max;

		sig_eff->Divide(sig_pass, sig_total, "cl=0.683 b(1,1) mode");

		Double_t
				*yArr = sig_eff->GetY(),
				*yErrLArr = sig_eff->GetEYlow(),
				*yErrHArr = sig_eff->GetEYhigh();

		Int_t n_bins = sig_eff->GetN();
		std::vector<Double_t>
				y(yArr, yArr + n_bins),
				eyl(yErrLArr, yErrLArr + n_bins),
				eyh(yErrHArr, yErrHArr + n_bins);

		Double_t
				weightedMean = Obj.WeightedAverageAsymmetric(y, eyl, eyh),
				weightedMeanErr = Obj.WAvgAsymmError(eyl, eyh);

		std::cout << "Weighted mean: " << weightedMean << " +/- " << weightedMeanErr << std::endl;

		TLine
				*lineAvg = new TLine(xMinRangeDisplay, weightedMean, xMaxRangeDisplay, weightedMean),
				*lineDown = new TLine(xMinRangeDisplay, weightedMean - weightedMeanErr, xMaxRangeDisplay, weightedMean - weightedMeanErr),
				*lineUp = new TLine(xMinRangeDisplay, weightedMean + weightedMeanErr, xMaxRangeDisplay, weightedMean + weightedMeanErr);

		lineAvg->SetLineColor(kRed);
		lineDown->SetLineColor(kRed);
		lineUp->SetLineColor(kRed);
		lineAvg->SetLineStyle(2); // np. przerywana
		lineDown->SetLineStyle(2);
		lineUp->SetLineStyle(2);
		lineAvg->SetLineWidth(2);
		lineDown->SetLineWidth(2);
		lineUp->SetLineWidth(2);

		c1->cd();
		paddown_c1->Draw();
		paddown_c1->cd();
		sig_eff->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);
		sig_eff->GetYaxis()->SetRangeUser(0.0, 1.0);
		sig_eff->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
		sig_eff->GetYaxis()->SetTitle("#varepsilon(#Deltat)");

		sig_eff->GetYaxis()->SetTitleSize(0.1);
		sig_eff->GetYaxis()->SetTitleOffset(0.5);
		sig_eff->GetXaxis()->SetTitleSize(0.1);
		sig_eff->GetYaxis()->SetLabelSize(0.05);
		sig_eff->GetXaxis()->SetLabelSize(0.08);

		sig_eff->Draw("APE");
		lineAvg->Draw("SAME");
		lineDown->Draw("SAME");
		lineUp->Draw("SAME");

		event.data->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);
		event.mc_sum->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);

		event.frac[0]->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);

		TRatioPlot *rp = new TRatioPlot(event.mc_sum, event.data, "diffsig");

		c1->cd();
		padup_c1->Draw();
		padup_c1->cd();

		rp->Draw();

		rp->SetSplitFraction(0.2);

		rp->GetLowerRefGraph()->SetMinimum(-5);
		rp->GetLowerRefGraph()->SetMaximum(5);
		rp->GetLowerRefGraph()->SetLineWidth(3);

		rp->SetLowBottomMargin(0.0);
		rp->SetLeftMargin(0.15);

		rp->GetLowerRefYaxis()->SetLabelSize(0.02);
		rp->GetLowerRefXaxis()->SetLabelSize(0.0);

		Double_t max_height = event.data->GetMaximum();

		rp->GetUpperRefYaxis()->SetRangeUser(0.0, 2 * max_height);
		rp->GetUpperRefYaxis()->SetTitle("Counts/2#tau_{S}");

		rp->GetLowerRefYaxis()->SetTitleSize(0.03);
		rp->GetLowerRefYaxis()->SetTitle("Residuals");

		rp->GetUpperPad()->cd();

		TLegend *legend_chann = new TLegend(0.6, 0.7, 0.9, 0.9);
		legend_chann->SetFillColor(kWhite);
		for (UInt_t i = 0; i < channNum; i++)
		{
			legend_chann->AddEntry(event.frac[i], channName[i], "l");
			event.frac[i]->Draw("HISTSAME");
		}

		legend_chann->AddEntry(event.frac[0], mcSumName, "l");
		legend_chann->AddEntry(event.data, dataName, "le");
		legend_chann->Draw();

		c1->Print(cpfit_dir + img_dir + "split_fit_with_corr" + ext_img);

		// Residuals graph
		TCanvas *c2 = new TCanvas("c2", "", 790, 790);
		TH1 *residuals_hist = new TH1D("Residuals hist", "", 11, -5., 5.);

		event.resi_vals = rp->GetLowerRefGraph()->GetY();

		for (Int_t i = 0; i < event.bin_number; i++)
		{
			residuals_hist->Fill(event.resi_vals[i]);
		}

		residuals_hist->GetYaxis()->SetRangeUser(0, 1.2 * residuals_hist->GetMaximum());
		residuals_hist->Fit("gaus");

		c2->cd();

		residuals_hist->SetXTitle("Residuals");
		residuals_hist->SetYTitle("Counts");
		residuals_hist->SetLineWidth(5);
		residuals_hist->Draw();
		residuals_hist->SetStats(1);
		gStyle->SetOptStat(1);
		gStyle->SetOptFit(1);

		residuals_hist->Draw();

		c2->Print(cpfit_dir + img_dir + "residuals_hist" + ext_img);

		properties["variables"]["CPFit"]["result"]["value"]["Re"] = par[0];
		properties["variables"]["CPFit"]["result"]["value"]["Im"] = par[1];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Signal"] = par[2];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["FarLeft"] = par[3];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["CloseLeft"] = par[4];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["CloseRight"] = par[5];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["FarRight"] = par[6];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Omegapi0"] = par[7];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Threepi0"] = par[8];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Semileptonic"] = par[9];
		properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Other"] = par[10];

		properties["variables"]["CPFit"]["result"]["error"]["Re"] = parErr[0];
		properties["variables"]["CPFit"]["result"]["error"]["Im"] = parErr[1];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Signal"] = parErr[2];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["FarLeft"] = parErr[3];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["CloseLeft"] = parErr[4];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["CloseRight"] = parErr[5];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["FarRight"] = parErr[6];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Omegapi0"] = parErr[7];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Threepi0"] = parErr[8];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Semileptonic"] = parErr[9];
		properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Other"] = parErr[10];

		properties["variables"]["CPFit"]["result"]["chi2"] = event.data->Chi2Test(event.mc_sum, "UW CHI2");
		properties["variables"]["CPFit"]["result"]["normChi2"] = event.data->Chi2Test(event.mc_sum, "UW CHI2/NDF");

		delete residuals_hist;
		delete c1;
		delete c2;

		delete rp;
	}
	else
	{
		TString
				mode = "FitResultErr",
				imageTitle = cpfit_dir + img_dir + "scan_of_errors_fit" + ext_img,
				realTitle = cpfit_dir + img_dir + "real_errors_vs_value_fit" + ext_img,
				imaginaryTitle = cpfit_dir + img_dir + "imaginary_errors_vs_value_fit" + ext_img;

		Double_t legendPos[4] = {0.2, 0.5, 0.7, 0.9};

		KLOE::MeasQualityGraph errorGraph(mode, numberOfPoints, cutLimit, errValue, xTitle, yTitle, yRightTitle, legendPos);

		errorGraph.DrawGraphs(imageTitle);

		mode = "errvsvalue";

		KLOE::MeasQualityGraph realGraph(mode, numberOfPoints, cutLimit, realValue, xTitle, yTitleReal, yRightTitleReal, legendPos);
		KLOE::MeasQualityGraph imaginaryGraph(mode, numberOfPoints, cutLimit, imaginaryValue, xTitle, yTitleImaginary, yRightTitleImaginary, legendPos);

		realGraph.DrawGraphs(realTitle);
		imaginaryGraph.DrawGraphs(imaginaryTitle);
	}

	properties["lastScript"] = "Final CP Parameters normalization";
	properties["lastUpdate"] = Obj.getCurrentTimestamp();

	std::ofstream outfile(propName);
	outfile << properties.dump(4);
	outfile.close();

	return 0;
}