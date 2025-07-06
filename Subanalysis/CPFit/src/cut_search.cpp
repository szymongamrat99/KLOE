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


#include "../inc/cpfit.hpp"

int cut_search(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
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

	const TString cpfit_res_dir = SystemPath::cpfit_dir + SystemPath::result_dir;

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

	// Initialization of cut formulas
	Int_t cuts_number = 6;

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

	cutStep = abs(cutLimits[1] - cutLimits[0]) / (Double_t)numberOfPoints;

	// Vectors needed for the initialization of graphs
	TString graphMode = "FitResultErr";
	std::vector<Double_t>
			cutLimit;
	std::vector<std::vector<Double_t>>
			errValue(2),
			realValue(2),
			imaginaryValue(2);

	// Values of cuts to be used in the analysis
	Double_t
			rho_cut = 0.0;

	std::vector<Double_t> par, parErr;

	Double_t
			sigmas = 1,
			sigma = (Double_t)properties["variables"]["RegenRejection"]["sigma"],
			ChHigher[2] = {
					properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["mean"][1],
					(Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["width"][1]},
			ChLower[2] = {properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["mean"][0], (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["width"][0]}, TrHigher[2] = {properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["mean"][1], (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["width"][1]}, TrLower[2] = {properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["mean"][0], (Double_t)properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["width"][0]}, OmegaCut1[2], OmegaCut2[2];

	// rho_cut = properties["variables"]["CPFit"]["cuts"]["omegaRho"];

	TrHigher[1] = 1.5;
	ChHigher[1] = 1.5;

	std::map<Int_t, std::vector<Event>> dataPoint;

	for (UInt_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		if (neutVars.done == 1 && doneOmega == 1)
		{
			velocity_kch = cVel * sqrt(pow(baseKin.Kchboost[0], 2) + pow(baseKin.Kchboost[1], 2) + pow(baseKin.Kchboost[2], 2)) / baseKin.Kchboost[3];

			velocity_kne = cVel * sqrt(pow(neutVars.Knerec[0], 2) + pow(neutVars.Knerec[1], 2) + pow(neutVars.Knerec[2], 2)) / neutVars.Knerec[3];

			tch_LAB = sqrt(pow(baseKin.Kchboost[6] - baseKin.ip[0], 2) + pow(baseKin.Kchboost[7] - baseKin.ip[1], 2) + pow(baseKin.Kchboost[8] - baseKin.ip[2], 2)) / velocity_kch;
			tne_LAB = sqrt(pow(neutVars.Knerec[6] - baseKin.ip[0], 2) + pow(neutVars.Knerec[7] - baseKin.ip[1], 2) + pow(neutVars.Knerec[8] - baseKin.ip[2], 2)) / velocity_kne;

			Kch_LAB[0] = baseKin.Kchboost[6] - baseKin.ip[0];
			Kch_LAB[1] = baseKin.Kchboost[7] - baseKin.ip[1];
			Kch_LAB[2] = baseKin.Kchboost[8] - baseKin.ip[2];
			Kch_LAB[3] = tch_LAB * cVel;

			Kchmom_LAB[0] = baseKin.Kchboost[0];
			Kchmom_LAB[1] = baseKin.Kchboost[1];
			Kchmom_LAB[2] = baseKin.Kchboost[2];
			Kchmom_LAB[3] = baseKin.Kchboost[3];

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

			// ----------------------------------------------------------------------------------------------------------------

			Event myEvent;

			// -----------------------------------------------------------------------------
			if (baseKin.mcflag == 1 && baseKin.mctruth_int != 0)
			{
				// Filling of reconstructed Delta T
				myEvent.DeltaTRec = baseKin.Dtboostlor;
				if (baseKin.mctruth_int == 1)
					myEvent.DeltaTGen = baseKin.Dtmc;

				myEvent.rho00 = rho_00_IP - meanRho00IP < stdDevRho00IP;
				myEvent.rhopm = rho_pm_IP - meanRhopmIP < stdDevRhopmIP;
				myEvent.z00 =
						abs(neu_vtx_avg[2] - baseKin.Kchboost[8]) < stdDevOmegaVtx[2];
				myEvent.zpm =
						abs(baseKin.Kchboost[8] - baseKin.bhabha_vtx[2]) < stdDevOmegaVtx[5];
				myEvent.totalOmegaFidVol =
						myEvent.rho00 &&
						myEvent.rhopm &&
						myEvent.z00 &&
						myEvent.zpm;

				myEvent.MinvOmega = Omegarec[5];

				myEvent.KinEnePi0 = Omegarec[6];

				myEvent.beamPipeCutNeg = radius[0];
				myEvent.beamPipeCutPos = radius_ch[0];

				myEvent.mctruth = baseKin.mctruth_int;
				myEvent.mcflag = 1;

				dataPoint[baseKin.mctruth_int].push_back(myEvent);
			}
			else if (baseKin.mcflag == 0)
			{
				// Filling of reconstructed Delta T
				myEvent.DeltaTRec = baseKin.Dtboostlor;

				myEvent.rho00 = rho_00_IP - meanRho00IP < stdDevRho00IP;
				myEvent.rhopm = rho_pm_IP - meanRhopmIP < stdDevRhopmIP;
				myEvent.z00 =
						abs(neu_vtx_avg[2] - baseKin.Kchboost[8]) < stdDevOmegaVtx[2];
				myEvent.zpm =
						abs(baseKin.Kchboost[8] - baseKin.bhabha_vtx[2]) < stdDevOmegaVtx[5];
				myEvent.totalOmegaFidVol =
						myEvent.rho00 &&
						myEvent.rhopm &&
						myEvent.z00 &&
						myEvent.zpm;

				myEvent.MinvOmega = Omegarec[5];

				myEvent.KinEnePi0 = Omegarec[6];

				myEvent.beamPipeCutNeg = radius[0];
				myEvent.beamPipeCutPos = radius_ch[0];

				myEvent.mctruth = 0;
				myEvent.mcflag = 0;

				dataPoint[baseKin.mctruth_int].push_back(myEvent);
			}
			// -----------------------------------------------------------------------------
		}
	}

	// Omega cuts
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
			KinEne[2] = {meanKinEne, 3 * stdKinEne},
			a = properties["variables"]["OmegaRec"]["combined"]["line"]["slope"],
			b = properties["variables"]["OmegaRec"]["combined"]["line"]["inter"],
			lineWidth = properties["variables"]["OmegaRec"]["combined"]["stdDev"]["value"];
	// -------------------------------------------------------------------------------

	Bool_t
			cut_regen_neg,
			cut_regen_pos,
			cut_line_down,
			cut_line_up,
			cut_KinEnePi0,
			cut_MinvOmega,
			cut_total;

	std::vector<std::vector<Double_t>>
			relativeErr(2),
			real(2),
			imaginary(2);

	Double_t k = cutLimits[0];

	for (Int_t i = 0; i < numberOfPoints; i++)
	{
		k += cutStep;

		cutLimit.push_back(k);

		std::cout << i << std::endl;

		for (const auto &pair : dataPoint)
			for (const auto &ev : pair.second)
			{
				if (ev.totalOmegaFidVol)
				{
					cut_MinvOmega = abs(ev.MinvOmega - meanInvMass) < k * stdInvMass;
					cut_KinEnePi0 = abs(ev.KinEnePi0 - meanKinEne) < k * stdKinEne;

					Double_t normFactor = k * lineWidth / sqrt(a * a + 1); // Norm of normal vector
					cut_line_up = (a * ev.KinEnePi0 + b + (1 - pow(a, 2)) * normFactor) > ev.MinvOmega;
					cut_line_down = (a * ev.KinEnePi0 + b - (1 - pow(a, 2)) * normFactor) < ev.MinvOmega;
				}
				else
				{
					cut_MinvOmega = 0;
					cut_KinEnePi0 = 0;
					cut_line_down = 0;
					cut_line_up = 0;
				}

				cut_regen_neg = abs(ev.beamPipeCutNeg - TrHigher[0]) < sigma * TrHigher[1];
				cut_regen_pos = abs(ev.beamPipeCutPos - ChHigher[0]) < sigma * ChHigher[1];

				cut_total = !(cut_regen_neg || cut_regen_pos) && !(cut_MinvOmega && cut_KinEnePi0 && cut_line_up && cut_line_down);

				if (cut_total)
				{
					if (ev.mcflag == 1)
					{
						if (ev.mctruth == 1)
						{
							event.time_diff_gen.push_back(ev.DeltaTGen);
							event.time_diff[0].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 3)
						{
							event.time_diff[1].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 4)
						{
							event.time_diff[2].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 5)
						{
							event.time_diff[3].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 6)
						{
							event.time_diff[4].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 7)
						{
							event.time_diff[5].push_back(ev.DeltaTRec);
						}
					}

					if (ev.mcflag == 0)
					{
						event.time_diff_data.push_back(ev.DeltaTRec);
					}
				}
			}

		cp_fit_func(event, relativeErr, real, imaginary, logger, Obj);

		event.ClearVectors();
	}

	std::vector<std::vector<Double_t>> cutLimitCombined;
	std::vector<std::vector<Double_t>>
			relativeErrCombined(2),
			realCombined(2),
			imaginaryCombined(2);

	for (Int_t i1 = 1; i1 <= numberOfPoints; i1++)
		for (Int_t i2 = 1; i2 <= numberOfPoints; i2++)
			for (Int_t i3 = 1; i3 <= numberOfPoints; i3++)
			{
				std::vector<Double_t> cutVal(3);
				cutVal[0] = i1 * cutStep;
				cutVal[1] = i2 * cutStep;
				cutVal[2] = i3 * cutStep;

				cutLimitCombined.push_back(cutVal);
			}

	/*for (Int_t i = 0; i < cutLimitCombined.size(); i++)
	{
		std::cout << cutLimitCombined[i][0] << " " << cutLimitCombined[i][1] << " " << cutLimitCombined[i][2] << std::endl;

		for (const auto &pair : dataPoint)
			for (const auto &ev : pair.second)
			{
				if (ev.totalOmegaFidVol)
				{
					cut_MinvOmega = abs(ev.MinvOmega - meanInvMass) < cutLimitCombined[i][0] * stdInvMass;
					cut_KinEnePi0 = abs(ev.KinEnePi0 - meanKinEne) < cutLimitCombined[i][1] * stdKinEne;

					Double_t normFactor = cutLimitCombined[i][2] * lineWidth / sqrt(a * a + 1); // Norm of normal vector
					cut_line_up = (a * ev.KinEnePi0 + b + (1 - pow(a, 2)) * normFactor) > ev.MinvOmega;
					cut_line_down = (a * ev.KinEnePi0 + b - (1 - pow(a, 2)) * normFactor) < ev.MinvOmega;
				}
				else
				{
					cut_MinvOmega = 0;
					cut_KinEnePi0 = 0;
					cut_line_down = 0;
					cut_line_up = 0;
				}

				cut_regen_neg = abs(ev.beamPipeCutNeg - TrHigher[0]) < sigma * TrHigher[1];
				cut_regen_pos = abs(ev.beamPipeCutPos - ChHigher[0]) < sigma * ChHigher[1];

				cut_total = !(cut_regen_neg && cut_regen_pos) && !(cut_MinvOmega && cut_KinEnePi0 && cut_line_up && cut_line_down);

				if (cut_total)
				{
					if (ev.mcflag == 1)
					{
						if (ev.mctruth == 1)
						{
							event.time_diff_gen.push_back(ev.DeltaTGen);
							event.time_diff[0].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 3)
						{
							event.time_diff[1].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 4)
						{
							event.time_diff[2].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 5)
						{
							event.time_diff[3].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 6)
						{
							event.time_diff[4].push_back(ev.DeltaTRec);
						}
						else if (ev.mctruth == 7)
						{
							event.time_diff[5].push_back(ev.DeltaTRec);
						}
					}

					if (ev.mcflag == 0)
					{
						event.time_diff_data.push_back(ev.DeltaTRec);
					}
				}
			}

		cp_fit_func(event, relativeErrCombined, realCombined, imaginaryCombined, logger, Obj);

		std::cout << relativeErrCombined[0][i] << " " << relativeErrCombined[1][i] << std::endl;
		std::cout << realCombined[0][i] << " " << realCombined[1][i] << std::endl;
		std::cout << imaginaryCombined[0][i] << " " << imaginaryCombined[1][i] << std::endl;

		event.ClearVectors();
	}*/

	TString xTitle = "Number of standard deviations [-]",
					yTitle = "Relative error of parameter",
					yRightTitle = "", //"|#sigma(Im(#varepsilon'/#varepsilon))/Im(#varepsilon'/#varepsilon)|",
			yTitleReal = "Re(#varepsilon'/#varepsilon)",
					yRightTitleReal = "#sigma(Re(#varepsilon'/#varepsilon))",
					yTitleImaginary = "Im(#varepsilon'/#varepsilon)",
					yRightTitleImaginary = "#sigma(Im(#varepsilon'/#varepsilon))",
					cutName = (std::string)properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["cutName"];

	TString
			modeGraph = "FitResultErr",
			imageTitle = SystemPath::cpfit_dir + SystemPath::img_dir + "scan_of_errors_fit" + ext_img,
			realTitle = SystemPath::cpfit_dir + SystemPath::img_dir + "real_errors_vs_value_fit" + ext_img,
			imaginaryTitle = SystemPath::cpfit_dir + SystemPath::img_dir + "imaginary_errors_vs_value_fit" + ext_img;

	Double_t legendPos[4] = {0.2, 0.5, 0.7, 0.9};

	KLOE::MeasQualityGraph errorGraph(modeGraph, numberOfPoints, cutLimit, relativeErr, xTitle, yTitle, yRightTitle, legendPos);

	errorGraph.DrawGraphs(imageTitle);

	modeGraph = "errvsvalue";

	KLOE::MeasQualityGraph realGraph(modeGraph, numberOfPoints, cutLimit, real, xTitle, yTitleReal, yRightTitleReal, legendPos);
	KLOE::MeasQualityGraph imaginaryGraph(modeGraph, numberOfPoints, cutLimit, imaginary, xTitle, yTitleImaginary, yRightTitleImaginary, legendPos);

	realGraph.DrawGraphs(realTitle);
	imaginaryGraph.DrawGraphs(imaginaryTitle);

	properties["lastScript"] = "Final CP Parameters normalization";
	properties["lastUpdate"] = Obj.getCurrentTimestamp();

	std::ofstream outfile(propName);
	outfile << properties.dump(4);
	outfile.close();

	return 0;
}