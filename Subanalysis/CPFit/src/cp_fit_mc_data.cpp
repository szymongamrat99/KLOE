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

int cp_fit_mc_data(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
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
	tree_omega->SetBranchAddress("anglePhiOmega", &anglePhiOmega);
	tree_omega->SetBranchAddress("chi2min", &chi2min);

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

	for (UInt_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		if (neutVars.done == 1)
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

			lorentz_transf(Phi_boost, Kch_LAB, Kch_CM);
			lorentz_transf(Phi_boost, Kne_LAB, Kne_CM);
			lorentz_transf(Phi_boost, Kchmom_LAB, Kchmom_CM);
			lorentz_transf(Phi_boost, Knemom_LAB, Knemom_CM);

			Kch_boost[0] = -Kchmom_CM[0] / Kchmom_CM[3];
			Kch_boost[1] = -Kchmom_CM[1] / Kchmom_CM[3];
			Kch_boost[2] = -Kchmom_CM[2] / Kchmom_CM[3];

			Kne_boost[0] = Kchmom_CM[0] / Kchmom_CM[3];
			Kne_boost[1] = Kchmom_CM[1] / Kchmom_CM[3];
			Kne_boost[2] = Kchmom_CM[2] / Kchmom_CM[3];

			lorentz_transf(Kch_boost, Kch_CM, Kch_CMCM);
			lorentz_transf(Kch_boost, Kne_CM, Kne_CMCM);

			baseKin.Dtboostlor = (Kch_CMCM[3] - Kne_CMCM[3]) / (cVel * tau_S_nonCPT);

			for (Int_t i = 0; i < 4; i++)
			{
				TRCV[i] = baseKin.cluster[3][neutVars.gtaken[i]] - (sqrt(pow(baseKin.cluster[0][neutVars.gtaken[i]] - neutVars.Knerec[6], 2) + pow(baseKin.cluster[1][neutVars.gtaken[i]] - neutVars.Knerec[7], 2) + pow(baseKin.cluster[2][neutVars.gtaken[i]] - neutVars.Knerec[8], 2)) / cVel) - tne_LAB;
			}

			trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

			Double_t radius[2] = {0., 0.},
							 radius_ch[2] = {0., 0.};

			Double_t sphere_bound = 10, bp_bound = 4.4;

			Double_t sigma = properties["variables"]["RegenRejection"]["sigma"];

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

			Int_t
					sigmas = 1.;

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
					stdDevRho00IP = sqrt((4-M_PI)/2)*(stdDevOmegaVtx[0] + stdDevOmegaVtx[1])/2.,
					stdDevRhopmIP = sqrt((4-M_PI)/2)*(stdDevOmegaVtx[3] + stdDevOmegaVtx[4])/2.;

			cond[0] = rho_00_IP < sigmas * stdDevRho00IP;
			cond[1] = rho_pm_IP < sigmas * stdDevRhopmIP;
			cond[2] = abs(neu_vtx_avg[2] - baseKin.Kchboost[8]) < sigmas * stdDevOmegaVtx[2];
			cond[3] = abs(baseKin.Kchboost[8] - baseKin.bhabha_vtx[2]) < sigmas * stdDevOmegaVtx[5];

			cond_tot = cond[0] && cond[1] && cond[2] && cond[3];
			// ----------------------------------------------------------------------------------------------------------------

			Double_t
					meanRadiusChHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["mean"][1],
					errorRadiusChHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["width"][1],
					meanRadiusChLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["mean"][0],
					errorRadiusChLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["width"][0],

					meanRadiusTriangleHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["mean"][1],
					errorRadiusTriangleHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["width"][1],
					meanRadiusTriangleLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["mean"][0],
					errorRadiusTriangleLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["width"][0];

			cuts[0] = abs(radius[0] - meanRadiusTriangleHigher) > errorRadiusTriangleHigher * sigma; // && radius[1] > 8;
			cuts[1] = 1;																																						 // abs(radius[1] - meanRadiusTriangleLower) > errorRadiusTriangleLower * sigma;     // && radius[1] <= 8;
			cuts[2] = abs(radius_ch[0] - meanRadiusChHigher) > errorRadiusChHigher * sigma;					 // && radius_ch[1] > 8;
			cuts[3] = 1;																																						 // abs(radius_ch[1] - meanRadiusChLower) > errorRadiusChLower * sigma;  // && radius_ch[1] <= 8;

			if (cond_tot)
			{
				cuts[4] = 1; //rho > 0.6;
				cuts[5] = 1; // chi2min > 4.0;
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

				if (cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5])
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

			if (baseKin.mcflag == 0 && cuts[0] && cuts[1] && cuts[2] && cuts[3] && (cuts[4] && cuts[5]))
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
								 step[num_of_vars] = {1E-5, 1E-5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

	Double_t
			limit_upper = properties["variables"]["CPFit"]["limitPer"]["upper"],
			limit_lower = properties["variables"]["CPFit"]["limitPer"]["lower"];

	minimum->SetVariable(0, "Real part", init_vars[0], step[0]);
	minimum->SetVariable(1, "Imaginary part", init_vars[1], step[1]);
	minimum->SetLimitedVariable(2, "Norm signal", init_vars[2], step[2], init_vars[2] - 1.0 * init_vars[2], init_vars[2] + 5.0 * init_vars[2]);

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

	Double_t par[2] = {minimum->X()[0], minimum->X()[1]};

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

				event.frac[i]->Fill(event.time_diff[i][j], event.interf_function(event.time_diff_gen[j], 0, par));
			}
			else if (i == 1)
			{
				if (event.time_diff[i][j] < event.left_x_split)
					event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[3]);
				else if (event.time_diff[i][j] > event.left_x_split && event.time_diff[i][j] < event.center_x_split)
					event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[4]);
				else if (event.time_diff[i][j] > event.center_x_split && event.time_diff[i][j] < event.right_x_split)
					event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[5]);
				else if (event.time_diff[i][j] > event.right_x_split)
					event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[6]);
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

	event.frac[0]->Scale(minimum->X()[2] * event.frac[0]->GetEntries() / event.frac[0]->Integral(0, nbins + 1));

	if (check_corr == true)
	{
		for (Int_t i = 0; i < nbins; i++)
		{
			event.frac[0]->SetBinContent(i + 1, event.frac[0]->GetBinContent(i + 1) * event.corr_vals[i]);
		}
	}

	event.frac[2]->Scale(minimum->X()[7] * event.frac[2]->GetEntries() / event.frac[2]->Integral(0, nbins + 1));
	event.frac[3]->Scale(minimum->X()[8] * event.frac[3]->GetEntries() / event.frac[3]->Integral(0, nbins + 1));
	event.frac[4]->Scale(minimum->X()[9] * event.frac[4]->GetEntries() / event.frac[4]->Integral(0, nbins + 1));
	event.frac[5]->Scale(minimum->X()[10] * event.frac[5]->GetEntries() / event.frac[5]->Integral(0, nbins + 1));

	for (UInt_t i = 0; i < channNum; i++)
	{
		event.mc_sum->Add(event.frac[i]);

		event.frac[i]->SetLineWidth(3);
		event.frac[i]->SetLineColor(channColor[i]);
	}

	properties["variables"]["CPFit"]["result"]["value"]["Re"] = minimum->X()[0];
	properties["variables"]["CPFit"]["result"]["value"]["Im"] = minimum->X()[1];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Signal"] = minimum->X()[2];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["FarLeft"] = minimum->X()[3];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["CloseLeft"] = minimum->X()[4];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["CloseRight"] = minimum->X()[5];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["FarRight"] = minimum->X()[6];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Omegapi0"] = minimum->X()[7];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Threepi0"] = minimum->X()[8];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Semileptonic"] = minimum->X()[9];
	properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Other"] = minimum->X()[10];

	properties["variables"]["CPFit"]["result"]["error"]["Re"] = minimum->Errors()[0];
	properties["variables"]["CPFit"]["result"]["error"]["Im"] = minimum->Errors()[1];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Signal"] = minimum->Errors()[2];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["FarLeft"] = minimum->Errors()[3];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["CloseLeft"] = minimum->Errors()[4];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["CloseRight"] = minimum->Errors()[5];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["FarRight"] = minimum->Errors()[6];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Omegapi0"] = minimum->Errors()[7];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Threepi0"] = minimum->Errors()[8];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Semileptonic"] = minimum->Errors()[9];
	properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Other"] = minimum->Errors()[10];

	properties["variables"]["CPFit"]["result"]["chi2"] = event.data->Chi2Test(event.mc_sum, "UW CHI2");
	properties["variables"]["CPFit"]["result"]["normChi2"] = event.data->Chi2Test(event.mc_sum, "UW CHI2/NDF");

	properties["lastUpdate"] = Obj.getCurrentDate();

	std::ofstream outfile(propName);
	outfile << properties.dump(4);
	outfile.close();

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

	gStyle->SetOptStat(0);

	TGraphAsymmErrors *sig_eff = new TGraphAsymmErrors();

	sig_eff->Divide(sig_pass, sig_total, "cl=0.683 b(1,1) mode");

	const Int_t
			n = sig_eff->GetN();

	Int_t
			n_mean = 0;

	Double_t ax[n], ay[n], mean = 0.;

	std::cout << "Tau_S" << " " << "Efficiency" << std::endl;

	for (Int_t k = 0; k < n; k++)
	{
		sig_eff->GetPoint(k, ax[k], ay[k]);
		if (abs(ax[k]) < 5.)
		{
			n_mean++;
			mean += ay[k];
			std::cout << ax[k] << " " << ay[k] << std::endl;
		}
	}

	mean = mean / (Float_t)n_mean;

	std::cout << mean << std::endl;

	c1->cd();
	paddown_c1->Draw();
	paddown_c1->cd();
	sig_eff->GetXaxis()->SetRangeUser(x_min, x_max);
	sig_eff->GetYaxis()->SetRangeUser(0.0, 1.0);
	sig_eff->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
	sig_eff->GetYaxis()->SetTitle("#varepsilon(#Deltat)");

	sig_eff->GetYaxis()->SetTitleSize(0.1);
	sig_eff->GetYaxis()->SetTitleOffset(0.5);
	sig_eff->GetXaxis()->SetTitleSize(0.1);
	sig_eff->GetYaxis()->SetLabelSize(0.05);
	sig_eff->GetXaxis()->SetLabelSize(0.08);

	sig_eff->Draw("APE");

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

	Double_t max_height = event.mc_sum->GetMaximum();

	rp->GetUpperRefYaxis()->SetRangeUser(0.0, 2 * max_height);
	rp->GetUpperRefYaxis()->SetTitle("Counts/2#tau_{S}");

	rp->GetLowerRefYaxis()->SetTitleSize(0.03);
	rp->GetLowerRefYaxis()->SetTitle("Residuals");

	rp->GetUpperPad()->cd();

	TLegend *legend_chann = new TLegend(0.6, 0.5, 0.9, 0.9);
	legend_chann->SetFillColor(kWhite);
	for (UInt_t i = 0; i < channNum; i++)
	{
		legend_chann->AddEntry(event.frac[i], channName[i], "l");
		event.frac[i]->Draw("HISTSAME");
	}

	legend_chann->AddEntry(event.mc_sum, mcSumName, "l");
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

	delete residuals_hist;
	delete c1;
	delete c2;

	delete rp;

	return 0;
}