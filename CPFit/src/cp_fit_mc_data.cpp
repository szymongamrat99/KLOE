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
#include "../../../Include/Codes/interference.h"
#include "../../../Include/Codes/kloe_class.h"
#include "../../../Include/Codes/lorentz_transf.h"

#define __FILENAME__ (__builtin_strrchr(__FILE__, '/') ? __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)

using namespace KLOE;

int cp_fit_mc_data(Int_t firstFile, Int_t lastFile, TString mode = "split", Bool_t check_corr = false, Int_t loopcount = 0, Int_t M = 0, Int_t range = 0)
{
	ErrorHandling::ErrorLogs errLogger;
	Controls::DataType data_type = Controls::DataType::MC_DATA;

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
			*file;
	TTree
			*tree,
			*tree_mctruth;
	// =============================================================================

	TString name = gen_vars_dir + root_files_dir + mctruth_filename + firstFile + "_" + lastFile + ext_root;

	file_mctruth = new TFile(name);
	tree_mctruth = (TTree *)file_mctruth->Get(gen_vars_tree);

	tree_mctruth->SetBranchAddress("mctruth", &baseKin.mctruth_int);

	Double_t *eff_vals;

	if (check_corr == true)
	{
		file_corr = new TFile("../Efficiency_analysis/correction_factor.root");
		TGraphAsymmErrors *eff_signal = (TGraphAsymmErrors *)file_corr->Get("correction_factor");
		eff_vals = eff_signal->GetY();

		delete eff_signal;
	}

	// ===========================================================================

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, firstFile, lastFile);

	chain->SetBranchAddress("mcflag", &baseKin.mcflag);

	chain->SetBranchAddress("Dtmc", &baseKin.Dtmc);

	chain->SetBranchAddress("Chi2", &baseKin.Chi2);
	chain->SetBranchAddress("minv4gam", &baseKin.minv4gam);
	chain->SetBranchAddress("Kchboost", baseKin.Kchboost);
	// chain->SetBranchAddress("Knereclor", neutVars.Knerec);
	chain->SetBranchAddress("Qmiss", &baseKin.Qmiss);

	Float_t
			PichFourMom[2][4];

	chain->SetBranchAddress("trk1", PichFourMom[0]);
	chain->SetBranchAddress("trk2", PichFourMom[1]);

	chain->SetBranchAddress("Xcl", baseKin.cluster[0]);
	chain->SetBranchAddress("Ycl", baseKin.cluster[1]);
	chain->SetBranchAddress("Zcl", baseKin.cluster[2]);
	chain->SetBranchAddress("Tcl", baseKin.cluster[3]);
	chain->SetBranchAddress("Enecl", baseKin.cluster[4]);

	chain->SetBranchAddress("ip", baseKin.ip);

	chain->SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
	chain->SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
	chain->SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
	chain->SetBranchAddress("Broots", &baseKin.phi_mom[3]);

	chain->SetBranchAddress("Bx", &baseKin.bhabha_vtx[0]);
	chain->SetBranchAddress("By", &baseKin.bhabha_vtx[1]);
	chain->SetBranchAddress("Bz", &baseKin.bhabha_vtx[2]);

	// ===========================================================================

	TString
			filename_trilateration = gen_vars_dir + root_files_dir + gen_vars_filename + firstFile + "_" + lastFile + ext_root,
			filename_trilateration_kin_fit = neutrec_dir + root_files_dir + neu_trilateration_kin_fit_filename + firstFile + "_" + lastFile + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root,
			filename_triangle = neutrec_dir + root_files_dir + neu_triangle_filename + firstFile + "_" + lastFile + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

	std::vector<TString> file_name, tree_name;

	file_name.push_back(filename_trilateration);
	file_name.push_back(filename_trilateration_kin_fit);
	file_name.push_back(filename_triangle);

	tree_name.push_back(neutrec_tri_tree);
	tree_name.push_back(neutrec_kin_fit_tree);
	tree_name.push_back(neutrec_triangle_tree);

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
	}
	catch (ErrorHandling::ErrorCodes err)
	{
		errLogger.setErrCount(err);
		errLogger.getErrLog(err);
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

	chain->AddFriend(tree_mctruth);
	chain->AddFriend(tree);

	UInt_t nentries = chain->GetEntries();

	Double_t x_min = -90.0, x_max = 90.0;
	UInt_t nbins = 1 + ((x_max - x_min) / 2.);

	Double_t split[3] = {-30.0, 0.0, 30.0};

	interference event(mode, check_corr, nbins, x_min, x_max, split);

	Double_t velocity_kch, velocity_kne, tch_LAB, tch_CM, tne_LAB, tne_CM, trcv_sum, TRCV[4];
	Float_t Kch_LAB[4], Kne_LAB[4], Kch_CM[4], Kne_CM[4], Kch_CMCM[4], Kne_CMCM[4], Kne_boost[3], Kch_boost[3], Phi_boost[3], Kchmom_LAB[4], Knemom_LAB[4], Kchmom_CM[4], Knemom_CM[4];

	TH1 *sig_total = new TH1D("sig_tot", ";#Delta t;Counts", nbins, x_min, x_max);
	TH1 *sig_pass = new TH1D("sig_pass", ";#Delta t;Counts", nbins, x_min, x_max);

	std::vector<Float_t> no_cuts_sig[2];

	for (UInt_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		if (neutVars.done == 1)
		{
			velocity_kch = cVel * sqrt(pow(baseKin.Kchboost[0], 2) + pow(baseKin.Kchboost[1], 2) + pow(baseKin.Kchboost[2], 2)) / baseKin.Kchboost[3];

			velocity_kne = cVel * sqrt(pow(neutVars.Knerec[0], 2) + pow(neutVars.Knerec[1], 2) + pow(neutVars.Knerec[2], 2)) / neutVars.Knerec[3];

			tch_LAB = sqrt(pow(baseKin.Kchboost[6] - baseKin.ip[0], 2) + pow(baseKin.Kchboost[7] - baseKin.ip[1], 2) + pow(baseKin.Kchboost[8] - baseKin.ip[2], 2)) / velocity_kch;
			tne_LAB = sqrt(pow(neutVars.Knerec[6] - baseKin.ip[0], 2) + pow(neutVars.Knerec[7] - baseKin.ip[1], 2) + pow(neutVars.Knerec[8] - baseKin.ip[2], 2)) / velocity_kne; // neutVars.Knerec[9];

			Kch_LAB[0] = baseKin.Kchboost[6] - baseKin.ip[0];
			Kch_LAB[1] = baseKin.Kchboost[7] - baseKin.ip[1];
			Kch_LAB[2] = baseKin.Kchboost[8] - baseKin.ip[2];
			Kch_LAB[3] = tch_LAB * cVel;

			Kchmom_LAB[0] = baseKin.Kchboost[0];
			Kchmom_LAB[1] = baseKin.Kchboost[1];
			Kchmom_LAB[2] = baseKin.Kchboost[2];
			Kchmom_LAB[3] = baseKin.Kchboost[3];

			Kne_LAB[0] = neutVars.Knerec[6] + 0.07 - baseKin.ip[0];
			Kne_LAB[1] = neutVars.Knerec[7] - 0.01 - baseKin.ip[1];
			Kne_LAB[2] = neutVars.Knerec[8] + 0.03 - baseKin.ip[2];
			Kne_LAB[3] = (tne_LAB - 1.3 * tau_S_nonCPT) * cVel;

			Knemom_LAB[0] = neutVars.Knerec[0] + 0.31;
			Knemom_LAB[1] = neutVars.Knerec[1] - 0.09;
			Knemom_LAB[2] = neutVars.Knerec[2] + 0.07;
			Knemom_LAB[3] = neutVars.Knerec[3] - 0.06;

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

			Float_t
					boost_vec_Kchboost[3] = {-baseKin.Kchboost[0] / baseKin.Kchboost[3],
																	 -baseKin.Kchboost[1] / baseKin.Kchboost[3],
																	 -baseKin.Kchboost[2] / baseKin.Kchboost[3]},
					PichFourMomKaonCM[2][4],
					anglePichKaonCM;

			lorentz_transf(boost_vec_Kchboost, PichFourMom[0], PichFourMomKaonCM[0]);
			lorentz_transf(boost_vec_Kchboost, PichFourMom[1], PichFourMomKaonCM[1]);

			TVector3
					pich1(PichFourMomKaonCM[0][0], PichFourMomKaonCM[0][1], PichFourMomKaonCM[0][2]),
					pich2(PichFourMomKaonCM[1][0], PichFourMomKaonCM[1][1], PichFourMomKaonCM[1][2]);

			anglePichKaonCM = pich1.Angle(pich2) * 180. / M_PI;

			Double_t radius[2] = {0., 0.},
							 radius_ch[2] = {0., 0.},
							 rho_pm_IP = 0.,
							 rho_00_IP = 0.;

			Double_t sphere_bound = 10, bp_bound = 4.4;

			Double_t sigma = 1.5;

			Bool_t cuts[4] = {false, false, false, false};

			rho_pm_IP = sqrt(pow(baseKin.Kchboost[6] - baseKin.bhabha_vtx[0], 2) + pow(baseKin.Kchboost[7] - baseKin.bhabha_vtx[1], 2));
			rho_00_IP = sqrt(pow(neutVars.Knerec[6] - baseKin.bhabha_vtx[0], 2) + pow(neutVars.Knerec[7] - baseKin.bhabha_vtx[1], 2));

			for (Int_t i = 0; i < 3; i++)
			{
				radius[0] += pow(neutVars.Knerec[6 + i], 2);
				radius_ch[0] += pow(baseKin.Kchboost[6 + i], 2);

				if (i < 2)
				{
					radius[1] += pow(neutVars.Knerec[6 + i], 2);
					radius_ch[1] += pow(baseKin.Kchboost[6 + i], 2);
				}
			}

			radius[0] = sqrt(radius[0]);
			radius[1] = sqrt(radius[1]);

			radius_ch[0] = sqrt(radius_ch[0]);
			radius_ch[1] = sqrt(radius_ch[1]);

			Int_t
					sigmas_neu = 1.,
					sigmas_ip = 3.;

			Bool_t
					cond[6],
					cond_tot;

			cond[0] = abs(baseKin.Kchboost[6] - neutVars.Knerec[6]) < sigmas_neu * 16.89;
			cond[1] = abs(baseKin.Kchboost[7] - neutVars.Knerec[7]) < sigmas_neu * 16.77;
			cond[2] = abs(baseKin.Kchboost[8] - neutVars.Knerec[8]) < sigmas_neu * 9.89;
			cond[3] = abs(baseKin.Kchboost[6] - baseKin.bhabha_vtx[0]) < sigmas_ip * 2.63;
			cond[4] = abs(baseKin.Kchboost[7] - baseKin.bhabha_vtx[1]) < sigmas_ip * 2.03;
			cond[5] = abs(baseKin.Kchboost[8] - baseKin.bhabha_vtx[2]) < sigmas_ip * 3.39;

			cond_tot = cond[0] && cond[1] && cond[2] && cond[3] && cond[4] && cond[5];

			cuts[0] = abs(radius[0] - 10.9554) > 1.19279 * sigma;
			cuts[1] = 1; // abs(radius[1] - 4.5) > 1.5 * sigma;
			cuts[2] = abs(radius_ch[0] - 10.5840) > 0.681656 * sigma;
			cuts[3] = 1; // abs(radius_ch[1] - 4.45238) > 1.29837 * sigma;

			if (cond_tot)
			{
				cuts[4] = 1;//sqrt(pow(rho_00_IP, 2) + pow(rho_pm_IP, 2)) > 0.5;
				cuts[5] = 1;//anglePichKaonCM > 179.5;
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

	ROOT::Math::Functor minimized_function(&event, &interference::interf_chi2, num_of_vars);

	minimum->SetFunction(minimized_function);

	const Double_t init_vars[num_of_vars] = {Re, M_PI * Im_nonCPT / 180., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
								 step[num_of_vars] = {1E-5, 1E-5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

	Double_t limit_span_signal = 0.3, limit_span = 0.3, limit_pars = 1.0;

	minimum->SetVariable(0, "Real part", init_vars[0], step[0]);			//, init_vars[0] - limit_pars*init_vars[0], init_vars[0] + limit_pars*init_vars[0]);
	minimum->SetVariable(1, "Imaginary part", init_vars[1], step[1]); //, init_vars[1] - limit_pars*init_vars[1], init_vars[1] + limit_pars*init_vars[1]);
	minimum->SetLimitedVariable(2, "Norm signal", init_vars[2], step[2], init_vars[2] - limit_span_signal * init_vars[2], init_vars[2] + limit_span_signal * init_vars[2]);

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

	minimum->SetLimitedVariable(7, "Norm omega", init_vars[7], step[7], init_vars[7] - limit_span * init_vars[7], init_vars[7] + limit_span * init_vars[7]);
	minimum->SetLimitedVariable(8, "Norm three", init_vars[8], step[8], init_vars[8] - limit_span * init_vars[8], init_vars[8] + limit_span * init_vars[8]);
	minimum->SetLimitedVariable(9, "Norm semi", init_vars[9], step[9], init_vars[9] - limit_span * init_vars[9], init_vars[9] + limit_span * init_vars[9]);
	minimum->SetLimitedVariable(10, "Norm other bcg", init_vars[10], step[10], init_vars[10] - limit_span * init_vars[10], init_vars[10] + limit_span * init_vars[10]);

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

	std::ofstream myfile;
	myfile.open(cpfit_res_dir + "results_split_fit.csv");
	myfile << "Parameter,Value,Error\n";
	myfile << "Real part," << minimum->X()[0] << "," << minimum->Errors()[0] << ",\n";
	myfile << "Imaginary part," << minimum->X()[1] << "," << minimum->Errors()[1] << ",\n";
	myfile << "Norm signal," << minimum->X()[2] << "," << minimum->Errors()[2] << ",\n";
	myfile << "Norm left DC wall," << minimum->X()[3] << "," << minimum->Errors()[3] << ",\n";
	myfile << "Norm left beam pipe," << minimum->X()[4] << "," << minimum->Errors()[4] << ",\n";
	myfile << "Norm right beam pipe," << minimum->X()[5] << "," << minimum->Errors()[5] << ",\n";
	myfile << "Norm right DC wall," << minimum->X()[6] << "," << minimum->Errors()[6] << ",\n";
	myfile << "Norm omega," << minimum->X()[7] << "," << minimum->Errors()[7] << ",\n";
	myfile << "Norm three," << minimum->X()[8] << "," << minimum->Errors()[8] << ",\n";
	myfile << "Norm semi," << minimum->X()[9] << "," << minimum->Errors()[9] << ",\n";
	myfile << "Norm other bcg," << minimum->X()[10] << "," << minimum->Errors()[10] << ",\n";
	myfile << "\u03C7\u00B2," << event.data->Chi2Test(event.mc_sum, "UW CHI2") << ",-,\n";
	myfile << "\u03C7\u00B2/" << (UInt_t)nbins << "," << event.data->Chi2Test(event.mc_sum, "UW CHI2/NDF") << ",-,\n";
	myfile.close();

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

	c1->Print(img_dir + "split_fit_with_corr" + ext_img);

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

	c2->Print(img_dir + "residuals_hist" + ext_img);

	delete residuals_hist;
	delete c1;
	delete c2;

	delete rp;

	return 0;
}