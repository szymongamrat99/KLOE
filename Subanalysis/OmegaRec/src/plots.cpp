#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <TTreeReader.h>
#include <TVirtualFitter.h>
#include <TLine.h>

#include <neutral_mom.h>
#include <pi0_photon_pair.h>
#include <triple_gaus.h>
#include <interference.h>
#include <kloe_class.h>

#include "../inc/omegarec.hpp"

int plots(TChain &chain, Short_t &loopcount, Short_t &numOfConstraints, Short_t &jmin, Short_t &jmax, Controls::DataType &dataType, KLOE::pm00 &Obj, ErrorHandling::ErrorLogs &logger)
{
	BaseKinematics baseKin;

	TFile
			*file_mctruth,
			*file_omega,
			*file_triangle;

	TTree
			*tree_mctruth,
			*tree_omega,
			*tree_triangle;

	TString
			omega_name = std::string(properties["variables"]["tree"]["filename"]["omegarec"]),
			tree_name = std::string(properties["variables"]["tree"]["treename"]["omegarec"]),
			mctruth_name = std::string(properties["variables"]["tree"]["filename"]["mctruth"]),
			mctruth_tree_name = std::string(properties["variables"]["tree"]["treename"]["mctruth"]);

	file_mctruth = new TFile(mctruth_name);
	tree_mctruth = (TTree *)file_mctruth->Get(mctruth_tree_name);

	file_omega = new TFile(omega_name);
	tree_omega = (TTree *)file_omega->Get(tree_name);

	Int_t
			doneOmega,
			g4takenomega[4];

	Float_t
			gammaomega[4][8],
			Omegapi0[6],
			Pi0[6],
			PichFourMom[2][4],
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
			angleOmegaPolar,
			angleKchPolar,
			neu_vtx_avg[3],
			Kchrec[9];

	tree_omega->SetBranchAddress("gamma1omega", gammaomega[0]);
	tree_omega->SetBranchAddress("gamma2omega", gammaomega[1]);
	tree_omega->SetBranchAddress("gamma3omega", gammaomega[2]);
	tree_omega->SetBranchAddress("gamma4omega", gammaomega[3]);

	tree_omega->SetBranchAddress("pich1", PichFourMom[0]);
	tree_omega->SetBranchAddress("pich2", PichFourMom[1]);

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
	tree_omega->SetBranchAddress("angleOmegaPolar", &angleOmegaPolar);
	tree_omega->SetBranchAddress("angleKchPolar", &angleKchPolar);

	tree_mctruth->SetBranchAddress("mctruth", &baseKin.mctruth_int);

	Float_t chi2min;

	// tree_omega->SetBranchAddress("chi2min", &chi2min);

	Float_t
			Knetri_kinfit[10],
			minv4gam,
			ip_kinfit_triangle[3],
			chi2min_triangle,
			gamma_kinfit_triangle[4][8];

	Int_t
			done_kinfit_triangle;

	Float_t
			ip_avg[3],
			Kchboost[9];

	chain.SetBranchAddress("nclu", &baseKin.nclu);
	chain.SetBranchAddress("Xcl", baseKin.cluster[0]);
	chain.SetBranchAddress("Ycl", baseKin.cluster[1]);
	chain.SetBranchAddress("Zcl", baseKin.cluster[2]);
	chain.SetBranchAddress("Tcl", baseKin.cluster[3]);
	chain.SetBranchAddress("Enecl", baseKin.cluster[4]);

	chain.SetBranchAddress("ncll", baseKin.ncll);

	chain.SetBranchAddress("Kchrec", Kchrec);
	chain.SetBranchAddress("Bx", &ip_avg[0]);
	chain.SetBranchAddress("By", &ip_avg[1]);
	chain.SetBranchAddress("Bz", &ip_avg[2]);
	chain.SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
	chain.SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
	chain.SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
	chain.SetBranchAddress("Broots", &baseKin.phi_mom[3]);

	chain.SetBranchAddress("Dtmc", &baseKin.Dtmc);
	chain.SetBranchAddress("mcflag", &baseKin.mcflag);

	////////////////////////////////////////////////////////

	chain.AddFriend(tree_omega, "omegatree");
	chain.AddFriend(tree_mctruth, "mctruthtree");

	/////////////////////////////////////////////////////////////////////

	const Int_t max_count = TMath::Factorial(4);
	Int_t count = 0, ind_gam[4], mc_ind[4] = {0, 1, 2, 3}, min_ind[max_count];
	Float_t clus_diff[max_count], clus_diff_min, neu_vtx[3], minv_pi0[2];

	Bool_t clus_time[max_count];

	Int_t nentries = (Int_t)chain.GetEntries();

	std::vector<TCanvas *> canva;
	TString canva_name = "";

	for (Int_t i = 0; i < 100; i++)
	{
		canva_name = "canva_" + std::to_string(i);
		canva.push_back(new TCanvas(canva_name, canva_name, 790, 790));
	};

	std::vector<TH1 *>
			hist[channNum],
			hist_control_chann[channNum],
			hist_pm_IP[channNum],
			hist_00_IP[channNum];
	TString hist_name = "";

	for (Int_t i = 0; i < channNum; i++)
		for (Int_t j = 0; j < 7; j++)
		{
			hist_name = "hist_" + std::to_string(i) + std::to_string(j);

			if (j < 2)
			{
				if (j == 0)
					hist[i].push_back(new TH1D(hist_name, "", 100.0, 0.0, 180.0));
				else if (j == 1)
					hist[i].push_back(new TH1D(hist_name, "", 100.0, 495.0, 500.0));

				hist_name = "hist_control_chann_" + std::to_string(i) + std::to_string(j);
				hist_control_chann[i].push_back(new TH1D(hist_name, "", 50.0, 0.0, 50.0));
				hist_name = "hist_pm_IP_" + std::to_string(i) + std::to_string(j);
				hist_pm_IP[i].push_back(new TH1D(hist_name, "", 50.0, 0.0, 50.0));
				hist_name = "hist_00_IP_" + std::to_string(i) + std::to_string(j);
				hist_00_IP[i].push_back(new TH1D(hist_name, "", 50.0, 0.0, 50.0));
			}
			else if (j == 2)
				hist[i].push_back(new TH1D(hist_name, "", 200.0, 300.0, 1000.0));
			else if (j == 3)
				hist[i].push_back(new TH1D(hist_name, "", 100.0, 0.0, 300.0));
			else
				hist[i].push_back(new TH1D(hist_name, "", 100.0, 0.0, 180.0));
		};

	std::vector<TCanvas *> canva2d;
	TString canva2d_name = "";

	for (Int_t i = 0; i < 3 * channNum; i++)
	{
		canva2d_name = "canva2d_" + std::to_string(i);
		canva2d.push_back(new TCanvas(canva2d_name, canva2d_name, 790, 790));
	};

	std::vector<TH2 *>
			hist2d,
			hist2d_00IP_pmIP[2];
	TString hist2d_name = "";

	for (Int_t i = 0; i < channNum; i++)
	{
		hist2d_name = channName[i];
		hist2d.push_back(new TH2D(hist2d_name, "", 200.0, 0.0, 300.0, 200.0, 600.0, 1000.0));
	};

	for (Int_t i = 0; i < channNum; i++)
	{
		hist2d_name = channName[i] + "100Ip";
		hist2d_00IP_pmIP[0].push_back(new TH2D(hist2d_name, "", 50.0, 0.0, 20.0, 50.0, 0.0, 20.0));
		hist2d_name = channName[i] + "200IP";
		hist2d_00IP_pmIP[1].push_back(new TH2D(hist2d_name, "", 100.0, -70.0, 70.0, 100.0, -70.0, 70.0));
	};

	std::vector<TCanvas *> canvas_cont;
	TString canvas_cont_name = "";

	for (Int_t i = 0; i < 10; i++)
	{
		canvas_cont_name = "canvas_cont_" + std::to_string(i);
		canvas_cont.push_back(new TCanvas(canvas_cont_name, canvas_cont_name, 790, 790));
	};

	std::vector<TH1 *> hist_control;
	TString hist_control_name = "";

	for (Int_t i = 0; i < 10; i++)
	{
		hist_control_name = "hist_control_" + std::to_string(i);

		if (i < 5)
			hist_control.push_back(new TH1D(hist_control_name, "", 20.0, -2.0, 2.0));
		else if (i == 5)
			hist_control.push_back(new TH1D(hist_control_name, "", 20.0, -5.0, 5.0));
		else if (i == 6)
			hist_control.push_back(new TH1D(hist_control_name, "", 100.0, 500.0, 1000.0));
		else if (i == 7)
			hist_control.push_back(new TH1D(hist_control_name, "", 30.0, 50.0, 250.0));
		else if (i >= 8)
			hist_control.push_back(new TH1D(hist_control_name, "", 30.0, 0.0, 3.0));
	};

	Double_t
			split[3] = {0., 0., 0.},
			par[2] = {Re, Im_nonCPT},
			M_omega_tmp[2] = {0.},
			M_omega_diff[2] = {0.};

	Int_t
			g4takenPi0[2][2];

	Float_t
			PhotonMom[4][8],
			PhotonMomPi0[2][2][4],
			Pi0Mom[2][4],
			OmegaMom[4];

	// Filling the histograms to get the limits for omega-pi0 decay
	chain.Draw("(omegatree.NeuVtxAvg[0] - Bx)>>hist_control_0", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("(omegatree.NeuVtxAvg[1] - By)>>hist_control_1", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("(omegatree.NeuVtxAvg[2] - Kchrec[8])>>hist_control_2", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("(Kchrec[6] - Bx)>>hist_control_3", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("(Kchrec[7] - By)>>hist_control_4", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("(Kchrec[8] - Bz)>>hist_control_5", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("omegatree.omega[5]>>hist_control_6", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("omegatree.omega[6]>>hist_control_7", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("sqrt(pow(omegatree.NeuVtxAvg[0] - Bx, 2) + pow(omegatree.NeuVtxAvg[1] - By, 2))>>hist_control_8", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");
	chain.Draw("sqrt(pow(Kchrec[6] - Bx, 2) + pow(Kchrec[7] - By, 2))>>hist_control_9", "omegatree.doneomega == 1 && mctruthtree.mctruth == 4");

	std::vector<TString> xTitleCont = {
			"x_{#pi^{0}#pi^{0}} - x_{IP} [cm]",
			"y_{#pi^{0}#pi^{0}} - y_{IP} [cm]",
			"z_{#pi^{0}#pi^{0}} - z_{#pi^{+}#pi^{-}} [cm]",
			"x_{#pi^{+}#pi^{-}} - x_{IP} [cm]",
			"y_{#pi^{+}#pi^{-}} - y_{IP} [cm]",
			"z_{#pi^{+}#pi^{-}} - z_{IP} [cm]",
			"m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]",
			"T_{#pi^{0}_{#omega}} [MeV]",
			"#sqrt{(x_{#pi^{0}#pi^{0}} - x_{IP})^{2} +(y_{#pi^{0}#pi^{0}} - y_{IP})^{2}} [cm]",
			"#sqrt{(x_{#pi^{+}#pi^{-}} - x_{IP})^{2} +(y_{#pi^{+}#pi^{-}} - y_{IP})^{2}} [cm]"};

	Float_t
			parameter[3];

	TString
			fit_stats[3];

	TPaveText *fit_text = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");

	Double_t
			meanInvMass = 0,
			meanInvMassErr = 0,
			stdInvMass = 0,
			stdInvMassErr = 0,
			meanKinEne = 0,
			meanKinEneErr = 0,
			stdKinEne = 0,
			stdKinEneErr = 0;

	std::vector<Double_t>
			stdDevOmegaVtx(6),
			stdDevOmegaVtxErr(6);

	for (Int_t i = 0; i < 10; i++)
	{
		TF1 *triple_fit;

		if (i < 8)
		{

			if (i < 5)
				triple_fit = new TF1("triple_gaus", triple_gaus, -2.0, 2.0, 9, 1);
			else if (i == 5)
				triple_fit = new TF1("triple_gaus", triple_gaus, -5.0, 5.0, 9, 1);
			else if (i == 6)
				triple_fit = new TF1("triple_gaus", triple_gaus, 500.0, 1000.0, 9, 1);
			else if (i == 7)
				triple_fit = new TF1("triple_gaus", triple_gaus, 80.0, 240.0, 9, 1);
			else
				triple_fit = new TF1("triple_gaus", triple_gaus, 80.0, 240.0, 9, 1);

			triple_fit->SetParNames("Norm1", "Avg1", "Std1", "Norm2", "Avg2", "Std2", "Norm3", "Avg3", "Std3");
			triple_fit->SetLineWidth(4);

			TFitResultPtr result;

			parameter[0] = hist_control[i]->GetEntries();
			parameter[1] = hist_control[i]->GetBinCenter(hist_control[i]->GetMaximumBin());
			parameter[2] = hist_control[i]->GetStdDev();

			triple_fit->SetParameters(parameter[0], parameter[1], parameter[2], 0.15 * parameter[0], parameter[1] - 1.5, parameter[2], parameter[0], parameter[1] + 1.0, parameter[2]);

			if (i >= 6 || i < 3)
			{
				triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
				triple_fit->FixParameter(3, 0.0);
				triple_fit->FixParameter(6, 0.0);

				triple_fit->SetParLimits(1, parameter[1] - 2.0, parameter[1] + 2.0);
				triple_fit->FixParameter(4, 0.0);
				triple_fit->FixParameter(7, 0.0);

				triple_fit->SetParLimits(2, 0.005 * parameter[2], 10.0 * parameter[2]);
				triple_fit->FixParameter(5, 1.0);
				triple_fit->FixParameter(8, 1.0);
			}
			else if (i >= 3 && i < 6)
			{
				triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
				triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
				triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

				triple_fit->SetParLimits(1, parameter[1] - 1.0, parameter[1] - 0.5);
				triple_fit->SetParLimits(4, parameter[1] - 0.5, parameter[1] + 0.5);
				triple_fit->SetParLimits(7, parameter[1] + 0.5, parameter[1] + 1.0);

				triple_fit->SetParLimits(2, 0.005 * parameter[2], 10.0 * parameter[2]);
				triple_fit->SetParLimits(5, 0.005 * parameter[2], 10.0 * parameter[2]);
				triple_fit->SetParLimits(8, 0.005 * parameter[2], 10.0 * parameter[2]);
			}

			TVirtualFitter::SetMaxIterations(1000000);

			result = hist_control[i]->Fit(triple_fit, "SL", "");

			if (result == 0)
			{
				if (i < 6)
				{
					stdDevOmegaVtx[i] = comb_std_dev(result->GetParams(), result->GetErrors());
					stdDevOmegaVtxErr[i] = comb_std_dev_err(result->GetParams(), result->GetErrors());
				}
				else if (i == 6)
				{
					meanInvMass = comb_mean(result->GetParams(), result->GetErrors());
					meanInvMassErr = comb_mean_err(result->GetParams(), result->GetErrors());
					stdInvMass = comb_std_dev(result->GetParams(), result->GetErrors());
					stdInvMassErr = comb_std_dev_err(result->GetParams(), result->GetErrors());
				}
				else if (i == 7)
				{
					meanKinEne = comb_mean(result->GetParams(), result->GetErrors());
					meanKinEneErr = comb_mean_err(result->GetParams(), result->GetErrors());
					stdKinEne = comb_std_dev(result->GetParams(), result->GetErrors());
					stdKinEneErr = comb_std_dev_err(result->GetParams(), result->GetErrors());
				}

				fit_stats[1] = Form("Mean = %.2f#pm%.2f", comb_mean(result->GetParams(), result->GetErrors()), comb_mean_err(result->GetParams(), result->GetErrors()));
				fit_stats[2] = Form("Std Dev = %.2f#pm%.2f", comb_std_dev(result->GetParams(), result->GetErrors()), comb_std_dev_err(result->GetParams(), result->GetErrors()));

				fit_text->AddText(fit_stats[1]);
				fit_text->AddText(fit_stats[2]);
			}
			else
			{
				stdDevOmegaVtx[i] = 0;
				stdDevOmegaVtxErr[i] = 0;
			}
		}

		canvas_cont[i]->cd();
		canvas_cont[i]->SetLogy(0);

		hist_control[i]->SetStats(0);
		hist_control[i]->SetLineColor(kBlue);
		hist_control[i]->SetLineWidth(3);

		hist_control[i]->GetXaxis()->SetTitle(xTitleCont[i]);
		hist_control[i]->GetYaxis()->SetTitle("Counts");
		hist_control[i]->GetYaxis()->SetRangeUser(0.0, 1.2 * hist_control[i]->GetMaximum());
		hist_control[i]->Draw();

		if(i < 8)
			fit_text->Draw();

		canvas_cont[i]->Print(img_dir + "OmegaRec/cont_plot_" + std::to_string(i) + ext_img);

		if(i < 8)
			fit_text->Clear();

		if(i < 8)
			delete triple_fit;
	}

	// -------------------------------------------------------------------------------------------------

	TH2 *energyVsMass = new TH2D("EnergyVsMass", ";T_{#pi^{0}_{#omega}} [MeV];m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]", 100, 0, 300, 100, 600, 1000);

	// Initialization of interference function
	KLOE::interference event("split", 0, 91, -90, 90, split);

	for (Int_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		for (Int_t j = 0; j < 4; j++)
			for (Int_t k = 0; k < 4; k++)
				PhotonMom[j][k] = gamma_kinfit_triangle[j][k];

		Int_t
				sigmas_neu = 1.,
				sigmas_ip = 3;

		Bool_t
				cond[6],
				cond_tot;

		cond[0] = rho_pm_IP < 0.5 * sigmas_ip;
		cond[1] = rho_00_IP < 0.6 * sigmas_ip;
		cond[2] = abs(Kchrec[8] - ip_avg[2]) < sigmas_ip * 1.06;
		cond[3] = abs(neu_vtx_avg[2] - Kchrec[8]) < sigmas_ip * 0.39;

		cond_tot = 1; // cond[0] && cond[1] && cond[2] && cond[3];

		if (cond_tot && doneOmega == 1)
		{
			if (baseKin.mctruth_int == 1)
			{
				hist[0][0]->Fill(angleKchPolar);
				hist[0][1]->Fill(Kchrec[5]);
				hist[0][2]->Fill(Omegarec[5]);
				hist[0][3]->Fill(Omegapi0[3] - Omegapi0[5]);
				hist[0][4]->Fill(angleOmegaPolar);
				hist[0][5]->Fill(anglePi0KaonCM);
				hist[0][6]->Fill(anglePichKaonCM);

				hist_control_chann[0][0]->Fill(rho_00);
				hist_control_chann[0][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[0][0]->Fill(rho_pm_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist_pm_IP[0][1]->Fill(Kchrec[8] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));

				hist_00_IP[0][0]->Fill(rho_00_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist_00_IP[0][1]->Fill(neu_vtx_avg[2] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));

				hist2d[0]->Fill(Omegapi0[3] - Omegapi0[5], Omegarec[5]);

				hist2d_00IP_pmIP[0][0]->Fill(rho_pm_IP, rho_00_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist2d_00IP_pmIP[1][0]->Fill(abs(Kchrec[8] - ip_avg[8]), abs(neu_vtx_avg[2] - ip_avg[8]), event.interf_function(baseKin.Dtmc, 0, par));
			}
			if (baseKin.mctruth_int == 3)
			{
				hist[1][0]->Fill(angleKchPolar);
				hist[1][1]->Fill(Kchrec[5]);
				hist[1][2]->Fill(Omegarec[5]);
				hist[1][3]->Fill(Omegapi0[3] - Omegapi0[5]);
				hist[1][4]->Fill(angleOmegaPolar);
				hist[1][5]->Fill(anglePi0KaonCM);
				hist[1][6]->Fill(anglePichKaonCM);

				hist_control_chann[1][0]->Fill(rho_00);
				hist_control_chann[1][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[1][0]->Fill(rho_pm_IP);
				hist_pm_IP[1][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[1][0]->Fill(rho_00_IP);
				hist_00_IP[1][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d_00IP_pmIP[0][1]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][1]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);

				hist2d[1]->Fill(Omegapi0[3] - Omegapi0[5], Omegarec[5]);
			}
			if (baseKin.mctruth_int == 4)
			{
				hist[2][0]->Fill(angleKchPolar);
				hist[2][1]->Fill(Kchrec[5]);
				hist[2][2]->Fill(Omegarec[5]);
				hist[2][3]->Fill(Omegapi0[3] - Omegapi0[5]);
				hist[2][4]->Fill(angleOmegaPolar);
				hist[2][5]->Fill(anglePi0KaonCM);
				hist[2][6]->Fill(anglePichKaonCM);

				hist_control_chann[2][0]->Fill(rho_00);
				hist_control_chann[2][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[2][0]->Fill(rho_pm_IP);
				hist_pm_IP[2][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[2][0]->Fill(rho_00_IP);
				hist_00_IP[2][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d_00IP_pmIP[0][2]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][2]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);

				hist2d[2]->Fill(Omegapi0[3] - Omegapi0[5], Omegarec[5]);

				// Fill histogram to get kin ene vs. mass width
				energyVsMass->Fill(Omegarec[6], Omegarec[5]);
			}
			if (baseKin.mctruth_int == 5)
			{
				hist[3][0]->Fill(angleKchPolar);
				hist[3][1]->Fill(Kchrec[5]);
				hist[3][2]->Fill(Omegarec[5]);
				hist[3][3]->Fill(Omegapi0[3] - Omegapi0[5]);
				hist[3][4]->Fill(angleOmegaPolar);
				hist[3][5]->Fill(anglePi0KaonCM);
				hist[3][6]->Fill(anglePichKaonCM);

				hist_control_chann[3][0]->Fill(rho_00);
				hist_control_chann[3][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[3][0]->Fill(rho_pm_IP);
				hist_pm_IP[3][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[3][0]->Fill(rho_00_IP);
				hist_00_IP[3][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d[3]->Fill(Omegapi0[3] - Omegapi0[5], Omegarec[5]);

				hist2d_00IP_pmIP[0][0]->Fill(rho_pm_IP, rho_00_IP);

				hist2d_00IP_pmIP[0][3]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][3]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);
			}
			if (baseKin.mctruth_int == 6)
			{
				hist[4][0]->Fill(angleKchPolar);
				hist[4][1]->Fill(Kchrec[5]);
				hist[4][2]->Fill(Omegarec[5]);
				hist[4][3]->Fill(Omegapi0[3] - Omegapi0[5]);
				hist[4][4]->Fill(angleOmegaPolar);
				hist[4][5]->Fill(anglePi0KaonCM);
				hist[4][6]->Fill(anglePichKaonCM);

				hist_control_chann[4][0]->Fill(rho_00);
				hist_control_chann[4][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[4][0]->Fill(rho_pm_IP);
				hist_pm_IP[4][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[4][0]->Fill(rho_00_IP);
				hist_00_IP[4][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d[4]->Fill(Omegapi0[3] - Omegapi0[5], Omegarec[5]);

				hist2d_00IP_pmIP[0][4]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][4]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);
			}
			if (baseKin.mctruth_int == 7)
			{
				hist[5][0]->Fill(angleKchPolar);
				hist[5][1]->Fill(Kchrec[5]);
				hist[5][2]->Fill(Omegarec[5]);
				hist[5][3]->Fill(Omegapi0[3] - Omegapi0[5]);
				hist[5][4]->Fill(angleOmegaPolar);
				hist[5][5]->Fill(anglePi0KaonCM);
				hist[5][6]->Fill(anglePichKaonCM);

				hist_control_chann[5][0]->Fill(rho_00);
				hist_control_chann[5][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[5][0]->Fill(rho_pm_IP);
				hist_pm_IP[5][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[5][0]->Fill(rho_00_IP);
				hist_00_IP[5][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d[5]->Fill(Omegapi0[3] - Omegapi0[5], Omegarec[5]);

				hist2d_00IP_pmIP[0][5]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][5]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);
			}
		}
	}

	hist_00_IP[0][0]->Scale(hist_00_IP[0][0]->GetEntries() / hist_00_IP[0][0]->Integral());
	hist_pm_IP[0][0]->Scale(hist_pm_IP[0][0]->GetEntries() / hist_pm_IP[0][0]->Integral());

	hist2d_00IP_pmIP[0][0]->Scale(hist2d_00IP_pmIP[0][0]->GetEntries() / hist2d_00IP_pmIP[0][0]->Integral());

	TLegend *legend_chann = new TLegend(0.1, 0.5, 0.4, 0.9);
	legend_chann->SetFillColor(kWhite);

	TLegend *legend_chi2 = new TLegend(0.6, 0.5, 0.9, 0.9);
	legend_chi2->SetFillColor(kWhite);

	TString xTitle[7] = {"#angle(#vec{p}_{K#rightarrow#pi^{+}#pi^{-}}, #hat{z}) [#circ]", "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]", "m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]", "T_{0} [MeV]", "#angle(#vec{p}_{#omega}, #hat{z}) [#circ]", "#angle(#pi^{0},#pi^{0}) [#circ]", "#angle(#pi^{+},#pi^{-}) [#circ]"};

	TString yTitle = "Counts";

	for (Int_t i = 0; i < 7; i++)
	{
		canva[i]->cd();
		canva[i]->SetLogy(1);
		for (Int_t j = 0; j < channNum; j++)
		{
			hist[j][i]->SetLineColor(channColor[j]);
			hist[j][i]->SetLineWidth(3);

			// legend_chann->AddEntry(hist[j][i], channName[j], "l");
			legend_chi2->AddEntry(hist[j][i], channName[j], "l");

			if (j == 0)
			{
				hist[j][i]->SetStats(0);
				hist[j][i]->GetXaxis()->SetTitle(xTitle[i]);
				hist[j][i]->GetYaxis()->SetTitle(yTitle);
				hist[j][i]->GetYaxis()->SetRangeUser(0.1, 10.0 * hist[3][i]->GetMaximum());
				hist[j][i]->Draw("HIST");
			}
			else
				hist[j][i]->Draw("SAME");
		}
		// if (i < 3)
		// 	legend_chann->Draw();
		if (i == 3)
			legend_chi2->Draw();

		canva[i]->Print(img_dir + "OmegaRec/hist" + std::to_string(i) + ext_img);
		// legend_chann->Clear();
		legend_chi2->Clear();
	}

	for (Int_t j = 0; j < channNum; j++)
	{
		canva2d[j]->cd();
		canva2d[j]->SetLogz(1);

		hist2d[j]->SetMarkerColor(channColor[j]);
		hist2d[j]->SetMarkerSize(2);

		hist2d[j]->GetXaxis()->SetTitle("T_{0} [MeV]");
		hist2d[j]->GetYaxis()->SetTitle("m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]");
		hist2d[j]->Draw("COLZ");

		canva2d[j]->Print(img_dir + "OmegaRec/angle_2d_kaon_" + channName[j] + ext_img);
	}

	TString xTitleControl[2] = {"#rho_{+-,00} [cm]", "z_{+-,00} [cm]"};

	TLegend *legend_chann_cont = new TLegend(0.6, 0.5, 0.9, 0.9);
	legend_chann->SetFillColor(kWhite);

	for (Int_t i = 0; i < 2; i++)
	{
		canva[i + 6]->cd();
		canva[i + 6]->SetLogy(1);
		for (Int_t j = 0; j < channNum; j++)
		{
			hist_control_chann[j][i]->SetLineColor(channColor[j]);
			hist_control_chann[j][i]->SetLineWidth(3);

			legend_chann_cont->AddEntry(hist_control_chann[j][i], channName[j], "l");
			legend_chi2->AddEntry(hist_control_chann[j][i], channName[j], "l");

			if (j == 3)
			{
				hist_control_chann[j][i]->SetStats(0);
				hist_control_chann[j][i]->GetXaxis()->SetTitle(xTitleControl[i]);
				hist_control_chann[j][i]->GetYaxis()->SetTitle(yTitle);
				hist_control_chann[j][i]->Draw("HIST");
				hist_control_chann[0][i]->Draw("SAME");
			}
			else if (j != 3 && j != 0)
				hist_control_chann[j][i]->Draw("SAME");
		}
		if (i < 3)
			legend_chann_cont->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i + 6]->Print(img_dir + "OmegaRec/hist_control_chann" + std::to_string(i) + ext_img);
		legend_chann_cont->Clear();
		legend_chi2->Clear();
	}

	TString xTitlepmIP[2] = {"#rho_{+-,IP} [cm]", "z_{+-,IP} [cm]"};

	for (Int_t i = 0; i < 2; i++)
	{
		canva[i + 6]->cd();
		canva[i + 6]->SetLogy(1);
		for (Int_t j = 0; j < channNum; j++)
		{
			hist_pm_IP[j][i]->SetLineColor(channColor[j]);
			hist_pm_IP[j][i]->SetLineWidth(3);

			legend_chann_cont->AddEntry(hist_pm_IP[j][i], channName[j], "l");
			legend_chi2->AddEntry(hist_pm_IP[j][i], channName[j], "l");

			if (j == 0)
			{
				hist_pm_IP[0][i]->SetStats(0);
				hist_pm_IP[0][i]->GetXaxis()->SetTitle(xTitlepmIP[i]);
				hist_pm_IP[j][i]->GetYaxis()->SetTitle(yTitle);
				hist_pm_IP[j][i]->GetYaxis()->SetRangeUser(1E-2, hist_pm_IP[3][i]->GetMaximum() * 100.);
				hist_pm_IP[j][i]->Draw("HIST");
			}
			else
				hist_pm_IP[j][i]->Draw("SAME");
		}
		if (i < 3)
			legend_chann_cont->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i + 6]->Print(img_dir + "OmegaRec/hist_pm_IP_" + std::to_string(i) + ext_img);
		legend_chann_cont->Clear();
		legend_chi2->Clear();
	}

	TString xTitle00IP[2] = {"#rho_{00,IP} [cm]", "z_{00,IP} [cm]"};

	for (Int_t i = 0; i < 2; i++)
	{
		canva[i + 8]->cd();
		canva[i + 8]->SetLogy(1);
		for (Int_t j = 0; j < channNum; j++)
		{
			hist_00_IP[j][i]->SetLineColor(channColor[j]);
			hist_00_IP[j][i]->SetLineWidth(3);

			legend_chann_cont->AddEntry(hist_00_IP[j][i], channName[j], "l");
			legend_chi2->AddEntry(hist_00_IP[j][i], channName[j], "l");

			if (j == 0)
			{
				hist_00_IP[j][i]->SetStats(0);
				hist_00_IP[j][i]->GetXaxis()->SetTitle(xTitle00IP[i]);
				hist_00_IP[j][i]->GetYaxis()->SetTitle(yTitle);
				hist_00_IP[j][i]->GetYaxis()->SetRangeUser(1E-2, hist_00_IP[3][i]->GetMaximum() * 100.);
				hist_00_IP[j][i]->Draw("HIST");
			}
			else
				hist_00_IP[j][i]->Draw("SAME");
		}
		if (i < 3)
			legend_chann_cont->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i + 8]->Print(img_dir + "OmegaRec/hist_00_IP_" + std::to_string(i) + ext_img);
		legend_chann_cont->Clear();
		legend_chi2->Clear();
	}

	for (Int_t j = 0; j < channNum; j++)
	{
		canva2d[j + channNum]->cd();
		canva2d[j + channNum]->SetLogz(1);

		hist2d_00IP_pmIP[0][j]->SetMarkerColor(channColor[j]);
		hist2d_00IP_pmIP[0][j]->SetMarkerSize(2);

		hist2d_00IP_pmIP[0][j]->GetXaxis()->SetTitle("#rho_{+-,IP} [cm]");
		hist2d_00IP_pmIP[0][j]->GetYaxis()->SetTitle("#rho_{00,IP} [cm]");
		hist2d_00IP_pmIP[0][j]->Draw("COLZ");

		canva2d[j + channNum]->Print(img_dir + "OmegaRec/rho_2d_" + channName[j] + ext_img);
	}

	for (Int_t j = 0; j < channNum; j++)
	{
		canva2d[j + 2 * channNum]->cd();
		canva2d[j + 2 * channNum]->SetLogz(1);

		hist2d_00IP_pmIP[1][j]->SetMarkerColor(channColor[j]);
		hist2d_00IP_pmIP[1][j]->SetMarkerSize(2);

		hist2d_00IP_pmIP[1][j]->GetXaxis()->SetTitle("z_{+-,IP} [cm]");
		hist2d_00IP_pmIP[1][j]->GetYaxis()->SetTitle("z_{00,IP} [cm]");
		hist2d_00IP_pmIP[1][j]->Draw("COLZ");

		canva2d[j + 2 * channNum]->Print(img_dir + "OmegaRec/z_2d_" + channName[j] + ext_img);
	}

	TF1 *fit2D = new TF1("fit2D", "[0]*x + [1]", 0, 300); // pol1: y = a*x + b

	energyVsMass->Fit(fit2D);

	Double_t a = fit2D->GetParameter(0);
	Double_t b = fit2D->GetParameter(1);

	Double_t theta = atan(a); // Get angle for rotation of points

	// Histogram of projection
	TH1 *hProj = new TH1D("hProj", ";Projection along fitted line [-];Counts", 50, 400, 800);
	TH2 *energyVsMassRot = new TH2D("Projection", "", 50, 400, 750, 50, 450, 650);

	for (int ix = 1; ix <= energyVsMass->GetNbinsX(); ix++)
	{
		for (int iy = 1; iy <= energyVsMass->GetNbinsY(); iy++)
		{
			int count = energyVsMass->GetBinContent(ix, iy);

			double x = energyVsMass->GetXaxis()->GetBinCenter(ix);
			double y = energyVsMass->GetYaxis()->GetBinCenter(iy);

			// Rotacja
			double x_rot = x * cos(-theta) - y * sin(-theta);
			double y_rot = x * sin(-theta) + y * cos(-theta);

			for (Int_t i = 0; i < count; i++)
				energyVsMassRot->Fill(x_rot, y_rot);
		}
	}

	hProj = energyVsMassRot->ProjectionY();

	// Dopasowanie Gaussa do szerokości
	TF1 *gaus = new TF1("gaus", "gaus", -100, 100);
	hProj->Fit(gaus);
	Double_t sigma = gaus->GetParameter(2);
	Double_t sigmaErr = gaus->GetParError(2);

	Double_t normFactor = 3 * sigma / sqrt(a * a + 1); // Norm of normal vector
	Double_t normalX = -a * normFactor;
	Double_t normalY = normFactor;

	TF1 *lineplus3sigma = new TF1("line_plus_3sigma", "[0]*x + [1]", 0, 300);		// pol1: y = a*x + b
	TF1 *lineminus3sigma = new TF1("line_minus_3sigma", "[0]*x + [1]", 0, 300); // pol1: y = a*x + b

	TLine *T0plus3sigma = new TLine(meanKinEne + stdKinEne, meanInvMass - stdInvMass, meanKinEne + stdKinEne, meanInvMass + stdInvMass);
	TLine *T0minus3sigma = new TLine(meanKinEne - stdKinEne, meanInvMass - stdInvMass, meanKinEne - stdKinEne, meanInvMass + stdInvMass);

	TLine *minvplus3sigma = new TLine(meanKinEne - stdKinEne, meanInvMass + stdInvMass, meanKinEne + stdKinEne, meanInvMass + stdInvMass);
	TLine *minvminus3sigma = new TLine(meanKinEne - stdKinEne, meanInvMass - stdInvMass, meanKinEne + stdKinEne, meanInvMass - stdInvMass);

	lineplus3sigma->SetParameter(0, a);
	lineminus3sigma->SetParameter(0, a);
	lineplus3sigma->SetParameter(1, b + (1 - pow(a, 2)) * normFactor);
	lineminus3sigma->SetParameter(1, b - (1 - pow(a, 2)) * normFactor);

	// Rysowanie
	TCanvas *c1 = new TCanvas("c1", "Rotacja TH2D", 1200, 600);
	c1->Divide(2, 1);
	c1->cd(1);

	gStyle->SetStatX(0.4); // prawa krawędź statboxa
	gStyle->SetStatY(0.9); // górna krawędź statboxa

	energyVsMass->Draw("COLZ");
	fit2D->SetLineColor(kBlack);
	lineplus3sigma->SetLineColor(kBlack);
	lineminus3sigma->SetLineColor(kBlack);
	fit2D->Draw("same");
	lineplus3sigma->Draw("same");
	lineminus3sigma->Draw("same");
	T0plus3sigma->Draw("same");
	T0minus3sigma->Draw("same");
	minvplus3sigma->Draw("same");
	minvminus3sigma->Draw("same");

	c1->cd(2);

	gStyle->SetStatX(0.4); // prawa krawędź statboxa
	gStyle->SetStatY(0.9); // górna krawędź statboxa

	hProj->GetXaxis()->SetTitle("Deviation from fitted line [-]");
	hProj->GetYaxis()->SetTitle("Counts");
	hProj->Draw();
	gaus->SetLineColor(kBlue);
	gaus->Draw("same");

	c1->Print(omegarec_dir + img_dir + "widthOfDistribution" + ext_img);

	// Addition of stdDevs to the properties file
	std::string decayType[2] = {"neutral", "charged"};
	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 3; j++)
		{
			properties["variables"]["OmegaRec"]["fiducialVolume"][decayType[i]]["stdDev"][j] = stdDevOmegaVtx[i * 3 + j];
			properties["variables"]["OmegaRec"]["fiducialVolume"][decayType[i]]["error"][j] = stdDevOmegaVtxErr[i * 3 + j];
		}

	properties["variables"]["OmegaRec"]["invMass"]["mean"]["value"] = meanInvMass;
	properties["variables"]["OmegaRec"]["invMass"]["mean"]["error"] = meanInvMassErr;
	properties["variables"]["OmegaRec"]["invMass"]["stdDev"]["value"] = stdInvMass;
	properties["variables"]["OmegaRec"]["invMass"]["stdDev"]["error"] = stdInvMassErr;

	properties["variables"]["OmegaRec"]["kinEne"]["mean"]["value"] = meanKinEne;
	properties["variables"]["OmegaRec"]["kinEne"]["mean"]["error"] = meanKinEneErr;
	properties["variables"]["OmegaRec"]["kinEne"]["stdDev"]["value"] = stdKinEne;
	properties["variables"]["OmegaRec"]["kinEne"]["stdDev"]["error"] = stdKinEneErr;

	properties["variables"]["OmegaRec"]["combined"]["stdDev"]["value"] = sigma;
	properties["variables"]["OmegaRec"]["combined"]["stdDev"]["error"] = sigmaErr;

	properties["variables"]["OmegaRec"]["combined"]["line"]["slope"] = a;
	properties["variables"]["OmegaRec"]["combined"]["line"]["inter"] = b;

	properties["lastScript"] = "Plots of Omega Reconstruction";
	properties["lastUpdate"] = Obj.getCurrentTimestamp();

	std::ofstream outfile(propName);
	outfile << properties.dump(4);
	outfile.close();

	return 0;
}
