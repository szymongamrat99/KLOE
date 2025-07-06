#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TPaveText.h>

#include <neutral_mom.h>
#include <pi0_photon_pair.h>
#include <triple_gaus.h>
#include <interference.h>
#include <kloe_class.h>

#include "../inc/omegarec.hpp"

int plots(int first_file, int last_file, int loopcount, int M, int range, Controls::DataType data_type)
{
	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

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
			mctruth_name = SystemPath::gen_vars_dir + SystemPath::root_files_dir + mctruth_filename + first_file + "_" + last_file + ext_root,
			omega_name = std::string(properties["variables"]["tree"]["filename"]["omegarec"]),
			tree_name = std::string(properties["variables"]["tree"]["treename"]["omegarec"]),
			filename_triangle = SystemPath::neutrec_dir + SystemPath::root_files_dir + neu_triangle_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

	file_mctruth = new TFile(mctruth_name);
	tree_mctruth = (TTree *)file_mctruth->Get(gen_vars_tree);

	tree_mctruth->SetBranchAddress("mctruth", &baseKin.mctruth_int);

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
	tree_omega->SetBranchAddress("anglePi0KaonCM", &anglePi0KaonCM),
	tree_omega->SetBranchAddress("anglePichKaonCM", &anglePichKaonCM),
	tree_omega->SetBranchAddress("anglePi0OmegaPhiCM", &anglePi0OmegaPhiCM),
	tree_omega->SetBranchAddress("anglePhiOmega", &anglePhiOmega);

	Float_t chi2min;

	tree_omega->SetBranchAddress("chi2min", &chi2min);

	// file_triangle = new TFile(filename_triangle);
	// tree_triangle = (TTree *)file_triangle->Get(neutrec_triangle_tree);

	Float_t
			Knetri_kinfit[10],
			minv4gam,
			ip_kinfit_triangle[3],
			chi2min_triangle,
			gamma_kinfit_triangle[4][8];

	Int_t
			done_kinfit_triangle;

	// tree_triangle->SetBranchAddress("fourgamma1triangle", gamma_kinfit_triangle[0]);
	// tree_triangle->SetBranchAddress("fourgamma2triangle", gamma_kinfit_triangle[1]);
	// tree_triangle->SetBranchAddress("fourgamma3triangle", gamma_kinfit_triangle[2]);
	// tree_triangle->SetBranchAddress("fourgamma4triangle", gamma_kinfit_triangle[3]);

	// tree_triangle->SetBranchAddress("fourKnetriangle", Knetri_kinfit);

	// tree_triangle->SetBranchAddress("chi2min", &chi2min_triangle);

	// tree_triangle->SetBranchAddress("minv4gam", &minv4gam);

	// tree_triangle->SetBranchAddress("iptriangle", ip_kinfit_triangle);
	// tree_triangle->SetBranchAddress("done_triangle", &done_kinfit_triangle);

	Float_t
			Kchrec[9],
			ip_avg[3],
			Kchboost[9];

	chain->SetBranchAddress("nclu", &baseKin.nclu);
	chain->SetBranchAddress("Xcl", baseKin.cluster[0]);
	chain->SetBranchAddress("Ycl", baseKin.cluster[1]);
	chain->SetBranchAddress("Zcl", baseKin.cluster[2]);
	chain->SetBranchAddress("Tcl", baseKin.cluster[3]);
	chain->SetBranchAddress("Enecl", baseKin.cluster[4]);

	chain->SetBranchAddress("ncll", baseKin.ncll);

	chain->SetBranchAddress("Kchrec", Kchrec);
	chain->SetBranchAddress("Kchboost", Kchboost);
	chain->SetBranchAddress("Bx", &ip_avg[0]);
	chain->SetBranchAddress("By", &ip_avg[1]);
	chain->SetBranchAddress("Bz", &ip_avg[2]);
	chain->SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
	chain->SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
	chain->SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
	chain->SetBranchAddress("Broots", &baseKin.phi_mom[3]);

	chain->SetBranchAddress("Dtmc", &baseKin.Dtmc);

	chain->SetBranchAddress("mctruth", &baseKin.mctruth);
	chain->SetBranchAddress("mcflag", &baseKin.mcflag);

	////////////////////////////////////////////////////////

	chain->AddFriend(tree_mctruth);
	chain->AddFriend(tree_omega);
	//chain->AddFriend(tree_triangle);

	/////////////////////////////////////////////////////////////////////

	const Int_t max_count = TMath::Factorial(4);
	Int_t count = 0, ind_gam[4], mc_ind[4] = {0, 1, 2, 3}, min_ind[max_count];
	Float_t clus_diff[max_count], clus_diff_min, neu_vtx[3], minv_pi0[2];

	Bool_t clus_time[max_count];

	Int_t nentries = (Int_t)chain->GetEntries();

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
				hist_control_chann[i].push_back(new TH1D(hist_name, "", 50.0, 0.0, 70.0));
				hist_name = "hist_pm_IP_" + std::to_string(i) + std::to_string(j);
				hist_pm_IP[i].push_back(new TH1D(hist_name, "", 50.0, 0.0, 70.0));
				hist_name = "hist_00_IP_" + std::to_string(i) + std::to_string(j);
				hist_00_IP[i].push_back(new TH1D(hist_name, "", 50.0, 0.0, 70.0));
			}
			else if (j == 2)
				hist[i].push_back(new TH1D(hist_name, "", 200.0, 500.0, 1000.0));
			else if (j == 3)
				hist[i].push_back(new TH1D(hist_name, "", 100.0, 0.0, 100.0));
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
		hist2d.push_back(new TH2D(hist2d_name, "", 50.0, 0.0, 100.0, 30.0, 0.0, 5.0));
	};

	for (Int_t i = 0; i < channNum; i++)
	{
		hist2d_name = channName[i] + "100Ip";
		hist2d_00IP_pmIP[0].push_back(new TH2D(hist2d_name, "", 50.0, 0.0, 5.0, 50.0, 0.0, 5.0));
		hist2d_name = channName[i] + "200IP";
		hist2d_00IP_pmIP[1].push_back(new TH2D(hist2d_name, "", 100.0, -70.0, 70.0, 100.0, -70.0, 70.0));
	};

	std::vector<TCanvas *> canvas_cont;
	TString canvas_cont_name = "";

	for (Int_t i = 0; i < 3; i++)
	{
		canvas_cont_name = "canvas_cont_" + std::to_string(i);
		canvas_cont.push_back(new TCanvas(canvas_cont_name, canvas_cont_name, 790, 790));
	};

	std::vector<TH1 *> hist_control;
	TString hist_control_name = "";

	for (Int_t i = 0; i < 3; i++)
	{
		hist_control_name = "hist_control_" + i;
		hist_control.push_back(new TH1D(hist_control_name, "", 20.0, -10.0, 10.0));
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

	// Initialization of interference function
	KLOE::interference event("split", 0, 91, -90, 90, split);

	TString
			name = SystemPath::efficiency_dir + SystemPath::root_files_dir + cut_vars_filename + first_file + "_" + last_file + ext_root;

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		for (Int_t j = 0; j < 4; j++)
			for (Int_t k = 0; k < 4; k++)
				PhotonMom[j][k] = gamma_kinfit_triangle[j][k];

		// Section for the histogram filling

		if (baseKin.mctruth_int == 4 && doneOmega == 1)
		{
			hist_control[0]->Fill(Kchrec[6] - neu_vtx_avg[0]);
			hist_control[1]->Fill(Kchrec[7] - neu_vtx_avg[1]);
			hist_control[2]->Fill(Kchrec[8] - neu_vtx_avg[2]);
		}

		Int_t
				sigmas_neu = 1.,
				sigmas_ip = 3.;

		Bool_t
				cond[6],
				cond_tot;

		Float_t
				resCh[3],
				resNeu[3],
				resComb[3];

		for(Int_t k = 0; k < 3; k++)
		{
			resCh[k] = properties["variables"]["Resolutions"]["vtxCharged"][k];
			resNeu[k] = properties["variables"]["Resolutions"]["vtxNeutral"]["triTriangle"][k];

			resComb[k] = resCh[k] + resNeu[k];
		}

		cond[0] = abs(Kchrec[6] - neu_vtx_avg[0]) < sigmas_neu * resComb[0];
		cond[1] = abs(Kchrec[7] - neu_vtx_avg[1]) < sigmas_neu * resComb[1];
		cond[2] = abs(Kchrec[8] - neu_vtx_avg[2]) < sigmas_neu * resComb[2];
		cond[3] = abs(Kchrec[6] - ip_avg[0]) < sigmas_ip * 2.63;
		cond[4] = abs(Kchrec[7] - ip_avg[1]) < sigmas_ip * 2.03;
		cond[5] = abs(Kchrec[8] - ip_avg[2]) < sigmas_ip * 3.39;

		cond_tot = cond[0] && cond[1] && cond[2];

		if (cond_tot && doneOmega == 1)
		{

			if (baseKin.mctruth_int == 1)
			{
				hist[0][0]->Fill(anglePi0OmegaPhiCM);
				hist[0][1]->Fill(Kchrec[5]);
				hist[0][2]->Fill(Omegarec[5]);
				hist[0][3]->Fill(chi2min);
				hist[0][4]->Fill(anglePhiOmega);
				hist[0][5]->Fill(anglePi0KaonCM);
				hist[0][6]->Fill(anglePichKaonCM);

				hist_control_chann[0][0]->Fill(rho_00);
				hist_control_chann[0][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[0][0]->Fill(rho_pm_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist_pm_IP[0][1]->Fill(Kchrec[8] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));

				hist_00_IP[0][0]->Fill(rho_00_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist_00_IP[0][1]->Fill(neu_vtx_avg[2] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));

				hist2d[0]->Fill(chi2min, rho);

				hist2d_00IP_pmIP[0][0]->Fill(rho_pm_IP, rho_00_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist2d_00IP_pmIP[1][0]->Fill(abs(Kchrec[8] - ip_avg[8]), abs(neu_vtx_avg[2] - ip_avg[8]), event.interf_function(baseKin.Dtmc, 0, par));
			}
			if (baseKin.mctruth_int == 3)
			{
				hist[1][0]->Fill(anglePi0OmegaPhiCM);
				hist[1][1]->Fill(Kchrec[5]);
				hist[1][2]->Fill(Omegarec[5]);
				hist[1][3]->Fill(chi2min);
				hist[1][4]->Fill(anglePhiOmega);
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

				hist2d[1]->Fill(chi2min, rho);
			}
			if (baseKin.mctruth_int == 4)
			{
				hist[2][0]->Fill(anglePi0OmegaPhiCM);
				hist[2][1]->Fill(Kchrec[5]);
				hist[2][2]->Fill(Omegarec[5]);
				hist[2][3]->Fill(chi2min);
				hist[2][4]->Fill(anglePhiOmega);
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

				hist2d[2]->Fill(chi2min, rho);
			}
			if (baseKin.mctruth_int == 5)
			{
				hist[3][0]->Fill(anglePi0OmegaPhiCM);
				hist[3][1]->Fill(Kchrec[5]);
				hist[3][2]->Fill(Omegarec[5]);
				hist[3][3]->Fill(chi2min);
				hist[3][4]->Fill(anglePhiOmega);
				hist[3][5]->Fill(anglePi0KaonCM);
				hist[3][6]->Fill(anglePichKaonCM);

				hist_control_chann[3][0]->Fill(rho_00);
				hist_control_chann[3][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[3][0]->Fill(rho_pm_IP);
				hist_pm_IP[3][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[3][0]->Fill(rho_00_IP);
				hist_00_IP[3][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d[3]->Fill(chi2min, rho);

				// hist2d_00IP_pmIP[0][0]->Fill(rho_pm_IP, rho_00_IP);

				hist2d_00IP_pmIP[0][3]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][3]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);
			}
			if (baseKin.mctruth_int == 6)
			{
				hist[4][0]->Fill(anglePi0OmegaPhiCM);
				hist[4][1]->Fill(Kchrec[5]);
				hist[4][2]->Fill(Omegarec[5]);
				hist[4][3]->Fill(chi2min);
				hist[4][4]->Fill(anglePhiOmega);
				hist[4][5]->Fill(anglePi0KaonCM);
				hist[4][6]->Fill(anglePichKaonCM);

				hist_control_chann[4][0]->Fill(rho_00);
				hist_control_chann[4][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[4][0]->Fill(rho_pm_IP);
				hist_pm_IP[4][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[4][0]->Fill(rho_00_IP);
				hist_00_IP[4][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d[4]->Fill(chi2min, rho);

				hist2d_00IP_pmIP[0][4]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][4]->Fill(Kchrec[8] - ip_avg[8], neu_vtx_avg[2] - ip_avg[8]);
			}
			if (baseKin.mctruth_int == 7)
			{
				hist[5][0]->Fill(anglePi0OmegaPhiCM);
				hist[5][1]->Fill(Kchrec[5]);
				hist[5][2]->Fill(Omegarec[5]);
				hist[5][3]->Fill(chi2min);
				hist[5][4]->Fill(anglePhiOmega);
				hist[5][5]->Fill(anglePi0KaonCM);
				hist[5][6]->Fill(anglePichKaonCM);

				hist_control_chann[5][0]->Fill(rho_00);
				hist_control_chann[5][1]->Fill(Kchrec[8] - neu_vtx_avg[2]);

				hist_pm_IP[5][0]->Fill(rho_pm_IP);
				hist_pm_IP[5][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[5][0]->Fill(rho_00_IP);
				hist_00_IP[5][1]->Fill(neu_vtx_avg[2] - ip_avg[8]);

				hist2d[5]->Fill(chi2min, rho);

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

	TString xTitle[7] = {"#angle(#pi^{0},#omega) [#circ]", "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]", "m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]", "#chi^{2}_{#omega#pi^{0}} [-]", "#angle(#phi,#omega) [#circ]", "#angle(#pi^{0},#pi^{0}) [#circ]", "#angle(#pi^{+},#pi^{-}) [#circ]"};

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
			// legend_chi2->AddEntry(hist[j][i], channName[j], "l");

			if (j == 0)
			{
				hist[j][i]->SetStats(0);
				hist[j][i]->GetXaxis()->SetTitle(xTitle[i]);
				hist[j][i]->GetYaxis()->SetTitle(yTitle);
				hist[j][i]->Draw("HIST");
			}
			else
				hist[j][i]->Draw("SAME");
		}
		// if (i < 3)
		// 	legend_chann->Draw();
		// else if (i == 3)
		// 	legend_chi2->Draw();

		canva[i]->Print(SystemPath::img_dir + "OmegaRec/hist" + std::to_string(i) + ext_img);
		// legend_chann->Clear();
		// legend_chi2->Clear();
	}

	for (Int_t j = 0; j < channNum; j++)
	{
		canva2d[j]->cd();
		canva2d[j]->SetLogz(1);

		hist2d[j]->SetMarkerColor(channColor[j]);
		hist2d[j]->SetMarkerSize(2);

		hist2d[j]->GetXaxis()->SetTitle("#chi^{2}_{#omega#pi^{0}} [-]");
		hist2d[j]->GetYaxis()->SetTitle("#sqrt{#rho_{+-}^{2} + #rho_{00}^{2}} [cm]");
		hist2d[j]->Draw("COLZ");

		canva2d[j]->Print(SystemPath::img_dir + "OmegaRec/angle_2d_kaon_cm_" + channName[j] + ext_img);
	}

	TString xTitleCont[3] = {"x_{#pi^{+}#pi^{-}} - x_{#pi^{0}#pi^{0}}", "y_{#pi^{+}#pi^{-}} - y_{#pi^{0}#pi^{0}}", "z_{#pi^{+}#pi^{-}} - z_{#pi^{0}#pi^{0}}"};

	TF1 *triple_fit = new TF1("triple_gaus", triple_gaus, -400.0, 400.0, 9, 1);
	triple_fit->SetParNames("Norm1", "Avg1", "Std1", "Norm2", "Avg2", "Std2", "Norm3", "Avg3", "Std3");
	triple_fit->SetLineWidth(4);

	TFitResultPtr result;

	Float_t
			parameter[3];

	TString
			fit_stats[3];

	TPaveText *fit_text = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");

	for (Int_t i = 0; i < 3; i++)
	{
		parameter[0] = hist_control[i]->GetEntries();
		parameter[1] = hist_control[i]->GetBinCenter(hist_control[i]->GetMaximumBin());
		parameter[2] = hist_control[i]->GetStdDev();

		triple_fit->SetParameters(0.15 * parameter[0], parameter[1], parameter[2], 0.15 * parameter[0], parameter[1] - 5.0, parameter[2], parameter[0], parameter[1] + 5.0, parameter[2]);

		triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
		triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
		triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

		triple_fit->SetParLimits(1, parameter[1] - 10.0, parameter[1] + 10.0);
		triple_fit->SetParLimits(4, parameter[1] - 10.0, parameter[1] + 10.0);
		triple_fit->SetParLimits(7, parameter[1] - 10.0, parameter[1] + 10.0);

		triple_fit->SetParLimits(2, 0.5 * parameter[2], 10.0 * parameter[2]);
		triple_fit->SetParLimits(5, 0.5 * parameter[2], 10.0 * parameter[2]);
		triple_fit->SetParLimits(8, 0.5 * parameter[2], 10.0 * parameter[2]);

		result = hist_control[i]->Fit(triple_fit, "SF");

		if(result == 0)
		{

			fit_stats[1] = Form("Mean = %.2f#pm%.2f", comb_mean(result->GetParams(), result->GetErrors()), comb_mean_err(result->GetParams(), result->GetErrors()));
			fit_stats[2] = Form("Std Dev = %.2f#pm%.2f", comb_std_dev(result->GetParams(), result->GetErrors()), comb_std_dev_err(result->GetParams(), result->GetErrors()));
		}

		fit_text->AddText(fit_stats[1]);
		fit_text->AddText(fit_stats[2]);

		canvas_cont[i]->cd();
		canvas_cont[i]->SetLogy(1);

		hist_control[i]->SetStats(0);
		hist_control[i]->SetLineColor(kBlue);
		hist_control[i]->SetLineWidth(3);

		hist_control[i]->GetXaxis()->SetTitle(xTitleCont[i]);
		hist_control[i]->GetYaxis()->SetTitle(yTitle);
		hist_control[i]->Draw();
		fit_text->Draw();

		canvas_cont[i]->Print(SystemPath::img_dir + "OmegaRec/cont_plot_" + std::to_string(i) + ext_img);

		fit_text->Clear();
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

			if (j == 1)
			{
				hist_control_chann[j][i]->SetStats(0);
				hist_control_chann[j][i]->GetXaxis()->SetTitle(xTitleControl[i]);
				hist_control_chann[j][i]->GetYaxis()->SetTitle(yTitle);
				hist_control_chann[j][i]->Draw("HIST");
				hist_control_chann[0][i]->Draw("SAME");
			}
			else if (j != 1 && j != 0)
				hist_control_chann[j][i]->Draw("SAME");
		}
		if (i < 3)
			legend_chann_cont->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i + 6]->Print(SystemPath::img_dir + "OmegaRec/hist_control_chann" + std::to_string(i) + ext_img);
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

			if (j == 1)
			{
				hist_pm_IP[j][i]->SetStats(0);
				hist_pm_IP[j][i]->GetXaxis()->SetTitle(xTitlepmIP[i]);
				hist_pm_IP[j][i]->GetYaxis()->SetTitle(yTitle);
				hist_pm_IP[j][i]->Draw("HIST");
				hist_pm_IP[0][i]->Draw("SAME");
			}
			else if (j != 1 && j != 0)
				hist_pm_IP[j][i]->Draw("SAME");
		}
		if (i < 3)
			legend_chann_cont->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i + 6]->Print(SystemPath::img_dir + "OmegaRec/hist_pm_IP_" + std::to_string(i) + ext_img);
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

			if (j == 1)
			{
				hist_00_IP[j][i]->SetStats(0);
				hist_00_IP[j][i]->GetXaxis()->SetTitle(xTitle00IP[i]);
				hist_00_IP[j][i]->GetYaxis()->SetTitle(yTitle);
				hist_00_IP[j][i]->Draw("HIST");
				hist_00_IP[0][i]->Draw("SAME");
			}
			else if (j != 1 && j != 0)
				hist_00_IP[j][i]->Draw("SAME");
		}
		if (i < 3)
			legend_chann_cont->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i + 8]->Print(SystemPath::img_dir + "OmegaRec/hist_00_IP_" + std::to_string(i) + ext_img);
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

		canva2d[j + channNum]->Print(SystemPath::img_dir + "OmegaRec/rho_2d_" + channName[j] + ext_img);
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

		canva2d[j + 2 * channNum]->Print(SystemPath::img_dir + "OmegaRec/z_2d_" + channName[j] + ext_img);
	}

	return 0;
}
