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
#include <lorentz_transf.h>

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
			mctruth_name = gen_vars_dir + root_files_dir + mctruth_filename + first_file + "_" + last_file + ext_root,
			omega_name = omegarec_dir + root_files_dir + omega_rec_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root,
			filename_triangle = neutrec_dir + root_files_dir + neu_triangle_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

	file_mctruth = new TFile(mctruth_name);
	tree_mctruth = (TTree *)file_mctruth->Get(gen_vars_tree);

	tree_mctruth->SetBranchAddress("mctruth", &baseKin.mctruth_int);

	file_omega = new TFile(omega_name);
	tree_omega = (TTree *)file_omega->Get(omegarec_tree);

	Float_t
			gamma_kinfit[4][8],
			Omegarec_kinfit[10],
			pi0_kinfit[2][10],
			ip_kinfit[3],
			chi2min;

	Int_t
			done_kinfit,
			g4taken_kinfit[4],
			bunchnum;

	tree_omega->SetBranchAddress("gamma1tri_kinfit", gamma_kinfit[0]);
	tree_omega->SetBranchAddress("gamma2tri_kinfit", gamma_kinfit[1]);
	// tree_omega->SetBranchAddress("gamma3tri_kinfit", gamma_kinfit[2]);
	// tree_omega->SetBranchAddress("gamma4tri_kinfit", gamma_kinfit[3]);

	tree_omega->SetBranchAddress("omegapi0tri_kinfit", pi0_kinfit[0]);
	tree_omega->SetBranchAddress("pi0_kinfit", pi0_kinfit[1]);

	tree_omega->SetBranchAddress("omega_kinfit", Omegarec_kinfit);

	tree_omega->SetBranchAddress("iptri_kinfit", ip_kinfit);
	tree_omega->SetBranchAddress("done4_kinfit", &done_kinfit);

	tree_omega->SetBranchAddress("g4takentri_kinfit", g4taken_kinfit);

	tree_omega->SetBranchAddress("bunchnum", &bunchnum);

	tree_omega->SetBranchAddress("chi2min", &chi2min);

	file_triangle = new TFile(filename_triangle);
	tree_triangle = (TTree *)file_triangle->Get(neutrec_triangle_tree);

	Float_t
			Knetri_kinfit[10],
			minv4gam,
			ip_kinfit_triangle[3],
			chi2min_triangle,
			gamma_kinfit_triangle[4][8];

	Int_t
			done_kinfit_triangle;

	tree_triangle->SetBranchAddress("fourgamma1triangle", gamma_kinfit_triangle[0]);
	tree_triangle->SetBranchAddress("fourgamma2triangle", gamma_kinfit_triangle[1]);
	tree_triangle->SetBranchAddress("fourgamma3triangle", gamma_kinfit_triangle[2]);
	tree_triangle->SetBranchAddress("fourgamma4triangle", gamma_kinfit_triangle[3]);

	tree_triangle->SetBranchAddress("fourKnetriangle", Knetri_kinfit);

	tree_triangle->SetBranchAddress("g4taken_triangle", g4taken_kinfit);

	tree_triangle->SetBranchAddress("chi2min", &chi2min_triangle);

	tree_triangle->SetBranchAddress("minv4gam", &minv4gam);

	tree_triangle->SetBranchAddress("iptriangle", ip_kinfit_triangle);
	tree_triangle->SetBranchAddress("done_triangle", &done_kinfit_triangle);

	Float_t
			PichFourMom[2][4],
			Kchrec[9],
			ip_avg[3],
			Kchboost[9];

	chain->SetBranchAddress("nclu", &baseKin.nclu);
	chain->SetBranchAddress("Xacl", baseKin.cluster[0]);
	chain->SetBranchAddress("Yacl", baseKin.cluster[1]);
	chain->SetBranchAddress("Zacl", baseKin.cluster[2]);
	chain->SetBranchAddress("Tcl", baseKin.cluster[3]);
	chain->SetBranchAddress("Enecl", baseKin.cluster[4]);

	chain->SetBranchAddress("ncll", baseKin.ncll);

	chain->SetBranchAddress("trk1", PichFourMom[0]);
	chain->SetBranchAddress("trk2", PichFourMom[1]);

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
	chain->AddFriend(tree_triangle);

	/////////////////////////////////////////////////////////////////////

	const Int_t max_count = TMath::Factorial(4);
	Int_t count = 0, ind_gam[4], mc_ind[4] = {0, 1, 2, 3}, min_ind[max_count];
	Float_t clus_diff[max_count], clus_diff_min, neu_vtx[3], minv_pi0[2], minv_omega;

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
		for (Int_t j = 0; j < 4; j++)
		{
			hist_name = "hist_" + std::to_string(i) + std::to_string(j);

			if (j < 2)
			{
				if (j == 0)
					hist[i].push_back(new TH1D(hist_name, "", 100.0, 0.0, 200.0));
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
			else
				hist[i].push_back(new TH1D(hist_name, "", 100.0, 0.0, 200.0));
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
		hist2d.push_back(new TH2D(hist2d_name, "", 100.0, 178.0, 182.0, 100.0, 0.0, 10.0));
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
			PhotonMom[4][4],
			PhotonMomPi0[2][2][4],
			Pi0Mom[2][4],
			OmegaMom[4];

	// Initialization of interference function
	KLOE::interference event("split", 0, 91, -90, 90, split);

	TString
			name = efficiency_dir + root_files_dir + cut_vars_filename + first_file + "_" + last_file + ext_root;

	Float_t
			rho_00,
			rho_pm_IP,
			rho_00_IP,
			rho,
			lengthKaon[4],
			lengthCh,
			anglePi0KaonCM,
			anglePichKaonCM,
			anglePi0OmegaPhiCM,
			anglePhiOmega;

	// Initialization of the cut vars TTree
	TFile *file = new TFile(name, "update");
	TTree *tree = new TTree(cut_vars_tree, "Cut variables to use");

	TBranch
			*b1 = tree->Branch("rho_00", &rho_00, "rho_00/F"),
			*b2 = tree->Branch("rho_00_IP", &rho_00_IP, "rho_00_IP/F"),
			*b3 = tree->Branch("rho_pm_IP", &rho_pm_IP, "rho_pm_IP/F"),
			*b4 = tree->Branch("rho", &rho, "rho/F"),
			*b5 = tree->Branch("anglePi0KaonCM", &anglePi0KaonCM, "anglePi0KaonCM/F"),
			*b6 = tree->Branch("anglePichKaonCM", &anglePichKaonCM, "anglePichKaonCM/F"),
			*b7 = tree->Branch("anglePi0OmegaPhiCM", &anglePi0OmegaPhiCM, "anglePi0OmegaPhiCM/F"),
			*b8 = tree->Branch("anglePhiOmega", &anglePhiOmega, "anglePhiOmega/F");

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		for (Int_t j = 0; j < 4; j++)
			for (Int_t k = 0; k < 4; k++)
				PhotonMom[j][k] = gamma_kinfit_triangle[j][k];

		// Pairing of photons to 2Pi0 - without Omega assumption!
		Pi0PhotonPair(g4taken_kinfit, PhotonMom, g4takenPi0, PhotonMomPi0, Pi0Mom, false, PichFourMom, OmegaMom);

		// Calculation of the "kaon" path from IP.
		rho_00 = sqrt(pow(Kchrec[6] - Knetri_kinfit[6], 2) + pow(Kchrec[7] - Knetri_kinfit[7], 2));
		rho_pm_IP = sqrt(pow(Kchrec[6] - ip_avg[0], 2) + pow(Kchrec[7] - ip_avg[1], 2));
		rho_00_IP = sqrt(pow(Knetri_kinfit[6] - ip_avg[0], 2) + pow(Knetri_kinfit[7] - ip_avg[1], 2));
		rho = sqrt(pow(rho_pm_IP, 2) + pow(rho_00_IP, 2));

		lengthCh = sqrt( pow(Kchrec[6] - ip_avg[0],2) +
										 pow(Kchrec[7] - ip_avg[1],2) +
										 pow(Kchrec[8] - ip_avg[2],2)	);

		for (Int_t ii = 0; ii < 4; ii++)
		{
			lengthKaon[ii] = sqrt( pow(gamma_kinfit_triangle[ii][4] - ip_avg[0],2) +
												  	pow(gamma_kinfit_triangle[ii][5] - ip_avg[1],2) +
												  	pow(gamma_kinfit_triangle[ii][6] - ip_avg[2],2)	) - (gamma_kinfit_triangle[ii][7] / cVel);
		}
		

		// Calculation of pi+pi-pi0 invariant mass - best option chosen by quadrature.
		for (Int_t j = 0; j < 2; j++)
		{
			M_omega_tmp[j] = sqrt(pow(PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[j][3], 2) -
														pow(PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[j][0], 2) -
														pow(PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[j][1], 2) -
														pow(PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[j][2], 2));
			M_omega_diff[j] = M_omega_tmp[j] - mOmega;
		}

		Float_t
				Pi0NonOmega[4];
		//
		if (std::isnan(M_omega_diff[0]) || std::isnan(M_omega_diff[1]))
			minv_omega = 999999.;
		else if (std::isinf(M_omega_diff[0]) || std::isinf(M_omega_diff[1]))
			minv_omega = 999999.;
		else
		{
			if (abs(M_omega_diff[0]) < abs(M_omega_diff[1]))
			{
				minv_omega = M_omega_tmp[0];
				OmegaMom[0] = PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[0][0];
				OmegaMom[1] = PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[0][1];
				OmegaMom[2] = PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[0][2];
				OmegaMom[3] = PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[0][3];

				Pi0NonOmega[0] = Pi0Mom[1][0];
				Pi0NonOmega[1] = Pi0Mom[1][1];
				Pi0NonOmega[2] = Pi0Mom[1][2];
				Pi0NonOmega[3] = Pi0Mom[1][3];
			}
			else
			{
				minv_omega = M_omega_tmp[1];
				OmegaMom[0] = PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[1][0];
				OmegaMom[1] = PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[1][1];
				OmegaMom[2] = PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[1][2];
				OmegaMom[3] = PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[1][3];

				Pi0NonOmega[0] = Pi0Mom[0][0];
				Pi0NonOmega[1] = Pi0Mom[0][1];
				Pi0NonOmega[2] = Pi0Mom[0][2];
				Pi0NonOmega[3] = Pi0Mom[0][3];
			}
		}
		//
		// Lorentz transformation of Pich to Kaon CM frame
		Float_t
				boost_vec_Kchboost[3] = {-Kchboost[0] / Kchboost[3],
																 -Kchboost[1] / Kchboost[3],
																 -Kchboost[2] / Kchboost[3]},
				PichFourMomKaonCM[2][4];

		lorentz_transf(boost_vec_Kchboost, PichFourMom[0], PichFourMomKaonCM[0]);
		lorentz_transf(boost_vec_Kchboost, PichFourMom[1], PichFourMomKaonCM[1]);

		TVector3
				pich1(PichFourMomKaonCM[0][0], PichFourMomKaonCM[0][1], PichFourMomKaonCM[0][2]),
				pich2(PichFourMomKaonCM[1][0], PichFourMomKaonCM[1][1], PichFourMomKaonCM[1][2]);

		anglePichKaonCM = pich1.Angle(pich2) * 180. / M_PI;

		// Lorentz transformation of Pi0 to Kaon CM frame
		Float_t
				boost_vec_Kne[3] = {-(baseKin.phi_mom[0] - Kchboost[0]) / (baseKin.phi_mom[3] - Kchboost[3]),
														-(baseKin.phi_mom[1] - Kchboost[1]) / (baseKin.phi_mom[3] - Kchboost[3]),
														-(baseKin.phi_mom[2] - Kchboost[2]) / (baseKin.phi_mom[3] - Kchboost[3])},
				Pi0KaonCM[2][4];

		lorentz_transf(boost_vec_Kne, Pi0Mom[0], Pi0KaonCM[0]);
		lorentz_transf(boost_vec_Kne, Pi0Mom[1], Pi0KaonCM[1]);

		TVector3
				pi01(Pi0KaonCM[0][0], Pi0KaonCM[0][1], Pi0KaonCM[0][2]),
				pi02(Pi0KaonCM[1][0], Pi0KaonCM[1][1], Pi0KaonCM[1][2]);

		anglePi0KaonCM = pi01.Angle(pi02) * 180. / M_PI;

		// Lorentz transformation of Pi0 to Kaon CM frame
		Float_t
				boost_vec_phi[3] = {-(baseKin.phi_mom[0]) / (baseKin.phi_mom[3]),
														-(baseKin.phi_mom[1]) / (baseKin.phi_mom[3]),
														-(baseKin.phi_mom[2]) / (baseKin.phi_mom[3])},
				Pi0NonOmegaCM[4],
				OmegaMomCM[4];

		lorentz_transf(boost_vec_phi, Pi0NonOmega, Pi0NonOmegaCM);
		lorentz_transf(boost_vec_phi, OmegaMom, OmegaMomCM);

		TVector3
				pi0CM(Pi0NonOmegaCM[0], Pi0NonOmegaCM[1], Pi0NonOmegaCM[2]),
				omegaCM(OmegaMomCM[0], OmegaMomCM[1], OmegaMomCM[2]);

		anglePi0OmegaPhiCM = pi0CM.Angle(omegaCM) * 180. / M_PI;

		// Angle between phi and omega

		TVector3
				phi(baseKin.phi_mom[0], baseKin.phi_mom[1], baseKin.phi_mom[2]),
				omega(OmegaMom[0], OmegaMom[1], OmegaMom[2]);

		anglePhiOmega = phi.Angle(omega) * 180. / M_PI;

		// Section for the histogram filling

		if(baseKin.mctruth_int == 4)
		{
			hist_control[0]->Fill(lengthKaon[0]);
			hist_control[0]->Fill(lengthKaon[1]);
			hist_control[0]->Fill(lengthKaon[2]);
			hist_control[0]->Fill(lengthKaon[3]);
			hist_control[1]->Fill(lengthCh);
			hist_control[2]->Fill(Kchrec[8] - Knetri_kinfit[8]);

			std::cout << lengthKaon[1] << std::endl;
		}

		Int_t
				sigmas_neu = 1.,
				sigmas_ip = 3.;

		Bool_t
				cond[6],
				cond_tot;

		cond[0] = abs(Kchrec[6] - Knetri_kinfit[6]) < sigmas_neu * 16.89;
		cond[1] = abs(Kchrec[7] - Knetri_kinfit[7]) < sigmas_neu * 16.77;
		cond[2] = abs(Kchrec[8] - Knetri_kinfit[8]) < sigmas_neu * 9.89;
		cond[3] = abs(Kchrec[6] - ip_avg[0]) < sigmas_ip * 2.63;
		cond[4] = abs(Kchrec[7] - ip_avg[1]) < sigmas_ip * 2.03;
		cond[5] = abs(Kchrec[8] - ip_avg[2]) < sigmas_ip * 3.39;

		cond_tot = 1;

		if (cond_tot)
		{

			if (baseKin.mctruth_int == 1)
			{
				hist[0][0]->Fill(anglePi0OmegaPhiCM);
				hist[0][1]->Fill(Kchrec[5]);
				hist[0][2]->Fill(minv_omega);
				hist[0][3]->Fill(anglePi0OmegaPhiCM);

				hist_control_chann[0][0]->Fill(rho_00);
				hist_control_chann[0][1]->Fill(Kchrec[8] - Knetri_kinfit[8]);

				hist_pm_IP[0][0]->Fill(rho_pm_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist_pm_IP[0][1]->Fill(Kchrec[8] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));

				hist_00_IP[0][0]->Fill(rho_00_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist_00_IP[0][1]->Fill(Knetri_kinfit[8] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));

				hist2d[0]->Fill(anglePichKaonCM, rho);

				hist2d_00IP_pmIP[0][0]->Fill(rho_pm_IP, rho_00_IP, event.interf_function(baseKin.Dtmc, 0, par));
				hist2d_00IP_pmIP[1][0]->Fill(Kchrec[8] - ip_avg[8], Knetri_kinfit[8] - ip_avg[8], event.interf_function(baseKin.Dtmc, 0, par));
			}
			if (baseKin.mctruth_int == 3)
			{
				hist[1][0]->Fill(anglePi0OmegaPhiCM);
				hist[1][1]->Fill(Kchrec[5]);
				hist[1][2]->Fill(minv_omega);
				hist[1][3]->Fill(anglePi0OmegaPhiCM);

				hist_control_chann[1][0]->Fill(rho_00);
				hist_control_chann[1][1]->Fill(Kchrec[8] - Knetri_kinfit[8]);

				hist_pm_IP[1][0]->Fill(rho_pm_IP);
				hist_pm_IP[1][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[1][0]->Fill(rho_00_IP);
				hist_00_IP[1][1]->Fill(Knetri_kinfit[8] - ip_avg[8]);

				hist2d_00IP_pmIP[0][1]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][1]->Fill(Kchrec[8] - ip_avg[8], Knetri_kinfit[8] - ip_avg[8]);

				hist2d[1]->Fill(anglePichKaonCM, rho);
			}
			if (baseKin.mctruth_int == 4)
			{
				hist[2][0]->Fill(anglePi0OmegaPhiCM);
				hist[2][1]->Fill(Kchrec[5]);
				hist[2][2]->Fill(minv_omega);
				hist[2][3]->Fill(anglePi0OmegaPhiCM);

				hist_control_chann[2][0]->Fill(rho_00);
				hist_control_chann[2][1]->Fill(Kchrec[8] - Knetri_kinfit[8]);

				hist_pm_IP[2][0]->Fill(rho_pm_IP);
				hist_pm_IP[2][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[2][0]->Fill(rho_00_IP);
				hist_00_IP[2][1]->Fill(Knetri_kinfit[8] - ip_avg[8]);

				hist2d_00IP_pmIP[0][2]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][2]->Fill(Kchrec[8] - ip_avg[8], Knetri_kinfit[8] - ip_avg[8]);

				hist2d[2]->Fill(anglePichKaonCM, rho);
			}
			if (baseKin.mctruth_int == 5)
			{
				hist[3][0]->Fill(anglePi0OmegaPhiCM);
				hist[3][1]->Fill(Kchrec[5]);
				hist[3][2]->Fill(minv_omega);
				hist[3][3]->Fill(anglePi0OmegaPhiCM);

				hist_control_chann[3][0]->Fill(rho_00);
				hist_control_chann[3][1]->Fill(Kchrec[8] - Knetri_kinfit[8]);

				hist_pm_IP[3][0]->Fill(rho_pm_IP);
				hist_pm_IP[3][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[3][0]->Fill(rho_00_IP);
				hist_00_IP[3][1]->Fill(Knetri_kinfit[8] - ip_avg[8]);

				hist2d[3]->Fill(anglePichKaonCM, rho);

				// hist2d_00IP_pmIP[0][0]->Fill(rho_pm_IP, rho_00_IP);

				hist2d_00IP_pmIP[0][3]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][3]->Fill(Kchrec[8] - ip_avg[8], Knetri_kinfit[8] - ip_avg[8]);
			}
			if (baseKin.mctruth_int == 6)
			{
				hist[4][0]->Fill(anglePi0OmegaPhiCM);
				hist[4][1]->Fill(Kchrec[5]);
				hist[4][2]->Fill(minv_omega);
				hist[4][3]->Fill(anglePi0OmegaPhiCM);

				hist_control_chann[4][0]->Fill(rho_00);
				hist_control_chann[4][1]->Fill(Kchrec[8] - Knetri_kinfit[8]);

				hist_pm_IP[4][0]->Fill(rho_pm_IP);
				hist_pm_IP[4][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[4][0]->Fill(rho_00_IP);
				hist_00_IP[4][1]->Fill(Knetri_kinfit[8] - ip_avg[8]);

				hist2d[4]->Fill(anglePichKaonCM, rho);

				hist2d_00IP_pmIP[0][4]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][4]->Fill(Kchrec[8] - ip_avg[8], Knetri_kinfit[8] - ip_avg[8]);
			}
			if (baseKin.mctruth_int == 7)
			{
				hist[5][0]->Fill(anglePi0OmegaPhiCM);
				hist[5][1]->Fill(Kchrec[5]);
				hist[5][2]->Fill(minv_omega);
				hist[5][3]->Fill(anglePi0OmegaPhiCM);

				hist_control_chann[5][0]->Fill(rho_00);
				hist_control_chann[5][1]->Fill(Kchrec[8] - Knetri_kinfit[8]);

				hist_pm_IP[5][0]->Fill(rho_pm_IP);
				hist_pm_IP[5][1]->Fill(Kchrec[8] - ip_avg[8]);

				hist_00_IP[5][0]->Fill(rho_00_IP);
				hist_00_IP[5][1]->Fill(Knetri_kinfit[8] - ip_avg[8]);

				hist2d[5]->Fill(anglePichKaonCM, rho);

				hist2d_00IP_pmIP[0][5]->Fill(rho_pm_IP, rho_00_IP);
				hist2d_00IP_pmIP[1][5]->Fill(Kchrec[8] - ip_avg[8], Knetri_kinfit[8] - ip_avg[8]);
			}
		}
		tree->Fill();
	}

	hist_00_IP[0][0]->Scale(hist_00_IP[0][0]->GetEntries() / hist_00_IP[0][0]->Integral());
	hist_pm_IP[0][0]->Scale(hist_pm_IP[0][0]->GetEntries() / hist_pm_IP[0][0]->Integral());

	hist2d_00IP_pmIP[0][0]->Scale(hist2d_00IP_pmIP[0][0]->GetEntries() / hist2d_00IP_pmIP[0][0]->Integral());

	TLegend *legend_chann = new TLegend(0.1, 0.5, 0.4, 0.9);
	legend_chann->SetFillColor(kWhite);

	TLegend *legend_chi2 = new TLegend(0.6, 0.5, 0.9, 0.9);
	legend_chi2->SetFillColor(kWhite);

	TString xTitle[4] = {"#angle(#pi^{0},#omega) [#circ]", "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]", "m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]", "#angle(#pi^{0},#omega) [#circ]"};

	TString yTitle = "Counts";

	for (Int_t i = 0; i < 4; i++)
	{
		canva[i]->cd();
		canva[i]->SetLogy(1);
		for (Int_t j = 0; j < channNum; j++)
		{
			hist[j][i]->SetLineColor(channColor[j]);
			hist[j][i]->SetLineWidth(3);

			legend_chann->AddEntry(hist[j][i], channName[j], "l");
			legend_chi2->AddEntry(hist[j][i], channName[j], "l");

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
		if (i < 3)
			legend_chann->Draw();
		else if (i == 3)
			legend_chi2->Draw();

		canva[i]->Print(img_dir + "hist" + std::to_string(i) + ext_img);
		legend_chann->Clear();
		legend_chi2->Clear();
	}

	for (Int_t j = 0; j < channNum; j++)
	{
		canva2d[j]->cd();
		canva2d[j]->SetLogz(1);

		hist2d[j]->SetMarkerColor(channColor[j]);
		hist2d[j]->SetMarkerSize(2);

		hist2d[j]->GetXaxis()->SetTitle("#angle(#pi^{+},#pi^{-}) [#circ]");
		hist2d[j]->GetYaxis()->SetTitle("#sqrt{#rho_{+-}^{2} + #rho_{00}^{2}} [cm]");
		hist2d[j]->Draw("COLZ");

		canva2d[j]->Print(img_dir + "angle_2d_kaon_cm_" + channName[j] + ext_img);
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

		triple_fit->SetParameters(0.15 * parameter[0], parameter[1], parameter[2], 0.15 * parameter[0], parameter[1] - 10.0, parameter[2], parameter[0], parameter[1] + 10.0, parameter[2]);

		triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
		triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
		triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

		triple_fit->SetParLimits(1, parameter[1] - 2.0, parameter[1] + 2.0);
		triple_fit->SetParLimits(4, parameter[1] - 2.0, parameter[1] + 2.0);
		triple_fit->SetParLimits(7, parameter[1] - 2.0, parameter[1] + 2.0);

		triple_fit->SetParLimits(2, 0.01 * parameter[2], 5.0 * parameter[2]);
		triple_fit->SetParLimits(5, 0.01 * parameter[2], 5.0 * parameter[2]);
		triple_fit->SetParLimits(8, 0.01 * parameter[2], 5.0 * parameter[2]);

		result = hist_control[i]->Fit(triple_fit, "SF");

		fit_stats[1] = Form("Mean = %.2f#pm%.2f", comb_mean(result->GetParams(), result->GetErrors()), comb_mean_err(result->GetParams(), result->GetErrors()));
		fit_stats[2] = Form("Std Dev = %.2f#pm%.2f", comb_std_dev(result->GetParams(), result->GetErrors()), comb_std_dev_err(result->GetParams(), result->GetErrors()));

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

		canvas_cont[i]->Print(img_dir + "cont_plot_" + std::to_string(i) + ext_img);

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

		canva[i + 6]->Print(img_dir + "hist_control_chann" + std::to_string(i) + ext_img);
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

		canva[i + 6]->Print(img_dir + "hist_pm_IP_" + std::to_string(i) + ext_img);
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

		canva[i + 8]->Print(img_dir + "hist_00_IP_" + std::to_string(i) + ext_img);
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

		canva2d[j + channNum]->Print(img_dir + "rho_2d_" + channName[j] + ext_img);
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

		canva2d[j + 2 * channNum]->Print(img_dir + "z_2d_" + channName[j] + ext_img);
	}

	tree->Print();

	file->Write();
	file->Close();
	delete file;

	return 0;
}
