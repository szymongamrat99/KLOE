#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "../../../Include/Codes/reconstructor.h"
#include "../../../Include/Codes/neu_triangle.h"
#include "../../../Include/Codes/plane_intersection.h"
#include "../inc/trilateration.hpp"

int triangle_neurec(int data_type, int first_file, int last_file, int good_clus) //	arguments are: 1. Number of points
																																	//                2. First analyzed file
																																	//                3. Last analyzed file

{
	ErrorHandling::ErrorLogs logger;
	std::ofstream LogFile;
	LogFile.open("TriangleNeuVtx.log");

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	try
	{
		TFile *file_gen = new TFile(gen_vars_dir + gen_vars_filename + first_file + "_" + last_file + ext_root);

		if (!file_gen)
		{
			throw ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
		}

		TTree *tree_gen = (TTree *)file_gen->Get(gen_vars_tree);

		TString name = neu_triangle_filename + std::to_string(first_file) + "_" + std::to_string(last_file) + ext_root;

		TFile *file = new TFile(name, "recreate");
		TTree *tree = new TTree(neutrec_triangle_tree, "Neu vtx rec with triangle");

		// Branches' addresses
		// Bhabha vars
		Float_t bhabha_mom[4], bhabha_vtx[3];

		chain->SetBranchAddress("Bpx", &bhabha_mom[0]);
		chain->SetBranchAddress("Bpy", &bhabha_mom[1]);
		chain->SetBranchAddress("Bpz", &bhabha_mom[2]);
		chain->SetBranchAddress("Broots", &bhabha_mom[3]);

		chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
		chain->SetBranchAddress("By", &bhabha_vtx[1]);
		chain->SetBranchAddress("Bz", &bhabha_vtx[2]);

		// Cluster vars
		Int_t nclu;
		UChar_t mctruth;
		Float_t cluster[5][500], Kchboost[9], Knereclor[9], Knemc[9], ip[3], Kchrec[9];

		chain->SetBranchAddress("nclu", &nclu);
		chain->SetBranchAddress("Xacl", cluster[0]);
		chain->SetBranchAddress("Yacl", cluster[1]);
		chain->SetBranchAddress("Zacl", cluster[2]);
		chain->SetBranchAddress("Tcl", cluster[3]);
		chain->SetBranchAddress("Enecl", cluster[4]);

		chain->SetBranchAddress("mctruth", &mctruth);

		chain->SetBranchAddress("ip", ip);
		chain->SetBranchAddress("Kchrec", Kchrec);
		chain->SetBranchAddress("Knemc", Knemc);

		Int_t good_clus_ind[4];
		Float_t pgammc[4][8];

		tree_gen->SetBranchAddress("pgammc1", pgammc[0]);
		tree_gen->SetBranchAddress("pgammc2", pgammc[1]);
		tree_gen->SetBranchAddress("pgammc3", pgammc[2]);
		tree_gen->SetBranchAddress("pgammc4", pgammc[3]);

		tree_gen->SetBranchAddress("clusindgood", good_clus_ind);

		chain->AddFriend(tree_gen);

		Int_t nentries = (Int_t)chain->GetEntries();

		// Variables
		Bool_t cond_ene, cond_clus[4], cond_sol[2][2], onesol;
		Int_t ind_gam[4], error = 0, index = 0;
		Float_t ene_sum, gamma_mom[2][6][4], gamma_len[2][6],
				kaon_mom[2][4], kaon_vel[2], kaon_len[2], path_diff[2],
				sol_tmp[4], total_err = 0., total_err_def = 0.,
										gamma_len_tmp[4], kaon_mom_tmp[4], kaon_inv_mass, ene_sum_tmp;
		Double_t solution[4];
		Float_t gammatri_tmp[4][8], Knetri_tmp[10];

		Float_t y_axis[3], neu_vtx[4] = {0.};

		Int_t done, fourg4taken[4], chosen;
		Float_t gammatri[4][8], Knetri[10], totalerr, ip_tri[3], sol1[4], sol2[4], sol1err, sol2err;

		TBranch *b_gamma1tri = tree->Branch("fourgamma1tri", gammatri[0], "fourgamma1tri[8]/F");
		TBranch *b_gamma2tri = tree->Branch("fourgamma2tri", gammatri[1], "fourgamma2tri[8]/F");
		TBranch *b_gamma3tri = tree->Branch("fourgamma3tri", gammatri[2], "fourgamma3tri[8]/F");
		TBranch *b_gamma4tri = tree->Branch("fourgamma4tri", gammatri[3], "fourgamma4tri[8]/F");

		TBranch *b_Knetri = tree->Branch("fourKnetri", Knetri, "fourKnetri[10]/F");
		TBranch *b_iptri = tree->Branch("iptri", ip_tri, "iptri[3]/F");

		TBranch *b_done = tree->Branch("done4", &done, "done4/I");

		TBranch *b_sol1 = tree->Branch("sol1", sol1, "sol1[4]/F");
		TBranch *b_sol2 = tree->Branch("sol2", sol2, "sol2[4]/F");

		TBranch *b_sol1err = tree->Branch("sol1err", &sol1err, "sol1err/F");
		TBranch *b_sol2err = tree->Branch("sol2err", &sol2err, "sol2err/F");

		TBranch *b_chosen = tree->Branch("chosen", &chosen, "chosen/I");

		TBranch *b_fourg4taken = tree->Branch("fourg4taken", fourg4taken, "fourg4taken[4]/I");

		Float_t distance = 0., distance_mc = 0., distance_diff = 0.;
		Float_t mc_dist[90], sigmas[90], Knerec[9], trc[4] = {0.};

		Float_t Clu5Vec[4][5];

		Float_t eqn_check[4], eqn_check_tmp[2][4], eqn_check_tmp_tot[2], kaon_inv_mass_tmp[2], vtxSigmaMin, TrcSumMin;

		for (Int_t i = 0; i < nentries; i++)
		{
			chain->GetEntry(i);

			for (Int_t j1 = 0; j1 < 4; j1++)
				for (Int_t j2 = 0; j2 < 8; j2++)
					gammatri[j1][j2] = -999.;

			for (Int_t j = 0; j < 10; j++)
				Knetri[j] = -999.;

			for (Int_t j = 0; j < 4; j++)
				fourg4taken[j] = -999;

			totalerr = -999.;

			done = 0;

			total_err_def = 999999.;

			vtxSigmaMin = 999999.;
			TrcSumMin = 999999.;

			if (nclu >= 4 && (mctruth == 1))
			{
				if (good_clus == 0)
				{
					for (Int_t j1 = 0; j1 < nclu - 3; j1++)
						for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
							for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
								for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
								{
									ene_sum = 0.;
									error = 2;
									total_err = 999999.;

									ind_gam[0] = j1;
									ind_gam[1] = j2;
									ind_gam[2] = j3;
									ind_gam[3] = j4;

									cond_ene = cluster[4][ind_gam[0]] > MIN_CLU_ENE && cluster[4][ind_gam[1]] > MIN_CLU_ENE &&
														 cluster[4][ind_gam[2]] > MIN_CLU_ENE && cluster[4][ind_gam[3]] > MIN_CLU_ENE;

									cond_clus[0] = cluster[3][ind_gam[0]] > 0 && cluster[0][ind_gam[0]] != 0 && cluster[1][ind_gam[0]] != 0 && cluster[2][ind_gam[0]] != 0;
									cond_clus[1] = cluster[3][ind_gam[1]] > 0 && cluster[0][ind_gam[1]] != 0 && cluster[1][ind_gam[1]] != 0 && cluster[2][ind_gam[1]] != 0;
									cond_clus[2] = cluster[3][ind_gam[2]] > 0 && cluster[0][ind_gam[2]] != 0 && cluster[1][ind_gam[2]] != 0 && cluster[2][ind_gam[2]] != 0;
									cond_clus[3] = cluster[3][ind_gam[3]] > 0 && cluster[0][ind_gam[3]] != 0 && cluster[1][ind_gam[3]] != 0 && cluster[2][ind_gam[3]] != 0;

									if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
									{
										for (Int_t k = 0; k < 4; k++)
										{
											Clu5Vec[k][0] = cluster[0][ind_gam[k]];
											Clu5Vec[k][1] = cluster[1][ind_gam[k]];
											Clu5Vec[k][2] = cluster[2][ind_gam[k]];
											Clu5Vec[k][3] = cluster[3][ind_gam[k]];
											Clu5Vec[k][4] = cluster[4][ind_gam[k]];
										}

										Knerec[0] = bhabha_mom[0] - Kchrec[0];
										Knerec[1] = bhabha_mom[1] - Kchrec[1];
										Knerec[2] = bhabha_mom[2] - Kchrec[2];
										Knerec[3] = bhabha_mom[3] - Kchrec[3];

										Float_t TrcSum = 0., vtxSigma = 0.;

										neu_triangle(&TrcSum, &vtxSigma, Clu5Vec, ip, bhabha_mom, Knerec, neu_vtx, trc);

										if (sqrt(pow(vtxSigma, 2) + pow(TrcSum, 2)) < sqrt(pow(vtxSigmaMin, 2) + pow(TrcSumMin, 2)))
										{
											vtxSigmaMin = vtxSigma;
											TrcSumMin = TrcSum;

											solution[0] = neu_vtx[0];
											solution[1] = neu_vtx[1];
											solution[2] = neu_vtx[2];
											solution[3] = neu_vtx[3];
										}
									}
								}
				}

				tree->Fill();

				std::cout << solution[0] << " " << solution[1] << " " << solution[2] << " " << solution[3] << std::endl;
				std::cout << Knemc[6] << " " << Knemc[7] << " " << Knemc[8] << std::endl
									<< std::endl;
			}
		}

		tree->Print();

		tree->Write();
		file->Close();
		delete file;
	}
	catch (ErrorHandling::ErrorCodes err)
	{
		logger.getErrLog(err, LogFile);
		return int(err);
	}

	return 0;
}
