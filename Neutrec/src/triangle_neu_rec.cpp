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

int triangle_neurec(int first_file, int last_file, int loopcount, int M, int range, Controls::DataType data_type, int good_clus)
{
	ErrorHandling::ErrorLogs errLogger;
  LogsHandling::Logs logger;
  std::ofstream LogFileMain, LogFileTriangle, LogFileTri, LogFileTriKinFit;
  std::ofstream ErrFileMain, ErrFileTriangle, ErrFileTri, ErrFileTriKinFit;

	TString
			filename_gen = gen_vars_dir + root_files_dir + gen_vars_filename + first_file + "_" + last_file + ext_root,
			filename_trilateration = neutrec_dir + root_files_dir + neu_trilateration_kin_fit_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	try
	{
		TFile *file_gen = new TFile(filename_gen);
		TFile *file_trilateration = new TFile(filename_trilateration);

		if (file_gen->IsZombie() || file_trilateration->IsZombie())
		{
			throw ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
		}

		TTree *tree_gen = (TTree *)file_gen->Get(gen_vars_tree);
		TTree *tree_trilateration = (TTree *)file_trilateration->Get(neutrec_kin_fit_tree);

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
		UChar_t mctruth, mcflag;
		Float_t cluster[5][500], Kchboost[9], Knereclor[9], Knemc[9], ip[3], Kchrec[9];

		chain->SetBranchAddress("nclu", &nclu);
		chain->SetBranchAddress("Xacl", cluster[0]);
		chain->SetBranchAddress("Yacl", cluster[1]);
		chain->SetBranchAddress("Zacl", cluster[2]);
		chain->SetBranchAddress("Tcl", cluster[3]);
		chain->SetBranchAddress("Enecl", cluster[4]);

		chain->SetBranchAddress("mctruth", &mctruth);
		chain->SetBranchAddress("mcflag", &mcflag);

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

		Int_t done_kinfit = 0, g4taken_kinfit[4] = {0};

		tree_trilateration->SetBranchAddress("done4_kinfit", &done_kinfit);
		tree_trilateration->SetBranchAddress("g4takentri_kinfit", g4taken_kinfit);

		// Adding friends to the chain
		chain->AddFriend(tree_gen);
		chain->AddFriend(tree_trilateration);

		Int_t nentries = (Int_t)chain->GetEntries();

		// Creation of the triangle method tree
		TString name = neutrec_dir + root_files_dir + neu_triangle_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

		TFile *file = new TFile(name, "recreate");
		TTree *tree = new TTree(neutrec_triangle_tree, "Neu vtx rec with triangle method");

		// Variables
		Bool_t cond_ene, cond_clus[4];
		Int_t ind_gam[4];
		Double_t solution[4];
		Float_t neu_vtx[4] = {0.};

		Int_t done, fourg4taken[4], chosen;
		Float_t gammatriangle[4][8], Knetriangle[10], totalerr, ip_triangle[3], sol1[4], sol2[4], sol1err, sol2err;

		TBranch *b_gamma1triangle = tree->Branch("fourgamma1triangle", gammatriangle[0], "fourgamma1triangle[8]/F");
		TBranch *b_gamma2triangle = tree->Branch("fourgamma2triangle", gammatriangle[1], "fourgamma2triangle[8]/F");
		TBranch *b_gamma3triangle = tree->Branch("fourgamma3triangle", gammatriangle[2], "fourgamma3triangle[8]/F");
		TBranch *b_gamma4triangle = tree->Branch("fourgamma4triangle", gammatriangle[3], "fourgamma4triangle[8]/F");

		TBranch *b_Knetriangle = tree->Branch("fourKnetriangle", Knetriangle, "fourKnetriangle[10]/F");
		TBranch *b_iptriangle = tree->Branch("iptriangle", ip_triangle, "iptriangle[3]/F");

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

		Bool_t
				data_flag = false;

		for (Int_t i = 0; i < nentries; i++)
		{
			chain->GetEntry(i);

			for (Int_t j1 = 0; j1 < 4; j1++)
				for (Int_t j2 = 0; j2 < 8; j2++)
					gammatriangle[j1][j2] = -999.;

			for (Int_t j = 0; j < 10; j++)
				Knetriangle[j] = -999.;

			for (Int_t j = 0; j < 4; j++)
				fourg4taken[j] = -999;

			totalerr = -999.;

			done = 0;

			vtxSigmaMin = 999999.;
			TrcSumMin = 999999.;

			// Setting the data type flags

			dataFlagSetter(data_type, data_flag, int(mcflag), int(mctruth));

			if (done_kinfit == 1 && data_flag)
			{
				if (good_clus == 0)
				{
					ind_gam[0] = g4taken_kinfit[0];
					ind_gam[1] = g4taken_kinfit[1];
					ind_gam[2] = g4taken_kinfit[2];
					ind_gam[3] = g4taken_kinfit[3];

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

						//! Using the charged part of the decay

						Knerec[0] = bhabha_mom[0] - Kchrec[0];
						Knerec[1] = bhabha_mom[1] - Kchrec[1];
						Knerec[2] = bhabha_mom[2] - Kchrec[2];
						Knerec[3] = bhabha_mom[3] - Kchrec[3];

						//!

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

				tree->Fill();
			}
		}

		tree->Print();

		tree->Write();
		file->Close();
		delete file;
	}
	catch (ErrorHandling::ErrorCodes err)
	{
		errLogger.getErrLog(err, ErrFileTriangle);
		ErrFileTriangle.close();

		errLogger.getErrLog(err);
		return int(err);
	}

	return 0;
}
