#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>

#include <reconstructor.h>
#include <neu_triangle.h>
#include <plane_intersection.h>
#include <neutral_mom.h>
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
		UChar_t mctruth, mcflag, g4taken[4], ncll[50];
		Float_t cluster[5][500], Kchboost[9], Knereclor[9], Knemc[9], ip[3];
		BaseKinematics baseKin;

		chain->SetBranchAddress("nclu", &nclu);
		chain->SetBranchAddress("Xacl", cluster[0]);
		chain->SetBranchAddress("Yacl", cluster[1]);
		chain->SetBranchAddress("Zacl", cluster[2]);
		chain->SetBranchAddress("Tcl", cluster[3]);
		chain->SetBranchAddress("Enecl", cluster[4]);

		chain->SetBranchAddress("mctruth", &mctruth);
		chain->SetBranchAddress("mcflag", &mcflag);

		chain->SetBranchAddress("ip", ip);
		chain->SetBranchAddress("Kchboost", Kchboost);
		chain->SetBranchAddress("Knereclor", Knereclor);
		chain->SetBranchAddress("g4taken", g4taken);
		chain->SetBranchAddress("ncll", baseKin.ncll);

		Int_t good_clus_ind[4];
		Float_t pgammc[4][8];

		tree_gen->SetBranchAddress("pgammc1", pgammc[0]);
		tree_gen->SetBranchAddress("pgammc2", pgammc[1]);
		tree_gen->SetBranchAddress("pgammc3", pgammc[2]);
		tree_gen->SetBranchAddress("pgammc4", pgammc[3]);

		tree_gen->SetBranchAddress("clusindgood", good_clus_ind);

		Int_t done_kinfit = 0, g4taken_kinfit[4] = {0};
		Float_t chi2min_tri, Knetri_kinfit[10];

		tree_trilateration->SetBranchAddress("done4_kinfit", &done_kinfit);
		tree_trilateration->SetBranchAddress("g4takentri_kinfit", g4taken_kinfit);
		tree_trilateration->SetBranchAddress("chi2min", &chi2min_tri);
		tree_trilateration->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

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
		Float_t gammatriangle[4][4], Knetriangle[10], totalerr, ip_triangle[3], chi2min, trcsum, trcfinal[4], minv4gam;

		TBranch *b_gamma1triangle = tree->Branch("fourgamma1triangle", gammatriangle[0], "fourgamma1triangle[4]/F");
		TBranch *b_gamma2triangle = tree->Branch("fourgamma2triangle", gammatriangle[1], "fourgamma2triangle[4]/F");
		TBranch *b_gamma3triangle = tree->Branch("fourgamma3triangle", gammatriangle[2], "fourgamma3triangle[4]/F");
		TBranch *b_gamma4triangle = tree->Branch("fourgamma4triangle", gammatriangle[3], "fourgamma4triangle[4]/F");

		TBranch *b_Knetriangle = tree->Branch("fourKnetriangle", Knetriangle, "fourKnetriangle[10]/F");
		TBranch *b_iptriangle = tree->Branch("iptriangle", ip_triangle, "iptriangle[3]/F");

		TBranch *b_done = tree->Branch("done_triangle", &done, "done_triangle/I");

		TBranch *b_chi2min = tree->Branch("chi2min", &chi2min, "chi2min/F");

		TBranch *b_trcsum = tree->Branch("trcsum", &trcsum, "trcsum/F");
		TBranch *b_trc = tree->Branch("trc", trcfinal, "trc/F");

		TBranch *b_minv4gam = tree->Branch("minv4gam", &minv4gam, "minv4gam/F");

		TBranch *b_fourg4taken = tree->Branch("g4taken_triangle", fourg4taken, "g4taken_triangle[4]/I");

		Float_t distance = 0., distance_mc = 0., distance_diff = 0.;
		Float_t mc_dist[90], sigmas[90], Knerec[9], trc[4] = {0.};

		Float_t Clu5Vec[4][5];

		Float_t eqn_check[4], eqn_check_tmp[2][4], eqn_check_tmp_tot[2], kaon_inv_mass_tmp[2], vtxSigmaMin, TrcSumMin;

		Bool_t
				data_flag = false;

		Int_t
				nbin = 100;

		Float_t
					x_lim[4][2] = {
												 {-50.0, 50.0},
												 {-50.0, 50.0},
												 {-50.0, 50.0},
												 {-10.0, 20.0}
					};

		std::vector<TH2*> test2d;
		TString test2d_name = "";
		
		for (Int_t i = 0; i < 4; i++)
		{
			test2d_name = "test2d_" + i;
			test2d.push_back(new TH2D(test2d_name, "", nbin, x_lim[i][0], x_lim[i][1], nbin, x_lim[i][0], x_lim[i][1]));
		};

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
			chi2min = chi2min_tri;

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

					cond_ene = cluster[4][baseKin.ncll[ind_gam[0]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[1]] - 1] > MIN_CLU_ENE &&
										 cluster[4][baseKin.ncll[ind_gam[2]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[3]] - 1] > MIN_CLU_ENE;

					cond_clus[0] = cluster[3][baseKin.ncll[ind_gam[0]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[0]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[0]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[0]] - 1] != 0;
					cond_clus[1] = cluster[3][baseKin.ncll[ind_gam[1]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[1]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[1]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[1]] - 1] != 0;
					cond_clus[2] = cluster[3][baseKin.ncll[ind_gam[2]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[2]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[2]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[2]] - 1] != 0;
					cond_clus[3] = cluster[3][baseKin.ncll[ind_gam[3]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[3]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[3]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[3]] - 1] != 0;

					if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
					{
						for (Int_t k = 0; k < 4; k++)
						{
							Clu5Vec[k][0] = cluster[0][baseKin.ncll[ind_gam[k]] - 1];
							Clu5Vec[k][1] = cluster[1][baseKin.ncll[ind_gam[k]] - 1];
							Clu5Vec[k][2] = cluster[2][baseKin.ncll[ind_gam[k]] - 1];
							Clu5Vec[k][3] = cluster[3][baseKin.ncll[ind_gam[k]] - 1];
							Clu5Vec[k][4] = cluster[4][baseKin.ncll[ind_gam[k]] - 1];
						}

						//! Using the charged part of the decay

						Knerec[0] = bhabha_mom[0] - Kchboost[0];
						Knerec[1] = bhabha_mom[1] - Kchboost[1];
						Knerec[2] = bhabha_mom[2] - Kchboost[2];
						Knerec[3] = bhabha_mom[3] - Kchboost[3];

						//!

						Float_t TrcSum = 0., vtxSigma = 0.;

						neu_triangle(&TrcSum, &vtxSigma, Clu5Vec, ip, bhabha_mom, Knerec, neu_vtx, trc);

						if (sqrt(pow(vtxSigma, 2) + pow(TrcSum, 2)) < sqrt(pow(vtxSigmaMin, 2) + pow(TrcSumMin, 2)))
						{
							vtxSigmaMin = vtxSigma;
							TrcSumMin = TrcSum;

							trcsum = TrcSumMin;

							done = 1;

							for (Int_t l = 0; l < 4; l++)
							{
								Knetriangle[l] = Knerec[l];

								// if (l == 3)
								// {
								// 	if (Knetri_kinfit[9] < 0.)
								// 		Knetriangle[6 + l] = -neu_vtx[l];
								// 	else if (Knetri_kinfit[9] >= 0.)
								// 		Knetriangle[6 + l] = neu_vtx[l];
								// }
								// else
								// {
									Knetriangle[6 + l] = neu_vtx[l];
								//}

								fourg4taken[l] = ind_gam[l];

								trcfinal[l] = trc[l];

								neutral_mom(cluster[0][baseKin.ncll[ind_gam[l]] - 1], cluster[1][baseKin.ncll[ind_gam[l]] - 1], cluster[2][baseKin.ncll[ind_gam[l]] - 1], cluster[4][baseKin.ncll[ind_gam[l]] - 1], neu_vtx, gammatriangle[l]);
							}

							minv4gam = sqrt(pow(gammatriangle[0][3] + gammatriangle[1][3] + gammatriangle[2][3] + gammatriangle[3][3], 2) -
															pow(gammatriangle[0][0] + gammatriangle[1][0] + gammatriangle[2][0] + gammatriangle[3][0], 2) -
															pow(gammatriangle[0][1] + gammatriangle[1][1] + gammatriangle[2][1] + gammatriangle[3][1], 2) -
															pow(gammatriangle[0][2] + gammatriangle[1][2] + gammatriangle[2][2] + gammatriangle[3][2], 2));

							Knetriangle[4] = 0.;
							for (Int_t l = 0; l < 3; l++)
							{
								Knetriangle[4] += pow(Knetriangle[l], 2);
								ip_triangle[l] = ip[l];
							}

							Knetriangle[5] = sqrt(pow(Knetriangle[3], 2) - Knetriangle[4]);
							Knetriangle[4] = sqrt(Knetriangle[4]);
						}
					}
				}

				Double_t vKne = cVel * (Knereclor[4]/Knereclor[3]);
				Double_t trec = sqrt(pow(Knereclor[6] - ip[0], 2) + 
														 pow(Knereclor[7] - ip[1], 2) +
														 pow(Knereclor[8] - ip[2], 2)) / vKne;

				if(mctruth == 1 || mctruth == 2)
				{
					test2d[0]->Fill(Knereclor[6], Knetriangle[6]);
					test2d[1]->Fill(Knereclor[7], Knetriangle[7]);
					test2d[2]->Fill(Knereclor[8], Knetriangle[8]);
					test2d[3]->Fill(trec, Knetriangle[9]);
				};
			}

			tree->Fill();
		}

		std::vector<TCanvas*> canva;
		TString canva_name = "";
		
		for (Int_t i = 0; i < 4; i++)
		{
			canva_name = "canva_" + std::to_string(i);
			canva.push_back(new TCanvas(canva_name, canva_name, 750, 750));

			test2d[i]->Draw("COLZ");

			canva[i]->Print(canva_name + ext_img);
		};

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
