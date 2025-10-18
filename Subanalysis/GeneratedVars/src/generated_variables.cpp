#include <cylinder_intersection.h>
#include "../inc/genvars.hpp"

Int_t GenVars(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	// Creation of filename for the analysis step
	std::string
			datestamp = Obj.getCurrentDate(),
			name = "";

	name = gen_vars_dir + root_files_dir + gen_vars_filename + datestamp + "_" + int(data_type) + ext_root;
	// -----------------------------------------------------------------------------------------

	TFile *file = new TFile(name.c_str(), "recreate");
	TTree *tree = new TTree(gen_vars_tree, "Gen vars for klspm00");

	// Branches' addresses
	// Bhabha vars
	Int_t ntmc, nvtxmc, nclu;
	UChar_t pidmcOld[200], vtxmcOld[200], motherOld[200], mctruth = 0, mcflag = 0;
	Float_t pos_mc[3][200], mom_mc[3][200], KnemcOld[9], cluster_rec[3][200], ipmcOld[3];

	chain.SetBranchAddress("ntmc", &ntmc);
	chain.SetBranchAddress("nvtxmc", &nvtxmc);

	chain.SetBranchAddress("pidmcOld", pidmcOld);
	chain.SetBranchAddress("vtxmcOld", vtxmcOld);
	chain.SetBranchAddress("motherOld", motherOld);

	chain.SetBranchAddress("xvmc", pos_mc[0]);
	chain.SetBranchAddress("yvmc", pos_mc[1]);
	chain.SetBranchAddress("zvmc", pos_mc[2]);

	chain.SetBranchAddress("pxmc", mom_mc[0]);
	chain.SetBranchAddress("pymc", mom_mc[1]);
	chain.SetBranchAddress("pzmc", mom_mc[2]);

	chain.SetBranchAddress("KnemcOld", KnemcOld);

	chain.SetBranchAddress("nclu", &nclu);
	chain.SetBranchAddress("Xcl", cluster_rec[0]);
	chain.SetBranchAddress("Ycl", cluster_rec[1]);
	chain.SetBranchAddress("Zcl", cluster_rec[2]);

	chain.SetBranchAddress("mctruth", &mctruth);
	chain.SetBranchAddress("mcflag", &mcflag);
	chain.SetBranchAddress("ipmcOld", ipmcOld);

	Int_t nentries = (Int_t)chain.GetEntries();

	Float_t pgammc[4][8], neu_vtx[3], cluster[3],
			Knemcnew[9], Kchmcnew[9], ipmcnew[3], Ks[9], Kl[9], trkMC[2][4];
	Int_t good_clus_ind[4], region[4];

	TBranch *b_pgammc1 = tree->Branch("pgammc1", pgammc[0], "pgammc1[8]/F");
	TBranch *b_pgammc2 = tree->Branch("pgammc2", pgammc[1], "pgammc2[8]/F");
	TBranch *b_pgammc3 = tree->Branch("pgammc3", pgammc[2], "pgammc3[8]/F");
	TBranch *b_pgammc4 = tree->Branch("pgammc4", pgammc[3], "pgammc4[8]/F");

	TBranch *b_trkMC1 = tree->Branch("trkMC1", trkMC[0], "trkMC1[4]/F");
	TBranch *b_trkMC2 = tree->Branch("trkMC2", trkMC[1], "trkMC2[4]/F");

	TBranch *b_Knemc = tree->Branch("Knemcnew", Knemcnew, "Knemcnew[9]/F");
	TBranch *b_Kchmc = tree->Branch("Kchmcnew", Kchmcnew, "Kchmcnew[9]/F");
	TBranch *b_ipmc = tree->Branch("ipmcnew", ipmcnew, "ipmcnew[3]/F");

	TBranch *b_clusindgood = tree->Branch("clusindgood", good_clus_ind, "clusindgood[4]/I");

	TBranch *b_region = tree->Branch("region", region, "region[4]/I");

	const Int_t max_count = TMath::Factorial(4);
	Int_t count = 0, ind_gam[4], mc_ind[4] = {0, 1, 2, 3}, min_ind[max_count];
	Float_t clus_diff[max_count], clus_diff_min;

	Bool_t clus_time[max_count];

	KLOE::CylinderIntersection CylIndObj;

	boost::progress_display show_progress(nentries);

	for (Int_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		clus_diff_min = 999999.;

		good_clus_ind[0] = 999;
		good_clus_ind[1] = 999;
		good_clus_ind[2] = 999;
		good_clus_ind[3] = 999;

		for (size_t j = 0; j < 2; j++)
			for (size_t k = 0; k < 4; k++)
			{
				trkMC[j][k] = 0;
			}

		if (mctruth == 1 || mctruth == 2 || mctruth == 3 || mctruth == 0)
		{
			for (Int_t j = 0; j < nvtxmc; j++)
			{
				if (motherOld[j] == 50)
				{
					ipmcnew[0] = pos_mc[0][j];
					ipmcnew[1] = pos_mc[1][j];
					ipmcnew[2] = pos_mc[2][j];
				}
			}

			for (Int_t j = 0; j < nvtxmc; j++)
			{
				if (motherOld[j] == 10)
				{
					Kl[6] = pos_mc[0][j];
					Kl[7] = pos_mc[1][j];
					Kl[8] = pos_mc[2][j];
				}
			}

			for (Int_t j = 0; j < nvtxmc; j++)
			{
				if (motherOld[j] == 16)
				{
					Ks[6] = pos_mc[0][j];
					Ks[7] = pos_mc[1][j];
					Ks[8] = pos_mc[2][j];
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if (pidmcOld[j] == 10)
				{
					Kl[0] = mom_mc[0][j];
					Kl[1] = mom_mc[1][j];
					Kl[2] = mom_mc[2][j];
					Kl[5] = PhysicsConstants::mK0;
					Kl[4] = pow(Kl[0], 2) + pow(Kl[1], 2) + pow(Kl[2], 2);
					Kl[3] = sqrt(Kl[4] + pow(Kl[5], 2));
					Kl[4] = sqrt(Kl[4]);
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if (pidmcOld[j] == 16)
				{
					Ks[0] = mom_mc[0][j];
					Ks[1] = mom_mc[1][j];
					Ks[2] = mom_mc[2][j];
					Ks[5] = PhysicsConstants::mK0;
					Ks[4] = pow(Ks[0], 2) + pow(Ks[1], 2) + pow(Ks[2], 2);
					Ks[3] = sqrt(Ks[4] + pow(Ks[5], 2));
					Ks[4] = sqrt(Ks[4]);
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if (motherOld[vtxmcOld[j] - 1] == 10)
				{
					if (pidmcOld[j] == 7)
					{
						for (Int_t k = 0; k < 9; k++)
						{
							Knemcnew[k] = Kl[k];
							Kchmcnew[k] = Ks[k];
						}
					}
					else if (pidmcOld[j] == 8 || pidmcOld[j] == 9)
					{
						for (Int_t k = 0; k < 9; k++)
						{
							Kchmcnew[k] = Kl[k];
							Knemcnew[k] = Ks[k];
						}
					}
				}

				if ((motherOld[vtxmcOld[j] - 1] == 10 || motherOld[vtxmcOld[j] - 1] == 16) && (pidmcOld[j] == 8 || pidmcOld[j] == 9))
				{
					if (trkMC[0][3] == 0)
					{
						trkMC[0][0] = mom_mc[0][j];
						trkMC[0][1] = mom_mc[1][j];
						trkMC[0][2] = mom_mc[2][j];
						trkMC[0][3] = sqrt(pow(trkMC[0][0], 2) +
															 pow(trkMC[0][1], 2) +
															 pow(trkMC[0][2], 2) +
															 pow(PhysicsConstants::mPiCh, 2));
					}
					else
					{
						trkMC[1][0] = mom_mc[0][j];
						trkMC[1][1] = mom_mc[1][j];
						trkMC[1][2] = mom_mc[2][j];
						trkMC[1][3] = sqrt(pow(trkMC[1][0], 2) +
															 pow(trkMC[1][1], 2) +
															 pow(trkMC[1][2], 2) +
															 pow(PhysicsConstants::mPiCh, 2));
					}
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if ((motherOld[vtxmcOld[j] - 1] == 7) && pidmcOld[j] == 1)
				{
					pgammc[count][0] = mom_mc[0][j];
					pgammc[count][1] = mom_mc[1][j];
					pgammc[count][2] = mom_mc[2][j];
					pgammc[count][3] = sqrt(pow(pgammc[count][0], 2) +
																	pow(pgammc[count][1], 2) +
																	pow(pgammc[count][2], 2));

					neu_vtx[0] = KnemcOld[6];
					neu_vtx[1] = KnemcOld[7];
					neu_vtx[2] = KnemcOld[8];

					region[count] = CylIndObj.inter_point(pgammc[count], neu_vtx, cluster);

					pgammc[count][4] = cluster[0];
					pgammc[count][5] = cluster[1];
					pgammc[count][6] = cluster[2];

					Float_t beta_c = PhysicsConstants::cVel * KnemcOld[4] / KnemcOld[3], length = sqrt(pow(KnemcOld[6] - ipmcOld[0], 2) + pow(KnemcOld[7] - ipmcOld[1], 2) + pow(KnemcOld[8] - ipmcOld[2], 2)), time_K = length / beta_c;

					Float_t length_clus = sqrt(pow(cluster[0] - KnemcOld[6], 2) + pow(cluster[1] - KnemcOld[7], 2) + pow(cluster[2] - KnemcOld[8], 2));

					pgammc[count][7] = time_K + (length_clus / PhysicsConstants::cVel);

					count++;
				}
			}

			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;

							for (Int_t k = 0; k < max_count; k++)
							{

								clus_diff[k] = sqrt(pow(cluster_rec[0][ind_gam[0]] - pgammc[mc_ind[0]][4], 2) +
																		pow(cluster_rec[1][ind_gam[0]] - pgammc[mc_ind[0]][5], 2) +
																		pow(cluster_rec[2][ind_gam[0]] - pgammc[mc_ind[0]][6], 2)) +
															 sqrt(pow(cluster_rec[0][ind_gam[1]] - pgammc[mc_ind[1]][4], 2) +
																		pow(cluster_rec[1][ind_gam[1]] - pgammc[mc_ind[1]][5], 2) +
																		pow(cluster_rec[2][ind_gam[1]] - pgammc[mc_ind[1]][6], 2)) +
															 sqrt(pow(cluster_rec[0][ind_gam[2]] - pgammc[mc_ind[2]][4], 2) +
																		pow(cluster_rec[1][ind_gam[2]] - pgammc[mc_ind[2]][5], 2) +
																		pow(cluster_rec[2][ind_gam[2]] - pgammc[mc_ind[2]][6], 2)) +
															 sqrt(pow(cluster_rec[0][ind_gam[3]] - pgammc[mc_ind[3]][4], 2) +
																		pow(cluster_rec[1][ind_gam[3]] - pgammc[mc_ind[3]][5], 2) +
																		pow(cluster_rec[2][ind_gam[3]] - pgammc[mc_ind[3]][6], 2));

								clus_time[k] = pgammc[mc_ind[0]][7] > 0. && pgammc[mc_ind[1]][7] > 0. && pgammc[mc_ind[2]][7] > 0. && pgammc[mc_ind[3]][7] > 0.;

								std::next_permutation(mc_ind, mc_ind + 4);
							}

							TMath::Sort(max_count, clus_diff, min_ind, kFALSE);

							if (clus_diff_min > clus_diff[min_ind[0]])
							{
								clus_diff_min = clus_diff[min_ind[0]];

								good_clus_ind[0] = ind_gam[0];
								good_clus_ind[1] = ind_gam[1];
								good_clus_ind[2] = ind_gam[2];
								good_clus_ind[3] = ind_gam[3];
							}
						}

			count = 0;
		}

		tree->Fill();

		++show_progress;
	}

	tree->Print();

	file->Write();
	file->Close();
	delete file;

	Utils::properties["variables"]["tree"]["filename"]["generatedvars"] = (std::string)name;
	Utils::properties["variables"]["tree"]["treename"]["generatedvars"] = (std::string)gen_vars_tree;

	Utils::properties["lastScript"] = "Generated variables for Monte Carlo.";
	Utils::properties["lastUpdate"] = Obj.getCurrentTimestamp();

	std::ofstream outfile(Paths::propName);
	outfile << Utils::properties.dump(4);
	outfile.close();

	return 0;
}
