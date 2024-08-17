#include <iostream>

#include <TTree.h>
#include <TFile.h>

#include "../../Include/Codes/chain_init.cpp"
#include "../../Include/Codes/cylinder_intersection.h"
#include "../../Include/const.h"

Int_t main(int argc, char *argv[])
{
	UInt_t first = atoi(argv[1]), last = atoi(argv[2]), nclus = atoi(argv[3]);

	TString first_s = argv[1], last_s = argv[2];

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first, last);

	TString name = "gen_vars_" + first_s + "_" + last_s + ".root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_gen_vars", "Gen vars for klspm00");

	// Branches' addresses
	// Bhabha vars
	Int_t ntmc, nvtxmc, nclu;
	UChar_t pidmc[200], vtxmc[200], mother[200], mctruth = 0, mcflag = 0;
	Float_t pos_mc[3][200], mom_mc[3][200], Knemc[9], cluster_rec[3][200], ipmc[3];

	chain->SetBranchAddress("ntmc", &ntmc);
	chain->SetBranchAddress("nvtxmc", &nvtxmc);

	chain->SetBranchAddress("pidmc", pidmc);
	chain->SetBranchAddress("vtxmc", vtxmc);
	chain->SetBranchAddress("mother", mother);

	chain->SetBranchAddress("xvmc", pos_mc[0]);
	chain->SetBranchAddress("yvmc", pos_mc[1]);
	chain->SetBranchAddress("zvmc", pos_mc[2]);

	chain->SetBranchAddress("pxmc", mom_mc[0]);
	chain->SetBranchAddress("pymc", mom_mc[1]);
	chain->SetBranchAddress("pzmc", mom_mc[2]);

	chain->SetBranchAddress("Knemc", Knemc);

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster_rec[0]);
	chain->SetBranchAddress("Ycl", cluster_rec[1]);
	chain->SetBranchAddress("Zcl", cluster_rec[2]);

	chain->SetBranchAddress("mctruth", &mctruth);
	chain->SetBranchAddress("mcflag", &mcflag);
	chain->SetBranchAddress("ipmc", ipmc);

	Int_t nentries = (Int_t)chain->GetEntries();

	Float_t pgammc[4][8], neu_vtx[3], cluster[3],
			Knemcnew[9], Kchmcnew[9], ipmcnew[3], Ks[9], Kl[9];
	Int_t good_clus_ind[4], region[4];

	TBranch *b_pgammc1 = tree->Branch("pgammc1", pgammc[0], "pgammc1[8]/F");
	TBranch *b_pgammc2 = tree->Branch("pgammc2", pgammc[1], "pgammc2[8]/F");
	TBranch *b_pgammc3 = tree->Branch("pgammc3", pgammc[2], "pgammc3[8]/F");
	TBranch *b_pgammc4 = tree->Branch("pgammc4", pgammc[3], "pgammc4[8]/F");

	TBranch *b_Knemc = tree->Branch("Knemcnew", Knemcnew, "Knemcnew[9]/F");
	TBranch *b_Kchmc = tree->Branch("Kchmcnew", Kchmcnew, "Kchmcnew[9]/F");
	TBranch *b_ipmc = tree->Branch("ipmcnew", ipmcnew, "ipmcnew[3]/F");

	TBranch *b_clusindgood = tree->Branch("clusindgood", good_clus_ind, "clusindgood[4]/I");

	TBranch *b_region = tree->Branch("region", region, "region[4]/I");

	const Int_t max_count = TMath::Factorial(nclus);
	Int_t count = 0, ind_gam[4], mc_ind[4] = {0, 1, 2, 3}, min_ind[max_count];
	Float_t clus_diff[max_count], clus_diff_min;

	Bool_t clus_time[max_count];

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		clus_diff_min = 999999.;

		good_clus_ind[0] = 999;
		good_clus_ind[1] = 999;
		good_clus_ind[2] = 999;
		good_clus_ind[3] = 999;

		if (mctruth == 1 || mctruth == 2)
		{
			for (Int_t j = 0; j < nvtxmc; j++)
			{
				if (mother[j] == 50)
				{
					ipmcnew[0] = pos_mc[0][j];
					ipmcnew[1] = pos_mc[1][j];
					ipmcnew[2] = pos_mc[2][j];
				}
			}

			for (Int_t j = 0; j < nvtxmc; j++)
			{
				if (mother[j] == 10)
				{
					Kl[6] = pos_mc[0][j];
					Kl[7] = pos_mc[1][j];
					Kl[8] = pos_mc[2][j];
				}
			}

			for (Int_t j = 0; j < nvtxmc; j++)
			{
				if (mother[j] == 16)
				{
					Ks[6] = pos_mc[0][j];
					Ks[7] = pos_mc[1][j];
					Ks[8] = pos_mc[2][j];
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if (pidmc[j] == 10)
				{
					Kl[0] = mom_mc[0][j];
					Kl[1] = mom_mc[1][j];
					Kl[2] = mom_mc[2][j];
					Kl[5] = mK0;
					Kl[4] = pow(Kl[0], 2) + pow(Kl[1], 2) + pow(Kl[2], 2);
					Kl[3] = sqrt(Kl[4] + pow(Kl[5], 2));
					Kl[4] = sqrt(Kl[4]);
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if (pidmc[j] == 16)
				{
					Ks[0] = mom_mc[0][j];
					Ks[1] = mom_mc[1][j];
					Ks[2] = mom_mc[2][j];
					Ks[5] = mK0;
					Ks[4] = pow(Ks[0], 2) + pow(Ks[1], 2) + pow(Ks[2], 2);
					Ks[3] = sqrt(Ks[4] + pow(Ks[5], 2));
					Ks[4] = sqrt(Ks[4]);
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if (mother[vtxmc[j] - 1] == 10)
				{
					if (pidmc[j] == 7)
					{
						for (Int_t k = 0; k < 9; k++)
						{
							Knemcnew[k] = Kl[k];
							Kchmcnew[k] = Ks[k];
						}
					}
					else if (pidmc[j] == 8 || pidmc[j] == 9)
					{
						for (Int_t k = 0; k < 9; k++)
						{
							Kchmcnew[k] = Kl[k];
							Knemcnew[k] = Ks[k];
						}
					}
				}
			}

			for (Int_t j = 0; j < ntmc; j++)
			{
				if ((mother[vtxmc[j] - 1] == 7) && pidmc[j] == 1)
				{
					pgammc[count][0] = mom_mc[0][j];
					pgammc[count][1] = mom_mc[1][j];
					pgammc[count][2] = mom_mc[2][j];
					pgammc[count][3] = sqrt(pow(pgammc[count][0], 2) +
																	pow(pgammc[count][1], 2) +
																	pow(pgammc[count][2], 2));

					neu_vtx[0] = Knemc[6];
					neu_vtx[1] = Knemc[7];
					neu_vtx[2] = Knemc[8];

					region[count] = inter_point(pgammc[count], neu_vtx, cluster);

					pgammc[count][4] = cluster[0];
					pgammc[count][5] = cluster[1];
					pgammc[count][6] = cluster[2];

					Float_t beta_c = cVel * Knemc[4] / Knemc[3], length = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2) + pow(Knemc[8] - ipmc[2], 2)), time_K = length / beta_c;

					Float_t length_clus = sqrt(pow(cluster[0] - Knemc[6], 2) + pow(cluster[1] - Knemc[7], 2) + pow(cluster[2] - Knemc[8], 2));

					pgammc[count][7] = time_K + (length_clus / cVel);

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

							if (clus_diff_min > clus_diff[min_ind[0]]) // && clus_diff[min_ind[0]] < 10)
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
	}

	tree->Print();

	file->Write();
	file->Close();
	delete file;

	return 0;
}
