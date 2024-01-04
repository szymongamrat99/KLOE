#include <iostream>

#include <TTree.h>
#include <TFile.h>

#include "chain_init.C"
#include "../../Include/Codes/line_emc_intersection.h"
#include "../../Include/const.h"

Int_t main(int argc, char *argv[])
{
	UInt_t first = atoi(argv[1]), last = atoi(argv[2]);

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
	Float_t pos_mc[3][200], mom_mc[3][200], Knemc[9], cluster_rec[3][200];

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

	Int_t nentries = (Int_t)chain->GetEntries();

	Float_t pgammc[4][7], neu_vtx[3], cluster[3];

	TBranch *b_pgammc1 = tree->Branch("pgammc1", pgammc[0], "pgammc1[7]/F");
	TBranch *b_pgammc2 = tree->Branch("pgammc2", pgammc[1], "pgammc2[7]/F");
	TBranch *b_pgammc3 = tree->Branch("pgammc3", pgammc[2], "pgammc3[7]/F");
	TBranch *b_pgammc4 = tree->Branch("pgammc4", pgammc[3], "pgammc4[7]/F");

	Int_t count = 0, ind_gam[4];
	const Int_t max_count = TMath::Factorial(4);
	Float_t clus_diff[max_count], clus_diff_min;

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		if (mctruth == 1 || mctruth == 2)
		{
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

					cluster_finder(neu_vtx, pgammc[count], cluster);

					pgammc[count][4] = cluster[0];
					pgammc[count][5] = cluster[1];
					pgammc[count][6] = cluster[2];

					count++;
				}
			}

			for(Int_t j1 = 0; j1 < nclu - 3; j1++)
				for(Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for(Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for(Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;


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
