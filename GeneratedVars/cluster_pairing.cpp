#include <iostream>

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include "chain_init.C"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/const.h"

Int_t clus[] = {0, 1, 2, 3};

Int_t main(int argc, char *argv[])
{
	const Int_t clus_num = atoi(argv[5]);

	UInt_t first = atoi(argv[1]), last = atoi(argv[2]);

	TString first_s = argv[1], last_s = argv[2], loopcount = argv[3], constraints = argv[4];

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first, last);

	TString name = "/internal/big_one/4/users/gamrat/scripts/Scripts/Neutrec/neuvtx_tri_rec_" + first_s + "_" + last_s + ".root";

	TFile *file_kinfit = new TFile(name);
	TTree *tree_kinfit = (TTree*)file_kinfit->Get("h_tri");

	name = "gen_vars_" + first_s + "_" + last_s + ".root";

	TFile *file = new TFile(name);
	TTree *tree = (TTree*)file->Get("h_gen_vars");

	// Branches' addresses
	// Bhabha vars
	Float_t gamma_kinfit[4][8];
	Int_t done_kinfit;

	tree_kinfit->SetBranchAddress("fourgamma1tri", gamma_kinfit[0]);
  tree_kinfit->SetBranchAddress("fourgamma2tri", gamma_kinfit[1]);
  tree_kinfit->SetBranchAddress("fourgamma3tri", gamma_kinfit[2]);
  tree_kinfit->SetBranchAddress("fourgamma4tri", gamma_kinfit[3]);

	tree_kinfit->SetBranchAddress("done4", &done_kinfit);


	Float_t pgammc[4][7];

	tree->SetBranchAddress("pgammc1", pgammc[0]);
	tree->SetBranchAddress("pgammc2", pgammc[1]);
	tree->SetBranchAddress("pgammc3", pgammc[2]);
	tree->SetBranchAddress("pgammc4", pgammc[3]);

	tree_kinfit->AddFriend(tree);

	TH1 *hist1 = new TH1F("hist1", "", 100, -100, 100);
	TH1 *hist2 = new TH1F("hist2", "", 100, -100, 100);
	TH1 *hist3 = new TH1F("hist3", "", 100, -100, 100);
	TH1 *hist4 = new TH1F("hist4", "", 100, -100, 100);

	Int_t nentries = (Int_t)tree_kinfit->GetEntries();

	const Int_t max_count = TMath::Factorial(clus_num);

	Int_t sort_ind[clus_num][max_count];
	Float_t chi2pair[clus_num][max_count];

	for (Int_t i = 0; i < nentries; i++)
	{
		tree_kinfit->GetEntry(i);

		if (done_kinfit == 1)
		{
			for(Int_t j = 0; j < max_count; j++)
			{
				for(Int_t k = 0; k < clus_num; k++)
				{
					chi2pair[k][j] = ( (pgammc[k][4] -gamma_kinfit[clus[k]][4])) ;
				}

				std::next_permutation(clus, clus + clus_num);
			}

			for(Int_t j = 0; j < clus_num; j++)
				TMath::Sort(max_count, chi2pair[j], sort_ind[j], kFALSE);

			hist1->Fill(chi2pair[0][sort_ind[0][0]]);
			hist2->Fill(chi2pair[1][sort_ind[1][0]]);
			hist3->Fill(chi2pair[2][sort_ind[2][0]]);
			hist4->Fill(chi2pair[3][sort_ind[3][0]]);
		}

	}

	TCanvas *c1 = new TCanvas("canva", "", 790, 790);

	hist1->Draw();
	hist2->SetLineColor(kRed);
	hist2->Draw("SAMES");
	hist3->SetLineColor(kBlack);
	hist3->Draw("SAMES");
	hist4->SetLineColor(kGreen);
	hist4->Draw("SAMES");

	c1->Print("test.png");

	delete file;
	delete file_kinfit;

	return 0;
}
