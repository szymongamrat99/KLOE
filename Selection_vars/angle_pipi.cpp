#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TError.h"

#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"
#include "../../Include/Codes/lorentz_transf.h"

#include "chain_init.C"

void selection_vars(UInt_t first_file, UInt_t last_file)
{
  TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	TString name = "selection_vars.root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_vars", "Neu vtx rec with trilateration kin fit");

	TFile *file_corr = new TFile("bunch_corr.root");
  TTree *tree_corr = (TTree *)file_corr->Get("h_bunch_corr");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_mom_err[4], bhabha_vtx[3];

	chain->SetBranchAddress("Bpx", &bhabha_mom[0]);
	chain->SetBranchAddress("Bpy", &bhabha_mom[1]);
	chain->SetBranchAddress("Bpz", &bhabha_mom[2]);
	chain->SetBranchAddress("Broots", &bhabha_mom[3]);

	chain->SetBranchAddress("Bwidpx", &bhabha_mom_err[0]);
	chain->SetBranchAddress("Bwidpy", &bhabha_mom_err[1]);
	chain->SetBranchAddress("Bwidpz", &bhabha_mom_err[2]);
	chain->SetBranchAddress("Brootserr", &bhabha_mom_err[3]);

	chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
	chain->SetBranchAddress("By", &bhabha_vtx[1]);
	chain->SetBranchAddress("Bz", &bhabha_vtx[2]);

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth, mcflag;
	Float_t cluster[5][500], Kchboost[9], Knerec[9], Knemc[9], ipmc[3], ip[3], Dtmc, bunch_corr;

	chain->SetBranchAddress("mctruth", &mctruth);
	chain->SetBranchAddress("mcflag", &mcflag);


}