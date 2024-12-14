#include <string.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <Math/Functor.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>

#include <reconstructor.h>
#include <const.h>
#include <uncertainties.h>
#include <charged_mom.h>
#include <neutral_mom.h>
#include <lorentz_transf.h>

void selection_vars(UInt_t firstFile, UInt_t lastFile)
{
	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, firstFile, lastFile);

	// =============================================================================
	BaseKinematics
			baseKin;
	NeutRec4
			neutVars;
	Int_t
			file_num;
	TFile
			*file,
			*fileTriangle,
			*fileMctruth,
			*fileOmega;
	TTree
			*tree,
			*treeTriangle,
			*treeMctruth,
			*treeOmega;
	// =============================================================================

	TString
			mctruthFileName = (std::string)properties["variables"]["tree"]["filename"]["mctruth"],
			mctruthTreeName = (std::string)properties["variables"]["tree"]["treename"]["mctruth"],
			fileName = (std::string)properties["variables"]["tree"]["filename"]["selectionvars"],
			treeName = (std::string)properties["variables"]["tree"]["treename"]["selectionvars"];

	// =============================================================================

	fileMctruth = new TFile(mctruthFileName);
	treeMctruth = (TTree *)fileMctruth->Get(mctruthTreeName);

	file = new TFile(fileName, "recreate");
	tree = new TTree(treeName, "Selection variables for the analysis");

	// =============================================================================

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

	TBranch
			*btCh = tree->Branch("Knetriangle");

	UInt_t nentries = chain->GetEntries();

	for (Int_t i = 0; i < nentries; i++)
	{

		tree->Fill();
	}

	file->Write();
	file->Close();
}