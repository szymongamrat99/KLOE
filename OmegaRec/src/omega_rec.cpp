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

#include "const.h"
#include "uncertainties.h"
#include "charged_mom.h"
#include "neutral_mom.h"
#include "lorentz_transf.h"
#include "plane_intersection.h"
#include "closest_approach.h"
#include "chi2_dist.h"
#include <pi0_photon_pair.h>

#include "../inc/omegarec.hpp"

using namespace std;

int omegarec(Int_t first_file, Int_t last_file, Controls::DataType data_type)
{

	gErrorIgnoreLevel = 6001;

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	TString name = "";

	name = omegarec_dir + root_files_dir + omega_rec_filename + first_file + "_" + last_file + "_" + int(data_type) + ext_root;

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree(omegarec_tree, "Omega reconstruction with kin fit");

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

	BaseKinematics baseKin;

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster[0]);
	chain->SetBranchAddress("Ycl", cluster[1]);
	chain->SetBranchAddress("Zcl", cluster[2]);
	chain->SetBranchAddress("Tcl", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);
	chain->SetBranchAddress("mcflag", &mcflag);
	chain->SetBranchAddress("ncll", baseKin.ncll);

	// Charged tracks momenta

	chain->SetBranchAddress("trk1", baseKin.trk[0]);
	chain->SetBranchAddress("trk2", baseKin.trk[1]);
	chain->SetBranchAddress("Kchrec", baseKin.Kchrec);

	// tree_corr->SetBranchAddress("bunchcorr", &bunch_corr);
	// chain->AddFriend(tree_corr);

	Int_t nentries = (Int_t)chain->GetEntries();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Int_t
			doneOmega,
			g4takenomega[4];

	Float_t
			gammaomega[4][8],
			Omegapi0[10],
			Pi0[10],
			Omegarec[10],
			lengthPhotonMin[4],
			lengthKch;

	TBranch *b_gamma1omega = tree->Branch("gamma1omega", gammaomega[0], "gamma1omega[8]/F");
	TBranch *b_gamma2omega = tree->Branch("gamma2omega", gammaomega[1], "gamma2omega[8]/F");
	TBranch *b_gamma3omega = tree->Branch("gamma3omega", gammaomega[2], "gamma1omega[8]/F");
	TBranch *b_gamma4omega = tree->Branch("gamma4omega", gammaomega[3], "gamma2omega[8]/F");

	TBranch *b_pi0omega = tree->Branch("omegapi0", Omegapi0, "omegapi0[10]/F");
	TBranch *b_pi0 = tree->Branch("pi0", Pi0, "pi0[10]/F");

	TBranch *b_omega = tree->Branch("omega", Omegarec, "omega[10]/F");

	TBranch *b_lengthphoton = tree->Branch("lengthphoton", lengthPhotonMin, "lengthphoton[4]/F");

	TBranch *b_lengthkch = tree->Branch("lengthKch", &lengthKch, "lengthKch/F");

	TBranch *b_done = tree->Branch("doneomega", &doneOmega, "doneomega/I");
	TBranch *b_g4takenomega = tree->Branch("g4takenomega", g4takenomega, "g4takenomega[4]/I");

	TH1 *chi2 = new TH1F("chi2", "", 100, -10.0, 30.0);

	TH1 *pi01 = new TH1F("pi01", "", 100, 600.0, 1000.0);
	TH1 *pi02 = new TH1F("pi02", "", 100, 0.0, 200.0);

	Bool_t data_flag;
	Bool_t cond_time_clus[2];

	Float_t lengthPhoton[4];

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		doneOmega = 0;
		lengthPhotonMin[0] = 1E6;
		lengthPhotonMin[1] = 1E6;
		lengthPhotonMin[2] = 1E6;
		lengthPhotonMin[3] = 1E6;

		dataFlagSetter(data_type, data_flag, int(mcflag), int(mctruth));

		if (nclu >= 4 && data_flag)
		{
			std::cout << 100 * i / (Float_t)nentries << "% done" << std::endl;

			lengthKch = sqrt(pow(baseKin.Kchrec[6] - bhabha_vtx[0], 2) +
											 pow(baseKin.Kchrec[7] - bhabha_vtx[1], 2) +
											 pow(baseKin.Kchrec[8] - bhabha_vtx[2], 2));

			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							Int_t ind_gam[4] = {j1, j2, j3, j4};

							for (Int_t k = 0; k < 4; k++)
							{
								lengthPhoton[k] = sqrt(pow(cluster[0][baseKin.ncll[ind_gam[k]] - 1] - bhabha_vtx[0], 2) +
																			 pow(cluster[1][baseKin.ncll[ind_gam[k]] - 1] - bhabha_vtx[1], 2) +
																			 pow(cluster[2][baseKin.ncll[ind_gam[k]] - 1] - bhabha_vtx[2], 2)) -
																	(cluster[3][baseKin.ncll[ind_gam[k]] - 1] / cVel);

								if (lengthPhoton[k] < lengthPhotonMin[k])
								{
									lengthPhotonMin[k] = lengthPhoton[k];
									g4takenomega[k] = ind_gam[k];
								}
							}
						}
		}
		else
		{
			// Clear of the variables		
			for (Int_t l2 = 0; l2 < 4; l2++)
			{
				Clear1DArray(8, gammaomega[l2]);
			}

			Clear1DArray(10, Omegapi0);
			Clear1DArray(10, Pi0);
			Clear1DArray(10, Omegarec);
		}

		tree->Fill();
	}

	tree->Print();

	file->Write();
	file->Close();
	delete file;

	return 0;
}
