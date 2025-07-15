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
#include <neu_triangle.h>

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
	chain->SetBranchAddress("TclOld", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);
	chain->SetBranchAddress("mcflag", &mcflag);
	chain->SetBranchAddress("ncll", baseKin.ncll);

	// Charged tracks momenta

	chain->SetBranchAddress("trk1", baseKin.trk[0]);
	chain->SetBranchAddress("trk2", baseKin.trk[1]);
	chain->SetBranchAddress("Kchrec", baseKin.Kchrec);
	chain->SetBranchAddress("Kchboost", baseKin.Kchboost);

	Int_t nentries = (Int_t)chain->GetEntries();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Int_t
			doneOmega,
			g4takenomega[4];

	Float_t
			gammaomega[4][8],
			Omegapi0[6],
			Pi0[6],
			PichFourMom[2][4],
			Omegarec[6],
			lengthPhotonMin[4],
			lengthKch,
			rho_00,
			rho_pm_IP,
			rho_00_IP,
			rho,
			anglePi0KaonCM,
			anglePichKaonCM,
			anglePi0OmegaPhiCM,
			anglePhiOmega,
			neu_vtx_avg[3];

	TBranch *b_gamma1omega = tree->Branch("gamma1omega", gammaomega[0], "gamma1omega[8]/F");
	TBranch *b_gamma2omega = tree->Branch("gamma2omega", gammaomega[1], "gamma2omega[8]/F");
	TBranch *b_gamma3omega = tree->Branch("gamma3omega", gammaomega[2], "gamma1omega[8]/F");
	TBranch *b_gamma4omega = tree->Branch("gamma4omega", gammaomega[3], "gamma2omega[8]/F");

	TBranch *b_pich1 = tree->Branch("pich1", PichFourMom[0], "pich1[4]/F");
	TBranch *b_pich2 = tree->Branch("pich2", PichFourMom[1], "pich2[4]/F");

	TBranch *b_pi0omega = tree->Branch("omegapi0", Omegapi0, "omegapi0[6]/F");
	TBranch *b_pi0 = tree->Branch("pi0", Pi0, "pi0[6]/F");

	TBranch *b_omega = tree->Branch("omega", Omegarec, "omega[6]/F");

	TBranch *b_neuvtxavg = tree->Branch("NeuVtxAvg", neu_vtx_avg, "NeuVtxAvg[3]/F");

	TBranch *b_lengthphoton = tree->Branch("lengthphoton", lengthPhotonMin, "lengthphoton[4]/F");

	TBranch *b_lengthkch = tree->Branch("lengthKch", &lengthKch, "lengthKch/F");

	TBranch *b_done = tree->Branch("doneomega", &doneOmega, "doneomega/I");
	TBranch *b_g4takenomega = tree->Branch("g4takenomega", g4takenomega, "g4takenomega[4]/I");

	TBranch
			*b1 = tree->Branch("rho_00", &rho_00, "rho_00/F"),
			*b2 = tree->Branch("rho_00_IP", &rho_00_IP, "rho_00_IP/F"),
			*b3 = tree->Branch("rho_pm_IP", &rho_pm_IP, "rho_pm_IP/F"),
			*b4 = tree->Branch("rho", &rho, "rho/F"),
			*b5 = tree->Branch("anglePi0KaonCM", &anglePi0KaonCM, "anglePi0KaonCM/F"),
			*b6 = tree->Branch("anglePichKaonCM", &anglePichKaonCM, "anglePichKaonCM/F"),
			*b7 = tree->Branch("anglePi0OmegaPhiCM", &anglePi0OmegaPhiCM, "anglePi0OmegaPhiCM/F"),
			*b8 = tree->Branch("anglePhiOmega", &anglePhiOmega, "anglePhiOmega/F");

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

		Float_t
				lengthPhotonMinAvg = 1E6;

		// Clear of the variables
		for (Int_t l2 = 0; l2 < 4; l2++)
		{
			Clear1DArray(8, gammaomega[l2]);
		}

		Clear1DArray(6, Omegapi0);
		Clear1DArray(6, Pi0);
		Clear1DArray(6, Omegarec);

		for (Int_t j = 0; j < 3; j++)
			neu_vtx_avg[j] = 0.;

		rho_00 = 999.;
		rho_00_IP = 999.;
		rho_pm_IP = 999.;
		rho = 999.;
		anglePhiOmega = 999.;
		anglePi0KaonCM = 999.;
		anglePi0OmegaPhiCM = 999.;
		anglePichKaonCM = 999.;

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

							Bool_t cond_ene = cluster[4][baseKin.ncll[ind_gam[0]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[1]] - 1] > MIN_CLU_ENE &&
																cluster[4][baseKin.ncll[ind_gam[2]] - 1] > MIN_CLU_ENE && cluster[4][baseKin.ncll[ind_gam[3]] - 1] > MIN_CLU_ENE;

							Bool_t cond_clus[4];

							cond_clus[0] = cluster[3][baseKin.ncll[ind_gam[0]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[0]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[0]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[0]] - 1] != 0;
							cond_clus[1] = cluster[3][baseKin.ncll[ind_gam[1]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[1]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[1]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[1]] - 1] != 0;
							cond_clus[2] = cluster[3][baseKin.ncll[ind_gam[2]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[2]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[2]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[2]] - 1] != 0;
							cond_clus[3] = cluster[3][baseKin.ncll[ind_gam[3]] - 1] > 0 && cluster[0][baseKin.ncll[ind_gam[3]] - 1] != 0 && cluster[1][baseKin.ncll[ind_gam[3]] - 1] != 0 && cluster[2][baseKin.ncll[ind_gam[3]] - 1] != 0;

							Float_t
									lengthPhotonAvg = 0.,
									totEnergy = 0.;

							for (Int_t k = 0; k < 4; k++)
							{
								lengthPhoton[k] = cluster[3][baseKin.ncll[ind_gam[k]] - 1] -
																	(sqrt(pow(cluster[0][baseKin.ncll[ind_gam[k]] - 1] - bhabha_vtx[0], 2) +
																				pow(cluster[1][baseKin.ncll[ind_gam[k]] - 1] - bhabha_vtx[1], 2) +
																				pow(cluster[2][baseKin.ncll[ind_gam[k]] - 1] - baseKin.Kchrec[8], 2)) /
																	 cVel);

								totEnergy += cluster[4][baseKin.ncll[ind_gam[k]] - 1];

								lengthPhotonAvg += pow(cluster[4][baseKin.ncll[ind_gam[k]] - 1] * lengthPhoton[k],2);
							}

							lengthPhotonAvg = sqrt(lengthPhotonAvg) / totEnergy;

							if (abs(lengthPhotonAvg) < abs(lengthPhotonMinAvg) && cond_ene && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
							{
								lengthPhotonMinAvg = lengthPhotonAvg;

								lengthPhotonMin[0] = lengthPhoton[0];
								lengthPhotonMin[1] = lengthPhoton[1];
								lengthPhotonMin[2] = lengthPhoton[2];
								lengthPhotonMin[3] = lengthPhoton[3];

								g4takenomega[0] = ind_gam[0];
								g4takenomega[1] = ind_gam[1];
								g4takenomega[2] = ind_gam[2];
								g4takenomega[3] = ind_gam[3];
							}
						}

			Float_t
					totEnergy = 0,
					photonVel[3] = {0.},
					kaonMom[3] = {0.},
					photonLength = 0.,
					IP_vtx[3] = {bhabha_vtx[0], bhabha_vtx[1], baseKin.Kchrec[8]};

			for (Int_t j = 0; j < 4; j++)
			{
				neutral_mom(cluster[0][baseKin.ncll[g4takenomega[j]] - 1], cluster[1][baseKin.ncll[g4takenomega[j]] - 1], cluster[2][baseKin.ncll[g4takenomega[j]] - 1], cluster[4][baseKin.ncll[g4takenomega[j]] - 1], IP_vtx, gammaomega[j]);

				for (Int_t k = 0; k < 4; k++)
				{
					gammaomega[j][k + 4] = cluster[k][baseKin.ncll[g4takenomega[j]] - 1];

					if (k < 3)
					{
						kaonMom[k] += gammaomega[j][k];
					}
				};
			}

			kaonMom[0] = bhabha_mom[0] - baseKin.Kchrec[0];
			kaonMom[1] = bhabha_mom[1] - baseKin.Kchrec[1];
			kaonMom[2] = bhabha_mom[2] - baseKin.Kchrec[2];
			kaonMom[3] = bhabha_mom[3] - baseKin.Kchrec[3];



			Float_t 
						kaonMomTot = sqrt(pow(kaonMom[0],2) + pow(kaonMom[1],2) + pow(kaonMom[2],2)),
						kaonVelTot = cVel *(kaonMomTot / kaonMom[3]);

			neu_vtx_avg[0] = lengthPhotonMinAvg * kaonVelTot * (kaonMom[0]/kaonMomTot) + IP_vtx[0];
			neu_vtx_avg[1] = lengthPhotonMinAvg * kaonVelTot * (kaonMom[1]/kaonMomTot) + IP_vtx[1];
			neu_vtx_avg[2] = lengthPhotonMinAvg * kaonVelTot * (kaonMom[2]/kaonMomTot) + IP_vtx[2];

			Int_t
					g4takenPi0[2][2];

			Float_t
					PhotonMomPi0[2][2][4],
					Pi0Mom[2][4],
					OmegaMom[4],
					M_omega_tmp[2],
					M_omega_diff[2],
					Pi0NonOmega[4];

			Pi0PhotonPair(g4takenomega, gammaomega, g4takenPi0, PhotonMomPi0, Pi0Mom, false, PichFourMom, OmegaMom);

			for (Int_t j = 0; j < 4; j++)
			{
				PichFourMom[0][j] = baseKin.trk[0][j];
				PichFourMom[1][j] = baseKin.trk[1][j];
			}

			// Calculation of pi+pi-pi0 invariant mass - best option chosen by quadrature.
			for (Int_t j = 0; j < 2; j++)
			{
				M_omega_tmp[j] = sqrt(pow(PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[j][3], 2) -
															pow(PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[j][0], 2) -
															pow(PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[j][1], 2) -
															pow(PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[j][2], 2));
				M_omega_diff[j] = M_omega_tmp[j] - mOmega;
			}

			//
			if (std::isnan(M_omega_diff[0]) || std::isnan(M_omega_diff[1]))
				Omegarec[5] = 999999.;
			else if (std::isinf(M_omega_diff[0]) || std::isinf(M_omega_diff[1]))
				Omegarec[5] = 999999.;
			else
			{
				doneOmega = 1;
				if (abs(M_omega_diff[0]) < abs(M_omega_diff[1]))
				{
					Omegarec[0] = PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[0][0];
					Omegarec[1] = PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[0][1];
					Omegarec[2] = PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[0][2];
					Omegarec[3] = PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[0][3];
					Omegarec[4] = sqrt(pow(Omegarec[0], 2) + pow(Omegarec[1], 2) + pow(Omegarec[2], 2));
					Omegarec[5] = M_omega_tmp[0];

					Omegapi0[0] = Pi0Mom[0][0];
					Omegapi0[1] = Pi0Mom[0][1];
					Omegapi0[2] = Pi0Mom[0][2];
					Omegapi0[3] = Pi0Mom[0][3];
					Omegapi0[4] = sqrt(pow(Pi0Mom[0][0], 2) + pow(Pi0Mom[0][1], 2) + pow(Pi0Mom[0][2], 2));
					Omegapi0[5] = sqrt(pow(Omegapi0[3], 2) - pow(Omegapi0[4], 2));

					Pi0[0] = Pi0Mom[1][0];
					Pi0[1] = Pi0Mom[1][1];
					Pi0[2] = Pi0Mom[1][2];
					Pi0[3] = Pi0Mom[1][3];
					Pi0[4] = sqrt(pow(Pi0Mom[1][0], 2) + pow(Pi0Mom[1][1], 2) + pow(Pi0Mom[1][2], 2));
					Pi0[5] = sqrt(pow(Pi0[3], 2) - pow(Pi0[4], 2));
				}
				else
				{
					Omegarec[0] = PichFourMom[0][0] + PichFourMom[1][0] + Pi0Mom[1][0];
					Omegarec[1] = PichFourMom[0][1] + PichFourMom[1][1] + Pi0Mom[1][1];
					Omegarec[2] = PichFourMom[0][2] + PichFourMom[1][2] + Pi0Mom[1][2];
					Omegarec[3] = PichFourMom[0][3] + PichFourMom[1][3] + Pi0Mom[1][3];
					Omegarec[4] = sqrt(pow(Omegarec[0], 2) + pow(Omegarec[1], 2) + pow(Omegarec[2], 2));
					Omegarec[5] = M_omega_tmp[1];

					Omegapi0[0] = Pi0Mom[1][0];
					Omegapi0[1] = Pi0Mom[1][1];
					Omegapi0[2] = Pi0Mom[1][2];
					Omegapi0[3] = Pi0Mom[1][3];
					Omegapi0[4] = sqrt(pow(Pi0Mom[1][0], 2) + pow(Pi0Mom[1][1], 2) + pow(Pi0Mom[1][2], 2));
					Omegapi0[5] = sqrt(pow(Omegapi0[3], 2) - pow(Omegapi0[4], 2));

					Pi0[0] = Pi0Mom[0][0];
					Pi0[1] = Pi0Mom[0][1];
					Pi0[2] = Pi0Mom[0][2];
					Pi0[3] = Pi0Mom[0][3];
					Pi0[4] = sqrt(pow(Pi0Mom[0][0], 2) + pow(Pi0Mom[0][1], 2) + pow(Pi0Mom[0][2], 2));
					Pi0[5] = sqrt(pow(Pi0[3], 2) - pow(Pi0[4], 2));
				}

				Float_t
						Kne[4] = {Pi0[0] + Omegapi0[0], Pi0[1] + Omegapi0[1], Pi0[2] + Omegapi0[2], Pi0[3] + Omegapi0[3]},
						boost_vec_Kchboost[3] = {-(bhabha_mom[0] - Kne[0]) / (bhabha_mom[3] - Kne[3]),
																		 -(bhabha_mom[1] - Kne[1]) / (bhabha_mom[3] - Kne[3]),
																		 -(bhabha_mom[2] - Kne[2]) / (bhabha_mom[3] - Kne[3])},
						PichFourMomKaonCM[2][4];

				lorentz_transf(boost_vec_Kchboost, PichFourMom[0], PichFourMomKaonCM[0]);
				lorentz_transf(boost_vec_Kchboost, PichFourMom[1], PichFourMomKaonCM[1]);

				TVector3
						pich1(PichFourMomKaonCM[0][0], PichFourMomKaonCM[0][1], PichFourMomKaonCM[0][2]),
						pich2(PichFourMomKaonCM[1][0], PichFourMomKaonCM[1][1], PichFourMomKaonCM[1][2]);

				anglePichKaonCM = pich1.Angle(pich2) * 180. / M_PI;

				// Lorentz transformation of Pi0 to Kaon CM frame
				Float_t
						boost_vec_Kne[3] = {-(Kne[0]) / (Kne[3]),
																-(Kne[1]) / (Kne[3]),
																-(Kne[2]) / (Kne[3])},
						Pi0KaonCM[2][4];

				lorentz_transf(boost_vec_Kne, Pi0Mom[0], Pi0KaonCM[0]);
				lorentz_transf(boost_vec_Kne, Pi0Mom[1], Pi0KaonCM[1]);

				TVector3
						pi01(Pi0KaonCM[0][0], Pi0KaonCM[0][1], Pi0KaonCM[0][2]),
						pi02(Pi0KaonCM[1][0], Pi0KaonCM[1][1], Pi0KaonCM[1][2]);

				anglePi0KaonCM = pi01.Angle(pi02) * 180. / M_PI;

				// Lorentz transformation of Pi0 to Kaon CM frame
				Float_t
						boost_vec_phi[3] = {-(bhabha_mom[0]) / (bhabha_mom[3]),
																-(bhabha_mom[1]) / (bhabha_mom[3]),
																-(bhabha_mom[2]) / (bhabha_mom[3])},
						Pi0NonOmegaCM[4],
						OmegaMomCM[4];

				lorentz_transf(boost_vec_phi, Pi0NonOmega, Pi0NonOmegaCM);
				lorentz_transf(boost_vec_phi, OmegaMom, OmegaMomCM);

				TVector3
						pi0CM(Pi0NonOmegaCM[0], Pi0NonOmegaCM[1], Pi0NonOmegaCM[2]),
						omegaCM(OmegaMomCM[0], OmegaMomCM[1], OmegaMomCM[2]);

				anglePi0OmegaPhiCM = pi0CM.Angle(omegaCM) * 180. / M_PI;

				// Angle between phi and omega

				TVector3
						phi(bhabha_mom[0], bhabha_mom[1], bhabha_mom[2]),
						omega(OmegaMom[0], OmegaMom[1], OmegaMom[2]);

				anglePhiOmega = phi.Angle(omega) * 180. / M_PI;

				rho_00 = 0.;
				rho_00_IP = sqrt(pow(neu_vtx_avg[0] - IP_vtx[0], 2) + pow(neu_vtx_avg[1] - IP_vtx[1], 2));
				rho_pm_IP = sqrt(pow(baseKin.Kchrec[6] - IP_vtx[0], 2) + pow(baseKin.Kchrec[7] - IP_vtx[1], 2));
				rho = sqrt(pow(rho_00_IP, 2) + pow(rho_pm_IP, 2));
			}
		}

		tree->Fill();
	}

	tree->Print();

	file->Write();
	file->Close();
	delete file;
	delete chain;

	return 0;
}
