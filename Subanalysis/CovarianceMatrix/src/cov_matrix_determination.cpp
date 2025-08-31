#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRatioPlot.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <ConfigManager.h>

#include "../inc/covmatrix.hpp"

int CovarianceMatrixDetermination(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	// =============================================================================
	BaseKinematics
			baseKin;

	ConfigManager &config = ConfigManager::getInstance();
	// =============================================================================

	std::vector<Float_t>
			*trk1KL = &baseKin.trkKL[0],
			*trk2KL = &baseKin.trkKL[1],
			*trk1KLmc = &baseKin.trkKLmc[0],
			*trk2KLmc = &baseKin.trkKLmc[1],
			*trk1KS = &baseKin.trkKS[0],
			*trk2KS = &baseKin.trkKS[1],
			*trk1KSmc = &baseKin.trkKSmc[0],
			*trk2KSmc = &baseKin.trkKSmc[1];
	;

	baseKin.trkKL[0].resize(4);
	baseKin.trkKL[1].resize(4);
	baseKin.trkKLmc[0].resize(4);
	baseKin.trkKLmc[1].resize(4);
	baseKin.trkKS[0].resize(4);
	baseKin.trkKS[1].resize(4);
	baseKin.trkKSmc[0].resize(4);
	baseKin.trkKSmc[1].resize(4);

	Int_t mcflag = 0;

	chain.SetBranchAddress("trk1KL", &trk1KL);
	chain.SetBranchAddress("trk2KL", &trk2KL);
	chain.SetBranchAddress("trk1KS", &trk1KS);
	chain.SetBranchAddress("trk2KS", &trk2KS);

	chain.SetBranchAddress("mcflag", &mcflag);
	chain.SetBranchAddress("mctruth", &baseKin.mctruth_int);

	chain.SetBranchAddress("trk1KLmc", &trk1KLmc);
	chain.SetBranchAddress("trk2KLmc", &trk2KLmc);
	chain.SetBranchAddress("trk1KSmc", &trk1KSmc);
	chain.SetBranchAddress("trk2KSmc", &trk2KSmc);

	const Int_t numberOfMomenta = 2;

	TVectorT<Double_t>
			momVecMC(numberOfMomenta * 3),
			momVecData(numberOfMomenta * 3);

	TMatrixT<Double_t>
			covMatrix(numberOfMomenta * 3, numberOfMomenta * 3);

	KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, momVecData, covMatrix);

	const Int_t nentries = chain.GetEntries();

	Double_t checkMCKL[2] = {0.}, checkMCKS[2] = {0.};

	std::string covMatrixType = config.getProperty<std::string>("flags.covMatrixType", "KL");

	for (Int_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		checkMCKL[0] = 999.;
		checkMCKL[1] = 999.;
		checkMCKS[0] = 999.;
		checkMCKS[1] = 999.;

		if (mcflag == 1 && baseKin.mctruth_int == 7)
		{
			if (covMatrixType == "KL" || covMatrixType == "MIXED")
			{
				checkMCKL[0] = sqrt(pow(trk1KLmc->at(0) - trk1KL->at(0), 2) +
														pow(trk1KLmc->at(1) - trk1KL->at(1), 2) +
														pow(trk1KLmc->at(2) - trk1KL->at(2), 2) +
														pow(trk2KLmc->at(0) - trk2KL->at(0), 2) +
														pow(trk2KLmc->at(1) - trk2KL->at(1), 2) +
														pow(trk2KLmc->at(2) - trk2KL->at(2), 2));

				checkMCKL[1] = sqrt(pow(trk1KLmc->at(0) - trk2KL->at(0), 2) +
														pow(trk1KLmc->at(1) - trk2KL->at(1), 2) +
														pow(trk1KLmc->at(2) - trk2KL->at(2), 2) +
														pow(trk2KLmc->at(0) - trk1KL->at(0), 2) +
														pow(trk2KLmc->at(1) - trk1KL->at(1), 2) +
														pow(trk2KLmc->at(2) - trk1KL->at(2), 2));

				for (Int_t j = 0; j < 3; j++)
				{
					if (checkMCKL[0] < checkMCKL[1])
					{
						momVecMC[j] = trk1KLmc->at(j);
						momVecData[j] = trk1KL->at(j);

						momVecMC[3 + j] = trk2KLmc->at(j);
						momVecData[3 + j] = trk2KL->at(j);
					}
					else
					{
						momVecMC[j] = trk1KLmc->at(j);
						momVecData[j] = trk2KL->at(j);

						momVecMC[3 + j] = trk2KLmc->at(j);
						momVecData[3 + j] = trk1KL->at(j);
					}
				}

				CovMatrixCalcObj.AddCovariancePoint(momVecMC, momVecData);
			}

			if (covMatrixType == "KS" || covMatrixType == "MIXED")
			{
				checkMCKS[0] = sqrt(pow(trk1KSmc->at(0) - trk1KS->at(0), 2) +
														pow(trk1KSmc->at(1) - trk1KS->at(1), 2) +
														pow(trk1KSmc->at(2) - trk1KS->at(2), 2) +
														pow(trk2KSmc->at(0) - trk2KS->at(0), 2) +
														pow(trk2KSmc->at(1) - trk2KS->at(1), 2) +
														pow(trk2KSmc->at(2) - trk2KS->at(2), 2));

				checkMCKS[1] = sqrt(pow(trk1KSmc->at(0) - trk2KS->at(0), 2) +
														pow(trk1KSmc->at(1) - trk2KS->at(1), 2) +
														pow(trk1KSmc->at(2) - trk2KS->at(2), 2) +
														pow(trk2KSmc->at(0) - trk1KS->at(0), 2) +
														pow(trk2KSmc->at(1) - trk1KS->at(1), 2) +
														pow(trk2KSmc->at(2) - trk1KS->at(2), 2));

				for (Int_t j = 0; j < 3; j++)
				{
					if (checkMCKS[0] < checkMCKS[1])
					{
						momVecMC[j] = trk1KSmc->at(j);
						momVecData[j] = trk1KS->at(j);

						momVecMC[3 + j] = trk2KSmc->at(j);
						momVecData[3 + j] = trk2KS->at(j);
					}
					else
					{
						momVecMC[j] = trk1KSmc->at(j);
						momVecData[j] = trk2KS->at(j);

						momVecMC[3 + j] = trk2KSmc->at(j);
						momVecData[3 + j] = trk1KS->at(j);
					}
				}

				CovMatrixCalcObj.AddCovariancePoint(momVecMC, momVecData);
			}
		}
	}

	CovMatrixCalcObj.GetCovMatrix();

	CovMatrixCalcObj.SaveCovMatrixToJSON("covarianceMatrixMCGenerated" + covMatrixType);

	// --- Phase 2: Uncertainty Calculation using Bootstrap ---
	const int num_bootstrap_samples = 100;
	std::cout << "\n--- Calculating Uncertainty with Bootstrap (" << num_bootstrap_samples << " samples) ---" << std::endl;
	CovMatrixCalcObj.CovMatrixUncertainty(num_bootstrap_samples);

	return 0;
}