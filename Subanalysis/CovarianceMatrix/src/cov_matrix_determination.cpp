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

#include "../inc/covmatrix.hpp"

int CovarianceMatrixDetermination(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	// =============================================================================
	BaseKinematics
			baseKin;
	// =============================================================================

	std::vector<Float_t>
			*trk1KL = &baseKin.trkKL[0],
			*trk2KL = &baseKin.trkKL[1],
			*trk1KLmc = &baseKin.trkKLmc[0],
			*trk2KLmc = &baseKin.trkKLmc[1];

	baseKin.trkKL[0].resize(4);
	baseKin.trkKL[1].resize(4);
	baseKin.trkKLmc[0].resize(4);
	baseKin.trkKLmc[1].resize(4);

	Int_t mcflag = 0;

	chain.SetBranchAddress("trk1KS", &trk1KL);
	chain.SetBranchAddress("trk2KS", &trk2KL);

	chain.SetBranchAddress("mcflag", &mcflag);
	chain.SetBranchAddress("mctruth", &baseKin.mctruth_int);

	chain.SetBranchAddress("trk1KSmc", &trk1KLmc);
	chain.SetBranchAddress("trk2KSmc", &trk2KLmc);

	const Int_t numberOfMomenta = 2;

	TVectorT<Double_t>
			momVecMC(numberOfMomenta * 3),
			momVecData(numberOfMomenta * 3);

	TMatrixT<Double_t>
			covMatrix(numberOfMomenta * 3, numberOfMomenta * 3);

	KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, momVecData, covMatrix);

	const Int_t nentries = chain.GetEntries();

	Double_t checkMC[2] = {0.};

	for (Int_t i = 0; i < nentries; i++)
	{
		chain.GetEntry(i);

		checkMC[0] = 999.;
		checkMC[1] = 999.;

		if (mcflag == 1 && baseKin.mctruth_int == 7)
		{
			checkMC[0] = sqrt(pow(trk1KLmc->at(0) - trk1KL->at(0), 2) +
												pow(trk1KLmc->at(1) - trk1KL->at(1), 2) +
												pow(trk1KLmc->at(2) - trk1KL->at(2), 2) +
												pow(trk2KLmc->at(0) - trk2KL->at(0), 2) +
												pow(trk2KLmc->at(1) - trk2KL->at(1), 2) +
												pow(trk2KLmc->at(2) - trk2KL->at(2), 2));

			checkMC[1] = sqrt(pow(trk1KLmc->at(0) - trk2KL->at(0), 2) +
												pow(trk1KLmc->at(1) - trk2KL->at(1), 2) +
												pow(trk1KLmc->at(2) - trk2KL->at(2), 2) +
												pow(trk2KLmc->at(0) - trk1KL->at(0), 2) +
												pow(trk2KLmc->at(1) - trk1KL->at(1), 2) +
												pow(trk2KLmc->at(2) - trk1KL->at(2), 2));

			for (Int_t j = 0; j < 3; j++)
			{
				if (checkMC[0] < checkMC[1])
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
	}

	CovMatrixCalcObj.GetCovMatrix();

	CovMatrixCalcObj.SaveCovMatrixToJSON("covarianceMatrixMCGenerated");

		// --- Phase 2: Uncertainty Calculation using Bootstrap ---
	const int num_bootstrap_samples = 1000;
	std::cout << "\n--- Calculating Uncertainty with Bootstrap (" << num_bootstrap_samples << " samples) ---" << std::endl;
	CovMatrixCalcObj.CovMatrixUncertainty(num_bootstrap_samples);

	return 0;
}