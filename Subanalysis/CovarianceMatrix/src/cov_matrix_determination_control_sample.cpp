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

int CovarianceMatrixDeterminationControlSample(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	// =============================================================================
	BaseKinematics
			baseKin;
	NeutRec4
			neutVars;
	TChain
			*chainDoublePiPi;
	// =============================================================================

	chainDoublePiPi = new TChain("h1");

	chainDoublePiPi->Add(charged_dir + root_files_dir + "2025-05-17/*.root");

	Float_t minDiff;
	Char_t vtxTwoTracks;
	Int_t mcflag_int, errflagks, errflagkl;

	std::vector<Float_t>
			*trk1KL = &baseKin.trkKL[0],
			*trk2KL = &baseKin.trkKL[1],
			*trk1KLTwoBody = &baseKin.trkKLTwoBody[0],
			*trk2KLTwoBody = &baseKin.trkKLTwoBody[1],
			*KchrecKL = &baseKin.KchrecKL,
			*KchrecKLTwoBody = &baseKin.KchrecKLTwoBody;

	baseKin.trkKL[0].resize(4);
	baseKin.trkKL[1].resize(4);
	baseKin.trkKLTwoBody[0].resize(4);
	baseKin.trkKLTwoBody[1].resize(4);

	chainDoublePiPi->SetBranchAddress("mcflag", &mcflag_int);
	chainDoublePiPi->SetBranchAddress("errflagks", &errflagks);
	chainDoublePiPi->SetBranchAddress("errflagkl", &errflagkl);

	chainDoublePiPi->SetBranchAddress("minDiff", &minDiff);

	chainDoublePiPi->SetBranchAddress("trk1KL", &trk1KL);
	chainDoublePiPi->SetBranchAddress("trk2KL", &trk2KL);

	chainDoublePiPi->SetBranchAddress("trk1TwoBody", &trk1KLTwoBody);
	chainDoublePiPi->SetBranchAddress("trk2TwoBody", &trk2KLTwoBody);

	chainDoublePiPi->SetBranchAddress("KchrecKL", &KchrecKL);
	chainDoublePiPi->SetBranchAddress("KchrecKLTwoBody", &KchrecKLTwoBody);

	const Int_t numberOfMomenta = 2;

	TVectorT<Double_t>
			momVecMC(numberOfMomenta * 3),
			momVecData(numberOfMomenta * 3);

	TMatrixT<Double_t>
			covMatrix(numberOfMomenta * 3, numberOfMomenta * 3);

	KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, momVecData, covMatrix);

	const Int_t nentries = chainDoublePiPi->GetEntries();

	for (Int_t i = 0; i < nentries; i++)
	{
		chainDoublePiPi->GetEntry(i);

		if (minDiff < 0.5 && abs(KchrecKL->at(5) - mK0) < 2.0)
		{

			for (Int_t j = 0; j < 3; j++)
			{
				momVecMC[j] = trk1KLTwoBody->at(j);
				momVecData[j] = trk1KL->at(j);

				momVecMC[3 + j] = trk2KLTwoBody->at(j);
				momVecData[3 + j] = trk2KL->at(j);
			}

			CovMatrixCalcObj.SetMCVector(momVecMC);
			CovMatrixCalcObj.SetDataVector(momVecData);

			CovMatrixCalcObj.CovCalcPart();
		}
	}

	CovMatrixCalcObj.GetCovMatrix();

	return 0;
}