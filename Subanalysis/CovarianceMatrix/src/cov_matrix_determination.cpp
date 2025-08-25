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
	NeutRec4
			neutVars;
	Int_t
			file_num;
	TFile
			*file_cs_pippim,
			*file_gen;
	TTree
			*tree_cs_pippim,
			*tree_gen;
	// =============================================================================

	TString
			filename_gen = std::string(properties["variables"]["tree"]["filename"]["generatedvars"]),
			treename_gen = std::string(properties["variables"]["tree"]["treename"]["generatedvars"]);

	try
	{

		file_gen = new TFile(filename_gen);

		if (file_gen->IsZombie())
			throw(ErrorHandling::ErrorCodes::FILE_NOT_EXIST);

		tree_gen = (TTree *)file_gen->Get(treename_gen);

		if (tree_gen->IsZombie())
			throw(ErrorHandling::ErrorCodes::TREE_NOT_EXIST);
	}
	catch (ErrorHandling::ErrorCodes err)
	{
		logger.getErrLog(err);
		return int(err);
	}

	chain.SetBranchAddress("trk1", baseKin.trk[0]);
	chain.SetBranchAddress("trk2", baseKin.trk[1]);

	chain.SetBranchAddress("mcflag", &baseKin.mcflag);
	chain.SetBranchAddress("mctruth", &baseKin.mctruth);

	Float_t trkMC[2][4];

	tree_gen->SetBranchAddress("trkMC1", trkMC[0]);
	tree_gen->SetBranchAddress("trkMC2", trkMC[1]);

	chain.AddFriend(tree_gen);

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

		checkMC[0] = sqrt(pow(trkMC[0][0] - baseKin.trk[0][0], 2) +
											pow(trkMC[0][1] - baseKin.trk[0][1], 2) +
											pow(trkMC[0][2] - baseKin.trk[0][2], 2) +
											pow(trkMC[1][0] - baseKin.trk[1][0], 2) +
											pow(trkMC[1][1] - baseKin.trk[1][1], 2) +
											pow(trkMC[1][2] - baseKin.trk[1][2], 2));

		checkMC[1] = sqrt(pow(trkMC[0][0] - baseKin.trk[1][0], 2) +
											pow(trkMC[0][1] - baseKin.trk[1][1], 2) +
											pow(trkMC[0][2] - baseKin.trk[1][2], 2) +
											pow(trkMC[1][0] - baseKin.trk[0][0], 2) +
											pow(trkMC[1][1] - baseKin.trk[0][1], 2) +
											pow(trkMC[1][2] - baseKin.trk[0][2], 2));

		for (Int_t j = 0; j < 3; j++)
		{
			if (baseKin.mcflag == 1 && baseKin.mctruth == 1)
			{
				if (checkMC[0] < checkMC[1])
				{
					momVecMC[j] = trkMC[0][j];
					momVecData[j] = baseKin.trk[0][j];

					momVecMC[3 + j] = trkMC[1][j];
					momVecData[3 + j] = baseKin.trk[1][j];
				}
				else
				{
					momVecMC[j] = trkMC[0][j];
					momVecData[j] = baseKin.trk[1][j];

					momVecMC[3 + j] = trkMC[1][j];
					momVecData[3 + j] = baseKin.trk[0][j];
				}

				CovMatrixCalcObj.SetMCVector(momVecMC);
				// CovMatrixCalcObj.SetDataVector(momVecData);

				// CovMatrixCalcObj.CovCalcPart();
			}
		}
	}

	CovMatrixCalcObj.GetCovMatrix();

	return 0;
}