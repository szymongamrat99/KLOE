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
	Int_t
			file_num;
	TFile
			*file_cs_pippim,
			*file_gen;
	TTree
			*tree_cs_pippim,
			*tree_gen;
	// =============================================================================

	Float_t KchrecKS[9], KchrecKL[9];
	Char_t vtxTwoTracks;

	chain.SetBranchAddress("Kchrec", KchrecKS);
	chain.SetBranchAddress("Kchreckl", KchrecKL);

	chain.SetBranchAddress("vtx_with_two_tracks", &vtxTwoTracks);

	chain.SetBranchAddress("mcflag", &baseKin.mcflag);
	chain.SetBranchAddress("mctruth", &baseKin.mctruth);

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

		if(vtxTwoTracks == true)
			std::cout << KchrecKS[5] << " " << KchrecKL[5] << std::endl;
	}

	return 0;
}