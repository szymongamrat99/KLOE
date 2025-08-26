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

	if (data_type == Controls::DataType::MC_ONLY)
		chainDoublePiPi->Add(charged_dir + root_files_dir + "2025-08-26/*.root");
	else if (data_type == Controls::DataType::DATA_ONLY)
		chainDoublePiPi->Add(charged_dir + root_files_dir + "2025-08-26/*.root");

	Float_t gamma;
	Char_t vtxTwoTracks;
	Int_t mcflag_int;

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

	baseKin.KchrecKL.resize(9);
	baseKin.KchrecKLTwoBody.resize(9);

	chainDoublePiPi->SetBranchAddress("mcflag", &mcflag_int);

	chainDoublePiPi->SetBranchAddress("gamma", &gamma);

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

	std::vector<TCanvas *> canvases_pi[2];
	TString canvases_pi_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			canvases_pi_name = "canvases_pi_" + std::to_string(i) + std::to_string(j);
			canvases_pi[i].push_back(new TCanvas(canvases_pi_name, canvases_pi_name, 790, 790));
		};

	std::vector<TH1 *> pi[2][2];
	TString pi_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 2; j++)
			for (Int_t k = 0; k < 4; k++)
			{
				pi_name = "pi_" + std::to_string(i) + std::to_string(j) + std::to_string(k);

				if (k != 3)
					pi[i][j].push_back(new TH1D(pi_name, "", 100, -400, 400));
				else
					pi[i][j].push_back(new TH1D(pi_name, "", 100, 100, 350));
			};

	TCanvas *canvaDiff = new TCanvas("canvaDiff", "canvaDiff", 790, 790);
	TH1 *histDiff = new TH1D("histDiff", "", 100, -0.1, 0.1);

	std::vector<TCanvas *> kaonCanva;
	TString kaonCanva_name = "";

	for (Int_t i = 0; i < 5; i++)
	{
		kaonCanva_name = "kaonCanva_" + std::to_string(i);
		kaonCanva.push_back(new TCanvas(kaonCanva_name, kaonCanva_name, 790, 790));
	};

	std::vector<TH1 *> KaonHist[2];
	TString KaonHist_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 5; j++)
		{
			KaonHist_name = "KaonHist_" + std::to_string(i) + std::to_string(j);
			if (j != 3 && j != 4)
				KaonHist[i].push_back(new TH1D(KaonHist_name, "", 100, -400.0, 400.0));
			else if (j == 3)
				KaonHist[i].push_back(new TH1D(KaonHist_name, "", 100, 490.0, 600.0));
			else
				KaonHist[i].push_back(new TH1D(KaonHist_name, "", 100, 490.0, 510.0));
		};

	KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, momVecData, covMatrix);

	const Int_t nentries = chainDoublePiPi->GetEntries();

	Int_t events[2] = {0, 0};
	Float_t efficiency[2];

	for (Int_t i = 0; i < nentries; i++)
	{
		chainDoublePiPi->GetEntry(i);

		histDiff->Fill(gamma);

		if (1)
		{
			KaonHist[0][0]->Fill(KchrecKL->at(0));
			KaonHist[0][1]->Fill(KchrecKL->at(1));
			KaonHist[0][2]->Fill(KchrecKL->at(2));
			KaonHist[0][3]->Fill(KchrecKL->at(3));
			KaonHist[0][4]->Fill(KchrecKL->at(5));

			KaonHist[1][0]->Fill(trk1KLTwoBody->at(0) + trk2KLTwoBody->at(0));
			KaonHist[1][1]->Fill(trk1KLTwoBody->at(1) + trk2KLTwoBody->at(1));
			KaonHist[1][2]->Fill(trk1KLTwoBody->at(2) + trk2KLTwoBody->at(2));
			KaonHist[1][3]->Fill(trk1KLTwoBody->at(3) + trk2KLTwoBody->at(3));
			KaonHist[1][4]->Fill(sqrt(pow(trk1KLTwoBody->at(3) + trk2KLTwoBody->at(3), 2) -
																pow(trk1KLTwoBody->at(0) + trk2KLTwoBody->at(0), 2) -
																pow(trk1KLTwoBody->at(1) + trk2KLTwoBody->at(1), 2) -
																pow(trk1KLTwoBody->at(2) + trk2KLTwoBody->at(2), 2)));

			events[0]++;
		}

		if (1)
		{
			events[1]++;
			for (Int_t j = 0; j < 4; j++)
			{
				pi[0][0][j]->Fill(trk1KL->at(j));
				pi[0][1][j]->Fill(trk2KL->at(j));

				pi[1][0][j]->Fill(trk1KLTwoBody->at(j));
				pi[1][1][j]->Fill(trk2KLTwoBody->at(j));

				if (j < 3)
				{
					momVecMC[j] = trk1KLTwoBody->at(j);
					momVecData[j] = trk1KL->at(j);

					momVecMC[3 + j] = trk2KLTwoBody->at(j);
					momVecData[3 + j] = trk2KL->at(j);
				}
			}

			CovMatrixCalcObj.AddCovariancePoint(momVecMC, momVecData);
		}
	}

	CovMatrixCalcObj.GetCovMatrix();

	if (data_type == Controls::DataType::MC_ONLY)
		CovMatrixCalcObj.SaveCovMatrixToJSON("covarianceMatrixMC");
	else if (data_type == Controls::DataType::DATA_ONLY)
		CovMatrixCalcObj.SaveCovMatrixToJSON("covarianceMatrix");

	// --- Phase 2: Uncertainty Calculation using Bootstrap ---
	const int num_bootstrap_samples = 100;
	std::cout << "\n--- Calculating Uncertainty with Bootstrap (" << num_bootstrap_samples << " samples) ---" << std::endl;
	CovMatrixCalcObj.CovMatrixUncertainty(num_bootstrap_samples);

	TString x_names_pi[4] = {"p^{#pi}_{x} [MeV/c]", "p^{#pi}_{y} [MeV/c]", "p^{#pi}_{z} [MeV/c]", "E^{#pi} [MeV]"};

	gStyle->SetOptStat(0);

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			TString name_canva = Form("%s_%d_%d.png", "Momentum", i, j);
			canvases_pi[i][j]->cd();

			TLegend *legendPi = new TLegend(0.2, 0.7, 0.5, 0.9);
			legendPi->AddEntry(pi[0][i][j], Form("Reconstructed charged pion %d", i), "l");
			legendPi->AddEntry(pi[1][i][j], Form("Charged pion %d from 2-body decay", i), "l");

			pi[0][i][j]->SetLineColor(kBlue);
			pi[1][i][j]->SetLineColor(kRed);

			pi[0][i][j]->GetXaxis()->SetTitle(x_names_pi[j]);
			pi[0][i][j]->GetYaxis()->SetTitle("Counts");

			Float_t maxY;

			if (pi[0][i][j]->GetMaximum() >= pi[1][i][j]->GetMaximum())
				maxY = pi[0][i][j]->GetMaximum();
			else
				maxY = pi[1][i][j]->GetMaximum();

			pi[0][i][j]->GetYaxis()->SetRangeUser(0.0, 1.2 * maxY);

			pi[0][i][j]->Draw("HIST");
			pi[1][i][j]->Draw("SAME");

			legendPi->Draw();

			canvases_pi[i][j]->Print(name_canva);
		}

	TString x_names_kaon[5] = {"p^{K2}_{x} [MeV/c]", "p^{K2}_{y} [MeV/c]", "p^{K2}_{z} [MeV/c]", "E^{K2} [MeV]", "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]"};

	for (Int_t i = 0; i < 5; i++)
	{
		TString name_canva = Form("%s_%d.png", "Kaon", i);
		kaonCanva[i]->cd();

		TLegend *legendPi = new TLegend(0.2, 0.7, 0.5, 0.9);
		legendPi->AddEntry(KaonHist[0][i], "Reconstructed charged K^{2}", "l");
		legendPi->AddEntry(KaonHist[1][i], "Corrected K^{2} based on 2-body decay", "l");

		KaonHist[0][i]->SetLineColor(kBlue);
		KaonHist[1][i]->SetLineColor(kRed);

		KaonHist[0][i]->GetXaxis()->SetTitle(x_names_kaon[i]);
		KaonHist[0][i]->GetYaxis()->SetTitle("Counts");

		Float_t maxY;

		maxY = KaonHist[0][i]->GetMaximum();

		KaonHist[0][i]->GetYaxis()->SetRangeUser(0.0, 1.2 * maxY);

		KaonHist[0][i]->Draw("HIST");
		KaonHist[1][i]->Draw("SAME");

		legendPi->Draw();

		kaonCanva[i]->Print(name_canva);
	}

	canvaDiff->cd();

	histDiff->SetLineColor(kBlue);
	histDiff->GetXaxis()->SetTitle("Correcting #gamma angle [^{#circ}]");
	histDiff->GetYaxis()->SetTitle("Counts");

	Float_t maxY = histDiff->GetMaximum();

	histDiff->GetYaxis()->SetRangeUser(0.0, 1.2 * maxY);

	histDiff->Draw("HIST");

	canvaDiff->Print("difference.png");

	return 0;
}