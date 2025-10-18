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

int CovarianceMatrixDeterminationControlSample(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	ConfigManager &config = ConfigManager::getInstance();

	// =============================================================================
	KLOE::BaseKinematics
			baseKin;
	NeutRec4
			neutVars;
	TChain
			*chainDoublePiPi;
	// =============================================================================

	std::string 
					filepathData = config.getProperty<std::string>("variables.tree.filepath.Data.fourPiTwoBody", ""),
					filepathMC = config.getProperty<std::string>("variables.tree.filepath.MC.fourPiTwoBody", "");

	chainDoublePiPi = new TChain("h1");

	if (data_type == Controls::DataType::MC_ONLY)
		chainDoublePiPi->Add(filepathMC.c_str());
	else if (data_type == Controls::DataType::DATA_ONLY)
		chainDoublePiPi->Add(filepathData.c_str());

	Float_t gammaKS, gammaKL;
	Char_t vtxTwoTracks;
	Int_t mcflag_int, mctruth_int;

	std::vector<Float_t>
			*trk1KL = &baseKin.trkKL[0],
			*trk2KL = &baseKin.trkKL[1],
			*trk1KS = &baseKin.trkKS[0],
			*trk2KS = &baseKin.trkKS[1],
			*trk1KLTwoBody = &baseKin.trkKLTwoBody[0],
			*trk2KLTwoBody = &baseKin.trkKLTwoBody[1],
			*trk1KSTwoBody = &baseKin.trkKSTwoBody[0],
			*trk2KSTwoBody = &baseKin.trkKSTwoBody[1],
			*KchrecKL = &baseKin.KchrecKL,
			*KchrecKLTwoBody = &baseKin.KchrecKLTwoBody,
			*KchboostKL = &baseKin.KchboostKL,
			*KchrecKS = &baseKin.KchrecKS,
			*KchrecKSTwoBody = &baseKin.KchrecKSTwoBody,
			*KchboostKS = &baseKin.KchboostKS;

	baseKin.trkKL[0].resize(4);
	baseKin.trkKL[1].resize(4);
	baseKin.trkKLTwoBody[0].resize(4);
	baseKin.trkKLTwoBody[1].resize(4);

	baseKin.trkKS[0].resize(4);
	baseKin.trkKS[1].resize(4);
	baseKin.trkKSTwoBody[0].resize(4);
	baseKin.trkKSTwoBody[1].resize(4);

	baseKin.KchrecKL.resize(9);
	baseKin.KchrecKLTwoBody.resize(9);

	baseKin.KchrecKS.resize(9);
	baseKin.KchrecKSTwoBody.resize(9);

	chainDoublePiPi->SetBranchAddress("mcflag", &mcflag_int);
	chainDoublePiPi->SetBranchAddress("mctruth", &mctruth_int);

	chainDoublePiPi->SetBranchAddress("gammaKS", &gammaKS);
	chainDoublePiPi->SetBranchAddress("gammaKL", &gammaKL);

	chainDoublePiPi->SetBranchAddress("trk1KL", &trk1KL);
	chainDoublePiPi->SetBranchAddress("trk2KL", &trk2KL);
	chainDoublePiPi->SetBranchAddress("trk1KS", &trk1KS);
	chainDoublePiPi->SetBranchAddress("trk2KS", &trk2KS);

	chainDoublePiPi->SetBranchAddress("trk1KLTwoBody", &trk1KLTwoBody);
	chainDoublePiPi->SetBranchAddress("trk2KLTwoBody", &trk2KLTwoBody);
	chainDoublePiPi->SetBranchAddress("trk1KSTwoBody", &trk1KSTwoBody);
	chainDoublePiPi->SetBranchAddress("trk2KSTwoBody", &trk2KSTwoBody);

	chainDoublePiPi->SetBranchAddress("KchrecKL", &KchrecKL);
	chainDoublePiPi->SetBranchAddress("KchboostKL", &KchboostKL);
	chainDoublePiPi->SetBranchAddress("KchrecKLTwoBody", &KchrecKLTwoBody);

	chainDoublePiPi->SetBranchAddress("KchrecKS", &KchrecKS);
	chainDoublePiPi->SetBranchAddress("KchboostKS", &KchboostKS);
	chainDoublePiPi->SetBranchAddress("KchrecKSTwoBody", &KchrecKSTwoBody);

	const Int_t numberOfMomenta = 2;

	TVectorT<Double_t>
			momVecMC(numberOfMomenta * 3),
			momVecData(numberOfMomenta * 3);

	TMatrixT<Double_t>
			covMatrix(numberOfMomenta * 3, numberOfMomenta * 3);

	std::vector<TCanvas *> canvases_piKL[2];
	TString canvases_piKL_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			canvases_piKL_name = "canvases_piKL_" + std::to_string(i) + std::to_string(j);
			canvases_piKL[i].push_back(new TCanvas(canvases_piKL_name, canvases_piKL_name, 790, 790));
		};

	std::vector<TH1 *> piKL[2][2];
	TString piKL_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 2; j++)
			for (Int_t k = 0; k < 4; k++)
			{
				piKL_name = "piKL_" + std::to_string(i) + std::to_string(j) + std::to_string(k);

				if (k != 3)
					piKL[i][j].push_back(new TH1D(piKL_name, "", 50, -400, 400));
				else
					piKL[i][j].push_back(new TH1D(piKL_name, "", 50, 100, 350));
			};

	TCanvas *canvaDiffKL = new TCanvas("canvaDiffKL", "canvaDiffKL", 790, 790);
	TH1 *histDiffKL = new TH1D("histDiffKL", "", 100, -0.1, 0.1);

	std::vector<TCanvas *> kaonCanvaKL;
	TString kaonCanvaKL_name = "";

	for (Int_t i = 0; i < 5; i++)
	{
		kaonCanvaKL_name = "kaonCanvaKL_" + std::to_string(i);
		kaonCanvaKL.push_back(new TCanvas(kaonCanvaKL_name, kaonCanvaKL_name, 790, 790));
	};

	std::vector<TH1 *> KaonHistKL[2];
	TString KaonHistKL_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 5; j++)
		{
			KaonHistKL_name = "KaonHistKL_" + std::to_string(i) + std::to_string(j);
			if (j != 3 && j != 4)
				KaonHistKL[i].push_back(new TH1D(KaonHistKL_name, "", 50, -400.0, 400.0));
			else if (j == 3)
				KaonHistKL[i].push_back(new TH1D(KaonHistKL_name, "", 50, 490.0, 600.0));
			else
				KaonHistKL[i].push_back(new TH1D(KaonHistKL_name, "", 50, 490.0, 510.0));
		};

	std::vector<TCanvas *> canvases_piKS[2];
	TString canvases_piKS_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			canvases_piKS_name = "canvases_piKS_" + std::to_string(i) + std::to_string(j);
			canvases_piKS[i].push_back(new TCanvas(canvases_piKS_name, canvases_piKS_name, 790, 790));
		};

	std::vector<TH1 *> piKS[2][2];
	TString piKS_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 2; j++)
			for (Int_t k = 0; k < 4; k++)
			{
				piKS_name = "piKS_" + std::to_string(i) + std::to_string(j) + std::to_string(k);

				if (k != 3)
					piKS[i][j].push_back(new TH1D(piKS_name, "", 50, -400, 400));
				else
					piKS[i][j].push_back(new TH1D(piKS_name, "", 50, 100, 350));
			};

	TCanvas *canvaDiffKS = new TCanvas("canvaDiffKS", "canvaDiffKS", 790, 790);
	TH1 *histDiffKS = new TH1D("histDiffKS", "", 100, -0.1, 0.1);

	std::vector<TCanvas *> kaonCanvaKS;
	TString kaonCanvaKS_name = "";

	for (Int_t i = 0; i < 5; i++)
	{
		kaonCanvaKS_name = "kaonCanvaKS_" + std::to_string(i);
		kaonCanvaKS.push_back(new TCanvas(kaonCanvaKS_name, kaonCanvaKS_name, 790, 790));
	};

	std::vector<TH1 *> KaonHistKS[2];
	TString KaonHistKS_name = "";

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 5; j++)
		{
			KaonHistKS_name = "KaonHistKS_" + std::to_string(i) + std::to_string(j);
			if (j != 3 && j != 4)
				KaonHistKS[i].push_back(new TH1D(KaonHistKS_name, "", 50, -400.0, 400.0));
			else if (j == 3)
				KaonHistKS[i].push_back(new TH1D(KaonHistKS_name, "", 50, 490.0, 600.0));
			else
				KaonHistKS[i].push_back(new TH1D(KaonHistKS_name, "", 50, 490.0, 510.0));
		};

	KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, momVecData, covMatrix);

	const Int_t nentries = chainDoublePiPi->GetEntries();

	Int_t events[2] = {0, 0};
	Float_t efficiency[2];

	std::string covMatrixType = config.getProperty<std::string>("flags.covMatrixType", "KL");

	for (Int_t i = 0; i < nentries; i++)
	{
		chainDoublePiPi->GetEntry(i);

		histDiffKL->Fill(gammaKL);
		histDiffKS->Fill(gammaKS);

		KaonHistKL[0][0]->Fill(KchrecKL->at(0));
		KaonHistKL[0][1]->Fill(KchrecKL->at(1));
		KaonHistKL[0][2]->Fill(KchrecKL->at(2));
		KaonHistKL[0][3]->Fill(KchrecKL->at(3));
		KaonHistKL[0][4]->Fill(KchrecKL->at(5));

		KaonHistKL[1][0]->Fill(KchrecKLTwoBody->at(0));
		KaonHistKL[1][1]->Fill(KchrecKLTwoBody->at(1));
		KaonHistKL[1][2]->Fill(KchrecKLTwoBody->at(2));
		KaonHistKL[1][3]->Fill(KchrecKLTwoBody->at(3));
		KaonHistKL[1][4]->Fill(KchrecKLTwoBody->at(5));

		KaonHistKS[0][0]->Fill(KchrecKS->at(0));
		KaonHistKS[0][1]->Fill(KchrecKS->at(1));
		KaonHistKS[0][2]->Fill(KchrecKS->at(2));
		KaonHistKS[0][3]->Fill(KchrecKS->at(3));
		KaonHistKS[0][4]->Fill(KchrecKS->at(5));

		KaonHistKS[1][0]->Fill(KchrecKSTwoBody->at(0));
		KaonHistKS[1][1]->Fill(KchrecKSTwoBody->at(1));
		KaonHistKS[1][2]->Fill(KchrecKSTwoBody->at(2));
		KaonHistKS[1][3]->Fill(KchrecKSTwoBody->at(3));
		KaonHistKS[1][4]->Fill(KchrecKSTwoBody->at(5));

		events[0]++;

		events[1]++;

		if (mcflag_int == 0 || (mctruth_int != 0 && mcflag_int == 1))
		{

			if (covMatrixType == "KL" || covMatrixType == "MIXED")
			{
				for (Int_t j = 0; j < 4; j++)
				{
					piKL[0][0][j]->Fill(trk1KL->at(j));
					piKL[0][1][j]->Fill(trk2KL->at(j));

					piKL[1][0][j]->Fill(trk1KLTwoBody->at(j));
					piKL[1][1][j]->Fill(trk2KLTwoBody->at(j));

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

			if (covMatrixType == "KS" || covMatrixType == "MIXED")
			{
				for (Int_t j = 0; j < 4; j++)
				{
					piKS[0][0][j]->Fill(trk1KS->at(j));
					piKS[0][1][j]->Fill(trk2KS->at(j));

					piKS[1][0][j]->Fill(trk1KSTwoBody->at(j));
					piKS[1][1][j]->Fill(trk2KSTwoBody->at(j));

					if (j < 3)
					{
						momVecMC[j] = trk1KSTwoBody->at(j);
						momVecData[j] = trk1KS->at(j);

						momVecMC[3 + j] = trk2KSTwoBody->at(j);
						momVecData[3 + j] = trk2KS->at(j);
					}
				}
				CovMatrixCalcObj.AddCovariancePoint(momVecMC, momVecData);
			}
		}
	}

	CovMatrixCalcObj.GetCovMatrix();

	if (data_type == Controls::DataType::MC_ONLY)
		CovMatrixCalcObj.SaveCovMatrixToJSON("covarianceMatrixMC" + covMatrixType);
	else if (data_type == Controls::DataType::DATA_ONLY)
		CovMatrixCalcObj.SaveCovMatrixToJSON("covarianceMatrix" + covMatrixType);

	// --- Phase 2: Uncertainty Calculation using Bootstrap ---
	const int num_bootstrap_samples = 10000;
	std::cout << "\n--- Calculating Uncertainty with Bootstrap (" << num_bootstrap_samples << " samples) ---" << std::endl;
	CovMatrixCalcObj.CovMatrixUncertainty(num_bootstrap_samples);

	TString x_names_pi[4] = {"p^{#pi}_{x} [MeV/c]", "p^{#pi}_{y} [MeV/c]", "p^{#pi}_{z} [MeV/c]", "E^{#pi} [MeV]"};

	gStyle->SetOptStat(0);

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			TString name_canva = Form("%s_%d_%d.png", "MomentumKL", i, j);
			canvases_piKL[i][j]->cd();

			TLegend *legendPi = new TLegend(0.15, 0.8, 0.6, 0.9);
			legendPi->AddEntry(piKL[0][i][j], Form("Reconstructed charged pion %d", i), "l");
			legendPi->AddEntry(piKL[1][i][j], Form("Charged pion %d from 2-body decay", i), "l");

			piKL[0][i][j]->SetLineColor(kBlue);
			piKL[1][i][j]->SetLineColor(kRed);

			piKL[0][i][j]->GetXaxis()->SetTitle(x_names_pi[j]);
			piKL[0][i][j]->GetYaxis()->SetTitle("Counts");

			Float_t maxY;

			if (piKL[0][i][j]->GetMaximum() >= piKL[1][i][j]->GetMaximum())
				maxY = piKL[0][i][j]->GetMaximum();
			else
				maxY = piKL[1][i][j]->GetMaximum();

			piKL[0][i][j]->GetYaxis()->SetRangeUser(0.0, 1.5 * maxY);

			piKL[0][i][j]->Draw("HIST");
			piKL[1][i][j]->Draw("SAME");

			legendPi->Draw();

			canvases_piKL[i][j]->Print(name_canva);
		}

	TString x_names_kaon[5] = {"p_{x} [MeV/c]", "p_{y} [MeV/c]", "p_{z} [MeV/c]", "E [MeV]", "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]"};

	for (Int_t i = 0; i < 5; i++)
	{
		TString name_canva = Form("%s_%d.png", "KaonKL", i);
		kaonCanvaKL[i]->cd();

		TLegend *legendPi = new TLegend(0.15, 0.8, 0.6, 0.9);
		legendPi->AddEntry(KaonHistKL[0][i], "Reconstructed charged K_{2}", "l");
		legendPi->AddEntry(KaonHistKL[1][i], "Corrected K_{2} based on 2-body decay", "l");

		KaonHistKL[0][i]->SetLineColor(kBlue);
		KaonHistKL[1][i]->SetLineColor(kRed);

		KaonHistKL[0][i]->GetXaxis()->SetTitle(x_names_kaon[i]);
		KaonHistKL[0][i]->GetYaxis()->SetTitle("Counts");

		Float_t maxY;

		maxY = KaonHistKL[0][i]->GetMaximum();

		KaonHistKL[0][i]->GetYaxis()->SetRangeUser(0.0, 1.5 * maxY);

		KaonHistKL[0][i]->Draw("HIST");
		KaonHistKL[1][i]->Draw("SAME");

		legendPi->Draw();

		kaonCanvaKL[i]->Print(name_canva);
	}

	canvaDiffKL->cd();

	histDiffKL->SetLineColor(kBlue);
	histDiffKL->GetXaxis()->SetTitle("Correcting #gamma angle [^{#circ}]");
	histDiffKL->GetYaxis()->SetTitle("Counts");

	Float_t maxY = histDiffKL->GetMaximum();

	histDiffKL->GetYaxis()->SetRangeUser(0.0, 1.2 * maxY);

	histDiffKL->Draw("HIST");

	canvaDiffKL->Print("differenceKL.png");

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			TString name_canva = Form("%s_%d_%d.png", "MomentumKS", i, j);
			canvases_piKS[i][j]->cd();

			TLegend *legendPi = new TLegend(0.15, 0.8, 0.6, 0.9);
			legendPi->AddEntry(piKS[0][i][j], Form("Reconstructed charged pion %d", i), "l");
			legendPi->AddEntry(piKS[1][i][j], Form("Charged pion %d from 2-body decay", i), "l");

			piKS[0][i][j]->SetLineColor(kBlue);
			piKS[1][i][j]->SetLineColor(kRed);

			piKS[0][i][j]->GetXaxis()->SetTitle(x_names_pi[j]);
			piKS[0][i][j]->GetYaxis()->SetTitle("Counts");

			Float_t maxY;

			if (piKS[0][i][j]->GetMaximum() >= piKS[1][i][j]->GetMaximum())
				maxY = piKS[0][i][j]->GetMaximum();
			else
				maxY = piKS[1][i][j]->GetMaximum();

			piKS[0][i][j]->GetYaxis()->SetRangeUser(0.0, 1.5 * maxY);

			piKS[0][i][j]->Draw("HIST");
			piKS[1][i][j]->Draw("SAME");

			legendPi->Draw();

			canvases_piKS[i][j]->Print(name_canva);
		}

	for (Int_t i = 0; i < 5; i++)
	{
		TString name_canva = Form("%s_%d.png", "KaonKS", i);
		kaonCanvaKS[i]->cd();

		TLegend *legendPi = new TLegend(0.15, 0.8, 0.6, 0.9);
		legendPi->AddEntry(KaonHistKS[0][i], "Reconstructed charged K_{1}", "l");
		legendPi->AddEntry(KaonHistKS[1][i], "Corrected K_{1} based on 2-body decay", "l");

		KaonHistKS[0][i]->SetLineColor(kBlue);
		KaonHistKS[1][i]->SetLineColor(kRed);

		KaonHistKS[0][i]->GetXaxis()->SetTitle(x_names_kaon[i]);
		KaonHistKS[0][i]->GetYaxis()->SetTitle("Counts");

		Float_t maxY;

		maxY = KaonHistKS[0][i]->GetMaximum();

		KaonHistKS[0][i]->GetYaxis()->SetRangeUser(0.0, 1.5 * maxY);

		KaonHistKS[0][i]->Draw("HIST");
		KaonHistKS[1][i]->Draw("SAME");

		legendPi->Draw();

		kaonCanvaKS[i]->Print(name_canva);
	}

	canvaDiffKS->cd();

	histDiffKS->SetLineColor(kBlue);
	histDiffKS->GetXaxis()->SetTitle("Correcting #gamma angle [^{#circ}]");
	histDiffKS->GetYaxis()->SetTitle("Counts");

	maxY = histDiffKS->GetMaximum();

	histDiffKS->GetYaxis()->SetRangeUser(0.0, 1.2 * maxY);

	histDiffKS->Draw("HIST");

	canvaDiffKS->Print("differenceKS.png");

	return 0;
}