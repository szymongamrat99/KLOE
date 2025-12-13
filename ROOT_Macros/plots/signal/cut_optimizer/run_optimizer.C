#include <iostream>
#include <TSystem.h>
#include <TROOT.h>
#include <TTime.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TLegend.h>

#include "src/SignalCutOptimizer.h"
#include <const.h>

TCanvas *CreateJointGraph(TString name = "default", TString title = "default", TString option = "111", std::vector<Double_t> limitValues = {}, std::vector<Double_t> significanceValues = {}, std::vector<Double_t> efficiencyValues = {}, std::vector<Double_t> purityValues = {});

void AdaptiveSampling(std::vector<Double_t> &limitValues, std::vector<Double_t> &significanceValues, std::vector<Double_t> &efficiencyValues, std::vector<Double_t> &purityValues, SignalCutOptimizer *optimizer, TChain *chain, Int_t chainEntries, Double_t minStep = 2.0, Double_t maxStep = 20.0, Double_t gradientThreshold = 0.05);

void run_optimizer(Double_t initialCutParam = 0.0, Double_t finalCutParam = 20.0, Double_t minStepSize = 2.0, Double_t maxStepSize = 20.0)
{

  // Set uniform style
  KLOE::setGlobalStyle();
  gStyle->SetPadTickY(0);

  ///////////////////////////////////////////////////
  //                 TChain setup                  //
  ///////////////////////////////////////////////////
  TChain *Schain = new TChain("h1");
  TChain *Bchain = new TChain("h1");

  for (Int_t i = 1; i <= 70; i++)
  {
    Schain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Signal_%d.root", i));
  }
  for (Int_t i = 1; i <= 24; i++)
  {
    Schain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Signal_%d.root", i));
  }
  for (Int_t i = 1; i <= 29; i++)
  {
    Schain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Signal_%d.root", i));
  }

  for (Int_t i = 1; i <= 98; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Regeneration_%d.root", i));
  }
  for (Int_t i = 1; i <= 54; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Regeneration_%d.root", i));
  }
  for (Int_t i = 1; i <= 53; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Regeneration_%d.root", i));
  }

  for (Int_t i = 1; i <= 136; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Semileptonic_%d.root", i));
  }
  for (Int_t i = 1; i <= 131; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Semileptonic_%d.root", i));
  }
  for (Int_t i = 1; i <= 78; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Semileptonic_%d.root", i));
  }

  for (Int_t i = 1; i <= 43; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_3pi0_%d.root", i));
  }
  for (Int_t i = 1; i <= 45; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_3pi0_%d.root", i));
  }
  for (Int_t i = 1; i <= 26; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_3pi0_%d.root", i));
  }

  for (Int_t i = 1; i <= 100; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Omega_%d.root", i));
  }
  for (Int_t i = 1; i <= 99; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Omega_%d.root", i));
  }
  for (Int_t i = 1; i <= 53; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Omega_%d.root", i));
  }

  for (Int_t i = 1; i <= 54; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Other_%d.root", i));
  }
  for (Int_t i = 1; i <= 55; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Other_%d.root", i));
  }
  for (Int_t i = 1; i <= 32; i++)
  {
    Bchain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Other_%d.root", i));
  }
  ///////////////////////////////////////////////////////////////////////////

  // Get the number of entries in the Signal chain
  Int_t signalEntries = Schain->GetEntries("mctruth != -1");

  // Get the number of entries in the Background chain
  Int_t backgroundEntries = Bchain->GetEntries("mctruth != -1");

  ///////////////////////////////////////////////////
  //          Analysis for cut optimization        //
  ///////////////////////////////////////////////////

  std::vector<Double_t> significanceValues,
      efficiencyValues,
      purityValues,
      limitValues;

  std::cout << "\n=== STARTING ADAPTIVE SAMPLING ===" << std::endl;

  Double_t currentLimitValue = initialCutParam;
  Double_t maxLimitValue = finalCutParam;
  Double_t currentStep = minStepSize;
  const Double_t minStep = minStepSize;
  const Double_t maxStep = maxStepSize;
  const Double_t gradientThreshold = 0.05;

  std::vector<Double_t> prevSignificance;

  while (currentLimitValue <= maxLimitValue)
  {
    Long64_t NCutSignal = 0;
    Long64_t NCutBackground = 0;

    Double_t significance = 0;
    Double_t purity = 0;
    Double_t efficiency = 0;

    Double_t limitValue = currentLimitValue;

    SignalCutOptimizer *signalOptimizer = new SignalCutOptimizer(Schain, limitValue);
    SignalCutOptimizer *bcgOptimizer = new SignalCutOptimizer(Bchain, limitValue);

    Schain->Process(signalOptimizer, "");
    Bchain->Process(bcgOptimizer, "");

    NCutSignal = signalOptimizer->GetNCut();
    NCutBackground = bcgOptimizer->GetNCut();

    significance = (NCutSignal + NCutBackground) > 0 ? NCutSignal / TMath::Sqrt(NCutSignal + NCutBackground) : 0;
    efficiency = signalEntries > 0 ? (Double_t)NCutSignal / signalEntries : 0;
    purity = (NCutSignal + NCutBackground) > 0 ? (Double_t)NCutSignal / (NCutSignal + NCutBackground) : 0;

    significanceValues.push_back(significance);
    efficiencyValues.push_back(efficiency);
    purityValues.push_back(purity);
    limitValues.push_back(limitValue);

    // Oblicz gradient dla adaptacji
    if (significanceValues.size() >= 2)
    {
      Double_t gradient = TMath::Abs(significanceValues.back() - significanceValues[significanceValues.size() - 2]) / (currentStep + 1e-6);

      if (gradient > gradientThreshold)
      {
        // Szybka zmiana - zmniejsz krok (gęstsze próbkowanie)
        currentStep = TMath::Max(minStep, currentStep * 0.75);
      }
      else
      {
        // Wolna zmiana - zwiększ krok (rzadsze próbkowanie)
        currentStep = TMath::Min(maxStep, currentStep * 1.25);
      }
    }

    std::cout << "Limit Value: " << limitValue << " (step: " << currentStep << "), S: " << significance << ", E: " << efficiency << ", P: " << purity << std::endl;

    signalOptimizer->Delete();
    bcgOptimizer->Delete();

    currentLimitValue += currentStep;
  }

  std::cout << "=== ADAPTIVE SAMPLING COMPLETED ===" << std::endl;
  std::cout << "Total points sampled: " << limitValues.size() << std::endl
            << std::endl;

  TCanvas *canva = CreateJointGraph("canva", "Cut Optimization Results (Adaptive Sampling)", "111", limitValues, significanceValues, efficiencyValues, purityValues);

  canva->SaveAs("cut_optimization_results.pdf");
}

TCanvas *CreateJointGraph(TString name = "default", TString title = "default", TString option = "111", std::vector<Double_t> limitValues = {}, std::vector<Double_t> significanceValues = {}, std::vector<Double_t> efficiencyValues = {}, std::vector<Double_t> purityValues = {})
{
  TCanvas *canva = new TCanvas("canva", "Cut Optimization Results", 1200, 800);
  canva->SetMargin(0.12, 0.25, 0.15, 0.1); // Zwiększ drugi parametr (prawy margines)

  // Create significance graph
  TGraph *significanceGraph = new TGraph(limitValues.size(), limitValues.data(), significanceValues.data());
  significanceGraph->SetTitle(";Cut Limit Value;Significance");
  significanceGraph->SetMarkerStyle(20);
  significanceGraph->SetMarkerColor(kBlack);
  significanceGraph->SetLineColor(kBlack);
  significanceGraph->SetLineWidth(2);
  significanceGraph->Draw("ALP");
  canva->Update();
  gPad->Update();

  significanceGraph->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  significanceGraph->GetYaxis()->SetTitleColor(kBlack);
  significanceGraph->GetYaxis()->SetLabelColor(kBlack);
  significanceGraph->GetYaxis()->SetTitleOffset(1.2);
  /////////////////////////////////

  if (limitValues.size() > 0)
  {
    Double_t xMin = *std::min_element(limitValues.begin(), limitValues.end());
    Double_t xMax = *std::max_element(limitValues.begin(), limitValues.end());

    // Dodaj margines 5%
    Double_t margin = (xMax - xMin) * 0.05;
    significanceGraph->GetXaxis()->SetRangeUser(xMin - margin, xMax + margin);
  }
  canva->Update();
  gPad->Update();

  // ==================== MAPOWANIE DANYCH NA PRAWĄ OŚ ====================
  // Pobierz zakresy osi Y (lewa)
  Double_t yLeftMin = significanceGraph->GetYaxis()->GetXmin();
  Double_t yLeftMax = significanceGraph->GetYaxis()->GetXmax();
  Double_t yRightMin = 0.0;
  Double_t yRightMax = 1.0;

  // Utwórz transformowane dane: mapuj (0,1) → (yLeftMin, yLeftMax)
  std::vector<Double_t> efficiencyTransformed(efficiencyValues.size());
  std::vector<Double_t> purityTransformed(purityValues.size());

  for (Int_t i = 0; i < efficiencyValues.size(); i++)
  {
    // Transformacja liniowa: (value - min_right) / (max_right - min_right) * (max_left - min_left) + min_left
    efficiencyTransformed[i] = (efficiencyValues[i] - yRightMin) / (yRightMax - yRightMin) * (yLeftMax - yLeftMin) + yLeftMin;
    purityTransformed[i] = (purityValues[i] - yRightMin) / (yRightMax - yRightMin) * (yLeftMax - yLeftMin) + yLeftMin;
  }

  // Wykres dla wydajności (efficiency)
  TGraph *efficiencyGraph = new TGraph(limitValues.size(), limitValues.data(), efficiencyTransformed.data());
  efficiencyGraph->SetMarkerStyle(21);
  efficiencyGraph->SetMarkerColor(kBlue);
  efficiencyGraph->SetLineColor(kBlue);
  efficiencyGraph->SetLineWidth(2);
  efficiencyGraph->Draw("LP SAME");

  // Wykres dla czystości (purity)
  TGraph *purityGraph = new TGraph(limitValues.size(), limitValues.data(), purityTransformed.data());
  purityGraph->SetMarkerStyle(22);
  purityGraph->SetMarkerColor(kRed);
  purityGraph->SetLineColor(kRed);
  purityGraph->SetLineWidth(2);
  purityGraph->Draw("LP SAME");
  canva->Update();

  // Utwórz drugą oś Y (po prawej stronie)
  TGaxis *axis2 = new TGaxis(canva->GetUxmax(), canva->GetUymin(),
                             canva->GetUxmax(), canva->GetUymax(),
                             0, 1, 510, "+L");
  axis2->SetTitle("Efficiency or Purity");
  axis2->SetTitleColor(kBlack);
  axis2->SetLabelColor(kBlack);
  axis2->SetTitleOffset(1.2);
  axis2->Draw();

    // ✅ LEGENDA W WSPÓŁRZĘDNYCH NDC (Normalized Device Coordinates)
  // Dla marginesu prawego 0.25: obszar od 0.75 do 1.0
  TLegend *legend = new TLegend(0.82, 0.15, 0.98, 0.35);  // NDC coordinates
  legend->AddEntry(significanceGraph, "S/#sqrt{S+B}", "LP");
  legend->AddEntry(efficiencyGraph, "Efficiency", "LP");
  legend->AddEntry(purityGraph, "Purity", "LP");
  legend->SetBorderSize(1);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.035);
  legend->Draw();

  canva->Update();

  return canva;
}

// Funkcja do adaptacyjnego próbkowania
void AdaptiveSampling(std::vector<Double_t> &limitValues, std::vector<Double_t> &significanceValues, std::vector<Double_t> &efficiencyValues, std::vector<Double_t> &purityValues, SignalCutOptimizer *optimizer, TChain *chain, Int_t chainEntries, Double_t minStep = 2.0, Double_t maxStep = 20.0, Double_t gradientThreshold = 0.05)
{
  // Guard clause - walidacja danych
  if (!chain || chainEntries == 0)
  {
    std::cerr << "ERROR: Invalid chain in AdaptiveSampling()" << std::endl;
    return;
  }

  std::vector<Double_t> tempLimitValues, tempSignValues, tempEffValues, tempPurValues;

  Double_t currentLimitValue = 0.0;
  Double_t maxLimitValue = 200.0;
  Double_t currentStep = maxStep;

  while (currentLimitValue <= maxLimitValue)
  {
    Long64_t NCutSignal = 0;
    Long64_t NCutBackground = 0;

    Double_t significance = 0;
    Double_t purity = 0;
    Double_t efficiency = 0;

    SignalCutOptimizer *signalOptimizer = new SignalCutOptimizer(chain, currentLimitValue);
    chain->Process(signalOptimizer, "", chainEntries);

    NCutSignal = signalOptimizer->GetNCut();
    significance = (NCutSignal + chainEntries > 0) ? NCutSignal / TMath::Sqrt(NCutSignal + chainEntries) : 0;
    efficiency = chainEntries > 0 ? (Double_t)NCutSignal / chainEntries : 0;

    tempLimitValues.push_back(currentLimitValue);
    tempSignValues.push_back(significance);
    tempEffValues.push_back(efficiency);

    signalOptimizer->Delete();

    // Adaptacyjna zmiana kroku na podstawie gradientu
    if (tempSignValues.size() >= 2)
    {
      Double_t gradient = TMath::Abs(tempSignValues.back() - tempSignValues[tempSignValues.size() - 2]) / (currentStep + 1e-6);

      if (gradient > gradientThreshold)
      {
        // Szybka zmiana - zmniejsz krok (gęstsze próbkowanie)
        currentStep = TMath::Max(minStep, currentStep * 0.8);
      }
      else
      {
        // Wolna zmiana - zwiększ krok (rzadsze próbkowanie)
        currentStep = TMath::Min(maxStep, currentStep * 1.2);
      }

      std::cout << "LimitValue: " << currentLimitValue << ", Step: " << currentStep << ", Gradient: " << gradient << std::endl;
    }

    currentLimitValue += currentStep;
  }

  // Kopiuj wyniki
  limitValues = tempLimitValues;
  significanceValues = tempSignValues;
  efficiencyValues = tempEffValues;
  purityValues = tempPurValues;
}