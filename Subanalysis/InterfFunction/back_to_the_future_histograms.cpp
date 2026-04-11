#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TF1.h>

#include <const.h>
#include <interf_function.h>

Bool_t NormalizeOversampledRegion(TH2 *hist2D, const TString &histName,
                                   Double_t oversampleMax = 20.0,
                                   Double_t neighborhoodMax = 300.0)
{
  if (!hist2D)
  {
    std::cerr << "ERROR: Null histogram passed to NormalizeOversampledRegion for " << histName << std::endl;
    return kFALSE;
  }

  const Int_t nx = hist2D->GetNbinsX();
  const Int_t ny = hist2D->GetNbinsY();

  Double_t oversampleCounts = 0.0;
  Double_t oversampleArea = 0.0;
  Double_t neighborhoodCounts = 0.0;
  Double_t neighborhoodArea = 0.0;

  for (Int_t ix = 1; ix <= nx; ++ix)
  {
    const Double_t xCenter = hist2D->GetXaxis()->GetBinCenter(ix);
    const Double_t dx = hist2D->GetXaxis()->GetBinWidth(ix);

    for (Int_t iy = 1; iy <= ny; ++iy)
    {
      const Double_t yCenter = hist2D->GetYaxis()->GetBinCenter(iy);
      const Double_t dy = hist2D->GetYaxis()->GetBinWidth(iy);
      const Double_t binArea = dx * dy;
      const Double_t content = hist2D->GetBinContent(ix, iy);

      // Oversampled region: [0, oversampleMax] x [0, oversampleMax]
      const Bool_t inOversample = (xCenter >= 0.0 && xCenter < oversampleMax && 
                                   yCenter >= 0.0 && yCenter < oversampleMax);

      // Neighborhood: [0, oversampleMax] x [oversampleMax, neighborhoodMax]
      //           or [oversampleMax, neighborhoodMax] x [0, oversampleMax]
      const Bool_t inHorizontalNeighbor = (xCenter >= 0.0 && xCenter < oversampleMax && 
                                           yCenter >= oversampleMax && yCenter <= neighborhoodMax);
      const Bool_t inVerticalNeighbor = (xCenter >= oversampleMax && xCenter <= neighborhoodMax && 
                                        yCenter >= 0.0 && yCenter < oversampleMax);

      if (inOversample)
      {
        oversampleCounts += content;
        oversampleArea += binArea;
      }

      if (inHorizontalNeighbor || inVerticalNeighbor)
      {
        neighborhoodCounts += content;
        neighborhoodArea += binArea;
      }
    }
  }

  if (oversampleArea <= 0.0 || neighborhoodArea <= 0.0 || oversampleCounts <= 0.0 || neighborhoodCounts <= 0.0)
  {
    std::cerr << "WARNING: Cannot normalize " << histName
              << " (oversampleArea=" << oversampleArea
              << ", neighborhoodArea=" << neighborhoodArea
              << ", oversampleCounts=" << oversampleCounts
              << ", neighborhoodCounts=" << neighborhoodCounts << ")" << std::endl;
    return kFALSE;
  }

  const Double_t oversampleDensity = oversampleCounts / oversampleArea;
  const Double_t neighborhoodDensity = neighborhoodCounts / neighborhoodArea;
  const Double_t scale = neighborhoodDensity / oversampleDensity;

  // Scale only the oversampled region
  for (Int_t ix = 1; ix <= nx; ++ix)
  {
    const Double_t xCenter = hist2D->GetXaxis()->GetBinCenter(ix);
    if (!(xCenter >= 0.0 && xCenter < oversampleMax))
      continue;

    for (Int_t iy = 1; iy <= ny; ++iy)
    {
      const Double_t yCenter = hist2D->GetYaxis()->GetBinCenter(iy);
      if (!(yCenter >= 0.0 && yCenter < oversampleMax))
        continue;

      const Double_t oldContent = hist2D->GetBinContent(ix, iy);
      const Double_t oldError = hist2D->GetBinError(ix, iy);
      hist2D->SetBinContent(ix, iy, oldContent * scale);
      hist2D->SetBinError(ix, iy, oldError * scale);
    }
  }

  std::cout << "Normalized oversampled region [0," << oversampleMax << "]x[0," << oversampleMax << "] for " << histName
            << " to neighborhood density: oversampleDensity=" << oversampleDensity
            << ", neighborhoodDensity=" << neighborhoodDensity
            << ", scale=" << scale << std::endl;

  return kTRUE;
}

struct Config
{
  Double_t tMinX = 0.0,
           tMinY = 0.0,
           tMaxX = 20.0,
           tMaxY = 300.0,
           resT1 = 1.0,
           resT2 = 1.0;

  Int_t binsX = (tMaxX - tMinX) / resT1,
        binsY = (tMaxY - tMinY) / resT2;

  Double_t
      physicsRe = PhysicsConstants::Re,
      physicsIm = PhysicsConstants::Im_nonCPT;
};

int main(int argc, char *argv[])
{
  KLOE::setGlobalStyle();

  TRandom3 *rand = new TRandom3(0);

  TString
      pm00_dir = Paths::initialanalysis_dir + Paths::root_files_dir + "2025-11-24/" + "mk0_initial_analysis_*_Signal*.root",
      pm00_dir_all_phys3 = Paths::initialanalysis_dir + Paths::root_files_dir + "2025-11-26/" + "mk0_initial_analysis_*_Signal*.root";

  TString pmpm_dir = Paths::initialanalysis_dir + Paths::root_files_dir + "2026-01-14/" + "mk0_initial_analysis_*.root";

  TChain *chain_pm00 = new TChain("h1");
  chain_pm00->Add(pm00_dir);
  chain_pm00->Add(pm00_dir_all_phys3);

  TChain *chain_pmpm = new TChain("h1");
  chain_pmpm->Add(pmpm_dir);

  Config config;

  Int_t binsX = (config.tMaxX - config.tMinX) / config.resT1, binsY = (config.tMaxY - config.tMinY) / config.resT2;

  TH2 *h_pm00_2D = new TH2F("h_pm00_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", binsX, config.tMinX, config.tMaxX, binsY, config.tMinY, config.tMaxY);
  TH2 *h_00pm_2D = new TH2F("h_00pm_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", binsX, config.tMinX, config.tMaxX, binsY, config.tMinY, config.tMaxY);
  TH2 *h_pmpm_2D = new TH2F("h_pmpm_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", binsX, config.tMinX, config.tMaxX, binsY, config.tMinY, config.tMaxY);

  Double_t resDelta = 1.5,
           deltaMin = -50.0,
           deltaMax = 50.0,
           binsDelta = 1 + (deltaMax - deltaMin) / resDelta;

  TH1 *h_delta = new TH1F("h_delta", "t1 vs t2 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", binsDelta, deltaMin, deltaMax);
  TH1 *h_delta_pmpm = new TH1F("h_delta_pmpm", "t1 vs t2 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", binsDelta, deltaMin, deltaMax);

  TTreeReader reader(chain_pm00);
  TTreeReader readerpmpm(chain_pmpm);

  std::cout << chain_pm00->GetEntries("mctruth == 1 || mctruth == 0") << " entries in pm00 chain." << std::endl;
  std::cout << chain_pmpm->GetEntries("mctruth == 7 || mctruth == 0") << " entries in pmpm chain." << std::endl;

  TTreeReaderValue<Float_t> tNeuMC(reader, "KaonNeTimeCMMC");
  TTreeReaderValue<Float_t> tChMC(reader, "KaonChTimeCMMC");
  TTreeReaderValue<Int_t> mctruth(reader, "mctruth");

  TTreeReaderValue<Float_t> tCh2MCpmpm(readerpmpm, "KaonNeTimeCMMC");
  TTreeReaderValue<Float_t> tCh1MCpmpm(readerpmpm, "KaonChTimeCMMC");
  TTreeReaderValue<Int_t> mctruthpmpm(readerpmpm, "mctruth");

  Float_t t1, t2;

  Double_t x[2], par[2];

  Printf("Do You want to use default physics constants values? (1 - yes, 0 - no): ");
  Int_t useDefault;
  std::cin >> useDefault;

  switch (useDefault)
  {
  case 0:
  {
    Printf("Set Re: ");
    std::cin >> par[0];
    Printf("Set Im (in rad): ");
    std::cin >> par[1];
    break;
  }
  case 1:
  {
    par[0] = PhysicsConstants::Re;
    par[1] = PhysicsConstants::Im_nonCPT;
    break;
  }
  default:
  {
    Printf("Wrong option selected. Exiting...");
    return 1;
  }
  }

  Double_t gammaS = 1.0; // Wartość gamma_S do ustawienia zakresów
  Double_t gammaL = PhysicsConstants::tau_S_nonCPT / PhysicsConstants::tau_L;

  auto double_exp = [gammaS, gammaL](Double_t t1, Double_t t2)
  {
    return exp(-gammaS * t1 - gammaL * t2) + exp(-gammaL * t1 - gammaS * t2);
  };

  while (reader.Next())
  {
    if (*mctruth != 1 && *mctruth != 0)
      continue;

    // Filling for pm00

    t1 = *tChMC;
    t2 = *tNeuMC;

    Double_t weight = 1.0;

    if (t1 < t2)
    {
      x[0] = t1;
      x[1] = t2;

      weight = interf_function_pm00(x, par) / double_exp(t1, t2);
    }
    else
    {
      x[0] = t2;
      x[1] = t1;

      weight = interf_function_00pm(x, par) / double_exp(t1, t2);
    }

    h_pm00_2D->Fill(t1, t2, weight);

    t1 = *tNeuMC;
    t2 = *tChMC;

    if (t1 < t2)
    {
      x[0] = t1;
      x[1] = t2;

      weight = interf_function_00pm(x, par) / double_exp(t1, t2);
    }
    else
    {
      x[0] = t2;
      x[1] = t1;

      weight = interf_function_pm00(x, par) / double_exp(t1, t2);
    }

    h_delta->Fill(t2 - t1, weight);
    h_00pm_2D->Fill(t1, t2, weight);
  }

  while (readerpmpm.Next())
  {
    if (*mctruthpmpm != 7 && *mctruthpmpm != 0)
      continue;

    // Filling for pmpm

    Double_t weight = 1.0;

    int wynik = rand->Binomial(1, 0.5);

    if (wynik)
    {
      t1 = *tCh1MCpmpm;
      t2 = *tCh2MCpmpm;
    }
    else
    {
      t1 = *tCh2MCpmpm;
      t2 = *tCh1MCpmpm;
    }

    if (t1 < t2)
    {
      x[0] = t1;
      x[1] = t2;

      weight = interf_function_pmpm(x, par) / double_exp(t1, t2);
    }
    else
    {
      x[0] = t2;
      x[1] = t1;

      weight = interf_function_pmpm(x, par) / double_exp(t1, t2);
    }

    h_delta_pmpm->Fill(t2 - t1, weight);
    h_pmpm_2D->Fill(t1, t2, weight);
  }

  // Normalization of 2D histograms
  h_00pm_2D->Scale(h_00pm_2D->GetEntries() / h_00pm_2D->Integral(0, h_00pm_2D->GetNbinsX() + 1, 0, h_00pm_2D->GetNbinsY() + 1));
  h_pm00_2D->Scale(h_pm00_2D->GetEntries() / h_pm00_2D->Integral(0, h_pm00_2D->GetNbinsX() + 1, 0, h_pm00_2D->GetNbinsY() + 1));
  h_pmpm_2D->Scale(h_pmpm_2D->GetEntries() / h_pmpm_2D->Integral(0, h_pmpm_2D->GetNbinsX() + 1, 0, h_pmpm_2D->GetNbinsY() + 1));

  NormalizeOversampledRegion(h_00pm_2D, "h_00pm_2D", 20.0, 300.0);
  NormalizeOversampledRegion(h_pm00_2D, "h_pm00_2D", 20.0, 300.0);
  NormalizeOversampledRegion(h_pmpm_2D, "h_pmpm_2D", 20.0, 300.0);

  Int_t binypmpmForpm00 = -1, binypmpmFor00pm = -1, binypm00 = -1, biny00pm = -1;

  Double_t integrationMax = config.tMaxY; // Set integration max to histogram max y

  Bool_t chooseIntegration = kTRUE;

  std::cout << "Do you want to set custom y integration limits? (1 - yes, 0 - no): ";
  std::cin >> chooseIntegration;

  if (!chooseIntegration)
    goto skipIntegrationInput;

  std::cout << "Set y val for integration - pmpm for 00pm: ";
  std::cin >> integrationMax;
  binypmpmFor00pm = h_pmpm_2D->GetYaxis()->FindBin(integrationMax);
  biny00pm = h_00pm_2D->GetYaxis()->FindBin(integrationMax);

  std::cout << "Set y val for integration - pmpm for pm00: ";
  std::cin >> integrationMax;
  binypmpmForpm00 = h_pmpm_2D->GetYaxis()->FindBin(integrationMax);

  std::cout << "Set y val for integration - pm00: ";
  std::cin >> integrationMax;
  binypm00 = h_pm00_2D->GetYaxis()->FindBin(integrationMax);

skipIntegrationInput:

  TH1 *h_pmpm_forpm00_projX = h_pmpm_2D->ProjectionX("h_pmpm_forpm00_projX", 0, binypmpmForpm00, "e");
  TH1 *h_pm00_projX = h_pm00_2D->ProjectionX("h_pm00_projX", 0, binypm00, "e");

  TH1 *h_pmpm_for00pm_projX = h_pmpm_2D->ProjectionX("h_pmpm_for00pm_projX", 0, binypmpmFor00pm, "e");
  TH1 *h_00pm_projX = h_00pm_2D->ProjectionX("h_00pm_projX", 0, biny00pm, "e");

  h_pmpm_forpm00_projX->SetTitle(";t_{1} [#tau_{S}];Counts");
  h_pm00_projX->SetTitle(";t_{1} [#tau_{S}];Counts");

  h_pmpm_for00pm_projX->SetTitle(";t_{1} [#tau_{S}];Counts");
  h_00pm_projX->SetTitle(";t_{1} [#tau_{S}];Counts");

  h_pmpm_forpm00_projX->Scale(1 / h_pmpm_forpm00_projX->Integral(-1, h_pmpm_forpm00_projX->GetNbinsX() + 1));
  h_pm00_projX->Scale(1 / h_pm00_projX->Integral(-1, h_pm00_projX->GetNbinsX() + 1));

  h_pmpm_for00pm_projX->Scale(1 / h_pmpm_for00pm_projX->Integral(-1, h_pmpm_for00pm_projX->GetNbinsX() + 1));
  h_00pm_projX->Scale(1 / h_00pm_projX->Integral(-1, h_00pm_projX->GetNbinsX() + 1));

  /////////////////////////////////////////////////////////////////////////
  // Cloning of histograms for BTTF analysis

  TH1 *h_bttf_RA = (TH1 *)h_pmpm_for00pm_projX->Clone("h_bttf_RA");
  TH1 *h_bttf_RB = (TH1 *)h_pmpm_forpm00_projX->Clone("h_bttf_RB");

  h_bttf_RA->SetTitle("R_{A}(t_{1});t_{1} [#tau_{S}];R_{A}(t_{1}) [-]");
  h_bttf_RB->SetTitle("R_{B}(t_{1});t_{1} [#tau_{S}];R_{B}(t_{1}) [-]");

  ////////////////////////////////////////////////////////////////////////
  h_bttf_RA->Divide(h_00pm_projX);
  h_bttf_RB->Divide(h_pm00_projX);

  TH1 *h_bttf_RARB = (TH1 *)h_bttf_RA->Clone("h_bttf_RARB");
  h_bttf_RARB->SetTitle("R_{C}(t_{1});t_{1} [#tau_{S}];R_{C}(t_{1}) [-]");
  h_bttf_RARB->Divide(h_bttf_RB);

  /////////////////////////////////////////////////////////////////////////

  TCanvas *c_delta_pmpm = new TCanvas("c_delta_pmpm", "", 800, 600);
  h_delta_pmpm->GetYaxis()->SetRangeUser(0, h_delta_pmpm->GetMaximum() * 1.1);
  h_delta_pmpm->GetXaxis()->SetTitle("t_{2} - t_{1} [#tau_{S}]");
  h_delta_pmpm->Draw("HIST");
  c_delta_pmpm->SaveAs(Paths::img_dir + "interf_func_bttf_histo_delta_pmpm.svg");

  TCanvas *c_delta = new TCanvas("c_delta", "", 800, 600);
  h_delta->GetYaxis()->SetRangeUser(0, h_delta->GetMaximum() * 1.1);
  h_delta->GetXaxis()->SetTitle("t_{ch} - t_{neu} [#tau_{S}]");
  h_delta->Draw("HIST");
  c_delta->SaveAs(Paths::img_dir + "interf_func_bttf_histo_delta.svg");

  // 2D

  TCanvas *c_00pm_2D = new TCanvas("c_00pm_2D", ";t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", 800, 600);
  c_00pm_2D->SetLogz(1);
  h_00pm_2D->Draw("COLZ");
  c_00pm_2D->SaveAs(Paths::img_dir + "interf_func_bttf_histo_00pm_2D.svg");

  TCanvas *c_pm00_2D = new TCanvas("c_pm00_2D", ";t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", 800, 600);
  c_pm00_2D->SetLogz(1);
  h_pm00_2D->Draw("COLZ");
  c_pm00_2D->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pm00_2D.svg");

  TCanvas *c_pmpm_2D = new TCanvas("c_pmpm_2D", ";t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", 800, 600);
  c_pmpm_2D->SetLogz(1);
  h_pmpm_2D->Draw("COLZ");
  c_pmpm_2D->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pmpm_2D.svg");

  // Projections

  TCanvas *c_pm00_projX = new TCanvas("c_pm00_projX", ";t_{1} [#tau_{S}]; Counts", 800, 600);
  h_pm00_projX->Draw("HIST");
  c_pm00_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pm00_projX.svg");

  TCanvas *c_00pm_projX = new TCanvas("c_00pm_projX", ";t_{1} [#tau_{S}]; Counts", 800, 600);
  h_00pm_projX->Draw("HIST");
  c_00pm_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_00pm_projX.svg");

  TCanvas *c_pmpm_for00pm_projX = new TCanvas("c_pmpm_for00pm_projX", "", 800, 600);
  h_pmpm_for00pm_projX->Draw("HIST");
  c_pmpm_for00pm_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pmpm_for00pm_projX.svg");

  TCanvas *c_pmpm_forpm00_projX = new TCanvas("c_pmpm_forpm00_projX", "", 800, 600);
  h_pmpm_forpm00_projX->Draw("HIST");
  c_pmpm_forpm00_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pmpm_forpm00_projX.svg");

  // BTTF canvases
  TCanvas *c_bttf_RA = new TCanvas("c_bttf_RA", "RA histogram;t_{1} [#tau_{S}];R_{A}(t_{1})", 800, 600);
  h_bttf_RA->Draw("PE1");
  c_bttf_RA->SaveAs(Paths::img_dir + "interf_func_bttf_histo_RA.svg");

  TCanvas *c_bttf_RB = new TCanvas("c_bttf_RB", "RB histogram;t_{1} [#tau_{S}];R_{B}(t_{1})", 800, 600);
  h_bttf_RB->Draw("PE1");
  c_bttf_RB->SaveAs(Paths::img_dir + "interf_func_bttf_histo_RB.svg");

  TCanvas *c_bttf_RARB = new TCanvas("c_bttf_RARB", "RC histogram;t_{1} [#tau_{S}];R_{C}(t_{1})", 800, 600);
  h_bttf_RARB->Draw("PE1");
  c_bttf_RARB->SaveAs(Paths::img_dir + "interf_func_bttf_histo_RARB.svg");

  return 0;
}