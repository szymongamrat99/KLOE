#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TF1.h>

#include <const.h>
#include <interf_function.h>

struct Config
{
  Double_t  tMinX = 0.0,
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

  Int_t binsX = (Config::tMaxX - Config::tMinX) / Config.resT1, binsY = (Config::tMaxY - Config::tMinY) / Config.resT2;

  TH2 *h_pm00_2D = new TH2F("h_pm00_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", binsX, tMinX, tMaxX, binsY, tMinY, tMaxY);
  TH2 *h_00pm_2D = new TH2F("h_00pm_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", binsX, tMinX, tMaxX, binsY, tMinY, tMaxY);
  TH2 *h_pmpm_2D = new TH2F("h_pmpm_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", binsX, tMinX, tMaxX, binsY, tMinY, tMaxY);

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

      weight = interf_function_pm00(x, par);
    }
    else
    {
      x[0] = t2;
      x[1] = t1;

      weight = interf_function_00pm(x, par);
    }

    h_pm00_2D->Fill(t1, t2, weight);

    t1 = *tNeuMC;
    t2 = *tChMC;

    if (t1 < t2)
    {
      x[0] = t1;
      x[1] = t2;

      weight = interf_function_00pm(x, par);
    }
    else
    {
      x[0] = t2;
      x[1] = t1;

      weight = interf_function_pm00(x, par);
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

      weight = interf_function_pmpm(x, par);
    }
    else
    {
      x[0] = t2;
      x[1] = t1;

      weight = interf_function_pmpm(x, par);
    }

    h_delta_pmpm->Fill(t2 - t1, weight);
    h_pmpm_2D->Fill(t1, t2, weight);
  }

  // Normalization of 2D histograms
  h_00pm_2D->Scale(h_00pm_2D->GetEntries() / h_00pm_2D->Integral(0, h_00pm_2D->GetNbinsX() + 1, 0, h_00pm_2D->GetNbinsY() + 1));
  h_pm00_2D->Scale(h_pm00_2D->GetEntries() / h_pm00_2D->Integral(0, h_pm00_2D->GetNbinsX() + 1, 0, h_pm00_2D->GetNbinsY() + 1));
  h_pmpm_2D->Scale(h_pmpm_2D->GetEntries() / h_pmpm_2D->Integral(0, h_pmpm_2D->GetNbinsX() + 1, 0, h_pmpm_2D->GetNbinsY() + 1));

  Int_t binypmpmForpm00 = -1, binypmpmFor00pm = -1, binypm00 = -1, biny00pm = -1;

  Double_t integrationMax = tMaxY; // Set integration max to histogram max y

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

  h_pmpm_forpm00_projX->Scale(1 / h_pmpm_forpm00_projX->GetMaximum());
  h_pm00_projX->Scale(1 / h_pm00_projX->GetMaximum());

  h_pmpm_for00pm_projX->Scale(1 / h_pmpm_for00pm_projX->GetMaximum());
  h_00pm_projX->Scale(1 / h_00pm_projX->GetMaximum());

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