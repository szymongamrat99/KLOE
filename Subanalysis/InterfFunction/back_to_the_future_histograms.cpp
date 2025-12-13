#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TF1.h>

#include <const.h>
#include <interf_function.h>

int main(int argc, char* argv[])
{
  KLOE::setGlobalStyle();

  TRandom3 *rand = new TRandom3(0);

  TString pm00_dir = Paths::initialanalysis_dir + Paths::root_files_dir + "2025-11-24/" + "mk0_initial_analysis_*_Signal*.root";

  TString pmpm_dir = Paths::initialanalysis_dir + Paths::root_files_dir + "2025-08-23/" + "mk0_initial_analysis_*.root";

  TChain *chain_pm00 = new TChain("h1");
  chain_pm00->Add(pm00_dir);

  TChain *chain_pmpm = new TChain("h1");
  chain_pmpm->Add(pmpm_dir);

  Int_t binsX = 20, bins = 20;

  TH1 *h_pm00 = new TH1F("h_pm00", ";t_{1} [#tau_{S}];Counts", binsX, 0, 20);

  TH2 *h_pm00_2D = new TH2F("h_pm00_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", bins, 0, 20, bins, 0, 20);

  TH1 *h_00pm = new TH1F("h_00pm", ";t_{1} [#tau_{S}];Counts", bins, 0, 20);

  TH2 *h_00pm_2D = new TH2F("h_00pm_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", bins, 0, 20, bins, 0, 20);
  TH1 *h_delta = new TH1F("h_delta", "t1 vs t2 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", bins, -20, 20);

  TH1 *h_pmpm = new TH1F("h_pmpm", ";t_{1} [#tau_{S}];Counts", bins, 0, 20);

  TH2 *h_pmpm_2D = new TH2F("h_pmpm_2D", ";t_{1} [#tau_{S}];t_{2} [#tau_{S}]", bins, 0, 20, bins, 0, 20);

  TH1 *h_delta_pmpm = new TH1F("h_delta_pmpm", "t1 vs t2 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", bins, -20, 20);

  TTreeReader reader(chain_pm00);
  TTreeReader readerpmpm(chain_pmpm);

  TTreeReaderValue<Float_t> tNeuMC(reader, "KaonNeTimeCMMC");
  TTreeReaderValue<Float_t> tChMC(reader, "KaonChTimeCMMC");
  TTreeReaderValue<Int_t> mctruth(reader, "mctruth");

  TTreeReaderValue<Float_t> tCh2MCpmpm(readerpmpm, "KaonNeTimeCMMC");
  TTreeReaderValue<Float_t> tCh1MCpmpm(readerpmpm, "KaonChTimeCMMC");
  TTreeReaderValue<Int_t> mctruthpmpm(readerpmpm, "mctruth");

  Float_t t1, t2;

  Double_t x[2], par[1];

  while (reader.Next())
  {
    if (*mctruth != 1)
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

    
    h_pm00->Fill(t1, weight);
    h_delta->Fill(t2 - t1, weight);
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

    h_00pm->Fill(t1, weight);
    h_delta->Fill(t2 - t1, weight);
    h_00pm_2D->Fill(t1, t2, weight);
  }

  while (readerpmpm.Next())
  {
    if (*mctruthpmpm != 7)
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

    h_pmpm->Fill(t1, weight);
    h_delta_pmpm->Fill(t2 - t1, weight);
    h_pmpm_2D->Fill(t1, t2, weight);
  }

  h_00pm->Sumw2();
  h_pm00->Sumw2();
  h_pmpm->Sumw2();

  h_00pm->Scale(1 / h_00pm->GetMaximum());
  h_pm00->Scale(1 / h_pm00->GetMaximum());
  h_pmpm->Scale(1 / h_pmpm->GetMaximum());

  TF1 *f1_00pm = new TF1("f1_00pm", "expo", 0, 20);
  f1_00pm->SetParameter(2, 0.25);
  f1_00pm->SetParameter(3, 10.0);
  f1_00pm->SetParameter(0, h_00pm->GetEntries());
  f1_00pm->SetParameters(1, h_00pm->GetMaximum());
  h_00pm->Fit(f1_00pm, "R");

  h_00pm_2D->Scale(h_00pm_2D->GetEntries() / h_00pm_2D->Integral(0, h_00pm_2D->GetNbinsX() + 1, 0, h_00pm_2D->GetNbinsY() + 1));
  h_pm00_2D->Scale(h_pm00_2D->GetEntries() / h_pm00_2D->Integral(0, h_pm00_2D->GetNbinsX() + 1, 0, h_pm00_2D->GetNbinsY() + 1));
  h_pmpm_2D->Scale(h_pmpm_2D->GetEntries() / h_pmpm_2D->Integral(0, h_pmpm_2D->GetNbinsX() + 1, 0, h_pmpm_2D->GetNbinsY() + 1));

  TH1 *h_pm00_projX = h_pm00_2D->ProjectionX("h_pm00_projX", 0, -1, "e");
  TH1 *h_pmpm_projX = h_pmpm_2D->ProjectionX("h_pmpm_projX", 0, -1, "e");
  TH1 *h_00pm_projX = h_00pm_2D->ProjectionX("h_00pm_projX", 0, -1, "e");

  h_pm00_projX->SetTitle(";t_{1} [#tau_{S}];Counts");
  h_00pm_projX->SetTitle(";t_{1} [#tau_{S}];Counts");
  h_pmpm_projX->SetTitle(";t_{1} [#tau_{S}];Counts");

  h_pm00_projX->Scale(1 / h_pm00_projX->GetMaximum());
  h_pmpm_projX->Scale(1 / h_pmpm_projX->GetMaximum());
  h_00pm_projX->Scale(1 / h_00pm_projX->GetMaximum());

  /////////////////////////////////////////////////////////////////////////
  // Cloning of histograms for BTTF analysis

  TH1 *h_bttf_RA = (TH1 *)h_pmpm->Clone("h_bttf_RA");
  TH1 *h_bttf_RB = (TH1 *)h_pmpm->Clone("h_bttf_RB");

  h_bttf_RA->SetTitle("R_{A}(t_{1});t_{1} [#tau_{S}];R_{A}(t_{1}) [-]");
  h_bttf_RB->SetTitle("R_{B}(t_{1});t_{1} [#tau_{S}];R_{B}(t_{1}) [-]");

  ////////////////////////////////////////////////////////////////////////
  h_bttf_RA->Divide(h_pmpm_projX, h_00pm_projX, 1, 1);
  h_bttf_RB->Divide(h_pmpm_projX, h_pm00_projX, 1, 1);

  TH1 *h_bttf_RARB = (TH1 *)h_bttf_RA->Clone("h_bttf_RARB");
  h_bttf_RARB->SetTitle("R_{C}(t_{1});t_{1} [#tau_{S}];R_{C}(t_{1}) [-]");
  h_bttf_RARB->Divide(h_bttf_RB);

  /////////////////////////////////////////////////////////////////////////

  TCanvas *c_pm00 = new TCanvas("c_pm00", "t1 distribution for K_{S}K_{L}#rightarrow#pi^{#pm}#pi^{0}#pi^{0}", 800, 600);
  h_pm00->Draw("HIST");
  c_pm00->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pm00.png");

  TCanvas *c_00pm = new TCanvas("c_00pm", "t1 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", 800, 600);
  // c_00pm->SetLogy(1);
  h_00pm->Draw("HIST");
  f1_00pm->Draw("SAME");
  c_00pm->SaveAs(Paths::img_dir + "interf_func_bttf_histo_00pm.png");

  TCanvas *c_delta = new TCanvas("c_delta", "t2 - t1 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", 800, 600);
  h_delta->Draw("HIST");
  c_delta->SaveAs(Paths::img_dir + "interf_func_bttf_histo_delta.png");

  TCanvas *c_pmpm = new TCanvas("c_pmpm", "t1 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", 800, 600);
  h_pmpm->Draw("HIST");
  c_pmpm->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pmpm.png");

  TCanvas *c_delta_pmpm = new TCanvas("c_delta_pmpm", "t2 - t1 distribution for K_{S}K_{L}#rightarrow#pi^{0}#pi^{0}#pi^{#pm}", 800, 600);
  h_delta_pmpm->Draw("HIST");
  c_delta_pmpm->SaveAs(Paths::img_dir + "interf_func_bttf_histo_delta_pmpm.png");

  // 2D

  TCanvas *c_00pm_2D = new TCanvas("c_00pm_2D", ";t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", 800, 600);
  c_00pm_2D->SetLogz(1);
  h_00pm_2D->Draw("COLZ");
  c_00pm_2D->SaveAs(Paths::img_dir + "interf_func_bttf_histo_00pm_2D.png");

  TCanvas *c_pm00_2D = new TCanvas("c_pm00_2D", ";t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", 800, 600);
  c_pm00_2D->SetLogz(1);
  h_pm00_2D->Draw("COLZ");
  c_pm00_2D->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pm00_2D.png");

  TCanvas *c_pmpm_2D = new TCanvas("c_pmpm_2D", ";t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", 800, 600);
  c_pmpm_2D->SetLogz(1);
  h_pmpm_2D->Draw("COLZ");
  c_pmpm_2D->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pmpm_2D.png");

  // Projections

  TCanvas *c_pm00_projX = new TCanvas("c_pm00_projX", ";t_{1} [#tau_{S}]; Counts", 800, 600);
  h_pm00_projX->Draw("HIST");
  c_pm00_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pm00_projX.png");

  TCanvas *c_00pm_projX = new TCanvas("c_00pm_projX", ";t_{1} [#tau_{S}]; Counts", 800, 600);
  h_00pm_projX->Draw("HIST");
  c_00pm_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_00pm_projX.png");

  TCanvas *c_pmpm_projX = new TCanvas("c_pmpm_projX", "", 800, 600);
  h_pmpm_projX->Draw("HIST");
  c_pmpm_projX->SaveAs(Paths::img_dir + "interf_func_bttf_histo_pmpm_projX.png");

  // BTTF canvases
  TCanvas *c_bttf_RA = new TCanvas("c_bttf_RA", "RA histogram;t_{1} [#tau_{S}];R_{A}(t_{1})", 800, 600);
  h_bttf_RA->Draw("PE1");
  c_bttf_RA->SaveAs(Paths::img_dir + "interf_func_bttf_histo_RA.png");

  TCanvas *c_bttf_RB = new TCanvas("c_bttf_RB", "RB histogram;t_{1} [#tau_{S}];R_{B}(t_{1})", 800, 600);
  h_bttf_RB->Draw("PE1");
  c_bttf_RB->SaveAs(Paths::img_dir + "interf_func_bttf_histo_RB.png");

  TCanvas *c_bttf_RARB = new TCanvas("c_bttf_RARB", "RC histogram;t_{1} [#tau_{S}];R_{C}(t_{1})", 800, 600);
  h_bttf_RARB->Draw("PE1");
  c_bttf_RARB->SaveAs(Paths::img_dir + "interf_func_bttf_histo_RARB.png");

  return 0;
}