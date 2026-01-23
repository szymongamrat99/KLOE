#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLegend.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TRandom3.h>

#include <boost/filesystem.hpp>

#include <interf_function.h>
#include <const.h>

void createDirIfNotExists(const TString& path)
{
  boost::filesystem::path dirPath(path.Data());

  if (!boost::filesystem::exists(dirPath))
    boost::filesystem::create_directories(dirPath);
}

int main(int argc, char *argv[])
{
  Float_t Tpmpm1 = 300,
          Tpmpm2 = 300,
          Tpm00 = 300,
          T00pm = 300;

  if (argc >= 5)
  {
    Tpmpm1 = atof(argv[1]);
    Tpmpm2 = atof(argv[2]);
    Tpm00 = atof(argv[3]);
    T00pm = atof(argv[4]);
  }

  KLOE::setGlobalStyle();

  std::cout << "Do You want to set custom parameters? (1 - yes, 0 - no): ";
  int customParams = 0;
  std::cin >> customParams;

  Double_t reParam = PhysicsConstants::Re;
  ;
  Double_t imParam = PhysicsConstants::Im_nonCPT;

  if (customParams)
  {
    std::cout << "Set Re(epsilon'/epsilon) (default " << PhysicsConstants::Re << "): ";
    std::cin >> reParam;

    std::cout << "Set Im(epsilon'/epsilon) (default " << PhysicsConstants::Im_nonCPT << "): ";
    std::cin >> imParam;
  }

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, 0.0, 300.0, 0.0, 300.0, 2);
  func_00pm->SetParameters(reParam, imParam);

  TF2 *func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, 0.0, 300.0, 0.0, 300.0, 2);
  func_pm00->SetParameters(reParam, imParam);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, 0.0, 300.0, 0.0, 300.0, 2);
  func_pmpm->SetParameters(reParam, imParam);

  // Sampling functions to create 2D histograms
  const Float_t sigmaT = 1.0; // Time resolution in tau_S units
  const Float_t t1Min = 0.0, t1Max = 300.0;
  const Float_t t2Min = 0.0, t2Max = 300.0;
  const Int_t nBinst1 = (t1Max - t1Min) / sigmaT, nBinst2 = (t2Max - t2Min) / sigmaT;
  Long64_t nSamplesTarget = 100000; // Number of samples for MC sampling
  Double_t maxIntegral = 300.0;

  std::cout << "Input max integral for sampling (default 300.0): ";
  std::cin >> maxIntegral;

  std::cout << "Set number of samples for MC sampling (default 1e5): ";
  std::cin >> nSamplesTarget;

  Long64_t nSamples = 0;

  TH2D *hist_00pm = new TH2D("hist_00pm", "I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);
  TH2D *hist_pm00 = new TH2D("hist_pm00", "I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);
  TH2D *hist_pmpm = new TH2D("hist_pmpm", "I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);

  // Find maximum values for accept-reject
  Double_t max_00pm = func_00pm->GetMaximum();
  Double_t max_pm00 = func_pm00->GetMaximum();
  Double_t max_pmpm = func_pmpm->GetMaximum();

  Int_t granularity = 1000;

  func_00pm->SetNpx(granularity);
  func_00pm->SetNpy(granularity);
  func_pm00->SetNpx(granularity);
  func_pm00->SetNpy(granularity);
  func_pmpm->SetNpx(granularity);
  func_pmpm->SetNpy(granularity);

  Double_t t1, t2;

      // Creation of the plot directory if it does not exist
  // Creation of the plot directory with tmp suffix
  TString plot_dir = "theoretical_plots/";
  TString dir_name = "integralLimit_" + TString::Format("%.0f", maxIntegral) + "_samples_" + TString::Format("%lld", nSamplesTarget);
  TString dir_tmp_path = Paths::img_dir + plot_dir + dir_name + "_tmp/";
  TString dir_final_path = Paths::img_dir + plot_dir + dir_name + "/";

  createDirIfNotExists(dir_tmp_path);

  // Sample func_00pm
  while (nSamples < nSamplesTarget)
  {
    func_00pm->GetRandom2(t1, t2);
    hist_00pm->Fill(t1, t2);

    func_pm00->GetRandom2(t1, t2);
    hist_pm00->Fill(t1, t2);

    func_pmpm->GetRandom2(t1, t2);
    hist_pmpm->Fill(t1, t2);

    nSamples++;
  }

  std::cout << "Sampling completed." << std::endl;
  std::cout << "hist_00pm entries: " << hist_00pm->GetEntries() << std::endl;
  std::cout << "hist_pm00 entries: " << hist_pm00->GetEntries() << std::endl;
  std::cout << "hist_pmpm entries: " << hist_pmpm->GetEntries() << std::endl;

  TCanvas *c_func_00pm = new TCanvas("c_func_00pm", "Interference function 00pm", 800, 600);
  c_func_00pm->SetLogz(1);
  hist_00pm->GetXaxis()->SetRangeUser(0, 20);
  hist_00pm->GetYaxis()->SetRangeUser(0, 20);
  hist_00pm->Draw("COLZ");
  c_func_00pm->SaveAs(dir_tmp_path + "interf_func_00pm.pdf");

  TCanvas *c_func_pm00 = new TCanvas("c_func_pm00", "Interference function pm00", 800, 600);
  c_func_pm00->SetLogz(1);
  hist_pm00->GetXaxis()->SetRangeUser(0, 20);
  hist_pm00->GetYaxis()->SetRangeUser(0, 20);
  hist_pm00->Draw("COLZ");
  c_func_pm00->SaveAs(dir_tmp_path + "interf_func_pm00.pdf");

  TCanvas *c_func_pmpm = new TCanvas("c_func_pmpm", "Interference function pmpm", 800, 600);
  c_func_pmpm->SetLogz(1);
  hist_pmpm->GetXaxis()->SetRangeUser(0, 20);
  hist_pmpm->GetYaxis()->SetRangeUser(0, 20);
  hist_pmpm->Draw("COLZ");
  c_func_pmpm->SaveAs(dir_tmp_path + "interf_func_pmpm.pdf");

  // Create 1D projections from 2D histograms with cuts
  // For projection we integrate over y-axis up to a certain limit (par[1] or par[0])

  auto create_projection_y_cut = [](TH2D *hist, Double_t t_max) -> TH1D *
  {
    TString name = TString::Format("%s_projx_tmax%.0f", hist->GetName(), t_max);
    return hist->ProjectionX(name, 0, hist->GetYaxis()->FindBin(t_max));
  };

  // Create projections for different T values
  TH1D *hist_00pm_1D = create_projection_y_cut(hist_00pm, maxIntegral);
  TH1D *hist_pm00_1D = create_projection_y_cut(hist_pm00, maxIntegral);
  TH1D *hist_pmpm_1D_for_RA = create_projection_y_cut(hist_pmpm, maxIntegral);
  TH1D *hist_pmpm_1D_for_RB = create_projection_y_cut(hist_pmpm, maxIntegral);

  hist_00pm_1D->Sumw2();
  hist_pm00_1D->Sumw2();
  hist_pmpm_1D_for_RA->Sumw2();
  hist_pmpm_1D_for_RB->Sumw2();

  hist_00pm_1D->SetTitle("Projection 00pm;t_{1} [#tau_{S}];Events");
  hist_pm00_1D->SetTitle("Projection pm00;t_{1} [#tau_{S}];Events");
  hist_pmpm_1D_for_RA->SetTitle("Projection pmpm (R_{A});t_{1} [#tau_{S}];Events");
  hist_pmpm_1D_for_RB->SetTitle("Projection pmpm (R_{B});t_{1} [#tau_{S}];Events");

  ////////////////////////////////////////////////////////

  // Rysowanie projekcji 1D
  TCanvas *c_pm00 = new TCanvas("c_pm00", "Projection of pm00", 800, 600);
  hist_pm00_1D->GetXaxis()->SetRangeUser(0, 20);
  hist_pm00_1D->Draw("PE1");
  c_pm00->SaveAs(dir_tmp_path + "interf_func_pm001D_draw.pdf");

  ///////////////////////////////////////////////////////////////////

  TCanvas *c_00pm = new TCanvas("c_00pm", "Projection of 00pm", 800, 600);
  hist_00pm_1D->GetXaxis()->SetRangeUser(0, 20);
  hist_00pm_1D->Draw("PE1");
  c_00pm->SaveAs(dir_tmp_path + "interf_func_00pm1D_draw.pdf");

  ////////////////////////////////////////////////////////////////////

  TCanvas *c_pmpm = new TCanvas("c_pmpm", "Projection of pmpm", 800, 600);
  hist_pmpm_1D_for_RA->GetXaxis()->SetRangeUser(0, 20);
  hist_pmpm_1D_for_RA->Draw("PE1");
  c_pmpm->SaveAs(dir_tmp_path + "interf_func_pmpm1D_draw.pdf");

  ////////////////////////////////////////////////////////////////////

  // Create ratio histograms R_A, R_B, R_C from projections
  TH1D *hist_RA = (TH1D *)hist_pmpm_1D_for_RA->Clone("hist_RA");
  hist_RA->SetTitle("R_{A}(t_{1});t_{1} [#tau_{S}];R_{A}(t_{1}) [-]");
  hist_RA->Divide(hist_00pm_1D);

  TH1D *hist_RB = (TH1D *)hist_pmpm_1D_for_RB->Clone("hist_RB");
  hist_RB->SetTitle("R_{B}(t_{1});t_{1} [#tau_{S}];R_{B}(t_{1}) [-]");
  hist_RB->Divide(hist_pm00_1D);

  TH1D *hist_RC = (TH1D *)hist_RA->Clone("hist_RC");
  hist_RC->SetTitle("R_{C}(t_{1});t_{1} [#tau_{S}];R_{C}(t_{1}) [-]");
  hist_RC->Divide(hist_RB);

  TF1 *const_line = new TF1("const_line", "1 + 6 * [0]", 0, 20);
  const_line->SetParameter(0, PhysicsConstants::Re);
  const_line->SetLineColor(kBlack);
  const_line->SetLineStyle(2);
  const_line->SetLineWidth(2);

  // Rysowanie R_A
  TCanvas *c = new TCanvas("c", "R_A from histogram ratio", 800, 600);
  hist_RA->GetXaxis()->SetRangeUser(0, 20);
  hist_RA->Draw("PE1");
  const_line->Draw("SAME");
  c->SaveAs(dir_tmp_path + "interf_func_RA_draw.pdf");

  // Rysowanie R_B
  TCanvas *cg = new TCanvas("cg", "R_B from histogram ratio", 800, 600);
  hist_RB->GetXaxis()->SetRangeUser(0, 20);
  hist_RB->Draw("PE1");
  const_line->Draw("SAME");
  cg->SaveAs(dir_tmp_path + "interf_func_RB_draw.pdf");

  // Rysowanie R_C
  TCanvas *c1 = new TCanvas("c1", "R_C from histogram ratio", 800, 600);
  hist_RC->GetXaxis()->SetRangeUser(0, 20);
  hist_RC->Draw("PE1");

  TF1 *const_line2 = new TF1("const_line2", "1 + 6 * [0]", 0, 20);
  const_line2->SetParameter(0, PhysicsConstants::Re);
  const_line2->SetLineColor(kBlack);
  const_line2->SetLineStyle(2);
  const_line2->SetLineWidth(2);

  TF1 *const_line3 = new TF1("const_line3", "1 - 6 * [0]", 0, 20);
  const_line3->SetParameter(0, PhysicsConstants::Re);
  const_line3->SetLineColor(kBlack);
  const_line3->SetLineStyle(2);
  const_line3->SetLineWidth(2);

  const_line2->Draw("SAME");
  const_line3->Draw("SAME");

  c1->SaveAs(dir_tmp_path + "interf_func_RC_double_draw.pdf");

  // Rename temporary directory to final name
  boost::filesystem::path tmp_path(dir_tmp_path.Data());
  boost::filesystem::path final_path(dir_final_path.Data());
  
  // Remove final directory if it exists
  if (boost::filesystem::exists(final_path))
  {
    boost::filesystem::remove_all(final_path);
  }
  
  // Rename tmp to final
  boost::filesystem::rename(tmp_path, final_path);
  
  std::cout << "Results saved to: " << dir_final_path << std::endl;
}