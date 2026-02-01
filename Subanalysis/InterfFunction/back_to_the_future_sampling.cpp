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
#include <TBranch.h>
#include <TTree.h>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <interf_function.h>
#include <const.h>

Double_t evTick = 10000.0;

void createDirIfNotExists(const TString &path)
{
  boost::filesystem::path dirPath(path.Data());

  if (!boost::filesystem::exists(dirPath))
    boost::filesystem::create_directories(dirPath);
}

void GenerateHitOrMissSample(TH2D *hist, TF2 *func, Long64_t nSamplesTarget, Double_t &t1Hit, Double_t &t2Hit, TTree *tree)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Double_t t1min, t1max, t2min, t2max;
  func->GetRange(t1min, t2min, t1max, t2max);
  Double_t maxFuncVal = func->GetMaximum();

  Long64_t nSamples = 0;
  TRandom3 randGen(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    Double_t t1 = randGen.Uniform(t1min, t1max);
    Double_t t2 = randGen.Uniform(t2min, t2max);

    // Evaluate function value at (t1, t2)
    Double_t funcValue = func->Eval(t1, t2);
    // Generate a random number for hit-or-miss
    Double_t randValue = randGen.Uniform(0, maxFuncVal);

    // Accept or reject the sample
    if (randValue <= funcValue)
    {
      hist->Fill(t1, t2);
      nSamples++;

      t1Hit = t1;
      t2Hit = t2;

      tree->Fill();

      if (nSamples % Int_t(evTick) == 0)
        ++progress;
    }
  }
}

double draw_from_mixed(double ts, double tl)
{
  if (gRandom->Uniform() < 0.5)
  {
    return gRandom->Exp(ts); // Losuj z KS
  }
  else
  {
    return gRandom->Exp(tl); // Losuj z KL
  }
}

void GenerateImportanceSample(TH2D *hist, TF2 *func, Long64_t nSamplesTarget)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Double_t t1min, t1max, t2min, t2max;
  func->GetRange(t1min, t2min, t1max, t2max);
  Double_t maxFuncVal = func->GetMaximum();

  Long64_t nSamples = 0;
  TRandom3 randGen(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    Double_t t1 = draw_from_mixed(1, PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT);
    Double_t t2 = draw_from_mixed(1, PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT);

    if (t1 < t1min || t1 > t1max || t2 < t2min || t2 > t2max)
      continue;

    auto pdf = [&](double t)
    {
      Double_t tau_L = PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT;

      return 0.5 * (exp(-t)) + 0.5 * (1.0 / tau_L * exp(-t / tau_L));
    };

    // Evaluate function value at (t1, t2)
    Double_t funcValue = func->Eval(t1, t2);
    Double_t weight = pdf(t1) * pdf(t2); // Jacobian for exponential sampling
    // Generate a random number for hit-or-miss
    Double_t M = 1000.0;

    Double_t randValue = randGen.Uniform(0, M);

    // Accept or reject the sample
    if (randValue <= funcValue / weight)
    {
      hist->Fill(t1, t2);
      nSamples++;

      if (nSamples % Int_t(evTick) == 0)
        ++progress;
    }
  }
}

void GenerateRandom2(TH2D *hist, TF2 *func, Long64_t nSamplesTarget)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Int_t granularity = 1000;

  func->SetNpx(granularity);
  func->SetNpy(granularity);

  Long64_t nSamples = 0;
  Double_t t1, t2;

  TRandom3 *randGen = new TRandom3(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    func->GetRandom2(t1, t2, randGen);
    hist->Fill(t1, t2); // Weight to account for equal probabilities
    nSamples++;

    if (nSamples % Int_t(evTick) == 0)
      ++progress;
  }
}

void GenerateRandom21D(TH1 *hist, TF2 *func, Long64_t nSamplesTarget)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Int_t granularity = 1000;

  func->SetNpx(granularity);
  func->SetNpy(granularity);

  Long64_t nSamples = 0;
  Double_t t1, t2;

  TRandom3 *randGen = new TRandom3(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    func->GetRandom2(t1, t2, randGen);
    hist->Fill(t1); // Weight to account for equal probabilities
    nSamples++;

    if (nSamples % Int_t(evTick) == 0)
      ++progress;
  }
}

static TRandom3 gRand(0);

void generate_t1_t2(double &t1, double &t2, bool equalProb)
{
  Double_t weightA = PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_kl_pippim,
           weightB = PhysicsConstants::br_ks_pippim * PhysicsConstants::br_kl_pi0pi0,
           totalWeight = weightA + weightB,
           probA = weightA / totalWeight;

  if (equalProb)
    probA = 0.5;

  Double_t u = gRand.Uniform();
  Double_t u1 = gRand.Uniform(), u2 = gRand.Uniform();

  if (u < probA)
  {
    t1 = -TMath::Log(1 - u1);
    t2 = -PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT * TMath::Log(1 - u2);
  }
  else
  {
    t1 = -PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT * TMath::Log(1 - u1);
    t2 = -TMath::Log(1 - u2);
  }
}

int main(int argc, char *argv[])
{
  KLOE::setGlobalStyle();

  std::cout << "Do You want to set custom parameters? (1 - yes, 0 - no): ";
  int customParams = 0;
  std::cin >> customParams;

  Double_t reParam = PhysicsConstants::Re;
  Double_t imParam = PhysicsConstants::Im_nonCPT;

  if (customParams)
  {
    std::cout << "Set Re(epsilon'/epsilon) (default " << PhysicsConstants::Re << "): ";
    std::cin >> reParam;

    std::cout << "Set Im(epsilon'/epsilon) (default " << PhysicsConstants::Im_nonCPT << "): ";
    std::cin >> imParam;
  }

  // Sampling functions to create 2D histograms
  const Float_t sigmaT = 1; // Time resolution in tau_S units
  const Float_t t1Min = 0.0, t1Max = 300.0;
  const Float_t t2Min = 0.0, t2Max = 300.0;
  const Int_t nBinst1 = (t1Max - t1Min) / sigmaT, nBinst2 = (t2Max - t2Min) / sigmaT;

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, t1Min, t1Max, t2Min, t2Max, 2);
  func_00pm->SetParameters(reParam, imParam);

  TF2 *func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, t1Min, t1Max, t2Min, t2Max, 2);
  func_pm00->SetParameters(reParam, imParam);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, t1Min, t1Max, t2Min, t2Max, 2);
  func_pmpm->SetParameters(reParam, imParam);

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

  // Creation of the plot directory if it does not exist
  // Creation of the plot directory with tmp suffix
  TString plot_dir = "theoretical_plots/";
  TString dir_name = "integralLimit_" + TString::Format("%.0f", maxIntegral) + "_samples_" + TString::Format("%lld", nSamplesTarget) + "_Re_" + TString::Format("%.5f", reParam) + "_Im_" + TString::Format("%.5f", imParam);
  TString dir_tmp_path = Paths::img_dir + plot_dir + dir_name + "_tmp/";
  TString dir_final_path = Paths::img_dir + plot_dir + dir_name + "/";
  createDirIfNotExists(dir_tmp_path);

  Double_t tne, tch, t1, t2;

  TFile *outputFile = new TFile(dir_tmp_path + "sampling_results_" + TString::Format("%.0f", maxIntegral) + ".root", "RECREATE");
  outputFile->cd();

  TTree *tree_00pm = new TTree("00pm", "Tree with sampled events");
  TTree *tree_pmpm = new TTree("pmpm", "Tree with sampled events");

  TBranch *branch_t1_00pm = tree_00pm->Branch("tne", &tne, "tne/D");
  TBranch *branch_t2_00pm = tree_00pm->Branch("tch", &tch, "tch/D");

  TBranch *branch_t1_pmpm = tree_pmpm->Branch("t1", &t1, "t1/D");
  TBranch *branch_t2_pmpm = tree_pmpm->Branch("t2", &t2, "t2/D");

  std::cout << "Which method of sampling do You want to use? (1 - Hit or Miss, 2 - GetRandom2, 3 - with weights): ";
  int samplingMethod = 3;
  std::cin >> samplingMethod;

  switch (samplingMethod)
  {
  case 1:
  {
    GenerateHitOrMissSample(hist_00pm, func_00pm, nSamplesTarget, tne, tch, tree_00pm);
    GenerateHitOrMissSample(hist_pm00, func_pm00, nSamplesTarget, tch, tne, tree_00pm);
    GenerateHitOrMissSample(hist_pmpm, func_pmpm, nSamplesTarget, t1, t2, tree_pmpm);
    break;
  }
  case 2:
  {
    GenerateRandom2(hist_00pm, func_00pm, nSamplesTarget);
    GenerateRandom2(hist_pm00, func_pm00, nSamplesTarget);
    GenerateRandom2(hist_pmpm, func_pmpm, nSamplesTarget);

    break;
  }
  case 3:
  {
    boost::progress_display progress00pm(nSamplesTarget / evTick);

    Long64_t nSamples = 0;
    while (nSamples < nSamplesTarget)
    {
      generate_t1_t2(tne, tch, true);

      if (tne <= maxIntegral && tch <= maxIntegral)
      {
        hist_00pm->Fill(tne, tch, func_00pm->Eval(tne, tch));
        hist_pm00->Fill(tch, tne, func_pm00->Eval(tch, tne));

        tree_00pm->Fill();

        nSamples++;

        if (nSamples % Int_t(evTick) == 0)
          ++progress00pm;
      }
    }


    boost::progress_display progresspmpm(nSamplesTarget / evTick);
    nSamples = 0;

    while (nSamples < nSamplesTarget)
    {
      generate_t1_t2(t1, t2, true);

      if (t1 <= maxIntegral && t2 <= maxIntegral)
      {
        hist_pmpm->Fill(t1, t2, func_pmpm->Eval(t1, t2));

        tree_pmpm->Fill();

        nSamples++;

        if (nSamples % Int_t(evTick) == 0)
          ++progresspmpm;
      }
    }
    break;
  }
  }

  TCanvas *c_func_00pm = new TCanvas("c_func_00pm", "Interference function 00pm", 800, 600);
  c_func_00pm->SetLogz(1);
  hist_00pm->Draw("COLZ");
  c_func_00pm->SaveAs(dir_tmp_path + "interf_func_00pm.pdf");

  TCanvas *c_func_pm00 = new TCanvas("c_func_pm00", "Interference function pm00", 800, 600);
  c_func_pm00->SetLogz(1);
  hist_pm00->Draw("COLZ");
  c_func_pm00->SaveAs(dir_tmp_path + "interf_func_pm00.pdf");

  TCanvas *c_func_pmpm = new TCanvas("c_func_pmpm", "Interference function pmpm", 800, 600);
  c_func_pmpm->SetLogz(1);
  hist_pmpm->Draw("COLZ");
  c_func_pmpm->SaveAs(dir_tmp_path + "interf_func_pmpm.pdf");

  // Create 1D projections from 2D histograms with cuts
  // For projection we integrate over y-axis up to a certain limit (par[1] or par[0])

  auto create_projection_y_cut = [](TH2D *hist, Double_t t_max) -> TH1D *
  {
    TRandom3 randGen(0);
    TString name = TString::Format("%s_projx_tmax%.0f_%d", hist->GetName(), t_max, randGen.Integer(400));
    return hist->ProjectionX(name, -1, hist->GetYaxis()->FindBin(t_max) + 1);
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
  hist_pm00_1D->Draw("PE1");
  c_pm00->SaveAs(dir_tmp_path + "interf_func_pm001D_draw.pdf");

  ///////////////////////////////////////////////////////////////////

  TCanvas *c_00pm = new TCanvas("c_00pm", "Projection of 00pm", 800, 600);
  hist_00pm_1D->Draw("PE1");
  c_00pm->SaveAs(dir_tmp_path + "interf_func_00pm1D_draw.pdf");

  ////////////////////////////////////////////////////////////////////

  TCanvas *c_pmpm = new TCanvas("c_pmpm", "Projection of pmpm", 800, 600);
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

  hist_00pm->Write();
  hist_pm00->Write();
  hist_pmpm->Write();

  hist_00pm_1D->Write();
  hist_pm00_1D->Write();
  hist_pmpm_1D_for_RA->Write();
  hist_pmpm_1D_for_RB->Write();

  hist_RA->Write();
  hist_RB->Write();
  hist_RC->Write();

  tree_00pm->Write();
  tree_pmpm->Write();

  outputFile->Close();

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