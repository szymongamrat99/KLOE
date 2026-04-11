#define signal_vs_bcg_v2_cxx
// The class definition in signal_vs_bcg_v2.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("signal_vs_bcg_v2.C")
// root> T->Process("signal_vs_bcg_v2.C","some options")
// root> T->Process("signal_vs_bcg_v2.C+")
//

#include "signal_vs_bcg_v2.h"
#include <kloe_class.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <TF1.h>
#include <TripleGaussFitter.h>
#include <DoublGaussFitter.h>
#include <NeutralReconstruction.h>
#include <triple_gaus.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <chi2_dist.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TRegexp.h>
#include <GeneratedVariables.h>

#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <TGraphAsymmErrors.h>

#include <interf_function.h>
#include <TEfficiency.h>
#include <const.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <fstream>
#include <ctime>

namespace KH = KLOE::Histograms;

std::map<TString, Float_t> channLumi = {
    {"Data", 78799},
    {"Signal", 1549784},
    {"Regeneration", 1574577},
    {"Omega", 1625861},
    {"3pi0", 230874},
    {"Semileptonic", 838182},
    {"Other", 343853}};

std::map<TString, Float_t> channFactor;
std::map<TString, Int_t> channEventsTotal, channEventsCut;

Int_t signal_num = 0, signal_tot = 0, tot_events = 0, bkg_tot = 0, signal_wo_err = 0;

std::map<TString, TCanvas *>
    canvas;

std::map<TString, std::map<TString, TCanvas *>>
    canvas2D;

std::map<TString, TCanvas *>
    canvasProfiles;

std::map<TString, TProfile *>
    profiles;

std::map<TString, std::map<TString, TH1 *>>
    histsReconstructed,
    histsFittedSignal;

std::map<TString, std::map<TString, TH2 *>>
    hists2DReconstructed,
    hists2DFittedSignal;

TCanvas *canvaEff, *canvaPurity;

TH1 *deltaTSignalTot, *deltaTTot;

KLOE::TripleGaussFitter *fitter;
KLOE::DoublGaussFitter *doubleFitter;

KLOE::pm00 Obj;
KLOE::NeutralReconstruction *neutReconst;

TEfficiency *efficiency;

TH1 *histCounts;

TF1 *chi2DistFunc;

struct ScenarioCounters
{
  ScenarioCounters() : signal_num(0), signal_tot(0), bkg_tot(0), passed_events(0), sel_mctruth_m1(0), sel_mctruth_0(0), sel_mctruth_1(0), sel_mctruth_0_before_cut(0), sel_mctruth_1_before_cut(0)
  {
    for (const auto &chann : KLOE::channName)
    {
      sel_mctruth_by_channel[(std::string)chann.second] = 0;
      sel_mctruth_by_channel_before_cut[(std::string)chann.second] = 0;
    }
  }

  Int_t signal_num = 0;
  Int_t signal_tot = 0;
  Int_t bkg_tot = 0;
  Int_t passed_events = 0;
  Int_t sel_mctruth_m1 = 0;
  Int_t sel_mctruth_0 = 0;
  Int_t sel_mctruth_1 = 0;
  Int_t sel_mctruth_0_before_cut = 0;
  Int_t sel_mctruth_1_before_cut = 0;

  std::unordered_map<std::string, Int_t> sel_mctruth_by_channel;
  std::unordered_map<std::string, Int_t> sel_mctruth_by_channel_before_cut;
};

struct ScenarioHistogramState
{
  TString folderPath;
  std::map<TString, std::map<TString, TH1 *>> histsReconstructed;
  std::map<TString, std::map<TString, TH1 *>> histsFittedSignal;
  std::map<TString, std::map<TString, TH2 *>> hists2DFittedSignal;
  std::map<TString, Int_t> channEventsTotal;
  std::map<TString, Int_t> channEventsCut;
  TH1 *deltaTSignalTot = nullptr;
  TH1 *deltaTTot = nullptr;
  TH1 *histCounts = nullptr;
  Int_t overflow = 0;
  Int_t cutPassed = 0;
  Int_t cutNPassed = 0;
};

std::vector<TString> g_activeScenarios;
std::map<TString, ScenarioCounters> g_scenarioCounters;
std::map<TString, ScenarioHistogramState> g_secondaryScenarioStates;
TString g_rootOutputFileName = "signal_vs_bcg_cache.root";

namespace
{
  void DrawEfficiencyWithFixedYAxis(TEfficiency *eff,
                                    TH1 *xRefHist,
                                    Double_t yMin,
                                    Double_t yMax,
                                    const TString &yTitle)
  {
    if (!eff || !xRefHist || !gPad)
      return;

    const Double_t xMin = xRefHist->GetXaxis()->GetXmin();
    const Double_t xMax = xRefHist->GetXaxis()->GetXmax();

    TH1 *frame = gPad->DrawFrame(xMin, yMin, xMax, yMax);
    if (frame)
    {
      frame->GetXaxis()->SetTitle(xRefHist->GetXaxis()->GetTitle());
      frame->GetYaxis()->SetTitle(yTitle);
    }

    eff->Draw("P SAME");

    gPad->Modified();
    gPad->Update();
  }
}

// Central place for tuning analysis cuts.
namespace CutDefs
{
  constexpr Double_t simonaBadClusSlope = -10.0 / 9.0;
  constexpr Double_t simonaBadClusDeltaTResMax = -3.0;
  constexpr Double_t simonaBadClusDeltaTShift = 2.0;

  constexpr Double_t simonaDeltaPhiCenter = 3.110;
  constexpr Double_t simonaDeltaPhiSigma = 0.135;
  constexpr Double_t simonaDeltaPhiNSigma = 2.0;
  constexpr Double_t simonaChi2Max = 30.0;

  constexpr Double_t oldCutsChi2Max = 40.0;
  constexpr Double_t oldCutsCombinedMassPi0Max = 35.0;
  constexpr Double_t oldCutsMassKchWindow = 1.2;
  constexpr Double_t oldCutsMassKneWindow = 76.0;
  constexpr Double_t oldCutsQmissMax = 3.75;
  constexpr Double_t oldCutsTrcSumMin = -1.0;
  constexpr Double_t oldCutsOpeningCosMin = -0.8;

  constexpr Double_t omegaRadiusLimit = 1.0;
  constexpr Double_t omegaZdistLimit = 0.6;
  constexpr Double_t omegaT0Center = 155.658;
  constexpr Double_t omegaT0Sigma = 5.691;
  constexpr Double_t omegaMassCenter = 782.994;
  constexpr Double_t omegaMassSigma = 5.620;
  constexpr Double_t omegaNSigma = 3.0;
  constexpr Double_t omegaLineA = 1.0;
  constexpr Double_t omegaLineB = 625.091;
  constexpr Double_t omegaLineBreal = 14.1421;

  constexpr Double_t kaonPathLimitNeutral = 50.0;
  constexpr Double_t kaonPathLimitCharged = 50.0;
  constexpr Double_t blobDeltaTMin = 75.0;

  constexpr Double_t zFidVol = 1.5;
  constexpr Double_t r00FidVol = 1.5;
  constexpr Double_t rpmFidVol = 2.0;
}

void SaveDeltaTComparisonPlot(const TString &outDir,
                              const TString &scenarioLabel,
                              std::map<TString, std::map<TString, TH1 *>> &histsFit,
                              std::map<TString, std::map<TString, TH1 *>> &histsRec)
{
  TH1 *hFit = histsFit["delta_t"]["Signal"];
  TH1 *hRec = histsRec["delta_t"]["Signal"];
  TH1 *hMC = histsFit["delta_t_MC"]["Signal"];

  if (!hFit || !hRec || !hMC)
    return;
  if (hFit->GetEntries() <= 0 || hRec->GetEntries() <= 0 || hMC->GetEntries() <= 0)
    return;

  TCanvas *cDt = new TCanvas(Form("c_delta_t_compare_%s", scenarioLabel.Data()),
                             Form("delta_t comparison (%s)", scenarioLabel.Data()),
                             750, 750);

  hFit->SetTitle(";#Delta t [#tau_{S}];Events");
  hFit->SetLineColor(kBlue + 1);
  hFit->SetLineWidth(3);
  hRec->SetLineColor(kRed + 1);
  hRec->SetLineWidth(3);
  hMC->SetLineColor(kGreen + 2);
  hMC->SetLineStyle(2);
  hMC->SetLineWidth(3);

  const Double_t yMax = TMath::Max(hFit->GetMaximum(), TMath::Max(hRec->GetMaximum(), hMC->GetMaximum()));
  hFit->GetYaxis()->SetRangeUser(0.0, yMax > 0.0 ? 1.25 * yMax : 1.0);

  hFit->Draw("HIST");
  hRec->Draw("HIST SAME");
  hMC->Draw("HIST SAME");

  TLegend *dtLegend = new TLegend(0.58, 0.72, 0.9, 0.9, "", "NDC");
  dtLegend->SetBorderSize(1);
  dtLegend->SetFillColor(kWhite);
  dtLegend->SetTextSize(0.032);
  dtLegend->AddEntry(hFit, "Fitted #Delta t", "l");
  dtLegend->AddEntry(hRec, "Reconstructed #Delta t", "l");
  dtLegend->AddEntry(hMC, "MC generated #Delta t", "l");
  dtLegend->Draw();

  cDt->SaveAs(Form("%s/delta_t_fitted_vs_reconstructed_vs_mc%s", outDir.Data(), Paths::ext_img.Data()));
  delete cDt;
}

Double_t ComputeSafeRatio(Int_t num, Int_t den)
{
  if (den <= 0)
    return 0.0;
  return static_cast<Double_t>(num) / static_cast<Double_t>(den);
}

TString ParseRootOutputFileName(const TString &option)
{
  TString outName = "signal_vs_bcg_cache.root";

  TString opt = option;
  TObjArray *tokens = opt.Tokenize(";");
  for (Int_t i = 0; i < tokens->GetEntries(); ++i)
  {
    TString token = ((TObjString *)tokens->At(i))->String();
    token = token.Strip(TString::kBoth);
    if (token.IsNull())
      continue;

    TString upper = token;
    upper.ToUpper();
    if (upper.BeginsWith("ROOTFILE=") || upper.BeginsWith("OUTPUT_ROOT=") || upper.BeginsWith("ROOT_OUT="))
    {
      const Ssiz_t eqPos = token.First('=');
      if (eqPos >= 0)
      {
        TString candidate = token(eqPos + 1, token.Length() - eqPos - 1);
        candidate = candidate.Strip(TString::kBoth);
        candidate.ReplaceAll("\"", "");
        if (!candidate.IsNull())
          outName = candidate;
      }
    }
  }
  delete tokens;

  if (!outName.EndsWith(".root"))
    outName += ".root";

  return outName;
}

TString ParseFolderAdditionalName(const TString &option)
{
  TString additionalName = "";

  TString opt = option;
  TObjArray *tokens = opt.Tokenize(";");
  for (Int_t i = 0; i < tokens->GetEntries(); ++i)
  {
    TString token = ((TObjString *)tokens->At(i))->String();
    token = token.Strip(TString::kBoth);
    if (token.IsNull())
      continue;

    TString upper = token;
    upper.ToUpper();
    if (upper.BeginsWith("ADDNAME") || upper.BeginsWith("FOLDER_NAME=") || upper.BeginsWith("DIR_NAME="))
    {
      const Ssiz_t eqPos = token.First('=');
      if (eqPos >= 0)
      {
        TString candidate = token(eqPos + 1, token.Length() - eqPos - 1);
        candidate = candidate.Strip(TString::kBoth);
        candidate.ReplaceAll("\"", "");
        if (!candidate.IsNull())
          additionalName = candidate;
      }
    }
  }
  delete tokens;

  if (!additionalName.IsNull())
    additionalName = "_" + additionalName;

  return additionalName;
}

void SaveHistMap1D(TDirectory *dir, const std::map<TString, std::map<TString, TH1 *>> &hMap)
{
  if (!dir)
    return;

  dir->cd();
  for (const auto &cfgEntry : hMap)
  {
    const TString &cfgName = cfgEntry.first;
    for (const auto &chEntry : cfgEntry.second)
    {
      const TString &chName = chEntry.first;
      TH1 *h = chEntry.second;
      if (!h)
        continue;
      h->Write(Form("%s__%s", cfgName.Data(), chName.Data()), TObject::kOverwrite);
    }
  }
}

void SaveHistMap2D(TDirectory *dir, const std::map<TString, std::map<TString, TH2 *>> &hMap)
{
  if (!dir)
    return;

  dir->cd();
  for (const auto &cfgEntry : hMap)
  {
    const TString &cfgName = cfgEntry.first;
    for (const auto &chEntry : cfgEntry.second)
    {
      const TString &chName = chEntry.first;
      TH2 *h = chEntry.second;
      if (!h)
        continue;
      h->Write(Form("%s__%s", cfgName.Data(), chName.Data()), TObject::kOverwrite);
    }
  }
}

std::vector<TString> ParseScenarioList(const TString &option)
{
  std::vector<TString> scenarios;
  Cuts cuts;

  TString opt = option;
  opt.ToUpper();

  if (opt.IsNull())
  {
    scenarios.push_back("NO_CUTS");
    return scenarios;
  }

  TObjArray *tokens = opt.Tokenize(";");
  for (Int_t i = 0; i < tokens->GetEntries(); ++i)
  {
    TString token = ((TObjString *)tokens->At(i))->String();
    token = token.Strip(TString::kBoth);
    if (!token.IsNull() && (token == "NO_CUTS" || cuts.Contains(token)))
    {
      scenarios.push_back(token);
    }
  }
  delete tokens;

  if (scenarios.empty())
    scenarios.push_back("NO_CUTS");

  return scenarios;
}

Bool_t PassScenario(const TString &scenario,
                    Bool_t shorterKaonPaths,
                    Bool_t oldChi2Cut,
                    Bool_t oldTrcSumCut,
                    Bool_t oldCombinedMassPi0Cut,
                    Bool_t oldMassKchCut,
                    Bool_t oldMassKneCut,
                    Bool_t oldQmissCut,
                    Bool_t oldOpeningAngleCut,
                    Bool_t oldOmegaGeometricalCut,
                    Bool_t oldOmegaFiducialVolume,
                    Bool_t simonaChi2Cut,
                    Bool_t badClusSimona,
                    Bool_t simonaKinCuts,
                    Bool_t simonaDeltaPhiCut,
                    Bool_t omegaMassT0Cut,
                    Bool_t blobCut,
                    Bool_t noBlobCut,
                    Bool_t newChi2Cut,
                    Bool_t newTrkAngleCut,
                    Bool_t newCombinedMassPi0Cut,
                    Bool_t newOmegaGeometricalCut,
                    Bool_t newOmegaMassT0Cut,
                    Bool_t newMassKchCut)
{
  if (scenario == "NO_CUTS")
    return kTRUE;
  if (scenario == "SHORTER_KAON_PATHS")
    return shorterKaonPaths;
  if (scenario == "OLD_CHI2_CUT")
    return oldChi2Cut;
  if (scenario == "OLD_TRCSUM_CUT")
    return oldTrcSumCut;
  if (scenario == "OLD_COMBINED_MASS_PI0_CUT")
    return oldCombinedMassPi0Cut;
  if (scenario == "OLD_MASS_KCH_CUT")
    return oldMassKchCut;
  if (scenario == "OLD_MASS_KNE_CUT")
    return oldMassKneCut;
  if (scenario == "OLD_QMISS_CUT")
    return oldQmissCut;
  if (scenario == "OLD_OPENING_ANGLE_CUT")
    return oldOpeningAngleCut;
  if (scenario == "OLD_OMEGA_GEOMETRICAL_CUT")
    return oldOmegaGeometricalCut;
  if (scenario == "OLD_OMEGA_FIDUCIAL_VOLUME")
    return oldOmegaFiducialVolume;
  if (scenario == "SIMONA_CHI2_CUT")
    return simonaChi2Cut;
  if (scenario == "BAD_CLUS_SIMONA")
    return badClusSimona;
  if (scenario == "SIMONA_KIN_CUTS")
    return simonaKinCuts;
  if (scenario == "SIMONA_ALL_CUTS")
    return simonaDeltaPhiCut;
  if (scenario == "OMEGA_MASS_T0_CUT")
    return omegaMassT0Cut;
  if (scenario == "BLOB")
    return blobCut;
  if (scenario == "NO_BLOB")
    return noBlobCut;
  if (scenario == "NEW_CHI2_CUT")
    return newChi2Cut;
  if (scenario == "NEW_TRK_ANGLE_CUT")
    return newTrkAngleCut;
  if (scenario == "NEW_COMBINED_MASS_PI0_CUT")
    return newCombinedMassPi0Cut;
  if (scenario == "NEW_OMEGA_GEOMETRICAL_CUT")
    return newOmegaGeometricalCut;
  if (scenario == "NEW_OMEGA_MASS_T0_CUT")
    return newOmegaMassT0Cut;
  if (scenario == "NEW_MASS_KCH_CUT")
    return newMassKchCut;
  return kFALSE;
}

Int_t InferChannelFromFileName(const TString &fileName)
{
  TString lowerFileName = fileName;
  // lowerFileName.ToLower();

  // Match explicit channel token in file name to avoid accidental matches.
  if (lowerFileName.Contains("_Signal_"))
    return 1;
  if (lowerFileName.Contains("_Regeneration_"))
    return 2;
  if (lowerFileName.Contains("_Omega_"))
    return 3;
  if (lowerFileName.Contains("_3pi0_"))
    return 4;
  if (lowerFileName.Contains("_Semileptonic_"))
    return 5;
  if (lowerFileName.Contains("_Other_"))
    return 6;

  return 0;
}

// Wczytaj konfiguracje histogramów
auto histogramConfigs1D = KLOE::Histograms::LoadHistogramConfigs1D(Paths::histogramConfig1DPath);
auto histogramConfigs2D = KLOE::Histograms::LoadHistogramConfigs2D(Paths::histogramConfig2DPath);
///////////////////////////////////////////////////////////////////////////////

TString g_additionalName = "";

void signal_vs_bcg_v2::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  ROOT::EnableImplicitMT();

  const char *datestamp = Obj.getCurrentTimestamp().c_str();

  // Logger setup
  ErrorHandling::ErrorLogs logger("log/");

  // Creation of output folder with cut names
  /////////////////////////////////////////
  TString option = GetOption();
  g_activeScenarios = ParseScenarioList(option);
  g_rootOutputFileName = ParseRootOutputFileName(option);
  g_scenarioCounters.clear();
  g_secondaryScenarioStates.clear();
  for (const auto &scenario : g_activeScenarios)
    g_scenarioCounters[scenario] = ScenarioCounters{};

  g_additionalName = ParseFolderAdditionalName(option);

  // Primary scenario keeps the legacy histogram/canvas workflow.
  fOption = g_activeScenarios.front();
  folderPath = Form("img/%s%s/%s", datestamp, g_additionalName.Data(), fOption.Data());
  FolderManagement(folderPath);
  std::cout << "ROOT output file: " << g_rootOutputFileName << std::endl;
  /////////////////////////////////////////

  chi2DistFunc = new TF1("chi2DistFunc", chi2dist, 0, 1E5, 2); // 2 parameters

  KLOE::setGlobalStyle();
  Utils::InitializeVariables(logger);

  fitter = new KLOE::TripleGaussFitter();
  doubleFitter = new KLOE::DoublGaussFitter();
  neutReconst = new KLOE::NeutralReconstruction(4, logger);

  // Create canvases
  for (const auto &config : histogramConfigs1D)
  {
    canvas[config.first] = new TCanvas(Form("c_%s", config.first.Data()),
                                       Form("Canvas for %s", config.first.Data()), 750, 750);
  }

  for (const auto &config : histogramConfigs2D)
  {
    canvasProfiles[config.first] = new TCanvas(Form("cProfiles_%s", config.first.Data()),
                                               Form("Canvas profiles for %s", config.first.Data()), 750, 750);
  }

  // Create canvases
  for (const auto &config : histogramConfigs2D)
  {
    for (const auto &name : KLOE::channName)
    {
      canvas2D[config.first][name.second] = new TCanvas(Form("c_%s_%s", config.first.Data(), name.second.Data()), Form("Canvas for %s (%s)", config.first.Data(), name.second.Data()), 750, 750);
    }
  }

  // Create histograms
  for (const auto &config : histogramConfigs1D)
  {
    for (const auto &name : KLOE::channName)
    {
      TString nameRec = Form("h_rec_%s_%s", config.first.Data(), name.second.Data());
      TString nameFit = Form("h_fit_%s_%s", config.first.Data(), name.second.Data());

      histsReconstructed[config.first][name.second] = KH::CreateHist1D(nameRec, config.second);
      histsFittedSignal[config.first][name.second] = KH::CreateHist1D(nameFit, config.second);
    }
  }

  // Create histograms for efficiency
  canvaEff = new TCanvas("Efficiency", "Efficiency", 800, 800);
  deltaTSignalTot = (TH1F *)histsReconstructed["delta_t_MC"]["Signal"]->Clone("EfficiencyHistTot");
  deltaTSignalTot->Reset();

  histCounts = (TH1F *)histsReconstructed["chi2_signalKinFit"]["Signal"]->Clone("CountsHist");
  histCounts->Reset();
  //////////////////////////////////////////////////////////////////////////////////////////////

  // Create histograms for Purity
  canvaPurity = new TCanvas("Purity", "Purity", 800, 800);
  deltaTTot = (TH1F *)histsFittedSignal["delta_t"]["Signal"]->Clone("PurityHistTot");
  deltaTTot->Reset();
  //////////////////////////////////////////////////////////////////////////////////////////////

  // Create histograms 2D
  for (const auto &config : histogramConfigs2D)
  {
    for (const auto &name : KLOE::channName)
    {
      TString nameRec = Form("h_rec2D_%s_%s", config.first.Data(), name.second.Data());
      TString nameFit = Form("h_fit2D_%s_%s", config.first.Data(), name.second.Data());

      hists2DReconstructed[config.first][name.second] = KH::CreateHist2D(nameRec, config.second);
      hists2DFittedSignal[config.first][name.second] = KH::CreateHist2D(nameFit, config.second);
    }
  }

  // Prepare independent histogram containers for non-primary scenarios.
  for (const auto &scenario : g_activeScenarios)
  {
    if (scenario == fOption)
      continue;

    ScenarioHistogramState state;
    state.folderPath = Form("img/%s%s/%s", datestamp, g_additionalName.Data(), scenario.Data());
    FolderManagement(state.folderPath);

    for (const auto &config : histogramConfigs1D)
    {
      for (const auto &name : KLOE::channName)
      {
        TH1 *srcRec = histsReconstructed[config.first][name.second];
        TH1 *srcFit = histsFittedSignal[config.first][name.second];

        state.histsReconstructed[config.first][name.second] =
            (TH1 *)srcRec->Clone(Form("%s__%s", srcRec->GetName(), scenario.Data()));
        state.histsReconstructed[config.first][name.second]->Reset();

        state.histsFittedSignal[config.first][name.second] =
            (TH1 *)srcFit->Clone(Form("%s__%s", srcFit->GetName(), scenario.Data()));
        state.histsFittedSignal[config.first][name.second]->Reset();
      }
    }

    for (const auto &config : histogramConfigs2D)
    {
      for (const auto &name : KLOE::channName)
      {
        TH2 *src2D = hists2DFittedSignal[config.first][name.second];
        state.hists2DFittedSignal[config.first][name.second] =
            (TH2 *)src2D->Clone(Form("%s__%s", src2D->GetName(), scenario.Data()));
        state.hists2DFittedSignal[config.first][name.second]->Reset();
      }
    }

    state.deltaTSignalTot = (TH1 *)deltaTSignalTot->Clone(Form("%s__%s", deltaTSignalTot->GetName(), scenario.Data()));
    state.deltaTSignalTot->Reset();
    state.deltaTTot = (TH1 *)deltaTTot->Clone(Form("%s__%s", deltaTTot->GetName(), scenario.Data()));
    state.deltaTTot->Reset();
    state.histCounts = (TH1 *)histCounts->Clone(Form("%s__%s", histCounts->GetName(), scenario.Data()));
    state.histCounts->Reset();

    for (const auto &lumi : channLumi)
    {
      state.channEventsTotal[lumi.first] = 0;
      state.channEventsCut[lumi.first] = 0;
    }

    g_secondaryScenarioStates[scenario] = std::move(state);
  }
}

void signal_vs_bcg_v2::SlaveBegin(TTree *tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  fChain = tree;

  TString maxKey;
  Int_t maxValue = 0;

  for (const auto &lumi : channLumi)
  {
    if (lumi.second > maxValue)
    {
      maxValue = lumi.second;
      maxKey = lumi.first;
    }
  }

  for (const auto &lumi : channLumi)
  {
    channFactor[lumi.first] = channLumi.at(maxKey) / lumi.second;

    channEventsTotal[lumi.first] = 0;
    channEventsCut[lumi.first] = 0;
  }

  // ✅ DEBUG - wyświetl dostępne kanały
  std::cout << "\n=== Available keys in KLOE::channName ===" << std::endl;
  for (const auto &ch : KLOE::channName)
  {
    std::cout << "  Key: " << ch.first << " => \"" << ch.second << "\"" << std::endl;
  }
  std::cout << "=========================================\n"
            << std::endl;
}

Int_t overflow = 0, cutPassed = 0, cutNPassed = 0;

Bool_t signal_vs_bcg_v2::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  GeneratedVariables genVarClassifier;

  Double_t dataPCA[2];

  TString fileNameTmp;

  Int_t mctruth_int = 0;

  fReader.SetLocalEntry(entry);

  mctruth_int = *mctruth;

  fileNameTmp = (TString)fChain->GetCurrentFile()->GetName();

  if (mctruth_int == 0 && *mcflag == 1)
  {
    const Int_t channelFromName = InferChannelFromFileName(fileNameTmp);
    if (channelFromName > 0)
    {
      mctruth_int = channelFromName;
    }

    mctruth_int = 1;
  }

  Float_t phiv1PhiCM, phivPlusPhiCM, phiv2PhiCM, phivMinusPhiCM, deltaPhiPhiCM;

  TVector3 boostPhi(-*Bpx / *Broots, -*Bpy / *Broots, -*Bpz / *Broots);

  TLorentzVector trk1LAB(trk1[0], trk1[1], trk1[2], trk1[3]), trk1PhiCM(0, 0, 0, 0);
  TLorentzVector trk2LAB(trk2[0], trk2[1], trk2[2], trk2[3]), trk2PhiCM(0, 0, 0, 0);

  Obj.lorentz_transf(boostPhi, trk1LAB, trk1PhiCM);
  Obj.lorentz_transf(boostPhi, trk2LAB, trk2PhiCM);

  phiv1PhiCM = trk1PhiCM.Phi();
  phiv2PhiCM = trk2PhiCM.Phi();

  Float_t KnerecPhotons[9];

  KnerecPhotons[0] = gammaMomTriangle1[0] + gammaMomTriangle2[0] +
                     gammaMomTriangle3[0] + gammaMomTriangle4[0];
  KnerecPhotons[1] = gammaMomTriangle1[1] + gammaMomTriangle2[1] +
                     gammaMomTriangle3[1] + gammaMomTriangle4[1];
  KnerecPhotons[2] = gammaMomTriangle1[2] + gammaMomTriangle2[2] +
                     gammaMomTriangle3[2] + gammaMomTriangle4[2];
  KnerecPhotons[3] = gammaMomTriangle1[3] + gammaMomTriangle2[3] +
                     gammaMomTriangle3[3] + gammaMomTriangle4[3];

  KnerecPhotons[4] = std::sqrt(std::pow(KnerecPhotons[0], 2) + std::pow(KnerecPhotons[1], 2) + std::pow(KnerecPhotons[2], 2));

  Float_t vKchFit = PhysicsConstants::cVel * KchboostFit[4] / KchboostFit[3],
          pathKchFit = std::sqrt(std::pow(KchboostFit[6] - ipFit[0], 2) +
                            std::pow(KchboostFit[7] - ipFit[1], 2) +
                            std::pow(KchboostFit[8] - ipFit[2], 2)),
          RKchFit = std::sqrt(std::pow(KchboostFit[6] - ipFit[0], 2) +
                         std::pow(KchboostFit[7] - ipFit[1], 2)),
          tKchFit = KchboostFit[9] / 0.0895,
          vKneFit = PhysicsConstants::cVel * KnereclorFit[4] / KnereclorFit[3],
          pathKneFit = std::sqrt(std::pow(KnerecFit[6] - ipFit[0], 2) +
                            std::pow(KnerecFit[7] - ipFit[1], 2) +
                            std::pow(KnerecFit[8] - ipFit[2], 2)),
          // std::pow(KnerecFit[8] - ipFit[2], 2)),
      RKneFit = std::sqrt(std::pow(KnerecFit[6] - ipFit[0], 2) +
                     std::pow(KnerecFit[7] - ipFit[1], 2)),
          tKneFit = KnereclorFit[9] / 0.0895,
          vKneMC = PhysicsConstants::cVel * Knemc[4] / Knemc[3],
          vKchMC = PhysicsConstants::cVel * Kchmc[4] / Kchmc[3],
          vKne = PhysicsConstants::cVel * Knerec[4] / Knerec[3],
          pathKne = std::sqrt(std::pow(Knerec[6] - ip[0], 2) +
                         std::pow(Knerec[7] - ip[1], 2) +
                         std::pow(Knerec[8] - ip[2], 2)),
          pathKch = std::sqrt(std::pow(Kchrec[6] - ip[0], 2) +
                         std::pow(Kchrec[7] - ip[1], 2) +
                         std::pow(Kchrec[8] - ip[2], 2)),
          RKne = std::sqrt(std::pow(Knerec[6] - ip[0], 2) +
                      std::pow(Knerec[7] - ip[1], 2)),
          RKch = std::sqrt(std::pow(Kchrec[6] - ip[0], 2) +
                      std::pow(Kchrec[7] - ip[1], 2)),
          tKne = pathKne / (vKne * 0.0895),
          pathKchMC = std::sqrt(std::pow(Kchmc[6] - ipmc[0], 2) +
                           std::pow(Kchmc[7] - ipmc[1], 2)),
          // std::pow(Kchmc[8] - ipmc[2], 2)),
      pathKneMC = std::sqrt(std::pow(Knemc[6] - ipmc[0], 2) +
                       std::pow(Knemc[7] - ipmc[1], 2));
  //  std::pow(Knemc[8] - ipmc[2], 2));

  std::vector<Float_t> kaonMom1 = {Kchboost[0],
                                   Kchboost[1],
                                   Kchboost[2],
                                   Kchboost[3]},
                       kaonMom2 = {KnerecPhotons[0],
                                   KnerecPhotons[1],
                                   KnerecPhotons[2],
                                   KnerecPhotons[3]},
                       kaonPos1 = {Kchboost[6],
                                   Kchboost[7],
                                   Kchboost[8],
                                   Kchboost[9]},
                       kaonPos2 = {Knerec[6],
                                   Knerec[7],
                                   Knerec[8],
                                   Knerec[9]},
                       kaonMomFit1 = {KchboostFit[0],
                                      KchboostFit[1],
                                      KchboostFit[2],
                                      KchboostFit[3]},
                       kaonMomFit2 = {KnerecFit[0],
                                      KnerecFit[1],
                                      KnerecFit[2],
                                      KnerecFit[3]},
                       kaonPosFit1 = {KchboostFit[6],
                                      KchboostFit[7],
                                      KchboostFit[8],
                                      KchboostFit[9]},
                       kaonPosFit2 = {KnerecFit[6],
                                      KnerecFit[7],
                                      KnerecFit[8],
                                      KnerecFit[9]},
                       ipVec = {ip[0], ip[1], ip[2]},
                       ipVecFit = {ipFit[0], ipFit[1], ipFit[2]};

  KLOE::KaonProperTimes propTimes = Obj.CalculateKaonProperTimes(kaonMom1,
                                                                 kaonPos1,
                                                                 kaonMom2,
                                                                 kaonPos2,
                                                                 ipVec);

  KLOE::KaonProperTimes propTimesFit = Obj.CalculateKaonProperTimes(kaonMomFit1,
                                                                    kaonPosFit1,
                                                                    kaonMomFit2,
                                                                    kaonPosFit2,
                                                                    ipVecFit);
  Float_t photon1path = std::sqrt(std::pow(photonFit1[4] - KnerecFit[6], 2) +
                             std::pow(photonFit1[5] - KnerecFit[7], 2) +
                             std::pow(photonFit1[6] - KnerecFit[8], 2)),
          photon2path = std::sqrt(std::pow(photonFit2[4] - KnerecFit[6], 2) +
                             std::pow(photonFit2[5] - KnerecFit[7], 2) +
                             std::pow(photonFit2[6] - KnerecFit[8], 2)),
          photon3path = std::sqrt(std::pow(photonFit3[4] - KnerecFit[6], 2) +
                             std::pow(photonFit3[5] - KnerecFit[7], 2) +
                             std::pow(photonFit3[6] - KnerecFit[8], 2)),
          photon4path = std::sqrt(std::pow(photonFit4[4] - KnerecFit[6], 2) +
                             std::pow(photonFit4[5] - KnerecFit[7], 2) +
                             std::pow(photonFit4[6] - KnerecFit[8], 2));

  Float_t trc1Fit = photonFit1[7] - photon1path / PhysicsConstants::cVel - tKneFit * 0.0895,
          trc2Fit = photonFit2[7] - photon2path / PhysicsConstants::cVel - tKneFit * 0.0895,
          trc3Fit = photonFit3[7] - photon3path / PhysicsConstants::cVel - tKneFit * 0.0895,
          trc4Fit = photonFit4[7] - photon4path / PhysicsConstants::cVel - tKneFit * 0.0895,
          TrcSumFit = trc1Fit + trc2Fit + trc3Fit + trc4Fit;

  Float_t deltaTfit = propTimesFit.deltaTimeCM,
          deltaT = propTimes.kaon1TimeCM - propTimes.kaon2TimeCM,
          deltaTMC = *KaonChTimeCMMC - *KaonNeTimeCMMC;

  Float_t Pi01FitMean = 134.83240924168854, Pi01FitSigma = 3.4027932219804007,
          Pi02FitMean = 134.87134080668457, Pi02FitSigma = 3.2275879781199186,
          rhoFit = -0.33996592037919965;

  Float_t deltaPi01Fit = (pi01Fit[5] - Pi01FitMean) / Pi01FitSigma,
          deltaPi02Fit = (pi02Fit[5] - Pi02FitMean) / Pi02FitSigma,
          rhoFactor = 1 / (1 - std::pow(rhoFit, 2)),
          pi0MassNorm = std::sqrt(rhoFactor * (std::pow(deltaPi01Fit, 2) + std::pow(deltaPi02Fit, 2) - 2 * rhoFit * deltaPi01Fit * deltaPi02Fit));

  // Rectangular cut on combined mass of two pi0 candidates
  Float_t u = ((pi01Fit[5] - Pi01FitMean) + (pi02Fit[5] - Pi02FitMean)) / std::sqrt(2),
          v = ((pi01Fit[5] - Pi01FitMean) - (pi02Fit[5] - Pi02FitMean)) / std::sqrt(2),
          varu = 0.5 * (std::pow(Pi01FitSigma, 2) + std::pow(Pi02FitSigma, 2) + 2 * rhoFit * Pi01FitSigma * Pi02FitSigma),
          varv = 0.5 * (std::pow(Pi01FitSigma, 2) + std::pow(Pi02FitSigma, 2) - 2 * rhoFit * Pi01FitSigma * Pi02FitSigma),
          sigmau = std::sqrt(varu),
          sigmav = std::sqrt(varv);

  Float_t combinedMassPi0Fit = std::sqrt(std::pow(pi01Fit[5] - PhysicsConstants::mPi0, 2) +
                                    std::pow(pi02Fit[5] - PhysicsConstants::mPi0, 2)),
          combinedMassPi0 = std::sqrt(std::pow(pi01[5] - PhysicsConstants::mPi0, 2) +
                                 std::pow(pi02[5] - PhysicsConstants::mPi0, 2));

  Float_t kaonChPath = std::sqrt(std::pow(Kchrec[6] - ip[0], 2) +
                            std::pow(Kchrec[7] - ip[1], 2) +
                            std::pow(Kchrec[8] - ip[2], 2)),
          kaonChMom[3] = {(Kchrec[6] - ip[0]) / kaonChPath * Kchboost[4],
                          (Kchrec[7] - ip[1]) / kaonChPath * Kchboost[4],
                          (Kchrec[8] - ip[2]) / kaonChPath * Kchboost[4]};

  TVector3 boostNeutralKaon(-Knerec[0] / Knerec[3],
                            -Knerec[1] / Knerec[3],
                            -Knerec[2] / Knerec[3]),
      boostChargedKaon(-kaonChMom[0] / Kchboost[3],
                       -kaonChMom[1] / Kchboost[3],
                       -kaonChMom[2] / Kchboost[3]);

  TLorentzVector pi1(trk1[0], trk1[1], trk1[2], trk1[3]);
  TLorentzVector pi2(trk2[0], trk2[1], trk2[2], trk2[3]);

  TLorentzVector pi01Vec(pi01[0], pi01[1], pi01[2], pi01[3]);
  TLorentzVector pi02Vec(pi02[0], pi02[1], pi02[2], pi02[3]);

  pi1.Boost(boostChargedKaon);
  pi2.Boost(boostChargedKaon);
  pi01Vec.Boost(boostNeutralKaon);
  pi02Vec.Boost(boostNeutralKaon);

  Float_t openingAngleCharged = pi1.Angle(pi2.Vect()) * 180.0 / TMath::Pi(),
          openingAngleNeutral = pi01Vec.Angle(pi02Vec.Vect()) * 180.0 / TMath::Pi();

  Float_t QmissFit = std::sqrt(std::pow(KchboostFit[0] - KchrecFit[0], 2) +
                          std::pow(KchboostFit[1] - KchrecFit[1], 2) +
                          std::pow(KchboostFit[2] - KchrecFit[2], 2) +
                          std::pow(KchboostFit[3] - KchrecFit[3], 2));

  Float_t weight = 1.0;

  Double_t gammaS = 1.0; // Wartość gamma_S do ustawienia zakresów
  Double_t gammaL = PhysicsConstants::tau_S_nonCPT / PhysicsConstants::tau_L;

  Double_t *x = new Double_t(*KaonChTimeCMMC - *KaonNeTimeCMMC),
           *par = nullptr;

  if ((mctruth_int == 1 || mctruth_int == 0) && *mcflag == 1)
  {
    weight = interf_function(x, par) / double_exponential(x, par);
  }

  TVector3 z_axis(0., 0., 1.),
      gamma1(gammaMomTriangle1[0], gammaMomTriangle1[1], gammaMomTriangle1[2]),
      gamma2(gammaMomTriangle2[0], gammaMomTriangle2[1], gammaMomTriangle2[2]),
      gamma3(gammaMomTriangle3[0], gammaMomTriangle3[1], gammaMomTriangle3[2]),
      gamma4(gammaMomTriangle4[0], gammaMomTriangle4[1], gammaMomTriangle4[2]);

  Float_t thetaGamma1 = gamma1.Angle(z_axis) * 180.0 / TMath::Pi(),
          thetaGamma2 = gamma2.Angle(z_axis) * 180.0 / TMath::Pi(),
          thetaGamma3 = gamma3.Angle(z_axis) * 180.0 / TMath::Pi(),
          thetaGamma4 = gamma4.Angle(z_axis) * 180.0 / TMath::Pi();

  Bool_t passSemi = std::abs(*Qmiss - 40.73) < 15 && openingAngleNeutral > 160. && openingAngleCharged > 160.;

  Double_t angleLower = 0., angleUpper = 180.;

  Bool_t passGamma = thetaGamma1 > angleLower && thetaGamma1 < angleUpper &&
                     thetaGamma2 > angleLower && thetaGamma2 < angleUpper &&
                     thetaGamma3 > angleLower && thetaGamma3 < angleUpper &&
                     thetaGamma4 > angleLower && thetaGamma4 < angleUpper;

  Float_t deltaPhi,
      deltaPhiMC,
      deltaPhiFit,
      deltaTheta;

  TVector3 trk1Vec(ParamSignalFit[23], ParamSignalFit[24], ParamSignalFit[25]),
      trk2Vec(ParamSignalFit[26], ParamSignalFit[27], ParamSignalFit[28]),
      trk1VecRec(trk1[0], trk1[1], trk1[2]),
      trk2VecRec(trk2[0], trk2[1], trk2[2]);

  Float_t Theta1Fit = trk1Vec.Theta(),
          Theta2Fit = trk2Vec.Theta(),
          Phi1Fit = trk1Vec.Phi(),
          Phi2Fit = trk2Vec.Phi(),
          Theta1Rec = trk1VecRec.Theta(),
          Theta2Rec = trk2VecRec.Theta(),
          Phi1Rec = trk1VecRec.Phi(),
          Phi2Rec = trk2VecRec.Phi();

  deltaPhiPhiCM = std::abs(phiv2PhiCM - phiv1PhiCM);
  deltaPhi = std::abs(Phi2Rec - Phi1Rec);
  deltaTheta = std::abs(Theta2Fit - Theta1Fit);
  deltaPhiFit = std::abs(Phi2Fit - Phi1Fit);

  // Analiza Simony ciecie na phi bad
  Bool_t condGeneral = (deltaTfit - deltaTMC)<CutDefs::simonaBadClusDeltaTResMax,
                                              condLowerLimit = (deltaTfit - deltaTMC)>
                           CutDefs::simonaBadClusSlope * deltaTMC,
         condUpperLimit = (deltaTfit - deltaTMC) < (CutDefs::simonaBadClusSlope * (deltaTMC - CutDefs::simonaBadClusDeltaTShift)),
         badClusSimona = condGeneral && condLowerLimit && condUpperLimit;
  ///////////////////////////////////////////////////////////////////////////////
  // Analiza Simony cięcie na 3 sigma mas
  Bool_t condMassKch = std::abs(Kchrec[5] - 497.605) < 3 * 0.879,
         condMassKne = std::abs(*minv4gam - 488.411) < 3 * 41.293,
         condMassPi01 = std::abs(pi01Fit[5] - 134.840) < 3 * 3.479,
         condMassPi02 = std::abs(pi02Fit[5] - 134.867) < 3 * 3.331;

  ///////////////////////////////////////////////////////////////////////////////

  Float_t radius00 = std::sqrt(std::pow(KnerecFit[6] - ipFit[0], 2) +
                          std::pow(KnerecFit[7] - ipFit[1], 2)),
          radiuspm = std::sqrt(std::pow(KchrecFit[6] - ipFit[0], 2) +
                          std::pow(KchrecFit[7] - ipFit[1], 2)),
          radius00Center = std::sqrt(std::pow(KnerecFit[6], 2) +
                                std::pow(KnerecFit[7], 2)),
          radiuspmCenter = std::sqrt(std::pow(KchrecFit[6], 2) +
                                std::pow(KchrecFit[7], 2)),
          zdist00 = std::abs(KnerecFit[8] - ipFit[2]),
          zdistpm = std::abs(KchrecFit[8] - ipFit[2]),
          path00 = std::sqrt(std::pow(radius00, 2) + std::pow(zdist00, 2)),
          pathpm = std::sqrt(std::pow(radiuspm, 2) + std::pow(zdistpm, 2)),
          path00MC = std::sqrt(std::pow(std::sqrt(std::pow(Knemc[6] - ipmc[0], 2) +
                                   std::pow(Knemc[7] - ipmc[1], 2)),
                              2) +
                          std::pow(std::abs(Knemc[8] - ipmc[2]), 2)),
          pathpmMC = std::sqrt(std::pow(std::sqrt(std::pow(Kchmc[6] - ipmc[0], 2) +
                                   std::pow(Kchmc[7] - ipmc[1], 2)),
                              2) +
                          std::pow(std::abs(Kchmc[8] - ipmc[2]), 2)),
          radius00MC = std::sqrt(std::pow(Knemc[6], 2) +
                            std::pow(Knemc[7], 2)),
          radiuspmMC = std::sqrt(std::pow(Kchmc[6], 2) +
                            std::pow(Kchmc[7], 2)),
          path00MCCenter = std::sqrt(std::pow(std::sqrt(std::pow(Knemc[6], 2) +
                                         std::pow(Knemc[7], 2)),
                                    2) +
                                std::pow(std::abs(Knemc[8]), 2)),
          pathpmMCCenter = std::sqrt(std::pow(std::sqrt(std::pow(Kchmc[6], 2) +
                                         std::pow(Kchmc[7], 2)),
                                    2) +
                                std::pow(std::abs(Kchmc[8]), 2)),

          zdist00MC = std::abs(Knemc[8] - ipmc[2]),
          zdistpmMC = std::abs(Kchmc[8] - ipmc[2]);

  // Calculating everything for the omega-pi0 rejection method - geometrically
  std::array<Float_t, 3> distNeutralCharged = {KchrecClosest[6] - Knerec[6],
                                               KchrecClosest[7] - Knerec[7],
                                               KchrecClosest[8] - Knerec[8]},
                         distNeutralIP = {Knerec[6] - *Bx,
                                          Knerec[7] - *By,
                                          Knerec[8] - KchrecClosest[8]},
                         distChargedIP = {KchrecClosest[6] - *Bx,
                                          KchrecClosest[7] - *By,
                                          KchrecClosest[8] - *Bz};

  Float_t
      rho_pm = std::sqrt(distChargedIP[0] * distChargedIP[0] + distChargedIP[1] * distChargedIP[1]),
      rho_00 = std::sqrt(distNeutralIP[0] * distNeutralIP[0] + distNeutralIP[1] * distNeutralIP[1]),
      rho = std::sqrt(std::pow(rho_pm, 2) + std::pow(rho_00, 2));
  //

  Float_t T0Omega = pi0OmegaFit1[3] - pi0OmegaFit1[5];

  TVector3 phiMeson = {ParamSignalFit[32], ParamSignalFit[33], ParamSignalFit[34]};
  TVector3 KchrecVec = {KchrecFit[0], KchrecFit[1], KchrecFit[2]};
  TVector3 KchrecVecMC = {Kchmc[0], Kchmc[1], Kchmc[2]};
  TVector3 trk1VecMC = {trk1MC[0], trk1MC[1], trk1MC[2]};
  TVector3 trk2VecMC = {trk2MC[0], trk2MC[1], trk2MC[2]};
  TVector3 pi01VecMom = {pi01Fit[0], pi01Fit[1], pi01Fit[2]};
  TVector3 pi02VecMom = {pi02Fit[0], pi02Fit[1], pi02Fit[2]};
  TVector3 KnerecVec = {KnerecFit[0], KnerecFit[1], KnerecFit[2]};

  Double_t phiTrk1Angle = trk1Vec.Angle(KchrecVec) * 180.0 / TMath::Pi(),
           phiTrk2Angle = trk2Vec.Angle(KchrecVec) * 180.0 / TMath::Pi(),
           phipi01Angle = pi01VecMom.Angle(KnerecVec) * 180.0 / TMath::Pi(),
           phipi02Angle = pi02VecMom.Angle(KnerecVec) * 180.0 / TMath::Pi(),
           phiTrk1AngleMC = trk1VecMC.Angle(KchrecVecMC) * 180.0 / TMath::Pi(),
           phiTrk2AngleMC = trk2VecMC.Angle(KchrecVecMC) * 180.0 / TMath::Pi();

  // mcflagCondition
  Bool_t mcflagCondition = (*mcflag == 1 && mctruth_int >= 0) || *mcflag == 0;

  // Geometrical omega-pi0 rejection cuts
  Bool_t
      fiducialVolume = std::sqrt(std::pow(distNeutralCharged[0], 2) + std::pow(distNeutralCharged[1], 2)) < 2.05 && std::abs(distNeutralCharged[2]) < 2.45,
      fiducialVolumeSimona = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

  // Simona Cuts
  Bool_t simonaChi2Cut = *Chi2SignalKinFit <= CutDefs::simonaChi2Max,
         simonaDeltaPhiCut = /*((fiducialVolume && (std::abs(cos(phiTrk2Angle * TMath::Pi() / 180.0)) < 0.8 && cos(phiTrk1Angle * TMath::Pi() / 180.0) < 0.8)) || !fiducialVolume) && simonaChi2Cut,*/ std::abs(deltaPhiFit - CutDefs::simonaDeltaPhiCenter) > CutDefs::simonaDeltaPhiNSigma * CutDefs::simonaDeltaPhiSigma && simonaChi2Cut,
         simonaKinCuts = condMassKch && condMassKne && condMassPi01 && condMassPi02 && simonaDeltaPhiCut,
         simonaPositionLimits = radius00 < CutDefs::omegaRadiusLimit && radiuspm < CutDefs::omegaRadiusLimit &&
                                zdist00 < CutDefs::omegaZdistLimit && zdistpm < CutDefs::omegaZdistLimit,
         omegaMassT0Cut = ((simonaPositionLimits &&
                            !(std::abs(T0Omega - CutDefs::omegaT0Center) < CutDefs::omegaNSigma * CutDefs::omegaT0Sigma &&
                              std::abs(omegaFit[5] - CutDefs::omegaMassCenter) < CutDefs::omegaNSigma * CutDefs::omegaMassSigma &&
                              omegaFit[5] < CutDefs::omegaLineA * T0Omega + CutDefs::omegaLineB + CutDefs::omegaLineBreal &&
                              omegaFit[5] > CutDefs::omegaLineA * T0Omega + CutDefs::omegaLineB - CutDefs::omegaLineBreal)) ||
                           !simonaPositionLimits) &&
                          simonaKinCuts;

  // Old cuts
  Bool_t oldChi2Cut = *Chi2SignalKinFit < CutDefs::oldCutsChi2Max,
         oldTrcSumCut = oldChi2Cut && *TrcSum > CutDefs::oldCutsTrcSumMin,
         oldCombinedMassPi0Cut = oldTrcSumCut && combinedMassPi0Fit < CutDefs::oldCutsCombinedMassPi0Max,
         oldMassKchCut = oldCombinedMassPi0Cut && std::abs(Kchrec[5] - PhysicsConstants::mK0) < CutDefs::oldCutsMassKchWindow,
         oldMassKneCut = oldMassKchCut && std::abs(*minv4gam - PhysicsConstants::mK0) < CutDefs::oldCutsMassKneWindow,
         oldQmissCut = oldMassKneCut && *Qmiss < CutDefs::oldCutsQmissMax,
         oldOpeningAngleCut = oldQmissCut && openingAngleCharged > acos(CutDefs::oldCutsOpeningCosMin),
         omegaPi0RejectionCut = ((rho > 0.8 && fiducialVolume) || !fiducialVolume) && oldOpeningAngleCut;

  

  // New cuts
  Bool_t newChi2Cut = *Chi2SignalKinFit < CutDefs::simonaChi2Max,
         newTrkAngleCut = newChi2Cut && ((fiducialVolumeSimona && (std::abs(cos(phiTrk2Angle * TMath::Pi() / 180.0)) < 0.8 && std::abs(cos(phiTrk1Angle * TMath::Pi() / 180.0)) < 0.8)) || !fiducialVolumeSimona),
         newKchCut = newTrkAngleCut && condMassKch,
         newCombinedMassPi0Cut = newKchCut && (std::abs(u) < 3.0 * sigmau && std::abs(v) < 3.0 * sigmav),
         newOmegaGeometricalCut = newCombinedMassPi0Cut && ((rho > 0.8 && fiducialVolume) || !fiducialVolume),
         newOmegaT0Cut = ((simonaPositionLimits &&
                           !(std::abs(T0Omega - CutDefs::omegaT0Center) < CutDefs::omegaNSigma * CutDefs::omegaT0Sigma &&
                             std::abs(omegaFit[5] - CutDefs::omegaMassCenter) < CutDefs::omegaNSigma * CutDefs::omegaMassSigma &&
                             omegaFit[5] < CutDefs::omegaLineA * T0Omega + CutDefs::omegaLineB + CutDefs::omegaLineBreal &&
                             omegaFit[5] > CutDefs::omegaLineA * T0Omega + CutDefs::omegaLineB - CutDefs::omegaLineBreal)) ||
                          !simonaPositionLimits) &&
                         newCombinedMassPi0Cut;

  // Additional cuts
  Bool_t shorterKaonPaths = pathKch < CutDefs::kaonPathLimitCharged && pathKne<CutDefs::kaonPathLimitNeutral,
                                                                               blobCut = *KaonNeTimeCMBoostTriFit - *KaonNeTimeCMBoostLor>
                                                                           CutDefs::blobDeltaTMin,
         noBlobCut = *KaonNeTimeCMBoostTriFit - *KaonNeTimeCMBoostLor <= CutDefs::blobDeltaTMin;

  if ((mctruth_int == 1) && *mcflag == 1)
  {
    deltaTSignalTot->Fill(deltaTMC, weight);
  }

  if (mctruth_int >= 0 && *mcflag == 1)
    channEventsTotal[KLOE::channName.at(mctruth_int)]++;

  for (const auto &scenario : g_activeScenarios)
  {
    auto &cnt = g_scenarioCounters[scenario];

    const Bool_t passScenario = PassScenario(scenario,
                                             shorterKaonPaths,
                                             oldChi2Cut,
                                             oldTrcSumCut,
                                             oldCombinedMassPi0Cut,
                                             oldMassKchCut,
                                             oldMassKneCut,
                                             oldQmissCut,
                                             oldOpeningAngleCut,
                                             omegaPi0RejectionCut,
                                             fiducialVolume,
                                             simonaChi2Cut,
                                             badClusSimona,
                                             simonaKinCuts,
                                             simonaDeltaPhiCut,
                                             omegaMassT0Cut,
                                             blobCut,
                                             noBlobCut,
                                             newChi2Cut,
                                             newTrkAngleCut,
                                             newCombinedMassPi0Cut,
                                             newOmegaGeometricalCut,
                                             newOmegaT0Cut,
                                             newKchCut);

    // For the requested preselection/selection/total efficiencies,
    // count MC truth classes after scenario cut and before mcflagCondition filtering.
    if (*mcflag == 1)
    {
      if (mctruth_int == 1 || mctruth_int == -1 || mctruth_int == 0)
        cnt.signal_tot++; // All signal events, along with errors

      if (mctruth_int == -1)
        cnt.sel_mctruth_m1++; // Error events
    }

    if (!mcflagCondition)
      continue;

    cnt.passed_events++; // All events, no matter the channel, before any cut, but without errors

    if (mctruth_int == 0 && *mcflag == 1)
      cnt.sel_mctruth_0_before_cut++; // Events which did not pass initial cuts

    if (mctruth_int == 1)
      cnt.sel_mctruth_1_before_cut++; // Events which passed initial cuts

    if (mctruth_int > 0)
    {
      cnt.sel_mctruth_by_channel_before_cut[(std::string)KLOE::channName.at(mctruth_int)]++; // Events in all channels before any cuts, without errors
    }

    if (!passScenario)
      continue;

    if (mctruth_int == 0)
      cnt.sel_mctruth_0++; // Events which did not pass initial cuts after new cut

    if (mctruth_int == 1)
      cnt.sel_mctruth_1++; // Events which passed initial cuts after new cut

    if (mctruth_int > 0)
    {
      cnt.sel_mctruth_by_channel[(std::string)KLOE::channName.at(mctruth_int)]++; // Events in all channels after cut, without errors
    }

    if (mctruth_int > 1)
      cnt.bkg_tot++; // Background which passed the cut
  }

  const Bool_t primaryPass = PassScenario(fOption,
                                          shorterKaonPaths,
                                          oldChi2Cut,
                                          oldTrcSumCut,
                                          oldCombinedMassPi0Cut,
                                          oldMassKchCut,
                                          oldMassKneCut,
                                          oldQmissCut,
                                          oldOpeningAngleCut,
                                          omegaPi0RejectionCut,
                                          fiducialVolume,
                                          simonaChi2Cut,
                                          badClusSimona,
                                          simonaKinCuts,
                                          simonaDeltaPhiCut,
                                          omegaMassT0Cut,
                                          blobCut,
                                          noBlobCut,
                                          newChi2Cut,
                                          newTrkAngleCut,
                                          newCombinedMassPi0Cut,
                                          newOmegaGeometricalCut,
                                          newOmegaT0Cut,
                                          newKchCut);

  Bool_t corrPosLimit = (radius00 < 1.5 && radiuspm < 1.5 && zdist00 < 1.0 && zdistpm < 1.0);
  Bool_t phivLimit = std::abs(deltaPhiFit - 3.110) < 2 * 0.135;
  Bool_t condResCorr = (corrPosLimit && !(cos(phiTrk1Angle * TMath::Pi() / 180.0) > 0.9 || cos(phiTrk2Angle * TMath::Pi() / 180.0) > 0.9)) || !corrPosLimit;

  auto fillAcceptedEvent = [&](std::map<TString, std::map<TString, TH1 *>> &targetHistsReconstructed,
                               std::map<TString, std::map<TString, TH1 *>> &targetHistsFitted,
                               std::map<TString, std::map<TString, TH2 *>> &targetHists2DFitted,
                               std::map<TString, Int_t> &targetEventsCut,
                               Int_t &targetCutPassed,
                               Int_t &targetCutNPassed,
                               Int_t &targetOverflow,
                               TH1 *targetHistCounts,
                               Bool_t addToPCA)
  {
    targetEventsCut[KLOE::channName.at(mctruth_int)]++;

    if ((*mcflag == 1 && mctruth_int >= 0) || *mcflag == 0)
    {
      targetHistsReconstructed["mass_Kch"][KLOE::channName.at(mctruth_int)]->Fill(Kchrec[5], weight);
      targetHistsReconstructed["mass_Kne"][KLOE::channName.at(mctruth_int)]->Fill(*minv4gam, weight);
      targetHistsReconstructed["mass_pi01"][KLOE::channName.at(mctruth_int)]->Fill(pi01[5], weight);
      targetHistsReconstructed["mass_pi02"][KLOE::channName.at(mctruth_int)]->Fill(pi02[5], weight);
      targetHistsReconstructed["time_neutral_MC"][KLOE::channName.at(mctruth_int)]->Fill(*TrcSum + *T0step1, weight);
      targetHistsReconstructed["combined_mass_pi0"][KLOE::channName.at(mctruth_int)]->Fill(combinedMassPi0Fit, weight);

      targetHistsFitted["mass_Kch"][KLOE::channName.at(mctruth_int)]->Fill(Kchrec[5], weight);
      targetHistsFitted["mass_Kne"][KLOE::channName.at(mctruth_int)]->Fill(*minv4gam, weight);
      targetHistsFitted["bestError"][KLOE::channName.at(mctruth_int)]->Fill(*bestError, weight);
      targetHistsFitted["mass_pi01"][KLOE::channName.at(mctruth_int)]->Fill(pi01Fit[5], weight);
      targetHistsFitted["mass_pi02"][KLOE::channName.at(mctruth_int)]->Fill(pi02Fit[5], weight);
      targetHistsFitted["chi2_signalKinFit"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2SignalKinFit / 10., weight);
      targetHistCounts->Fill(*Chi2SignalKinFit / 10.);
      targetHistsFitted["chi2_trilaterationKinFit"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2TriKinFit / 5., weight);
      targetHistsFitted["chi2_omegaKinFit"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2OmegaKinFit / 8., weight);
      targetHistsFitted["prob_signal"][KLOE::channName.at(mctruth_int)]->Fill(TMath::Prob(*Chi2SignalKinFit, 10), weight);
      targetHistsFitted["combined_mass_pi0"][KLOE::channName.at(mctruth_int)]->Fill(combinedMassPi0Fit, weight);

      for (Int_t i = 0; i < 39; i++)
      {
        targetHistsFitted["pull" + std::to_string(i + 1)][KLOE::channName.at(mctruth_int)]->Fill(pullsSignalFit[i], weight);
      }

      targetHistsFitted["time_neutral_MC"][KLOE::channName.at(mctruth_int)]->Fill(*TrcSum, weight);
      targetHistsFitted["openingAngleCharged"][KLOE::channName.at(mctruth_int)]->Fill(openingAngleCharged, weight);
      targetHistsFitted["openingAngleNeutral"][KLOE::channName.at(mctruth_int)]->Fill(openingAngleNeutral, weight);
      targetHistsFitted["Qmiss"][KLOE::channName.at(mctruth_int)]->Fill(*Qmiss, weight);
      targetHistsReconstructed["delta_t"][KLOE::channName.at(mctruth_int)]->Fill(deltaT, weight);
      targetHistsFitted["delta_t"][KLOE::channName.at(mctruth_int)]->Fill(deltaTfit, weight);
      targetHistsFitted["delta_t_MC"][KLOE::channName.at(mctruth_int)]->Fill(deltaTMC, weight);
      targetHistsFitted["deltaPhiv"][KLOE::channName.at(mctruth_int)]->Fill(deltaPhi, weight);
      targetHistsFitted["deltaPhivFit"][KLOE::channName.at(mctruth_int)]->Fill(deltaPhiFit, weight);
      targetHistsFitted["deltaTheta"][KLOE::channName.at(mctruth_int)]->Fill(deltaTheta, weight);
      targetHistsFitted["TransvRadius"][KLOE::channName.at(mctruth_int)]->Fill(path00MCCenter, weight);
      targetHistsFitted["T0Omega"][KLOE::channName.at(mctruth_int)]->Fill(T0Omega, weight);
      targetHistsFitted["mass_omega"][KLOE::channName.at(mctruth_int)]->Fill(omegaFit[5], weight);
      targetHistsFitted["dist_z_Neu"][KLOE::channName.at(mctruth_int)]->Fill(KnerecFit[8] - ipFit[2], weight);
      targetHistsFitted["dist_z_Ch"][KLOE::channName.at(mctruth_int)]->Fill(KchrecFit[8] - ipFit[2], weight);

      targetHistsFitted["dist_ch_neu_closest_x"][KLOE::channName.at(mctruth_int)]->Fill(distNeutralCharged[0], weight);
      targetHistsFitted["dist_ch_neu_closest_y"][KLOE::channName.at(mctruth_int)]->Fill(distNeutralCharged[1], weight);
      targetHistsFitted["dist_ch_neu_closest_z"][KLOE::channName.at(mctruth_int)]->Fill(distNeutralCharged[2], weight);

      targetHistsFitted["rho_omega"][KLOE::channName.at(mctruth_int)]->Fill(rho, weight);

      if (*muonAlertPlus > 0 || *muonAlertMinus > 0)
        targetHistsFitted["muon_alert"][KLOE::channName.at(mctruth_int)]->Fill(1., weight);
      else
        targetHistsFitted["muon_alert"][KLOE::channName.at(mctruth_int)]->Fill(0., weight);

      targetHists2DFitted["T0_omega_vs_mass_omega"][KLOE::channName.at(mctruth_int)]->Fill(T0Omega, omegaFit[5], weight);
      targetHists2DFitted["delta_t_mc_vs_delta_t_res"][KLOE::channName.at(mctruth_int)]->Fill(deltaTMC, deltaTfit - deltaTMC, weight);
      targetHists2DFitted["delta_t_vs_delta_t_mc"][KLOE::channName.at(mctruth_int)]->Fill(deltaTMC, deltaT - deltaTMC, weight);
      targetHists2DFitted["delta_t_fit_vs_delta_t_mc"][KLOE::channName.at(mctruth_int)]->Fill(deltaTMC, deltaTfit, weight);
      targetHists2DFitted["t_ch_fit_vs_t_ch_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMMC, *KaonChTimeCMSignalFit - *KaonChTimeCMMC, weight);
      targetHists2DFitted["t_neu_fit_vs_t_neu_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMSignalFit - *KaonNeTimeCMMC, weight);
      targetHists2DFitted["t_ch_rec_vs_t_ch_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMMC, propTimes.kaon1TimeCM - *KaonChTimeCMMC, weight);
      targetHists2DFitted["t_neu_rec_vs_t_neu_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, propTimes.kaon2TimeCM - *KaonNeTimeCMMC, weight);
      targetHists2DFitted["t_ch_fit_MC_vs_t_neu_fit_MC"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMSignalFit - *KaonChTimeCMMC, *KaonNeTimeCMSignalFit - *KaonNeTimeCMMC, weight);

      if (*KaonNeTimeCMSignalFit <= *KaonNeTimeCMMC - (KLOE::T0 / 2.))
        targetHists2DFitted["t_ch_fit_vs_t_neu_fit"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMSignalFit, *KaonNeTimeCMSignalFit, weight);

      targetHists2DFitted["chi2_signalKinFit_vs_chi2_trilaterationKinFit"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2SignalKinFit / 10., *Chi2TriKinFit / 5., weight);
      targetHists2DFitted["chi2_signalKinFit_vs_chi2_omegaKinFit"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2SignalKinFit / 10., *Chi2OmegaKinFit / 8., weight);
      targetHists2DFitted["t00_tri_vs_t00_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMBoostTriFit, weight);
      targetHists2DFitted["t00_triangle_vs_t00_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMBoostLor, weight);
      targetHists2DFitted["t00_triangle_vs_t00_tri"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMBoostLor, *KaonNeTimeCMBoostTriFit, weight);
      targetHists2DFitted["Rkne_vs_Rkch"][KLOE::channName.at(mctruth_int)]->Fill(RKne, RKch, weight);
      targetHists2DFitted["chi2_signalKinFit_vs_mass_Kch"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2SignalKinFit, Kchrec[5] - PhysicsConstants::mK0, weight);
      targetHists2DFitted["chi2_signalKinFit_vs_mass_Kne"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2SignalKinFit, Knerec[5] - PhysicsConstants::mK0, weight);
      targetHists2DFitted["chi2_signalKinFit_vs_delta_t"][KLOE::channName.at(mctruth_int)]->Fill(*Chi2SignalKinFit / 10., deltaTfit, weight);
      targetHists2DFitted["TransvRadius_vs_delta_t"][KLOE::channName.at(mctruth_int)]->Fill(pathpm, deltaTfit, weight);
      targetHists2DFitted["Qmiss_vs_deltaPhi"][KLOE::channName.at(mctruth_int)]->Fill(*Qmiss, deltaPhi, weight);
      targetHists2DFitted["mass_pi01_vs_mass_pi02"][KLOE::channName.at(mctruth_int)]->Fill(pi01Fit[5], pi02Fit[5], weight);
      targetHists2DFitted["mass_Kch_vs_mass_Kne"][KLOE::channName.at(mctruth_int)]->Fill(Kchrec[5], *minv4gam, weight);

      Double_t ppmMC = std::sqrt(std::pow(Kchmc[0], 2) + std::pow(Kchmc[1], 2) + std::pow(Kchmc[2], 2));
      Double_t p00MC = std::sqrt(std::pow(Knemc[0], 2) + std::pow(Knemc[1], 2) + std::pow(Knemc[2], 2));
      Double_t ppm = std::sqrt(std::pow(KchrecFit[0], 2) + std::pow(KchrecFit[1], 2) + std::pow(KchrecFit[2], 2));
      Double_t p00 = std::sqrt(std::pow(KnerecFit[0], 2) + std::pow(KnerecFit[1], 2) + std::pow(KnerecFit[2], 2));

      targetHists2DFitted["pt_ch_fit_vs_pt_ch_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMMC, ppm - ppmMC, weight);
      targetHists2DFitted["pt_neu_fit_vs_pt_neu_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, p00 - p00MC, weight);
      targetHists2DFitted["rt_ch_fit_vs_rt_ch_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMMC, radiuspm - radiuspmMC, weight);
      targetHists2DFitted["rt_neu_fit_vs_rt_neu_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, radius00 - radius00MC, weight);
      targetHists2DFitted["z_ch_fit_vs_z_ch_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonChTimeCMMC, zdistpm - zdistpmMC, weight);
      targetHists2DFitted["z_neu_fit_vs_z_neu_mc"][KLOE::channName.at(mctruth_int)]->Fill(*KaonNeTimeCMMC, zdist00 - zdist00MC, weight);
      targetHists2DFitted["Phi_trk1_angle_vs_Phi_trk2_angle"][KLOE::channName.at(mctruth_int)]->Fill(cos(phiTrk1Angle * TMath::Pi() / 180.0), cos(phiTrk2Angle * TMath::Pi() / 180.0), weight);
      targetHists2DFitted["Phi_pi01_angle_vs_Phi_pi02_angle"][KLOE::channName.at(mctruth_int)]->Fill(cos(phipi01Angle * TMath::Pi() / 180.0), cos(phipi02Angle * TMath::Pi() / 180.0), weight);
      targetHists2DFitted["Energy_trk1_vs_Energy_trk2"][KLOE::channName.at(mctruth_int)]->Fill(trk1Fit[3], trk2Fit[3], weight);
      targetHists2DFitted["Energy_pi01_vs_Energy_pi02"][KLOE::channName.at(mctruth_int)]->Fill(pi01Fit[3], pi02Fit[3], weight);

      targetHists2DFitted["dist_ch_neu_closest_x_vs_dist_ch_neu_closest_y"][KLOE::channName.at(mctruth_int)]->Fill(distNeutralCharged[0], distNeutralCharged[1], weight);

      targetHists2DFitted["rho_pm_vs_rho_00"][KLOE::channName.at(mctruth_int)]->Fill(rho_pm, rho_00, weight);

      if (addToPCA)
      {
        dataPCA[0] = T0Omega;
        dataPCA[1] = omegaFit[5];
        AddDataToPCA(dataPCA);
      }
    }
  };

  if (primaryPass && mcflagCondition)
  {
    if ((mctruth_int == 1) && *mcflag == 1)
      signal_num++;

    if (mctruth_int > 1 && *mcflag == 1)
      bkg_tot++;

    fillAcceptedEvent(histsReconstructed,
                      histsFittedSignal,
                      hists2DFittedSignal,
                      channEventsCut,
                      cutPassed,
                      cutNPassed,
                      overflow,
                      histCounts,
                      kTRUE);
  }

  // Fill independent histograms for non-primary scenarios.
  for (const auto &scenario : g_activeScenarios)
  {
    if (scenario == fOption)
      continue;

    auto stateIt = g_secondaryScenarioStates.find(scenario);
    if (stateIt == g_secondaryScenarioStates.end())
      continue;

    ScenarioHistogramState &state = stateIt->second;

    if ((mctruth_int == 1) && *mcflag == 1)
      state.deltaTSignalTot->Fill(deltaTMC, weight);

    if (mctruth_int >= 0 && *mcflag == 1)
      state.channEventsTotal[KLOE::channName.at(mctruth_int)]++;

    const Bool_t passScenario = PassScenario(scenario,
                                             shorterKaonPaths,
                                             oldChi2Cut,
                                             oldTrcSumCut,
                                             oldCombinedMassPi0Cut,
                                             oldMassKchCut,
                                             oldMassKneCut,
                                             oldQmissCut,
                                             oldOpeningAngleCut,
                                             omegaPi0RejectionCut,
                                             fiducialVolume,
                                             simonaChi2Cut,
                                             badClusSimona,
                                             simonaKinCuts,
                                             simonaDeltaPhiCut,
                                             omegaMassT0Cut,
                                             blobCut,
                                             noBlobCut,
                                             newChi2Cut,
                                             newTrkAngleCut,
                                             newCombinedMassPi0Cut,
                                             newOmegaGeometricalCut,
                                             newOmegaT0Cut,
                                             newKchCut);

    if (!passScenario || !mcflagCondition)
      continue;

    fillAcceptedEvent(state.histsReconstructed,
                      state.histsFittedSignal,
                      state.hists2DFittedSignal,
                      state.channEventsCut,
                      state.cutPassed,
                      state.cutNPassed,
                      state.overflow,
                      state.histCounts,
                      kFALSE);
  }

  return kTRUE;
}

void signal_vs_bcg_v2::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void signal_vs_bcg_v2::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file

  gErrorIgnoreLevel = kFatal;

  // MC sum

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *rootOutFile = TFile::Open(g_rootOutputFileName, "RECREATE");
  if (!rootOutFile || rootOutFile->IsZombie())
  {
    std::cerr << "ERROR: Cannot open ROOT output file: " << g_rootOutputFileName << std::endl;
    rootOutFile = nullptr;
  }

  std::map<TString, std::map<TString, Int_t>> integrals;
  std::map<TString, Int_t> maxNum;
  std::map<TString, TString> maxChannel;

  std::map<TString, Float_t> channEventsCutCorr, channEventsTotalCorr;

  for (const auto &config : histogramConfigs1D)
  {
    maxNum[config.first] = 0;
    for (const auto &channelType : KLOE::channName)
    {
      if (histsFittedSignal[config.first][channelType.second]->GetEntries() > 0. && channelType.second != "MC sum")
      {
        histsReconstructed[config.first][channelType.second]->Scale(histsReconstructed[config.first][channelType.second]->GetEntries() / histsReconstructed[config.first][channelType.second]->Integral(0, histsReconstructed[config.first][channelType.second]->GetNbinsX() + 1));

        histsReconstructed[config.first][channelType.second]->Scale(channFactor[channelType.second]);

        histsFittedSignal[config.first][channelType.second]->Scale(histsFittedSignal[config.first][channelType.second]->GetEntries() / histsFittedSignal[config.first][channelType.second]->Integral(0, histsFittedSignal[config.first][channelType.second]->GetNbinsX() + 1));

        histsFittedSignal[config.first][channelType.second]->Scale(channFactor[channelType.second]);

        channEventsCutCorr[channelType.second] = channEventsCut[channelType.second] * channFactor[channelType.second];
        channEventsTotalCorr[channelType.second] = channEventsTotal[channelType.second] * channFactor[channelType.second];
      }

      integrals[config.first][channelType.second] = histsFittedSignal[config.first][channelType.second]->GetMaximum();

      if (integrals[config.first][channelType.second] > maxNum[config.first])
      {
        maxNum[config.first] = integrals[config.first][channelType.second];
        maxChannel[config.first] = channelType.second;
      }
    }
  }

  // Efficiency calculation
  std::map<TString, Float_t> channEffAna;

  Int_t tot_events_after = 0;

  for (const auto &lumi : channEventsTotal)
  {
    if (lumi.first != "Data" && lumi.first != "MC sum")
    {
      tot_events += channEventsCutCorr[lumi.first];
    }

    channEffAna[lumi.first] = CalculateEfficiency(channEventsCutCorr[lumi.first], channEventsTotalCorr[lumi.first]);
  }

  Double_t eff = CalculateEfficiency(signal_num, signal_tot);
  Double_t purity = CalculatePurity(channEventsCutCorr["Signal"], tot_events);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TString resultsFitter = "";
  TPaveText *textresultsFitter = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
  textresultsFitter->SetTextSize(0.03);

  for (const auto &config : histogramConfigs1D)
  {
    for (const auto &channelType : KLOE::channName)
    {
      if (channelType.second != "Data" && channelType.second != "MC sum")
      {
        histsReconstructed[config.first][KLOE::channName.at(8)]->Add(histsReconstructed[config.first][channelType.second]);
        histsFittedSignal[config.first][KLOE::channName.at(8)]->Add(histsFittedSignal[config.first][channelType.second]);
      }
    }
  }

  Double_t meanMass = 0, sigmaMass = 0, meanT0 = 0, sigmaT0 = 0;

  TF1 *chi2 = new TF1("chi2", chi2dist, 0, 1E5, 2);

  for (const auto &config : histogramConfigs1D)
  {
    Bool_t fitDoubleGaus = 0, //(config.first == "mass_Kch" || config.first == "mass_Kne" || config.first == "mass_pi01" || config.first == "mass_pi02" || config.first == "time_neutral_MC"),
        fitOmegaGaus = (config.first == "T0Omega" || config.first == "mass_omega"),
           fitSignalBadClus = (config.first == "deltaPhivFit" && fOption == "BAD_CLUS_SIMONA");
    Bool_t chi2Fit = 0; //(config.first == "chi2_signalKinFit");

    Bool_t logCond = (config.first == "Qmiss" || config.first == "prob_signal" || config.first == "Energy_Kne" || config.first == "bestError"), // || config.first == "chi2_signalKinFit"), // || config.first == "TransvRadius"),
        logCondX = 0;                                                                                                                           //(config.first == "chi2_signalKinFit");

    std::vector<TString> labels;

    // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
    TLegend *legend = new TLegend(0.6, 0.65, 0.9, 0.9, "", "NDC");
    TLegend *deltaTCompLegend = nullptr;

    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.03);

    canvas[config.first]->cd();
    if (logCond)
      canvas[config.first]->SetLogy(1);
    else
      canvas[config.first]->SetLogy(0);

    if (logCondX)
      canvas[config.first]->SetLogx(1);
    else
      canvas[config.first]->SetLogx(0);

    for (const auto &channelType : KLOE::channName)
    {
      if (channelType.second == "MC sum")
        continue;

      histsReconstructed[config.first][channelType.second]->SetLineColor(KLOE::channColor.at(channelType.second) + 1);
      histsFittedSignal[config.first][channelType.second]->SetLineColor(KLOE::channColor.at(channelType.second));

      histsFittedSignal[config.first][channelType.second]->SetTitle("");

      if (logCond)
        histsFittedSignal[config.first][channelType.second]->GetYaxis()->SetRangeUser(10, maxNum[config.first] * 5);
      else
        histsFittedSignal[config.first][channelType.second]->GetYaxis()->SetRangeUser(0, maxNum[config.first] * 1.2);

      if (channelType.second == "Data")
      {
        histsFittedSignal[config.first][channelType.second]->Draw("PE");
        legend->AddEntry(histsFittedSignal[config.first][channelType.second], channelType.second, "pe");
      }
      else
      {
        if ((channelType.second == "Signal" && fitDoubleGaus) || (channelType.second == "Omega" && fitOmegaGaus))
        {
          textresultsFitter->Clear();

          if (config.first == "T0Omega")
          {
            doubleFitter->FitHistogram(histsFittedSignal[config.first][channelType.second], KLOE::DoublGaussFitter::FitType::kDoubleGauss);
            meanT0 = doubleFitter->GetLastResults().mean;
            sigmaT0 = doubleFitter->GetLastResults().coreSigma;
          }
          else if (config.first == "mass_omega")
          {
            doubleFitter->FitHistogram(histsFittedSignal[config.first][channelType.second], KLOE::DoublGaussFitter::FitType::kSingleGauss);
            meanMass = doubleFitter->GetLastResults().mean;
            sigmaMass = doubleFitter->GetLastResults().coreSigma;
          }
          else
          {
            doubleFitter->FitHistogram(histsFittedSignal[config.first][channelType.second], KLOE::DoublGaussFitter::FitType::kDoubleGauss);
          }

          doubleFitter->DrawFitOnCurrentPad(true, false);

          textresultsFitter->AddText(Form("Fit results:"));
          textresultsFitter->AddText(Form("Chi2 / ndof: %.3f / %d", doubleFitter->GetLastResults().chi2, doubleFitter->GetLastResults().ndf));
          textresultsFitter->AddText(Form("Mean: %.3f #pm %.3f", doubleFitter->GetLastResults().mean, doubleFitter->GetLastResults().meanErr));
          textresultsFitter->AddText(Form("Sigma: %.3f #pm %.3f", doubleFitter->GetLastResults().coreSigma, doubleFitter->GetLastResults().coreSigmaErr));
          textresultsFitter->Draw();
        }

        if ((channelType.second == "Signal" && (fitSignalBadClus || chi2Fit)))
        {
          textresultsFitter->Clear();

          Double_t mean = histsFittedSignal[config.first][channelType.second]->GetMean();
          Double_t rms = histsFittedSignal[config.first][channelType.second]->GetRMS();
          Double_t maxBin = histsFittedSignal[config.first][channelType.second]->GetMaximum();
          Double_t integral = histsFittedSignal[config.first][channelType.second]->Integral();
          Double_t binWidth = histsFittedSignal[config.first][channelType.second]->GetBinWidth(1);

          Int_t dof = 10;

          // Użyj "gaus" i oblicz poprawną amplitudę
          TF1 *gausFit = new TF1("gausFit", "gausn", mean - 10 * rms, mean + 10 * rms);

          TF1 *chi2FitDist = new TF1("chi2FitDist", chi2dist, 0, 1E5, 2);

          // Parametry dla "gaus": [0] = amplituda, [1] = mean, [2] = sigma
          // Amplituda = Integral / (sigma * std::sqrt(2*pi))
          Double_t amplitude = maxBin; // * binWidth / (rms * TMath::Sqrt(2 * TMath::Pi()));

          gausFit->SetParameter(0, amplitude);
          gausFit->SetParameter(1, mean);
          gausFit->SetParameter(2, 0.005 * rms);

          // Opcjonalnie: ustaw granice parametrów
          gausFit->SetParLimits(0, 0, 100 * amplitude);             // Amplituda > 0
          gausFit->SetParLimits(1, mean - 1 * rms, mean + 1 * rms); // Amplituda > 0
          gausFit->SetParLimits(2, 0.001 * rms, 0.01 * rms);        // Sigma rozumna

          chi2FitDist->SetParameter(0, integral); // degrees of freedom
          chi2FitDist->SetParameter(1, dof);

          TFitResultPtr result;

          if (config.first == "T0Omega")
          {
            result = histsFittedSignal[config.first][channelType.second]->Fit(gausFit, "RSQ", "", mean - rms, mean + rms);

            meanT0 = result->Parameter(1);
            sigmaT0 = result->Parameter(2);
          }
          else if (config.first == "mass_omega")
          {
            result = histsFittedSignal[config.first][channelType.second]->Fit(gausFit, "RSQ", "", mean - rms, mean + rms);

            meanMass = result->Parameter(1);
            sigmaMass = result->Parameter(2);
          }
          else if (chi2Fit)
          {
            result = histsFittedSignal[config.first][channelType.second]->Fit(chi2FitDist, "RS");
          }
          else
          {
            result = histsFittedSignal[config.first][channelType.second]->Fit(gausFit, "RS");
          }

          histsFittedSignal[config.first][channelType.second]->Sumw2(false);

          textresultsFitter->AddText(Form("Fit results:"));
          textresultsFitter->AddText(Form("Chi2/NDF: %.3f / %d", result->Chi2(), result->Ndf()));
          textresultsFitter->AddText(Form("Mean: %.3f #pm %.3f", result->Parameter(1), result->ParError(1)));
          textresultsFitter->AddText(Form("Sigma: %.3f #pm %.3f", result->Parameter(2), result->ParError(2)));
          textresultsFitter->Draw();
        }

        if (config.first == "delta_t" && channelType.second == "Signal")
        {
          histsFittedSignal["delta_t"][channelType.second]->SetLineColor(kBlue + 1);
          histsFittedSignal["delta_t"][channelType.second]->SetLineWidth(3);
          histsFittedSignal["delta_t"][channelType.second]->Draw("HIST SAME");
          histsReconstructed["delta_t"][channelType.second]->SetLineColor(kRed + 1);
          histsReconstructed["delta_t"][channelType.second]->SetLineWidth(3);
          histsReconstructed["delta_t"][channelType.second]->Draw("HIST SAME");
          histsFittedSignal["delta_t_MC"][channelType.second]->SetLineColor(kGreen + 2);
          histsFittedSignal["delta_t_MC"][channelType.second]->SetLineStyle(2);
          histsFittedSignal["delta_t_MC"][channelType.second]->SetLineWidth(3);
          histsFittedSignal["delta_t_MC"][channelType.second]->Draw("HIST SAME");

          deltaTCompLegend = new TLegend(0.18, 0.72, 0.52, 0.9, "", "NDC");
          deltaTCompLegend->SetBorderSize(1);
          deltaTCompLegend->SetFillColor(kWhite);
          deltaTCompLegend->SetTextSize(0.03);
          deltaTCompLegend->AddEntry(histsFittedSignal["delta_t"][channelType.second], "Fitted #Delta t", "l");
          deltaTCompLegend->AddEntry(histsReconstructed["delta_t"][channelType.second], "Reconstructed #Delta t", "l");
          deltaTCompLegend->AddEntry(histsFittedSignal["delta_t_MC"][channelType.second], "MC generated #Delta t", "l");
        }
        else
        {
          for (Int_t bin = 1; bin <= histsFittedSignal[config.first][channelType.second]->GetNbinsX(); bin++)
          {
            histsFittedSignal[config.first][channelType.second]->SetBinError(bin, 0);
          }
          histsFittedSignal[config.first][channelType.second]->Draw("HIST SAME");
        }

        gPad->Update();
        legend->AddEntry(histsFittedSignal[config.first][channelType.second], channelType.second, "l");
      }
    }

    // meanT0 = 132.215;
    // sigmaT0 = 31.848;

    // meanMass = 776.571;
    // sigmaMass = 28.026;

    // if (config.first == "chi2_signalKinFit")
    // {
    //   chi2->SetParameter(0, histsFittedSignal[config.first]["Signal"]->Integral(1,30));
    //   chi2->SetParameter(1, 10); // degrees of freedom
    //   chi2->SetLineColor(KLOE::channColor.at("Signal"));
    //   histsFittedSignal[config.first]["Data"]->GetYaxis()->SetRangeUser(0, chi2->GetMaximum(10.0, 10.0) * 1.2);
    //   chi2->Draw("SAME");
    // }

    labels.push_back(Form("Events: %d", tot_events));

    // labels.push_back(Form("Chi2 (Data vs. MC Sum): %.3f", histsFittedSignal[config.first]["MC sum"]->Chi2Test(histsFittedSignal[config.first]["Data"], "WU CHI2/NDF")));

    DrawLabelOnHisto(labels);

    gPad->Update();

    if (!fitOmegaGaus && !fitSignalBadClus && !(config.first == "TransvRadius") && channEventsTotal["Omega"] > 0)
      legend->Draw();

    if (deltaTCompLegend)
      deltaTCompLegend->Draw();

    // Sprawdź czy JAKIKOLWIEK histogram ma wpisy
    Bool_t hasEntries = kFALSE;

    // Sprawdź wszystkie kanały
    for (const auto &name : KLOE::channName)
    {
      if (histsFittedSignal[config.first][name.second]->GetEntries() > 0)
      {
        hasEntries = kTRUE;
        break;
      }
    }

    // Jeśli żaden histogram nie ma wpisów - pomiń ten canvas
    if (!hasEntries)
    {
      continue;
    }

    canvas[config.first]->SaveAs(Form("%s/%s_comparison%s", folderPath.Data(), config.first.Data(), Paths::ext_img.Data()));
  }

  Double_t sigmas[2], means[2] = {154.061, 779.152}, vLong[2], vTransv[2];

  meanT0 = means[0];
  meanMass = means[1];

  sigmaT0 = 5.691;
  sigmaMass = 5.620;

  // Calculation of PCA components
  // WidthOfCorrelatedHist(means, sigmas, vLong, vTransv);

  // ===== SLOPE z PCA (z skalowanych wektorów) =====
  Double_t a = 1;                       // GetPCASlope(); // Pobierz slope bezpośrednio ze skalowanych wektorów
  Double_t b = means[1] - a * means[0]; // y - a*x w centrum masowym

  std::cout << "\n=== SLOPE FROM PCA (via GetPCASlope) ===" << std::endl;
  std::cout << "Slope a = " << a << std::endl;
  std::cout << "Intercept b = " << b << std::endl;

  // ===== HARDCUT Y limits =====
  // Dolna linia na y = 620, górna na y = 640
  Double_t yLowerLimit = 620.0;
  Double_t yUpperLimit = 640.0;

  // Średnia wartość między limitami to efektywny "center"
  Double_t yCenterLimit = (yLowerLimit + yUpperLimit) / 2.0;
  Double_t yHalfWidth = (yUpperLimit - yLowerLimit) / 2.0;

  // Przelicz to na Breal (odległość prostopadła od linii)
  Double_t normFactor = TMath::Sqrt(1 + a * a);
  Double_t Breal = yHalfWidth * normFactor; // Przekształć różnicę Y na odległość prostopadłą

  std::cout << "\n=== HARDCUT Y LIMITS ===" << std::endl;
  std::cout << "Y lower limit: " << yLowerLimit << std::endl;
  std::cout << "Y upper limit: " << yUpperLimit << std::endl;
  std::cout << "Y center: " << yCenterLimit << std::endl;
  std::cout << "Y half-width: " << yHalfWidth << std::endl;
  std::cout << "Breal (perpendicular distance): " << Breal << std::endl;
  std::cout << "=======================\n"
            << std::endl;

  std::cout << "Trend line: y = " << a << " * x + " << b << std::endl;

  for (const auto &config : histogramConfigs2D)
  {

    Bool_t withProfile = (config.first == "t_ch_fit_vs_t_ch_mc" || config.first == "t_neu_fit_vs_t_neu_mc" || config.first == "t_ch_rec_vs_t_ch_mc" || config.first == "t_neu_rec_vs_t_neu_mc" || config.first == "delta_t_mc_vs_delta_t_res" || config.first == "delta_t_vs_delta_t_mc");

    TH2 *h2D = hists2DFittedSignal[config.first]["Signal"];

    TCanvas *c;

    if (!withProfile)
      goto skip_profile;

    // Stwórz canvas z oboma profilami
    c = CreateCanvasWithProfiles(h2D,
                                 config.first + "_with_profiles",
                                 kTRUE,  // Rysuj mean profile
                                 kTRUE); // Rysuj sigma profile

    c->SaveAs(folderPath + "/" + config.first + "_with_profiles" + Paths::ext_img.Data());
    canvasProfiles[config.first] = c;

    continue;

  skip_profile: // Go to here if no profiles needed

    for (const auto &channelType : KLOE::channName)
    {
      canvas2D[config.first][channelType.second]->cd();
      canvas2D[config.first][channelType.second]->SetLogz(1);

      hists2DFittedSignal[config.first][channelType.second]->Draw("COLZ");

      if (config.first == "T0_omega_vs_mass_omega") // && channelType.second == "Signal")
      {
        // Sprawdź zakres histogramu
        Double_t xMin = hists2DFittedSignal[config.first][channelType.second]->GetXaxis()->GetXmin();
        Double_t xMax = hists2DFittedSignal[config.first][channelType.second]->GetXaxis()->GetXmax();
        Double_t yMin = hists2DFittedSignal[config.first][channelType.second]->GetYaxis()->GetXmin();
        Double_t yMax = hists2DFittedSignal[config.first][channelType.second]->GetYaxis()->GetXmax();

        TLine *pc1_line_cut = new TLine(xMin,
                                        a * xMin + b + Breal,
                                        xMax,
                                        a * xMax + b + Breal);

        TLine *pc2_line_cut = new TLine(xMin,
                                        a * xMin + b - Breal,
                                        xMax,
                                        a * xMax + b - Breal);

        TLine *T0_left_cut = new TLine(meanT0 - 3 * sigmaT0, meanMass - 3 * sigmaMass,
                                       meanT0 - 3 * sigmaT0, meanMass + 3 * sigmaMass);
        TLine *T0_right_cut = new TLine(meanT0 + 3 * sigmaT0, meanMass - 3 * sigmaMass,
                                        meanT0 + 3 * sigmaT0, meanMass + 3 * sigmaMass);

        TLine *Mass_left_cut = new TLine(meanT0 - 3 * sigmaT0, meanMass - 3 * sigmaMass,
                                         meanT0 + 3 * sigmaT0, meanMass - 3 * sigmaMass);
        TLine *Mass_right_cut = new TLine(meanT0 - 3 * sigmaT0, meanMass + 3 * sigmaMass,
                                          meanT0 + 3 * sigmaT0, meanMass + 3 * sigmaMass);

        std::cout << meanT0 << " " << sigmaT0 << " " << meanMass << " " << sigmaMass << std::endl;

        pc1_line_cut->SetLineColor(kBlack);
        pc1_line_cut->SetLineWidth(3);
        pc1_line_cut->Draw("SAME"); // DODAJ "SAME"!

        pc2_line_cut->SetLineColor(kBlack);
        pc2_line_cut->SetLineWidth(3);
        pc2_line_cut->Draw("SAME"); // DODAJ "SAME"!

        T0_left_cut->SetLineColor(kBlack);
        T0_left_cut->SetLineWidth(3);
        T0_left_cut->Draw("SAME");

        T0_right_cut->SetLineColor(kBlack);
        T0_right_cut->SetLineWidth(3);
        T0_right_cut->Draw("SAME");

        Mass_left_cut->SetLineColor(kBlack);
        Mass_left_cut->SetLineWidth(3);
        Mass_left_cut->Draw("SAME");

        Mass_right_cut->SetLineColor(kBlack);
        Mass_right_cut->SetLineWidth(3);
        Mass_right_cut->Draw("SAME");

        canvas2D[config.first][channelType.second]->GetListOfPrimitives()->Add(pc1_line_cut);
        canvas2D[config.first][channelType.second]->GetListOfPrimitives()->Add(pc2_line_cut);
        canvas2D[config.first][channelType.second]->Modified();
        canvas2D[config.first][channelType.second]->Update();
      }

      // Sprawdź czy JAKIKOLWIEK histogram ma wpisy
      Bool_t hasEntries = kFALSE;

      // Sprawdź wszystkie kanały

      if (hists2DFittedSignal[config.first][channelType.second]->GetEntries() > 0)
      {
        hasEntries = kTRUE;
      }

      // Jeśli żaden histogram nie ma wpisów - pomiń ten canvas
      if (!hasEntries)
      {
        continue;
      }

      canvas2D[config.first][channelType.second]->SaveAs(Form("%s/%s_%s_2D%s", folderPath.Data(), config.first.Data(), channelType.second.Data(), Paths::ext_img.Data()));
    }
  }

  Double_t scale = deltaTSignalTot->GetEntries() / deltaTSignalTot->Integral(0, deltaTSignalTot->GetNbinsX() + 1);

  deltaTSignalTot->Scale(scale);

  efficiency = new TEfficiency(*histsFittedSignal["delta_t_MC"]["Signal"], *deltaTSignalTot);

  efficiency->SetUseWeightedEvents(kTRUE);

  efficiency->SetStatisticOption(TEfficiency::kBUniform);

  canvaEff->cd();
  DrawEfficiencyWithFixedYAxis(efficiency, deltaTSignalTot, 0.0, 1.0, "Efficiency");

  // Dodaj box z metrykami
  // TPaveText *metricsBox = new TPaveText(0.15, 0.15, 0.45, 0.35, "NDC");
  TPaveText *metricsBox = new TPaveText(0.15, 0.85, 0.45, 0.9, "NDC");
  metricsBox->SetFillColor(kWhite);
  metricsBox->SetBorderSize(1);
  metricsBox->SetTextAlign(12);
  metricsBox->SetTextSize(0.03);

  metricsBox->AddText(Form("Efficiency: %.2f%%", channEffAna["Signal"] * 100));

  metricsBox->Draw();

  canvaEff->SaveAs(folderPath + "/efficiency_delta_t" + Paths::ext_img.Data());
  if (rootOutFile)
  {
    TDirectory *primaryDir = rootOutFile->GetDirectory("primary");
    if (!primaryDir)
      primaryDir = rootOutFile->mkdir("primary");
    if (primaryDir)
    {
      primaryDir->cd();
      efficiency->Write("efficiency_delta_t_obj", TObject::kOverwrite);
      canvaEff->Write("efficiency_delta_t_canvas", TObject::kOverwrite);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////

  scale = deltaTTot->GetEntries() / deltaTTot->Integral(0, deltaTTot->GetNbinsX() + 1);

  deltaTTot->Scale(scale);

  efficiency = new TEfficiency(*histsFittedSignal["delta_t"]["Signal"], *histsFittedSignal["delta_t"]["MC sum"]);
  efficiency->SetUseWeightedEvents(kTRUE);

  efficiency->SetStatisticOption(TEfficiency::kBUniform);

  canvaPurity->cd();
  DrawEfficiencyWithFixedYAxis(efficiency, histsFittedSignal["delta_t"]["MC sum"], 0.0, 1.0, "Purity");

  // Dodaj box z metrykami
  // TPaveText *metricsBox = new TPaveText(0.15, 0.15, 0.45, 0.35, "NDC");
  metricsBox = new TPaveText(0.15, 0.85, 0.45, 0.9, "NDC");
  metricsBox->SetFillColor(kWhite);
  metricsBox->SetBorderSize(1);
  metricsBox->SetTextAlign(12);
  metricsBox->SetTextSize(0.03);

  metricsBox->AddText(Form("Purity: %.2f%%", purity * 100));

  metricsBox->Draw();

  canvaPurity->SaveAs(folderPath + "/purity_delta_t" + Paths::ext_img.Data());
  if (rootOutFile)
  {
    TDirectory *primaryDir = rootOutFile->GetDirectory("primary");
    if (!primaryDir)
      primaryDir = rootOutFile->mkdir("primary");
    if (primaryDir)
    {
      primaryDir->cd();
      efficiency->Write("purity_delta_t_obj", TObject::kOverwrite);
      canvaPurity->Write("purity_delta_t_canvas", TObject::kOverwrite);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////

  // === Analiza czystości w podzbiorach |ΔT| ===
  std::vector<PuritySubset> purity_subsets = CalculatePurityInSubsets(
      histsFittedSignal["delta_t"]["Signal"], // Histogram sygnału
      histsFittedSignal["delta_t"]["MC sum"], // Histogram całkowitego
      30.0,                                   // max limit [ns]
      30                                      // liczba podzbiorów
  );

  TCanvas *canvas_purity_subsets = DrawPuritySubsets(purity_subsets, "delta_t");
  if (canvas_purity_subsets)
  {
    canvas_purity_subsets->SaveAs(folderPath + "/purity_subsets_delta_t" + Paths::ext_img.Data());
    if (rootOutFile)
    {
      TDirectory *primaryDir = rootOutFile->GetDirectory("primary");
      if (!primaryDir)
        primaryDir = rootOutFile->mkdir("primary");
      if (primaryDir)
      {
        primaryDir->cd();
        canvas_purity_subsets->Write("purity_subsets_delta_t_canvas", TObject::kOverwrite);
      }
    }
  }

  // -- REPORT OF THE CUT STATISTICS IN THE ANALYSIS --

  std::cout << "Total Efficiency signal: " << 100 * eff << " % (" << signal_num << "/" << signal_tot << ")" << std::endl;

  for (const auto &entry : channEffAna)
  {
    std::cout << "Channel: " << entry.first << ", Efficiency: " << entry.second * 100 << " % (" << channEventsCutCorr[entry.first] << "/" << channEventsTotalCorr[entry.first] << "), (" << channEventsCut[entry.first] << "/" << channEventsTotal[entry.first] << ")" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Background events: " << bkg_tot << std::endl;

  std::cout << "Purity: " << 100 * purity << " % (" << channEventsCutCorr["Signal"] << "/" << tot_events << ")" << std::endl;

  std::cout << "\n=== Single-pass multi-scenario summary ===" << std::endl;
  for (const auto &entry : g_scenarioCounters)
  {
    const TString &scenario = entry.first;
    const ScenarioCounters &cnt = entry.second;

    const Double_t effScenario = (cnt.signal_tot > 0)
                                     ? static_cast<Double_t>(cnt.sel_mctruth_by_channel.at("Signal")) / cnt.sel_mctruth_by_channel_before_cut.at("Signal")
                                     : 0.0; // Not taking into account 0 mctruth
    const Int_t totalScenario = cnt.sel_mctruth_by_channel.at("Signal") + cnt.bkg_tot;
    const Double_t purityScenario = (totalScenario > 0)
                                        ? static_cast<Double_t>(cnt.sel_mctruth_by_channel.at("Signal")) / totalScenario
                                        : 0.0;

    const Int_t denomPreselection = cnt.signal_tot;
    const Int_t numerPreselection = cnt.sel_mctruth_0_before_cut + cnt.sel_mctruth_1_before_cut;
    const Double_t effPreselection = ComputeSafeRatio(numerPreselection, denomPreselection);
    const Double_t effSelection = ComputeSafeRatio(cnt.sel_mctruth_1 + cnt.sel_mctruth_0, cnt.sel_mctruth_0_before_cut + cnt.sel_mctruth_1_before_cut); // Taking into account 0 mctruth (thus not taking into account initial cuts!)
    const Double_t effTotalTruth = ComputeSafeRatio(cnt.sel_mctruth_1 + cnt.sel_mctruth_0, denomPreselection);

    std::cout << "Scenario: " << scenario
              << ", Efficiency: " << 100.0 * effScenario << " % (" << cnt.sel_mctruth_by_channel.at("Signal") << "/" << cnt.sel_mctruth_by_channel_before_cut.at("Signal") << ")"
              << ", Purity: " << 100.0 * purityScenario << " % (" << cnt.sel_mctruth_by_channel.at("Signal") << "/" << totalScenario << ")"
              << ", Preselection eff: " << 100.0 * effPreselection << " % ((mctruth 0+1)/(mctruth -1+0+1) = " << numerPreselection << "/" << denomPreselection << ")"
              << ", Selection eff: " << 100.0 * effSelection << " % ((mctruth 1)/(mctruth 0+1) = " << cnt.sel_mctruth_1 + cnt.sel_mctruth_0 << "/" << (numerPreselection) << ")"
              << ", Total eff: " << 100.0 * effTotalTruth << " % ((mctruth 1)/(mctruth -1+0+1) = " << cnt.sel_mctruth_1 + cnt.sel_mctruth_0 << "/" << denomPreselection << ")"
              << ", All events without errors before cuts (bkg + signal): " << cnt.passed_events
              << std::endl;
  }

  // Save scenario summary reports to log/YYYY-MM-DD with timestamp.
  {
    std::time_t now = std::time(nullptr);
    std::tm *lt = std::localtime(&now);
    char dateBuf[11];  // YYYY-MM-DD
    char stampBuf[20]; // YYYY-MM-DD_HH-MM-SS
    std::strftime(dateBuf, sizeof(dateBuf), "%Y-%m-%d", lt);
    std::strftime(stampBuf, sizeof(stampBuf), "%Y-%m-%d_%H-%M-%S", lt);

    TString logDir = Form("log/%s", dateBuf);
    system(Form("mkdir -p %s", logDir.Data()));

    TString csvPath = Form("%s/scenario_report_%s.csv", logDir.Data(), stampBuf);
    TString txtPath = Form("%s/scenario_report_%s.txt", logDir.Data(), stampBuf);

    std::ofstream csv(csvPath.Data());
    std::ofstream txt(txtPath.Data());

    if (csv.is_open())
    {
      csv << "scenario,signal_num_before_cuts,signal_m1_0_1,bkg_tot,total_selected,purity,all_events_before_cut,mctruth_m1_total,mctruth_0_total,mctruth_1_total,mctruth_0_selected,mctruth_1_selected,preselection_eff,selection_eff,efficiency,total_eff\n";
      for (const auto &entry : g_scenarioCounters)
      {
        const TString &scenario = entry.first;
        const ScenarioCounters &cnt = entry.second;
        const Double_t effScenario = (cnt.signal_tot > 0)
                                         ? static_cast<Double_t>(cnt.sel_mctruth_by_channel.at("Signal")) / cnt.sel_mctruth_by_channel_before_cut.at("Signal")
                                         : 0.0; // Not taking into account 0 mctruth
        const Int_t totalScenario = cnt.sel_mctruth_by_channel.at("Signal") + cnt.bkg_tot;
        const Double_t purityScenario = (totalScenario > 0)
                                            ? static_cast<Double_t>(cnt.sel_mctruth_by_channel.at("Signal")) / totalScenario
                                            : 0.0;

        const Int_t denomPreselection = cnt.signal_tot;
        const Int_t numerPreselection = cnt.sel_mctruth_0_before_cut + cnt.sel_mctruth_1_before_cut;
        const Double_t effPreselection = ComputeSafeRatio(numerPreselection, denomPreselection);
        const Double_t effSelection = ComputeSafeRatio(cnt.sel_mctruth_1 + cnt.sel_mctruth_0, cnt.sel_mctruth_0_before_cut + cnt.sel_mctruth_1_before_cut); // Taking into account 0 mctruth (thus not taking into account initial cuts!)
        const Double_t effTotalTruth = ComputeSafeRatio(cnt.sel_mctruth_1 + cnt.sel_mctruth_0, denomPreselection);

        csv
            << scenario.Data() << ","
            << cnt.sel_mctruth_1_before_cut + cnt.sel_mctruth_0_before_cut << ","
            << cnt.signal_tot << ","
            << cnt.bkg_tot << ","
            << totalScenario << ","
            << purityScenario << ","
            << cnt.passed_events << ","
            << cnt.sel_mctruth_m1 << ","
            << cnt.sel_mctruth_0_before_cut << ","
            << cnt.sel_mctruth_1_before_cut << ","
            << cnt.sel_mctruth_0 << ","
            << cnt.sel_mctruth_1 << ","
            << effPreselection << ","
            << effSelection << ","
            << effScenario << ","
            << effTotalTruth << "\n";
      }
    }

    if (txt.is_open())
    {
      txt << "Scenario report generated at: " << stampBuf << "\n\n";
      for (const auto &entry : g_scenarioCounters)
      {
        const TString &scenario = entry.first;
        const ScenarioCounters &cnt = entry.second;
        const Double_t effScenario = (cnt.signal_tot > 0)
                                         ? static_cast<Double_t>(cnt.sel_mctruth_by_channel.at("Signal")) / cnt.sel_mctruth_by_channel_before_cut.at("Signal")
                                         : 0.0; // Not taking into account 0 mctruth
        const Int_t totalScenario = cnt.sel_mctruth_by_channel.at("Signal") + cnt.bkg_tot;
        const Double_t purityScenario = (totalScenario > 0)
                                            ? static_cast<Double_t>(cnt.sel_mctruth_by_channel.at("Signal")) / totalScenario
                                            : 0.0;

        const Int_t denomPreselection = cnt.signal_tot;
        const Int_t numerPreselection = cnt.sel_mctruth_0_before_cut + cnt.sel_mctruth_1_before_cut;
        const Double_t effPreselection = ComputeSafeRatio(numerPreselection, denomPreselection);
        const Double_t effSelection = ComputeSafeRatio(cnt.sel_mctruth_1 + cnt.sel_mctruth_0, cnt.sel_mctruth_0_before_cut + cnt.sel_mctruth_1_before_cut); // Taking into account 0 mctruth (thus not taking into account initial cuts!)
        const Double_t effTotalTruth = ComputeSafeRatio(cnt.sel_mctruth_1 + cnt.sel_mctruth_0, denomPreselection);

        txt << "Scenario: " << scenario << "\n"
            << "  signal_num_before_cuts: " << cnt.sel_mctruth_1_before_cut + cnt.sel_mctruth_0_before_cut << "\n"
            << "  signal_m1_0_1: " << cnt.signal_tot << "\n"
            << "  bkg_tot: " << cnt.bkg_tot << "\n"
            << "  total_selected: " << totalScenario << "\n"
            << "  purity: " << purityScenario << "\n"
            << "  mctruth_-1_selected: " << cnt.sel_mctruth_m1 << "\n"
            << "  mctruth_0_selected: " << cnt.sel_mctruth_0 << "\n"
            << "  mctruth_1_selected: " << cnt.sel_mctruth_1 << "\n"
            << "  preselection_eff ((mctruth 0+1)/(mctruth -1+0+1)): " << effPreselection << "\n"
            << "  selection_eff ((mctruth 1)/(mctruth 0+1)): " << effSelection << "\n"
            << "  efficiency: " << effScenario << "\n"
            << "  total_eff ((mctruth 1)/(mctruth -1+0+1)): " << effTotalTruth << "\n"
            << "  all_events_before_cut: " << cnt.passed_events << "\n\n";
      }
    }

    std::cout << "Saved scenario reports: " << csvPath << " and " << txtPath << std::endl;
  }

  // Render separate output plots for all non-primary scenarios.
  for (auto &entry : g_secondaryScenarioStates)
  {
    const TString &scenario = entry.first;
    ScenarioHistogramState &state = entry.second;

    std::cout << "Rendering plots for scenario: " << scenario << std::endl;

    // Build MC sum for secondary scenarios, analogicznie do scenariusza glownego.
    for (const auto &config : histogramConfigs1D)
    {
      TH1 *mcSumRec = state.histsReconstructed[config.first][KLOE::channName.at(8)];
      TH1 *mcSumFit = state.histsFittedSignal[config.first][KLOE::channName.at(8)];

      if (mcSumRec)
        mcSumRec->Reset();
      if (mcSumFit)
        mcSumFit->Reset();

      for (const auto &channelType : KLOE::channName)
      {
        if (channelType.second == "Data" || channelType.second == "MC sum")
          continue;

        TH1 *hRec = state.histsReconstructed[config.first][channelType.second];
        TH1 *hFit = state.histsFittedSignal[config.first][channelType.second];

        if (mcSumRec && hRec)
          mcSumRec->Add(hRec);
        if (mcSumFit && hFit)
          mcSumFit->Add(hFit);
      }
    }

    for (const auto &config : histogramConfigs1D)
    {
      TCanvas *c = new TCanvas(Form("c_%s_%s", config.first.Data(), scenario.Data()),
                               Form("Canvas for %s (%s)", config.first.Data(), scenario.Data()),
                               750, 750);

      Double_t maxVal = 0.0;
      for (const auto &channelType : KLOE::channName)
      {
        TH1 *h = state.histsFittedSignal[config.first][channelType.second];
        if (h && h->GetMaximum() > maxVal)
          maxVal = h->GetMaximum();
      }

      Bool_t drawn = kFALSE;
      TLegend *legend = new TLegend(0.6, 0.65, 0.9, 0.9, "", "NDC");
      legend->SetBorderSize(1);
      legend->SetFillColor(kWhite);
      legend->SetTextSize(0.03);

      for (const auto &channelType : KLOE::channName)
      {
        if (channelType.second == "MC sum")
          continue;

        TH1 *h = state.histsFittedSignal[config.first][channelType.second];
        if (!h || h->GetEntries() <= 0)
          continue;

        h->SetLineColor(KLOE::channColor.at(channelType.second));
        h->SetTitle("");
        h->GetYaxis()->SetRangeUser(0, maxVal > 0.0 ? maxVal * 1.2 : 1.0);

        if (!drawn)
        {
          if (channelType.second == "Data")
            h->Draw("PE");
          else
            h->Draw("HIST");
          drawn = kTRUE;
        }
        else
        {
          if (channelType.second == "Data")
            h->Draw("PE SAME");
          else
            h->Draw("HIST SAME");
        }

        legend->AddEntry(h, channelType.second, channelType.second == "Data" ? "pe" : "l");
      }

      if (drawn)
      {
        legend->Draw();
        c->SaveAs(Form("%s/%s_comparison%s", state.folderPath.Data(), config.first.Data(), Paths::ext_img.Data()));
      }

      delete c;
    }

    for (const auto &config : histogramConfigs2D)
    {
      for (const auto &channelType : KLOE::channName)
      {
        TH2 *h2 = state.hists2DFittedSignal[config.first][channelType.second];
        if (!h2 || h2->GetEntries() <= 0)
          continue;

        TCanvas *c2 = new TCanvas(Form("c2_%s_%s_%s", config.first.Data(), channelType.second.Data(), scenario.Data()),
                                  Form("Canvas2D for %s %s (%s)", config.first.Data(), channelType.second.Data(), scenario.Data()),
                                  750, 750);
        c2->SetLogz(1);
        h2->Draw("COLZ");
        c2->SaveAs(Form("%s/%s_%s_2D%s", state.folderPath.Data(), config.first.Data(), channelType.second.Data(), Paths::ext_img.Data()));
        delete c2;
      }
    }

    if (state.deltaTSignalTot && state.deltaTSignalTot->Integral(0, state.deltaTSignalTot->GetNbinsX() + 1) > 0.0)
    {
      TH1 *numEff = state.histsFittedSignal["delta_t_MC"]["Signal"];
      if (numEff && numEff->Integral(0, numEff->GetNbinsX() + 1) > 0.0)
      {
        TCanvas *cEff = new TCanvas(Form("Efficiency_%s", scenario.Data()), "Efficiency", 800, 800);
        TEfficiency *effSc = new TEfficiency(*numEff, *state.deltaTSignalTot);
        effSc->SetUseWeightedEvents(kTRUE);
        effSc->SetStatisticOption(TEfficiency::kBUniform);
        cEff->cd();
        DrawEfficiencyWithFixedYAxis(effSc, state.deltaTSignalTot, 0.0, 1.0, "Efficiency");
        cEff->SaveAs(state.folderPath + "/efficiency_delta_t" + Paths::ext_img.Data());
        if (rootOutFile)
        {
          TString scenDirName = TString("scenario_") + scenario;
          TDirectory *scDir = rootOutFile->GetDirectory(scenDirName);
          if (!scDir)
            scDir = rootOutFile->mkdir(scenDirName);
          if (scDir)
          {
            scDir->cd();
            effSc->Write("efficiency_delta_t_obj", TObject::kOverwrite);
            cEff->Write("efficiency_delta_t_canvas", TObject::kOverwrite);
          }
        }
        delete cEff;
      }
    }

    TH1 *numPurity = state.histsFittedSignal["delta_t"]["Signal"];
    TH1 *denPurity = state.histsFittedSignal["delta_t"]["MC sum"];
    if (numPurity && denPurity && denPurity->Integral(0, denPurity->GetNbinsX() + 1) > 0.0)
    {
      TCanvas *cPurity = new TCanvas(Form("Purity_%s", scenario.Data()), "Purity", 800, 800);
      TEfficiency *purSc = new TEfficiency(*numPurity, *denPurity);
      purSc->SetUseWeightedEvents(kTRUE);
      purSc->SetStatisticOption(TEfficiency::kBUniform);
      cPurity->cd();
      DrawEfficiencyWithFixedYAxis(purSc, denPurity, 0.0, 1.0, "Purity");
      cPurity->SaveAs(state.folderPath + "/purity_delta_t" + Paths::ext_img.Data());
      if (rootOutFile)
      {
        TString scenDirName = TString("scenario_") + scenario;
        TDirectory *scDir = rootOutFile->GetDirectory(scenDirName);
        if (!scDir)
          scDir = rootOutFile->mkdir(scenDirName);
        if (scDir)
        {
          scDir->cd();
          purSc->Write("purity_delta_t_obj", TObject::kOverwrite);
          cPurity->Write("purity_delta_t_canvas", TObject::kOverwrite);
        }
      }
      delete cPurity;
    }

    // Purity in |Delta t| subsets for each active scenario.
    if (numPurity && denPurity)
    {
      std::vector<PuritySubset> puritySubsetsSc = CalculatePurityInSubsets(
          numPurity,
          denPurity,
          20.0,
          40);

      TCanvas *canvasPuritySubsetsSc = DrawPuritySubsets(puritySubsetsSc, Form("delta_t_%s", scenario.Data()));
      if (canvasPuritySubsetsSc)
      {
        canvasPuritySubsetsSc->SaveAs(state.folderPath + "/purity_subsets_delta_t" + Paths::ext_img.Data());
        if (rootOutFile)
        {
          TString scenDirName = TString("scenario_") + scenario;
          TDirectory *scDir = rootOutFile->GetDirectory(scenDirName);
          if (!scDir)
            scDir = rootOutFile->mkdir(scenDirName);
          if (scDir)
          {
            scDir->cd();
            canvasPuritySubsetsSc->Write("purity_subsets_delta_t_canvas", TObject::kOverwrite);
          }
        }
      }
    }

    SaveDeltaTComparisonPlot(state.folderPath, scenario, state.histsFittedSignal, state.histsReconstructed);

    if (rootOutFile)
    {
      TString scenDirName = TString("scenario_") + scenario;
      TDirectory *scDir = rootOutFile->GetDirectory(scenDirName);
      if (!scDir)
        scDir = rootOutFile->mkdir(scenDirName);
      if (scDir)
      {
        TDirectory *hRecDir = scDir->GetDirectory("hists_reconstructed");
        if (!hRecDir)
          hRecDir = scDir->mkdir("hists_reconstructed");
        SaveHistMap1D(hRecDir, state.histsReconstructed);

        TDirectory *hFitDir = scDir->GetDirectory("hists_fitted");
        if (!hFitDir)
          hFitDir = scDir->mkdir("hists_fitted");
        SaveHistMap1D(hFitDir, state.histsFittedSignal);

        TDirectory *h2Dir = scDir->GetDirectory("hists2d_fitted");
        if (!h2Dir)
          h2Dir = scDir->mkdir("hists2d_fitted");
        SaveHistMap2D(h2Dir, state.hists2DFittedSignal);

        scDir->cd();
        if (state.deltaTSignalTot)
          state.deltaTSignalTot->Write("deltaTSignalTot", TObject::kOverwrite);
        if (state.deltaTTot)
          state.deltaTTot->Write("deltaTTot", TObject::kOverwrite);
        if (state.histCounts)
          state.histCounts->Write("histCounts", TObject::kOverwrite);
      }
    }
  }

  SaveDeltaTComparisonPlot(folderPath, fOption, histsFittedSignal, histsReconstructed);

  if (rootOutFile)
  {
    TDirectory *primaryDir = rootOutFile->GetDirectory("primary");
    if (!primaryDir)
      primaryDir = rootOutFile->mkdir("primary");
    if (primaryDir)
    {
      TDirectory *hRecDir = primaryDir->GetDirectory("hists_reconstructed");
      if (!hRecDir)
        hRecDir = primaryDir->mkdir("hists_reconstructed");
      SaveHistMap1D(hRecDir, histsReconstructed);

      TDirectory *hFitDir = primaryDir->GetDirectory("hists_fitted");
      if (!hFitDir)
        hFitDir = primaryDir->mkdir("hists_fitted");
      SaveHistMap1D(hFitDir, histsFittedSignal);

      TDirectory *h2Dir = primaryDir->GetDirectory("hists2d_fitted");
      if (!h2Dir)
        h2Dir = primaryDir->mkdir("hists2d_fitted");
      SaveHistMap2D(h2Dir, hists2DFittedSignal);

      primaryDir->cd();
      if (deltaTSignalTot)
        deltaTSignalTot->Write("deltaTSignalTot", TObject::kOverwrite);
      if (deltaTTot)
        deltaTTot->Write("deltaTTot", TObject::kOverwrite);
      if (histCounts)
        histCounts->Write("histCounts", TObject::kOverwrite);
      if (canvaEff)
        canvaEff->Write("efficiency_canvas_latest", TObject::kOverwrite);
      if (canvaPurity)
        canvaPurity->Write("purity_canvas_latest", TObject::kOverwrite);
    }

    rootOutFile->Write();
    rootOutFile->Close();
    std::cout << "Saved ROOT cache file: " << g_rootOutputFileName << std::endl;
  }
}

Double_t signal_vs_bcg_v2::CalculatePurity(Int_t signal, Int_t total) const
{
  if (total == 0)
    return 0.0; // Unikaj dzielenia przez zero
  return static_cast<Double_t>(signal) / total;
}

Double_t signal_vs_bcg_v2::CalculateEfficiency(Int_t signal, Int_t total) const
{
  if (total == 0)
    return 0.0; // Unikaj dzielenia przez zero

  return static_cast<Double_t>(signal) / total;
}

void signal_vs_bcg_v2::FolderManagement(TString folderPath) const
{
  // Usuń folder, jeśli istnieje
  system(Form("rm -rf %s", folderPath.Data()));

  // Utwórz folder
  system(Form("mkdir -p %s", folderPath.Data()));
}

CutCount signal_vs_bcg_v2::CountEventsAroundCut(TH1 *hist, Double_t cutValue)
{
  CutCount result;
  result.totalEvents = hist->GetEntries();
  result.eventsBelow = 0;
  result.eventsAbove = 0;

  if (result.totalEvents == 0)
  {
    std::cerr << "WARNING: Histogram is empty!" << std::endl;
    return result;
  }

  // Znajdź bin odpowiadający wartości cięcia
  Int_t cutBin = hist->FindBin(cutValue);

  // Policz zdarzenia poniżej cięcia
  for (Int_t bin = 0; bin <= cutBin; bin++)
  {
    result.eventsBelow += hist->GetBinContent(bin);
  }

  // Policz zdarzenia powyżej (włącznie z binem cięcia)
  for (Int_t bin = cutBin + 1; bin <= hist->GetNbinsX() + 1; bin++)
  {
    result.eventsAbove += hist->GetBinContent(bin);
  }

  return result;
}

// Przykład użycia:
void signal_vs_bcg_v2::PrintCutCount(const TString &histName, Double_t cutValue, const CutCount &count)
{
  std::cout << "Histogram: " << histName << ", Cut at: " << cutValue << std::endl;
  std::cout << "  Events below: " << count.eventsBelow
            << " (" << (100.0 * count.eventsBelow / count.totalEvents) << "%)" << std::endl;
  std::cout << "  Events above: " << count.eventsAbove
            << " (" << (100.0 * count.eventsAbove / count.totalEvents) << "%)" << std::endl;
  std::cout << "  Total: " << count.totalEvents << std::endl;
}

void signal_vs_bcg_v2::DrawLabelOnHisto(std::vector<TString> labels)
{
  for (const auto &label : labels)
  {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(11); // Align left
    latex.DrawLatex(0.15, 0.92, label + "\n");
  }
}

void signal_vs_bcg_v2::AddDataToPCA(Double_t *data)
{
  pca->AddRow(data);
}

Double_t signal_vs_bcg_v2::GetPCASlope()
{
  pca->MakePrincipals();
  const TVectorD *principalValues = pca->GetEigenValues();
  const TMatrixD *principalVectors = pca->GetEigenVectors();

  Double_t v1x = (*principalVectors)(0, 0);
  Double_t v1y = (*principalVectors)(1, 0);

  Double_t lambda1 = (*principalValues)[0];

  // Skaluj wektor przez std::sqrt(lambda) aby uzyskać prawdziwy slope
  Double_t scaledV1x = v1x * TMath::Sqrt(lambda1);
  Double_t scaledV1y = v1y * TMath::Sqrt(lambda1);

  if (TMath::Abs(scaledV1x) < 1e-10)
  {
    std::cout << "WARNING: scaledV1x near zero!" << std::endl;
    return -10.0 / 9.0; // fallback
  }

  Double_t slope = scaledV1y / scaledV1x;

  std::cout << "GetPCASlope: v1x=" << v1x << ", v1y=" << v1y
            << ", lambda1=" << lambda1 << std::endl;
  std::cout << "  scaledV1x=" << scaledV1x << ", scaledV1y=" << scaledV1y
            << ", slope=" << slope << std::endl;

  return slope;
}

void signal_vs_bcg_v2::WidthOfCorrelatedHist(Double_t *means, Double_t *sigmas, Double_t *vLong, Double_t *vTransv)
{
  pca->MakePrincipals();
  const TVectorD *principalValues = pca->GetEigenValues();   // Wariancje λ
  const TMatrixD *principalVectors = pca->GetEigenVectors(); // Wektory v (znormalizowane)
  const TVectorD *dataMean = pca->GetMeanValues();
  const TVectorD *dataSigma = pca->GetSigmas();

  means[0] = (*dataMean)[0]; // Mean x
  means[1] = (*dataMean)[1]; // Mean y

  // Znormalizowane wektory własne
  Double_t v1x = (*principalVectors)(0, 0); // Pierwsza os (długa)
  Double_t v1y = (*principalVectors)(1, 0);
  Double_t v2x = (*principalVectors)(0, 1); // Druga os (poprzeczna)
  Double_t v2y = (*principalVectors)(1, 1);

  // Wariancje (eigenvalues)
  Double_t lambda1 = (*principalValues)[0]; // Większa wariancja
  Double_t lambda2 = (*principalValues)[1]; // Mniejsza wariancja

  // WŁAŚCIWE użycie wariancji - skalowanie wektora:
  Double_t scaledV1x = v1x * TMath::Sqrt(lambda1);
  Double_t scaledV1y = v1y * TMath::Sqrt(lambda1);

  Double_t slope = scaledV1y / scaledV1x; // Rzeczywisty slope w przestrzeni danych

  // Normalizuj wektory do długości 1
  Double_t normLong = TMath::Sqrt(scaledV1x * scaledV1x + scaledV1y * scaledV1y);
  vLong[0] = scaledV1x / normLong;
  vLong[1] = scaledV1y / normLong;

  // Wektor prostopadły
  vTransv[0] = -scaledV1y / normLong;
  vTransv[1] = scaledV1x / normLong;

  // Sigmy (odchylenia standardowe)
  Double_t sigmaLong = TMath::Sqrt(lambda1);   // Wzdłuż głównej osi
  Double_t sigmaTransv = TMath::Sqrt(lambda2); // Prostopadle do osi

  sigmas[0] = sigmaTransv; // Sigma poprzeczna (szerokość)
  sigmas[1] = sigmaLong;   // Sigma wzdłużna

  std::cout << "\n=== PCA RESULTS WITH VARIANCES ===" << std::endl;
  std::cout << "Lambda1 (Long variance): " << lambda1 << std::endl;
  std::cout << "Lambda2 (Transv variance): " << lambda2 << std::endl;
  std::cout << "v1x (raw): " << v1x << ", v1y (raw): " << v1y << std::endl;
  std::cout << "scaledV1x: " << scaledV1x << ", scaledV1y: " << scaledV1y << std::endl;
  std::cout << "Trend line slope (from PCA): " << slope << std::endl;
  std::cout << "Mean X: " << means[0] << ", Mean Y: " << means[1] << std::endl;
  std::cout << "Sigma Transverse: " << sigmas[0] << std::endl;
  std::cout << "Sigma Longitudinal: " << sigmas[1] << std::endl;
  std::cout << "Long vector (normalized): (" << vLong[0] << ", " << vLong[1] << ")" << std::endl;
  std::cout << "Transv vector (normalized): (" << vTransv[0] << ", " << vTransv[1] << ")" << std::endl;
  std::cout << "===================================\n"
            << std::endl;
}

// Funkcja do tworzenia wykresu RMS vs X
TGraphErrors *signal_vs_bcg_v2::CreateRMSProfile(TH2 *h2D, const char *name, const char *title)
{
  Int_t nbinsX = h2D->GetNbinsX();

  std::vector<Double_t> xPoints, yPoints, xErrors, yErrors;

  // Iteruj po binach X
  for (Int_t binX = 1; binX <= nbinsX; binX++)
  {
    // Weź projekcję Y dla danego binu X (slice)
    TH1D *projY = h2D->ProjectionY(Form("_py_%d", binX), binX, binX);

    // Pomiń puste biny lub z małą statystyką
    if (projY->GetEntries() < 10)
    { // Wymaga minimum 10 wpisów do sensownego fitu
      delete projY;
      continue;
    }

    // Środek binu X
    Double_t x = h2D->GetXaxis()->GetBinCenter(binX);
    Double_t xErr = h2D->GetXaxis()->GetBinWidth(binX) / 2.0;

    // Parametry wstępne dla fitu Gaussa
    Double_t mean = projY->GetMean();
    Double_t rms = projY->GetRMS();
    Double_t maxVal = projY->GetMaximum();

    // Definicja funkcji Gaussa
    TF1 *gausFit = new TF1(Form("gaus_%d", binX), "gaus", mean - 3 * rms, mean + 3 * rms);
    gausFit->SetParameters(maxVal, mean, rms);

    // Fit z opcją "Q" (quiet), "R" (range), "S" (return fit result)
    TFitResultPtr fitResult = projY->Fit(gausFit, "QRS");

    Double_t sigma = 0.0;
    Double_t sigmaError = 0.0;

    // Sprawdź czy fit się powiódł
    if (fitResult->IsValid() && fitResult->Status() == 0)
    {
      sigma = fitResult->Parameter(2);     // Sigma z fitu
      sigmaError = fitResult->ParError(2); // Błąd sigmy

      // Dodatkowa walidacja: sigma musi być sensowna
      if (sigma <= 0 || sigma > 10 * rms || sigmaError / sigma > 0.5)
      {
        // Jeśli fit jest zły, użyj RMS jako fallback
        sigma = rms;
        sigmaError = projY->GetRMSError();
      }
    }
    else
    {
      // Jeśli fit się nie powiódł, użyj RMS
      sigma = rms;
      sigmaError = projY->GetRMSError();
    }

    xPoints.push_back(x);
    yPoints.push_back(sigma);
    xErrors.push_back(xErr);
    yErrors.push_back(sigmaError);

    delete gausFit;
    delete projY;
  }

  // Stwórz TGraphErrors
  TGraphErrors *graph = new TGraphErrors(xPoints.size(),
                                         xPoints.data(),
                                         yPoints.data(),
                                         xErrors.data(),
                                         yErrors.data());

  graph->SetName(name);
  graph->SetTitle(title);
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);

  return graph;
}

TCanvas *signal_vs_bcg_v2::CreateCanvasWithProfiles(TH2 *h2D, const TString &name,
                                                    Bool_t drawMeanProfile,
                                                    Bool_t drawSigmaProfile)
{
  // Oblicz liczbę padów
  Int_t nPads = 1; // Zawsze histogram 2D
  if (drawMeanProfile)
    nPads++;
  if (drawSigmaProfile)
    nPads++;

  // Oblicz wysokość canvasu
  Int_t canvasHeight = 400 + (nPads - 1) * 250; // 400 dla 2D, 250 dla każdego profilu

  // Stwórz canvas
  TCanvas *c = new TCanvas(name, name, 800, canvasHeight);

  // Proporcje padów (od DOŁU do GÓRY)
  Double_t padBottom = 0.0;
  Double_t padTop = 0.0;

  TString yTitleh2D = h2D->GetYaxis()->GetTitle();
  TRegexp pattern("\\[.*\\]"); // Wzorzec dla [cokolwiek]

  TString match = yTitleh2D(pattern);

  TString yTitle2DReplaced = yTitleh2D.ReplaceAll(match, "").ReplaceAll("(", "").ReplaceAll(")", "").Strip(),
          yTitleProf = Form("#mu(%s) %s", yTitle2DReplaced.Data(), match.Data()),
          yTitleProfSigma = Form("#sigma(%s) %s", yTitle2DReplaced.Data(), match.Data());

  // === PAD 3 (najniższy): Sigma z fitu Gaussa ===
  if (drawSigmaProfile)
  {
    padBottom = 0.0;
    padTop = 0.25; // 25% wysokości

    TPad *pad3 = new TPad("pad3", "pad3", 0.0, padBottom, 1.0, padTop);
    pad3->SetTopMargin(0.02);
    pad3->SetBottomMargin(0.25); // Miejsce na etykiety
    pad3->SetRightMargin(0.15);
    pad3->Draw();
    pad3->cd();

    TGraphErrors *grRMS = CreateRMSProfile(h2D,
                                           Form("rms_%s", name.Data()),
                                           "");

    grRMS->SetMarkerStyle(20);
    grRMS->SetMarkerColor(kBlue);
    grRMS->SetLineColor(kBlue);
    grRMS->SetMarkerSize(0.8);

    grRMS->GetXaxis()->SetTitle(h2D->GetXaxis()->GetTitle());
    grRMS->GetYaxis()->SetTitle(yTitleProfSigma);

    grRMS->GetXaxis()->SetTitleSize(0.08);
    grRMS->GetXaxis()->SetLabelSize(0.08);
    grRMS->GetYaxis()->SetTitleSize(0.08);
    grRMS->GetYaxis()->SetLabelSize(0.08);
    grRMS->GetYaxis()->SetTitleOffset(0.5);

    grRMS->GetXaxis()->SetLimits(h2D->GetXaxis()->GetXmin(),
                                 h2D->GetXaxis()->GetXmax());

    grRMS->Draw("AP");

    c->cd(); // Wróć do głównego canvasu
  }

  // === PAD 2 (środkowy): TProfile (średnia) ===
  if (drawMeanProfile)
  {
    if (drawSigmaProfile)
    {
      padBottom = 0.25; // Zaczyna się gdzie kończy pad3
      padTop = 0.5;     // 25% wysokości
    }
    else
    {
      padBottom = 0.0;
      padTop = 0.3; // 30% jeśli nie ma sigma profile
    }

    TPad *pad2 = new TPad("pad2", "pad2", 0.0, padBottom, 1.0, padTop);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(drawSigmaProfile ? 0.02 : 0.25);
    pad2->SetRightMargin(0.15);
    pad2->Draw();
    pad2->cd();

    TProfile *prof = h2D->ProfileX("_pfx", 1, -1, "");
    prof->SetTitle("");
    prof->SetMarkerStyle(20);
    prof->SetMarkerColor(kRed);
    prof->SetLineColor(kRed);
    prof->SetMarkerSize(0.8);

    prof->GetYaxis()->SetTitle(yTitleProf);
    prof->GetYaxis()->SetTitleSize(0.08);
    prof->GetYaxis()->SetLabelSize(0.08);
    prof->GetYaxis()->SetTitleOffset(0.5);

    if (drawSigmaProfile)
    {
      prof->GetXaxis()->SetLabelSize(0);
      prof->GetXaxis()->SetTitleSize(0);
    }
    else
    {
      prof->GetXaxis()->SetTitle(h2D->GetXaxis()->GetTitle());
      prof->GetXaxis()->SetTitleSize(0.08);
      prof->GetXaxis()->SetLabelSize(0.08);
    }

    prof->GetXaxis()->SetRangeUser(h2D->GetXaxis()->GetXmin(),
                                   h2D->GetXaxis()->GetXmax());

    prof->GetYaxis()->SetRangeUser(h2D->GetYaxis()->GetXmin(),
                                   h2D->GetYaxis()->GetXmax());
    prof->Draw("E1");

    TLine *line = new TLine(h2D->GetXaxis()->GetXmin(), 0,
                            h2D->GetXaxis()->GetXmax(), 0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->Draw();

    c->cd(); // Wróć do głównego canvasu
  }

  // === PAD 1 (najwyższy): Histogram 2D ===
  if (drawMeanProfile && drawSigmaProfile)
  {
    padBottom = 0.5; // Zaczyna się gdzie kończy pad2
    padTop = 1.0;    // 50% wysokości
  }
  else if (drawMeanProfile || drawSigmaProfile)
  {
    padBottom = drawMeanProfile ? 0.3 : 0.25;
    padTop = 1.0; // 70% lub 75%
  }
  else
  {
    padBottom = 0.0;
    padTop = 1.0; // 100%
  }

  TPad *pad1 = new TPad("pad1", "pad1", 0.0, padBottom, 1.0, padTop);
  pad1->SetBottomMargin(nPads > 1 ? 0.02 : 0.12);
  pad1->SetRightMargin(0.15);
  pad1->SetLogz(1);
  pad1->Draw();
  pad1->cd();

  h2D->Draw("COLZ");
  if (nPads > 1)
  {
    h2D->GetXaxis()->SetLabelSize(0);
    h2D->GetXaxis()->SetTitleSize(0);
  }

  c->cd();
  return c;
}

// Funkcja główna - oblicza czystość w podzbiorach |ΔT| < limit
std::vector<PuritySubset> signal_vs_bcg_v2::CalculatePurityInSubsets(
    TH1 *hist_signal,
    TH1 *hist_total,
    Double_t max_limit,
    Int_t n_subsets)
{
  if (!hist_signal || !hist_total)
  {
    std::cerr << "ERROR: Invalid histograms in CalculatePurityInSubsets" << std::endl;
    return {};
  }

  const Double_t signalIntegral = hist_signal->Integral(0, hist_signal->GetNbinsX() + 1);
  const Double_t totalIntegral = hist_total->Integral(0, hist_total->GetNbinsX() + 1);
  if (signalIntegral <= 0.0 || totalIntegral <= 0.0)
  {
    std::cerr << "ERROR: Empty histogram integrals in CalculatePurityInSubsets (signal="
              << signalIntegral << ", total=" << totalIntegral << ")" << std::endl;
    return {};
  }

  std::vector<PuritySubset> results;

  // Oblicz limity dla podzbiorów (np. 2, 4, 6, 8, 10 ps)
  Double_t limit_step = max_limit / n_subsets;

  for (Int_t i = 1; i <= n_subsets; i++)
  {
    Double_t current_limit = i * limit_step;

    // Oblicz czystość w symetrycznym, kumulacyjnym przedziale |ΔT| < limit.
    Int_t bin_low = hist_signal->FindBin(-current_limit);
    Int_t bin_high = hist_signal->FindBin(current_limit);

    Double_t signal_count_d = hist_signal->Integral(bin_low, bin_high);
    Double_t total_count_d = hist_total->Integral(bin_low, bin_high);

    Int_t signal_count = static_cast<Int_t>(signal_count_d);
    Int_t total_count = static_cast<Int_t>(total_count_d);

    Double_t purity = 0.0;
    Double_t purity_error = 0.0;

    if (total_count_d > 0.0)
    {
      purity = signal_count_d / total_count_d;

      if (purity > 1.0)
      {
        std::cout << "WARNING: purity > 1 for |DeltaT| < " << current_limit
                  << " (signal=" << signal_count_d << ", total=" << total_count_d << ")" << std::endl;
      }

      // Błąd Bayesowski dla czystości (binomial)
      // σ_purity = std::sqrt(p*(1-p)/N)
      if (total_count_d > 1.0)
      {
        purity_error = TMath::Sqrt(TMath::Abs(purity * (1.0 - purity)) / total_count_d);
      }
    }

    PuritySubset subset;
    subset.deltaT_limit = current_limit;
    subset.deltaT_limit_error = limit_step / 2.0;
    subset.signal_events = signal_count;
    subset.total_events = total_count;
    subset.purity = purity;
    subset.purity_error = purity_error;

    results.push_back(subset);

    std::cout << "Subset |ΔT| < " << current_limit << " ps: "
              << "Signal=" << signal_count << ", Total=" << total_count
              << ", Purity=" << (purity * 100.0) << "% ± "
              << (purity_error * 100.0) << "%" << std::endl;
  }

  return results;
}

// Funkcja do rysowania wykresu czystości vs |ΔT| limit
TCanvas *signal_vs_bcg_v2::DrawPuritySubsets(
    const std::vector<PuritySubset> &subsets,
    const TString &name)
{
  if (subsets.empty())
  {
    std::cerr << "ERROR: Empty subsets in DrawPuritySubsets" << std::endl;
    return nullptr;
  }

  // Przygotuj dane do TGraphErrors
  std::vector<Double_t> limits, purities, limit_errors, purity_errors;

  for (const auto &subset : subsets)
  {
    limits.push_back(subset.deltaT_limit);
    purities.push_back(subset.purity);
    limit_errors.push_back(subset.deltaT_limit_error);
    purity_errors.push_back(subset.purity_error);
  }

  // Stwórz TGraphErrors
  TGraphErrors *graph_purity = new TGraphErrors(
      limits.size(),
      limits.data(),
      purities.data(),
      limit_errors.data(),
      purity_errors.data());

  graph_purity->SetName(Form("purity_%s", name.Data()));
  graph_purity->SetTitle("Purity in |#Deltat| regions");
  graph_purity->SetMarkerStyle(20);
  graph_purity->SetMarkerColor(kBlue);
  graph_purity->SetMarkerSize(1.0);
  graph_purity->SetLineColor(kBlue);
  graph_purity->SetLineWidth(2);

  // Sformułuj osie
  graph_purity->GetXaxis()->SetTitle("|#Deltat| regions [#tau_{S}]");
  graph_purity->GetYaxis()->SetTitle("Purity");
  graph_purity->GetXaxis()->SetTitleSize(0.05);
  graph_purity->GetYaxis()->SetTitleSize(0.05);
  graph_purity->GetXaxis()->SetLabelSize(0.04);
  graph_purity->GetYaxis()->SetLabelSize(0.04);

  // Stwórz canvas
  TCanvas *canvas = new TCanvas(Form("c_purity_%s", name.Data()),
                                Form("Purity subsets: %s", name.Data()),
                                800, 600);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);

  // Rysuj
  graph_purity->Draw("APL");
  Double_t maxPurity = 0.0;
  for (const auto &subset : subsets)
  {
    if (subset.purity > maxPurity)
      maxPurity = subset.purity;
  }
  graph_purity->GetYaxis()->SetRangeUser(0.0, TMath::Max(1.0, 1.1 * maxPurity));
  gPad->Update();

  // Dodaj siatkę
  gPad->SetGrid(1, 1);

  // Dodaj linię referencyjną na 100%
  TLine *line_max = new TLine(
      limits.front(), 1.0,
      limits.back(), 1.0);
  line_max->SetLineStyle(2);
  line_max->SetLineColor(kGray);
  line_max->Draw("SAME");

  // Dodaj tabelkę z liczbami
  // TPaveText *stats_box = new TPaveText(0.15, 0.15, 0.55, 0.5, "NDC");
  // stats_box->SetFillColor(kWhite);
  // stats_box->SetBorderSize(1);
  // stats_box->SetTextAlign(12);
  // stats_box->SetTextSize(0.03);
  // stats_box->SetTextFont(42);

  // stats_box->AddText("Purity per |#Delta t| subset:");
  // stats_box->AddText("");

  // for (size_t i = 0; i < subsets.size(); i++)
  // {
  //   stats_box->AddText(Form(
  //       "|#Delta t| < %.1f ps: %.2f%% (%d/%d)",
  //       subsets[i].deltaT_limit,
  //       subsets[i].purity * 100.0,
  //       subsets[i].signal_events,
  //       subsets[i].total_events));
  // }

  // stats_box->Draw();
  // gPad->Update();

  return canvas;
}

// Dodatkowo: Funkcja do porównania purity w podzbiorach symetrycznych
std::map<Double_t, Double_t> signal_vs_bcg_v2::ComparePuritySymmetric(
    TH1 *hist_signal,
    TH1 *hist_total,
    const std::vector<Double_t> &limits)
{
  std::map<Double_t, Double_t> purity_map;

  for (Double_t limit : limits)
  {
    Int_t bin_low = hist_signal->FindBin(-limit);
    Int_t bin_high = hist_signal->FindBin(limit);

    Int_t signal_count = static_cast<Int_t>(hist_signal->Integral(bin_low, bin_high));
    Int_t total_count = static_cast<Int_t>(hist_total->Integral(bin_low, bin_high));

    Double_t purity = (total_count > 0)
                          ? static_cast<Double_t>(signal_count) / total_count
                          : 0.0;

    purity_map[limit] = purity;
  }

  return purity_map;
}