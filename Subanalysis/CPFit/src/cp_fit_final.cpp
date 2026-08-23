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
#include <TBufferJSON.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <boost/progress.hpp>
#include <TF1.h>
#include <set>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <BRCorrectionFactors.h>

#include "../inc/fit_setter.h"
#include "../../RegenerationAnalysis/inc/regen_frac_fit.h"

#include "../inc/cpfit.hpp"

int cp_fit_final(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj, ConfigWatcher &cfgWatcher)
{
  gErrorIgnoreLevel = kBreak;

  int mcFiles = 0, dataFiles = 0;
  double lumiFactor = 1.0;

  mcFiles = Utils::CountFiles(&chain, true);
  dataFiles = Utils::CountFiles(&chain, false);

  lumiFactor = 1.07906;//(dataFiles > 0) ? static_cast<double>(dataFiles) / static_cast<double>(mcFiles) : 1.0;

  std::cout << "INFO: Number of MC files: " << mcFiles << std::endl;
  std::cout << "INFO: Number of Data files: " << dataFiles << std::endl;
  std::cout << "INFO: Luminosity factor: " << lumiFactor << std::endl;

  // =============================================================================
  KLOE::BaseKinematics
      baseKin;
  Int_t
      file_num;
  TFile
      *file_corr;
  // =============================================================================
  const TString cpfit_res_dir = Paths::cpfit_dir + Paths::result_dir;

  Double_t *eff_vals;

  if (check_corr == true)
  {
    file_corr = new TFile("../Efficiency_analysis/correction_factor.root");
    TGraphAsymmErrors *eff_signal = (TGraphAsymmErrors *)file_corr->Get("correction_factor");
    eff_vals = eff_signal->GetY();

    delete eff_signal;
  }

  // =============================================================================
  // Preparation of the Three pi0 hist object to calculate regeneration weights
  TFile *fFileWeightsThreePi0 = TFile::Open("/data/ssd/gamrat/python-kloe-analysis/scripts/results/control_sample_corr_factors/control_sample_corr_factors.root", "READ");

  TCanvas *c = (TCanvas *)fFileWeightsThreePi0->Get("three_pi0/correction_factors/cCorr");
  TH1 *threePi0WeightsHist = (TH1 *)c->FindObject("correction_factors")->Clone();

  // Calculation of average weighting factor for three pi0 events
  TFitResultPtr threePi0WeightResult = threePi0WeightsHist->Fit("pol0", "QS");
  double avgThreePi0Weight = 0.9375;//threePi0WeightResult->Parameter(0);
  std::cout << "INFO: Average weighting factor for three pi0 events: " << avgThreePi0Weight << std::endl;

  // Preparation of the Semileptonic hist object to calculate regeneration weights
  TFile *fFileWeightsSemileptonic = TFile::Open("/data/ssd/gamrat/python-kloe-analysis/scripts/results/control_sample_corr_factors/control_sample_corr_factors.root", "READ");

  TCanvas *ctest = (TCanvas *)fFileWeightsSemileptonic->Get("semileptonic/correction_factors/cCorr");
  TH1 *semileptonicWeightsHist = (TH1 *)ctest->FindObject("correction_factors")->Clone();

  // Calculation of average weighting factor for semileptonic events
  TFitResultPtr semileptonicWeightResult = semileptonicWeightsHist->Fit("pol0", "QS");
  double avgSemileptonicWeight = 0.9272;//semileptonicWeightResult->Parameter(0);
  std::cout << "INFO: Average weighting factor for semileptonic events: " << avgSemileptonicWeight << std::endl;

  // Preparation of signal weights
  TFile *fFileWeightsSignal = TFile::Open("/data/ssd/gamrat/python-kloe-analysis/scripts/results/control_sample_corr_factors/combined_signal_correction_factors.root", "READ");

  TCanvas *cSignal = (TCanvas*)fFileWeightsSignal->Get("CorrFactors");
  TH1* signalWeightsHist = (TH1*)cSignal->FindObject("CombinedCorrFactor")->Clone();
  

  // Preparation of the RegenerationFractionFit object to calculate regeneration weights
  std::string weightsFilePath = (std::string)Paths::regen_analysis_dir + (std::string)Paths::result_dir + "regeneration_analysis_results.root";
  RegenerationFractionFit regen_weights(weightsFilePath);

  // ===========================================================================

  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> mcflag(reader, "mcflag");
  TTreeReaderValue<Int_t> mctruth(reader, "mctruth");

  TTreeReaderValue<Double_t> Chi2SignalKinFit(reader, "Chi2SignalKinFit");
  TTreeReaderValue<Double_t> minv4gam(reader, "minv4gam");

  TTreeReaderValue<Double_t> Qmiss(reader, "Qmiss");

  TTreeReaderValue<Double_t> Bx(reader, "Bx");
  TTreeReaderValue<Double_t> By(reader, "By");
  TTreeReaderValue<Double_t> Bz(reader, "Bz");

  TTreeReaderValue<Double_t> KaonChTimeCMMC(reader, "KaonChTimeCMMC");
  TTreeReaderValue<Double_t> KaonNeTimeCMMC(reader, "KaonNeTimeCMMC");
  TTreeReaderValue<Double_t> KaonChTimeCMSignalFit(reader, "KaonChTimeCMSignalFit");
  TTreeReaderValue<Double_t> KaonNeTimeCMSignalFit(reader, "KaonNeTimeCMSignalFit");
  TTreeReaderValue<Double_t> KaonChTimeCMBoostTriFit(reader, "KaonChTimeCMBoostTriFit");
  TTreeReaderValue<Double_t> KaonNeTimeCMBoostTriFit(reader, "KaonNeTimeCMBoostTriFit");


  TTreeReaderValue<Double_t> bestErrorSixGamma(reader, "bestErrorSixGamma");

  TTreeReaderArray<Double_t> Knerec(reader, "Knerec");
  TTreeReaderArray<Double_t> KnerecSix(reader, "KnerecSix");
  TTreeReaderArray<Double_t> KchrecClosest(reader, "KchrecClosest");

  TTreeReaderArray<Double_t> KnerecFit(reader, "KnerecFit");
  TTreeReaderArray<Double_t> KchrecFit(reader, "KchrecFit");
  TTreeReaderArray<Double_t> ipFit(reader, "ipFit");

  TTreeReaderArray<Double_t> Kchrec(reader, "Kchrec");
  TTreeReaderArray<Double_t> ip(reader, "ip");

  TTreeReaderArray<Double_t> trk1Fit(reader, "trk1Fit");
  TTreeReaderArray<Double_t> trk2Fit(reader, "trk2Fit");

  TTreeReaderArray<Double_t> pi01Fit(reader, "pi01Fit");
  TTreeReaderArray<Double_t> pi02Fit(reader, "pi02Fit");

  // ===========================================================================

  ///////////////////////////////////////
  FitSetter setter((std::string)Paths::cpfit_dir + "config/config.json");
  FitConfig cfg = setter.getFitConfig();
  Double_t
      x_min = cfg.deltaTConfig.xRange[0],
      x_max = cfg.deltaTConfig.xRange[1],
      res_deltaT = cfg.deltaTConfig.resolution;
  UInt_t
      nbins = Int_t((x_max - x_min) / res_deltaT);

  std::cout << "INFO: Fit settings loaded successfully." << std::endl;

  if (nbins % 2 == 0)
    nbins += 1; // Ensure an odd number of bins for symmetry around zero

  // Splits related to regeneration
  // Single mode: jeden parametr A_regen bez parameter_ranges → granice z regenerationSplit
  // Split mode:  4 parametry A_regen_far_left/near_left/near_right/far_right → granice z ich parameter_ranges
  bool regenSplitMode = cfg.parameters.count("A_regen_far_left") &&
                        cfg.parameters.at("A_regen_far_left").enabled;

  Double_t left_split, center_split, right_split;
  if (regenSplitMode)
  {
    left_split = cfg.parameters.at("A_regen_far_left").parameter_ranges[0][1];
    center_split = cfg.parameters.at("A_regen_near_left").parameter_ranges[0][1];
    right_split = cfg.parameters.at("A_regen_near_right").parameter_ranges[0][1];
  }
  else
  {
    left_split = cfg.regenerationSplit[0];
    center_split = cfg.regenerationSplit[1];
    right_split = cfg.regenerationSplit[2];
  }

  Double_t split[3] = {left_split, center_split, right_split};

  Bool_t regenerationExclusionFlag = cfg.regenerationExclusionFlag;

  std::cout << "INFO: Fit settings loaded successfully." << std::endl;

  // Pobranie zakresów wyświetlania (używane później w rysowaniu)
  const Double_t xMinRangeDisplay = cfg.deltaTConfig.xRangeDisplay[0];
  const Double_t xMaxRangeDisplay = cfg.deltaTConfig.xRangeDisplay[1];

  // Profile of relative errors
  TH2 *hEffErrorProfile = new TH2D("hEffErrorProfile", "Regen corr factor relative errors vs. deltaT;#Deltat [#tau_{S}];Relative error on weight [-]", nbins, x_min, x_max, 50, 0, 3);
  TH2 *hEffWeightVsRelErr = new TH2D("hEffWeightVsRelErr", "Regen corr factor weights vs. relative errors;Weight [-];Relative error on weight [-]", 50, 0, 7.0, 50, 0, 3);

  KLOE::interference event(mode, check_corr, nbins, x_min, x_max, split);

  // -------------------------------------------------------------------------
  // Zewnętrzne stałe skalujące histogram regeneracji (per split).
  // Kolejność: [0]=far_left, [1]=near_left, [2]=near_right, [3]=far_right.
  // Dostępne w całej funkcji – używane zarówno w chi2 jak i histogramach
  // prezentacyjnych (spójna propagacja błędów).
  std::array<KLOE::interference::RegenSplitScaling, 4> regenScaling = {{
      {1.0, 0.0}, // [0] far_left
      {1.0, 0.0}, // [1] near_left
      {1.0, 0.0}, // [2] near_right
      {1.0, 0.0}, // [3] far_right
  }};
  event.SetRegenScaling(regenScaling);
  // -------------------------------------------------------------------------
  // Get settings

  UInt_t num_of_vars = cfg.getNumOfEnabledParameters(); // Default number of variables (can be adjusted based on mode)
  ROOT::Math::Minimizer *minimum =
      ROOT::Math::Factory::CreateMinimizer(cfg.minimizer_type, cfg.minimizer_algorithm);

  // set tolerance , etc...
  minimum->SetMaxFunctionCalls(cfg.max_function_calls); // for Minuit/Minuit2
  minimum->SetTolerance(cfg.tolerance);
  minimum->SetPrintLevel(cfg.print_level);
  minimum->SetStrategy(cfg.strategy);

  std::cout << "INFO: Fit settings loaded successfully." << std::endl;

  std::map<TString, Int_t> param_index_map;
  Int_t par_num = 0; // Licznik parametrów faktycznie przekazanych do minimizera

  for (auto it = cfg.parameters.begin(); it != cfg.parameters.end(); ++it)
  {
    const TString &pName = it->first;
    const Parameter &p = it->second;

    if (!p.enabled)
    {
      // Parametr całkowicie pominięty w minimizerze
      event.SetParamIndex(pName, -1);
      param_index_map[pName] = -1;
      std::cout << "INFO: Parameter " << pName << " is DISABLED (omitted from fit)." << std::endl;
      continue;
    }

    // Rejestrujemy indeks w minimizerze
    event.SetParamIndex(pName, par_num);
    param_index_map[pName] = par_num;

    std::cout << "INFO: Parameter " << pName << " is ENABLED with initial value " << p.initial_value
              << ", step " << p.step
              << (p.fixed ? ", FIXED." : (p.limit_lower > 0.0 || p.limit_upper > 0.0 ? ", LIMITED." : ", UNLIMITED."))
              << std::endl;

    if (p.fixed)
    {
      // Parametr jest w fitowaniu, ale ma stałą wartość (nie zmienia się)
      minimum->SetFixedVariable(par_num, pName.Data(), p.initial_value);
    }
    else
    {
      auto lims = p.limits();
      if (p.limit_lower > 0.0 || p.limit_upper > 0.0)
      {
        minimum->SetLimitedVariable(par_num, pName.Data(), p.initial_value, p.step, lims[0], lims[1]);
      }
      else
      {
        minimum->SetVariable(par_num, pName.Data(), p.initial_value, p.step);
      }
    }
    par_num++; // Zwiększamy tylko jeśli parametr został dodany do minimum
  }

  // Czyścimy i wypełniamy powiązania kanał -> parametry
  event.channel_to_indices.clear();
  for (auto const &ch : KLOE::channName)
  {
    TString chName = ch.second;
    std::vector<std::pair<Double_t, Int_t>> range_based;

    for (auto const &pEntry : cfg.parameters)
    {
      const auto &p = pEntry.second;
      if (!p.enabled)
        continue; // ← FIX 1: pomiń wyłączone
      if (std::find(p.channels.begin(), p.channels.end(), std::string(chName.Data())) == p.channels.end())
        continue;

      Int_t idx = param_index_map[pEntry.first];
      if (!p.parameter_ranges.empty())
        range_based.push_back({p.parameter_ranges[0][0], idx}); // ← zbieramy z rangą
      else
        event.channel_to_indices[chName].push_back(idx);
    }

    // FIX 2: sortowanie po dolnej granicy zakresu → gwarantuje kolejność far_left, near_left, near_right, far_right
    if (!range_based.empty())
    {
      std::sort(range_based.begin(), range_based.end());
      for (const auto &pr : range_based)
        event.channel_to_indices[chName].push_back(pr.second);
    }
  }

  ROOT::Math::Functor minimized_function(&event, &KLOE::interference::interf_chi2, num_of_vars);
  minimum->SetFunction(minimized_function);
  ////////////////////////////////////

  TH1 *sig_total = new TH1D("sig_tot", ";#Delta t;Counts", nbins, x_min, x_max);
  TH1 *sig_pass = new TH1D("sig_pass", ";#Delta t;Counts", nbins, x_min, x_max);

  std::vector<Double_t> no_cuts_sig[2];

  Long64_t nentries = reader.GetEntries(true);

  boost::progress_display display(nentries);

  double u0 = -92.9524, v0 = 7.1644, su = 33.7633, sv = 42.0163, rho = 0.556;
  double n_sigma_cut = 2.5;

  // Definicja lambdy
  auto ellipse_cut = [=](double u, double v)
  {
    double du = u - u0;
    double dv = v - v0;
    double inv_rho2 = 1.0 / (1.0 - rho * rho);

    double dist2 = inv_rho2 * ((du * du) / (su * su) +
                               (dv * dv) / (sv * sv) -
                               2.0 * rho * (du * dv) / (su * sv));

    return std::sqrt(std::max(0.0, dist2)) > n_sigma_cut;
  };

  Long64_t event_data = reader.GetEntries("mcflag == 0"),
           event_mc = reader.GetEntries("mcflag == 1 && mctruth != -1 && mctruth != 0");

  // Luminosity scaling factor
  Double_t lumi_scaling_factor = event_data / static_cast<Double_t>(event_mc);
  // ---

  struct RegenKinVars
  {
    Double_t rhoCharged, rhoNeutral, rCharged, rNeutral, dt;
  };
  std::vector<RegenKinVars> regen_kin_vars;

  reader.SetEntry(0); // Ensure we start from the first entry

  while (reader.Next())
  {
    // Cut setting
    std::array<Double_t, 3> distNeutralCharged = {KchrecFit[6] - KnerecFit[6],
                                                  KchrecFit[7] - KnerecFit[7],
                                                  KchrecFit[8] - KnerecFit[8]},
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

    Double_t radius00 = std::sqrt(std::pow(Knerec[6] - *Bx, 2) +
                                  std::pow(Knerec[7] - *By, 2)),
             radiuspm = std::sqrt(std::pow(KchrecClosest[6] - *Bx, 2) +
                                  std::pow(KchrecClosest[7] - *By, 2)),
             zdist00 = std::abs(Knerec[8] - KchrecClosest[8]),
             zdistpm = std::abs(KchrecClosest[8] - *Bz);

    Double_t fiducialVolume = std::sqrt(std::pow(distNeutralCharged[0], 2) + std::pow(distNeutralCharged[1], 2)) < 2.05 && std::abs(distNeutralCharged[2]) < 2.45,
             fiducialVolumeClose = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

    TVector3 KchrecVec = {KchrecFit[0], KchrecFit[1], KchrecFit[2]};
    TVector3 trk1VecFit = {trk1Fit[0], trk1Fit[1], trk1Fit[2]};
    TVector3 trk2VecFit = {trk2Fit[0], trk2Fit[1], trk2Fit[2]};

    Double_t phiTrk1Angle = cos(trk1VecFit.Angle(KchrecVec)),
             phiTrk2Angle = cos(trk2VecFit.Angle(KchrecVec));

    Bool_t global_cut = ((fiducialVolume && (abs(phiTrk1Angle) < 0.8 || abs(phiTrk2Angle) < 0.8)) || !fiducialVolume) && ((fiducialVolumeClose && rho > 1.5) || !fiducialVolumeClose) && ellipse_cut(*minv4gam - PhysicsConstants::mK0, KnerecSix[5] - PhysicsConstants::mK0);

    // Regeneration cut
    TVector3 ipVector, neutralVtxVector, chargedVtxVector;

    // Set up geometry vectors
    ipVector.SetXYZ(0, 0, 0); // ip[0], ip[1], ip[2]);
    neutralVtxVector.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);
    chargedVtxVector.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);

    // Calculate radius vectors
    TVector3 ipToCharged = chargedVtxVector - ipVector;
    TVector3 ipToNeutral = neutralVtxVector - ipVector;

    // Cylindrical radius (rho) to charged and neutral vertices
    double rhoCharged = ipToCharged.Perp();
    double rhoNeutral = ipToNeutral.Perp();

    // Spherical radius to charged and neutral vertices
    double rCharged = ipToCharged.Mag();
    double rNeutral = ipToNeutral.Mag();

    std::map<HistogramType, double> radiusMap = {
        {HistogramType::RHO_CHARGED, rhoCharged},
        {HistogramType::RHO_NEUTRAL, rhoNeutral},
        {HistogramType::R_CHARGED, rCharged},
        {HistogramType::R_NEUTRAL, rNeutral}};

    Double_t ch_Spherical_Mean = 10.4941, ch_Spherical_Sigma = 0.957544;
    Double_t ne_Spherical_Mean = 10.3769, ne_Spherical_Sigma = 1.23898;
    Double_t ch_Cylindrical_Mean = 4.84397, ch_Cylindrical_Sigma = 0.877508;
    Double_t ne_Cylindrical_Mean = 4.63652, ne_Cylindrical_Sigma = 0.703466;

    Double_t beamPipeBound = 5; // Przykładowa granica między obszarem cylindra a sferą (do dostosowania)

    Bool_t regeneration_cut = true; // Domyślnie przepuszczamy wszystkie zdarzenia, jeśli flaga jest wyłączona

    if (regenerationExclusionFlag)
    {
      Bool_t regeneration_cut_ch = false, regeneration_cut_ne = false;

      if (rhoCharged > beamPipeBound)
        regeneration_cut_ch = std::abs(rCharged - ch_Spherical_Mean) < 1.5 * ch_Spherical_Sigma;

      if (rhoNeutral > beamPipeBound)
        regeneration_cut_ne = std::abs(rNeutral - ne_Spherical_Mean) < 1.5 * ne_Spherical_Sigma;

      regeneration_cut = !(regeneration_cut_ch || regeneration_cut_ne);
    }

    baseKin.Dtmc = *KaonChTimeCMMC - *KaonNeTimeCMMC;
    baseKin.Dtboostlor = *KaonChTimeCMSignalFit - *KaonNeTimeCMSignalFit;

    double Dttri = *KaonChTimeCMBoostTriFit - *KaonNeTimeCMBoostTriFit;

    if (*mcflag == 1)
    {
      if (*mctruth == 1 || *mctruth == 0)
      {
        no_cuts_sig[0].push_back(baseKin.Dtmc);
        no_cuts_sig[1].push_back(baseKin.Dtboostlor);
      }

      if (global_cut)
      {
        if (*mctruth == 1)
        {
          event.time_diff_gen.push_back(baseKin.Dtmc);
          event.time_diff["Signal"].push_back(baseKin.Dtboostlor);
          double weight = signalWeightsHist->Interpolate(Dttri);

          // if (weight <= 0 || std::isnan(weight) || std::isinf(weight))
            weight = 1;

          event.signal_event_weights.push_back(weight);

          sig_pass->Fill(baseKin.Dtboostlor);
        }

        if (*mctruth == 2)
        {
          event.time_diff["Regeneration"].push_back(baseKin.Dtboostlor);

          if (cfg.regenShapeCorrection.enabled)
            event.regen_event_weights.push_back(regen_weights.GetContinuousRegenerationWeight(baseKin.Dtboostlor, radiusMap));
        }

        if (*mctruth == 3)
        {
          event.time_diff["Omega"].push_back(baseKin.Dtboostlor);
        }

        if (*mctruth == 4)
        {
          event.time_diff["3pi0"].push_back(baseKin.Dtboostlor);
          event.three_pi0_weights.push_back(1.0);//threePi0WeightsHist->Interpolate(Kchrec[5]));
        }

        if (*mctruth == 5)
        {
          event.time_diff["Semileptonic"].push_back(baseKin.Dtboostlor);
          event.semileptonic_weights.push_back(1.0);//semileptonicWeightsHist->Interpolate(*minv4gam));
        }

        if (*mctruth == 6)
        {
          event.time_diff["Other"].push_back(baseKin.Dtboostlor);
        }
      }
    }

    if (*mcflag == 0 && global_cut)
    {
      event.time_diff["Data"].push_back(baseKin.Dtboostlor);
    }

    ++display;
  }

  minimum->Minimize();

  std::vector<Double_t> par(num_of_vars), parErr(num_of_vars);

  // Wypełnij par/parErr tylko dla aktywnych parametrów
  for (const auto &pm : param_index_map)
  {
    if (pm.second < 0)
      continue;
    par[pm.second] = minimum->X()[pm.second];
    parErr[pm.second] = minimum->Errors()[pm.second];
  }

  std::cout << "INFO: Parameters" << std::endl;
  std::cout << "------------------------" << std::endl;
  for (const auto &pm : param_index_map)
  {
    if (pm.second < 0)
      continue;
    std::cout << pm.first << ": " << par[pm.second] << " ± " << parErr[pm.second] << std::endl;
  }

  Double_t sum_of_events = 0.;
  std::map<TString, Double_t> fractions;

  for (const auto &c : KLOE::channName)
    sum_of_events += event.time_diff[c.second].size();

  for (const auto &c : KLOE::channName)
    fractions[c.second] = 100 * event.time_diff[c.second].size() / sum_of_events;
  std::ofstream myfile_num;
  myfile_num.open(cpfit_res_dir + "num_of_events.csv");
  myfile_num << "Channel,Number of events,Fraction\n";
  myfile_num << "Signal," << event.time_diff["Signal"].size() << "," << fractions["Signal"] << "%,\n";
  myfile_num << "Regeneration," << event.time_diff["Regeneration"].size() << "," << fractions["Regeneration"] << "%,\n";
  myfile_num << "Omega," << event.time_diff["Omega"].size() << "," << fractions["Omega"] << "%,\n";
  myfile_num << "Three," << event.time_diff["3pi0"].size() << "," << fractions["3pi0"] << "%,\n";
  myfile_num << "Semi," << event.time_diff["Semileptonic"].size() << "," << fractions["Semileptonic"] << "%,\n";
  myfile_num << "Other bcg," << event.time_diff["Other"].size() << "," << fractions["Other"] << "%,\n";
  myfile_num.close();

  for (UInt_t i = 0; i < no_cuts_sig[1].size(); i++)
  {
    sig_total->Fill(no_cuts_sig[1][i]);
  }

  std::cout << "Integral of Regeneration: " << event.getFracHistogram("Regeneration")->Integral(0, nbins + 1) << std::endl;
  std::cout << "Integral of Omega: " << event.getFracHistogram("Omega")->Integral(0, nbins + 1) << std::endl;
  std::cout << "Integral of 3pi0: " << event.getFracHistogram("3pi0")->Integral(0, nbins + 1) << std::endl;
  std::cout << "Integral of Semileptonic: " << event.getFracHistogram("Semileptonic")->Integral(0, nbins + 1) << std::endl;
  std::cout << "Integral of Other: " << event.getFracHistogram("Other")->Integral(0, nbins + 1) << std::endl;
  std::cout << "Integral of Signal: " << event.getFracHistogram("Signal")->Integral(0, nbins + 1) << std::endl;
  std::cout << "Integral of MC sum: " << event.getFracHistogram("MC sum")->Integral(0, nbins + 1) << std::endl;

  // ─── Purity & Efficiency calculations ────────────────────────────────────────

  const Double_t range_lo = -20.0, range_hi = 20.0;
  Int_t bin_lo = event.getFracHistogram("Signal")->FindBin(range_lo);
  Int_t bin_hi = event.getFracHistogram("Signal")->FindBin(range_hi);

  // Całkowita wydajność
  Double_t pass_total = sig_pass->Integral(0, nbins + 1);
  Double_t total_total = sig_total->Integral(0, nbins + 1);
  Double_t total_eff = (total_total > 0.0) ? pass_total / total_total : 0.0;

  // Całkowita czystość
  Double_t sig_total_integral = event.getFracHistogram("Signal")->Integral(0, nbins + 1);
  Double_t mc_sum_total_integral = event.getFracHistogram("MC sum")->Integral(0, nbins + 1);
  Double_t total_purity = (mc_sum_total_integral > 0.0) ? sig_total_integral / mc_sum_total_integral : 0.0;

  // Czystość w [-20, 20]
  Double_t sig_range_integral = event.getFracHistogram("Signal")->Integral(bin_lo, bin_hi);
  Double_t mc_sum_range_integral = event.getFracHistogram("MC sum")->Integral(bin_lo, bin_hi);
  Double_t purity_range = (mc_sum_range_integral > 0.0) ? sig_range_integral / mc_sum_range_integral : 0.0;

  // Wydajność w [-20, 20]
  Int_t bin_lo_raw = sig_total->FindBin(range_lo);
  Int_t bin_hi_raw = sig_total->FindBin(range_hi);
  Double_t pass_range = sig_pass->Integral(bin_lo_raw, bin_hi_raw);
  Double_t total_range = sig_total->Integral(bin_lo_raw, bin_hi_raw);
  Double_t eff_range = (total_range > 0.0) ? pass_range / total_range : 0.0;

  // std::cout << "Branching Ratios:" << std::endl;
  // for (const auto &br : BRCF.BranchingRatios)
  // {
  //   std::cout << "  " << br.first << ": (" << br.second * 100 << " ± " << BRCF.BranchingRatios_err[br.first] * 100 << ") %" << std::endl;
  // }

  // std::cout << std::endl;
  // std::cout << "Old Branching Ratios:" << std::endl;
  // for (const auto &br : BRCF.BranchingRatios_old)
  // {
  //   std::cout << "  " << br.first << ": " << br.second * 100 << " %" << std::endl;
  // }
  // std::cout << std::endl;

  // std::cout << "BR Correction factors: " << std::endl;
  // for (const auto &factor : BRCF.BRcorrectionFactors)
  // {
  //   std::cout << "  " << factor.first << ": " << factor.second << " ± " << BRCF.BRcorrectionFactors_err[factor.first] << std::endl;
  // }
  std::cout << std::endl;

  std::cout << "=== Purity & Efficiency ===" << std::endl;
  std::cout << "Total efficiency:              " << total_eff * 100.0 << " %" << std::endl;
  std::cout << "Efficiency in [-20, 20] tau_S: " << eff_range * 100.0 << " %" << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << "Total purity:                  " << total_purity * 100.0 << " %" << std::endl;
  std::cout << "Purity in [-20, 20] tau_S:     " << purity_range * 100.0 << " %" << std::endl;
  std::cout << "===========================" << std::endl;

  TCanvas *c1 = new TCanvas("c1", "", 790, 1200);

  c1->SetBottomMargin(0.5);
  c1->Draw();

  TPad *padup_c1 = new TPad("pad_up_c1", "", 0.0, 0.3, 1.0, 1.0);
  TPad *paddown_c1 = new TPad("pad_down_c1", "", 0.0, 0.0, 1.0, 0.3);
  paddown_c1->SetBottomMargin(0.3);
  paddown_c1->SetRightMargin(0.10);

  gStyle->SetOptStat(0);

  TGraphAsymmErrors *sig_eff = new TGraphAsymmErrors();

  sig_eff->Divide(sig_pass, sig_total, "cl=0.683 b(1,1) mode");

  Double_t
      *yArr = sig_eff->GetY(),
      *yErrLArr = sig_eff->GetEYlow(),
      *yErrHArr = sig_eff->GetEYhigh();

  Int_t n_bins = sig_eff->GetN();
  std::vector<Double_t>
      y(yArr, yArr + n_bins),
      eyl(yErrLArr, yErrLArr + n_bins),
      eyh(yErrHArr, yErrHArr + n_bins);

  Double_t
      weightedMean = Obj.WeightedAverageAsymmetric(y, eyl, eyh),
      weightedMeanErr = Obj.WAvgAsymmError(eyl, eyh);

  TLine
      *lineAvg = new TLine(xMinRangeDisplay, weightedMean, xMaxRangeDisplay, weightedMean),
      *lineDown = new TLine(xMinRangeDisplay, weightedMean - weightedMeanErr, xMaxRangeDisplay, weightedMean - weightedMeanErr),
      *lineUp = new TLine(xMinRangeDisplay, weightedMean + weightedMeanErr, xMaxRangeDisplay, weightedMean + weightedMeanErr);

  lineAvg->SetLineColor(kRed);
  lineDown->SetLineColor(kRed);
  lineUp->SetLineColor(kRed);
  lineAvg->SetLineStyle(2); // np. przerywana
  lineDown->SetLineStyle(2);
  lineUp->SetLineStyle(2);
  lineAvg->SetLineWidth(2);
  lineDown->SetLineWidth(2);
  lineUp->SetLineWidth(2);

  c1->cd();
  paddown_c1->Draw();
  paddown_c1->cd();
  sig_eff->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);
  sig_eff->GetYaxis()->SetRangeUser(0.0, 1.0);
  sig_eff->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
  sig_eff->GetYaxis()->SetTitle("#varepsilon(#Deltat)");

  sig_eff->GetYaxis()->SetTitleSize(0.1);
  sig_eff->GetYaxis()->SetTitleOffset(0.5);
  sig_eff->GetXaxis()->SetTitleSize(0.1);
  sig_eff->GetYaxis()->SetLabelSize(0.05);
  sig_eff->GetXaxis()->SetLabelSize(0.08);

  sig_eff->Draw("APE");
  lineAvg->Draw("SAME");
  lineDown->Draw("SAME");
  lineUp->Draw("SAME");

  event.getFracHistogram("Data")->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);
  event.getFracHistogram("MC sum")->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);
  event.getFracHistogram("Signal")->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);

  // Set marker style for Data
  event.getFracHistogram("Data")->SetMarkerStyle(20);
  event.getFracHistogram("Data")->SetMarkerSize(0.8);
  event.getFracHistogram("Data")->SetLineColor(kBlack);

  // Calculate chi2 between Data and MC sum
  Double_t chi2 = event.getFracHistogram("Data")->Chi2Test(event.getFracHistogram("MC sum"), "WW CHI2/NDF");
  std::cout << "Chi2/NDF between Data and MC sum: " << chi2 << " (NDF = " << event.getFracHistogram("Data")->GetNbinsX() - 1 << ")" << std::endl;

  // Get correlations between parameters
  for (const auto &p1 : param_index_map)
  {
    if (p1.second < 0)
      continue; // Skip disabled parameters
    for (const auto &p2 : param_index_map)
    {
      if (p2.second < 0)
        continue; // Skip disabled parameters
      Double_t corr = minimum->Correlation(p1.second, p2.second);
      std::cout << "Correlation between " << p1.first << " and " << p2.first << ": " << corr << std::endl;
    }
  }

  TRatioPlot *rp = new TRatioPlot(event.getFracHistogram("MC sum"), event.getFracHistogram("Data"), "diffsig");

  c1->cd();
  padup_c1->Draw();
  padup_c1->cd();

  event.getFracHistogram("MC sum")->SetDrawOption("HIST E2");

  rp->Draw();

  rp->SetSplitFraction(0.2);

  rp->GetLowerRefGraph()->SetMinimum(-5);
  rp->GetLowerRefGraph()->SetMaximum(5);
  rp->GetLowerRefGraph()->SetLineWidth(3);

  rp->SetLowBottomMargin(0.0);
  rp->SetLeftMargin(0.15);

  rp->GetLowerRefYaxis()->SetLabelSize(0.02);
  rp->GetLowerRefXaxis()->SetLabelSize(0.0);

  Double_t max_height = event.getFracHistogram("Data")->GetMaximum();

  rp->GetUpperRefYaxis()->SetRangeUser(0.0, 2 * max_height);
  rp->GetUpperRefYaxis()->SetTitle(Form("Counts / %.2f #tau_{S}", res_deltaT));

  rp->GetLowerRefYaxis()->SetTitleSize(0.03);
  rp->GetLowerRefYaxis()->SetTitle("Residuals");

  rp->GetUpperPad()->cd();

  // Draw Signal as separate line
  event.getFracHistogram("Signal")->Draw("HIST SAME");

  // Draw other background channels
  for (auto const &name : KLOE::channName)
  {
    if (name.second == "Data" || name.second == "MC sum" || name.second == "Signal")
      continue;

    event.getFracHistogram(name.second)->Draw("HIST SAME");
  }

  TLegend *legend_chann = new TLegend(0.58, 0.48, 0.88, 0.88);
  legend_chann->SetFillColor(kWhite);
  legend_chann->AddEntry(event.getFracHistogram("Data"), KLOE::channTitle.at("Data"), "pe");

  for (auto const &name : KLOE::channName)
  {
    if (name.second == "Data" || name.second == "pi+pi-pi+pi-")
      continue;

    legend_chann->AddEntry(event.getFracHistogram(name.second), KLOE::channTitle.at(name.second), "l");
  }

  legend_chann->Draw();

  c1->Print(Paths::cpfit_dir + Paths::img_dir + "split_fit_with_corr" + Paths::ext_img);

  // Residuals graph
  TCanvas *c2 = new TCanvas("c2", "", 790, 790);
  TH1 *residuals_hist = new TH1D("Residuals hist", "", 51, -5., 5.);

  event.resi_vals = rp->GetLowerRefGraph()->GetY();

  for (Int_t i = 0; i < event.getBinNumber(); i++)
  {
    residuals_hist->Fill(event.resi_vals[i]);
  }

  residuals_hist->GetYaxis()->SetRangeUser(0, 1.2 * residuals_hist->GetMaximum());
  residuals_hist->Fit("gaus");

  c2->cd();

  residuals_hist->SetXTitle("Residuals");
  residuals_hist->SetYTitle("Counts");
  residuals_hist->SetLineWidth(5);
  residuals_hist->Draw();
  residuals_hist->SetStats(1);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  residuals_hist->Draw();

  c2->Print(Paths::cpfit_dir + Paths::img_dir + "residuals_hist" + Paths::ext_img);

  TCanvas *c3 = new TCanvas("c3", "", 790, 790);
  hEffWeightVsRelErr->Draw("COLZ");
  c3->Print(Paths::cpfit_dir + Paths::img_dir + "eff_weight_vs_relerr" + Paths::ext_img);

  TCanvas *c4 = new TCanvas("c4", "", 790, 790);
  hEffErrorProfile->GetXaxis()->SetRangeUser(xMinRangeDisplay, xMaxRangeDisplay);
  hEffErrorProfile->Draw("COLZ");
  c4->Print(Paths::cpfit_dir + Paths::img_dir + "eff_error_profile" + Paths::ext_img);

  Utils::properties["variables"]["CPFit"]["result"]["value"]["Re"] = par[0];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Im"] = par[1];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Signal"] = par[2];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["FarLeft"] = par[3];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["CloseLeft"] = par[4];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["CloseRight"] = par[5];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Regeneration"]["FarRight"] = par[6];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Omegapi0"] = par[7];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Threepi0"] = par[8];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Semileptonic"] = par[9];
  Utils::properties["variables"]["CPFit"]["result"]["value"]["Norm"]["Other"] = par[10];

  Utils::properties["variables"]["CPFit"]["result"]["error"]["PhysicsConstants::Re"] = parErr[0];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Im"] = parErr[1];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Signal"] = parErr[2];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["FarLeft"] = parErr[3];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["CloseLeft"] = parErr[4];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["CloseRight"] = parErr[5];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Regeneration"]["FarRight"] = parErr[6];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Omegapi0"] = parErr[7];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Threepi0"] = parErr[8];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Semileptonic"] = parErr[9];
  Utils::properties["variables"]["CPFit"]["result"]["error"]["Norm"]["Other"] = parErr[10];

  Utils::properties["variables"]["CPFit"]["result"]["chi2"] = event.getFracHistogram("Data")->Chi2Test(event.getFracHistogram("MC sum"), "UW CHI2");
  Utils::properties["variables"]["CPFit"]["result"]["normChi2"] = event.getFracHistogram("Data")->Chi2Test(event.getFracHistogram("MC sum"), "UW CHI2/NDF");

  // delete residuals_hist;
  // delete c1;
  // delete c2;

  // delete rp;

  Utils::properties["lastScript"] = "Final CP Parameters normalization";
  Utils::properties["lastUpdate"] = Obj.getCurrentTimestamp();

  std::ofstream outfile(Paths::propName);
  outfile << Utils::properties.dump(4);
  outfile.close();

  return 0;
}