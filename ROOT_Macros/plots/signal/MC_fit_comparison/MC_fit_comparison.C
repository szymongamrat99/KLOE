#define MC_fit_comparison_cxx
// The class definition in MC_fit_comparison.h has been generated automatically
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
// root> T->Process("MC_fit_comparison.C")
// root> T->Process("MC_fit_comparison.C","some options")
// root> T->Process("MC_fit_comparison.C+")
//

#include "MC_fit_comparison.h"
#include <kloe_class.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TripleGaussFitter.h>
#include <DoublGaussFitter.h>
#include <triple_gaus.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TRegexp.h>
#include <interf_function.h>
#include <charged_mom.h>
#include <StatisticalCutter.h>

#include <TFitResult.h>
#include <TFitResultPtr.h>

namespace KH = KLOE::Histograms;

Int_t signal_num = 0, signal_tot = 0, signal_tot_err = 0, tot_events = 0, bkg_tot = 0;

std::map<TString, TCanvas *>
    canvas;

std::map<TString, TCanvas *>
    canvasProfiles;

std::map<TString, TProfile *>
    profiles;

std::map<TString, TCanvas *>
    canvas2D;

std::map<TString, TH1 *>
    histsReconstructed,
    histsFittedSignal;

std::map<TString, TH2 *>
    hists2DReconstructed,
    hists2DFittedSignal;

TCanvas *canvaEff;

TH1 *deltaTSignalTot;

KLOE::TripleGaussFitter *fitter;
KLOE::DoublGaussFitter *doubleFitter;

KLOE::pm00 Obj;
KLOE::ChargedVtxRec<> chargedVtxRec;

// Wczytaj konfiguracje histogramów
auto histogramConfigs1D = KLOE::Histograms::LoadHistogramConfigs1D(Paths::histogramConfig1DPath);
auto histogramConfigs2D = KLOE::Histograms::LoadHistogramConfigs2D(Paths::histogramConfig2DPath);
///////////////////////////////////////////////////////////////////////////////

Int_t nbins = 121;

Float_t ndofSignal = 10., ndofOmega = 8.;

Float_t Chi2SignalReduced;

Float_t deltaPhi,
    deltaPhiMC,
    deltaPhiFit,
    deltaTheta,
    RtCh,
    RtNeu,
    RtChRec,
    RtNeuRec,
    ZChRec,
    ZNeuRec;

Float_t T0Omega = 0;

Float_t a = 1,
        b = 625.091,
        Breal = 14.1421;

// Definitions of cuts

StatisticalCutter *cutter;

std::map<std::string, std::function<Float_t()>> cutValues, cutLimits;
std::map<std::string, Float_t> centralValues;

std::vector<size_t> cutIndices;

void MC_fit_comparison::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  // Creation of output folder with cut names
  /////////////////////////////////////////
  TString option = GetOption();

  if (option.IsNull())
    option = "NO_CUTS";

  InitializeCutSelector(option);
  cutIndices = GetCutIndicesForOption(option);

  TString folderName = SanitizeFolderName(option);
  folderPath = "img/" + folderName;
  FolderManagement(folderPath);

  KLOE::setGlobalStyle();

  fitter = new KLOE::TripleGaussFitter();
  doubleFitter = new KLOE::DoublGaussFitter();

  canvaEff = new TCanvas("Efficiency", "Efficiency", 800, 800);
  deltaTSignalTot = new TH1D("EfficiencyHistTot", "Efficiency Histogram Total; #Deltat [#tau_{S}]; Efficiency [-]", nbins, -30, 30);

  // Create canvases
  for (const auto &histName : histogramConfigs1D)
  {
    canvas[histName.first] = new TCanvas(Form("c_%s", histName.first.Data()),
                                         Form("Canvas for %s", histName.first.Data()), 750, 750);
  }

  // Create canvases
  for (const auto &histName : histogramConfigs2D)
  {
    canvas2D[histName.first] = new TCanvas(Form("c_%s", histName.first.Data()), Form("Canvas for %s", histName.first.Data()), 750, 750);
  }

  for (const auto &config : histogramConfigs2D)
  {
    canvasProfiles[config.first] = new TCanvas(Form("cProfiles_%s", config.first.Data()),
                                               Form("Canvas profiles for %s", config.first.Data()), 750, 750);
  }

  // Create histograms
  for (const auto &histName : histogramConfigs1D)
  {

    TString nameRec = Form("h_rec_%s", histName.first.Data());
    TString nameFit = Form("h_fit_%s", histName.first.Data());

    histsReconstructed[histName.first] = KH::CreateHist1D(nameRec, histName.second);
    histsFittedSignal[histName.first] = KH::CreateHist1D(nameFit, histName.second);
  }

  // Create histograms 2D
  for (const auto &histName : histogramConfigs2D)
  {
    TString nameRec = Form("h_rec2D_%s", histName.first.Data());
    TString nameFit = Form("h_fit2D_%s", histName.first.Data());

    hists2DReconstructed[histName.first] = KH::CreateHist2D(nameRec, histName.second);
    hists2DFittedSignal[histName.first] = KH::CreateHist2D(nameFit, histName.second);
  }
}

void MC_fit_comparison::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

Int_t numberOfAllGood = 0,
      numberOfAtLeastOneBad = 0;

Bool_t MC_fit_comparison::Process(Long64_t entry)
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

  TString option = GetOption();

  fReader.SetLocalEntry(entry);

  Int_t mctruth_int = *mctruth;

  if (mctruth_int == 0 && *mcflag == 1)
  {
    mctruth_int = 1;
  }

  Float_t weight = 1.0;

  if ((mctruth_int == 1 || mctruth_int == 0) && *mcflag == 1)
    weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

  Float_t vKchFit = PhysicsConstants::cVel * KchboostFit[4] / KchboostFit[3],
          pathKchFit = sqrt(pow(KchboostFit[6] - 0, 2) +
                            pow(KchboostFit[7] - 0, 2) +
                            pow(KchboostFit[8] - 0, 2)),
          tKchFit = KchboostFit[9] / (0.0895),
          vKneFit = PhysicsConstants::cVel * KnereclorFit[4] / KnereclorFit[3],
          pathKneFit = sqrt(pow(KnereclorFit[6] - 0, 2) +
                            pow(KnereclorFit[7] - 0, 2) +
                            pow(KnereclorFit[8] - 0, 2)),
          tKneFit = KnerecFit[9] / (0.0895),
          vKneMC = PhysicsConstants::cVel * Knemc[4] / Knemc[3],
          vKchMC = PhysicsConstants::cVel * Kchmc[4] / Kchmc[3],
          vKne = PhysicsConstants::cVel * Knerec[4] / Knerec[3],
          pathKne = sqrt(pow(Knerec[6], 2) +
                         pow(Knerec[7], 2) +
                         pow(Knerec[8], 2)),
          pathKch = sqrt(pow(Kchrec[6], 2) +
                         pow(Kchrec[7], 2) +
                         pow(Kchrec[8], 2)),
          tKne = pathKne / (vKne * 0.0895);

  Float_t combinedMassPi0Fit = sqrt(pow(pi01Fit[5] - PhysicsConstants::mPi0, 2) +
                                    pow(pi02Fit[5] - PhysicsConstants::mPi0, 2)),
          combinedMassPi0 = sqrt(pow(pi01[5] - PhysicsConstants::mPi0, 2) +
                                 pow(pi02[5] - PhysicsConstants::mPi0, 2));

  Float_t photon1path = sqrt(pow(photonFit1[4] - KnerecFit[6], 2) +
                             pow(photonFit1[5] - KnerecFit[7], 2) +
                             pow(photonFit1[6] - KnerecFit[8], 2)),
          photon2path = sqrt(pow(photonFit2[4] - KnerecFit[6], 2) +
                             pow(photonFit2[5] - KnerecFit[7], 2) +
                             pow(photonFit2[6] - KnerecFit[8], 2)),
          photon3path = sqrt(pow(photonFit3[4] - KnerecFit[6], 2) +
                             pow(photonFit3[5] - KnerecFit[7], 2) +
                             pow(photonFit3[6] - KnerecFit[8], 2)),
          photon4path = sqrt(pow(photonFit4[4] - KnerecFit[6], 2) +
                             pow(photonFit4[5] - KnerecFit[7], 2) +
                             pow(photonFit4[6] - KnerecFit[8], 2));

  Float_t trc1Fit = photonFit1[7] - photon1path / PhysicsConstants::cVel - tKneFit * 0.0895,
          trc2Fit = photonFit2[7] - photon2path / PhysicsConstants::cVel - tKneFit * 0.0895,
          trc3Fit = photonFit3[7] - photon3path / PhysicsConstants::cVel - tKneFit * 0.0895,
          trc4Fit = photonFit4[7] - photon4path / PhysicsConstants::cVel - tKneFit * 0.0895,
          TrcSumFit = trc1Fit + trc2Fit + trc3Fit + trc4Fit;

  std::vector<Float_t> KchboostMom = {Kchboost[0], Kchboost[1], Kchboost[2], Kchboost[3]},
                       KchrecMom = {Kchrec[0], Kchrec[1], Kchrec[2], Kchrec[3]},
                       KchboostFitMom = {KchrecFit[0], KchrecFit[1], KchrecFit[2], KchrecFit[3]},
                       KnerecMom = {Knerec[0], Knerec[1], Knerec[2], Knerec[3]},
                       KnemcMom = {Knemc[0], Knemc[1], Knemc[2], Knemc[3]},
                       KnereclorMom = {*Bpx - Kchboost[0], *Bpy - Kchboost[1], *Bpz - Kchboost[2], *Broots - Kchboost[3]},
                       KnerecFitMom = {KnerecFit[0], KnerecFit[1], KnerecFit[2], KnerecFit[3]},
                       KchPos = {Kchboost[6], Kchboost[7], Kchboost[8]},
                       KnePos = {Knerec[6], Knerec[7], Knerec[8]},
                       KnePosMC = {Knemc[6], Knemc[7], Knemc[8]},
                       KchPosFit = {KchboostFit[6], KchboostFit[7], KchboostFit[8]},
                       KnePosFit = {KnerecFit[6], KnerecFit[7], KnerecFit[8]},
                       ipPos = {ip[0], ip[1], ip[2]},
                       ipPosFit = {ipFit[0], ipFit[1], ipFit[2]},
                       ipBhaPos = {*Bx, *By, *Bz},
                       ipCenter = {0., 0., 0.};

  KLOE::KaonProperTimes timesBoostLor = Obj.CalculateKaonProperTimes(KchboostMom,
                                                                     KchPos,
                                                                     KnereclorMom,
                                                                     KnePos,
                                                                     ipPos);

  KLOE::KaonProperTimes timesRecRec = Obj.CalculateKaonProperTimes(KchrecMom,
                                                                   KchPos,
                                                                   KnerecMom,
                                                                   KnePos,
                                                                   ipPos);

  KLOE::KaonProperTimes timesSignalFit = Obj.CalculateKaonProperTimes(KchboostFitMom,
                                                                      KchPosFit,
                                                                      KnerecFitMom,
                                                                      KnePosFit,
                                                                      ipPosFit);

  Float_t Omega1MassTmp = sqrt(pow(trk1[3] + trk2[3] + pi0Omega1[3], 2) -
                               pow(trk1[0] + trk2[0] + pi0Omega1[0], 2) -
                               pow(trk1[1] + trk2[1] + pi0Omega1[1], 2) -
                               pow(trk1[2] + trk2[2] + pi0Omega1[2], 2)),
          Omega2MassTmp = sqrt(pow(trk1[3] + trk2[3] + pi0Omega2[3], 2) -
                               pow(trk1[0] + trk2[0] + pi0Omega2[0], 2) -
                               pow(trk1[1] + trk2[1] + pi0Omega2[1], 2) -
                               pow(trk1[2] + trk2[2] + pi0Omega2[2], 2)),
          Omega1ErrTmp = abs(Omega1MassTmp - PhysicsConstants::mOmega),
          Omega2ErrTmp = abs(Omega2MassTmp - PhysicsConstants::mOmega);

  Float_t radius00 = sqrt(pow(KnerecFit[6] - ipFit[0], 2) +
                          pow(KnerecFit[7] - ipFit[1], 2)),
          radiuspm = sqrt(pow(KchrecFit[6] - ipFit[0], 2) +
                          pow(KchrecFit[7] - ipFit[1], 2)),
          zdist00 = abs(KnerecFit[8] - ipFit[2]),
          zdistpm = abs(KchrecFit[8] - ipFit[2]),
          path00 = sqrt(pow(Knerec[6] - *Bx, 2) +
                        pow(Knerec[7] - *By, 2) +
                        pow(Knerec[8] - *Bz, 2)),
          pathpm = sqrt(pow(Kchrec[6] - *Bx, 2) +
                        pow(Kchrec[7] - *By, 2) +
                        pow(Kchrec[8] - *Bz, 2)),
          path00Center = sqrt(pow(Knerec[6], 2) +
                              pow(Knerec[7], 2) +
                              pow(Knerec[8], 2)),
          pathpmCenter = sqrt(pow(Kchrec[6], 2) +
                              pow(Kchrec[7], 2) +
                              pow(Kchrec[8], 2)),
          path00FitCenter = sqrt(pow(KnerecFit[6], 2) +
                                 pow(KnerecFit[7], 2) +
                                 pow(KnerecFit[8], 2)),
          pathpmFitCenter = sqrt(pow(KchrecFit[6], 2) +
                                 pow(KchrecFit[7], 2) +
                                 pow(KchrecFit[8], 2)),
          path00MC = sqrt(pow(Knemc[6] - ipmc[0], 2) +
                          pow(Knemc[7] - ipmc[1], 2) +
                          pow(Knemc[8] - ipmc[2], 2)),
          pathpmMC = sqrt(pow(Kchmc[6] - ipmc[0], 2) +
                          pow(Kchmc[7] - ipmc[1], 2) +
                          pow(Kchmc[8] - ipmc[2], 2)),
          radius00MC = sqrt(pow(Knemc[6] - ipmc[0], 2) +
                            pow(Knemc[7] - ipmc[1], 2)),
          radiuspmMC = sqrt(pow(Kchmc[6] - ipmc[0], 2) +
                            pow(Kchmc[7] - ipmc[1], 2)),
          path00MCCenter = sqrt(pow(Knemc[6], 2) +
                                pow(Knemc[7], 2) +
                                pow(Knemc[8], 2)),
          pathpmMCCenter = sqrt(pow(Kchmc[6], 2) +
                                pow(Kchmc[7], 2) +
                                pow(Kchmc[8], 2)),
          radius00MCCenter = sqrt(pow(Knemc[6], 2) +
                                  pow(Knemc[7], 2)),
          radiuspmMCCenter = sqrt(pow(Kchmc[6], 2) +
                                  pow(Kchmc[7], 2)),
          radiusLimit = 1,
          zdistLimit = 0.6;

  Bool_t isInsideFiducialVolume = (radius00 < radiusLimit) && (zdist00 < zdistLimit) &&
                                  (radiuspm < radiusLimit) && (zdistpm < zdistLimit);

  T0Omega = pi0OmegaFit1[3] - pi0OmegaFit1[5];

  Float_t deltaTfit = timesSignalFit.deltaTimeCM,
          deltaT = timesRecRec.deltaTimeCM,
          deltaTMC = *KaonChTimeCMMC - *KaonNeTimeCMMC;

  Chi2SignalReduced = *Chi2SignalKinFit / ndofSignal;

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

  deltaPhi = abs(Phi2Rec - Phi1Rec);
  deltaTheta = abs(Theta2Fit - Theta1Fit);
  deltaPhiFit = abs(Phi2Fit - Phi1Fit);

  // Analiza Simony cięcie na 3 sigma mas
  Bool_t condMassKch = abs(Kchrec[5] - 497.605) < 3 * 0.879,
         condMassKne = abs(*minv4gam - 488.411) < 3 * 41.293,
         condMassPi01 = abs(pi01Fit[5] - 134.840) < 3 * 3.479,
         condMassPi02 = abs(pi02Fit[5] - 134.867) < 3 * 3.331;
  ///////////////////////////////////////////////////////////////////////////////

  Bool_t simonaCuts = abs(deltaPhiFit - 3.110) > 2 * 0.135 && *Chi2SignalKinFit < 30.,
         simonaKinCuts = condMassKch && condMassKne && condMassPi01 && condMassPi02 && simonaCuts;

  if (mctruth_int == 0 || mctruth_int == -1 || mctruth_int == 1)
    signal_tot_err++;

  if ((mctruth_int == 0 || mctruth_int == 1))
    signal_tot++;

  Double_t limitRadiusChMC = 25.,
           limitRadiusNeMC = 25.;

  RtNeu = sqrt(pow(Knemc[6] - ipmc[0], 2) +
               pow(Knemc[7] - ipmc[1], 2));
  RtCh = sqrt(pow(Kchmc[6] - ipmc[0], 2) +
              pow(Kchmc[7] - ipmc[1], 2));

  RtChRec = sqrt(pow(KchrecFit[6] - ipFit[0], 2) +
                 pow(KchrecFit[7] - ipFit[1], 2));
  RtNeuRec = sqrt(pow(KnerecFit[6] - ipFit[0], 2) +
                  pow(KnerecFit[7] - ipFit[1], 2));

  ZChRec = abs(KchrecFit[8] - ipFit[2]);
  ZNeuRec = abs(KnerecFit[8] - ipFit[2]);

  Bool_t lastSimonaCut = (isInsideFiducialVolume && !(cutter->PassCut(12) && cutter->PassCut(13) && cutter->PassCut(14) && cutter->PassCut(15))) || (!isInsideFiducialVolume);

  // Sprawdź czy zdarzenie przechodzi wszystkie aktywne cięcia (z obsługą grup background rejection)
  bool passesCuts = cutter->PassCuts(cutIndices) && lastSimonaCut;

  // Update statistics for all events
  cutter->UpdateStats(mctruth_int);

  if (!passesCuts)
    return kTRUE;

  if (mctruth_int == 1)
  {
    Int_t mctruth_tmp = mctruth_int;

    if (mctruth_int == 0)
      mctruth_tmp = 1;

    if (*goodClustersTriKinFitSize >= 3)
      numberOfAtLeastOneBad++;
    if (*goodClustersTriKinFitSize >= 4)
      numberOfAllGood++;

    signal_num++;

    histsReconstructed["TransvRadius"]->Fill(pathpmMCCenter, weight);
    histsFittedSignal["TransvRadius"]->Fill(path00MCCenter, weight);

    // Fill histograms for reconstructed variables
    histsReconstructed["px_Kch"]->Fill(Kchrec[0] - Kchmc[0], weight);
    histsReconstructed["py_Kch"]->Fill(Kchrec[1] - Kchmc[1], weight);
    histsReconstructed["pz_Kch"]->Fill(Kchrec[2] - Kchmc[2], weight);
    histsReconstructed["Energy_Kch"]->Fill(Kchrec[3] - Kchmc[3], weight);
    histsReconstructed["mass_Kch"]->Fill(Kchrec[5], weight);

    histsReconstructed["px_Kne"]->Fill(Knerec[0] - Knemc[0], weight);
    histsReconstructed["py_Kne"]->Fill(Knerec[1] - Knemc[1], weight);
    histsReconstructed["pz_Kne"]->Fill(Knerec[2] - Knemc[2], weight);
    histsReconstructed["Energy_Kne"]->Fill(Knerec[3] - Knemc[3], weight);
    histsReconstructed["mass_Kne"]->Fill(*minv4gam, weight);

    histsReconstructed["mass_pi01"]->Fill(pi01[5], weight);
    histsReconstructed["mass_pi02"]->Fill(pi02[5], weight);

    histsReconstructed["px_Phi"]->Fill(*Bpx - Kchmc[0] - Knemc[0], weight);
    histsReconstructed["py_Phi"]->Fill(*Bpy - Kchmc[1] - Knemc[1], weight);
    histsReconstructed["pz_Phi"]->Fill(*Bpz - Kchmc[2] - Knemc[2], weight);
    histsReconstructed["Energy_Phi"]->Fill(*Broots - Kchmc[3] - Knemc[3], weight);

    if (abs(Omega1ErrTmp) > abs(Omega2ErrTmp))
      histsReconstructed["mass_omega"]->Fill(Omega2MassTmp, weight);
    else
      histsReconstructed["mass_omega"]->Fill(Omega1MassTmp, weight);

    if (abs(Omega1ErrTmp) > abs(Omega2ErrTmp))
      histsReconstructed["mass_omega_rec"]->Fill(Omega2MassTmp, weight);
    else
      histsReconstructed["mass_omega_rec"]->Fill(Omega1MassTmp, weight);

    histsReconstructed["vtxNeu_x"]->Fill(Knerec[6] - Knemc[6], weight);
    histsReconstructed["vtxNeu_y"]->Fill(Knerec[7] - Knemc[7], weight);
    histsReconstructed["vtxNeu_z"]->Fill(Knerec[8] - Knemc[8], weight);

    histsReconstructed["vtxCh_x"]->Fill(Kchrec[6] - Kchmc[6], weight);
    histsReconstructed["vtxCh_y"]->Fill(Kchrec[7] - Kchmc[7], weight);
    histsReconstructed["vtxCh_z"]->Fill(Kchrec[8] - Kchmc[8], weight);

    histsReconstructed["phi_vtx_x"]->Fill(ip[0] - ipmc[0], weight);
    histsReconstructed["phi_vtx_y"]->Fill(ip[1] - ipmc[1], weight);
    histsReconstructed["phi_vtx_z"]->Fill(ip[2] - ipmc[2], weight);

    histsReconstructed["vtxNeu_x_Fit"]->Fill(Knerec[6] - Knemc[6], weight);
    histsReconstructed["vtxNeu_y_Fit"]->Fill(Knerec[7] - Knemc[7], weight);
    histsReconstructed["vtxNeu_z_Fit"]->Fill(Knerec[8] - Knemc[8], weight);

    // histsReconstructed["vKne"]->Fill(vKne - vKneMC, weight);

    histsReconstructed["time_neutral_MC"]->Fill(*TrcSum, weight);

    histsReconstructed["delta_t_boostlor_res"]->Fill(timesBoostLor.deltaTimeCM - deltaTMC, weight);

    histsReconstructed["delta_t_rec_res"]->Fill(timesRecRec.deltaTimeCM - deltaTMC, weight);

    histsReconstructed["combined_mass_pi0"]->Fill(combinedMassPi0, weight);

    histsReconstructed["T0Omega"]->Fill(T0Omega, weight);

    // Decide which reconstructed track corresponds to which MC particle

    Double_t errorPart10, errorPart20, errorPart11, errorPart21, error1021, error1120;

    errorPart10 = pow(abs(*CurvSmeared1) - abs(CurvMC[0]), 2) +
                  pow(abs(*PhivSmeared1) - abs(PhivMC[0]), 2) +
                  pow(abs(*CotvSmeared1) - abs(CotvMC[0]), 2);

    errorPart20 = pow(abs(*CurvSmeared2) - abs(CurvMC[0]), 2) +
                  pow(abs(*PhivSmeared2) - abs(PhivMC[0]), 2) +
                  pow(abs(*CotvSmeared2) - abs(CotvMC[0]), 2);

    errorPart11 = pow(abs(*CurvSmeared1) - abs(CurvMC[1]), 2) +
                  pow(abs(*PhivSmeared1) - abs(PhivMC[1]), 2) +
                  pow(abs(*CotvSmeared1) - abs(CotvMC[1]), 2);

    errorPart21 = pow(abs(*CurvSmeared2) - abs(CurvMC[1]), 2) +
                  pow(abs(*PhivSmeared2) - abs(PhivMC[1]), 2) +
                  pow(abs(*CotvSmeared2) - abs(CotvMC[1]), 2);

    error1021 = sqrt(errorPart10 + errorPart21);
    error1120 = sqrt(errorPart11 + errorPart20);

    if (error1021 < error1120)
    {
      histsReconstructed["curv1"]->Fill(abs(*CurvSmeared1) - abs(CurvMC[0]), weight);
      histsReconstructed["phiv1"]->Fill(abs(*PhivSmeared1) - abs(PhivMC[0]), weight);
      histsReconstructed["cotv1"]->Fill(abs(*CotvSmeared1) - abs(CotvMC[0]), weight);

      histsReconstructed["curv2"]->Fill(abs(*CurvSmeared2) - abs(CurvMC[1]), weight);
      histsReconstructed["phiv2"]->Fill(abs(*PhivSmeared2) - abs(PhivMC[1]), weight);
      histsReconstructed["cotv2"]->Fill(abs(*CotvSmeared2) - abs(CotvMC[1]), weight);
    }
    else
    {
      histsReconstructed["curv1"]->Fill(abs(*CurvSmeared1) - abs(CurvMC[1]), weight);
      histsReconstructed["phiv1"]->Fill(abs(*PhivSmeared1) - abs(PhivMC[1]), weight);
      histsReconstructed["cotv1"]->Fill(abs(*CotvSmeared1) - abs(CotvMC[1]), weight);

      histsReconstructed["curv2"]->Fill(abs(*CurvSmeared2) - abs(CurvMC[0]), weight);
      histsReconstructed["phiv2"]->Fill(abs(*PhivSmeared2) - abs(PhivMC[0]), weight);
      histsReconstructed["cotv2"]->Fill(abs(*CotvSmeared2) - abs(CotvMC[0]), weight);
    }

    errorPart10 = pow(trk1[0] - trk1MC[0], 2) +
                  pow(trk1[1] - trk1MC[1], 2) +
                  pow(trk1[2] - trk1MC[2], 2);

    errorPart20 = pow(trk2[0] - trk1MC[0], 2) +
                  pow(trk2[1] - trk1MC[1], 2) +
                  pow(trk2[2] - trk1MC[2], 2);

    errorPart11 = pow(trk1[0] - trk2MC[0], 2) +
                  pow(trk1[1] - trk2MC[1], 2) +
                  pow(trk1[2] - trk2MC[2], 2);

    errorPart21 = pow(trk2[0] - trk2MC[0], 2) +
                  pow(trk2[1] - trk2MC[1], 2) +
                  pow(trk2[2] - trk2MC[2], 2);

    error1021 = sqrt(errorPart10 + errorPart21);
    error1120 = sqrt(errorPart11 + errorPart20);

    if (error1021 < error1120)
    {
      histsReconstructed["trk1x"]->Fill(trk1[0] - trk1MC[0], weight);
      histsReconstructed["trk1y"]->Fill(trk1[1] - trk1MC[1], weight);
      histsReconstructed["trk1z"]->Fill(trk1[2] - trk1MC[2], weight);

      histsReconstructed["trk2x"]->Fill(trk2[0] - trk2MC[0], weight);
      histsReconstructed["trk2y"]->Fill(trk2[1] - trk2MC[1], weight);
      histsReconstructed["trk2z"]->Fill(trk2[2] - trk2MC[2], weight);
    }
    else
    {
      histsReconstructed["trk1x"]->Fill(trk1[0] - trk2MC[0], weight);
      histsReconstructed["trk1y"]->Fill(trk1[1] - trk2MC[1], weight);
      histsReconstructed["trk1z"]->Fill(trk1[2] - trk2MC[2], weight);

      histsReconstructed["trk2x"]->Fill(trk2[0] - trk1MC[0], weight);
      histsReconstructed["trk2y"]->Fill(trk2[1] - trk1MC[1], weight);
      histsReconstructed["trk2z"]->Fill(trk2[2] - trk1MC[2], weight);
    }

    // Fitted signal variables
    histsFittedSignal["px_Kch"]->Fill(KchrecFit[0] - Kchmc[0], weight);
    histsFittedSignal["py_Kch"]->Fill(KchrecFit[1] - Kchmc[1], weight);
    histsFittedSignal["pz_Kch"]->Fill(KchrecFit[2] - Kchmc[2], weight);
    histsFittedSignal["Energy_Kch"]->Fill(KchrecFit[3] - Kchmc[3], weight);
    histsFittedSignal["mass_Kch"]->Fill(KchrecFit[5], weight);

    histsFittedSignal["px_Kne"]->Fill(KnerecFit[0] - Knemc[0], weight);
    histsFittedSignal["py_Kne"]->Fill(KnerecFit[1] - Knemc[1], weight);
    histsFittedSignal["pz_Kne"]->Fill(KnerecFit[2] - Knemc[2], weight);
    histsFittedSignal["Energy_Kne"]->Fill(KnerecFit[3] - Knemc[3], weight);
    histsFittedSignal["mass_Kne"]->Fill(KnerecFit[5], weight);

    histsFittedSignal["mass_pi01"]->Fill(pi01Fit[5], weight);
    histsFittedSignal["mass_pi02"]->Fill(pi02Fit[5], weight);

    histsFittedSignal["mass_omega"]->Fill(omegaFit[5], weight);

    histsReconstructed["chi2_signalKinFit"]->Fill(*Chi2SignalKinFit / 10., weight);
    histsReconstructed["chi2_trilaterationKinFit"]->Fill(*Chi2TriKinFit / 5., weight);
    histsReconstructed["prob_signal"]->Fill(TMath::Prob(*Chi2SignalKinFit, 10), weight);

    histsFittedSignal["vtxNeu_x_Fit"]->Fill(KnerecFit[6] - Knemc[6], weight);
    histsFittedSignal["vtxNeu_y_Fit"]->Fill(KnerecFit[7] - Knemc[7], weight);
    histsFittedSignal["vtxNeu_z_Fit"]->Fill(KnerecFit[8] - Knemc[8], weight);

    histsFittedSignal["px_Phi"]->Fill(ParamSignalFit[32] - Kchmc[0] - Knemc[0], weight);
    histsFittedSignal["py_Phi"]->Fill(ParamSignalFit[33] - Kchmc[1] - Knemc[1], weight);
    histsFittedSignal["pz_Phi"]->Fill(ParamSignalFit[34] - Kchmc[2] - Knemc[2], weight);
    Double_t energyPhiFit = sqrt(pow(ParamSignalFit[32], 2) +
                                 pow(ParamSignalFit[33], 2) +
                                 pow(ParamSignalFit[34], 2) +
                                 pow(PhysicsConstants::mPhi, 2));
    histsFittedSignal["Energy_Phi"]->Fill(energyPhiFit - Kchmc[3] - Knemc[3], weight);

    histsFittedSignal["phi_vtx_x"]->Fill(ipFit[0] - ipmc[0], weight);
    histsFittedSignal["phi_vtx_y"]->Fill(ipFit[1] - ipmc[1], weight);
    histsFittedSignal["phi_vtx_z"]->Fill(ipFit[2] - ipmc[2], weight);

    // histsFittedSignal["vKne"]->Fill(vKneFit - vKneMC, weight);

    histsFittedSignal["delta_t_boostlor_res"]->Fill(deltaTfit - deltaTMC, weight);

    histsFittedSignal["delta_t_rec_res"]->Fill(deltaTfit - deltaTMC, weight);

    histsFittedSignal["combined_mass_pi0"]->Fill(combinedMassPi0Fit, weight);

    for (Int_t i = 0; i < 39; i++)
    {
      histsFittedSignal["pull" + std::to_string(i + 1)]->Fill(pullsSignalFit[i], weight);
    }

    histsFittedSignal["time_neutral_MC"]->Fill(TrcSumFit, weight);

    for (Int_t i = 0; i < 38; i++)
    {
      Double_t chi2inp = pow(ParamSignalFit[i] - ParamSignal[i], 2) / pow(ErrorsSignal[i], 2);
      histsReconstructed["Chi2Comp" + std::to_string(i + 1)]->Fill(chi2inp, weight);
    }

    hists2DFittedSignal["px_ch_res_vs_px_neu_res"]->Fill(ParamSignalFit[32] - Kchmc[0] - Knemc[0], KchboostFit[0] - Kchmc[0], weight);
    hists2DFittedSignal["py_ch_res_vs_py_neu_res"]->Fill(ParamSignalFit[33] - Kchmc[1] - Knemc[1], KchboostFit[1] - Kchmc[1], weight);
    hists2DFittedSignal["pz_ch_res_vs_pz_neu_res"]->Fill(ParamSignalFit[34] - Kchmc[2] - Knemc[2], KchboostFit[2] - Kchmc[2], weight);
    hists2DFittedSignal["E_ch_res_vs_E_neu_res"]->Fill(ParamSignalFit[35] - Kchmc[3] - Knemc[3], KchboostFit[3] - Kchmc[3], weight);

    hists2DFittedSignal["T0_omega_vs_mass_omega"]->Fill(T0Omega, omegaFit[5], weight);

    hists2DFittedSignal["delta_t_mc_vs_delta_t_res"]->Fill(deltaTMC, deltaTfit - deltaTMC, weight);
    hists2DFittedSignal["delta_t_vs_delta_t_mc"]->Fill(deltaTMC, deltaT - deltaTMC, weight);

    hists2DFittedSignal["delta_t_fit_vs_delta_t_mc"]->Fill(deltaTMC, deltaTfit, weight);

    hists2DFittedSignal["t_ch_fit_vs_t_ch_mc"]->Fill(*KaonChTimeCMMC, *KaonChTimeCMSignalFit - *KaonChTimeCMMC, weight);
    hists2DFittedSignal["t_neu_fit_vs_t_neu_mc"]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMSignalFit - *KaonNeTimeCMMC, weight);

    hists2DFittedSignal["t_ch_rec_vs_t_ch_mc"]->Fill(*KaonChTimeCMMC, timesRecRec.kaon1TimeCM - *KaonChTimeCMMC, weight);
    hists2DFittedSignal["t_neu_rec_vs_t_neu_mc"]->Fill(*KaonNeTimeCMMC, timesRecRec.kaon2TimeCM - *KaonNeTimeCMMC, weight);
  }

  return kTRUE;
}

void MC_fit_comparison::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void MC_fit_comparison::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Bool_t recHasEntries = kFALSE,
         fitHasEntries = kFALSE;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for (const auto &config : histogramConfigs1D)
  {
    if (histsFittedSignal[config.first]->GetEntries() > 0)
    {
      histsReconstructed[config.first]->Scale(histsReconstructed[config.first]->GetEntries() / histsReconstructed[config.first]->Integral(0, histsReconstructed[config.first]->GetNbinsX() + 1));
      histsFittedSignal[config.first]->Scale(histsFittedSignal[config.first]->GetEntries() / histsFittedSignal[config.first]->Integral(0, histsFittedSignal[config.first]->GetNbinsX() + 1));
    }
  }

  for (const auto &histName : histogramConfigs1D)
  {
    canvas[histName.first]->cd();

    recHasEntries = histsReconstructed[histName.first]->GetEntries() > 0;
    fitHasEntries = histsFittedSignal[histName.first]->GetEntries() > 0;

    // Jeśli żaden histogram nie ma wpisów - pomiń ten canvas
    if (!fitHasEntries && !recHasEntries)
    {
      continue;
    }

    TRegexp pattern("vtxNeu.*Fit"), pattern1("Chi2Comp.*"), pattern2("trk.*"), patternpull("pull.*"), patternPhi(".*Phi");

    if (histName.first.Contains(pattern) || histName.first.Contains(pattern1))
    {
      histsReconstructed[histName.first]->GetYaxis()->SetRangeUser(1, std::max(histsReconstructed[histName.first]->GetMaximum(), histsFittedSignal[histName.first]->GetMaximum()) * 100);
      canvas[histName.first]->SetLogy(1);
    }
    else
    {
      histsReconstructed[histName.first]->GetYaxis()->SetRangeUser(0, std::max(histsReconstructed[histName.first]->GetMaximum(), histsFittedSignal[histName.first]->GetMaximum()) * 1.2);
      canvas[histName.first]->SetLogy(0);
    }

    histsReconstructed[histName.first]->SetTitle("");

    histsReconstructed[histName.first]->SetLineColor(kBlue);

    Bool_t fitcond = histName.first == "curv1" || histName.first == "phiv1" || histName.first == "cotv1" ||
                     histName.first == "curv2" || histName.first == "phiv2" || histName.first == "cotv2" ||
                     histName.first == "vtxNeu_x" || histName.first == "vtxNeu_y" || histName.first == "vtxNeu_z" || histName.first == "delta_t_boostlor_res" || histName.first == "delta_t_rec_res" || histName.first == "phi_vtx_x" || histName.first == "phi_vtx_y" || histName.first == "phi_vtx_z" || histName.first == "vtxCh_x" || histName.first == "vtxCh_y" || histName.first == "vtxCh_z" || histName.first.Contains(pattern2) || histName.first.Contains(patternpull) || histName.first.Contains(patternPhi);

    if (fitcond && fitHasEntries)
    {
      Bool_t fitSuccess = false;

      std::cout << "\n=== Fitting delta_t with Double Gaussian ===" << std::endl;

      // Użyj starej, sprawdzonej metody TF1 dla delta_t
      fitSuccess = doubleFitter->FitHistogram(histsFittedSignal[histName.first], KLOE::DoublGaussFitter::FitType::kDoubleGauss);

      if (fitSuccess)
      {
        // Pobierz wyniki fitu
        const KLOE::DoublGaussFitter::FitResult &result = doubleFitter->GetLastResults();

        // Wypisz wyniki do konsoli
        std::cout << "\n=== Fit results for " << histName.first << " ===" << std::endl;
        std::cout << "Status: " << result.status << ", Converged: " << result.converged << std::endl;
        std::cout << "Chi2/NDF = " << result.chi2 << "/" << result.ndf
                  << " = " << doubleFitter->GetChi2NDF() << std::endl;
        std::cout << "Combined mean = " << result.mean
                  << " ± " << result.meanErr << std::endl;
        std::cout << "Combined sigma = " << result.coreSigma
                  << " ± " << result.coreSigmaErr << std::endl;

        // Wypisz parametry poszczególnych gaussów
        for (int i = 0; i < 2; i++)
        {
          int idx = i * 2;
          std::cout << "  Gauss " << (i + 1) << ": A=" << result.parameters[idx]
                    << ", μ=" << result.parameters[idx + 1]
                    << ", σ=" << result.parameters[idx + 2] << std::endl;
        }
      }
      else
      {
        std::cerr << "WARNING: Fit failed for " << histName.first << std::endl;
      }
    }

    if (recHasEntries)
    {
      histsReconstructed[histName.first]->Draw("HIST");
      histsFittedSignal[histName.first]->SetLineColor(kRed);
      histsFittedSignal[histName.first]->Draw("HIST SAME");
    }
    else if (fitHasEntries)
    {
      histsFittedSignal[histName.first]->SetLineColor(kRed);
      histsFittedSignal[histName.first]->Draw("HIST");
      histsReconstructed[histName.first]->Draw("HIST SAME");
    }

    // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
    TLegend *legend;

    if (!fitcond && histName.first != "TransvRadius" && histName.first != "Chi2SignalKinFit" && histName.first != "Chi2TriKinFit" && histName.first != "prob_signal")
    {
      legend = new TLegend(0.15, 0.75, 0.4, 0.9, "", "NDC");
      if (recHasEntries)
        legend->AddEntry(histsReconstructed[histName.first], "Reconstructed", "l");
      if (fitHasEntries)
        legend->AddEntry(histsFittedSignal[histName.first], "Fitted Signal", "l");
    }
    else if (histName.first == "TransvRadius")
    {
      legend = new TLegend(0.65, 0.75, 0.9, 0.9, "", "NDC");
      legend->AddEntry(histsReconstructed[histName.first], "Charged", "l");
      legend->AddEntry(histsFittedSignal[histName.first], "Neutral", "l");
    }
    else
    {
      legend = new TLegend(0.15, 0.7, 0.45, 0.9, "", "NDC");
      if (recHasEntries)
        legend->AddEntry(histsReconstructed[histName.first], "Reconstructed", "l");
      if (fitHasEntries)
        legend->AddEntry(histsFittedSignal[histName.first], "Fitted Data", "l");

      // Sprawdź czy fit się udał przed użyciem wyników
      if (doubleFitter->IsLastFitSuccessful())
      {
        const KLOE::DoublGaussFitter::FitResult &r = doubleFitter->GetLastResults();
        legend->AddEntry((TObject *)0, Form("#chi^{2}/NDF = %.2f", doubleFitter->GetChi2NDF()), "");
        legend->AddEntry((TObject *)0, Form("#mu = %.3f #pm %.3f", r.mean, r.meanErr), "");
        legend->AddEntry((TObject *)0, Form("#sigma = %.3f #pm %.3f", r.coreSigma, r.coreSigmaErr), "");

        // Narysuj fit z komponentami
        doubleFitter->DrawFitOnCurrentPad(true, false); // drawComponents=true, showStats=false
      }
      else
      {
        legend->AddEntry((TObject *)0, "Fit FAILED", "");
      }
    }

    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.03);
    legend->Draw();

    canvas[histName.first]->Update();

    canvas[histName.first]->SaveAs(Form(folderPath + "/%s_comparison%s", histName.first.Data(), Paths::ext_img.Data()));
  }

  for (const auto &config : histogramConfigs2D)
  {

    Bool_t withProfile = (config.first == "t_ch_fit_vs_t_ch_mc" || config.first == "t_neu_fit_vs_t_neu_mc" || config.first == "t_ch_rec_vs_t_ch_mc" || config.first == "t_neu_rec_vs_t_neu_mc" || config.first == "delta_t_mc_vs_delta_t_res" || config.first == "delta_t_vs_delta_t_mc");

    TH2 *h2D = hists2DFittedSignal[config.first];

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

    canvas2D[config.first]->cd();
    canvas2D[config.first]->SetLogz(1);

    hists2DFittedSignal[config.first]->Draw("COLZ");

    // Sprawdź czy JAKIKOLWIEK histogram ma wpisy
    Bool_t hasEntries = kFALSE;

    // Sprawdź wszystkie kanały

    if (hists2DFittedSignal[config.first]->GetEntries() > 0)
    {
      hasEntries = kTRUE;
    }

    // Jeśli żaden histogram nie ma wpisów - pomiń ten canvas
    if (!hasEntries)
    {
      continue;
    }

    canvas2D[config.first]->SaveAs(Form("%s/%s_2D%s", folderPath.Data(), config.first.Data(), Paths::ext_img.Data()));
  }

  Float_t fullyGoodPer = CalculateEfficiency(numberOfAllGood, numberOfAllGood + numberOfAtLeastOneBad) * 100,
          preselectionEff = CalculateEfficiency(cutter->GetTotalSignalExcludingMinus1(), cutter->GetTotalSignal()) * 100;

  std::cout << "How many reconstructed events had all good clusters? " << numberOfAllGood << std::endl;
  std::cout << "How many reconstructed events had at least one bad cluster? " << numberOfAtLeastOneBad << std::endl;
  std::cout << "Percentage of fully good events: " << fullyGoodPer << " %" << std::endl
            << std::endl;

  std::cout << std::endl;

  std::cout << "Signal events without errors (entire space): " << signal_tot << " over " << signal_tot_err << " (" << preselectionEff << " %) --> Efficiency of preselection" << std::endl;

  // Generuj raport z cięć
  std::cout << "Generating cuts report..." << std::endl;
  try
  {
    cutter->GenerateReport(Paths::reportConfigPath, folderPath.Data());
    std::cout << "Report generated successfully in: " << folderPath << std::endl;
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error generating report: " << e.what() << std::endl;
  }

  std::cout << std::endl;
  std::cout << "=== MC Fit Comparison Terminated ===" << std::endl;
}

Double_t MC_fit_comparison::CalculatePurity(Int_t signal, Int_t total) const
{
  if (total == 0)
    return 0.0; // Unikaj dzielenia przez zero
  return static_cast<Double_t>(signal) / total;
}

Double_t MC_fit_comparison::CalculateEfficiency(Int_t signal, Int_t total) const
{
  if (total == 0)
    return 0.0; // Unikaj dzielenia przez zero

  return static_cast<Double_t>(signal) / total;
}

void MC_fit_comparison::FolderManagement(TString folderPath) const
{
  // Usuń folder, jeśli istnieje
  system(Form("rm -rf %s", folderPath.Data()));

  // Utwórz folder
  system(Form("mkdir -p %s", folderPath.Data()));
}

void MC_fit_comparison::InitializeCutSelector(const TString &option)
{
  cutter = new StatisticalCutter(Paths::cutlimitsName, 1, KLOE::HypothesisCode::SIMONA_ANALYSIS);

  cutValues = {
      {"RtCh", [this]()
       { return RtCh; }},
      {"RtNeu", [this]()
       { return RtNeu; }},
      {"Chi2SignalReduced", [this]()
       { return Chi2SignalReduced; }},
      {"DeltaPhivFit", [this]()
       { return deltaPhiFit; }},
      {"InvMassKch", [this]()
       { return Kchrec[5]; }},
      {"InvMassKne", [this]()
       { return (*minv4gam); }},
      {"InvMassPi01", [this]()
       { return pi01Fit[5]; }},
      {"InvMassPi02", [this]()
       { return pi02Fit[5]; }},
      {"Pi01OmegaKineticEnergy", [this]()
       { return T0Omega; }},
      {"MassOmega", [this]()
       { return omegaFit[5]; }},
      {"T0OmegaUpperLimit", [this]()
       { return a * T0Omega + b + Breal; }},
      {"T0OmegaLowerLimit", [this]()
       { return a * T0Omega + b - Breal; }},
      {"RtChOmega", [this]()
       { return RtChRec; }},
      {"RtNeuOmega", [this]()
       { return RtNeuRec; }},
      {"ZChOmega", [this]()
       { return ZChRec; }},
      {"ZNeuOmega", [this]()
       { return ZNeuRec; }}};

  centralValues = {
      {"DeltaPhivFit", 3.130},
      {"InvMassKch", 497.605},
      {"InvMassKne", 489.467},
      {"InvMassPi01", 134.954},
      {"InvMassPi02", 134.841},
      {"Pi01OmegaKineticEnergy", 155.658},
      {"MassOmega", 782.994}};

  cutLimits = {
      {"T0OmegaUpperLimit", [this]()
       { return omegaFit[5]; }},
      {"T0OmegaLowerLimit", [this]()
       { return omegaFit[5]; }}};

  cutter->SetTree(fChain);

  // Rejestruj gettery tylko dla rzeczywistych (nie-syntetycznych) cięć
  // Syntetyczne cięcia są obsługiwane automatycznie na podstawie swoich członków
  for (const auto &cutPair : cutValues)
  {
    cutter->RegisterVariableGetter(cutPair.first, cutPair.second);
  }

  for (const auto &centralPair : centralValues)
  {
    Float_t centralValue = centralPair.second;
    cutter->RegisterCentralValueGetter(centralPair.first, [centralValue]()
                                       { return centralValue; });
  }

  for (const auto &limitPair : cutLimits)
  {
    cutter->RegisterCutValueGetter(limitPair.first, limitPair.second);
  }

  const auto &cuts = cutter->GetCuts();
  for (size_t i = 0; i < cuts.size(); ++i)
  {
    cutNameToIndex[cuts[i].cutId] = i;
  }

  std::cout << "StatisticalCutter initialized with " << cuts.size() << " cuts:" << std::endl;
  auto visibleCuts = cutter->GetVisibleCuts();
  std::cout << "  Visible cuts: " << visibleCuts.size() << std::endl;
  for (size_t i : visibleCuts)
  {
    const auto &cut = cuts[i];
    if (cut.isSyntheticGroup)
    {
      std::cout << "    [" << i << "] " << cut.cutId << " (SYNTHETIC GROUP with " << cut.groupMembers.size() << " members)" << std::endl;
    }
    else
    {
      std::cout << "    [" << i << "] " << cut.cutId << std::endl;
    }
  }

  // Pokaż ukryte cięcia (członkowie grup)
  auto hiddenCount = cuts.size() - visibleCuts.size();
  if (hiddenCount > 0)
  {
    std::cout << "  Hidden cuts (group members): " << hiddenCount << std::endl;
    for (size_t i = 0; i < cuts.size(); ++i)
    {
      if (cutter->IsCutInGroup(i))
      {
        std::cout << "    [" << i << "] " << cuts[i].cutId << " (member of group)" << std::endl;
      }
    }
  }
}

std::vector<size_t> MC_fit_comparison::GetCutIndicesForOption(const TString &option)
{
  std::vector<size_t> indices;

  if (option.IsNull() || option == "NO_CUTS")
  {
    std::cout << "No cuts applied" << std::endl;
    cutter->ClearActiveCuts();
    return indices;
  }

  TObjArray *tokens = option.Tokenize("+");

  for (Int_t i = 0; i < tokens->GetEntriesFast(); ++i)
  {
    TString token = ((TObjString *)tokens->At(i))->GetString();

    auto it = cutNameToIndex.find(std::string(token.Data()));
    if (it != cutNameToIndex.end())
    {
      indices.push_back(it->second);
    }
    else if (token != "")
    {
      std::cerr << "Warning: Cut '" << token << "' not found in loaded cuts" << std::endl;
    }
  }

  delete tokens;

  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());

  std::cout << "Option '" << option << "' uses " << indices.size() << " cuts:" << std::endl;
  const auto &cuts = cutter->GetCuts();
  for (size_t idx : indices)
  {
    if (idx < cuts.size())
      std::cout << "  - " << cuts[idx].cutId << std::endl;
  }

  // Ustawić aktywne cięcia w cutter
  cutter->SetActiveCuts(indices);

  return indices;
}

TString MC_fit_comparison::SanitizeFolderName(const TString &option)
{
  TString sanitized = option;

  sanitized.ReplaceAll("+", "_");
  sanitized.ReplaceAll(" ", "_");
  sanitized.ReplaceAll("-", "_");

  sanitized.ToLower();

  return sanitized;
}

// Funkcja do tworzenia wykresu RMS vs X
TGraphErrors *MC_fit_comparison::CreateRMSProfile(TH2 *h2D, const char *name, const char *title)
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

TCanvas *MC_fit_comparison::CreateCanvasWithProfiles(TH2 *h2D, const TString &name,
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