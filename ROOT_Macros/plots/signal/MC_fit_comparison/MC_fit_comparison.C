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
    deltaTheta;

// Definitions of cuts

StatisticalCutter *cutter;

std::map<std::string, std::function<Float_t()>> cutValues;
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

  Float_t weight = 1.0;

  if ((*mctruth == 1 || *mctruth == 0) && *mcflag == 1)
    weight = 1; // interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

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

  Float_t T0Omega = 0;

  if (abs(Omega1ErrTmp) > abs(Omega2ErrTmp))
    T0Omega = Omega2MassTmp - (trk1[3] + trk2[3]) - pi0Omega2[5];
  else
    T0Omega = Omega1MassTmp - (trk1[3] + trk2[3]) - pi0Omega1[5];

  Float_t deltaTfit = timesSignalFit.deltaTimeCM,
          deltaT = timesBoostLor.deltaTimeCM,
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

  ///////////////////////////////////////////////////////////////////////////////
  // Analiza Simony cięcie na 3 sigma mas
  Bool_t condMassKch = abs(Kchrec[5] - 497.606) < 3 * 1.099,
         condMassKne = abs(*minv4gam - 489.718) < 3 * 38.863,
         condMassPi01 = abs(pi01Fit[5] - 134.894) < 3 * 6.291,
         condMassPi02 = abs(pi02Fit[5] - 134.745) < 3 * 6.772;
  ///////////////////////////////////////////////////////////////////////////////

  Bool_t simonaCuts = abs(deltaPhiFit - 3.132) > 2 * 0.167 && *Chi2SignalKinFit < 30.,
         simonaKinCuts = condMassKch && condMassKne && condMassPi01 && condMassPi02 && simonaCuts;

  if (*mctruth == 0 || *mctruth == -1 || *mctruth == 1)
    signal_tot_err++;

  if ((*mctruth == 0 || *mctruth == 1))
    signal_tot++;

  Double_t limitRadiusChMC = 25.,
           limitRadiusNeMC = 25.;

  Double_t
      pathKneMC = sqrt(pow(Knemc[6] - ipmc[0], 2) +
                       pow(Knemc[7] - ipmc[1], 2)),
      pathKchMC = sqrt(pow(Kchmc[6] - ipmc[0], 2) +
                       pow(Kchmc[7] - ipmc[1], 2));

  bool passesCuts = true;
  if (!cutIndices.empty())
  {
    for (size_t cutIdx : cutIndices)
    {
      if (!cutter->PassCut(cutIdx))
      {
        passesCuts = false;
        break;
      }
    }
  }

  // Update statistics for all events
  cutter->UpdateStats(*mctruth);

  if (!passesCuts)
    return kTRUE;

  if (*mctruth == 1)
  {
    Int_t mctruth_tmp = *mctruth;

    if (*mctruth == 0)
      mctruth_tmp = 1;

    if (*goodClustersTriKinFitSize < 4)
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

    errorPart10 = pow(abs(trk1[0]) - abs(CurvMC[0]), 2) +
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

    if (fitcond)
    {
      Bool_t fitSuccess = false;

      // ✅ SPECJALNE TRAKTOWANIE DLA delta_t
      if (1) // histName.first == "delta_t")
      {
        std::cout << "\n=== Fitting delta_t with TF1 (proven method) ===" << std::endl;

        // Użyj starej, sprawdzonej metody TF1 dla delta_t
        fitSuccess = doubleFitter->FitHistogram(histsFittedSignal[histName.first], KLOE::DoublGaussFitter::FitType::kDoubleGauss);
      }
      else
      {
        // ✅ Dla innych zmiennych użyj automatycznych parametrów
        fitSuccess = doubleFitter->FitHistogram(histsReconstructed[histName.first]);
      }

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

    // histsReconstructed[histName.first]->GetYaxis()->SetRangeUser(0.0, 1.5 * std::max(histsReconstructed[histName.first]->GetMaximum(), histsFittedSignal[histName.first]->GetMaximum()));

    histsReconstructed[histName.first]->Draw();
    histsFittedSignal[histName.first]->SetLineColor(kRed);
    histsFittedSignal[histName.first]->Draw("SAME");

    // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
    TLegend *legend;

    if (!fitcond && histName.first != "TransvRadius")
    {
      legend = new TLegend(0.15, 0.75, 0.4, 0.9, "", "NDC");
      legend->AddEntry(histsReconstructed[histName.first], "Reconstructed", "l");
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
      legend->AddEntry(histsReconstructed[histName.first], "Reconstructed", "l");
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

    canvas[histName.first]->SaveAs(Form(folderPath + "/%s_comparison.png", histName.first.Data()));
  }

  for (const auto &config : histogramConfigs2D)
  {

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

    canvas2D[config.first]->SaveAs(Form("%s/%s_2D.png", folderPath.Data(), config.first.Data()));
  }

  Float_t fullyGoodPer = CalculateEfficiency(numberOfAllGood, numberOfAllGood + numberOfAtLeastOneBad) * 100,
          preselectionEff = CalculateEfficiency(signal_tot, signal_tot_err) * 100,
          analysisEff = CalculateEfficiency(signal_num, signal_tot) * 100,
          totalEff = CalculateEfficiency(signal_num, signal_tot_err) * 100;

  std::cout << "How many reconstructed events had all good clusters? " << numberOfAllGood << std::endl;
  std::cout << "How many reconstructed events had at least one bad cluster? " << numberOfAtLeastOneBad << std::endl;
  std::cout << "Percentage of fully good events: " << fullyGoodPer << " %" << std::endl
            << std::endl;

  std::cout << "Signal events without errors: " << signal_tot << " over " << signal_tot_err << " (" << preselectionEff << " %) --> Efficiency of preselection" << std::endl;
  std::cout << "Signal events after cuts: " << signal_num << " over " << signal_tot << " (" << analysisEff << " %) --> Efficiency of analysis" << std::endl;

  std::cout << "Total Efficiency: " << totalEff << " %" << std::endl;

  for (const auto &cutId : cutIndices)
  {
    const auto efficiency = cutter->GetEfficiency(cutId);

    std::cout << "Applied cut: " << cutId << " with eff " << efficiency << std::endl;
  }

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
       { return pi02Fit[5]; }}};

  centralValues = {
      {"DeltaPhivFit", 3.130},
      {"InvMassKch", 497.605},
      {"InvMassKne", 489.467},
      {"InvMassPi01", 134.954},
      {"InvMassPi02", 134.841}};

  cutter->SetTree(fChain);

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

  const auto &cuts = cutter->GetCuts();
  for (size_t i = 0; i < cuts.size(); ++i)
  {
    cutNameToIndex[cuts[i].cutId] = i;
  }

  std::cout << "StatisticalCutter initialized with " << cuts.size() << " cuts:" << std::endl;
  for (size_t i = 0; i < cuts.size(); ++i)
  {
    std::cout << "  [" << i << "] " << cuts[i].cutId << std::endl;
  }
}

std::vector<size_t> MC_fit_comparison::GetCutIndicesForOption(const TString &option)
{
  std::vector<size_t> indices;

  if (option.IsNull() || option == "NO_CUTS")
  {
    std::cout << "No cuts applied" << std::endl;
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