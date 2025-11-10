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

#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <TGraphAsymmErrors.h>

#include <interf_function.h>
#include <TEfficiency.h>

namespace KH = KLOE::Histograms;

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

TCanvas *canvaEff;

TH1 *deltaTSignalTot;

KLOE::TripleGaussFitter *fitter;
KLOE::DoublGaussFitter *doubleFitter;

KLOE::pm00 Obj;
KLOE::NeutralReconstruction *neutReconst;

TEfficiency *efficiency;

TH1 *histCounts;

TF1 *chi2DistFunc;

// Wczytaj konfiguracje histogramów
auto histogramConfigs1D = KLOE::Histograms::LoadHistogramConfigs1D(Paths::histogramConfig1DPath);
auto histogramConfigs2D = KLOE::Histograms::LoadHistogramConfigs2D(Paths::histogramConfig2DPath);
///////////////////////////////////////////////////////////////////////////////

void signal_vs_bcg_v2::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  // Creation of output folder with cut names
  /////////////////////////////////////////
  Cuts cuts;
  TString option = GetOption();
  fOption = option;

  if (option.IsNull() || !cuts.Contains(option))
    option = "NO_CUTS";
  folderPath = "img/" + option;
  FolderManagement(folderPath);
  /////////////////////////////////////////

  chi2DistFunc = new TF1("chi2DistFunc", chi2dist, 0, 1E5, 2); // 2 parameters

  KLOE::setGlobalStyle();

  fitter = new KLOE::TripleGaussFitter();
  doubleFitter = new KLOE::DoublGaussFitter();
  neutReconst = new KLOE::NeutralReconstruction(4);

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
  deltaTSignalTot = (TH1F *)histsReconstructed["delta_t"]["Signal"]->Clone("EfficiencyHistTot");
  deltaTSignalTot->Reset();

  histCounts = (TH1F *)histsReconstructed["chi2_signalKinFit"]["Signal"]->Clone("CountsHist");
  histCounts->Reset();
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
}

void signal_vs_bcg_v2::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
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

  Double_t dataPCA[2];

  fReader.SetLocalEntry(entry);

  Float_t phiv1PhiCM, phivPlusPhiCM, phiv2PhiCM, phivMinusPhiCM, deltaPhiPhiCM;

  TVector3 boostPhi(-*Bpx / *Broots, -*Bpy / *Broots, -*Bpz / *Broots);

  TLorentzVector trk1LAB(trk1[0], trk1[1], trk1[2], trk1[3]), trk1PhiCM(0, 0, 0, 0);
  TLorentzVector trk2LAB(trk2[0], trk2[1], trk2[2], trk2[3]), trk2PhiCM(0, 0, 0, 0);

  Obj.lorentz_transf(boostPhi, trk1LAB, trk1PhiCM);
  Obj.lorentz_transf(boostPhi, trk2LAB, trk2PhiCM);

  phiv1PhiCM = trk1PhiCM.Phi();
  phiv2PhiCM = trk2PhiCM.Phi();

  // Float_t Knerec[9];

  // Knerec[0] = gammaMomTriangle1[0] + gammaMomTriangle2[0] +
  //             gammaMomTriangle3[0] + gammaMomTriangle4[0];
  // Knerec[1] = gammaMomTriangle1[1] + gammaMomTriangle2[1] +
  //             gammaMomTriangle3[1] + gammaMomTriangle4[1];
  // Knerec[2] = gammaMomTriangle1[2] + gammaMomTriangle2[2] +
  //             gammaMomTriangle3[2] + gammaMomTriangle4[2];
  // Knerec[3] = gammaMomTriangle1[3] + gammaMomTriangle2[3] +
  //             gammaMomTriangle3[3] + gammaMomTriangle4[3];

  // Knerec[4] = sqrt(pow(Knerec[0], 2) + pow(Knerec[1], 2) + pow(Knerec[2], 2));

  Float_t vKchFit = PhysicsConstants::cVel * KchboostFit[4] / KchboostFit[3],
          pathKchFit = sqrt(pow(KchboostFit[6] - ipFit[0], 2) +
                            pow(KchboostFit[7] - ipFit[1], 2) +
                            pow(KchboostFit[8] - ipFit[2], 2)),
          RKchFit = sqrt(pow(KchboostFit[6] - ipFit[0], 2) +
                         pow(KchboostFit[7] - ipFit[1], 2)),
          tKchFit = KchboostFit[9] / 0.0895,
          vKneFit = PhysicsConstants::cVel * KnereclorFit[4] / KnereclorFit[3],
          pathKneFit = sqrt(pow(KnerecFit[6] - ipFit[0], 2) +
                            pow(KnerecFit[7] - ipFit[1], 2) +
                            pow(KnerecFit[8] - ipFit[2], 2)),
          // pow(KnerecFit[8] - ipFit[2], 2)),
      RKneFit = sqrt(pow(KnerecFit[6] - ipFit[0], 2) +
                     pow(KnerecFit[7] - ipFit[1], 2)),
          tKneFit = KnereclorFit[9] / 0.0895,
          vKneMC = PhysicsConstants::cVel * Knemc[4] / Knemc[3],
          vKchMC = PhysicsConstants::cVel * Kchmc[4] / Kchmc[3],
          vKne = PhysicsConstants::cVel * Knerec[4] / Knerec[3],
          pathKne = sqrt(pow(Knerec[6] - ip[0], 2) +
                         pow(Knerec[7] - ip[1], 2) +
                         pow(Knerec[8] - ip[2], 2)),
          RKne = sqrt(pow(Knerec[6] - ip[0], 2) +
                      pow(Knerec[7] - ip[1], 2)),
          RKch = sqrt(pow(Kchrec[6] - ip[0], 2) +
                      pow(Kchrec[7] - ip[1], 2)),
          tKne = pathKne / (vKne * 0.0895),
          pathKchMC = sqrt(pow(Kchmc[6] - ipmc[0], 2) +
                           pow(Kchmc[7] - ipmc[1], 2)),
          // pow(Kchmc[8] - ipmc[2], 2)),
      pathKneMC = sqrt(pow(Knemc[6] - ipmc[0], 2) +
                       pow(Knemc[7] - ipmc[1], 2));
  //  pow(Knemc[8] - ipmc[2], 2));

  std::vector<Float_t> kaonMom1 = {Kchboost[0],
                                   Kchboost[1],
                                   Kchboost[2],
                                   Kchboost[3]},
                       kaonMom2 = {Knerec[0],
                                   Knerec[1],
                                   Knerec[2],
                                   Knerec[3]},
                       kaonPos1 = {Kchboost[6],
                                   Kchboost[7],
                                   Kchboost[8],
                                   Kchboost[9]},
                       kaonPos2 = {Knerec[6],
                                   Knerec[7],
                                   Knerec[8],
                                   Knerec[9]},
                       ipVec = {ip[0], ip[1], ip[2]};

  KLOE::KaonProperTimes propTimes = Obj.CalculateKaonProperTimes(kaonMom1,
                                                                 kaonPos1,
                                                                 kaonMom2,
                                                                 kaonPos2,
                                                                 ipVec);

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

  Float_t deltaTfit = *KaonChTimeCMSignalFit - *KaonNeTimeCMSignalFit,
          deltaT = propTimes.deltaTimeCM,
          deltaTMC = *KaonChTimeCMMC - *KaonNeTimeCMMC;

  Float_t combinedMassPi0Fit = sqrt(pow(pi01Fit[5] - PhysicsConstants::mPi0, 2) +
                                    pow(pi02Fit[5] - PhysicsConstants::mPi0, 2)),
          combinedMassPi0 = sqrt(pow(pi01[5] - PhysicsConstants::mPi0, 2) +
                                 pow(pi02[5] - PhysicsConstants::mPi0, 2));

  Float_t kaonChPath = sqrt(pow(Kchrec[6] - ip[0], 2) +
                            pow(Kchrec[7] - ip[1], 2) +
                            pow(Kchrec[8] - ip[2], 2)),
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
          openingAngleNeutral = pi01Vec.Angle(pi02Vec.Vect()) * 180.0 / TMath::Pi(),
          acosCutAngle = acos(-0.8) * 180.0 / TMath::Pi();

  Float_t QmissFit = sqrt(pow(KchboostFit[0] - KchrecFit[0], 2) +
                          pow(KchboostFit[1] - KchrecFit[1], 2) +
                          pow(KchboostFit[2] - KchrecFit[2], 2) +
                          pow(KchboostFit[3] - KchrecFit[3], 2));

  Float_t weight = 1.0;

  if ((*mctruth == 1 || *mctruth == 0) && *mcflag == 1)
    weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

  TVector3 z_axis(0., 0., 1.),
      gamma1(gammaMomTriangle1[0], gammaMomTriangle1[1], gammaMomTriangle1[2]),
      gamma2(gammaMomTriangle2[0], gammaMomTriangle2[1], gammaMomTriangle2[2]),
      gamma3(gammaMomTriangle3[0], gammaMomTriangle3[1], gammaMomTriangle3[2]),
      gamma4(gammaMomTriangle4[0], gammaMomTriangle4[1], gammaMomTriangle4[2]);

  Float_t thetaGamma1 = gamma1.Angle(z_axis) * 180.0 / TMath::Pi(),
          thetaGamma2 = gamma2.Angle(z_axis) * 180.0 / TMath::Pi(),
          thetaGamma3 = gamma3.Angle(z_axis) * 180.0 / TMath::Pi(),
          thetaGamma4 = gamma4.Angle(z_axis) * 180.0 / TMath::Pi();

  Bool_t passSemi = abs(*Qmiss - 40.73) < 15 && openingAngleNeutral > 160. && openingAngleCharged > 160.;

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

  deltaPhiPhiCM = abs(phiv2PhiCM - phiv1PhiCM);
  deltaPhi = abs(Phi2Rec - Phi1Rec);
  deltaTheta = abs(Theta2Fit - Theta1Fit);
  deltaPhiFit = abs(Phi2Fit - Phi1Fit);

  // Analiza Simony ciecie na phi bad
  const Double_t SLOPE = -10.0 / 9.0;

  Bool_t condGeneral = (deltaTfit - deltaTMC)<-5.0,
                                              condLowerLimit = (deltaTfit - deltaTMC)> SLOPE *
                       deltaTMC,
         condUpperLimit = (deltaTfit - deltaTMC) < (SLOPE * (deltaTMC - 2.0)),
         badClusSimona = condGeneral && *Chi2SignalKinFit < 30.;
  ///////////////////////////////////////////////////////////////////////////////
  // Analiza Simony cięcie na 3 sigma mas
  Bool_t condMassKch = abs(Kchrec[5] - 497.613) < 3 * 1.074,
         condMassKne = abs(*minv4gam - 489.274) < 3 * 41.319,
         condMassPi01 = abs(pi01Fit[5] - 134.964) < 3 * 4.214,
         condMassPi02 = abs(pi02Fit[5] - 134.805) < 3 * 4.057;

  ///////////////////////////////////////////////////////////////////////////////

  Float_t radius00 = sqrt(pow(KnerecFit[6] - ipFit[0], 2) +
                          pow(KnerecFit[7] - ipFit[1], 2)),
          radiuspm = sqrt(pow(Kchrec[6] - ip[0], 2) +
                          pow(Kchrec[7] - ip[1], 2)),
          zdist00 = abs(KnerecFit[8] - ip[2]),
          zdistpm = abs(Kchrec[8] - ip[2]),
          path00 = sqrt(pow(radius00, 2) + pow(zdist00, 2)),
          pathpm = sqrt(pow(radiuspm, 2) + pow(zdistpm, 2)),
          path00MC = sqrt(pow(sqrt(pow(Knemc[6] - ipmc[0], 2) +
                                   pow(Knemc[7] - ipmc[1], 2)),
                              2) +
                          pow(abs(Knemc[8] - ipmc[2]), 2)),
          pathpmMC = sqrt(pow(sqrt(pow(Kchmc[6] - ipmc[0], 2) +
                                   pow(Kchmc[7] - ipmc[1], 2)),
                              2) +
                          pow(abs(Kchmc[8] - ipmc[2]), 2)),
          path00MCCenter = sqrt(pow(sqrt(pow(Knemc[6], 2) +
                                         pow(Knemc[7], 2)),
                                    2) +
                                pow(abs(Knemc[8]), 2)),
          pathpmMCCenter = sqrt(pow(sqrt(pow(Kchmc[6], 2) +
                                         pow(Kchmc[7], 2)),
                                    2) +
                                pow(abs(Kchmc[8]), 2)),
          radiusLimit = 1,
          zdistLimit = 0.6;

  Float_t T0Omega = pi0OmegaFit1[3] - pi0OmegaFit1[5];

  Double_t numSigmaSimona = 3;

  Float_t limitRadiusNeMC = 25.,
          limitRadiusChMC = 25.;

  Bool_t mcflagCondition = (*mcflag == 1 && *mctruth >= 0) || *mcflag == 0,
         condAnalysisOld = *Chi2SignalKinFit < 40. && combinedMassPi0Fit < 35. && *TrcSum > -1 && abs(Kchrec[5] - PhysicsConstants::mK0) < 1.2 && abs(*minv4gam - PhysicsConstants::mK0) < 76. && *Qmiss < 3.75 && openingAngleCharged > acosCutAngle,
         simonaCuts = abs(deltaPhiFit - 3.132) > 2 * 0.167 && *Chi2SignalKinFit < 30.,
         simonaKinCuts = condMassKch && condMassKne && condMassPi01 && condMassPi02 && simonaCuts,
         shorterKaonPaths = pathKchMC < limitRadiusChMC && pathKneMC < limitRadiusNeMC,
         blobCut = *KaonNeTimeCMBoostTriFit - *KaonNeTimeCMBoostLor > 75.,
         noBlobCut = *KaonNeTimeCMBoostTriFit - *KaonNeTimeCMBoostLor <= 75.,
         simonaChi2Cut = *Chi2SignalKinFit <= 30,
         simonaPositionLimits = radius00 < radiusLimit && radiuspm < radiusLimit &&
                                zdist00 < zdistLimit && zdistpm < zdistLimit,
         omegaMassT0Cut = ((simonaPositionLimits && !(abs(T0Omega - 132.215) < numSigmaSimona * 31.848 && abs(omega[5] - 776.571) < numSigmaSimona * 28.026 && omega[5] < T0Omega + 649.516 && omega[5] > T0Omega + 640.176)) || !simonaPositionLimits) && simonaKinCuts;

  if ((*mctruth == 1 || *mctruth == -1 || *mctruth == 0) && *mcflag == 1)
    signal_tot++;

  if ((*mctruth == 1 || *mctruth == 0) && *mcflag == 1)
    signal_wo_err++;

  if ((*mctruth == 1) && *mcflag == 1)
  {
    deltaTSignalTot->Fill(deltaTMC, weight);
  }

  // Option to use in analysis
  TString option = GetOption();
  option.ToUpper();
  ////////////////////////////

  if (option.Contains("SHORTER_KAON_PATHS"))
    if (!shorterKaonPaths)
      return kTRUE;

  if (option.Contains("OLD_CUTS"))
    if (!condAnalysisOld)
      return kTRUE;

  if (option.Contains("SIMONA_CHI2_CUT"))
    if (!simonaChi2Cut)
      return kTRUE;

  if (option.Contains("BAD_CLUS_SIMONA"))
    if (!badClusSimona)
      return kTRUE;

  if (option.Contains("SIMONA_KIN_CUTS"))
    if (!simonaKinCuts)
      return kTRUE;

  if (option.Contains("SIMONA_ALL_CUTS"))
    if (!simonaCuts)
      return kTRUE;

  if (option.Contains("OMEGA_MASS_T0_CUT"))
    if (!omegaMassT0Cut)
      return kTRUE;

  if (option == "BLOB")
    if (!blobCut)
      return kTRUE;

  if (option == "NO_BLOB")
    if (!noBlobCut)
      return kTRUE;

  if (mcflagCondition)
  {
    Int_t mctruth_tmp = *mctruth;

    if (mctruth_tmp == 0 && *mcflag == 1)
      mctruth_tmp = 1;

    if ((mctruth_tmp == 1) && *mcflag == 1)
      signal_num++;

    if (*mctruth > 1 && *mcflag == 1)
      bkg_tot++;

    if ((*mcflag == 1 && *mctruth >= 0) || *mcflag == 0)
    {
      if (simonaChi2Cut)
        cutPassed++;
      else if (!simonaChi2Cut && *Chi2SignalKinFit < 100.)
        cutNPassed++;
      else if (!simonaChi2Cut && *Chi2SignalKinFit >= 100.)
        overflow++;

      // Fill histograms for reconstructed variables
      histsReconstructed["mass_Kch"][KLOE::channName.at(mctruth_tmp)]->Fill(Kchrec[5], weight);

      histsReconstructed["mass_Kne"][KLOE::channName.at(mctruth_tmp)]->Fill(*minv4gam, weight);

      histsReconstructed["mass_pi01"][KLOE::channName.at(mctruth_tmp)]->Fill(pi01[5], weight);
      histsReconstructed["mass_pi02"][KLOE::channName.at(mctruth_tmp)]->Fill(pi02[5], weight);

      histsReconstructed["time_neutral_MC"][KLOE::channName.at(mctruth_tmp)]->Fill(*TrcSum, weight);

      histsReconstructed["combined_mass_pi0"][KLOE::channName.at(mctruth_tmp)]->Fill(combinedMassPi0Fit, weight);

      // Fitted signal variables
      histsFittedSignal["mass_Kch"][KLOE::channName.at(mctruth_tmp)]->Fill(Kchrec[5], weight);

      histsFittedSignal["mass_Kne"][KLOE::channName.at(mctruth_tmp)]->Fill(*minv4gam, weight);

      histsFittedSignal["mass_pi01"][KLOE::channName.at(mctruth_tmp)]->Fill(pi0OmegaFit1[5], weight);
      histsFittedSignal["mass_pi02"][KLOE::channName.at(mctruth_tmp)]->Fill(pi0OmegaFit2[5], weight);

      histsFittedSignal["chi2_signalKinFit"][KLOE::channName.at(mctruth_tmp)]->Fill(*Chi2SignalKinFit / 10., weight);
      histCounts->Fill(*Chi2SignalKinFit / 10.);

      histsFittedSignal["chi2_trilaterationKinFit"][KLOE::channName.at(mctruth_tmp)]->Fill(*Chi2OmegaKinFit / 8., weight);
      histsFittedSignal["prob_signal"][KLOE::channName.at(mctruth_tmp)]->Fill(TMath::Prob(*Chi2SignalKinFit, 10), weight);

      histsFittedSignal["combined_mass_pi0"][KLOE::channName.at(mctruth_tmp)]->Fill(combinedMassPi0Fit, weight);

      for (Int_t i = 0; i < 39; i++)
      {
        histsFittedSignal["pull" + std::to_string(i + 1)][KLOE::channName.at(mctruth_tmp)]->Fill(pullsSignalFit[i], weight);
      }

      histsFittedSignal["time_neutral_MC"][KLOE::channName.at(mctruth_tmp)]->Fill(*TrcSum, weight);

      histsFittedSignal["openingAngleCharged"][KLOE::channName.at(mctruth_tmp)]->Fill(openingAngleCharged, weight);
      histsFittedSignal["openingAngleNeutral"][KLOE::channName.at(mctruth_tmp)]->Fill(openingAngleNeutral, weight);

      histsFittedSignal["Qmiss"][KLOE::channName.at(mctruth_tmp)]->Fill(*Qmiss, weight);

      histsReconstructed["delta_t"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaT, weight);
      histsFittedSignal["delta_t"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaTfit, weight);

      histsFittedSignal["delta_t_MC"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaTMC, weight);

      histsFittedSignal["deltaPhiv"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaPhi, weight);
      histsFittedSignal["deltaPhivFit"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaPhiFit, weight);

      histsFittedSignal["deltaTheta"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaTheta, weight);

      histsFittedSignal["TransvRadius"][KLOE::channName.at(mctruth_tmp)]->Fill(path00MCCenter, weight);

      histsFittedSignal["T0Omega"][KLOE::channName.at(mctruth_tmp)]->Fill(T0Omega, weight);

      histsFittedSignal["mass_omega"][KLOE::channName.at(mctruth_tmp)]->Fill(omegaFit[5], weight);

      histsFittedSignal["dist_z_Neu"][KLOE::channName.at(mctruth_tmp)]->Fill(KnerecFit[8] - ipFit[2], weight);
      histsFittedSignal["dist_z_Ch"][KLOE::channName.at(mctruth_tmp)]->Fill(KchrecFit[8] - ipFit[2], weight);


      if(*muonAlertPlus > 0 || *muonAlertMinus > 0)
        histsFittedSignal["muon_alert"][KLOE::channName.at(mctruth_tmp)]->Fill(1., weight);
      else
        histsFittedSignal["muon_alert"][KLOE::channName.at(mctruth_tmp)]->Fill(0., weight);

      hists2DFittedSignal["T0_omega_vs_mass_omega"][KLOE::channName.at(mctruth_tmp)]->Fill(T0Omega, omegaFit[5], weight);

      hists2DFittedSignal["delta_t_mc_vs_delta_t_res"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaTMC, deltaTfit - deltaTMC, weight);
      hists2DFittedSignal["delta_t_vs_delta_t_mc"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaTMC, deltaT - deltaTMC, weight);

      hists2DFittedSignal["delta_t_fit_vs_delta_t_mc"][KLOE::channName.at(mctruth_tmp)]->Fill(deltaTMC, deltaTfit, weight);

      hists2DFittedSignal["t_ch_fit_vs_t_ch_mc"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonChTimeCMMC, *KaonChTimeCMSignalFit - *KaonChTimeCMMC, weight);
      hists2DFittedSignal["t_neu_fit_vs_t_neu_mc"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMSignalFit - *KaonNeTimeCMMC, weight);

      hists2DFittedSignal["t_ch_fit_MC_vs_t_neu_fit_MC"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonChTimeCMSignalFit - *KaonChTimeCMMC, *KaonNeTimeCMSignalFit - *KaonNeTimeCMMC, weight);

      if (*KaonNeTimeCMSignalFit <= *KaonNeTimeCMMC - (KLOE::T0 / 2.))
        hists2DFittedSignal["t_ch_fit_vs_t_neu_fit"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonChTimeCMSignalFit, *KaonNeTimeCMSignalFit, weight);

      hists2DFittedSignal["chi2_signalKinFit_vs_chi2_trilaterationKinFit"][KLOE::channName.at(mctruth_tmp)]->Fill(*Chi2SignalKinFit, *Chi2OmegaKinFit, weight);

      hists2DFittedSignal["t00_tri_vs_t00_mc"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMBoostTriFit, weight);

      hists2DFittedSignal["t00_triangle_vs_t00_mc"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonNeTimeCMMC, *KaonNeTimeCMBoostLor, weight);

      hists2DFittedSignal["t00_triangle_vs_t00_tri"][KLOE::channName.at(mctruth_tmp)]->Fill(*KaonNeTimeCMBoostLor, *KaonNeTimeCMBoostTriFit, weight);

      hists2DFittedSignal["Rkne_vs_Rkch"][KLOE::channName.at(mctruth_tmp)]->Fill(RKne, RKch, weight);

      hists2DFittedSignal["chi2_signalKinFit_vs_mass_Kch"][KLOE::channName.at(mctruth_tmp)]->Fill(*Chi2SignalKinFit, Kchrec[5] - PhysicsConstants::mK0, weight);

      hists2DFittedSignal["chi2_signalKinFit_vs_mass_Kne"][KLOE::channName.at(mctruth_tmp)]->Fill(*Chi2SignalKinFit, Knerec[5] - PhysicsConstants::mK0, weight);

      hists2DFittedSignal["chi2_signalKinFit_vs_delta_t"][KLOE::channName.at(mctruth_tmp)]->Fill(*Chi2SignalKinFit, deltaTfit, weight);

      hists2DFittedSignal["TransvRadius_vs_delta_t"][KLOE::channName.at(mctruth_tmp)]->Fill(path00MCCenter, deltaTMC, weight);

      hists2DFittedSignal["Qmiss_vs_deltaPhi"][KLOE::channName.at(mctruth_tmp)]->Fill(*Qmiss, deltaPhi, weight);

      dataPCA[0] = T0Omega;
      dataPCA[1] = omega[5];

      AddDataToPCA(dataPCA);
    }
  }

  return kTRUE;
}

void signal_vs_bcg_v2::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  std::cout << "Overflow entries in chi2_signalKinFit: " << overflow << std::endl;

  std::cout << "Number of events passing cut: " << cutPassed << std::endl;
  std::cout << "Number of events NOT passing cut: " << cutNPassed << std::endl;
  std::cout << "Overflows NOT passing cut: " << overflow << std::endl;
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

  std::map<TString, std::map<TString, Int_t>> integrals;
  std::map<TString, Int_t> maxNum;
  std::map<TString, TString> maxChannel;

  for (const auto &config : histogramConfigs1D)
  {
    maxNum[config.first] = 0;
    for (const auto &channelType : KLOE::channName)
    {
      if (histsFittedSignal[config.first][channelType.second]->GetEntries() > 0. && channelType.second != "Data" && channelType.second != "MC sum" && channelType.second == "Signal")
      {
        histsReconstructed[config.first][channelType.second]->Scale(histsReconstructed[config.first][channelType.second]->GetEntries() / histsReconstructed[config.first][channelType.second]->Integral(0, histsReconstructed[config.first][channelType.second]->GetNbinsX() + 1));

        histsFittedSignal[config.first][channelType.second]->Scale(histsFittedSignal[config.first][channelType.second]->GetEntries() / histsFittedSignal[config.first][channelType.second]->Integral(0, histsFittedSignal[config.first][channelType.second]->GetNbinsX() + 1));
      }

      integrals[config.first][channelType.second] = histsFittedSignal[config.first][channelType.second]->GetMaximum();

      if (integrals[config.first][channelType.second] > maxNum[config.first])
      {
        maxNum[config.first] = integrals[config.first][channelType.second];
        maxChannel[config.first] = channelType.second;
      }
    }
  }

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

    Double_t scalefactorReconstructed = 1.0;
    Double_t scalefactorFittedSignal = 1.0;

    if (histsReconstructed[config.first]["MC sum"]->GetEntries() > 0 && histsReconstructed[config.first]["Data"]->GetEntries() > 0)
    {
      scalefactorReconstructed = histsReconstructed[config.first]["Data"]->GetEntries() / histsReconstructed[config.first]["MC sum"]->GetEntries();

      histsReconstructed[config.first]["MC sum"]->Scale(scalefactorReconstructed);
    }

    if (histsFittedSignal[config.first]["MC sum"]->GetEntries() > 0 && histsFittedSignal[config.first]["Data"]->GetEntries() > 0)
    {
      scalefactorFittedSignal = histsFittedSignal[config.first]["Data"]->GetEntries() / histsFittedSignal[config.first]["MC sum"]->GetEntries();
      histsFittedSignal[config.first]["MC sum"]->Scale(scalefactorFittedSignal);
    }

    for (const auto &channelType : KLOE::channName)
    {
      if (channelType.second != "Data" && channelType.second != "MC sum")
      {
        histsReconstructed[config.first][channelType.second]->Scale(scalefactorReconstructed);
        histsFittedSignal[config.first][channelType.second]->Scale(scalefactorFittedSignal);
      }
    }
  }

  Double_t meanMass = 0, sigmaMass = 0, meanT0 = 0, sigmaT0 = 0;

  TF1 *chi2 = new TF1("chi2", chi2dist, 0, 1E5, 2);

  for (const auto &config : histogramConfigs1D)
  {
    Bool_t fitDoubleGaus = (config.first == "mass_Kch" || config.first == "mass_Kne" || config.first == "mass_pi01" || config.first == "mass_pi02"),
           fitOmegaGaus = (config.first == "T0Omega" || config.first == "mass_omega"),
           fitSignalBadClus = \(config.first == "deltaPhivFit" && fOption == "BAD_CLUS_SIMONA");
    Bool_t chi2Fit = 0; //(config.first == "chi2_signalKinFit");

    Bool_t logCond = (config.first == "time_neutral_MC" || config.first == "Qmiss" || config.first == "prob_signal"), // || config.first == "chi2_signalKinFit"), // || config.first == "TransvRadius"),
        logCondX = 0;                                                                                                 //(config.first == "chi2_signalKinFit");

    std::vector<TString> labels;

    // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
    TLegend *legend = new TLegend(0.6, 0.65, 0.9, 0.9, "", "NDC");

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
      // if (channelType.second == "MC sum")
      //   continue;

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
          doubleFitter->FitHistogram(histsFittedSignal[config.first][channelType.second], KLOE::DoublGaussFitter::FitType::kSingleGauss);
          doubleFitter->DrawFitOnCurrentPad(true, false);

          if (config.first == "T0Omega")
          {
            meanT0 = doubleFitter->GetLastResults().mean;
            sigmaT0 = doubleFitter->GetLastResults().coreSigma;
          }
          else if (config.first == "mass_omega")
          {
            meanMass = doubleFitter->GetLastResults().mean;
            sigmaMass = doubleFitter->GetLastResults().coreSigma;
          }

          textresultsFitter->AddText(Form("Fit results:"));
          textresultsFitter->AddText(Form("Chi2 / ndof: %.3f / %d", doubleFitter->GetLastResults().chi2, doubleFitter->GetLastResults().ndf));
          textresultsFitter->AddText(Form("Mean: %.3f #pm %.3f", doubleFitter->GetLastResults().mean, doubleFitter->GetLastResults().meanErr));
          textresultsFitter->AddText(Form("Sigma: %.3f #pm %.3f", doubleFitter->GetLastResults().coreSigma, doubleFitter->GetLastResults().coreSigmaErr));
          textresultsFitter->Draw();
        }

        if (0) //(channelType.second == "Signal" && (fitSignalBadClus || chi2Fit)))
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
          // Amplituda = Integral / (sigma * sqrt(2*pi))
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
          histsFittedSignal["delta_t"][channelType.second]->SetLineColor(kBlue);
          histsFittedSignal["delta_t"][channelType.second]->Draw("HIST SAME");
          histsReconstructed["delta_t"][channelType.second]->SetLineColor(kRed);
          histsReconstructed["delta_t"][channelType.second]->Draw("HIST SAME");
          histsFittedSignal["delta_t_MC"][channelType.second]->SetLineColor(kBlack);
          histsFittedSignal["delta_t_MC"][channelType.second]->Draw("HIST SAME");
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

    labels.push_back(Form("Events: %d", (Int_t)histsFittedSignal[config.first]["MC sum"]->GetEntries()));

    // labels.push_back(Form("Chi2 (Data vs. MC Sum): %.3f", histsFittedSignal[config.first]["MC sum"]->Chi2Test(histsFittedSignal[config.first]["Data"], "WU CHI2/NDF")));

    DrawLabelOnHisto(labels);

    gPad->Update();

    // if (!fitTripleGaus && !fitOmegaGaus && !fitSignalBadClus && !(config.first == "TransvRadius"))
    //   legend->Draw();

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

    canvas[config.first]->SaveAs(Form("%s/%s_comparison.png", folderPath.Data(), config.first.Data()));
  }

  Double_t sigmas[2], means[2], vLong[2], vTransv[2];
  // Calculation of PCA components
  WidthOfCorrelatedHist(means, sigmas, vLong, vTransv);

  Double_t numSigmas = 3.0;
  Double_t vTransvSigma[2] = {numSigmas * sigmas[1] * vTransv[0],
                              numSigmas * sigmas[1] * vTransv[1]},
           point1[2] = {means[0] - vLong[0] * 200,
                        means[1] - vLong[1] * 200},
           point2[2] = {means[0] + vLong[0] * 200,
                        means[1] + vLong[1] * 200},
           a = (point1[1] - point2[1]) / (point1[0] - point2[0]),
           b = point1[1] - a * point1[0],
           B = numSigmas * sigmas[1] * (vTransv[1] + a * vTransv[0]),
           Breal = numSigmas * sigmas[1] * TMath::Sqrt(1 + a * a);

  std::cout << "Trend line: y = " << a << " * x + " << b << std::endl;
  std::cout << "B: " << B << std::endl;
  std::cout << "Breal: " << Breal << std::endl;

  for (const auto &config : histogramConfigs2D)
  {
    if ((config.first == "t_ch_fit_vs_t_ch_mc" || config.first == "t_neu_fit_vs_t_neu_mc" || config.first == "delta_t_fit_vs_delta_t_mc"))
    {
      TH2 *h2D = hists2DFittedSignal[config.first]["Signal"];

      // Stwórz canvas z oboma profilami
      TCanvas *c = CreateCanvasWithProfiles(h2D,
                                            config.first + "_with_profiles",
                                            kTRUE,  // Rysuj mean profile
                                            kTRUE); // Rysuj sigma profile

      c->SaveAs(folderPath + "/" + config.first + "_with_profiles.png");
      canvasProfiles[config.first] = c;
    }

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

        TLine *pc1_line_cut = new TLine((means[0] + vTransvSigma[0]) - vLong[0] * 100,
                                        (means[1] + vTransvSigma[1]) - vLong[1] * 100,
                                        (means[0] + vTransvSigma[0]) + vLong[0] * 100,
                                        (means[1] + vTransvSigma[1]) + vLong[1] * 100);

        TLine *pc2_line_cut = new TLine((means[0] - vTransvSigma[0]) - vLong[0] * 100,
                                        (means[1] - vTransvSigma[1]) - vLong[1] * 100,
                                        (means[0] - vTransvSigma[0]) + vLong[0] * 100,
                                        (means[1] - vTransvSigma[1]) + vLong[1] * 100);

        TLine *T0_left_cut = new TLine(meanT0 - 3 * sigmaT0, meanMass - 3 * sigmaMass,
                                       meanT0 - 3 * sigmaT0, meanMass + 3 * sigmaMass);
        TLine *T0_right_cut = new TLine(meanT0 + 3 * sigmaT0, meanMass - 3 * sigmaMass,
                                        meanT0 + 3 * sigmaT0, meanMass + 3 * sigmaMass);

        TLine *Mass_left_cut = new TLine(meanT0 - 3 * sigmaT0, meanMass - 3 * sigmaMass,
                                         meanT0 + 3 * sigmaT0, meanMass - 3 * sigmaMass);
        TLine *Mass_right_cut = new TLine(meanT0 - 3 * sigmaT0, meanMass + 3 * sigmaMass,
                                          meanT0 + 3 * sigmaT0, meanMass + 3 * sigmaMass);

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

      canvas2D[config.first][channelType.second]->SaveAs(Form("%s/%s_%s_2D.png", folderPath.Data(), config.first.Data(), channelType.second.Data()));
    }
  }

  Double_t scale = deltaTSignalTot->GetEntries() / deltaTSignalTot->Integral(0, deltaTSignalTot->GetNbinsX() + 1);

  deltaTSignalTot->Scale(scale);

  efficiency = new TEfficiency(*histsFittedSignal["delta_t_MC"]["Signal"], *deltaTSignalTot);

  efficiency->SetUseWeightedEvents(kTRUE);

  efficiency->SetStatisticOption(TEfficiency::kBUniform);

  // TF1 *f_const = new TF1("f_const", "pol0", -30., 30.);

  // f_const->SetParameter(0, 0.1);

  // TFitResultPtr fitResult = efficiency->Fit(f_const, "SR");

  // if (fitResult->IsValid())
  // {
  //   std::cout << "Fit results for efficiency:" << std::endl;
  //   fitResult->Print("V");

  //   fitResult->GetErrors();
  // }
  // else
  // {
  //   std::cout << "Fit failed!" << std::endl;
  // }

  // Przy założeniu, że eff i f_const są już zdefiniowane i fit został wykonany

  // // 1. Stworzenie histogramu dla Pulls
  // TH1F *h_pulls = new TH1F("h_pulls", "Pulls Distribution; (Eps_i - Eps_fit) / sigma_i; Entries", 50, -5.0, 5.0);

  // // 2. Pobranie dopasowanej stałej i jej błędu
  // Double_t eff_fit = f_const->GetParameter(0);
  // Double_t eff_fit_error = f_const->GetParError(0);

  // // 3. Iteracja po binach TEfficiency (iterujemy po histogramie mianownika)
  // TH1 *h_total = (TH1 *)efficiency->GetTotalHistogram();
  // for (int i = 1; i <= h_total->GetNbinsX(); ++i)
  // {
  //   // Sprawdzenie, czy bin jest w zakresie fitowania (opcjonalnie, ale zalecane)
  //   Double_t x_center = h_total->GetBinCenter(i);

  //   // Wydobycie wartości wydajności i jej błędu z TEfficiency
  //   Double_t epsilon_i = efficiency->GetEfficiency(i);
  //   // Użyj większego z błędów dla konserwatywnej analizy (upper/lower error)
  //   Double_t error_i = efficiency->GetEfficiencyErrorUp(i) > efficiency->GetEfficiencyErrorLow(i) ? efficiency->GetEfficiencyErrorUp(i) : efficiency->GetEfficiencyErrorLow(i);

  //   // Tylko dla binów z niezerową statystyką i błędem
  //   if (efficiency->GetTotalHistogram()->GetBinContent(i) > 0 && error_i > 0)
  //   {
  //     Double_t pull = (epsilon_i - eff_fit) / error_i;
  //     h_pulls->Fill(pull);
  //   }
  // }

  // // 4. Dopasowanie Gaussa do rozkładu Pulls (Test Płaskości)
  // // Jeśli średnia jest 0, a sigma 1, wydajność jest stała.
  // TF1 *f_gaus = new TF1("f_gaus", "gaus", -5, 5);
  // f_gaus->SetParameters(h_pulls->GetMaximum(), 0.0, 1.0); // Wstępne parametry
  // h_pulls->Fit(f_gaus, "R");

  // // Odczytaj średnią (p1) i odchylenie standardowe (p2) Gaussa
  // Double_t mean_pull = f_gaus->GetParameter(1);
  // Double_t sigma_pull = f_gaus->GetParameter(2);

  // Double_t mean_pull_err = f_gaus->GetParError(1);
  // Double_t sigma_pull_err = f_gaus->GetParError(2);

  tot_events = signal_num + bkg_tot;

  Double_t eff = CalculateEfficiency(signal_num, signal_tot);
  Double_t effAna = CalculateEfficiency(signal_num, signal_wo_err);
  Double_t purity = CalculatePurity(signal_num, tot_events);
  // Double_t metric = sigma_pull / (pow(eff, 2) * purity);

  // std::cout << "Pulls Fit Results: Mean = " << mean_pull << " +- " << mean_pull_err << ", Sigma = " << sigma_pull << " +- " << sigma_pull_err << std::endl;

  canvaEff->cd();

  efficiency->Draw();
  gPad->Update();

  efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.0);

  // Dodaj box z metrykami
  // TPaveText *metricsBox = new TPaveText(0.15, 0.15, 0.45, 0.35, "NDC");
  TPaveText *metricsBox = new TPaveText(0.15, 0.85, 0.45, 0.9, "NDC");
  metricsBox->SetFillColor(kWhite);
  metricsBox->SetBorderSize(1);
  metricsBox->SetTextAlign(12);
  metricsBox->SetTextSize(0.03);

  metricsBox->AddText(Form("Efficiency: %.2f%%", effAna * 100));
  // metricsBox->AddText(Form("Purity: %.2f%%", purity * 100));
  // metricsBox->AddText(Form("#sigma_{pull}: %.3f #pm %.3f", sigma_pull, sigma_pull_err));
  // // metricsBox->AddText(Form("Metric: %.4f", metric));
  // metricsBox->AddText(Form("#varepsilon_{fit}: (%.3f #pm %.3f)%%", eff_fit * 100, eff_fit_error * 100));

  metricsBox->Draw();

  canvaEff->SaveAs(folderPath + "/efficiency_delta_t.png");

  std::cout << "Signal events: " << signal_num << std::endl;
  std::cout << "Signal total events: " << signal_tot << std::endl;
  std::cout << "Background events: " << bkg_tot << std::endl;

  std::cout << "Total Efficiency signal: " << 100 * eff << " % (" << signal_num << "/" << signal_tot << ")" << std::endl;

  std::cout << "Analysis Efficiency signal: " << 100 * effAna << " % (" << signal_num << "/" << signal_wo_err << ")" << std::endl;

  std::cout << "Purity: " << 100 * purity << " % (" << signal_num << "/" << tot_events << ")" << std::endl;

  // std::cout << "Metryka: " << sigma_pull / (pow(eff, 2) * purity) << std::endl;

  std::cout << "Topology eff signal: " << 100 * signal_wo_err / (Float_t)signal_tot << " % (" << signal_wo_err << "/" << signal_tot << ")" << std::endl;
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
void signal_vs_bcg_v2::WidthOfCorrelatedHist(Double_t *means, Double_t *sigmas, Double_t *vLong, Double_t *vTransv)
{
  pca->MakePrincipals();
  const TVectorD *principalValues = pca->GetEigenValues();
  const TMatrixD *principalVectors = pca->GetEigenVectors();
  const TVectorD *dataMean = pca->GetMeanValues();
  const TVectorD *dataSigma = pca->GetSigmas();

  means[0] = (*dataMean)[0];              // Mean x
  means[1] = (*dataMean)[1];              // Mean y
  vLong[0] = (*principalVectors)(0, 0);   // Longitudinal vector
  vLong[1] = (*principalVectors)(1, 0);   // Longitudinal vector
  vTransv[0] = (*principalVectors)(0, 1); // Transverse vector
  vTransv[1] = (*principalVectors)(1, 1); // Transverse vector

  Double_t sigmaTransv = (*dataSigma)(1); //(*principalValues)[1] * (pow(vTransv[0], 2) * pow((*dataSigma)[0], 2) + pow(vTransv[1], 2) * pow((*dataSigma)[1], 2));

  Double_t sigmaLong = (*dataSigma)(0); //(*principalValues)[0] * (pow(vLong[0], 2) * pow((*dataSigma)[0], 2) + pow(vLong[1], 2) * pow((*dataSigma)[1], 2));

  sigmas[0] = TMath::Sqrt(sigmaTransv); // Long Sigma
  sigmas[1] = TMath::Sqrt(sigmaLong);   // Transv Sigma

  std::cout << "PCA Results:" << std::endl;
  std::cout << "Mean X: " << means[0] << ", Mean Y: " << means[1] << std::endl;
  std::cout << "Sigma Long: " << sigmas[0] << ", Sigma Transv: " << sigmas[1] << std::endl;
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