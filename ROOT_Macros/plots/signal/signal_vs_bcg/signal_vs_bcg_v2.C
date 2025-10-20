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
#include <triple_gaus.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <TGraphAsymmErrors.h>

#include <interf_function.h>
#include <TEfficiency.h>

namespace KH = KLOE::Histograms;

Int_t signal_num = 0, signal_tot = 0, tot_events = 0, bkg_tot = 0;

std::map<TString, TCanvas *>
    canvas;

std::map<TString, std::map<TString, TCanvas *>>
    canvas2D;

std::map<TString, std::map<TString, TH1 *>>
    histsReconstructed,
    histsFittedSignal;

std::map<TString, std::map<TString, TH2 *>>
    hists2DReconstructed,
    hists2DFittedSignal;

TCanvas *canvaEff;

TH1 *deltaTSignalTot;

KLOE::TripleGaussFitter *fitter;

KLOE::pm00 Obj;

Int_t nbins = 121;

void signal_vs_bcg_v2::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  KLOE::setGlobalStyle();

  fitter = new KLOE::TripleGaussFitter();

  canvaEff = new TCanvas("Efficiency", "Efficiency", 800, 800);
  deltaTSignalTot = new TH1D("EfficiencyHistTot", "Efficiency Histogram Total; #Deltat [#tau_{S}]; Efficiency [-]", nbins, -30, 30);

  // Create canvases
  for (const auto &histName : KH::varNames)
  {
    canvas[histName] = new TCanvas(Form("c_%s", histName.Data()),
                                   Form("Canvas for %s", histName.Data()), 750, 750);
  }

  // Create canvases
  for (const auto &histName : KH::histConfigs2D)
  {
    for (const auto &name : KLOE::channName)
    {
      canvas2D[histName.first][name.second] = new TCanvas(Form("c_%s_%s", histName.first.Data(), name.second.Data()), Form("Canvas for %s (%s)", histName.first.Data(), name.second.Data()), 750, 750);
    }
  }

  // Create histograms
  for (const auto &histName : KH::varNames)
  {
    for (const auto &name : KLOE::channName)
    {
      TString nameRec = Form("h_rec_%s_%s", histName.Data(), name.second.Data());
      TString nameFit = Form("h_fit_%s_%s", histName.Data(), name.second.Data());

      histsReconstructed[histName][name.second] = KH::CreateHist1D(histName, name.second, nameRec);
      histsFittedSignal[histName][name.second] = KH::CreateHist1D(histName, name.second, nameFit);
    }
  }

  // Create histograms 2D
  for (const auto &histName : KH::histConfigs2D)
  {
    for (const auto &name : KLOE::channName)
    {
      TString nameRec = Form("h_rec2D_%s_%s", histName.first.Data(), name.second.Data());
      TString nameFit = Form("h_fit2D_%s_%s", histName.first.Data(), name.second.Data());

      hists2DReconstructed[histName.first][name.second] = KH::CreateHist2D(histName.first, nameRec);
      hists2DFittedSignal[histName.first][name.second] = KH::CreateHist2D(histName.first, nameFit);
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

  fReader.SetLocalEntry(entry);

  Float_t Knerec[9];

  Knerec[0] = gammaMomTriangle1[0] + gammaMomTriangle2[0] +
              gammaMomTriangle3[0] + gammaMomTriangle4[0];
  Knerec[1] = gammaMomTriangle1[1] + gammaMomTriangle2[1] +
              gammaMomTriangle3[1] + gammaMomTriangle4[1];
  Knerec[2] = gammaMomTriangle1[2] + gammaMomTriangle2[2] +
              gammaMomTriangle3[2] + gammaMomTriangle4[2];
  Knerec[3] = gammaMomTriangle1[3] + gammaMomTriangle2[3] +
              gammaMomTriangle3[3] + gammaMomTriangle4[3];

  Knerec[4] = sqrt(pow(Knerec[0], 2) + pow(Knerec[1], 2) + pow(Knerec[2], 2));

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
          RKneFit = sqrt(pow(KnerecFit[6] - ipFit[0], 2) +
                         pow(KnerecFit[7] - ipFit[1], 2)),
          tKneFit = KnereclorFit[9] / 0.0895,
          vKneMC = 0, // PhysicsConstants::cVel * Knemc[4] / Knemc[3],
      vKne = PhysicsConstants::cVel * Knerec[4] / Knerec[3],
          pathKne = sqrt(pow(Knerec[6] - ip[0], 2) +
                         pow(Knerec[7] - ip[1], 2) +
                         pow(Knerec[8] - ip[2], 2)),
          tKne = pathKne / (vKne * 0.0895),
          pathKchMC = sqrt(pow(Kchmc[6] - ipmc[0], 2) +
                           pow(Kchmc[7] - ipmc[1], 2) +
                           pow(Kchmc[8] - ipmc[2], 2)),
          pathKneMC = sqrt(pow(Knemc[6] - ipmc[0], 2) +
                           pow(Knemc[7] - ipmc[1], 2) +
                           pow(Knemc[8] - ipmc[2], 2));

  std::vector<Float_t> kaonMom1Fit = {KchboostFit[0],
                                      KchboostFit[1],
                                      KchboostFit[2],
                                      KchboostFit[3]},
                       kaonMom2Fit = {KnereclorFit[0],
                                      KnereclorFit[1],
                                      KnereclorFit[2],
                                      KnereclorFit[3]},
                       kaonPos1Fit = {KchboostFit[6],
                                      KchboostFit[7],
                                      KchboostFit[8],
                                      KchboostFit[9]},
                       kaonPos2Fit = {KnereclorFit[6],
                                      KnereclorFit[7],
                                      KnereclorFit[8],
                                      KnereclorFit[9]},
                       ipFitVec = {ipFit[0], ipFit[1], ipFit[2]};

  KLOE::KaonProperTimes propTimesFit = Obj.CalculateKaonProperTimes(kaonMom1Fit,
                                                                    kaonPos1Fit,
                                                                    kaonMom2Fit,
                                                                    kaonPos2Fit,
                                                                    ipFitVec);

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

  Float_t deltaTfit = propTimesFit.deltaTimeCM,
          deltaT = *KaonChTimeCMBoostLor - *KaonNeTimeCMBoostLor,
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

  if (*mctruth == 1)
    weight = 1.0;//interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

  if ((*mctruth == 1 || *mctruth == -1 || *mctruth == 0) && *mcflag == 1)
    signal_tot++;

  if (*mctruth == 1)
    deltaTSignalTot->Fill(deltaTfit, weight);

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
      deltaPhiFit;

  if (*PhivSmeared1 > *PhivSmeared2)
    deltaPhi = *PhivSmeared1 - *PhivSmeared2;
  else if (*PhivSmeared1 < *PhivSmeared2)
    deltaPhi = *PhivSmeared2 - *PhivSmeared1;

  if (ParamSignalFit[24] > ParamSignalFit[27])
    deltaPhiFit = ParamSignalFit[24] - ParamSignalFit[27];
  else if (ParamSignalFit[24] < ParamSignalFit[27])
    deltaPhiFit = ParamSignalFit[27] - ParamSignalFit[24];

  // Analiza Simony ciecie na phi bad
  const Double_t SLOPE = -10.0 / 9.0;

  Bool_t condGeneral = (deltaTfit - deltaTMC)<2.0,
                                              condLowerLimit = (deltaTfit - deltaTMC)> SLOPE *
                       deltaTMC,
         condUpperLimit = (deltaTfit - deltaTMC) < (SLOPE * (deltaTMC - 2.0)),
         condTotalSimona = condGeneral && condLowerLimit && condUpperLimit;
  ///////////////////////////////////////////////////////////////////////////////
  // Analiza Simony cięcie na 3 sigma mas
  Bool_t condMassKch = abs(Kchrec[5] - 0.937 - PhysicsConstants::mK0) < 3 * 1.583,
         condMassKne = abs(*minv4gam + 13.106 - PhysicsConstants::mK0) < 3 * 28.991,
         condMassPi01 = abs(pi01Fit[5] + 0.255 - PhysicsConstants::mPi0) < 3 * 4.383,
         condMassPi02 = abs(pi02Fit[5] + 0.063 - PhysicsConstants::mPi0) < 3 * 4.061;

  ///////////////////////////////////////////////////////////////////////////////

  if (*mctruth >= 0 /*&& combinedMassPi0Fit < 10. && *Chi2SignalKinFit < 30.*/ /*&& *TrcSum > -1 && abs(Kchrec[5] - PhysicsConstants::mK0) < 1.2 && abs(*minv4gam - PhysicsConstants::mK0) < 76. &&  *Qmiss < 3.75 && openingAngleCharged > acosCutAngle && abs(deltaPhi - 3.09) > 2 * 0.087 && condMassKch && condMassKne && condMassPi01 && condMassPi02 && pathKchMC < 30.0 && pathKneMC < 30.0*/)
  {
    if ((*mctruth == 1) && *mcflag == 1)
      signal_num++;

    if (*mctruth > 1 && *mcflag == 1)
      bkg_tot++;

    if ((*mcflag == 1 && *mctruth >= 1) || *mcflag == 0)
    {
      // Fill histograms for reconstructed variables
      histsReconstructed["mass_Kch"][KLOE::channName.at(*mctruth)]->Fill(Kchrec[5] - PhysicsConstants::mK0, weight);

      histsReconstructed["mass_Kne"][KLOE::channName.at(*mctruth)]->Fill(*minv4gam - PhysicsConstants::mK0, weight);

      histsReconstructed["mass_pi01"][KLOE::channName.at(*mctruth)]->Fill(pi01[5] - PhysicsConstants::mPi0, weight);
      histsReconstructed["mass_pi02"][KLOE::channName.at(*mctruth)]->Fill(pi02[5] - PhysicsConstants::mPi0, weight);

      histsReconstructed["time_neutral_MC"][KLOE::channName.at(*mctruth)]->Fill(*TrcSum, weight);

      histsReconstructed["combined_mass_pi0"][KLOE::channName.at(*mctruth)]->Fill(combinedMassPi0Fit, weight);

      // Fitted signal variables
      histsFittedSignal["mass_Kch"][KLOE::channName.at(*mctruth)]->Fill(Kchrec[5] - PhysicsConstants::mOmega, weight);

      histsFittedSignal["mass_Kne"][KLOE::channName.at(*mctruth)]->Fill(*minv4gam - PhysicsConstants::mK0, weight);

      histsFittedSignal["mass_pi01"][KLOE::channName.at(*mctruth)]->Fill(pi01Fit[5] - PhysicsConstants::mPi0, weight);
      histsFittedSignal["mass_pi02"][KLOE::channName.at(*mctruth)]->Fill(pi02Fit[5] - PhysicsConstants::mPi0, weight);

      histsFittedSignal["chi2_signalKinFit"][KLOE::channName.at(*mctruth)]->Fill(*Chi2SignalKinFit, weight);
      histsFittedSignal["chi2_trilaterationKinFit"][KLOE::channName.at(*mctruth)]->Fill(*Chi2TriKinFit, weight);
      histsFittedSignal["prob_signal"][KLOE::channName.at(*mctruth)]->Fill(TMath::Prob(*Chi2SignalKinFit, 10), weight);

      histsFittedSignal["combined_mass_pi0"][KLOE::channName.at(*mctruth)]->Fill(combinedMassPi0Fit, weight);

      for (Int_t i = 0; i < 36; i++)
      {
        histsFittedSignal["pull" + std::to_string(i + 1)][KLOE::channName.at(*mctruth)]->Fill(pullsSignalFit[i], weight);
      }

      histsFittedSignal["time_neutral_MC"][KLOE::channName.at(*mctruth)]->Fill(*TrcSum, weight);

      histsFittedSignal["openingAngleCharged"][KLOE::channName.at(*mctruth)]->Fill(openingAngleCharged, weight);
      histsFittedSignal["openingAngleNeutral"][KLOE::channName.at(*mctruth)]->Fill(openingAngleNeutral, weight);

      histsFittedSignal["Qmiss"][KLOE::channName.at(*mctruth)]->Fill(*Qmiss, weight);

      histsFittedSignal["delta_t"][KLOE::channName.at(*mctruth)]->Fill(deltaTfit, weight);

      if (condTotalSimona)
      {
        histsFittedSignal["deltaPhiv"][KLOE::channName.at(*mctruth)]->Fill(deltaPhi, weight);

        histsFittedSignal["deltaPhivFit"][KLOE::channName.at(*mctruth)]->Fill(deltaPhiFit, weight);
      }

      histsFittedSignal["nev"][KLOE::channName.at(*mctruth)]->Fill(*nev, weight);
      histsFittedSignal["nrun"][KLOE::channName.at(*mctruth)]->Fill(*nrun, weight);

      if (abs(pathKchMC - pathKneMC) > 5)
        histsFittedSignal["TransvRadius"][KLOE::channName.at(*mctruth)]->Fill(pathKneFit, weight);

      hists2DFittedSignal["delta_t_fit_vs_delta_t_mc"][KLOE::channName.at(*mctruth)]->Fill(deltaTMC, deltaTfit - deltaTMC, weight);
      hists2DFittedSignal["delta_t_vs_delta_t_mc"][KLOE::channName.at(*mctruth)]->Fill(deltaTMC, deltaT - deltaTMC, weight);

      hists2DFittedSignal["chi2_signalKinFit_vs_chi2_trilaterationKinFit"][KLOE::channName.at(*mctruth)]->Fill(*Chi2SignalKinFit, *Chi2TriKinFit, weight);
    }
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
  // the results graphically or save the results to file.

  // MC sum

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::map<TString, std::map<TString, Int_t>> integrals;
  std::map<TString, Int_t> maxNum;
  std::map<TString, TString> maxChannel;

  for (const auto &histName : KH::varNames)
  {
    maxNum[histName] = 0;
    for (const auto &channelType : KLOE::channName)
    {
      if (histsFittedSignal[histName][channelType.second]->GetEntries() > 0.)
        histsFittedSignal[histName][channelType.second]->Scale(histsFittedSignal[histName][channelType.second]->GetEntries() / histsFittedSignal[histName][channelType.second]->Integral());

      integrals[histName][channelType.second] = histsFittedSignal[histName][channelType.second]->GetMaximum();

      if (integrals[histName][channelType.second] > maxNum[histName])
      {
        maxNum[histName] = integrals[histName][channelType.second];
        maxChannel[histName] = channelType.second;
      }
    }
  }

  for (const auto &histName : KH::varNames)
  {
    for (const auto &channelType : KLOE::channName)
    {
      if (channelType.second != "Data" && channelType.second != "MC sum")
      {
        histsReconstructed[histName][KLOE::channName.at(8)]->Add(histsReconstructed[histName][channelType.second]);
        histsFittedSignal[histName][KLOE::channName.at(8)]->Add(histsFittedSignal[histName][channelType.second]);
      }
    }

    Double_t scalefactorReconstructed = 1.0;
    Double_t scalefactorFittedSignal = 1.0;

    if (histsReconstructed[histName]["MC sum"]->GetEntries() > 0 && histsReconstructed[histName]["Data"]->GetEntries() > 0)
    {
      scalefactorReconstructed = histsReconstructed[histName]["Data"]->GetEntries() / histsReconstructed[histName]["MC sum"]->GetEntries();

      histsReconstructed[histName]["MC sum"]->Scale(scalefactorReconstructed);
    }

    if (histsFittedSignal[histName]["MC sum"]->GetEntries() > 0 && histsFittedSignal[histName]["Data"]->GetEntries() > 0)
    {
      scalefactorFittedSignal = histsFittedSignal[histName]["Data"]->GetEntries() / histsFittedSignal[histName]["MC sum"]->GetEntries();
      histsFittedSignal[histName]["MC sum"]->Scale(scalefactorFittedSignal);
    }

    for (const auto &channelType : KLOE::channName)
    {
      if (channelType.second != "Data" && channelType.second != "MC sum")
      {
        histsReconstructed[histName][channelType.second]->Scale(scalefactorReconstructed);
        histsFittedSignal[histName][channelType.second]->Scale(scalefactorFittedSignal);
      }
    }
  }

  for (const auto &histName : KH::varNames)
  {
    // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
    TLegend *legend = new TLegend(0.6, 0.65, 0.9, 0.9, "", "NDC");

    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.03);

    canvas[histName]->cd();
    canvas[histName]->SetLogy(0);
    for (const auto &channelType : KLOE::channName)
    {
      histsFittedSignal[histName][channelType.second]->SetLineColor(KLOE::channColor.at(channelType.second));

      histsFittedSignal[histName][channelType.second]->SetTitle("");

      histsFittedSignal[histName][channelType.second]->GetYaxis()->SetRangeUser(0, maxNum[histName] * 1.2);

      if (channelType.second == "Data")
      {
        histsFittedSignal[histName][channelType.second]->Draw("PE");
        legend->AddEntry(histsFittedSignal[histName][channelType.second], channelType.second, "pe");
      }
      else
      {
        histsFittedSignal[histName][channelType.second]->Draw("HIST SAME");
        legend->AddEntry(histsFittedSignal[histName][channelType.second], channelType.second, "l");

        if (channelType.second == "Signal" && (histName == "mass_Kch" || histName == "mass_Kne" || histName == "mass_pi01" || histName == "mass_pi02"))
        {

          fitter->FitHistogram(histsFittedSignal[histName][channelType.second]);
          fitter->DrawFitOnCurrentPad();
        }
      }
    }
    gPad->Update();
    legend->Draw();

    // Sprawdź czy JAKIKOLWIEK histogram ma wpisy
    Bool_t hasEntries = kFALSE;

    // Sprawdź wszystkie kanały
    for (const auto &name : KLOE::channName)
    {
      if (histsFittedSignal[histName][name.second]->GetEntries() > 0)
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

    canvas[histName]->SaveAs(Form("img/%s_comparison.png", histName.Data()));
  }

  for (const auto &histName : KH::histConfigs2D)
  {
    for (const auto &channelType : KLOE::channName)
    {
      canvas2D[histName.first][channelType.second]->cd();
      canvas2D[histName.first][channelType.second]->SetLogz(1);

      hists2DFittedSignal[histName.first][channelType.second]->Draw("COLZ");

      // Sprawdź czy JAKIKOLWIEK histogram ma wpisy
      Bool_t hasEntries = kFALSE;

      // Sprawdź wszystkie kanały

      if (hists2DFittedSignal[histName.first][channelType.second]->GetEntries() > 0)
      {
        hasEntries = kTRUE;
      }

      // Jeśli żaden histogram nie ma wpisów - pomiń ten canvas
      if (!hasEntries)
      {
        continue;
      }

      canvas2D[histName.first][channelType.second]->SaveAs(Form("img/%s_%s_2D.png", histName.first.Data(), channelType.second.Data()));
    }
  }

  deltaTSignalTot->Scale(deltaTSignalTot->GetEntries() / deltaTSignalTot->Integral());

  TEfficiency *efficiency = new TEfficiency(*histsFittedSignal["delta_t"]["Signal"], *deltaTSignalTot);

  efficiency->SetStatisticOption(TEfficiency::kBUniform);

  TF1 *f_const = new TF1("f_const", "pol0", -300., 300.);

  f_const->SetParameter(0, 0.1);

  TFitResultPtr fitResult = efficiency->Fit(f_const, "SR");

  if (fitResult->IsValid())
  {
    std::cout << "Fit results for efficiency:" << std::endl;
    fitResult->Print("V");

    fitResult->GetErrors();
  }
  else
  {
    std::cout << "Fit failed!" << std::endl;
  }

  // Przy założeniu, że eff i f_const są już zdefiniowane i fit został wykonany

  // 1. Stworzenie histogramu dla Pulls
  TH1F *h_pulls = new TH1F("h_pulls", "Pulls Distribution; (Eps_i - Eps_fit) / sigma_i; Entries", 50, -5.0, 5.0);

  // 2. Pobranie dopasowanej stałej i jej błędu
  Double_t eff_fit = f_const->GetParameter(0);

  // 3. Iteracja po binach TEfficiency (iterujemy po histogramie mianownika)
  TH1 *h_total = (TH1 *)efficiency->GetTotalHistogram();
  for (int i = 1; i <= h_total->GetNbinsX(); ++i)
  {
    // Sprawdzenie, czy bin jest w zakresie fitowania (opcjonalnie, ale zalecane)
    Double_t x_center = h_total->GetBinCenter(i);

    // Wydobycie wartości wydajności i jej błędu z TEfficiency
    Double_t epsilon_i = efficiency->GetEfficiency(i);
    // Użyj większego z błędów dla konserwatywnej analizy (upper/lower error)
    Double_t error_i = efficiency->GetEfficiencyErrorUp(i) > efficiency->GetEfficiencyErrorLow(i) ? efficiency->GetEfficiencyErrorUp(i) : efficiency->GetEfficiencyErrorLow(i);

    // Tylko dla binów z niezerową statystyką i błędem
    if (efficiency->GetTotalHistogram()->GetBinContent(i) > 0 && error_i > 0)
    {
      Double_t pull = (epsilon_i - eff_fit) / error_i;
      h_pulls->Fill(pull);
    }
  }

  // 4. Dopasowanie Gaussa do rozkładu Pulls (Test Płaskości)
  // Jeśli średnia jest 0, a sigma 1, wydajność jest stała.
  TF1 *f_gaus = new TF1("f_gaus", "gaus", -5, 5);
  f_gaus->SetParameters(h_pulls->GetMaximum(), 0.0, 1.0); // Wstępne parametry
  h_pulls->Fit(f_gaus, "R");

  // Odczytaj średnią (p1) i odchylenie standardowe (p2) Gaussa
  Double_t mean_pull = f_gaus->GetParameter(1);
  Double_t sigma_pull = f_gaus->GetParameter(2);

  Double_t mean_pull_err = f_gaus->GetParError(1);
  Double_t sigma_pull_err = f_gaus->GetParError(2);

  std::cout << "Pulls Fit Results: Mean = " << mean_pull << " +- " << mean_pull_err << ", Sigma = " << sigma_pull << " +- " << sigma_pull_err << std::endl;

  canvaEff->cd();

  efficiency->Draw();
  gPad->Update();

  efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 0.5);

  canvaEff->SaveAs("img/efficiency.png");

  tot_events = signal_num + bkg_tot;

  std::cout << "Signal events: " << signal_num << std::endl;
  std::cout << "Signal total events: " << signal_tot << std::endl;
  std::cout << "Background events: " << bkg_tot << std::endl;

  std::cout << "Efficiency signal: " << 100 * signal_num / (Float_t)signal_tot << " % (" << signal_num << "/" << signal_tot << ")" << std::endl;

  std::cout << "Purity: " << 100 * signal_num / (Float_t)tot_events << " % (" << signal_num << "/" << tot_events << ")" << std::endl;

  Double_t eff = (Double_t)signal_num / (Double_t)signal_tot,
           purity = (Double_t)signal_num / (Double_t)tot_events;

  std::cout << "Metryka: " << sigma_pull / (pow(eff, 2) * purity) << std::endl;
}