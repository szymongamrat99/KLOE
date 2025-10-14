#define signal_vs_bcg_cxx
// The class definition in signal_vs_bcg.h has been generated automatically
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
// root> T->Process("signal_vs_bcg.C")
// root> T->Process("signal_vs_bcg.C","some options")
// root> T->Process("signal_vs_bcg.C+")
//

#include "signal_vs_bcg.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
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

Int_t signal_num = 0, signal_tot = 0, tot_events = 0, bkg_tot = 0;

std::map<Int_t, TString> channelTypes = {
    {0, "Data"},
    {1, "Signal"},
    {2, "Regeneration"},
    {3, "Omega-pi0"},
    {4, "3pi0"},
    {5, "Semileptonic"},
    {6, "Other"},
    {7, "pi+pi-pi+pi-"}};

std::map<TString, Color_t> channelColors = {
    {"Data", kBlack},
    {"Signal", kRed},
    {"Regeneration", kGreen},
    {"Omega-pi0", kViolet},
    {"3pi0", kCyan},
    {"Semileptonic", kBlue},
    {"Other", kGreen - 1},
    {"pi+pi-pi+pi-", kYellow}};

std::vector<TString> histList = {"px_Pi1", "py_Pi1", "pz_Pi1", "Energy_Pi1",
                                 "px_Pi2", "py_Pi2", "pz_Pi2", "Energy_Pi2",
                                 "px_Kch", "py_Kch", "pz_Kch", "Energy_Kch",
                                 "px_Kne", "py_Kne", "pz_Kne", "Energy_Kne",
                                 "px_phi", "py_phi", "pz_phi", "Energy_phi",
                                 "mass_Kch", "mass_Kne", "mass_phi", "chi2_signalKinFit",
                                 "chi2_trilaterationKinFit", "curv1", "phiv1", "cotv1",
                                 "curv2", "phiv2", "cotv2", "vtxNeu_x", "vtxNeu_y", "vtxNeu_z",
                                 "vtxNeu_x_Fit", "vtxNeu_y_Fit", "vtxNeu_z_Fit",
                                 "mass_pi01", "mass_pi02",
                                 "time_neutral_MC", "prob_signal", "delta_t",
                                 "combined_mass_pi0",
                                 "pull1", "pull2", "pull3", "pull4", "pull5", "phi_vtx_x", "phi_vtx_y", "phi_vtx_z", "vKne", "openingAngleCharged", "openingAngleNeutral",
                                 "Qmiss", "deltaPhiv", "deltaPhivFit", "nev", "nrun"}; // List of histograms to be created

std::map<TString, std::vector<TString>> histTitles = {
    {"px_Pi1", {"p_{x} - p_{x}^{MC} [MeV/c]"}},
    {"py_Pi1", {"p_{y} - p_{y}^{MC} [MeV/c]"}},
    {"pz_Pi1", {"p_{z} - p_{z}^{MC} [MeV/c]"}},
    {"Energy_Pi1", {"E - E^{MC} [MeV]"}},
    {"px_Pi2", {"p_{x} - p_{x}^{MC} [MeV/c]"}},
    {"py_Pi2", {"p_{y} - p_{y}^{MC} [MeV/c]"}},
    {"pz_Pi2", {"p_{z} - p_{z}^{MC} [MeV/c]"}},
    {"Energy_Pi2", {"E - E^{MC} [MeV]"}},
    {"px_Kch", {"p_{x} - p_{x}^{MC} [MeV/c]"}},
    {"py_Kch", {"p_{y} - p_{y}^{MC} [MeV/c]"}},
    {"pz_Kch", {"p_{z} - p_{z}^{MC} [MeV/c]"}},
    {"Energy_Kch", {"E - E^{MC} [MeV]"}},
    {"px_Kne", {"p_{x} - p_{x}^{MC} [MeV/c]"}},
    {"py_Kne", {"p_{y} - p_{y}^{MC} [MeV/c]"}},
    {"pz_Kne", {"p_{z} - p_{z}^{MC} [MeV/c]"}},
    {"Energy_Kne", {"E - E^{MC} [MeV]"}},
    {"px_phi", {"p_{x} - p_{x}^{MC} [MeV/c]"}},
    {"py_phi", {"p_{y} - p_{y}^{MC} [MeV/c]"}},
    {"pz_phi", {"p_{z} - p_{z}^{MC} [MeV/c]"}},
    {"Energy_phi", {"E - E^{MC} [MeV]"}},
    {"mass_Kch", {"m_{#pi^{+}#pi^{-}} - m_{K^{0}} [MeV]"}},
    {"mass_Kne", {"m^{inv}_{4#gamma} - m_{K^{0}} [MeV]"}},
    {"mass_phi", {"m^{inv}_{#phi} - m_{#phi} [MeV]"}},
    {"chi2_signalKinFit", {"#chi^{2} of signal kinematic fit"}},
    {"chi2_trilaterationKinFit", {"#chi^{2} of trilateration kinematic fit"}},
    {"curv1", {"Curvature_{1} - Curvature_{1}^{MC} [1/cm]"}},
    {"phiv1", {"#phi_{1} - #phi_{1}^{MC} [rad]"}},
    {"cotv1", {"cot(#theta_{1}) - cot(#theta_{1}^{MC})"}},
    {"curv2", {"Curvature_{2} - Curvature_{2}^{MC} [1/cm]"}},
    {"phiv2", {"#phi_{2} - #phi_{2}^{MC} [rad]"}},
    {"cotv2", {"cot(#theta_{2}) - cot(#theta_{2}^{MC})"}},
    {"vtxNeu_x", {"x_{neu} - x_{neu}^{MC} [cm]"}},
    {"vtxNeu_y", {"y_{neu} - y_{neu}^{MC} [cm]"}},
    {"vtxNeu_z", {"z_{neu} - z_{neu}^{MC} [cm]"}},
    {"vtxNeu_x_Fit", {"x_{neu} - x_{neu}^{MC} [cm]"}},
    {"vtxNeu_y_Fit", {"y_{neu} - y_{neu}^{MC} [cm]"}},
    {"vtxNeu_z_Fit", {"z_{neu} - z_{neu}^{MC} [cm]"}},
    {"mass_pi01", {"m_{#gamma#gamma} - m_{#pi^{0}} [MeV]"}},
    {"mass_pi02", {"m_{#gamma#gamma} - m_{#pi^{0}} [MeV]"}},
    {"time_neutral_MC", {"#sum_{i} T_{cl,i} - t_{K_{ne}} - t_{#gamma,i} [ns]"}},
    {"prob_signal", {"Probability of signal"}},
    {"delta_t", {"#Deltat - #Deltat^{MC} [#tau_{S}]"}},
    {"combined_mass_pi0", {"#sqrt{(m_{#gamma#gamma,1} - m_{#pi^{0}})^2 + (m_{#gamma#gamma,2} - m_{#pi^{0}})^2} [MeV/c^{2}]"}},
    {"pull1", {"Pull_{1} [MeV]"}},
    {"pull2", {"Pull_{2} [MeV]"}},
    {"pull3", {"Pull_{3} [MeV]"}},
    {"pull4", {"Pull_{4} [MeV]"}},
    {"pull5", {"Pull_{5} [MeV]"}},
    {"phi_vtx_x", {"#phi_{vtx,x} - #phi_{vtx,x}^{MC} [cm]"}},
    {"phi_vtx_y", {"#phi_{vtx,y} - #phi_{vtx,y}^{MC} [cm]"}},
    {"phi_vtx_z", {"#phi_{vtx,z} - #phi_{vtx,z}^{MC} [cm]"}},
    {"vKne", {"v_{K_{ne}} - v_{K_{ne}}^{MC} [cm/ns]"}},
    {"openingAngleCharged", {"Opening angle between #pi^{+} and #pi^{-} [deg]"}},
    {"openingAngleNeutral", {"Opening angle between neutral pions [deg]"}},
    {"Qmiss", {"Missing energy Q_{miss} [MeV]"}},
    {"deltaPhiv", {"#Delta#phi_{+-} [rad]"}},
    {"deltaPhivFit", {"#Delta#phi_{+-}^{fit} [rad]"}},
    {"nev", {"Event number"}},
    {"nrun", {"Run number"}}};

std::map<TString, std::vector<Double_t>>
    histLimits = {{"px_Pi1", {-200, 200}},
                  {"py_Pi1", {-200, 200}},
                  {"pz_Pi1", {-250, 250}},
                  {"Energy_Pi1", {-20, 20}},
                  {"px_Pi2", {-200, 200}},
                  {"py_Pi2", {-200, 200}},
                  {"pz_Pi2", {-250, 250}},
                  {"Energy_Pi2", {-20, 20}},
                  {"px_Kch", {-200, 200}},
                  {"py_Kch", {-200, 200}},
                  {"pz_Kch", {-250, 250}},
                  {"Energy_Kch", {-20, 20}},
                  {"px_Kne", {-200, 200}},
                  {"py_Kne", {-200, 200}},
                  {"pz_Kne", {-250, 250}},
                  {"Energy_Kne", {-20, 20}},
                  {"px_phi", {-200, 200}},
                  {"py_phi", {-200, 200}},
                  {"pz_phi", {-250, 250}},
                  {"Energy_phi", {-20, 20}},
                  {"mass_Kch", {-20, 20}},
                  {"mass_Kne", {-200, 200}},
                  {"mass_phi", {-20, 20}},
                  {"chi2_signalKinFit", {0, 10}},
                  {"chi2_trilaterationKinFit", {0, 100}},
                  {"curv1", {-0.5, 0.5}},
                  {"phiv1", {-0.1, 0.1}},
                  {"cotv1", {-0.1, 0.1}},
                  {"curv2", {-0.5, 0.5}},
                  {"phiv2", {-0.1, 0.1}},
                  {"cotv2", {-0.1, 0.1}},
                  {"vtxNeu_x", {-10, 10}},
                  {"vtxNeu_y", {-10, 10}},
                  {"vtxNeu_z", {-10, 10}},
                  {"vtxNeu_x_Fit", {-10, 10}},
                  {"vtxNeu_y_Fit", {-10, 10}},
                  {"vtxNeu_z_Fit", {-10, 10}},
                  {"mass_pi01", {-100, 100}},
                  {"mass_pi02", {-100, 100}},
                  {"time_neutral_MC", {-5, 2}},
                  {"prob_signal", {0, 1}},
                  {"delta_t", {-20, 20}},
                  {"combined_mass_pi0", {-100, 100}},
                  {"pull1", {-5, 5}},
                  {"pull2", {-5, 5}},
                  {"pull3", {-5, 5}},
                  {"pull4", {-5, 5}},
                  {"pull5", {-5, 5}},
                  {"phi_vtx_x", {-1, 1}},
                  {"phi_vtx_y", {-0.2, 0.2}},
                  {"phi_vtx_z", {-10, 10}},
                  {"vKne", {-2, 2}},
                  {"openingAngleCharged", {0, 182}},
                  {"openingAngleNeutral", {150, 190}},
                  {"Qmiss", {0., 300}},
                  {"deltaPhiv", {2, 4}},
                  {"deltaPhivFit", {2, 4}},
                  {"nev", {0, 60e3}},
                  {"nrun", {30000, 43000}}};

std::map<TString, TCanvas *> canvas;

std::map<TString, std::map<TString, TH1 *>>
    histsReconstructed,
    histsFittedSignal;

TCanvas *canvaEff;

TH1 *deltaTSignalTot;

KLOE::TripleGaussFitter *fitter;

Int_t nbins = 81;

void signal_vs_bcg::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  fitter = new KLOE::TripleGaussFitter();

  canvaEff = new TCanvas("Efficiency", "Efficiency", 800, 800);
  deltaTSignalTot = new TH1D("EfficiencyHistTot", "Efficiency Histogram Total; #Delta t; Efficiency", nbins, -20, 20);

  // Create canvases
  for (const auto &histName : histList)
  {
    canvas[histName] = new TCanvas(Form("c_%s", histName.Data()),
                                   Form("Canvas for %s", histName.Data()), 750, 750);
  }

  // Create histograms
  for (const auto &histName : histList)
  {

    for (Int_t i = 0; i < channelTypes.size(); i++)
    {
      histsReconstructed[histName][channelTypes[i]] = new TH1F(Form("h_reconstructed_%s_%s", histName.Data(), channelTypes[i].Data()), Form("Reconstructed %s; %s; Counts", histName.Data(), histName.Data()), nbins, histLimits[histName][0], histLimits[histName][1]);
      histsFittedSignal[histName][channelTypes[i]] = new TH1F(Form("h_fittedSignal_%s_%s", histName.Data(), channelTypes[i].Data()), Form("Fitted Signal %s; %s; Counts", histName.Data(), histName.Data()), nbins, histLimits[histName][0], histLimits[histName][1]);
    }
  }
}

void signal_vs_bcg::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

Bool_t signal_vs_bcg::Process(Long64_t entry)
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

  Float_t vKchFit = cVel * KchboostFit[4] / KchboostFit[3],
          pathKchFit = sqrt(pow(KchboostFit[6] - ipFit[0], 2) +
                            pow(KchboostFit[7] - ipFit[1], 2) +
                            pow(KchboostFit[8] - ipFit[2], 2)),
          tKchFit = pathKchFit / (vKchFit * 0.0895),
          vKneFit = cVel * KnereclorFit[4] / KnereclorFit[3],
          pathKneFit = sqrt(pow(ParamSignalFit[33] - ipFit[0], 2) +
                            pow(ParamSignalFit[34] - ipFit[1], 2) +
                            pow(ParamSignalFit[35] - ipFit[2], 2)),
          tKneFit = pathKneFit / (vKneFit * 0.0895),
          vKneMC = cVel * Knemc[4] / Knemc[3],
          vKne = cVel * Knerec[4] / Knerec[3],
          pathKne = sqrt(pow(KneTriangle[6] - ip[0], 2) +
                         pow(KneTriangle[7] - ip[1], 2) +
                         pow(KneTriangle[8] - ip[2], 2)),
          tKne = pathKne / (vKne * 0.0895);

  Float_t photon1path = sqrt(pow(ParamSignalFit[0] - KnerecFit[6], 2) +
                             pow(ParamSignalFit[1] - KnerecFit[7], 2) +
                             pow(ParamSignalFit[2] - KnerecFit[8], 2)),
          photon2path = sqrt(pow(ParamSignalFit[5] - KnerecFit[6], 2) +
                             pow(ParamSignalFit[6] - KnerecFit[7], 2) +
                             pow(ParamSignalFit[7] - KnerecFit[8], 2)),
          photon3path = sqrt(pow(ParamSignalFit[10] - KnerecFit[6], 2) +
                             pow(ParamSignalFit[11] - KnerecFit[7], 2) +
                             pow(ParamSignalFit[12] - KnerecFit[8], 2)),
          photon4path = sqrt(pow(ParamSignalFit[15] - KnerecFit[6], 2) +
                             pow(ParamSignalFit[16] - KnerecFit[7], 2) +
                             pow(ParamSignalFit[17] - KnerecFit[8], 2));

  Float_t trc1Fit = ParamSignalFit[3] - photon1path / cVel - tKneFit * 0.0895,
          trc2Fit = ParamSignalFit[8] - photon2path / cVel - tKneFit * 0.0895,
          trc3Fit = ParamSignalFit[13] - photon3path / cVel - tKneFit * 0.0895,
          trc4Fit = ParamSignalFit[18] - photon4path / cVel - tKneFit * 0.0895,
          TrcSumFit = trc1Fit + trc2Fit + trc3Fit + trc4Fit;

  Float_t deltaTfit = tKchFit - tKneFit,
          deltaT = *KaonChTimeLAB - tKne,
          deltaTMC = *KaonChTimeLABMC - *KaonNeTimeLABMC;

  Float_t combinedMassPi0Fit = sqrt(pow(pi01Fit[5] - mPi0, 2) +
                                    pow(pi02Fit[5] - mPi0, 2)),
          combinedMassPi0 = sqrt(pow(pi01[5] - mPi0, 2) +
                                 pow(pi02[5] - mPi0, 2));

  Float_t kaonChPath = sqrt(pow(Kchrec[6] - ip[0], 2) +
                            pow(Kchrec[7] - ip[1], 2) +
                            pow(Kchrec[8] - ip[2], 2)),
          kaonChMom[3] = {(Kchrec[6] - ip[0]) / kaonChPath * Kchboost[4],
                          (Kchrec[7] - ip[1]) / kaonChPath * Kchboost[4],
                          (Kchrec[8] - ip[2]) / kaonChPath * Kchboost[4]};

  TVector3 boostNeutralKaon(-KneTriangle[0] / KneTriangle[3],
                            -KneTriangle[1] / KneTriangle[3],
                            -KneTriangle[2] / KneTriangle[3]),
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
    // weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

  if ((*mctruth == 1 || *mctruth == 0 || *mctruth == -1) && *mcflag == 1)
    signal_tot++;

  if (*mctruth == 1)
    deltaTSignalTot->Fill(tKchFit - tKneFit, weight);

  // if (*mctruth == 5 && *mcflag == 1)
  //   signal_tot++;

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

  if (*mctruth > 0 /*&& *Chi2SignalKinFit < 30. && abs(Kchrec[5] - mK0) < 1.2 && abs(*minv4gam - mK0) < 90. && combinedMassPi0Fit < 12. && *Qmiss < 3.75*/)
  {
    if ((*mctruth == 1 || *mctruth == 0) && *mcflag == 1)
      signal_num++;

    if (*mctruth > 1 && *mcflag == 1)
      bkg_tot++;

    if (*mcflag == 1 && *mctruth >= 1)
    {

      // Fill histograms for reconstructed variables
      histsReconstructed["mass_Kch"][channelTypes[*mctruth]]->Fill(Kchrec[5] - mK0, weight);

      histsReconstructed["mass_Kne"][channelTypes[*mctruth]]->Fill(*minv4gam - mK0, weight);

      histsReconstructed["mass_pi01"][channelTypes[*mctruth]]->Fill(pi01[5] - mPi0, weight);
      histsReconstructed["mass_pi02"][channelTypes[*mctruth]]->Fill(pi02[5] - mPi0, weight);

      histsReconstructed["time_neutral_MC"][channelTypes[*mctruth]]->Fill(*TrcSum, weight);

      histsReconstructed["combined_mass_pi0"][channelTypes[*mctruth]]->Fill(combinedMassPi0, weight);

      // Fitted signal variables
      histsFittedSignal["mass_Kch"][channelTypes[*mctruth]]->Fill(KchrecFit[5] - mK0, weight);

      histsFittedSignal["mass_Kne"][channelTypes[*mctruth]]->Fill(KnerecFit[5] - mK0, weight);

      histsFittedSignal["mass_pi01"][channelTypes[*mctruth]]->Fill(pi01Fit[5] - mPi0, weight);
      histsFittedSignal["mass_pi02"][channelTypes[*mctruth]]->Fill(pi02Fit[5] - mPi0, weight);

      histsFittedSignal["chi2_signalKinFit"][channelTypes[*mctruth]]->Fill(*Chi2SignalKinFit / 10., weight);
      histsFittedSignal["chi2_trilaterationKinFit"][channelTypes[*mctruth]]->Fill(*Chi2TriKinFit, weight);
      histsFittedSignal["prob_signal"][channelTypes[*mctruth]]->Fill(TMath::Prob(*Chi2SignalKinFit, 10), weight);

      histsFittedSignal["combined_mass_pi0"][channelTypes[*mctruth]]->Fill(combinedMassPi0Fit, weight);

      histsFittedSignal["pull1"][channelTypes[*mctruth]]->Fill(pullsSignalFit[0], weight);
      histsFittedSignal["pull2"][channelTypes[*mctruth]]->Fill(pullsSignalFit[1], weight);
      histsFittedSignal["pull3"][channelTypes[*mctruth]]->Fill(pullsSignalFit[2], weight);
      histsFittedSignal["pull4"][channelTypes[*mctruth]]->Fill(pullsSignalFit[3], weight);
      histsFittedSignal["pull5"][channelTypes[*mctruth]]->Fill(pullsSignalFit[4], weight);

      histsFittedSignal["time_neutral_MC"][channelTypes[*mctruth]]->Fill(TrcSumFit, weight);

      histsFittedSignal["openingAngleCharged"][channelTypes[*mctruth]]->Fill(openingAngleCharged, weight);
      histsFittedSignal["openingAngleNeutral"][channelTypes[*mctruth]]->Fill(openingAngleNeutral, weight);

      histsFittedSignal["Qmiss"][channelTypes[*mctruth]]->Fill(*Qmiss, weight);

      histsFittedSignal["delta_t"][channelTypes[*mctruth]]->Fill(tKchFit - tKneFit, weight);

      histsFittedSignal["deltaPhiv"][channelTypes[*mctruth]]->Fill(*PhivSmeared1 - *PhivSmeared2, weight);

      histsFittedSignal["deltaPhivFit"][channelTypes[*mctruth]]->Fill(ParamSignalFit[24] - ParamSignalFit[27], weight);

      histsFittedSignal["nev"][channelTypes[*mctruth]]->Fill(*nev, weight);
      histsFittedSignal["nrun"][channelTypes[*mctruth]]->Fill(*nrun, weight);
    }
  }

  return kTRUE;
}

void signal_vs_bcg::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void signal_vs_bcg::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::map<TString, std::map<TString, Int_t>> integrals;
  std::map<TString, Int_t> maxNum;
  std::map<TString, TString> maxChannel;

  for (const auto &histName : histList)
  {
    maxNum[histName] = 0;
    for (const auto &channelType : channelTypes)
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

  for (const auto &histName : histList)
  {
    canvas[histName]->cd();
    canvas[histName]->SetLogy(0);
    for (const auto &channelType : channelTypes)
    {
      histsFittedSignal[histName][channelType.second]->SetLineColor(channelColors[channelType.second]);

      if (channelType.second == maxChannel[histName])
      {
        histsFittedSignal[histName][channelType.second]->SetTitle("");

        histsFittedSignal[histName][channelType.second]->GetXaxis()->SetTitle(histTitles[histName][0]);

        histsFittedSignal[histName][maxChannel[histName]]->GetYaxis()->SetRangeUser(0, maxNum[histName] * 1.2);
        histsFittedSignal[histName][channelType.second]->Draw("HIST");
      }
      else
      {
        histsFittedSignal[histName][channelType.second]->Draw("HIST SAME");
      }

      // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
      // TLegend *legend = new TLegend(0.15, 0.75, 0.4, 0.9, "", "NDC");

      // legend->SetBorderSize(1);
      // legend->SetFillColor(kWhite);
      // legend->SetFillStyle(1001);
      // legend->SetTextSize(0.03);
      // legend->Draw();
    }
    gPad->Update();
    canvas[histName]->SaveAs(Form("img/%s_comparison.png", histName.Data()));
  }

  deltaTSignalTot->Scale(deltaTSignalTot->GetEntries() / deltaTSignalTot->Integral());

  std::cout << histsFittedSignal["delta_t"]["Signal"]->Integral() << " " << deltaTSignalTot->Integral() << std::endl;

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

  efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.0);

  canvaEff->SaveAs("img/efficiency.png");

  tot_events = signal_num + bkg_tot;

  std::cout << "Signal events: " << signal_num << std::endl;
  std::cout << "Background events: " << bkg_tot << std::endl;

  std::cout << "Efficiency signal: " << 100 * signal_num / (Float_t)signal_tot << " % (" << signal_num << "/" << signal_tot << ")" << std::endl;

  std::cout << "Purity: " << 100 * signal_num / (Float_t)tot_events << " % (" << signal_num << "/" << tot_events << ")" << std::endl;

  Double_t eff = (Double_t)signal_num / (Double_t)signal_tot,
           purity = (Double_t)signal_num / (Double_t)tot_events;

  std::cout << "Metryka: " << sigma_pull / (pow(eff, 2) * purity) << std::endl;
}