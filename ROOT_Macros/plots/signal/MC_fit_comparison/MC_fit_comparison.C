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
#include <triple_gaus.h>
#include <TPaveText.h>
#include <TLegend.h>

#include <TFitResult.h>
#include <TFitResultPtr.h>

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
                                 "pull1", "pull2", "pull3", "pull4", "pull5", "phi_vtx_x", "phi_vtx_y", "phi_vtx_z", "vKne",
                                 "goodClusNumTriKinFit", "pathKne", "pathKch"}; // List of histograms to be created

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
    {"goodClusNumTriKinFit", {"Number of good clusters used in trilateration kinematic fit"}},
    {"pathKne", {"Path length [cm]"}},
    {"pathKch", {"Path length of K_{ch} [cm]"}}};

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
                  {"chi2_signalKinFit", {0, 50}},
                  {"chi2_trilaterationKinFit", {0, 50}},
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
                  {"delta_t", {-15, 15}},
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
                  {"goodClusNumTriKinFit", {-1, 6}},
                  {"pathKne", {0, 50.}},
                  {"pathKch", {0, 50.}}};

std::map<TString, TCanvas *> canvas;

std::map<TString, TH1 *>
    histsReconstructed,
    histsFittedSignal;

KLOE::TripleGaussFitter *fitter;

void MC_fit_comparison::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  fitter = new KLOE::TripleGaussFitter();

  // Create canvases
  for (const auto &histName : histList)
  {
    canvas[histName] = new TCanvas(Form("c_%s", histName.Data()),
                                   Form("Canvas for %s", histName.Data()), 750, 750);
  }

  // Create histograms
  for (const auto &histName : histList)
  {

    histsReconstructed[histName] = new TH1F(Form("h_reconstructed_%s", histName.Data()),
                                            Form("Reconstructed %s; %s; Counts", histName.Data(), histName.Data()),
                                            100, histLimits[histName][0], histLimits[histName][1]);
    histsFittedSignal[histName] = new TH1F(Form("h_fittedSignal_%s", histName.Data()),
                                           Form("Fitted Signal %s; %s; Counts", histName.Data(), histName.Data()),
                                           100, histLimits[histName][0], histLimits[histName][1]);
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
          tKchFit = KchboostFit[9] / (0.0895),
          vKneFit = PhysicsConstants::cVel * KnereclorFit[4] / KnereclorFit[3],
          pathKneFit = sqrt(pow(KnereclorFit[6] - ipFit[0], 2) +
                            pow(KnereclorFit[7] - ipFit[1], 2) +
                            pow(KnereclorFit[8] - ipFit[2], 2)),
          tKneFit = KnereclorFit[9] / (0.0895),
          vKneMC = PhysicsConstants::cVel * Knemc[4] / Knemc[3],
          vKne = PhysicsConstants::cVel * Knerec[4] / Knerec[3],
          pathKne = sqrt(pow(Knerec[6] - ip[0], 2) +
                         pow(Knerec[7] - ip[1], 2) +
                         pow(Knerec[8] - ip[2], 2)),
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

  Float_t deltaTfit = tKchFit - tKneFit,
          deltaT = *KaonChTimeLABBoostLor - tKne,
          deltaTMC = *KaonChTimeLABMC - *KaonNeTimeLABMC;

  Float_t deltaPhi = *PhivSmeared1 - *PhivSmeared2;

  if (*mctruth == 1)
  {
    if (*goodClustersTriKinFitSize < 4)
      numberOfAtLeastOneBad++;
    if (*goodClustersTriKinFitSize >= 4)
      numberOfAllGood++;

    // Fill histograms for reconstructed variables
    histsReconstructed["px_Kch"]->Fill(Kchrec[0] - Kchmc[0]);
    histsReconstructed["py_Kch"]->Fill(Kchrec[1] - Kchmc[1]);
    histsReconstructed["pz_Kch"]->Fill(Kchrec[2] - Kchmc[2]);
    histsReconstructed["Energy_Kch"]->Fill(Kchrec[3] - Kchmc[3]);
    histsReconstructed["mass_Kch"]->Fill(Kchrec[5] - PhysicsConstants::mK0);

    histsReconstructed["px_Kne"]->Fill(Knerec[0] - Knemc[0]);
    histsReconstructed["py_Kne"]->Fill(Knerec[1] - Knemc[1]);
    histsReconstructed["pz_Kne"]->Fill(Knerec[2] - Knemc[2]);
    histsReconstructed["Energy_Kne"]->Fill(Knerec[3] - Knemc[3]);
    histsReconstructed["mass_Kne"]->Fill(*minv4gam - PhysicsConstants::mK0);

    histsReconstructed["mass_pi01"]->Fill(pi01[5] - PhysicsConstants::mPi0);
    histsReconstructed["mass_pi02"]->Fill(pi02[5] - PhysicsConstants::mPi0);

    histsReconstructed["vtxNeu_x"]->Fill(Knerec[6] - Knemc[6]);
    histsReconstructed["vtxNeu_y"]->Fill(Knerec[7] - Knemc[7]);
    histsReconstructed["vtxNeu_z"]->Fill(Knerec[8] - Knemc[8]);

    histsReconstructed["phi_vtx_x"]->Fill(*Bx - ipmc[0]);
    histsReconstructed["phi_vtx_y"]->Fill(*By - ipmc[1]);
    histsReconstructed["phi_vtx_z"]->Fill(*Bz - ipmc[2]);

    histsReconstructed["vtxNeu_x_Fit"]->Fill(Knerec[6] - Knemc[6]);
    histsReconstructed["vtxNeu_y_Fit"]->Fill(Knerec[7] - Knemc[7]);
    histsReconstructed["vtxNeu_z_Fit"]->Fill(Knerec[8] - Knemc[8]);

    histsReconstructed["vKne"]->Fill(vKne - vKneMC);

    histsReconstructed["time_neutral_MC"]->Fill(*TrcSum);

    histsReconstructed["delta_t"]->Fill(deltaT - deltaTMC);

    histsReconstructed["combined_mass_pi0"]->Fill(combinedMassPi0);

    histsReconstructed["pathKne"]->Fill(pathKchFit);

    // Decide which reconstructed track corresponds to which MC particle

    Double_t error1 = sqrt(pow(*CurvSmeared1 - CurvMC[0], 2) +
                           pow(*PhivSmeared1 - PhivMC[0], 2) +
                           pow(*CotvSmeared1 - CotvMC[0], 2)),
             error2 = sqrt(pow(*CurvSmeared1 - CurvMC[1], 2) +
                           pow(*PhivSmeared1 - PhivMC[1], 2) +
                           pow(*CotvSmeared1 - CotvMC[1], 2));

    if (error1 < error2)
    {
      histsReconstructed["curv1"]->Fill(*CurvSmeared1 - CurvMC[0]);
      histsReconstructed["phiv1"]->Fill(*PhivSmeared1 - PhivMC[0]);
      histsReconstructed["cotv1"]->Fill(*CotvSmeared1 - CotvMC[0]);

      histsReconstructed["curv2"]->Fill(*CurvSmeared2 - CurvMC[1]);
      histsReconstructed["phiv2"]->Fill(*PhivSmeared2 - PhivMC[1]);
      histsReconstructed["cotv2"]->Fill(*CotvSmeared2 - CotvMC[1]);
    }
    else
    {
      histsReconstructed["curv1"]->Fill(*CurvSmeared1 - CurvMC[1]);
      histsReconstructed["phiv1"]->Fill(*PhivSmeared1 - PhivMC[1]);
      histsReconstructed["cotv1"]->Fill(*CotvSmeared1 - CotvMC[1]);

      histsReconstructed["curv2"]->Fill(*CurvSmeared2 - CurvMC[0]);
      histsReconstructed["phiv2"]->Fill(*PhivSmeared2 - PhivMC[0]);
      histsReconstructed["cotv2"]->Fill(*CotvSmeared2 - CotvMC[0]);
    }

    // Fitted signal variables
    histsFittedSignal["px_Kch"]->Fill(KchrecFit[0] - Kchmc[0]);
    histsFittedSignal["py_Kch"]->Fill(KchrecFit[1] - Kchmc[1]);
    histsFittedSignal["pz_Kch"]->Fill(KchrecFit[2] - Kchmc[2]);
    histsFittedSignal["Energy_Kch"]->Fill(KchrecFit[3] - Kchmc[3]);
    histsFittedSignal["mass_Kch"]->Fill(KchrecFit[5] - PhysicsConstants::mK0);

    histsFittedSignal["px_Kne"]->Fill(KnerecFit[0] - Knemc[0]);
    histsFittedSignal["py_Kne"]->Fill(KnerecFit[1] - Knemc[1]);
    histsFittedSignal["pz_Kne"]->Fill(KnerecFit[2] - Knemc[2]);
    histsFittedSignal["Energy_Kne"]->Fill(KnerecFit[3] - Knemc[3]);
    histsFittedSignal["mass_Kne"]->Fill(KnerecFit[5] - PhysicsConstants::mK0);

    histsFittedSignal["mass_pi01"]->Fill(pi01Fit[5] - PhysicsConstants::mPi0);
    histsFittedSignal["mass_pi02"]->Fill(pi02Fit[5] - PhysicsConstants::mPi0);

    histsFittedSignal["chi2_signalKinFit"]->Fill(*Chi2SignalKinFit);
    histsFittedSignal["chi2_trilaterationKinFit"]->Fill(*Chi2TriKinFit);
    histsFittedSignal["prob_signal"]->Fill(TMath::Prob(*Chi2SignalKinFit, 10));

    histsFittedSignal["vtxNeu_x_Fit"]->Fill(KnerecFit[6] - Knemc[6]);
    histsFittedSignal["vtxNeu_y_Fit"]->Fill(KnerecFit[7] - Knemc[7]);
    histsFittedSignal["vtxNeu_z_Fit"]->Fill(KnerecFit[8] - Knemc[8]);

    histsFittedSignal["phi_vtx_x"]->Fill(ipFit[0] - ipmc[0]);
    histsFittedSignal["phi_vtx_y"]->Fill(ipFit[1] - ipmc[1]);
    histsFittedSignal["phi_vtx_z"]->Fill(ipFit[2] - ipmc[2]);

    histsFittedSignal["vKne"]->Fill(vKneFit - vKneMC);

    histsFittedSignal["delta_t"]->Fill(deltaTfit - deltaTMC);

    histsFittedSignal["combined_mass_pi0"]->Fill(combinedMassPi0Fit);

    histsFittedSignal["pull1"]->Fill(pullsSignalFit[0]);
    histsFittedSignal["pull2"]->Fill(pullsSignalFit[1]);
    histsFittedSignal["pull3"]->Fill(pullsSignalFit[2]);
    histsFittedSignal["pull4"]->Fill(pullsSignalFit[3]);
    histsFittedSignal["pull5"]->Fill(pullsSignalFit[4]);

    histsFittedSignal["time_neutral_MC"]->Fill(TrcSumFit);

    histsFittedSignal["goodClusNumTriKinFit"]->Fill(*goodClustersTriKinFitSize);

    histsFittedSignal["pathKne"]->Fill(pathKneFit);
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

  for (const auto &histName : histList)
  {
    canvas[histName]->cd();
    canvas[histName]->SetLogy(0);

    histsReconstructed[histName]->SetTitle("");

    histsReconstructed[histName]->GetYaxis()->SetRangeUser(0, std::max(histsReconstructed[histName]->GetMaximum(), histsFittedSignal[histName]->GetMaximum()) * 1.2);
    histsReconstructed[histName]->SetLineColor(kBlue);

    Bool_t fitcond = histName == "curv1" || histName == "phiv1" || histName == "cotv1" ||
                     histName == "curv2" || histName == "phiv2" || histName == "cotv2" ||
                     histName == "vtxNeu_x" || histName == "vtxNeu_y" || histName == "vtxNeu_z" || histName == "delta_t";

    Bool_t fitcondGaus = histName == "phi_vtx_x" || histName == "phi_vtx_y" || histName == "phi_vtx_z";

    if (fitcond)
    {
      fitter->FitHistogram(histsFittedSignal[histName]);
      KLOE::TripleGaussFitter::FitResult result = fitter->GetLastResults();
    }

    if (fitcondGaus)
    {
      histsFittedSignal[histName]->Fit("gaus");
    }

    histsReconstructed[histName]->GetXaxis()->SetTitle(histTitles[histName][0]);
    histsReconstructed[histName]->GetYaxis()->SetRangeUser(0.0, 1.5 * std::max(histsReconstructed[histName]->GetMaximum(), histsFittedSignal[histName]->GetMaximum()));

    histsReconstructed[histName]->Draw();
    histsFittedSignal[histName]->SetLineColor(kRed);
    histsFittedSignal[histName]->Draw("SAME");

    // DODAJ WŁASNĄ LEGENDĘ W LEWYM GÓRNYM ROGU:
    TLegend *legend = new TLegend(0.15, 0.75, 0.4, 0.9, "", "NDC");

    if (!fitcond)
    {
      legend->AddEntry(histsReconstructed[histName], "Reconstructed", "l");
      legend->AddEntry(histsFittedSignal[histName], "Fitted Signal", "l");
    }
    else
    {
      legend->AddEntry(histsReconstructed[histName], "Reconstructed", "l");
      legend->AddEntry(histsFittedSignal[histName], "3x Gauss", "l");
      fitter->DrawFitOnCurrentPad();
    }

    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.03);
    legend->Draw();

    canvas[histName]->Update();

    canvas[histName]->SaveAs(Form("img/%s_comparison.png", histName.Data()));
  }

  std::cout << "How many reconstructed events had all good clusters? " << numberOfAllGood << std::endl;
  std::cout << "How many reconstructed events had at least one bad cluster? " << numberOfAtLeastOneBad << std::endl;
  std::cout << "Percentage of fully good events: " << (Float_t)numberOfAllGood / (numberOfAllGood + numberOfAtLeastOneBad) * 100 << " %" << std::endl;
}