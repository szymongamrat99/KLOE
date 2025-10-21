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

KLOE::pm00 Obj;

Int_t nbins = 121;

void MC_fit_comparison::Begin(TTree * /*tree*/)
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
    canvas2D[histName.first] = new TCanvas(Form("c_%s", histName.first.Data()), Form("Canvas for %s", histName.first.Data()), 750, 750);
  }

  // Create histograms
  for (const auto &histName : KH::varNames)
  {

    TString nameRec = Form("h_rec_%s", histName.Data());
    TString nameFit = Form("h_fit_%s", histName.Data());

    histsReconstructed[histName] = KH::CreateHist1D(histName, "signal", nameRec);
    histsFittedSignal[histName] = KH::CreateHist1D(histName, "signal", nameFit);
  }

  // Create histograms 2D
  for (const auto &histName : KH::histConfigs2D)
  {
    TString nameRec = Form("h_rec2D_%s", histName.first.Data());
    TString nameFit = Form("h_fit2D_%s", histName.first.Data());

    hists2DReconstructed[histName.first] = KH::CreateHist2D(histName.first, nameRec);
    hists2DFittedSignal[histName.first] = KH::CreateHist2D(histName.first, nameFit);
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

  std::vector<Float_t> KchboostMom = {Kchboost[0], Kchboost[1], Kchboost[2], Kchboost[3]},
                       KnereclorMom = {*Bpx - Kchboost[0], *Bpy - Kchboost[1], *Bpz - Kchboost[2], *Broots - Kchboost[3]},
                       KchPos = {Kchboost[6], Kchboost[7], Kchboost[8]},
                       KnePos = {Knerec[6], Knerec[7], Knerec[8]},
                       ipPos = {ip[0], ip[1], ip[2]};

  KLOE::KaonProperTimes timesBoostLor = Obj.CalculateKaonProperTimes(KchboostMom,
                                                                     KchPos,
                                                                     KnereclorMom,
                                                                     KnePos,
                                                                     ipPos);

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
          radiusLimit = 1,
          zdistLimit = 0.6;

  Bool_t isInsideFiducialVolume = (radius00 < radiusLimit) && (zdist00 < zdistLimit) &&
                                  (radiuspm < radiusLimit) && (zdistpm < zdistLimit);

  Float_t T0Omega = 0;

  if(abs(Omega1ErrTmp) > abs(Omega2ErrTmp))
    T0Omega = Omega2MassTmp - (trk1[3] + trk2[3]) - pi0Omega2[5];
  else
    T0Omega = Omega1MassTmp - (trk1[3] + trk2[3]) - pi0Omega1[5];

  Float_t deltaTfit = *KaonChTimeCMSignalFit - *KaonNeTimeCMSignalFit,
          deltaT = timesBoostLor.deltaTimeCM,
          deltaTMC = *KaonChTimeCMMC - *KaonNeTimeCMMC;

  Float_t deltaPhi = *PhivSmeared1 - *PhivSmeared2;

  if (*mctruth == 0 || *mctruth == -1 || *mctruth == 1)
    signal_tot_err++;

  if (*mctruth == 0 || *mctruth == 1)
    signal_tot++;

  if (*mctruth == 1 && *Chi2SignalKinFit < 30.)// && isInsideFiducialVolume)
  {
    if (*goodClustersTriKinFitSize < 4)
      numberOfAtLeastOneBad++;
    if (*goodClustersTriKinFitSize >= 4)
      numberOfAllGood++;

    signal_num++;

    // Fill histograms for reconstructed variables
    histsReconstructed["px_Kch"]->Fill(Kchboost[0] - Kchmc[0]);
    histsReconstructed["py_Kch"]->Fill(Kchboost[1] - Kchmc[1]);
    histsReconstructed["pz_Kch"]->Fill(Kchboost[2] - Kchmc[2]);
    histsReconstructed["Energy_Kch"]->Fill(Kchboost[3] - Kchmc[3]);
    histsReconstructed["mass_Kch"]->Fill(Kchrec[5] - PhysicsConstants::mK0);

    histsReconstructed["px_Kne"]->Fill(Knerec[0] - Knemc[0]);
    histsReconstructed["py_Kne"]->Fill(Knerec[1] - Knemc[1]);
    histsReconstructed["pz_Kne"]->Fill(Knerec[2] - Knemc[2]);
    histsReconstructed["Energy_Kne"]->Fill(Knerec[3] - Knemc[3]);
    histsReconstructed["mass_Kne"]->Fill(*minv4gam - PhysicsConstants::mK0);

    histsReconstructed["mass_pi01"]->Fill(pi01[5] - PhysicsConstants::mPi0);
    histsReconstructed["mass_pi02"]->Fill(pi02[5] - PhysicsConstants::mPi0);

    if(abs(Omega1ErrTmp) > abs(Omega2ErrTmp))
      histsReconstructed["mass_omega"]->Fill(Omega2MassTmp - PhysicsConstants::mOmega);
    else
      histsReconstructed["mass_omega"]->Fill(Omega1MassTmp - PhysicsConstants::mOmega);

    if(abs(Omega1ErrTmp) > abs(Omega2ErrTmp))
      histsReconstructed["mass_omega_rec"]->Fill(Omega2MassTmp - PhysicsConstants::mOmega);
    else
      histsReconstructed["mass_omega_rec"]->Fill(Omega1MassTmp - PhysicsConstants::mOmega);

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

    histsReconstructed["T0Omega"]->Fill(T0Omega);

    // histsReconstructed["pathKne"]->Fill(pathKchFit);

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

    histsFittedSignal["mass_omega"]->Fill(omegaFit[5] - PhysicsConstants::mOmega);

    histsFittedSignal["chi2_signalKinFit"]->Fill(*Chi2SignalKinFit);
    histsFittedSignal["chi2_trilaterationKinFit"]->Fill(*Chi2TriKinFit);
    histsFittedSignal["prob_signal"]->Fill(TMath::Prob(*Chi2SignalKinFit, 10));

    histsFittedSignal["vtxNeu_x_Fit"]->Fill(KnereclorFit[6] - Knemc[6]);
    histsFittedSignal["vtxNeu_y_Fit"]->Fill(KnereclorFit[7] - Knemc[7]);
    histsFittedSignal["vtxNeu_z_Fit"]->Fill(KnereclorFit[8] - Knemc[8]);

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

    // histsFittedSignal["goodClusNumTriKinFit"]->Fill(*goodClustersTriKinFitSize);

    // histsFittedSignal["pathKne"]->Fill(pathKneFit);
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

  for (const auto &histName : KH::varNames)
  {
    canvas[histName]->cd();
    canvas[histName]->SetLogy(0);

    histsReconstructed[histName]->SetTitle("");

    histsReconstructed[histName]->GetYaxis()->SetRangeUser(0, std::max(histsReconstructed[histName]->GetMaximum(), histsFittedSignal[histName]->GetMaximum()) * 1.2);
    histsReconstructed[histName]->SetLineColor(kBlue);

    Bool_t fitcond = histName == "curv1" || histName == "phiv1" || histName == "cotv1" ||
                     histName == "curv2" || histName == "phiv2" || histName == "cotv2" ||
                     histName == "vtxNeu_x" || histName == "vtxNeu_y" || histName == "vtxNeu_z" || histName == "delta_t" || histName == "phi_vtx_x" || histName == "phi_vtx_y" || histName == "phi_vtx_z";
    ;

    if (fitcond)
    {
      Bool_t fitSuccess = false;

      // ✅ SPECJALNE TRAKTOWANIE DLA delta_t
      if (histName == "delta_t")
      {
        std::cout << "\n=== Fitting delta_t with TF1 (proven method) ===" << std::endl;

        // Użyj starej, sprawdzonej metody TF1 dla delta_t
        fitter->UseRooFit(false); // TF1 - działa lepiej dla delta_t
        fitter->SetVerbose(false);
        fitSuccess = fitter->FitHistogram(histsFittedSignal[histName]);
      }
      else
      {
        // ✅ Dla innych zmiennych użyj automatycznych parametrów
        fitter->UseRooFit(false); // Możesz zmienić na false żeby użyć TF1
        fitter->SetVerbose(false);
        fitSuccess = fitter->FitHistogram(histsFittedSignal[histName]);
      }

      if (fitSuccess)
      {
        // Pobierz wyniki fitu
        const KLOE::TripleGaussFitter::FitResult &result = fitter->GetLastResults();

        // Wypisz wyniki do konsoli
        std::cout << "\n=== Fit results for " << histName << " ===" << std::endl;
        std::cout << "Status: " << result.status << ", Converged: " << result.converged << std::endl;
        std::cout << "Chi2/NDF = " << result.chi2 << "/" << result.ndf
                  << " = " << fitter->GetChi2NDF() << std::endl;
        std::cout << "Combined mean = " << result.combinedMean
                  << " ± " << result.combinedMeanErr << std::endl;
        std::cout << "Combined sigma = " << result.combinedSigma
                  << " ± " << result.combinedSigmaErr << std::endl;

        // Wypisz parametry poszczególnych gaussów
        for (int i = 0; i < 3; i++)
        {
          int idx = i * 3;
          std::cout << "  Gauss " << (i + 1) << ": A=" << result.parameters[idx]
                    << ", μ=" << result.parameters[idx + 1]
                    << ", σ=" << result.parameters[idx + 2] << std::endl;
        }
      }
      else
      {
        std::cerr << "WARNING: Fit failed for " << histName << std::endl;
      }
    }

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
      legend->AddEntry(histsFittedSignal[histName], "Fitted Data", "l");

      // Sprawdź czy fit się udał przed użyciem wyników
      if (fitter->IsLastFitSuccessful())
      {
        const KLOE::TripleGaussFitter::FitResult &r = fitter->GetLastResults();
        legend->AddEntry((TObject *)0, Form("#chi^{2}/NDF = %.2f", fitter->GetChi2NDF()), "");
        legend->AddEntry((TObject *)0, Form("#mu = %.3f #pm %.3f", r.combinedMean, r.combinedMeanErr), "");
        legend->AddEntry((TObject *)0, Form("#sigma = %.3f #pm %.3f", r.combinedSigma, r.combinedSigmaErr), "");

        // Narysuj fit z komponentami
        fitter->DrawFitOnCurrentPad(true, false); // drawComponents=true, showStats=false
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

    canvas[histName]->Update();

    canvas[histName]->SaveAs(Form("img/%s_comparison.png", histName.Data()));
  }

  std::cout << "How many reconstructed events had all good clusters? " << numberOfAllGood << std::endl;
  std::cout << "How many reconstructed events had at least one bad cluster? " << numberOfAtLeastOneBad << std::endl;
  std::cout << "Percentage of fully good events: " << (Float_t)numberOfAllGood / (numberOfAllGood + numberOfAtLeastOneBad) * 100 << " %" << std::endl
            << std::endl;

  std::cout << "Signal events without errors: " << signal_tot << " over " << signal_tot_err << " (" << (Float_t)signal_tot / signal_tot_err * 100 << " %) --> Efficiency of preselection" << std::endl;
  std::cout << "Signal events after cuts: " << signal_num << " over " << signal_tot << " (" << (Float_t)signal_num / signal_tot * 100 << " %) --> Efficiency of analysis" << std::endl;
}