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
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>

std::vector<TString> histList = {"px_Pi1", "py_Pi1", "pz_Pi1", "Energy_Pi1",
                                 "px_Pi2", "py_Pi2", "pz_Pi2", "Energy_Pi2",
                                 "px_Kch", "py_Kch", "pz_Kch", "Energy_Kch",
                                 "px_Kne", "py_Kne", "pz_Kne", "Energy_Kne",
                                 "px_phi", "py_phi", "pz_phi", "Energy_phi",
                                 "mass_Kch", "mass_Kne", "mass_phi", "chi2_signalKinFit",
                                 "chi2_trilaterationKinFit", "curv1", "phiv1", "cotv1",
                                 "curv2", "phiv2", "cotv2", "vtxNeu_x", "vtxNeu_y", "vtxNeu_z",
                                 "vtxNeu_x_Fit", "vtxNeu_y_Fit", "vtxNeu_z_Fit",
                                 "mass_pi01", "mass_pi02"}; // List of histograms to be created

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
                  {"mass_Kne", {-20, 20}},
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
                  {"mass_pi02", {-100, 100}}}; // Histogram limits

std::map<TString, TCanvas *> canvas;

std::map<TString, TH1 *>
    histsReconstructed,
    histsFittedSignal;

void MC_fit_comparison::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  std::cout << option << std::endl;

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

  if (*mctruth == 1)
  {
    // Fill histograms for reconstructed variables
    histsReconstructed["px_Kch"]->Fill(Kchrec[0] - Kchmc[0]);
    histsReconstructed["py_Kch"]->Fill(Kchrec[1] - Kchmc[1]);
    histsReconstructed["pz_Kch"]->Fill(Kchrec[2] - Kchmc[2]);
    histsReconstructed["Energy_Kch"]->Fill(Kchrec[3] - Kchmc[3]);
    histsReconstructed["mass_Kch"]->Fill(Kchrec[5] - mK0);

    histsReconstructed["px_Kne"]->Fill(KneTriangle[0] - Knemc[0]);
    histsReconstructed["py_Kne"]->Fill(KneTriangle[1] - Knemc[1]);
    histsReconstructed["pz_Kne"]->Fill(KneTriangle[2] - Knemc[2]);
    histsReconstructed["Energy_Kne"]->Fill(KneTriangle[3] - Knemc[3]);
    histsReconstructed["mass_Kne"]->Fill(*minv4gam - mK0);

    histsReconstructed["mass_pi01"]->Fill(pi01[5] - mPi0);
    histsReconstructed["mass_pi02"]->Fill(pi02[5] - mPi0);

    histsReconstructed["vtxNeu_x"]->Fill(Kchrec[6] - Kchmc[6]);
    histsReconstructed["vtxNeu_y"]->Fill(Kchrec[7] - Kchmc[7]);
    histsReconstructed["vtxNeu_z"]->Fill(Kchrec[8] - Kchmc[8]);

    histsReconstructed["vtxNeu_x_Fit"]->Fill(Kchrec[6] - Kchmc[6]);
    histsReconstructed["vtxNeu_y_Fit"]->Fill(Kchrec[7] - Kchmc[7]);
    histsReconstructed["vtxNeu_z_Fit"]->Fill(Kchrec[8] - Kchmc[8]);

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
    histsFittedSignal["px_Kch"]->Fill(Kchboost[0] - Kchmc[0]);
    histsFittedSignal["py_Kch"]->Fill(Kchboost[1] - Kchmc[1]);
    histsFittedSignal["pz_Kch"]->Fill(Kchboost[2] - Kchmc[2]);
    histsFittedSignal["Energy_Kch"]->Fill(Kchboost[3] - Kchmc[3]);
    histsFittedSignal["mass_Kch"]->Fill(Kchboost[5] - mK0);

    histsFittedSignal["px_Kne"]->Fill(KnerecFit[0] - Knemc[0]);
    histsFittedSignal["py_Kne"]->Fill(KnerecFit[1] - Knemc[1]);
    histsFittedSignal["pz_Kne"]->Fill(KnerecFit[2] - Knemc[2]);
    histsFittedSignal["Energy_Kne"]->Fill(KnerecFit[3] - Knemc[3]);
    histsFittedSignal["mass_Kne"]->Fill(KnerecFit[5] - mK0);

    histsFittedSignal["mass_pi01"]->Fill(pi01Fit[5] - mPi0);
    histsFittedSignal["mass_pi02"]->Fill(pi02Fit[5] - mPi0);

    histsFittedSignal["chi2_signalKinFit"]->Fill(*Chi2SignalKinFit);
    histsFittedSignal["chi2_trilaterationKinFit"]->Fill(*Chi2TriKinFit);

    histsFittedSignal["vtxNeu_x_Fit"]->Fill(KchrecFit[6] - Kchmc[6]);
    histsFittedSignal["vtxNeu_y_Fit"]->Fill(KchrecFit[7] - Kchmc[7]);
    histsFittedSignal["vtxNeu_z_Fit"]->Fill(KchrecFit[8] - Kchmc[8]);
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

  for (const auto &histName : histList)
  {
    canvas[histName]->cd();
    histsReconstructed[histName]->GetYaxis()->SetRangeUser(0, std::max(histsReconstructed[histName]->GetMaximum(), histsFittedSignal[histName]->GetMaximum()) * 1.2);
    histsReconstructed[histName]->SetLineColor(kBlue);

    if (histName == "curv1" || histName == "phiv1" || histName == "cotv1" ||
        histName == "curv2" || histName == "phiv2" || histName == "cotv2" ||
        histName == "vtxNeu_x" || histName == "vtxNeu_y" || histName == "vtxNeu_z")
    {

      histsReconstructed[histName]->Fit("gaus");
    }

    histsReconstructed[histName]->Draw();
    histsFittedSignal[histName]->SetLineColor(kRed);
    histsFittedSignal[histName]->Draw("HIST SAME");
    canvas[histName]->BuildLegend();
    canvas[histName]->Update();

    canvas[histName]->SaveAs(Form("img/%s_comparison.png", histName.Data()));
  }
}