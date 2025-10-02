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
                                 "mass_Kch", "mass_Kne", "mass_phi"};

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
      if (histName.Contains("px") || histName.Contains("py"))
      {

         histsReconstructed[histName] = new TH1F(Form("h_reconstructed_%s", histName.Data()),
                                                 Form("Reconstructed %s; %s; Counts", histName.Data(), histName.Data()),
                                                 100, -200, 200);
         histsFittedSignal[histName] = new TH1F(Form("h_fittedSignal_%s", histName.Data()),
                                                Form("Fitted Signal %s; %s; Counts", histName.Data(), histName.Data()),
                                                100, -200, 200);
      }

      if (histName.Contains("pz"))
      {

         histsReconstructed[histName] = new TH1F(Form("h_reconstructed_%s", histName.Data()),
                                                 Form("Reconstructed %s; %s; Counts", histName.Data(), histName.Data()),
                                                 100, -250, 250);
         histsFittedSignal[histName] = new TH1F(Form("h_fittedSignal_%s", histName.Data()),
                                                Form("Fitted Signal %s; %s; Counts", histName.Data(), histName.Data()),
                                                100, -250, 250);
      }

      if (histName.Contains("Energy"))
      {

         histsReconstructed[histName] = new TH1F(Form("h_reconstructed_%s", histName.Data()),
                                                 Form("Reconstructed %s; %s; Counts", histName.Data(), histName.Data()),
                                                 50, -20, 20);
         histsFittedSignal[histName] = new TH1F(Form("h_fittedSignal_%s", histName.Data()),
                                                Form("Fitted Signal %s; %s; Counts", histName.Data(), histName.Data()),
                                                50, -20, 20);
      }

      if (histName.Contains("mass"))
      {

         histsReconstructed[histName] = new TH1F(Form("h_reconstructed_%s", histName.Data()),
                                                 Form("Reconstructed %s; %s; Counts", histName.Data(), histName.Data()),
                                                 50, -20, 20);
         histsFittedSignal[histName] = new TH1F(Form("h_fittedSignal_%s", histName.Data()),
                                                Form("Fitted Signal %s; %s; Counts", histName.Data(), histName.Data()),
                                                50, -20, 20);
      }
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
      histsReconstructed["px_Kch"]->Fill(Kchboost[0] - Kchmc[0]);
      histsReconstructed["py_Kch"]->Fill(Kchboost[1] - Kchmc[1]);
      histsReconstructed["pz_Kch"]->Fill(Kchboost[2] - Kchmc[2]);
      histsReconstructed["Energy_Kch"]->Fill(Kchboost[3] - Kchmc[3]);
      histsReconstructed["mass_Kch"]->Fill(Kchrec[5] - mK0);

      histsReconstructed["px_Kne"]->Fill(KneTriangle[0] - Knemc[0]);
      histsReconstructed["py_Kne"]->Fill(KneTriangle[1] - Knemc[1]);
      histsReconstructed["pz_Kne"]->Fill(KneTriangle[2] - Knemc[2]);
      histsReconstructed["Energy_Kne"]->Fill(KneTriangle[3] - Knemc[3]);
      histsReconstructed["mass_Kne"]->Fill(*minv4gam - mK0);

      // Fitted signal variables
      histsFittedSignal["px_Kch"]->Fill(KchboostFit[0] - Kchboost[0]);
      histsFittedSignal["py_Kch"]->Fill(KchboostFit[1] - Kchboost[1]);
      histsFittedSignal["pz_Kch"]->Fill(KchboostFit[2] - Kchboost[2]);
      histsFittedSignal["Energy_Kch"]->Fill(KchboostFit[3] - Kchboost[3]);
      histsFittedSignal["mass_Kch"]->Fill(KchrecFit[5] - mK0);

      histsFittedSignal["px_Kne"]->Fill(KnereclorFit[0] - KneTriangle[0]);
      histsFittedSignal["py_Kne"]->Fill(KnereclorFit[1] - KneTriangle[1]);
      histsFittedSignal["pz_Kne"]->Fill(KnereclorFit[2] - KneTriangle[2]);
      histsFittedSignal["Energy_Kne"]->Fill(KnereclorFit[3] - KneTriangle[3]);
      histsFittedSignal["mass_Kne"]->Fill(KnerecFit[5] - mK0);
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
      histsReconstructed[histName]->SetLineColor(kBlue);
      histsReconstructed[histName]->Draw("HIST");
      histsFittedSignal[histName]->SetLineColor(kRed);
      histsFittedSignal[histName]->Draw("HIST SAME");
      canvas[histName]->BuildLegend();
      canvas[histName]->Update();

      canvas[histName]->SaveAs(Form("img/%s_comparison.png", histName.Data()));
   }
}