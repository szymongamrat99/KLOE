#define PlotsPiPi_cxx
// The class definition in PlotsPiPi.h has been generated automatically
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
// root> T->Process("PlotsPiPi.C")
// root> T->Process("PlotsPiPi.C","some options")
// root> T->Process("PlotsPiPi.C+")
//

#include "PlotsPiPi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <const.h>

TCanvas *c1, *c2;

std::vector<TH1 *> histKS;
std::vector<TH1 *> histKL;
TLegend legend(0.7, 0.5, 0.9, 0.9);

void PlotsPiPi::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   c1 = new TCanvas("", "", 790, 790);
   c2 = new TCanvas("", "", 790, 790);

   TString histKS_name = "";

   for (Int_t i = 0; i < channNum; i++)
   {
      histKS_name = "histKS_" + std::to_string(i);
      histKS.push_back(new TH1D(histKS_name, "", 100, -400, 400));
   };

   TString histKL_name = "";

   for (Int_t i = 0; i < channNum; i++)
   {
      histKL_name = "histKL_" + std::to_string(i);
      histKL.push_back(new TH1D(histKL_name, "", 100, -400, 400));
   };

   TString option = GetOption();
}

void PlotsPiPi::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t PlotsPiPi::Process(Long64_t entry)
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

   switch (*mctruth)
   {
   case 1:
   {
      histKS[0]->Fill(KchrecKS[2]);
      histKL[0]->Fill(KchrecKL[2]);
      break;
   }
   case 2:
   {
      histKS[1]->Fill(KchrecKS[2]);
      histKL[1]->Fill(KchrecKL[2]);
      break;
   }
   case 3:
   {
      histKS[2]->Fill(KchrecKS[2]);
      histKL[2]->Fill(KchrecKL[2]);
      break;
   }
   case 4:
   {
      histKS[3]->Fill(KchrecKS[2]);
      histKL[3]->Fill(KchrecKL[2]);
      break;
   }
   case 5:
   {
      histKS[4]->Fill(KchrecKS[2]);
      histKL[4]->Fill(KchrecKL[2]);
      break;
   }
   case 6:
   {
      histKS[5]->Fill(KchrecKS[2]);
      histKL[5]->Fill(KchrecKL[2]);
      break;
   }
   case 7:
   {
      histKS[6]->Fill(KchrecKS[2]);
      histKL[6]->Fill(KchrecKL[2]);
      break;
   }
   }

   return kTRUE;
}

void PlotsPiPi::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void PlotsPiPi::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   gStyle->SetOptStat(0);

   Float_t
      ymaxtrue = histKS[6]->GetMaximum(),
      ymax1 = ymaxtrue * 1.1;

   TLine *line1 = new TLine(497.611, 0, 497.611, ymax1 - 0.05*ymaxtrue);
   line1->SetLineColor(kBlack); // opcjonalnie kolor
   line1->SetLineWidth(2);    // opcjonalnie grubość

   for (Int_t i = 0; i < channNum; i++)
   {
      histKS[i]->SetLineColor(channColor[i]);
      histKL[i]->SetLineColor(channColor[i]);

      legend.AddEntry(histKS[i], channName[i], "l");
   }

   histKS[6]->GetXaxis()->SetTitle("m^{inv,K_{S}}_{#pi^{+}#pi^{-}} [MeV/c^{2}]");
   histKS[6]->GetYaxis()->SetTitle("Counts");
   histKS[6]->GetYaxis()->SetRangeUser(0, ymax1);

   histKL[6]->GetXaxis()->SetTitle("m^{inv,K_{L}}_{#pi^{+}#pi^{-}} [MeV/c^{2}]");
   histKL[6]->GetYaxis()->SetTitle("Counts");
   histKL[6]->GetYaxis()->SetRangeUser(0, ymax1);

   c1->cd();
   histKS[6]->Draw();
   histKS[0]->Draw("SAME");
   histKS[1]->Draw("SAME");
   histKS[2]->Draw("SAME");
   histKS[3]->Draw("SAME");
   histKS[4]->Draw("SAME");
   histKS[5]->Draw("SAME");
   // line1->Draw();
   legend.Draw();

   c1->SaveAs("KS_inv_mass.png");

   c2->cd();
   histKL[6]->Draw();
   histKL[0]->Draw("SAME");
   histKL[1]->Draw("SAME");
   histKL[2]->Draw("SAME");
   histKL[3]->Draw("SAME");
   histKL[4]->Draw("SAME");
   histKL[5]->Draw("SAME");
   // line1->Draw();
   legend.Draw();

   c2->SaveAs("KL_inv_mass.png");
}