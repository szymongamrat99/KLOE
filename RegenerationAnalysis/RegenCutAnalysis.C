#define RegenCutAnalysis_cxx
// The class definition in RegenCutAnalysis.h has been generated automatically
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
// root> T->Process("RegenCutAnalysis.C")
// root> T->Process("RegenCutAnalysis.C","some options")
// root> T->Process("RegenCutAnalysis.C+")
//

#include <iostream>
#include <fstream>

#include "RegenCutAnalysis.h"

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TEfficiency.h>
#include <TEventList.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>
#include <TFitResult.h>

using namespace std;

/* Global creation of canvas and histograms */

Double_t height = 750, width = 750;
TString dir = "img/", ext = ".svg";

TCanvas *canvas[100];
TH1 *h_radius[2], *h_radius_tri[2], *h_radius_ch[2], *h_dt;

TEventList *elist = 0;
Bool_t useList, fillList;

void RegenCutAnalysis::Begin(TTree *tree)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Int_t done4 = 0, g4takentri_kinfit[4], nentries = 0;
Float_t fourKnetri[10] = {0.};

void RegenCutAnalysis::SlaveBegin(TTree */*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TFile file_tri("../Neutrec/root_files/neuvtx_tri_kin_fit_1_56_10_9_1_2_bunch_corr_automatic.root");
   tree_tri = (TTree *)file_tri.Get("h_tri_kin_fit");

   tree_tri->SetBranchAddress("fourKnetri_kinfit", fourKnetri);
   tree_tri->SetBranchAddress("done4_kinfit", &done4);

   nentries = tree_tri->GetEntries();

   TString hist_name, canvas_name;

   for (Int_t i = 0; i < 100; i++)
   {
      canvas_name = "Canva_" + to_string(i);
      canvas[i] = new TCanvas(canvas_name, canvas_name, width, height);
   }

   for (Int_t i = 0; i < 2; i++)
   {
      hist_name = "Hist_radius_" + to_string(i);
      h_radius[i] = new TH1D(hist_name, hist_name, 100, 0, 50.0);

      hist_name = "Hist_radius_tri_" + to_string(i);
      h_radius_tri[i] = new TH1D(hist_name, hist_name, 100, 0, 50.0);

      hist_name = "Hist_radius_ch_" + to_string(i);
      h_radius_ch[i] = new TH1D(hist_name, hist_name, 100, 0, 50.0);
   }

   h_dt = new TH1D("Delta T", "Delta T", 200, -90.0, 90.0);

   Double_t radius[2] = {0., 0.};

   Double_t sphere_bound = 10, bp_bound = 4.4;

   Double_t sigma = 1.5;

   Bool_t cuts[4] = {false, false, false, false};

   for(Int_t i = 0; i < nentries; i++)
   {
      tree_tri->GetEntry(i);


   }

   TString option = GetOption();

   printf("Starting analysis with process option: %s", option.Data());

   // process cases with event list
   fillList = kFALSE;
   useList = kFALSE;
}

Bool_t RegenCutAnalysis::Process(Long64_t entry)
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

   Double_t radius[2] = {0., 0.}, radius_ch[2] = {0., 0.};

   Double_t sphere_bound = 10, bp_bound = 4.4;

   Double_t sigma = 1.5;

   Bool_t cuts[4] = {false, false, false, false};

   fReader.SetLocalEntry(entry);

   for (Int_t i = 0; i < 3; i++)
   {
      radius[0] += pow(fourKnetri[6 + i], 2);
      radius_ch[0] += pow(Kchboost[6 + i], 2);

      if (i < 2)
      {
         radius[1] += pow(fourKnetri[6 + i], 2);
         radius_ch[1] += pow(Kchboost[6 + i], 2);
      }
   }

   radius[0] = sqrt(radius[0]);
   radius[1] = sqrt(radius[1]);

   radius_ch[0] = sqrt(radius_ch[0]);
   radius_ch[1] = sqrt(radius_ch[1]);

   if (*mctruth != 0 && *mctruth != 2)
   {
      if (abs(Knereclor[8]) < 10)
      {
         h_radius[0]->Fill(radius[0]);
      }
      else
      {
         h_radius[1]->Fill(radius[1]);
      }

      if (abs(Kchboost[8]) < 10)
      {
         h_radius_ch[0]->Fill(radius_ch[0]);
      }
      else
      {
         h_radius_ch[1]->Fill(radius_ch[1]);
      }

      cuts[0] = abs(radius[0] - 11.1134) > 1.48753 * sigma;
      cuts[1] = abs(radius[1] - 4.9646) > 0.652839 * sigma;
      cuts[2] = abs(radius_ch[0] - 10.5195) > 0.982962 * sigma;
      cuts[3] = abs(radius_ch[1] - 4.85379) > 1.06893 * sigma;

      if (cuts[0] && cuts[1] && cuts[2] && cuts[3])
         h_dt->Fill(*Dtboostlor);
   }

   return kTRUE;
}

TFitResultPtr r;

void RegenCutAnalysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   TString img_name;

   for (Int_t i = 0; i < 2; i++)
   {
      img_name = dir + "radius_" + i + ext;
      canvas[i]->cd();

      h_radius[i]->Draw();

      if (i == 0)
      {
         r = h_radius[i]->Fit("gaus", "SM", "", 9., 13.);
         cout << "Center of radius 0: (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
         cout << "Std dev of radius 0: (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
      }
      else if (i == 1)
      {
         r = h_radius[i]->Fit("gaus", "SM", "", 3., 6.);
         cout << "Center of radius 1: (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
         cout << "Std dev of radius 1: (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
      }

      canvas[i]->Print(img_name);
   }

   for (Int_t i = 2; i < 4; i++)
   {
      img_name = dir + "radius_ch" + (i - 2) + ext;
      canvas[i]->cd();

      if (i == 2)
      {
         r = h_radius_ch[i - 2]->Fit("gaus", "SM", "", 9., 12.);
         cout << "Center of radius ch 0: (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
         cout << "Std dev of radius ch 0: (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
      }
      else if (i == 3)
      {
         r = h_radius_ch[i - 2]->Fit("gaus", "SM", "", 3., 6.);
         cout << "Center of radius ch 1: (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
         cout << "Std dev of radius ch 1: (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
      }

      h_radius_ch[i - 2]->Draw();
      canvas[i]->Print(img_name);
   }

   img_name = dir + "deltaT" + ext;
   canvas[4]->cd();

   h_dt->Draw();
   canvas[4]->Print(img_name);
}

void RegenCutAnalysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
}