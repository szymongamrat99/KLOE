#define neutral_vtx_resolutions_cxx
// The class definition in neutral_vtx_resolutions.h has been generated automatically
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
// root> T->Process("neutral_vtx_resolutions.C")
// root> T->Process("neutral_vtx_resolutions.C","some options")
// root> T->Process("neutral_vtx_resolutions.C+")
//

#include <string>

#include "neutral_vtx_resolutions.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFitResult.h>

TH2 *sigmas_std[5], *sigmas_tri[5], *sigmas_tri_kin_fit[5];
TH1 *neu_std_hist[5], *neu_tri_hist[5], *neu_tri_kin_fit_hist[5];
TH1 *res_std_hist[5], *res_tri_hist[5], *res_tri_kin_fit_hist[5];

TH2 *tchvstne[6];

TCanvas *canvas[30];

const UInt_t number_of_points = 10;

void neutral_vtx_resolutions::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TString name[4] = {"x", "y", "z", "length"};

   for (Int_t i = 0; i < 4; i++)
   {
      if (i < 3)
      {
         neu_std_hist[i] = new TH1F(name[i] + " std coor", "", 100, -200, 200);
         neu_tri_hist[i] = new TH1F(name[i] + " tri coor", "", 100, -200, 200);
         neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit coor", "", 100, -165, 165);
         res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 200);
         res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 200);
         res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 165);
      }
      else if (i == 3)
      {
         neu_std_hist[i] = new TH1F(name[i] + " std", "", 100, 0, 50);
         neu_tri_hist[i] = new TH1F(name[i] + " tri", "", 100, 0, 50);
         neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit", "", 100, 0, 50);
         res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 50);
         res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 50);
         res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 50);
      }

      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 50, 100, -5, 5);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 50, 100, -20, 20);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 50, 100, -20, 20);
   }

   for(Int_t i = 0; i < 6; i++)
   {

   }

   for (Int_t i = 0; i < 30; i++)
      canvas[i] = new TCanvas(("Canvas" + std::to_string(i)).c_str(), "", 750, 750);
}

void neutral_vtx_resolutions::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t neutral_vtx_resolutions::Process(Long64_t entry)
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

   Float_t radius_std, length_std, radius_mc, length_mc;

   fReader.SetLocalEntry(entry);

   if (*mctruth == 1)
   {
      length_mc = sqrt(pow(KnemcOld[6] - ipmcOld[0], 2) + pow(KnemcOld[7] - ipmcOld[1], 2) + pow(KnemcOld[8] - ipmcOld[2], 2));

      length_std = sqrt(pow(Knerec[6] - *Bx, 2) + pow(Knerec[7] - *By, 2) + pow(Knerec[8] - *Bz, 2));

      sigmas_std[0]->Fill(KnemcOld[6], Knerec[6] - KnemcOld[6]);
      sigmas_std[1]->Fill(KnemcOld[7], Knerec[7] - KnemcOld[7]);
      sigmas_std[2]->Fill(KnemcOld[8], Knerec[8] - KnemcOld[8]);
      sigmas_std[3]->Fill(length_mc, length_std - length_mc);
   }

   return kTRUE;
}

void neutral_vtx_resolutions::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void neutral_vtx_resolutions::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TFitResultPtr r;
   Float_t parameter[3];

   TF1 *func = new TF1("f1", "[0]*exp(-0.5*pow((x-[1])/[2],2))", -200, 200);
   func->SetParNames("Constant", "Mean", "Sigma");

   Color_t res_color[4] = {kRed, kBlue, kGreen, kBlack};

   for (Int_t j = 0; j < 4; j++)
   {
      for (Int_t i = 1; i <= number_of_points; i++)
      {
         parameter[0] = sigmas_std[j]->ProjectionY("_py", i, i)->GetMaximum();
         parameter[1] = sigmas_std[j]->ProjectionY("_py", i, i)->GetBinCenter(sigmas_std[j]->ProjectionY("_py", i, i)->GetMaximumBin());
         parameter[2] = 1.;

         func->SetParameters(parameter[0], parameter[1], parameter[2]);

         r = sigmas_std[j]->ProjectionY("_py", i, i)->Fit("f1", "S");
         res_std_hist[j]->SetBinContent(i, r->Parameter(2));
         res_std_hist[j]->SetBinError(i, r->ParError(2));
      }

      res_std_hist[j]->SetLineColor(res_color[j]);
   }

   TString title_x[2] = {"V_{gen,i}^{K_{ne}} [cm]",
                         "|V^{K_{ne}}_{gen} - V^{#Phi}_{gen}| [cm]"};

   TString title_y[2] = {"#sigma(V_{rec,i}^{K_{ne}} - V_{gen,i}^{K_{ne}}) [cm]",
                         "#sigma(V^{K_{ne}}_{rec} - V^{K_{ne}}_{gen}) [cm]"};

   gStyle->SetOptStat(0);

   canvas[0]->cd();
   res_std_hist[0]->SetXTitle(title_x[0]);
   res_std_hist[0]->SetYTitle(title_y[0]);
   res_std_hist[0]->GetYaxis()->SetRangeUser(0.2, 2.5);
   res_std_hist[0]->Draw("PE1");
   res_std_hist[1]->Draw("PE1SAME");
   res_std_hist[2]->Draw("PE1SAME");
   canvas[0]->Print(("sigmas_std" + std::to_string(0) + ".png").c_str());

   canvas[1]->cd();
   res_std_hist[3]->SetXTitle(title_x[1]);
   res_std_hist[3]->SetYTitle(title_y[1]);
   res_std_hist[3]->GetYaxis()->SetRangeUser(0.5, 1.2);
   res_std_hist[3]->Draw("PE1");
   canvas[1]->Print(("sigmas_std" + std::to_string(2) + ".png").c_str());
}