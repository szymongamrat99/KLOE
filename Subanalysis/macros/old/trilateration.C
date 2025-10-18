#define trilateration_cxx
// The class definition in trilateration.h has been generated automatically
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
// root> T->Process("trilateration.C")
// root> T->Process("trilateration.C","some options")
// root> T->Process("trilateration.C+")
//


#include "trilateration.h"
#include <TH2.h>
#include <TStyle.h>

#include "../Include/const.h"

TCanvas *canva1d, *canva2d[KLOE::channNum], *canvaproj;

TH1 *hist1d[2][KLOE::channNum];
TH2 *hist2d[2][KLOE::channNum];

TH1 *hist;

TLegend *legend;

Int_t entries[KLOE::channNum-2] = {0}, entries_sel[KLOE::channNum-2] = {0}, pred_entries_sel[KLOE::channNum-2] = {0};

void trilateration::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   canva1d = new TCanvas("canva1d", "canva1d", 750, 750);
   canvaproj = new TCanvas("canvaproj", "canvaproj", 750, 750); 
   for(Int_t i = 0; i < KLOE::channNum-2; i++) canva2d[i] = new TCanvas(channName[i] + " canva2d", channName[i] + " canva2d", 750, 750);

   for(Int_t i = 0; i < 2; i++)
      for(Int_t j = 0; j < KLOE::channNum; j++)
         hist1d[i][j] = new TH1F(channName[j] + i, "", 100, -50, 300);
   for(Int_t i = 0; i < 2; i++)
      for(Int_t j = 0; j < KLOE::channNum; j++)
         hist2d[i][j] = new TH2F(channName[j] + " 2d" + i, "", 50, 0, 60, 100, 0, 300);

   hist = new TH1F("Std Devs", "Std Devs", 13, 0, 60);

   legend = new TLegend(0.65,0.65,0.9,0.9);

   for(Int_t i = 0; i < KLOE::channNum-2; i++)
   {
      if(i == 8) legend->AddEntry(hist1d[1][i], channName[i], "PE1");
      else legend->AddEntry(hist1d[1][i], channName[i], "L");
   }

   TString option = GetOption();
}

void trilateration::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t trilateration::Process(Long64_t entry)
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

   Float_t distance_rec3, distance_rec2, distance_mc;

   fReader.SetLocalEntry(entry);

   if(*mctruth == 1) entries[0]++;
   if(*mctruth == 3) entries[1]++;
   if(*mctruth == 4) entries[2]++;
   if(*mctruth == 5) entries[3]++;
   if(*mctruth == 6) entries[4]++;
   if(*mctruth_pipi == 1) entries[5]++;
   if(*mctruth == 7) entries[6]++;

   distance_mc = sqrt(pow(Knemc_three[6] - *Bx, 2) + pow(Knemc_three[7] - *By, 2) + pow(Knemc_three[8] - *Bz, 2));

   if(*done4 == 1 && abs(*Qmiss_inv - 104.6) < 15 && abs(*anglepipi_CM_kch - 145) < 25 )
   {
      distance_rec2 = sqrt(pow(fourKnetri[6] - *Bx, 2) + pow(fourKnetri[7] - *By, 2) + pow(fourKnetri[8] - *Bz, 2));

      if(*mctruth == 1)
      {
         hist1d[0][0]->Fill(distance_rec2);
         hist2d[0][0]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[0]++;
      }
      if(*mctruth == 3)
      {
         hist1d[0][1]->Fill(distance_rec2);
         hist2d[0][1]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[1]++;
      }
      if(*mctruth == 4)
      {
         hist1d[0][2]->Fill(distance_rec2);
         hist2d[0][2]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[2]++;
      }
      if(*mctruth == 5)
      {
         hist1d[0][3]->Fill(distance_rec2);
         hist2d[0][3]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[3]++;
      }
      if(*mctruth == 6 )
      {
         hist1d[0][4]->Fill(distance_rec2);
         hist2d[0][4]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[4]++;
      }
      if(*mctruth_pipi == 1)
      {
         hist1d[0][5]->Fill(distance_rec2);
         hist2d[0][5]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[5]++;
      }
      if(*mctruth == 7)
      {
         hist1d[0][6]->Fill(distance_rec2);
         hist2d[0][6]->Fill(distance_mc, distance_rec2 - distance_mc);
//         entries_sel[6]++;
      }
   }

   if(*done == 1)
   {
      distance_rec3 = sqrt(pow(Knetri[6] - *Bx, 2) + pow(Knetri[7] - *By, 2) + pow(Knetri[8] - *Bz, 2));

      if(*mctruth == 1)
      {
         hist1d[1][0]->Fill(distance_rec3 - distance_mc);
         hist2d[1][0]->Fill(distance_mc, distance_rec3);
	 entries_sel[0]++;
      }
      if(*mctruth == 3)
      {
         hist1d[1][1]->Fill(distance_rec3 - distance_mc);
         hist2d[1][1]->Fill(distance_mc, distance_rec3);
	 entries_sel[1]++;

      }
      if(*mctruth == 4)
      {
         hist1d[1][2]->Fill(distance_rec3 - distance_mc);
         hist2d[1][2]->Fill(distance_mc, distance_rec3);
	 entries_sel[2]++;

      }
      if(*mctruth == 5)
      {
         hist1d[1][3]->Fill(distance_rec3 - distance_mc);
         hist2d[1][3]->Fill(distance_mc, distance_rec3);
	 entries_sel[3]++;

      }
      if(*mctruth == 6)
      {
         hist1d[1][4]->Fill(distance_rec3 - distance_mc);
         hist2d[1][4]->Fill(distance_mc, distance_rec3);
	 entries_sel[4]++;

      }
      if(*mctruth_pipi == 1)
      {
         hist1d[1][5]->Fill(distance_rec3 - distance_mc);
         hist2d[1][5]->Fill(distance_mc, distance_rec3);
	 entries_sel[5]++;
      }
      if(*mctruth == 7)
      {
         hist1d[1][6]->Fill(distance_rec3 - distance_mc);
         hist2d[1][6]->Fill(distance_mc, distance_rec3);
	 entries_sel[6]++;
      }
   }

   return kTRUE;
}

void trilateration::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void trilateration::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   gStyle->SetOptStat(0);

	canva1d->SetLeftMargin(0.15);
   canva1d->SetBottomMargin(0.15);
   canva1d->cd();
   canva1d->SetLogy(1);

   Int_t index_sort[KLOE::channNum];
   Double_t max_counts[KLOE::channNum];
   TString title[100];

   title[0] = "|V_{K^{0}}^{rec} - V_{#phi}^{rec}|" + units[0];
   title[1] = "|V_{K^{0}}^{rec} - V_{#phi}^{rec}| - |V_{K^{0}}^{gen} - V_{#phi}^{gen}|" + units[0];
   title[2] = "Counts";
   title[3] = "|V_{K^{0}}^{gen} - V_{#phi}^{gen}|" + units[0];
   title[4] = "|V_{K^{0}}^{rec} - V_{#phi}^{rec}| - |V_{K^{0}}^{gen} - V_{#phi}^{gen}|"; 

   Float_t x1d[2], y1d[2], x2d[2], y2d[2];

   for(Int_t i = 0; i < KLOE::channNum; i++) max_counts[i] = hist1d[0][i]->GetMaximum();

   TMath::Sort(KLOE::channNum, max_counts, index_sort);

   x1d[0] = -50;
   x1d[1] = 300;
   y1d[0] = 0.1;
   y1d[1] = 1E6;//max_counts[index_sort[0]] + 200.;

   x2d[0] = 0;
   x2d[1] = 60;
   y2d[0] = 0;
   y2d[1] = 300;

   for(Int_t i = 0; i < KLOE::channNum-2; i++)
   {
      hist1d[1][i]->SetLineColor(channColor[i]);

         if(i == 0)
         {
            hist1d[1][index_sort[i]]->GetXaxis()->SetMaxDigits(3);
            hist1d[1][index_sort[i]]->GetYaxis()->SetMaxDigits(3);

            hist1d[1][index_sort[i]]->GetXaxis()->SetTitle(title[1]);
            hist1d[1][index_sort[i]]->GetYaxis()->SetTitle(title[2]);

            hist1d[1][index_sort[i]]->GetXaxis()->CenterTitle(1);
            hist1d[1][index_sort[i]]->GetYaxis()->CenterTitle(1);

            hist1d[1][index_sort[i]]->GetXaxis()->SetRangeUser(x1d[0], x1d[1]);
            hist1d[1][index_sort[i]]->GetYaxis()->SetRangeUser(y1d[0], y1d[1]);
            hist1d[1][index_sort[i]]->Draw();
         }
         else
         {
            hist1d[1][index_sort[i]]->Draw("SAMES");
         }

   }

   legend->Draw();

   for(Int_t i = 0; i < KLOE::channNum-2; i++)
   {
      canva2d[i]->SetLeftMargin(0.15);
      canva2d[i]->SetRightMargin(0.15);
      canva2d[i]->SetBottomMargin(0.15);
      canva2d[i]->cd();

      hist2d[1][i]->GetXaxis()->SetMaxDigits(3);
      hist2d[1][i]->GetYaxis()->SetMaxDigits(3);
      hist2d[1][i]->GetZaxis()->SetMaxDigits(3);

      hist2d[1][i]->GetXaxis()->SetTitle(title[3]);
      hist2d[1][i]->GetYaxis()->SetTitle(title[0]);

      hist2d[1][i]->GetXaxis()->CenterTitle(1);
      hist2d[1][i]->GetYaxis()->CenterTitle(1);

      hist2d[1][i]->GetXaxis()->SetRangeUser(x2d[0], x2d[1]);
      hist2d[1][i]->GetYaxis()->SetRangeUser(y2d[0], y2d[1]);
      hist2d[1][i]->Draw("COLZ");
   }

   TFitResultPtr r;

   canvaproj->cd();

   Float_t parameter[3];

   TF1* func = new TF1("f1","[0]*exp(-0.5*pow((x-[1])/[2],2))",-30,200);
   func->SetParNames("Constant", "Mean", "Sigma");

   for(Int_t i = 1; i <= 13; i++)
   {
      parameter[0] = hist2d[1][3]->ProjectionY("_py", i, i)->GetMaximum();   
      parameter[1] = hist2d[1][3]->ProjectionY("_py", i, i)->GetBinCenter(hist2d[1][3]->ProjectionY("_py", i, i)->GetMaximumBin());
      parameter[2] = 1.;

      func->SetParameters(parameter[0], parameter[1], parameter[2]);

      r = hist2d[1][3]->ProjectionY("_py", i, i)->Fit("f1","S");
      hist->SetBinContent(i, r->Parameter(2));
      hist->SetBinError(i, r->ParError(2));
   }

   hist->GetYaxis()->SetRangeUser(0,30);
   hist->GetXaxis()->CenterTitle(1);
   hist->GetYaxis()->CenterTitle(1);
   hist->GetXaxis()->SetTitle(title[3]);
   hist->SetTitle("");
   hist->GetYaxis()->SetTitle("#sigma(" + title[4] + ")" + units[0]);
   hist->Draw("PE1");

   cout << "Purity: " << entries_sel[3]/(Float_t)(entries_sel[0] + entries_sel[1] + entries_sel[2] + entries_sel[3] + entries_sel[4] + entries_sel[5] + entries_sel[6]) << endl;
   cout << "Efficiency: " << entries_sel[3]/(Float_t)(entries[3]) << endl;
   cout << "Predicted efficiency: " << entries_sel[3]/(Float_t)(entries[3]) << endl;

   canva1d->Print("tri_error.png");
   for(Int_t i = 0; i < KLOE::channNum - 2; i++) canva2d[i]->Print(("tri_2d_plane_corr_" + to_string(i) + ".png").c_str());
   canvaproj->Print("sigma_tri.png");
}
