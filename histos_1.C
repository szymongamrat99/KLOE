#define histos_1_cxx
// The class definition in histos_1.h has been generated automatically
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
// root> T->Process("histos_1.C")
// root> T->Process("histos_1.C","some options")
// root> T->Process("histos_1.C+")
//


#include "histos_1.h"
#include <TH2.h>
#include <TStyle.h>

#include "../Include/const.h"

TCanvas *canva1d, *canva2d[chann_num], *canvaproj[chann_num];

TH1 *hist1d[2][chann_num];
TH2 *hist2d[2][chann_num];

TLegend *legend;

Int_t entries[chann_num-2] = {0}, entries_sel[chann_num-2] = {0}, pred_entries_sel[chann_num-2] = {0};
void histos_1::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   canva1d = new TCanvas("canva1d", "canva1d", 750, 750);
   for(Int_t i = 0; i < chann_num-2; i++) canva2d[i] = new TCanvas(chann_name[i] + " canva2d", chann_name[i] + " canva2d", 750, 750);
   for(Int_t i = 0; i < chann_num-2; i++) canvaproj[i] = new TCanvas(chann_name[i] + " canvaproj", chann_name[i] + " canvaproj", 750, 750);

   for(Int_t i = 0; i < 2; i++)
      for(Int_t j = 0; j < chann_num; j++)
         hist1d[i][j] = new TH1F(chann_name[j] + i, "", 500, 0, 1000);
   for(Int_t i = 0; i < 2; i++)
      for(Int_t j = 0; j < chann_num; j++)
         hist2d[i][j] = new TH2F(chann_name[j] + " 2d" + i, "", 13, 0, 60, 100, -100, 300);


   legend = new TLegend(0.65,0.65,0.9,0.9);
   //legend = new TLegend(0.15,0.65,0.4,0.9);

   for(Int_t i = 0; i < chann_num-2; i++)
   {
      if(i == 7) legend->AddEntry(hist1d[0][i], chann_name[i], "PE1");
      else legend->AddEntry(hist1d[0][i], chann_name[i], "L");
   }


   TString option = GetOption();
}

void histos_1::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t histos_1::Process(Long64_t entry)
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

   if(*mctruth == 1) entries[0]++;
   if(*mctruth == 3) entries[1]++;
   if(*mctruth == 4) entries[2]++;
   if(*mctruth == 5) entries[3]++;
   if(*mctruth == 6) entries[4]++;
   if(*mctruth_pipi == 1) entries[5]++;
   if(*mctruth == 7) entries[6]++;

   if(1)//abs(*Qmiss_inv - 104.6) < 15 && abs(*anglepipi_CM_kch - 145) < 25 && *done4 == 1 && abs(Kchrec1[5] - m_k0) < 76)
   {
   	   if(*mctruth == 1)
	   {
		   hist1d[0][0]->Fill(Kchrec2[5]);
		   entries_sel[0]++;
	   }
	   if(*mctruth == 3)
	   {
		   hist1d[0][1]->Fill(Kchrec2[5]);
		   entries_sel[1]++;
	   }
	   if(*mctruth == 4)
	   {
		   hist1d[0][2]->Fill(Kchrec2[5]);
		   entries_sel[2]++;
	   }
	   if(*mctruth == 5)
	   {
		   hist1d[0][3]->Fill(Kchrec2[5]);
		   entries_sel[3]++;
	   }
	   if(*mctruth == 6)
	   {
		   hist1d[0][4]->Fill(Kchrec2[5]);
		   entries_sel[4]++;
	   }
	   if(*mctruth_pipi == 1)
	   {
		   hist1d[0][5]->Fill(Kchrec2[5]);
		   entries_sel[5]++;
	   }
	   if(*mctruth == 7)
	   {
		   hist1d[0][6]->Fill(Kchrec2[5]);
		   entries_sel[6]++;
	   }
   }
   return kTRUE;
}

void histos_1::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void histos_1::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   gStyle->SetOptStat(0);

   canva1d->SetLeftMargin(0.15);
   canva1d->SetBottomMargin(0.15);
   canva1d->cd();
   canva1d->SetLogy(1);

   Int_t index_sort[chann_num];
   Double_t max_counts[chann_num];
   TString title[100];

   title[0] = "|V_{K^{0}}^{rec} - V_{#phi}^{rec}|" + units[0];
   title[1] = "|V_{K^{0}}^{rec} - V_{#phi}^{rec}|" + units[0];
   title[2] = "Counts";
   title[3] = "|V_{K^{0}}^{gen} - V_{#phi}^{gen}|" + units[0];
   title[4] = "m^{inv}_{#pi^{+}#pi^{-}}" + units[3];
   title[5] = "m^{inv}_{#pi^{0}#pi^{0},standard}" + units[3];
   title[6] = "m^{inv}_{#pi^{0}#pi^{0},trilateration}" + units[3];
   title[7] = "m^{inv}_{3#pi^{0},trilateration}" + units[3];
   title[8] = "Q_{miss}" + units[5];
   title[9] = "#alpha^{CM}_{#pi^{+},#pi^{-}} [#circ]";
   title[10] ="m^{inv}_{#pi^{+}#pi^{-},K_{1}}" + units[3];  
   title[11] ="m^{inv}_{#pi^{+}#pi^{-},K_{2}}" + units[3];  

   Float_t x1d[2], y1d[2], x2d[2], y2d[2]; 

   for(Int_t i = 0; i < chann_num; i++) max_counts[i] = hist1d[0][i]->GetMaximum();

   TMath::Sort(chann_num, max_counts, index_sort);

   x1d[0] = 0;
   x1d[1] = 1000;
   y1d[0] = 0.1;
   y1d[1] = 1E5;//max_counts[index_sort[0]] + 200.;

   x2d[0] = 0;
   x2d[1] = 60;
   y2d[0] = -100;
   y2d[1] = 300;

   for(Int_t i = 0; i < chann_num-2; i++)
   {
      hist1d[0][i]->SetLineColor(chann_color[i]);

         if(i == 0)
         {
            hist1d[0][index_sort[i]]->GetXaxis()->SetMaxDigits(3);
            hist1d[0][index_sort[i]]->GetYaxis()->SetMaxDigits(3);

            hist1d[0][index_sort[i]]->GetXaxis()->SetTitle(title[11]);
            hist1d[0][index_sort[i]]->GetYaxis()->SetTitle(title[2]);

            hist1d[0][index_sort[i]]->GetXaxis()->CenterTitle(1);
            hist1d[0][index_sort[i]]->GetYaxis()->CenterTitle(1);

            hist1d[0][index_sort[i]]->GetXaxis()->SetRangeUser(x1d[0], x1d[1]);
            hist1d[0][index_sort[i]]->GetYaxis()->SetRangeUser(y1d[0], y1d[1]);
            hist1d[0][index_sort[i]]->Draw();
         }
         else
         {
            hist1d[0][index_sort[i]]->Draw("SAMES");
         }
      
   }

   legend->Draw();

   for(Int_t i = 0; i < chann_num-2; i++)
   {
      canva2d[i]->SetLeftMargin(0.15);
      canva2d[i]->SetRightMargin(0.15);
      canva2d[i]->SetBottomMargin(0.15);
      canva2d[i]->cd();

      hist2d[0][i]->GetXaxis()->SetMaxDigits(3);
      hist2d[0][i]->GetYaxis()->SetMaxDigits(3);
      hist2d[0][i]->GetZaxis()->SetMaxDigits(3);

      hist2d[0][i]->GetXaxis()->SetTitle(title[3]);
      hist2d[0][i]->GetYaxis()->SetTitle(title[1]);

      hist2d[0][i]->GetXaxis()->CenterTitle(1);
      hist2d[0][i]->GetYaxis()->CenterTitle(1);

      hist2d[0][i]->GetXaxis()->SetRangeUser(x2d[0], x2d[1]);
      hist2d[0][i]->GetYaxis()->SetRangeUser(y2d[0], y2d[1]);
      hist2d[0][i]->Draw("COLZ");
   }

   canva1d->Print("kchrec2_minv.png");

   cout << "Purity: " << entries_sel[4]/(Float_t)(entries_sel[0] + entries_sel[1] + entries_sel[2] + entries_sel[3] + entries_sel[4] + entries_sel[5] + entries_sel[6]) << endl;
   cout << "Efficiency: " << entries_sel[4]/(Float_t)(entries[4]) << endl;
   cout << "Predicted efficiency: " << entries_sel[4]/(Float_t)(entries[4]) << endl;
}
