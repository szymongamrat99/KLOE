#define histos_cxx
// The class definition in histos.h has been generated automatically
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
// root> T->Process("histos.C")
// root> T->Process("histos.C","some options")
// root> T->Process("histos.C+")
//


#include "histos.h"
#include <TH2.h>
#include <TStyle.h>

#include "../Include/const.h"

TCanvas *canva1d_semi[10], *canva1d_three[10], *canva1d_pipi[10], *canva2d[10][KLOE::channNum], *canvaproj[10][KLOE::channNum];
TCanvas *canva1d_seminocuts, *canva1d_threenocuts, *canva1d_pipinocuts;
TCanvas *canva1d_semiwithcuts, *canva1d_threewithcuts, *canva1d_pipiwithcuts;
TCanvas *canva_efficiency, *canva_eff_signal;


TH1 *hist_semi[10][KLOE::channNum], *hist_three[10][KLOE::channNum], *hist_pipi[10][KLOE::channNum];
TH1 *hist_signal, *hist_signal_nocuts;
TH1 *hist_semi_nocuts, *hist_three_nocuts, *hist_pipi_nocuts;
TH1 *hist_semi_withcuts, *hist_three_withcuts, *hist_pipi_withcuts;
TH2 *hist2d[10][KLOE::channNum];

TLegend *legend[10];

Int_t entries[KLOE::channNum] = {0}, entries_sel[KLOE::channNum] = {0}, pred_entries_sel[KLOE::channNum] = {0};
void histos::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   for(Int_t i = 0; i < 10; i++) canva1d_semi[i] = new TCanvas(("canva1d_semi" + to_string(i)).c_str(), "", 750, 750);
   for(Int_t i = 0; i < 10; i++) canva1d_three[i] = new TCanvas(("canva1d_three" + to_string(i)).c_str(), "", 750, 750);
   for(Int_t i = 0; i < 10; i++) canva1d_pipi[i] = new TCanvas(("canva1d_pipi" + to_string(i)).c_str(), "", 750, 750);

   canva_efficiency = new TCanvas("Canva efficiency", "", 750, 750);
   canva_eff_signal = new TCanvas("Canva eff signal", "", 750, 750);
   canva1d_seminocuts = new TCanvas("Canva semi", "", 750, 750);
   canva1d_threenocuts = new TCanvas("Canva three", "", 750, 750);
   canva1d_pipinocuts = new TCanvas("Canva pipi", "", 750, 750);
   canva1d_semiwithcuts = new TCanvas("Canva semi cuts", "", 750, 750);
   canva1d_threewithcuts = new TCanvas("Canva three cuts", "", 750, 750);
   canva1d_pipiwithcuts = new TCanvas("Canva pipi cuts", "", 750, 750);

  // for(Int_t i = 0; i < 10; i++)
  //    for(Int_t j = 0; j < KLOE::channNum-2; j++) canva2d[i][j] = new TCanvas((channName[j] + " canva2d" + to_string(i)).c_str(), "", 750, 750);
  // for(Int_t i = 0; i < 10; i++)
  //    for(Int_t j = 0; j < KLOE::channNum-2; j++) canvaproj[i][j] = new TCanvas((channName[j] + " canvaproj" + to_string(i)).c_str(), "", 750, 750);

   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[0][j] = new TH1F(channName[j] + "0semi", "", 201, -100, 100);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[1][j] = new TH1F(channName[j] + "1semi", "", 100, 0, 50);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[2][j] = new TH1F(channName[j] + "2semi", "", 100, 0, 50);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[3][j] = new TH1F(channName[j] + "3semi", "", 200, -10, 10);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[4][j] = new TH1F(channName[j] + "4semi", "", 200, -10, 10);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[5][j] = new TH1F(channName[j] + "5semi", "", 100, -100, 5);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[6][j] = new TH1F(channName[j] + "6semi", "", 200, 300, 1000);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_semi[7][j] = new TH1F(channName[j] + "7semi", "", 200, 300, 1000);

   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[0][j] = new TH1F(channName[j] + "0three", "", 201, -100, 100);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[1][j] = new TH1F(channName[j] + "1three", "", 100, 0, 50);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[2][j] = new TH1F(channName[j] + "2three", "", 100, 0, 50);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[3][j] = new TH1F(channName[j] + "3three", "", 200, -10, 10);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[4][j] = new TH1F(channName[j] + "4three", "", 200, -10, 10);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[5][j] = new TH1F(channName[j] + "5three", "", 30, -1.0, 0.0);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[6][j] = new TH1F(channName[j] + "6three", "", 50, PhysicsConstants::mK0 - 5, PhysicsConstants::mK0 + 5);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_three[7][j] = new TH1F(channName[j] + "7three", "", 50, PhysicsConstants::mK0 - 5, PhysicsConstants::mK0 + 5);

   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[0][j] = new TH1F(channName[j] + "0pipi", "", 201, -100, 100);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[1][j] = new TH1F(channName[j] + "1pipi", "", 100, 0, 50);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[2][j] = new TH1F(channName[j] + "2pipi", "", 100, 0, 50);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[3][j] = new TH1F(channName[j] + "3pipi", "", 200, -10, 10);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[4][j] = new TH1F(channName[j] + "4pipi", "", 200, -10, 10);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[5][j] = new TH1F(channName[j] + "5pipi", "", 100, -100, 5);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[6][j] = new TH1F(channName[j] + "6pipi", "", 200, 300, 1000);
   for(Int_t j = 0; j < KLOE::channNum; j++) hist_pipi[7][j] = new TH1F(channName[j] + "7pipi", "", 200, 300, 1000);

   hist_signal = new TH1F("Signal histo", "", 100, -100, 10);
   hist_signal_nocuts = new TH1F("Signal histo nocuts", "", 201, -100, 100);
   hist_semi_nocuts = new TH1F("semi histo", "", 201, -100, 100);
   hist_three_nocuts = new TH1F("three histo", "", 201, -100, 100);
   hist_pipi_nocuts = new TH1F("pipi histo", "", 201, -100, 100);
   hist_semi_withcuts = new TH1F("semi histo cuts", "", 201, -100, 100);
   hist_three_withcuts = new TH1F("three histo cuts", "", 201, -100, 100);
   hist_pipi_withcuts = new TH1F("pipi histo cuts", "", 201, -100, 100);

//   for(Int_t i = 0; i < 10; i++)
//      for(Int_t j = 0; j < KLOE::channNum; j++)
//         hist2d[i][j] = new TH2F(channName[j] + " 2d" + to_string(i), "", 13, 0, 60, 100, -100, 300);


   for(Int_t i = 0; i < 10; i++)
   {
      //legend[i] = new TLegend(0.65,0.65,0.9,0.9);
      legend[i] = new TLegend(0.65,0.25,0.9,0.5);

      for(Int_t j = 0; j < KLOE::channNum; j++)
      {
         if(j == 6) legend[i]->AddEntry(hist_semi[i][j], channName[j], "PE1");
         else legend[i]->AddEntry(hist_semi[i][j], channName[j], "L");
      }
   }


   TString option = GetOption();
}

void histos::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t histos::Process(Long64_t entry)
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

   Float_t TRCV[4], k_path00, k_pathpm, k_betapm, k_beta00, trcv_sum, trcv_sum_signal;
   Float_t tpm, t00, DeltaT;

   Bool_t cuts_signal_charged_cs[3], cuts_signal_neutral_cs[2], cuts_semi[2], cuts_three, cuts_pipi;
   Bool_t tot_cuts_semi, tot_cuts_neutral, tot_cuts_charged;

   fReader.SetLocalEntry(entry);

   k_beta00 = sqrt(pow(fourKnetri[0],2) + pow(fourKnetri[1],2) + pow(fourKnetri[2],2))/fourKnetri[3];
   k_path00 = sqrt(pow(fourKnetri[6] - *Bx,2) + pow(fourKnetri[7] - *By,2) + pow(fourKnetri[8] - *Bz,2));

   for(Int_t i = 0; i < 4; i++) TRCV[i] = TclOld[fourg4taken[i]] - (sqrt(pow(Xcl[fourg4taken[i]] - fourKnetri[6],2) + pow(Ycl[fourg4taken[i]] - fourKnetri[7],2) + pow(Zcl[fourg4taken[i]] - fourKnetri[8],2))/PhysicsConstants::cVel) - (k_path00/(k_beta00*PhysicsConstants::cVel));

   trcv_sum = TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3];
	
   k_betapm = sqrt(pow(Kchboost[0],2) + pow(Kchboost[1],2) + pow(Kchboost[2],2))/Kchboost[3];

   k_pathpm = sqrt(pow(Kchboost[6] - *Bx,2) + pow(Kchboost[7] - *By,2) + pow(Kchboost[8] - *Bz,2));

   //Calculation of time difference

   t00 = fourKnetri[9]/PhysicsConstants::tau_S_nonCPT;
   tpm = k_pathpm/(k_betapm*PhysicsConstants::cVel*PhysicsConstants::tau_S_nonCPT); 

   DeltaT = tpm - t00;

   cuts_semi[0] = abs(*Qmiss_inv - 104.6) < 15;
   cuts_semi[1] = abs(*anglepipi_CM_kch - 145) < 25;

   cuts_signal_neutral_cs[0] = (abs(fourKnetri[5] - PhysicsConstants::mK0) < 30);
   cuts_signal_neutral_cs[1] = (trcv_sum > -5);

   cuts_signal_charged_cs[0] = (abs(Kchrec[5] - PhysicsConstants::mK0) < 1.2);
   cuts_signal_charged_cs[1] = (*Qmiss_inv < 3.75);
   cuts_signal_charged_cs[2] = (cos(M_PI**anglepipi_CM_kch/180.) < -0.8);

   tot_cuts_semi = cuts_semi[0] && cuts_semi[1];
   tot_cuts_neutral = cuts_signal_neutral_cs[0] && cuts_signal_neutral_cs[1];
   tot_cuts_charged = cuts_signal_charged_cs[0] && cuts_signal_charged_cs[1] && cuts_signal_charged_cs[2];

   trcv_sum_signal = trcv[ncll[g4taken[0]-1]-1] + trcv[ncll[g4taken[1]-1]-1] + trcv[ncll[g4taken[2]-1]-1] + trcv[ncll[g4taken[3]-1]-1];

   //Signal selection

   if(*mctruth == 1)
	{
		hist_signal_nocuts->Fill(*Dtrec);
   }

   if(abs(*minv4gam - PhysicsConstants::mK0) < 150 && trcv_sum_signal > -1 && abs(Kchrec[5] - PhysicsConstants::mK0) < 2 && *Qmiss_inv < 3.75 && cos(M_PI**anglepipi_CM_kch/180.) < -0.8)
   {
   	if(*mctruth == 1)
	   {
		   hist_signal->Fill(*Dtrec);
	   }
   }

   //Semi selection

   if(tot_cuts_semi && *done4 == 1)
   {
      hist_semi_nocuts->Fill(DeltaT, ((br_kl_piele + br_kl_pimu)/(br_kl_piele + br_kl_pimu + br_ks_piele + br_ks_pimu))*br_ks_pi0pi0 + ((br_ks_piele + br_ks_pimu)/(br_kl_piele + br_kl_pimu + br_ks_piele + br_ks_pimu))*br_kl_pi0pi0);


      if(tot_cuts_neutral)
      {
          hist_semi_withcuts->Fill(DeltaT);

         if(*mctruth == 1)
         {
            hist_semi[0][0]->Fill(DeltaT);
            hist_semi[1][0]->Fill(k_pathpm);
            hist_semi[2][0]->Fill(k_path00);
            hist_semi[3][0]->Fill(fourKnetri[8]);
            hist_semi[4][0]->Fill(Kchboost[8]);
            hist_semi[5][0]->Fill(trcv_sum);
            hist_semi[6][0]->Fill(Kchrec[5]);
            hist_semi[7][0]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 3)
         {
            hist_semi[0][1]->Fill(DeltaT);
            hist_semi[1][1]->Fill(k_pathpm);
            hist_semi[2][1]->Fill(k_path00);
            hist_semi[3][1]->Fill(fourKnetri[8]);
            hist_semi[4][1]->Fill(Kchboost[8]);
            hist_semi[5][1]->Fill(trcv_sum);
            hist_semi[6][1]->Fill(Kchrec[5]);
            hist_semi[7][1]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 4)
         {
            hist_semi[0][2]->Fill(DeltaT);
            hist_semi[1][2]->Fill(k_pathpm);
            hist_semi[2][2]->Fill(k_path00);
            hist_semi[3][2]->Fill(fourKnetri[8]);
            hist_semi[4][2]->Fill(Kchboost[8]);
            hist_semi[5][2]->Fill(trcv_sum);
            hist_semi[6][2]->Fill(Kchrec[5]);
            hist_semi[7][2]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 5)
         {
            hist_semi[0][3]->Fill(DeltaT);
            hist_semi[1][3]->Fill(k_pathpm);
            hist_semi[2][3]->Fill(k_path00);
            hist_semi[3][3]->Fill(fourKnetri[8]);
            hist_semi[4][3]->Fill(Kchboost[8]);
            hist_semi[5][3]->Fill(trcv_sum);
            hist_semi[6][3]->Fill(Kchrec[5]);
            hist_semi[7][3]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 6)
         {
            hist_semi[0][4]->Fill(DeltaT);
            hist_semi[1][4]->Fill(k_pathpm);
            hist_semi[2][4]->Fill(k_path00);
            hist_semi[3][4]->Fill(fourKnetri[8]);
            hist_semi[4][4]->Fill(Kchboost[8]);
            hist_semi[5][4]->Fill(trcv_sum);
            hist_semi[6][4]->Fill(Kchrec[5]);
            hist_semi[7][4]->Fill(fourKnetri[5]);
         }
         if(*mctruth_pipi == 1)
         {
            hist_semi[0][5]->Fill(DeltaT);
            hist_semi[1][5]->Fill(k_pathpm);
            hist_semi[2][5]->Fill(k_path00);
            hist_semi[3][5]->Fill(fourKnetri[8]);
            hist_semi[4][5]->Fill(Kchboost[8]);
            hist_semi[5][5]->Fill(trcv_sum);
            hist_semi[6][5]->Fill(Kchrec[5]);
            hist_semi[7][5]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 7)
         {
            hist_semi[0][6]->Fill(DeltaT);
            hist_semi[1][6]->Fill(k_pathpm);
            hist_semi[2][6]->Fill(k_path00);
            hist_semi[3][6]->Fill(fourKnetri[8]);
            hist_semi[4][6]->Fill(Kchboost[8]);
            hist_semi[5][6]->Fill(trcv_sum);
            hist_semi[6][6]->Fill(Kchrec[5]);
            hist_semi[7][6]->Fill(fourKnetri[5]);
         }
	 if(*mcflag == 0)
         {
            hist_semi[0][7]->Fill(DeltaT);
            hist_semi[1][7]->Fill(k_pathpm);
            hist_semi[2][7]->Fill(k_path00);
            hist_semi[3][7]->Fill(fourKnetri[8]);
            hist_semi[4][7]->Fill(Kchboost[8]);
            hist_semi[5][7]->Fill(trcv_sum);
            hist_semi[6][7]->Fill(Kchrec[5]);
            hist_semi[7][7]->Fill(fourKnetri[5]);
	 }
      }
   }

   //Three selection
   if(*done == 1 && trcv_sum > -10 && *done4 == 1)
   {
      hist_three_nocuts->Fill(DeltaT, br_ks_pippim);

      if(tot_cuts_charged)
      {
         hist_three_withcuts->Fill(DeltaT);

         if(*mctruth == 1)
         {
            hist_three[0][0]->Fill(DeltaT);
            hist_three[1][0]->Fill(k_pathpm);
            hist_three[2][0]->Fill(k_path00);
            hist_three[3][0]->Fill(fourKnetri[8]);
            hist_three[4][0]->Fill(Kchboost[8]);
            hist_three[5][0]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][0]->Fill(Kchrec[5]);
            hist_three[7][0]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 3)
         {
            hist_three[0][1]->Fill(DeltaT);
            hist_three[1][1]->Fill(k_pathpm);
            hist_three[2][1]->Fill(k_path00);
            hist_three[3][1]->Fill(fourKnetri[8]);
            hist_three[4][1]->Fill(Kchboost[8]);
            hist_three[5][1]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][1]->Fill(Kchrec[5]);
            hist_three[7][1]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 4)
         {
            hist_three[0][2]->Fill(DeltaT);
            hist_three[1][2]->Fill(k_pathpm);
            hist_three[2][2]->Fill(k_path00);
            hist_three[3][2]->Fill(fourKnetri[8]);
            hist_three[4][2]->Fill(Kchboost[8]);
            hist_three[5][2]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][2]->Fill(Kchrec[5]);
            hist_three[7][2]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 5)
         {
            hist_three[0][3]->Fill(DeltaT);
            hist_three[1][3]->Fill(k_pathpm);
            hist_three[2][3]->Fill(k_path00);
            hist_three[3][3]->Fill(fourKnetri[8]);
            hist_three[4][3]->Fill(Kchboost[8]);
            hist_three[5][3]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][3]->Fill(Kchrec[5]);
            hist_three[7][3]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 6)
         {
            hist_three[0][4]->Fill(DeltaT);
            hist_three[1][4]->Fill(k_pathpm);
            hist_three[2][4]->Fill(k_path00);
            hist_three[3][4]->Fill(fourKnetri[8]);
            hist_three[4][4]->Fill(Kchboost[8]);
            hist_three[5][4]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][4]->Fill(Kchrec[5]);
            hist_three[7][4]->Fill(fourKnetri[5]);
         }
         if(*mctruth_pipi == 1)
         {
            hist_three[0][5]->Fill(DeltaT);
            hist_three[1][5]->Fill(k_pathpm);
            hist_three[2][5]->Fill(k_path00);
            hist_three[3][5]->Fill(fourKnetri[8]);
            hist_three[4][5]->Fill(Kchboost[8]);
            hist_three[5][5]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][5]->Fill(Kchrec[5]);
            hist_three[7][5]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 7)
         {
            hist_three[0][6]->Fill(DeltaT);
            hist_three[1][6]->Fill(k_pathpm);
            hist_three[2][6]->Fill(k_path00);
            hist_three[3][6]->Fill(fourKnetri[8]);
            hist_three[4][6]->Fill(Kchboost[8]);
            hist_three[5][6]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][6]->Fill(Kchrec[5]);
            hist_three[7][6]->Fill(fourKnetri[5]);
         }
	 if(*mcflag == 0)
         {
            hist_three[0][7]->Fill(DeltaT);
            hist_three[1][7]->Fill(k_pathpm);
            hist_three[2][7]->Fill(k_path00);
            hist_three[3][7]->Fill(fourKnetri[8]);
            hist_three[4][7]->Fill(Kchboost[8]);
            hist_three[5][7]->Fill(cos(M_PI**anglepipi_CM_kch/180.));
            hist_three[6][7]->Fill(Kchrec[5]);
            hist_three[7][7]->Fill(fourKnetri[5]);
	 }
      }
   }

   //Pipi selection
   if(donepipi[0] == 1 && donepipi[1] == 1 /*&& abs(Kchrec1[5] - PhysicsConstants::mK0) < 5*/ && *done4 == 1)
   {
      hist_pipi_nocuts->Fill(DeltaT, br_ks_pippim + br_kl_pippim);

      if(1)//tot_cuts_charged)
      {
         hist_pipi_withcuts->Fill(DeltaT);

         if(*mctruth == 1)
         {
            hist_pipi[0][0]->Fill(DeltaT);
            hist_pipi[1][0]->Fill(k_pathpm);
            hist_pipi[2][0]->Fill(k_path00);
            hist_pipi[3][0]->Fill(fourKnetri[8]);
            hist_pipi[4][0]->Fill(Kchboost[8]);
            hist_pipi[5][0]->Fill(trcv_sum);
            hist_pipi[6][0]->Fill(Kchrec[5]);
            hist_pipi[7][0]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 3)
         {
            hist_pipi[0][1]->Fill(DeltaT);
            hist_pipi[1][1]->Fill(k_pathpm);
            hist_pipi[2][1]->Fill(k_path00);
            hist_pipi[3][1]->Fill(fourKnetri[8]);
            hist_pipi[4][1]->Fill(Kchboost[8]);
            hist_pipi[5][1]->Fill(trcv_sum);
            hist_pipi[6][1]->Fill(Kchrec[5]);
            hist_pipi[7][1]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 4)
         {
            hist_pipi[0][2]->Fill(DeltaT);
            hist_pipi[1][2]->Fill(k_pathpm);
            hist_pipi[2][2]->Fill(k_path00);
            hist_pipi[3][2]->Fill(fourKnetri[8]);
            hist_pipi[4][2]->Fill(Kchboost[8]);
            hist_pipi[5][2]->Fill(trcv_sum);
            hist_pipi[6][2]->Fill(Kchrec[5]);
            hist_pipi[7][2]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 5)
         {
            hist_pipi[0][3]->Fill(DeltaT);
            hist_pipi[1][3]->Fill(k_pathpm);
            hist_pipi[2][3]->Fill(k_path00);
            hist_pipi[3][3]->Fill(fourKnetri[8]);
            hist_pipi[4][3]->Fill(Kchboost[8]);
            hist_pipi[5][3]->Fill(trcv_sum);
            hist_pipi[6][3]->Fill(Kchrec[5]);
            hist_pipi[7][3]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 6)
         {
            hist_pipi[0][4]->Fill(DeltaT);
            hist_pipi[1][4]->Fill(k_pathpm);
            hist_pipi[2][4]->Fill(k_path00);
            hist_pipi[3][4]->Fill(fourKnetri[8]);
            hist_pipi[4][4]->Fill(Kchboost[8]);
            hist_pipi[5][4]->Fill(trcv_sum);
            hist_pipi[6][4]->Fill(Kchrec[5]);
            hist_pipi[7][4]->Fill(fourKnetri[5]);
         }
         if(*mctruth_pipi == 1)
         {
            hist_pipi[0][5]->Fill(DeltaT);
            hist_pipi[1][5]->Fill(k_pathpm);
            hist_pipi[2][5]->Fill(k_path00);
            hist_pipi[3][5]->Fill(fourKnetri[8]);
            hist_pipi[4][5]->Fill(Kchboost[8]);
            hist_pipi[5][5]->Fill(trcv_sum);
            hist_pipi[6][5]->Fill(Kchrec[5]);
            hist_pipi[7][5]->Fill(fourKnetri[5]);
         }
         if(*mctruth == 7)
         {
            hist_pipi[0][6]->Fill(DeltaT);
            hist_pipi[1][6]->Fill(k_pathpm);
            hist_pipi[2][6]->Fill(k_path00);
            hist_pipi[3][6]->Fill(fourKnetri[8]);
            hist_pipi[4][6]->Fill(Kchboost[8]);
            hist_pipi[5][6]->Fill(trcv_sum);
            hist_pipi[6][6]->Fill(Kchrec[5]);
            hist_pipi[7][6]->Fill(fourKnetri[5]);
         }
	 if(*mcflag == 0)
         {
            hist_pipi[0][7]->Fill(DeltaT);
            hist_pipi[1][7]->Fill(k_pathpm);
            hist_pipi[2][7]->Fill(k_path00);
            hist_pipi[3][7]->Fill(fourKnetri[8]);
            hist_pipi[4][7]->Fill(Kchboost[8]);
            hist_pipi[5][7]->Fill(trcv_sum);
            hist_pipi[6][7]->Fill(Kchrec[5]);
            hist_pipi[7][7]->Fill(fourKnetri[5]);
	 }
      }
   }
   return kTRUE;
}

void histos::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void histos::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   //canva1d_semi->SetLogy(1);
   //

   TObjArray *mc_semi[10];
   TObjArray *mc_three[10];
   TObjArray *mc_pipi[10];

   for(Int_t i = 0; i < 8; i++)
   {
	mc_three[i] = new TObjArray(KLOE::channNum-2);
	mc_pipi[i] = new TObjArray(KLOE::channNum-2);
	mc_semi[i] = new TObjArray(KLOE::channNum-2);
   }

   for(Int_t i = 0; i < 8; i++)
	   for(Int_t j = 0; j < KLOE::channNum - 2; j++) mc_semi[i]->Add(hist_semi[i][j]);

   for(Int_t i = 0; i < 8; i++)
	   for(Int_t j = 0; j < KLOE::channNum - 2; j++) mc_three[i]->Add(hist_three[i][j]);

   for(Int_t i = 0; i < 8; i++)
	   for(Int_t j = 0; j < KLOE::channNum - 2; j++) mc_pipi[i]->Add(hist_pipi[i][j]);

   TFractionFitter* fit_semi[10];
   TFractionFitter* fit_three[10];
   TFractionFitter* fit_pipi[10];

   for(Int_t i = 0; i < 8; i++)
   {
	fit_semi[i] = new TFractionFitter(hist_semi[i][KLOE::channNum-1-1],mc_semi[i]);
	fit_three[i] = new TFractionFitter(hist_three[i][KLOE::channNum-1-1],mc_three[i]);
	fit_pipi[i] = new TFractionFitter(hist_pipi[i][KLOE::channNum-1-1],mc_pipi[i]);

   	for(Int_t j = 0; j < KLOE::channNum - 2; j++)
   	{
		fit_semi[i]->Constrain(j,0.0,1.0);
		fit_three[i]->Constrain(j,0.0,1.0);
		fit_pipi[i]->Constrain(j,0.0,1.0);
	}

	fit_semi[i]->Fit();
	fit_three[i]->Fit();
	fit_pipi[i]->Fit();

	hist_semi[i][KLOE::channNum - 1] = (TH1F*) fit_semi[i]->GetPlot();
	hist_three[i][KLOE::channNum - 1] = (TH1F*) fit_three[i]->GetPlot();
	hist_pipi[i][KLOE::channNum - 1] = (TH1F*) fit_pipi[i]->GetPlot();
   }

   gStyle->SetStatX(0.9);
   gStyle->SetStatY(0.9);
   gStyle->SetStatH(0.15);
   gStyle->SetStatW(0.25);

   Int_t index_sort_semi[10][KLOE::channNum], index_sort_three[10][KLOE::channNum], index_sort_pipi[10][KLOE::channNum];
   Double_t max_counts_semi[10][KLOE::channNum], max_counts_three[10][KLOE::channNum], max_counts_pipi[10][KLOE::channNum];
   TString title_semi[100], title_three[100], title_pipi[100];

   title_semi[0] = "#Deltat [#tau_{S}]";
   title_semi[1] = "l_{K#rightarrow#pi^{+}#pi^{-}} [cm]";
   title_semi[2] = "l_{K#rightarrow#pi^{0}#pi^{0}} [cm]";
   title_semi[3] = "z_{neu} [cm]";
   title_semi[4] = "z_{ch} [cm]";
   title_semi[5] = "#sum_{i=1}^{4}t_{i,r} [ns]";
   title_semi[6] = "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]";
   title_semi[7] = "m^{inv}_{4#gamma} [MeV/c^{2}]";

   title_three[0] = "#Deltat [#tau_{S}]";
   title_three[1] = "l_{K#rightarrow#pi^{+}#pi^{-}} [cm]";
   title_three[2] = "l_{K#rightarrow#pi^{0}#pi^{0}} [cm]";
   title_three[3] = "z_{neu} [cm]";
   title_three[4] = "z_{ch} [cm]";
   title_three[5] = "cos(#alpha^{CM}_{#pi^{+},#pi^{-}})";
   title_three[6] = "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]";
   title_three[7] = "m^{inv}_{4#gamma} [MeV/c^{2}]";

   title_pipi[0] = "#Deltat [#tau_{S}]";
   title_pipi[1] = "l_{K#rightarrow#pi^{+}#pi^{-}} [cm]";
   title_pipi[2] = "l_{K#rightarrow#pi^{0}#pi^{0}} [cm]";
   title_pipi[3] = "z_{neu} [cm]";
   title_pipi[4] = "z_{ch} [cm]";
   title_pipi[5] = "#sum_{i=1}^{4}t_{i,r} [ns]";
   title_pipi[6] = "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]";
   title_pipi[7] = "m^{inv}_{4#gamma} [MeV/c^{2}]";

   Float_t y1d_semi[10][2], y1d_three[10][2], y1d_pipi[10][2], x2d[2], y2d[2]; 

   for(Int_t i = 0; i < 8; i++)
      for(Int_t j = 0; j < KLOE::channNum; j++) max_counts_semi[i][j] = hist_semi[i][j]->GetMaximum();

   for(Int_t i = 0; i < 8; i++)
      for(Int_t j = 0; j < KLOE::channNum; j++) max_counts_three[i][j] = hist_three[i][j]->GetMaximum();
   
   for(Int_t i = 0; i < 8; i++)
      for(Int_t j = 0; j < KLOE::channNum; j++) max_counts_pipi[i][j] = hist_pipi[i][j]->GetMaximum();

   for(Int_t i = 0; i < 8; i++)
   {
      TMath::Sort(KLOE::channNum, max_counts_semi[i], index_sort_semi[i]);
      y1d_semi[i][0] = 0;
      y1d_semi[i][1] = max_counts_semi[i][index_sort_semi[i][0]] + 200.;
   }

   for(Int_t i = 0; i < 8; i++)
   {
      TMath::Sort(KLOE::channNum, max_counts_three[i], index_sort_three[i]);
      y1d_three[i][0] = 0;
      y1d_three[i][1] = max_counts_three[i][index_sort_three[i][0]] + 200.;
   }

   for(Int_t i = 0; i < 8; i++)
   {
      TMath::Sort(KLOE::channNum, max_counts_pipi[i], index_sort_pipi[i]);
      y1d_pipi[i][0] = 0;
      y1d_pipi[i][1] = max_counts_pipi[i][index_sort_pipi[i][0]] + 200.;
   }

   for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < KLOE::channNum; j++)
      {
         canva1d_semi[i]->SetLeftMargin(0.15);
         canva1d_semi[i]->SetBottomMargin(0.15);
         canva1d_semi[i]->cd();

         hist_semi[i][j]->SetLineColor(channColor[j]);

            if(j == 0)
            {
               hist_semi[i][index_sort_semi[i][j]]->GetXaxis()->SetMaxDigits(3);
               hist_semi[i][index_sort_semi[i][j]]->GetYaxis()->SetMaxDigits(3);

               hist_semi[i][index_sort_semi[i][j]]->GetXaxis()->SetTitle(title_semi[i]);
               hist_semi[i][index_sort_semi[i][j]]->GetYaxis()->SetTitle("Counts");

               hist_semi[i][index_sort_semi[i][j]]->GetXaxis()->CenterTitle(1);
               hist_semi[i][index_sort_semi[i][j]]->GetYaxis()->CenterTitle(1);

               hist_semi[i][index_sort_semi[i][j]]->GetYaxis()->SetRangeUser(y1d_semi[i][0], y1d_semi[i][1]);
               hist_semi[i][index_sort_semi[i][j]]->Draw();
            }
            else
            {
               hist_semi[i][index_sort_semi[i][j]]->Draw("SAME");
            }
         
      }
      legend[i]->Draw();
   }

   for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < KLOE::channNum; j++)
      {
         canva1d_three[i]->SetLeftMargin(0.15);
         canva1d_three[i]->SetBottomMargin(0.15);
         canva1d_three[i]->cd();

         hist_three[i][j]->SetLineColor(channColor[j]);

            if(j == 0)
            {
               hist_three[i][index_sort_three[i][j]]->GetXaxis()->SetMaxDigits(3);
               hist_three[i][index_sort_three[i][j]]->GetYaxis()->SetMaxDigits(3);

               hist_three[i][index_sort_three[i][j]]->GetXaxis()->SetTitle(title_three[i]);
               hist_three[i][index_sort_three[i][j]]->GetYaxis()->SetTitle("Counts");

               hist_three[i][index_sort_three[i][j]]->GetXaxis()->CenterTitle(1);
               hist_three[i][index_sort_three[i][j]]->GetYaxis()->CenterTitle(1);

               hist_three[i][index_sort_three[i][j]]->GetYaxis()->SetRangeUser(y1d_three[i][0], y1d_three[i][1]);
               hist_three[i][index_sort_three[i][j]]->Draw();
            }
            else
            {
               hist_three[i][index_sort_three[i][j]]->Draw("SAME");
            }
         
      }
      legend[i]->Draw();
   }

   for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < KLOE::channNum; j++)
      {
         canva1d_pipi[i]->SetLeftMargin(0.15);
         canva1d_pipi[i]->SetBottomMargin(0.15);
         canva1d_pipi[i]->cd();

         hist_pipi[i][j]->SetLineColor(channColor[j]);

            if(j == 0)
            {
               hist_pipi[i][index_sort_pipi[i][j]]->GetXaxis()->SetMaxDigits(3);
               hist_pipi[i][index_sort_pipi[i][j]]->GetYaxis()->SetMaxDigits(3);

               hist_pipi[i][index_sort_pipi[i][j]]->GetXaxis()->SetTitle(title_pipi[i]);
               hist_pipi[i][index_sort_pipi[i][j]]->GetYaxis()->SetTitle("Counts");

               hist_pipi[i][index_sort_pipi[i][j]]->GetXaxis()->CenterTitle(1);
               hist_pipi[i][index_sort_pipi[i][j]]->GetYaxis()->CenterTitle(1);

               hist_pipi[i][index_sort_pipi[i][j]]->GetYaxis()->SetRangeUser(y1d_pipi[i][0], y1d_pipi[i][1]);
               hist_pipi[i][index_sort_pipi[i][j]]->Draw();
            }
            else
            {
               hist_pipi[i][index_sort_pipi[i][j]]->Draw("SAME");
            }
         
      }
      legend[i]->Draw();
   }

   /*for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < KLOE::channNum-2; j++)
      {
         canva2d[i][j]->SetLeftMargin(0.15);
         canva2d[i][j]->SetRightMargin(0.15);
         canva2d[i][j]->SetBottomMargin(0.15);
         canva2d[i][j]->cd();

         hist2d[i][j]->GetXaxis()->SetMaxDigits(3);
         hist2d[i][j]->GetYaxis()->SetMaxDigits(3);
         hist2d[i][j]->GetZaxis()->SetMaxDigits(3);

         hist2d[i][j]->GetXaxis()->SetTitle(title[3]);
         hist2d[i][j]->GetYaxis()->SetTitle(title[1]);

         hist2d[i][j]->GetXaxis()->CenterTitle(1);
         hist2d[i][j]->GetYaxis()->CenterTitle(1);

         //hist2d[0][j]->GetXaxis()->SetRangeUser(x2d[0], x2d[1]);
         //hist2d[0][j]->GetYaxis()->SetRangeUser(y2d[0], y2d[1]);
         hist2d[i][j]->Draw("COLZ");
      }
   }*/

   canva1d_seminocuts->SetLeftMargin(0.15);
   canva1d_seminocuts->SetBottomMargin(0.15);
   canva1d_seminocuts->cd();

   hist_semi_nocuts->Draw("HIST");

   canva1d_threenocuts->SetLeftMargin(0.15);
   canva1d_threenocuts->SetBottomMargin(0.15);
   canva1d_threenocuts->cd();

   hist_three_nocuts->Draw("HIST");

   canva1d_pipinocuts->SetLeftMargin(0.15);
   canva1d_pipinocuts->SetBottomMargin(0.15);
   canva1d_pipinocuts->cd();

   hist_pipi_nocuts->Draw("HIST");

   //hist_semi_nocuts->Add(hist_three_nocuts);
   //hist_semi_nocuts->Add(hist_pipi_nocuts);

   //hist_semi_withcuts->Add(hist_three_withcuts);
   //hist_semi_withcuts->Add(hist_pipi_withcuts);
 
   hist_semi_withcuts->Divide(hist_semi_nocuts);

   auto gr = new TGraphAsymmErrors(hist_three_withcuts, hist_three_nocuts);

   canva_efficiency->SetLeftMargin(0.15);
   canva_efficiency->SetBottomMargin(0.15);
   canva_efficiency->cd();

   gr->GetXaxis()->SetLimits(-100.0, 100.0);
   gr->GetHistogram()->SetMaximum(1.0);
   gr->GetHistogram()->SetMinimum(0.0);
   //gr->Draw("AP");

   hist_semi_withcuts->Draw("HIST");

   auto gr1 = new TGraphAsymmErrors();
   gr1->Divide(hist_signal, hist_signal_nocuts);

   canva_eff_signal->SetLeftMargin(0.15);
   canva_eff_signal->SetBottomMargin(0.15);
   canva_eff_signal->cd();

   //gr1->GetHistogram()->SetMaximum(0.05);
   //gr1->GetHistogram()->SetMinimum(0.0);
   hist_signal->SetLineColor(kRed);
   hist_signal->Draw("HIST");

   canva1d_semi[0]->Print("plots/Semi/deltat_semi_signal_cuts.png");
   canva1d_semi[1]->Print("plots/Semi/kpathch_semi_signal_cuts.png");
   canva1d_semi[2]->Print("plots/Semi/kpathneu_semi_signal_cuts.png");
   canva1d_semi[3]->Print("plots/Semi/zcoorneu_semi_signal_cuts.png");
   canva1d_semi[4]->Print("plots/Semi/zcoorch_semi_signal_cuts.png");
   canva1d_semi[5]->Print("plots/Semi/trcv_semi_signal_cuts.png");
   canva1d_semi[6]->Print("plots/Semi/minvch_semi_signal_cuts.png");
   canva1d_semi[7]->Print("plots/Semi/minvneu_semi_signal_cuts.png");

   canva1d_three[0]->Print("plots/Three/deltat_three_signal_cuts.png");
   canva1d_three[1]->Print("plots/Three/kpathch_three_signal_cuts.png");
   canva1d_three[2]->Print("plots/Three/kpathneu_three_signal_cuts.png");
   canva1d_three[3]->Print("plots/Three/zcoorneu_three_signal_cuts.png");
   canva1d_three[4]->Print("plots/Three/zcoorch_three_signal_cuts.png");
   canva1d_three[5]->Print("plots/Three/trcv_three_signal_cuts.png");
   canva1d_three[6]->Print("plots/Three/minvch_three_signal_cuts.png");
   canva1d_three[7]->Print("plots/Three/minvneu_three_signal_cuts.png");

   canva1d_pipi[0]->Print("plots/Pipi/deltat_pipi_signal_cuts.png");
   canva1d_pipi[1]->Print("plots/Pipi/kpathch_pipi_signal_cuts.png");
   canva1d_pipi[2]->Print("plots/Pipi/kpathneu_pipi_signal_cuts.png");
   canva1d_pipi[3]->Print("plots/Pipi/zcoorneu_pipi_signal_cuts.png");
   canva1d_pipi[4]->Print("plots/Pipi/zcoorch_pipi_signal_cuts.png");
   canva1d_pipi[5]->Print("plots/Pipi/trcv_pipi_signal_cuts.png");
   canva1d_pipi[6]->Print("plots/Pipi/minvch_pipi_signal_cuts.png");
   canva1d_pipi[7]->Print("plots/Pipi/minvneu_pipi_signal_cuts.png");

   canva1d_seminocuts->Print("plots/Semi/deltat_semi_nocuts.png");
   canva1d_threenocuts->Print("plots/Three/deltat_three_nocuts.png");
   canva1d_pipinocuts->Print("plots/Pipi/deltat_pipi_nocuts.png");
   canva_efficiency->Print("efficiency_mc_cs.png");

   canva_eff_signal->Print("efficiency_mc.png");
}
