#define full_analysis_cxx
// The class definition in full_analysis.h has been generated automatically
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
// root> T->Process("full_analysis.C")
// root> T->Process("full_analysis.C","some options")
// root> T->Process("full_analysis.C+")
//

#include "full_analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TFractionFitter.h>
#include <TMath.h>
#include <TLine.h>
#include <TBox.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include "../Include/const.h"
#include "../Include/Codes/interf_function.h"

TCanvas *canva1d_semi[10], *canva1d_three[10], *canva1d_pipi[10], *canva2d[10][channNum], *canvaproj[10][channNum];
TCanvas *canva1d_seminocuts, *canva1d_threenocuts, *canva1d_pipinocuts;
TCanvas *canva1d_semiwithcuts, *canva1d_threewithcuts, *canva1d_pipiwithcuts;
TCanvas *canva_efficiency, *canva_eff_signal;

TH1 *hist_semi[10][channNum], *hist_three[10][channNum], *hist_pipi[10][channNum];
TH1 *hist_signal, *hist_signal_nocuts;
TH1 *hist_semi_nocuts, *hist_three_nocuts, *hist_pipi_nocuts;
TH1 *hist_semi_withcuts, *hist_three_withcuts, *hist_pipi_withcuts;
TH2 *hist2d[10][channNum];

TLegend *legend[10];

TH1 *signal_before, *signal_after;

TEfficiency *pEff_signal;
TEfficiency *pEff_signal_tri;
TEfficiency *pEff_semi;
TEfficiency *pEff_three;
TEfficiency *pEff_pipi;

TEfficiency *pEff_three_length;
TEfficiency *pEff_pipi_length;

TEfficiency *pEff_semimc;
TEfficiency *pEff_threemc;
TEfficiency *pEff_pipimc;

Int_t entries[channNum] = {0}, entries_sel[3][channNum] = {};

Double_t tot_br_semi = ((br_kl_piele + br_kl_pimu)/(br_kl_piele + br_kl_pimu + br_ks_piele + br_ks_pimu))*br_ks_pi0pi0 + ((br_ks_piele + br_ks_pimu)/(br_kl_piele + br_kl_pimu + br_ks_piele + br_ks_pimu))*br_kl_pi0pi0;
Double_t tot_br_pipi = ((br_kl_pippim)/(br_kl_pippim + br_ks_pippim))*br_ks_pippim + ((br_ks_pippim)/(br_kl_pippim + br_ks_pippim))*br_kl_pippim;

void full_analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   for(Int_t i = 0; i < 10; i++) canva1d_semi[i] = new TCanvas(("canva1d_semi" + std::to_string(i)).c_str(), "", 750, 750);
   for(Int_t i = 0; i < 10; i++) canva1d_three[i] = new TCanvas(("canva1d_three" + std::to_string(i)).c_str(), "", 750, 750);
   for(Int_t i = 0; i < 10; i++) canva1d_pipi[i] = new TCanvas(("canva1d_pipi" + std::to_string(i)).c_str(), "", 750, 750);

   canva_efficiency = new TCanvas("Canva efficiency", "", 750, 750);
   canva_eff_signal = new TCanvas("Canva eff signal", "", 750, 750);
   canva1d_seminocuts = new TCanvas("Canva semi", "", 750, 750);
   canva1d_threenocuts = new TCanvas("Canva three", "", 750, 750);
   canva1d_pipinocuts = new TCanvas("Canva pipi", "", 750, 750);
   canva1d_semiwithcuts = new TCanvas("Canva semi cuts", "", 750, 750);
   canva1d_threewithcuts = new TCanvas("Canva three cuts", "", 750, 750);
   canva1d_pipiwithcuts = new TCanvas("Canva pipi cuts", "", 750, 750);

  // for(Int_t i = 0; i < 10; i++)
  //    for(Int_t j = 0; j < channNum-2; j++) canva2d[i][j] = new TCanvas((channName[j] + " canva2d" + to_string(i)).c_str(), "", 750, 750);
  // for(Int_t i = 0; i < 10; i++)
  //    for(Int_t j = 0; j < channNum-2; j++) canvaproj[i][j] = new TCanvas((channName[j] + " canvaproj" + to_string(i)).c_str(), "", 750, 750);

   for(Int_t j = 0; j < channNum; j++) hist_semi[0][j] = new TH1F(channName[j] + "0semi", "", 91, -90, 90);
   for(Int_t j = 0; j < channNum; j++) hist_semi[1][j] = new TH1F(channName[j] + "1semi", "", 101, 0, 100);
   for(Int_t j = 0; j < channNum; j++) hist_semi[2][j] = new TH1F(channName[j] + "2semi", "", 101, 0, 100);
   for(Int_t j = 0; j < channNum; j++) hist_semi[3][j] = new TH1F(channName[j] + "3semi", "", 300, -100, 100);
   for(Int_t j = 0; j < channNum; j++) hist_semi[4][j] = new TH1F(channName[j] + "4semi", "", 150, 0, 180);
   for(Int_t j = 0; j < channNum; j++) hist_semi[5][j] = new TH1F(channName[j] + "5semi", "", 150, 0, 220);
   for(Int_t j = 0; j < channNum; j++) hist_semi[6][j] = new TH1F(channName[j] + "6semi", "", 200, 300, 1000);
   for(Int_t j = 0; j < channNum; j++) hist_semi[7][j] = new TH1F(channName[j] + "7semi", "", 200, 300, 1000);

   for(Int_t j = 0; j < channNum; j++) hist_three[0][j] = new TH1F(channName[j] + "0three", "", 91, -90, 90);
   for(Int_t j = 0; j < channNum; j++) hist_three[1][j] = new TH1F(channName[j] + "1three", "", 100, 0, 50);
   for(Int_t j = 0; j < channNum; j++) hist_three[2][j] = new TH1F(channName[j] + "2three", "", 100, 0, 50);
   for(Int_t j = 0; j < channNum; j++) hist_three[3][j] = new TH1F(channName[j] + "3three", "", 200, -10, 10);
   for(Int_t j = 0; j < channNum; j++) hist_three[4][j] = new TH1F(channName[j] + "4three", "", 200, 0, 7000);
   for(Int_t j = 0; j < channNum; j++) hist_three[5][j] = new TH1F(channName[j] + "5three", "", 30, -1.0, 0.0);
   for(Int_t j = 0; j < channNum; j++) hist_three[6][j] = new TH1F(channName[j] + "6three", "", 200, 300, 1000);
   for(Int_t j = 0; j < channNum; j++) hist_three[7][j] = new TH1F(channName[j] + "7three", "", 200, 300, 1000);

   for(Int_t j = 0; j < channNum; j++) hist_pipi[0][j] = new TH1F(channName[j] + "0pipi", "", 91, -90, 90);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[1][j] = new TH1F(channName[j] + "1pipi", "", 100, 0, 50);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[2][j] = new TH1F(channName[j] + "2pipi", "", 100, 0, 50);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[3][j] = new TH1F(channName[j] + "3pipi", "", 200, -10, 10);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[4][j] = new TH1F(channName[j] + "4pipi", "", 200, -10, 10);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[5][j] = new TH1F(channName[j] + "5pipi", "", 100, -100, 5);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[6][j] = new TH1F(channName[j] + "6pipi", "", 12, mK0 - 20, mK0 + 20);
   for(Int_t j = 0; j < channNum; j++) hist_pipi[7][j] = new TH1F(channName[j] + "7pipi", "", 200, 300, 1000);

   hist_signal = new TH1F("Signal histo", "", 100, -100, 10);
   hist_signal_nocuts = new TH1F("Signal histo nocuts", "", 91, -90, 90);
   hist_semi_nocuts = new TH1F("semi histo", "", 91, -90, 90);
   hist_three_nocuts = new TH1F("three histo", "", 91, -90, 90);
   hist_pipi_nocuts = new TH1F("pipi histo", "", 91, -90, 90);
   hist_semi_withcuts = new TH1F("semi histo cuts", "", 91, -90, 90);
   hist_three_withcuts = new TH1F("three histo cuts", "", 91, -90, 90);
   hist_pipi_withcuts = new TH1F("pipi histo cuts", "", 91, -90, 90);

   //   for(Int_t i = 0; i < 10; i++)
   //      for(Int_t j = 0; j < channNum; j++)
   //         hist2d[i][j] = new TH2F(channName[j] + " 2d" + to_string(i), "", 13, 0, 60, 100, -100, 300);

   pEff_signal = new TEfficiency("eff_signal",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);
   pEff_signal_tri = new TEfficiency("eff_signal_tri",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);
   pEff_semi = new TEfficiency("eff_semi",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);
   pEff_three = new TEfficiency("eff_three",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);
   pEff_pipi = new TEfficiency("eff_pipi",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);

   pEff_three_length = new TEfficiency("eff_three_length",";Length of path of Kaon [cm];Efficiency",100,0.,50.);
   pEff_pipi_length = new TEfficiency("eff_pipi_length",";Length of path of Kaon [cm];Efficiency",100,0.,50.);

   pEff_semimc = new TEfficiency("eff_semimc",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);
   pEff_threemc = new TEfficiency("eff_threemc",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);
   pEff_pipimc = new TEfficiency("eff_pipimc",";#Deltat [#tau_{S}];Efficiency",91,-90.,90.);

   signal_before = new TH1F("signal before", ";#Deltat [#tau_{S}];Counts", 21, -20, 20);
   signal_after = new TH1F("signal_after", ";#Deltat [#tau_{S}];Counts", 21, -20, 20);

   TString option = GetOption();
}

void full_analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t full_analysis::Process(Long64_t entry)
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

   //Lorentz vectors
   TLorentzVector mom_kaon_pm, mom_kaon_00_std, mom_kaon_00_tri, mom_kaon_s, mom_kaon_l, 
                  mom_kaon_pm_alt, mom_kaon_00_std_alt, mom_kaon_00_tri_alt, mom_kaon_s_alt, mom_kaon_l_alt, 
                  mom_kne_mc, mom_kch_mc, momPhi, momPhi_mc;
   TLorentzVector pos_kaon_pm, pos_kaon_00_std, pos_kaon_00_tri, pos_kaon_s, pos_kaon_l, pos_phi, pos_kne_mc, pos_kch_mc;

   TVector3 phi_boost, kaon_00_std_boost, kaon_00_tri_boost, kaon_s_boost, kaon_l_boost, kaon_pm_boost;

   //Kaon vars
   Float_t k_pathpm, k_betapm;
   Float_t k_path00_std, k_beta00_std;
   Float_t k_path00_tri, k_beta00_tri;
   Float_t k_pathks, k_betaks;
   Float_t k_betakl, k_pathkl;

   //trcv
   Float_t TRCV[4], trcv_sum_tri, trcv_sum_std;

   //Boolean cuts
   Bool_t cuts_signal_charged_cs[3], cuts_signal_neutral_cs[2], cuts_semi[2], cuts_three, cuts_pipi;
   Bool_t tot_cuts_semi, tot_cuts_neutral, tot_cuts_charged;

   //Auxilliary vars
   Float_t phi_vel, kaon_pm_vel, kaon_s_vel, kaon_l_vel, kaon_00_std_vel, kaon_00_tri_vel;

   Float_t tpm, t00_std, t00_tri, ts, tl, tpm_mc, t00_mc;

   Float_t DeltaT_signal, DeltaT_control, DeltaT_pipi, DeltaT_mc;

   Float_t trcv_sum, trcv_sum_signal;

   fReader.SetLocalEntry(entry);

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //Momentum and position of Phi in LAB
   momPhi_mc(0) = Kchmc[0] + Knemc[0];
   momPhi_mc(1) = Kchmc[1] + Knemc[1];
   momPhi_mc(2) = Kchmc[2] + Knemc[2];
   momPhi_mc(3) = Kchmc[3] + Knemc[3];

   //Momentum and position of Phi in LAB
   momPhi(0) = *Bpx;
   momPhi(1) = *Bpy;
   momPhi(2) = *Bpz;
   momPhi(3) = *Broots;

   phi_vel = cVel*(sqrt(pow(momPhi(0),2) + pow(momPhi(1),2) + pow(momPhi(2),2))/momPhi(3));

   pos_phi(0) = 0.;
   pos_phi(1) = 0.;
   pos_phi(2) = 0.;
   pos_phi(3) = cVel*(sqrt(pow(pos_phi(0),2) + pow(pos_phi(1),2) + pow(pos_phi(2),2))/phi_vel);

   //Momentum and position of charged kaon MC
   mom_kch_mc(0) = Kchmc[0];
   mom_kch_mc(1) = Kchmc[1];
   mom_kch_mc(2) = Kchmc[2];
   mom_kch_mc(3) = Kchmc[3];

   kaon_pm_vel = cVel*(sqrt(pow(mom_kch_mc(0),2) + pow(mom_kch_mc(1),2) + pow(mom_kch_mc(2),2))/mom_kch_mc(3));

   pos_kch_mc(0) = Kchmc[6] - ipmc[0];
   pos_kch_mc(1) = Kchmc[7] - ipmc[1];
   pos_kch_mc(2) = Kchmc[8] - ipmc[2];
   pos_kch_mc(3) = cVel*(sqrt(pow(pos_kch_mc(0),2) + pow(pos_kch_mc(1),2) + pow(pos_kch_mc(2),2))/kaon_pm_vel);

   //Momentum and position of neutral kaon MC
   mom_kne_mc(0) = Knemc[0];
   mom_kne_mc(1) = Knemc[1];
   mom_kne_mc(2) = Knemc[2];
   mom_kne_mc(3) = Knemc[3];

   kaon_pm_vel = cVel*(sqrt(pow(mom_kne_mc(0),2) + pow(mom_kne_mc(1),2) + pow(mom_kne_mc(2),2))/mom_kne_mc(3));

   pos_kne_mc(0) = Knemc[6] - ipmc[0];
   pos_kne_mc(1) = Knemc[7] - ipmc[1];
   pos_kne_mc(2) = Knemc[8] - ipmc[2];
   pos_kne_mc(3) = cVel*(sqrt(pow(pos_kne_mc(0),2) + pow(pos_kne_mc(1),2) + pow(pos_kne_mc(2),2))/kaon_pm_vel);

   //Momentum and position of charged kaon
   mom_kaon_pm(0) = Kchboost[0];
   mom_kaon_pm(1) = Kchboost[1];
   mom_kaon_pm(2) = Kchboost[2];
   mom_kaon_pm(3) = Kchboost[3];

   kaon_pm_vel = cVel*(sqrt(pow(mom_kaon_pm(0),2) + pow(mom_kaon_pm(1),2) + pow(mom_kaon_pm(2),2))/mom_kaon_pm(3));

   pos_kaon_pm(0) = Kchboost[6] - *Bx;
   pos_kaon_pm(1) = Kchboost[7] - *By;
   pos_kaon_pm(2) = Kchboost[8] - *Bz;
   pos_kaon_pm(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_pm(0),2) + pow(pos_kaon_pm(1),2) + pow(pos_kaon_pm(2),2))/kaon_pm_vel);

   //Momentum and position of charged kaon
   mom_kaon_pm_alt(0) = sqrt(pow(mom_kaon_pm(3),2) - pow(mK0,2))*(pos_kaon_pm(0) - pos_phi(0))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(1) = sqrt(pow(mom_kaon_pm(3),2) - pow(mK0,2))*(pos_kaon_pm(1) - pos_phi(1))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(2) = sqrt(pow(mom_kaon_pm(3),2) - pow(mK0,2))*(pos_kaon_pm(2) - pos_phi(2))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(3) = mom_kaon_pm(3);

   if(sqrt(pow(Kchrec1[6] - *Bx,2) + pow(Kchrec1[7] - *By,2) + pow(Kchrec1[8] - *Bz,2)) < sqrt(pow(Kchrec2[6] - *Bx,2) + pow(Kchrec2[7] - *By,2) + pow(Kchrec2[8] - *Bz,2)))
   {
      //Momentum and position of charged kaon_s
      mom_kaon_s(0) = Kchrec1[0];
      mom_kaon_s(1) = Kchrec1[1];
      mom_kaon_s(2) = Kchrec1[2];
      mom_kaon_s(3) = Kchrec1[3];

      kaon_s_vel = cVel*(sqrt(pow(mom_kaon_s(0),2) + pow(mom_kaon_s(1),2) + pow(mom_kaon_s(2),2))/mom_kaon_s(3));

      pos_kaon_s(0) = Kchrec1[6] - *Bx;
      pos_kaon_s(1) = Kchrec1[7] - *By;
      pos_kaon_s(2) = Kchrec1[8] - *Bz;
      pos_kaon_s(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_s(0),2) + pow(pos_kaon_s(1),2) + pow(pos_kaon_s(2),2))/kaon_s_vel);

      //Momentum and position of charged kaon_s alternative
      mom_kaon_s_alt(0) = sqrt(pow(mom_kaon_s(3),2) - pow(mK0,2))*(pos_kaon_s(0) - pos_phi(0))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));
      mom_kaon_s_alt(1) = sqrt(pow(mom_kaon_s(3),2) - pow(mK0,2))*(pos_kaon_s(1) - pos_phi(1))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(2) = sqrt(pow(mom_kaon_s(3),2) - pow(mK0,2))*(pos_kaon_s(2) - pos_phi(2))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(3) = mom_kaon_s(3);

      //Momentum and position of charged kaon_l
      mom_kaon_l(0) = Kchrec2[0];
      mom_kaon_l(1) = Kchrec2[1];
      mom_kaon_l(2) = Kchrec2[2];
      mom_kaon_l(3) = Kchrec2[3];

      kaon_l_vel = cVel*(sqrt(pow(mom_kaon_l(0),2) + pow(mom_kaon_l(1),2) + pow(mom_kaon_l(2),2))/mom_kaon_l(3));

      pos_kaon_l(0) = Kchrec2[6] - *Bx;
      pos_kaon_l(1) = Kchrec2[7] - *By;
      pos_kaon_l(2) = Kchrec2[8] - *Bz;
      pos_kaon_l(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_l(0),2) + pow(pos_kaon_l(1),2) + pow(pos_kaon_l(2),2))/kaon_l_vel);

      //Momentum and position of charged kaon_l alternative
      mom_kaon_l_alt(0) = sqrt(pow(mom_kaon_l(3),2) - pow(mK0,2))*(pos_kaon_l(0) - pos_phi(0))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));
      mom_kaon_l_alt(1) = sqrt(pow(mom_kaon_l(3),2) - pow(mK0,2))*(pos_kaon_l(1) - pos_phi(1))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(2) = sqrt(pow(mom_kaon_l(3),2) - pow(mK0,2))*(pos_kaon_l(2) - pos_phi(2))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(3) = mom_kaon_l(3);
   }
   else
   {
      //Momentum and position of charged kaon_s
      mom_kaon_s(0) = Kchrec2[0];
      mom_kaon_s(1) = Kchrec2[1];
      mom_kaon_s(2) = Kchrec2[2];
      mom_kaon_s(3) = Kchrec2[3];

      kaon_s_vel = cVel*(sqrt(pow(mom_kaon_s(0),2) + pow(mom_kaon_s(1),2) + pow(mom_kaon_s(2),2))/mom_kaon_s(3));

      pos_kaon_s(0) = Kchrec2[6] - *Bx;
      pos_kaon_s(1) = Kchrec2[7] - *By;
      pos_kaon_s(2) = Kchrec2[8] - *Bz;
      pos_kaon_s(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_s(0),2) + pow(pos_kaon_s(1),2) + pow(pos_kaon_s(2),2))/kaon_s_vel);

      //Momentum and position of charged kaon_s alternative
      mom_kaon_s_alt(0) = sqrt(pow(mom_kaon_s(3),2) - pow(mK0,2))*(pos_kaon_s(0) - pos_phi(0))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));
      mom_kaon_s_alt(1) = sqrt(pow(mom_kaon_s(3),2) - pow(mK0,2))*(pos_kaon_s(1) - pos_phi(1))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(2) = sqrt(pow(mom_kaon_s(3),2) - pow(mK0,2))*(pos_kaon_s(2) - pos_phi(2))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(3) = mom_kaon_s(3);

      //Momentum and position of charged kaon_l
      mom_kaon_l(0) = Kchrec1[0];
      mom_kaon_l(1) = Kchrec1[1];
      mom_kaon_l(2) = Kchrec1[2];
      mom_kaon_l(3) = Kchrec1[3];

      kaon_l_vel = cVel*(sqrt(pow(mom_kaon_l(0),2) + pow(mom_kaon_l(1),2) + pow(mom_kaon_l(2),2))/mom_kaon_l(3));

      pos_kaon_l(0) = Kchrec1[6] - *Bx;
      pos_kaon_l(1) = Kchrec1[7] - *By;
      pos_kaon_l(2) = Kchrec1[8] - *Bz;
      pos_kaon_l(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_l(0),2) + pow(pos_kaon_l(1),2) + pow(pos_kaon_l(2),2))/kaon_l_vel);

      //Momentum and position of charged kaon_l alternative
      mom_kaon_l_alt(0) = sqrt(pow(mom_kaon_l(3),2) - pow(mK0,2))*(pos_kaon_l(0) - pos_phi(0))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));
      mom_kaon_l_alt(1) = sqrt(pow(mom_kaon_l(3),2) - pow(mK0,2))*(pos_kaon_l(1) - pos_phi(1))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(2) = sqrt(pow(mom_kaon_l(3),2) - pow(mK0,2))*(pos_kaon_l(2) - pos_phi(2))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(3) = mom_kaon_l(3);
   }

   //Momentum and position of neutral standard kaon
   mom_kaon_00_std(0) = Knereclor[0];
   mom_kaon_00_std(1) = Knereclor[1];
   mom_kaon_00_std(2) = Knereclor[2];
   mom_kaon_00_std(3) = Knereclor[3];

   kaon_00_std_vel = cVel*(sqrt(pow(mom_kaon_00_std(0),2) + pow(mom_kaon_00_std(1),2) + pow(mom_kaon_00_std(2),2))/mom_kaon_00_std(3));

   pos_kaon_00_std(0) = Knereclor[6] - *Bx;
   pos_kaon_00_std(1) = Knereclor[7] - *By;
   pos_kaon_00_std(2) = Knereclor[8] - *Bz;
   pos_kaon_00_std(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_00_std(0),2) + pow(pos_kaon_00_std(1),2) + pow(pos_kaon_00_std(2),2))/kaon_00_std_vel);

   //Momentum and position of neutral standard kaon
   mom_kaon_00_std_alt(0) = sqrt(pow(mom_kaon_00_std(3),2) - pow(mK0,2))*(pos_kaon_00_std(0) - pos_phi(0))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));
   mom_kaon_00_std_alt(1) = sqrt(pow(mom_kaon_00_std(3),2) - pow(mK0,2))*(pos_kaon_00_std(1) - pos_phi(1))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));;
   mom_kaon_00_std_alt(2) = sqrt(pow(mom_kaon_00_std(3),2) - pow(mK0,2))*(pos_kaon_00_std(2) - pos_phi(2))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));;
   mom_kaon_00_std_alt(3) = mom_kaon_00_std(3);

   //Momentum and position of neutral trilateration kaon
   mom_kaon_00_tri(0) = fourKnetri[0];
   mom_kaon_00_tri(1) = fourKnetri[1];
   mom_kaon_00_tri(2) = fourKnetri[2];
   mom_kaon_00_tri(3) = fourKnetri[3];

   kaon_00_tri_vel = cVel*mom_kaon_00_tri.BoostVector().Mag();

   pos_kaon_00_tri(0) = fourKnetri[6] - *Bx;
   pos_kaon_00_tri(1) = fourKnetri[7] - *By;
   pos_kaon_00_tri(2) = fourKnetri[8] - *Bz;
   pos_kaon_00_tri(3) = cVel*(pos_phi(3) + sqrt(pow(pos_kaon_00_tri(0),2) + pow(pos_kaon_00_tri(1),2) + pow(pos_kaon_00_tri(2),2))/kaon_00_tri_vel);

   //Momentum and position of neutral standard kaon
   mom_kaon_00_tri_alt(0) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(mK0,2))*(pos_kaon_00_tri(0) - pos_phi(0))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));
   mom_kaon_00_tri_alt(1) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(mK0,2))*(pos_kaon_00_tri(1) - pos_phi(1))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));;
   mom_kaon_00_tri_alt(2) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(mK0,2))*(pos_kaon_00_tri(2) - pos_phi(2))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));;
   mom_kaon_00_tri_alt(3) = mom_kaon_00_tri(3);
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Lorentz transformation to CM of Phi

   k_pathpm = sqrt(pow(pos_kaon_pm(0),2) + pow(pos_kaon_pm(1),2) + pow(pos_kaon_pm(2),2));
   k_betapm = sqrt(pow(mom_kaon_pm(0),2) + pow(mom_kaon_pm(1),2) + pow(mom_kaon_pm(2),2))/mom_kaon_pm(3);
   k_path00_tri = sqrt(pow(pos_kaon_00_tri(0),2) + pow(pos_kaon_00_tri(1),2) + pow(pos_kaon_00_tri(2),2));
   k_beta00_tri = sqrt(pow(mom_kaon_00_tri(0),2) + pow(mom_kaon_00_tri(1),2) + pow(mom_kaon_00_tri(2),2))/mom_kaon_00_tri(3);

   phi_boost = momPhi.BoostVector();

   pos_phi.Boost(-phi_boost);
   pos_kaon_pm.Boost(-phi_boost);
   pos_kaon_00_std.Boost(-phi_boost);
   pos_kaon_00_tri.Boost(-phi_boost);
   pos_kaon_s.Boost(-phi_boost);
   pos_kaon_l.Boost(-phi_boost);

   mom_kaon_pm.Boost(-phi_boost);
   mom_kaon_00_std.Boost(-phi_boost);
   mom_kaon_00_tri.Boost(-phi_boost);
   mom_kaon_s.Boost(-phi_boost);
   mom_kaon_l.Boost(-phi_boost);

   mom_kaon_pm_alt.Boost(-phi_boost);
   mom_kaon_00_std_alt.Boost(-phi_boost);
   mom_kaon_00_tri_alt.Boost(-phi_boost);
   mom_kaon_s_alt.Boost(-phi_boost);
   mom_kaon_l_alt.Boost(-phi_boost);

   momPhi.Boost(-phi_boost);

   k_pathpm = sqrt(pow(pos_kaon_pm(0),2) + pow(pos_kaon_pm(1),2) + pow(pos_kaon_pm(2),2));
   k_betapm = sqrt(pow(mom_kaon_pm(0),2) + pow(mom_kaon_pm(1),2) + pow(mom_kaon_pm(2),2))/mom_kaon_pm(3);

   phi_boost = momPhi_mc.BoostVector();
   mom_kne_mc.Boost(-phi_boost);
   mom_kch_mc.Boost(-phi_boost);
   pos_kne_mc.Boost(-phi_boost);
   pos_kch_mc.Boost(-phi_boost);

   kaon_pm_boost = mom_kaon_pm.BoostVector();
   kaon_00_std_boost = mom_kaon_00_std.BoostVector();
   kaon_00_tri_boost = mom_kaon_00_tri.BoostVector();
   kaon_s_boost = mom_kaon_s.BoostVector();
   kaon_l_boost = mom_kaon_l.BoostVector();

   pos_kaon_pm.Boost(-kaon_pm_boost);
   pos_kaon_00_std.Boost(-kaon_00_std_boost);
   pos_kaon_00_tri.Boost(-kaon_00_tri_boost);
   pos_kaon_s.Boost(-kaon_s_boost);
   pos_kaon_l.Boost(-kaon_l_boost);

   mom_kaon_pm.Boost(-kaon_pm_boost);
   mom_kaon_00_std.Boost(-kaon_00_std_boost);
   mom_kaon_00_tri.Boost(-kaon_00_tri_boost);
   mom_kaon_s.Boost(-kaon_s_boost);
   mom_kaon_l.Boost(-kaon_l_boost);

   kaon_pm_boost = mom_kch_mc.BoostVector();
   kaon_00_std_boost = mom_kne_mc.BoostVector();

   pos_kch_mc.Boost(-kaon_pm_boost);
   pos_kne_mc.Boost(-kaon_00_std_boost);

   //k_pathpm = sqrt(pow(pos_kaon_pm(0),2) + pow(pos_kaon_pm(1),2) + pow(pos_kaon_pm(2),2));
   //k_betapm = sqrt(pow(mom_kaon_pm(0),2) + pow(mom_kaon_pm(1),2) + pow(mom_kaon_pm(2),2))/mom_kaon_pm(3);
   //k_path00_tri = sqrt(pow(pos_kaon_00_std(0),2) + pow(pos_kaon_00_std(1),2) + pow(pos_kaon_00_std(2),2));
   //k_beta00_tri = sqrt(pow(mom_kaon_00_std(0),2) + pow(mom_kaon_00_std(1),2) + pow(mom_kaon_00_std(2),2))/mom_kaon_00_std(3);

   

   // Calculation of time difference

   tpm = pos_kaon_pm(3)/(tau_S_nonCPT*cVel);
   t00_std = pos_kaon_00_std(3)/(tau_S_nonCPT*cVel);
   t00_tri = pos_kaon_00_tri(3)/(tau_S_nonCPT*cVel);
   ts = pos_kaon_s(3)/(tau_S_nonCPT*cVel);
   tl = pos_kaon_l(3)/(tau_S_nonCPT*cVel);
   tpm_mc = pos_kch_mc(3)/(tau_S_nonCPT*cVel);
   t00_mc = pos_kne_mc(3)/(tau_S_nonCPT*cVel);  

   DeltaT_signal = (tpm - t00_std);
   DeltaT_control = (tpm - t00_tri);
   DeltaT_pipi = (tl - ts);
   DeltaT_mc = (tpm_mc - t00_mc);

   ///////////////////////////////////////////////////////////////////////////////////////////////

   //for(Int_t i = 0; i < 4; i++) TRCV[i] = TclOld[fourg4taken[i]] - (sqrt(pow(Xcl[fourg4taken[i]] - fourKnetri[6],2) + pow(Ycl[fourg4taken[i]] - fourKnetri[7],2) + pow(Zcl[fourg4taken[i]] - fourKnetri[8],2))/cVel) - (k_path00/(k_beta00*cVel));

   //trcv_sum = TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3];

   
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   //Signal truth check
   ///////////////////////////////////////////////////////////////////////////////////////////////////

   for(Int_t i = 0; i < 4; i++) TRCV[i] = TclOld[fourg4taken[i]] - (sqrt(pow(Xcl[fourg4taken[i]] - fourKnetri[6],2) + pow(Ycl[fourg4taken[i]] - fourKnetri[7],2) + pow(Zcl[fourg4taken[i]] - fourKnetri[8],2))/cVel) - (k_path00_tri/(k_beta00_tri*cVel));
   trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

   if(*mctruth == 1)
	{
      trcv_sum_signal = trcv[g4taken[0]-1] + trcv[g4taken[1]-1] + trcv[g4taken[2]-1] + trcv[g4taken[3]-1];
      pEff_signal->FillWeighted(trcv_sum_signal > -1 && abs(*minv4gam - mK0) < 76 && abs(Kchrec[5] - mK0) < 1.2 && *Qmiss_inv < 3.75 && cos(M_PI*(*anglepipi_CM_kch/180.)) < -0.8, interf_function(DeltaT_mc,1,0), DeltaT_signal);
      pEff_signal_tri->FillWeighted(trcv_sum > -1 && abs(fourKnetri[5] - mK0) < 76 && abs(Kchrec[5] - mK0) < 1.2 && *Qmiss_inv < 3.75 && cos(M_PI*(*anglepipi_CM_kch/180.)) < -0.8, interf_function(DeltaT_mc,1,0), DeltaT_signal);

      signal_before->Fill(DeltaT_signal, interf_function(DeltaT_mc,1,0));

      if(trcv_sum_signal > -1 && abs(*minv4gam - mK0) < 76 && abs(Kchrec[5] - mK0) < 1.2 && *Qmiss_inv < 3.75 && cos(M_PI*(*anglepipi_CM_kch/180.)) < -0.8)
      {
         signal_after->Fill(DeltaT_signal, interf_function(DeltaT_mc,1,0));
      }

   }

   cuts_semi[0] = abs(*Qmiss_inv - 71.13) < 25;
   cuts_semi[1] = abs(*anglepipi_CM_kch - 145.8) < 10;

   cuts_signal_neutral_cs[0] = (abs(fourKnetri[5] - mK0) < 76);
   cuts_signal_neutral_cs[1] = (trcv_sum > -1.0);

   cuts_signal_charged_cs[0] = (abs(Kchrec[5] - mK0) < 1.2);
   cuts_signal_charged_cs[1] = (*Qmiss_inv < 3.75);
   cuts_signal_charged_cs[2] = (cos(M_PI*(*anglepipi_CM_kch)/180.) < -0.8);

   tot_cuts_semi = cuts_semi[0] && cuts_semi[1];
   tot_cuts_neutral = cuts_signal_neutral_cs[0] && cuts_signal_neutral_cs[1];
   tot_cuts_charged = cuts_signal_charged_cs[0] && cuts_signal_charged_cs[1] && cuts_signal_charged_cs[2];

   ///////////////////////////////////////////////////////////////////////////////////////////////////
   // Control samples selection
   ///////////////////////////////////////////////////////////////////////////////////////////////////

   if(*mctruth == 1) entries[0]++;
   if(*mctruth == 3) entries[1]++;
   if(*mctruth == 4) entries[2]++;
   if(*mctruth == 5) entries[3]++;
   if(*mctruth == 6) entries[4]++;
   if(*mctruth_pipi == 1) entries[5]++;
   if(*mctruth == 7) entries[6]++;
   if(*mcflag == 0) entries[7]++;

   if(tot_cuts_semi && *done4 == 1)
   {
      if(1/*tot_cuts_neutral*/)
      {
         if(*mctruth == 1)
         {
            hist_semi[0][0]->Fill(DeltaT_control,tot_br_semi*interf_function(DeltaT_mc,1,0));
            hist_semi[1][0]->Fill(tpm);
            hist_semi[2][0]->Fill(t00_tri);
            hist_semi[3][0]->Fill(trcv_sum);
            hist_semi[4][0]->Fill(*anglepipi_CM_kch);
            hist_semi[5][0]->Fill(*Qmiss_inv);
            hist_semi[6][0]->Fill(Kchrec[5]);
            hist_semi[7][0]->Fill(fourKnetri[5]);
	         entries_sel[0][0]++;
         }
         if(*mctruth == 3)
         {
            hist_semi[0][1]->Fill(DeltaT_control,tot_br_semi);
            hist_semi[1][1]->Fill(tpm);
            hist_semi[2][1]->Fill(t00_tri);
            hist_semi[3][1]->Fill(trcv_sum);
            hist_semi[4][1]->Fill(*anglepipi_CM_kch);
            hist_semi[5][1]->Fill(*Qmiss_inv);
            hist_semi[6][1]->Fill(Kchrec[5]);
            hist_semi[7][1]->Fill(fourKnetri[5]);
	         entries_sel[0][1]++;
         }
         if(*mctruth == 4)
         {
            hist_semi[0][2]->Fill(DeltaT_control,tot_br_semi);
            hist_semi[1][2]->Fill(tpm);
            hist_semi[2][2]->Fill(t00_tri);
            hist_semi[3][2]->Fill(trcv_sum);
            hist_semi[4][2]->Fill(*anglepipi_CM_kch);
            hist_semi[5][2]->Fill(*Qmiss_inv);
            hist_semi[6][2]->Fill(Kchrec[5]);
            hist_semi[7][2]->Fill(fourKnetri[5]);
	         entries_sel[0][2]++;
         }
         if(*mctruth == 5)
         {
            hist_semi[0][3]->Fill(DeltaT_control,tot_br_semi);
            hist_semi[1][3]->Fill(tpm);
            hist_semi[2][3]->Fill(t00_tri);
            hist_semi[3][3]->Fill(trcv_sum);
            hist_semi[4][3]->Fill(*anglepipi_CM_kch);
            hist_semi[5][3]->Fill(*Qmiss_inv);
            hist_semi[6][3]->Fill(Kchrec[5]);
            hist_semi[7][3]->Fill(fourKnetri[5]);
	         entries_sel[0][3]++;
         }
         if(*mctruth == 6)
         {
            hist_semi[0][4]->Fill(DeltaT_control,tot_br_semi);
            hist_semi[1][4]->Fill(tpm);
            hist_semi[2][4]->Fill(t00_tri);
            hist_semi[3][4]->Fill(trcv_sum);
            hist_semi[4][4]->Fill(*anglepipi_CM_kch);
            hist_semi[5][4]->Fill(*Qmiss_inv);
            hist_semi[6][4]->Fill(Kchrec[5]);
            hist_semi[7][4]->Fill(fourKnetri[5]);
	         entries_sel[0][4]++;
         }
         if(*mctruth_pipi == 1)
         {
            hist_semi[0][5]->Fill(DeltaT_control,tot_br_semi);
            hist_semi[1][5]->Fill(tpm);
            hist_semi[2][5]->Fill(t00_tri);
            hist_semi[3][5]->Fill(trcv_sum);
            hist_semi[4][5]->Fill(*anglepipi_CM_kch);
            hist_semi[5][5]->Fill(*Qmiss_inv);
            hist_semi[6][5]->Fill(Kchrec[5]);
            hist_semi[7][5]->Fill(fourKnetri[5]);
	         entries_sel[0][5]++;
         }
         if(*mctruth == 7)
         {
            hist_semi[0][6]->Fill(DeltaT_control,tot_br_semi);
            hist_semi[1][6]->Fill(tpm);
            hist_semi[2][6]->Fill(t00_tri);
            hist_semi[3][6]->Fill(trcv_sum);
            hist_semi[4][6]->Fill(*anglepipi_CM_kch);
            hist_semi[5][6]->Fill(*Qmiss_inv);
            hist_semi[6][6]->Fill(Kchrec[5]);
            hist_semi[7][6]->Fill(fourKnetri[5]);
	         entries_sel[0][6]++;
         }
         if(*mcflag == 0)
               {

                  hist_semi[0][7]->Fill(DeltaT_control,tot_br_semi);
                  hist_semi[1][7]->Fill(tpm);
                  hist_semi[2][7]->Fill(t00_tri);
                  hist_semi[3][7]->Fill(trcv_sum);
                  hist_semi[4][7]->Fill(*anglepipi_CM_kch);
                  hist_semi[5][7]->Fill(*Qmiss_inv);
                  hist_semi[6][7]->Fill(Kchrec[5]);
                  hist_semi[7][7]->Fill(fourKnetri[5]);
         }
      }
   }

   //Three selection
   if(*done == 1 && *totalerr < 2000 && *done4 == 1)
   {
      if(1)//tot_cuts_charged)
      {
         if(*mctruth == 1)
         {
            hist_three[0][0]->Fill(DeltaT_control);
            hist_three[1][0]->Fill(tpm);
            hist_three[2][0]->Fill(t00_tri);
            hist_three[3][0]->Fill(fourKnetri[8]);
            hist_three[4][0]->Fill(*totalerr);
            hist_three[5][0]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][0]->Fill(Kchrec[5]);
            hist_three[7][0]->Fill(fourKnetri[5]);
            entries_sel[1][0]++;
         }
         if(*mctruth == 3)
         {
            hist_three[0][1]->Fill(DeltaT_control);
            hist_three[1][1]->Fill(tpm);
            hist_three[2][1]->Fill(t00_tri);
            hist_three[3][1]->Fill(fourKnetri[8]);
            hist_three[4][1]->Fill(*totalerr);
            hist_three[5][1]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][1]->Fill(Kchrec[5]);
            hist_three[7][1]->Fill(fourKnetri[5]);
            entries_sel[1][1]++;
         }
         if(*mctruth == 4)
         {
            hist_three[0][2]->Fill(DeltaT_control);
            hist_three[1][2]->Fill(tpm);
            hist_three[2][2]->Fill(t00_tri);
            hist_three[3][2]->Fill(fourKnetri[8]);
            hist_three[4][2]->Fill(*totalerr);
            hist_three[5][2]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][2]->Fill(Kchrec[5]);
            hist_three[7][2]->Fill(fourKnetri[5]);
            entries_sel[1][2]++;
         }
         if(*mctruth == 5)
         {
            hist_three[0][3]->Fill(DeltaT_control);
            hist_three[1][3]->Fill(tpm);
            hist_three[2][3]->Fill(t00_tri);
            hist_three[3][3]->Fill(fourKnetri[8]);
            hist_three[4][3]->Fill(*totalerr);
            hist_three[5][3]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][3]->Fill(Kchrec[5]);
            hist_three[7][3]->Fill(fourKnetri[5]);
            entries_sel[1][3]++;
         }
         if(*mctruth == 6)
         {
            hist_three[0][4]->Fill(DeltaT_control);
            hist_three[1][4]->Fill(tpm);
            hist_three[2][4]->Fill(t00_tri);
            hist_three[3][4]->Fill(fourKnetri[8]);
            hist_three[4][4]->Fill(*totalerr);
            hist_three[5][4]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][4]->Fill(Kchrec[5]);
            hist_three[7][4]->Fill(fourKnetri[5]);
            entries_sel[1][4]++;
         }
         if(*mctruth_pipi == 1)
         {
            hist_three[0][5]->Fill(DeltaT_control);
            hist_three[1][5]->Fill(tpm);
            hist_three[2][5]->Fill(t00_tri);
            hist_three[3][5]->Fill(fourKnetri[8]);
            hist_three[4][5]->Fill(*totalerr);
            hist_three[5][5]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][5]->Fill(Kchrec[5]);
            hist_three[7][5]->Fill(fourKnetri[5]);
            entries_sel[1][5]++;
         }
         if(*mctruth == 7)
         {
            hist_three[0][6]->Fill(DeltaT_control);
            hist_three[1][6]->Fill(tpm);
            hist_three[2][6]->Fill(t00_tri);
            hist_three[3][6]->Fill(fourKnetri[8]);
            hist_three[4][6]->Fill(*totalerr);
            hist_three[5][6]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][6]->Fill(Kchrec[5]);
            hist_three[7][6]->Fill(fourKnetri[5]);
            entries_sel[1][6]++;
         }
	 if(*mcflag == 0)
         {
            hist_three[0][7]->Fill(DeltaT_control);
            hist_three[1][7]->Fill(tpm);
            hist_three[2][7]->Fill(t00_tri);
            hist_three[3][7]->Fill(fourKnetri[8]);
            hist_three[4][7]->Fill(*totalerr);
            hist_three[5][7]->Fill(cos(M_PI*(*anglepipi_CM_kch)/180.));
            hist_three[6][7]->Fill(Kchrec[5]);
            hist_three[7][7]->Fill(fourKnetri[5]);
	 }
      }
   }

   //Pipi selection
   if(donepipi[0] == 1 && donepipi[1] == 1)// && abs(Kchrec1[5] - mK0) < 2)
   {
      if(1)//tot_cuts_charged)
      {

         if(*mctruth == 1)
         {
            hist_pipi[0][0]->Fill(DeltaT_pipi);
            hist_pipi[1][0]->Fill(ts);
            hist_pipi[2][0]->Fill(tl);
            hist_pipi[3][0]->Fill(fourKnetri[8]);
            hist_pipi[4][0]->Fill(Kchboost[8]);
            hist_pipi[5][0]->Fill(trcv_sum);
            hist_pipi[6][0]->Fill(Kchrec1[5]);
            hist_pipi[7][0]->Fill(Kchrec2[5]);
            entries_sel[2][0]++;
         }
         if(*mctruth == 3)
         {
            hist_pipi[0][1]->Fill(DeltaT_pipi);
            hist_pipi[1][1]->Fill(ts);
            hist_pipi[2][1]->Fill(tl);
            hist_pipi[3][1]->Fill(fourKnetri[8]);
            hist_pipi[4][1]->Fill(Kchboost[8]);
            hist_pipi[5][1]->Fill(trcv_sum);
            hist_pipi[6][1]->Fill(Kchrec1[5]);
            hist_pipi[7][1]->Fill(Kchrec2[5]);
            entries_sel[2][1]++;
         }
         if(*mctruth == 4)
         {
            hist_pipi[0][2]->Fill(DeltaT_pipi);
            hist_pipi[1][2]->Fill(ts);
            hist_pipi[2][2]->Fill(tl);
            hist_pipi[3][2]->Fill(fourKnetri[8]);
            hist_pipi[4][2]->Fill(Kchboost[8]);
            hist_pipi[5][2]->Fill(trcv_sum);
            hist_pipi[6][2]->Fill(Kchrec1[5]);
            hist_pipi[7][2]->Fill(Kchrec2[5]);
            entries_sel[2][2]++;
         }
         if(*mctruth == 5)
         {
            hist_pipi[0][3]->Fill(DeltaT_pipi);
            hist_pipi[1][3]->Fill(ts);
            hist_pipi[2][3]->Fill(tl);
            hist_pipi[3][3]->Fill(fourKnetri[8]);
            hist_pipi[4][3]->Fill(Kchboost[8]);
            hist_pipi[5][3]->Fill(trcv_sum);
            hist_pipi[6][3]->Fill(Kchrec1[5]);
            hist_pipi[7][3]->Fill(Kchrec2[5]);
            entries_sel[2][3]++;
         }
         if(*mctruth == 6)
         {
            hist_pipi[0][4]->Fill(DeltaT_pipi);
            hist_pipi[1][4]->Fill(ts);
            hist_pipi[2][4]->Fill(tl);
            hist_pipi[3][4]->Fill(fourKnetri[8]);
            hist_pipi[4][4]->Fill(Kchboost[8]);
            hist_pipi[5][4]->Fill(trcv_sum);
            hist_pipi[6][4]->Fill(Kchrec1[5]);
            hist_pipi[7][4]->Fill(Kchrec2[5]);
            entries_sel[2][4]++;
         }
         if(*mctruth_pipi == 1)
         {
            hist_pipi[0][5]->Fill(DeltaT_pipi);
            hist_pipi[1][5]->Fill(ts);
            hist_pipi[2][5]->Fill(tl);
            hist_pipi[3][5]->Fill(fourKnetri[8]);
            hist_pipi[4][5]->Fill(Kchboost[8]);
            hist_pipi[5][5]->Fill(trcv_sum);
            hist_pipi[6][5]->Fill(Kchrec1[5]);
            hist_pipi[7][5]->Fill(Kchrec2[5]);
            entries_sel[2][5]++;
         }
         if(*mctruth == 7)
         {
            hist_pipi[0][6]->Fill(DeltaT_pipi);
            hist_pipi[1][6]->Fill(ts);
            hist_pipi[2][6]->Fill(tl);
            hist_pipi[3][6]->Fill(fourKnetri[8]);
            hist_pipi[4][6]->Fill(Kchboost[8]);
            hist_pipi[5][6]->Fill(trcv_sum);
            hist_pipi[6][6]->Fill(Kchrec1[5]);
            hist_pipi[7][6]->Fill(Kchrec2[5]);
            entries_sel[2][6]++;
         }
	 if(*mcflag == 0)
         {
            hist_pipi[0][7]->Fill(DeltaT_pipi);
            hist_pipi[1][7]->Fill(ts);
            hist_pipi[2][7]->Fill(tl);
            hist_pipi[3][7]->Fill(fourKnetri[8]);
            hist_pipi[4][7]->Fill(Kchboost[8]);
            hist_pipi[5][7]->Fill(trcv_sum);
            hist_pipi[6][7]->Fill(Kchrec1[5]);
            hist_pipi[7][7]->Fill(Kchrec2[5]);
	 }
      }
   }

   if(tot_cuts_semi && *mcflag == 1 && *mctruth != 0 && *done4 == 1) pEff_semimc->FillWeighted(tot_cuts_neutral, tot_br_semi, DeltaT_control);
   if(*done == 1 && *mcflag == 1 && *mctruth != 0 && *done4 == 1) pEff_threemc->FillWeighted(tot_cuts_charged, br_ks_pippim, DeltaT_control);
   if(donepipi[0] == 1 && donepipi[1] == 1 && abs(Kchrec1[5] - mK0) < 2 && *mcflag == 1 && *mctruth != 0) pEff_pipimc->FillWeighted(tot_cuts_charged, tot_br_pipi, DeltaT_pipi);

   if(tot_cuts_semi && *mcflag == 0 && *done4 == 1) pEff_semi->FillWeighted(tot_cuts_neutral, tot_br_semi, DeltaT_control);
   if(*done == 1 && *mcflag == 0 && *done4 == 1) pEff_three->FillWeighted(tot_cuts_charged, br_ks_pippim, DeltaT_control);
   if(donepipi[0] == 1 && donepipi[1] == 1 && abs(Kchrec1[5] - mK0) < 2 && *mcflag == 0) pEff_pipi->FillWeighted(tot_cuts_charged, tot_br_pipi, DeltaT_pipi);

   if(*done == 1 && *mcflag == 0 && *done4 == 1) pEff_three_length->FillWeighted(tot_cuts_charged, br_ks_pippim, sqrt(pow(Kchrec[6] - *Bx,2) + pow(Kchrec[7] - *By,2) + pow(Kchrec[8] - *Bz,2)));
   if(donepipi[0] == 1 && donepipi[1] == 1 && abs(Kchrec1[5] - mK0) < 2 && *mcflag == 0) pEff_pipi_length->FillWeighted(tot_cuts_charged, tot_br_pipi, sqrt(pow(Kchrec2[6] - *Bx,2) + pow(Kchrec2[7] - *By,2) + pow(Kchrec2[8] - *Bz,2)));

   return kTRUE;
}

void full_analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void full_analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TObjArray *mc_semi;
   TObjArray *mc_three;
   TObjArray *mc_pipi;


   hist_semi[0][0]->Scale(tot_br_semi*(Float_t)entries_sel[0][0]/hist_semi[0][0]->Integral(0,hist_semi[0][0]->GetNbinsX()+1));
   mc_three = new TObjArray(channNum-2);
   mc_pipi = new TObjArray(channNum-2);
   mc_semi = new TObjArray(channNum-3);
   
   std::cout << hist_semi[0][5]->GetEntries() << std::endl;
   
	mc_semi->Add(hist_semi[0][0]);
   mc_semi->Add(hist_semi[0][1]);
   mc_semi->Add(hist_semi[0][2]);
   mc_semi->Add(hist_semi[0][3]);
   mc_semi->Add(hist_semi[0][4]);
   mc_semi->Add(hist_semi[0][6]);

   /*for(Int_t i = 0; i < 8; i++)
	   for(Int_t j = 0; j < channNum - 2; j++) mc_three[i]->Add(hist_three[i][j]);

   for(Int_t i = 0; i < 8; i++)
	   for(Int_t j = 0; j < channNum - 2; j++) mc_pipi[i]->Add(hist_pipi[i][j]);*/

   TFractionFitter* fit_semi;
   //TFractionFitter* fit_three[10];
   //TFractionFitter* fit_pipi[10];

   Int_t bins_data, bin_width_data;
   Double_t scaling_factor, factor_error, int_initial[channNum], int_data, int_initial_tot;

	fit_semi = new TFractionFitter(hist_semi[0][channNum-2],mc_semi);
	//fit_three[i] = new TFractionFitter(hist_three[i][channNum-2],mc_three[i]);
	//fit_pipi[i] = new TFractionFitter(hist_pipi[i][channNum-2],mc_pipi[i]);

	fit_semi->Constrain(0,0.0,1.0);
	fit_semi->Constrain(1,0.0,1.0);
	fit_semi->Constrain(2,0.0,1.0);
	fit_semi->Constrain(3,0.0,1.0);
	fit_semi->Constrain(4,0.0,1.0);
	fit_semi->Constrain(5,0.0,1.0);
	//fit_semi->Constrain(6,0.0,1.0);
		//fit_three[i]->Constrain(j,0.0,1.0);
		//fit_pipi[i]->Constrain(j,0.0,1.0);

	//fit_semi->Fit();
	//fit_three[i]->Fit();
	//fit_pipi[i]->Fit();

	//hist_semi[0][channNum - 1] = (TH1F*) fit_semi->GetPlot();
	//hist_three[i][channNum - 1] = (TH1F*) fit_three[i]->GetPlot();
	//hist_pipi[i][channNum - 1] = (TH1F*) fit_pipi[i]->GetPlot();
	
	/*bins_data = hist_semi[i][channNum-2]->GetNbinsX();
	bin_width_data = hist_semi[i][channNum-2]->GetBinWidth(1);
	int_initial[0] = hist_semi[i][0]->Integral(1,bins_data);
	int_initial[1] = hist_semi[i][1]->Integral(1,bins_data);
	int_initial[2] = hist_semi[i][2]->Integral(1,bins_data);
	int_initial[3] = hist_semi[i][3]->Integral(1,bins_data);
	int_initial[4] = hist_semi[i][4]->Integral(1,bins_data);
	int_initial[5] = hist_semi[i][5]->Integral(1,bins_data);
	int_initial[6] = hist_semi[i][6]->Integral(1,bins_data);
	int_data = hist_semi[i][channNum-2]->Integral(1,bins_data);*/

	//int_initial_tot = int_initial[0] + int_initial[1] + int_initial[2] + int_initial[3] + int_initial[4] + int_initial[5] + int_initial[6];

   	/*for(Int_t j = 0; j < channNum - 2; j++)
   	{
		fit_semi[i]->GetResult(j,scaling_factor,factor_error);

		cout << int_data*scaling_factor/int_initial_tot << endl;
		cout << scaling_factor << endl << endl;

		hist_semi[i][j]->Scale(scaling_factor*int_data/int_initial[j]);

		//fit_three[i]->Constrain(j,0.0,1.0);
		//fit_pipi[i]->Constrain(j,0.0,1.0);
		//

		hist_semi[i][channNum-1]->Add(hist_semi[i][j]);
	}*/



   for(Int_t i = 0; i < 10; i++)
   {
      legend[i] = new TLegend(0.65,0.7,0.9,0.9);
      //legend[i] = new TLegend(0.15,0.7,0.4,0.9);

      for(Int_t j = 0; j < channNum-2; j++)
      {
         if(j == channNum - 2) legend[i]->AddEntry(hist_semi[i][j], channName[j], "PE1");
         else legend[i]->AddEntry(hist_semi[i][j], channName[j], "L");
      }
   }
   
   gStyle->SetOptStat(0);
   gStyle->SetStatX(0.9);
   gStyle->SetStatY(0.9);
   gStyle->SetStatH(0.15);
   gStyle->SetStatW(0.25);

   Int_t index_sort_semi[10][channNum], index_sort_three[10][channNum], index_sort_pipi[10][channNum];
   Double_t max_counts_semi[10][channNum], max_counts_three[10][channNum], max_counts_pipi[10][channNum];
   TString title_semi[100], title_three[100], title_pipi[100];

   title_semi[0] = "#Deltat [#tau_{S}]";
   title_semi[1] = "t_{K#rightarrow#pi^{+}#pi^{-}} [#tau_{S}]";
   title_semi[2] = "t_{K#rightarrow#pi^{0}#pi^{0}} [#tau_{S}]";
   title_semi[3] = "t_{sum,r} [ns]";
   title_semi[4] = "#alpha_{#pi^{+},#pi^{-}} [#circ]";
   title_semi[5] = "Q_{miss} [MeV]";
   title_semi[6] = "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]";
   title_semi[7] = "m^{inv}_{4#gamma} [MeV/c^{2}]";

   title_three[0] = "#Deltat [#tau_{S}]";
   title_three[1] = "t_{K#rightarrow#pi^{+}#pi^{-}} [#tau_{S}]";
   title_three[2] = "t_{K#rightarrow#pi^{0}#pi^{0}} [#tau_{S}]";
   title_three[3] = "z_{neu} [cm]";
   title_three[4] = "R [cm]";
   title_three[5] = "cos(#alpha^{CM}_{#pi^{+},#pi^{-}})";
   title_three[6] = "m^{inv}_{#pi^{+}#pi^{-}} [MeV/c^{2}]";
   title_three[7] = "m^{inv}_{4#gamma} [MeV/c^{2}]";

   title_pipi[0] = "#Deltat [#tau_{S}]";
   title_pipi[1] = "t_{K_{S}#rightarrow#pi^{+}#pi^{-}} [#tau_{S}]";
   title_pipi[2] = "t_{K_{L}#rightarrow#pi^{+}#pi^{-}} [#tau_{S}]";
   title_pipi[3] = "z_{neu} [cm]";
   title_pipi[4] = "z_{ch} [cm]";
   title_pipi[5] = "#sum_{i=1}^{4}t_{i,r} [ns]";
   title_pipi[6] = "m^{inv}_{#pi^{+}#pi^{-},1} [MeV/c^{2}]";
   title_pipi[7] = "m^{inv}_{#pi^{+}#pi^{-},2} [MeV/c^{2}]";

   Float_t y1d_semi[10][2], y1d_three[10][2], y1d_pipi[10][2], x2d[2], y2d[2]; 

   for(Int_t i = 0; i < 8; i++)
      for(Int_t j = 0; j < channNum; j++) max_counts_semi[i][j] = hist_semi[i][j]->GetMaximum();

   for(Int_t i = 0; i < 8; i++)
      for(Int_t j = 0; j < channNum; j++) max_counts_three[i][j] = hist_three[i][j]->GetMaximum();
   
   for(Int_t i = 0; i < 8; i++)
      for(Int_t j = 0; j < channNum; j++) max_counts_pipi[i][j] = hist_pipi[i][j]->GetMaximum();

   for(Int_t i = 0; i < 8; i++)
   {
      TMath::Sort(channNum-2, max_counts_semi[i], index_sort_semi[i]);
      y1d_semi[i][0] = 0.1;
      y1d_semi[i][1] = 1E10;//max_counts_semi[i][index_sort_semi[i][0]] + 200.;
   }

   for(Int_t i = 0; i < 8; i++)
   {
      TMath::Sort(channNum-2, max_counts_three[i], index_sort_three[i]);
      y1d_three[i][0] = 0.1;
      y1d_three[i][1] = 1E10;//max_counts_three[i][index_sort_three[i][0]] + 200.;
   }

   for(Int_t i = 0; i < 8; i++)
   {
      TMath::Sort(channNum-2, max_counts_pipi[i], index_sort_pipi[i]);
      y1d_pipi[i][0] = 0.1;
      y1d_pipi[i][1] = 1E10;//max_counts_pipi[i][index_sort_pipi[i][0]] + 200.;
   }

   Float_t cut_value = 2; 

   std::cout << hist_semi[4][4]->GetBinCenter(hist_semi[4][4]->GetMaximumBin()) << std::endl;

   TLine *line_down = new TLine(mK0-2,0,mK0-2,1E10);
   TLine *line_up = new TLine(mK0+2,0,mK0+2,1E10);
   line_up->SetLineWidth(3);
   line_down->SetLineWidth(3);
   line_up->SetLineStyle(10);
   line_down->SetLineStyle(10);

   TBox *box_down = new TBox(mK0-20,0,mK0-2,1E10);
   TBox *box_up = new TBox(mK0+2,0,mK0+20,1E10);

   box_down->SetFillStyle(3544);
   box_up->SetFillStyle(3544);

   box_down->SetFillColor(kBlack);
   box_up->SetFillColor(kBlack);

   for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < channNum-2; j++)
      {
         canva1d_semi[i]->SetLeftMargin(0.15);
         canva1d_semi[i]->SetBottomMargin(0.15);
	      canva1d_semi[i]->SetLogy(1);
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
	       if(index_sort_semi[i][j] == channNum - 2) hist_semi[i][index_sort_semi[i][j]]->Draw("PE1");
	       else hist_semi[i][index_sort_semi[i][j]]->Draw("HIST");
            }
            else
            {
	       if(index_sort_semi[i][j] == channNum - 2) hist_semi[i][index_sort_semi[i][j]]->Draw("PE1SAME");
	       else hist_semi[i][index_sort_semi[i][j]]->Draw("HISTSAME");
            }
         
      }

	//box_up->Draw();
	//box_down->Draw();

	//line_up->Draw();
	//line_down->Draw();
      legend[i]->Draw();
   }

   for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < channNum-2; j++)
      {
         canva1d_three[i]->SetLeftMargin(0.15);
         canva1d_three[i]->SetBottomMargin(0.15);
	 canva1d_three[i]->SetLogy(1);
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
	       if(index_sort_three[i][j] == channNum - 2) hist_three[i][index_sort_three[i][j]]->Draw("PE1");
	       else hist_three[i][index_sort_three[i][j]]->Draw("HIST");
            }
            else
            {
	       if(index_sort_three[i][j] == channNum - 2) hist_three[i][index_sort_three[i][j]]->Draw("PE1SAME");
	       else hist_three[i][index_sort_three[i][j]]->Draw("HISTSAME");
            }
         
      }
      //box_down->Draw();
      //line_down->Draw();
      legend[i]->Draw();
   }

   for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < channNum-2; j++)
      {
         canva1d_pipi[i]->SetLeftMargin(0.15);
         canva1d_pipi[i]->SetBottomMargin(0.15);
	 canva1d_pipi[i]->SetLogy(1);
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

	       if(index_sort_pipi[i][j] == channNum - 2) hist_pipi[i][index_sort_pipi[i][j]]->Draw("PE1");
	       else hist_pipi[i][index_sort_pipi[i][j]]->Draw("HIST");
            }
            else
            {
	       if(index_sort_pipi[i][j] == channNum - 2) hist_pipi[i][index_sort_pipi[i][j]]->Draw("PE1SAME");
	       else hist_pipi[i][index_sort_pipi[i][j]]->Draw("HISTSAME");
            }
         
      }

	//box_up->Draw();
	//box_down->Draw();

	//line_up->Draw();
	//line_down->Draw();
      legend[i]->Draw();
   }

   /*for(Int_t i = 0; i < 8; i++)
   {
      for(Int_t j = 0; j < channNum-2; j++)
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

   TLegend *legend_nocuts[3];

   for(Int_t i = 0; i < 3; i++)
   {
      legend_nocuts[i] = new TLegend(0.45,0.8,0.9,0.9);
      //legend_nocuts[i] = new TLegend(0.15,0.7,0.4,0.9);
   }

   legend_nocuts[0]->AddEntry(hist_semi[0][channNum-1], "Expected #Deltat distribution for K#rightarrow#pi^{0}#pi^{0} from MC", "L");
   legend_nocuts[0]->AddEntry(hist_semi[0][channNum-2], "Expected #Deltat distribution for K#rightarrow#pi^{0}#pi^{0} from DATA", "PE1");
   legend_nocuts[1]->AddEntry(hist_semi_nocuts, "Expected #Deltat distribution for K_{S}#rightarrow#pi^{+}#pi^{-}", "L");
   legend_nocuts[2]->AddEntry(hist_semi_nocuts, "Expected #Deltat distribution for K#rightarrow#pi^{+}#pi^{-}", "L");

   canva1d_seminocuts->SetLeftMargin(0.15);
   canva1d_seminocuts->SetBottomMargin(0.15);
   canva1d_seminocuts->cd();

   hist_semi[0][channNum-1]->GetYaxis()->SetRangeUser(0,1.3*hist_semi[0][channNum-1]->GetMaximum());
   hist_semi[0][channNum-1]->SetLineColor(kGray+1);
   hist_semi[0][channNum-1]->GetYaxis()->SetTitle("Counts");
   hist_semi[0][channNum-1]->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");

   hist_semi[0][channNum-2]->SetMarkerStyle(8);
   hist_semi[0][channNum-2]->SetMarkerColor(kBlack);
   hist_semi[0][channNum-2]->SetLineColor(kBlack);

   hist_semi[0][channNum-1]->Draw("HIST");
   hist_semi[0][channNum-2]->Draw("PE1SAME");
   legend_nocuts[0]->Draw();

   /*canva1d_threenocuts->SetLeftMargin(0.15);
   canva1d_threenocuts->SetBottomMargin(0.15);
   canva1d_threenocuts->cd();

   hist_three_nocuts->Draw("HIST");

   canva1d_pipinocuts->SetLeftMargin(0.15);
   canva1d_pipinocuts->SetBottomMargin(0.15);
   canva1d_pipinocuts->cd();

   hist_pipi_nocuts->Draw("HIST");*/

   //hist_semi_nocuts->Add(hist_three_nocuts);
   //hist_semi_nocuts->Add(hist_pipi_nocuts);

   //hist_semi_withcuts->Add(hist_three_withcuts);
   //hist_semi_withcuts->Add(hist_pipi_withcuts);
 
   //hist_semi_withcuts->Divide(hist_semi_nocuts);

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

   canva1d_seminocuts->Print("plots/Semi/deltat_semi_fit_exp.png");
   canva1d_threenocuts->Print("plots/Three/deltat_three_nocuts.png");
   canva1d_pipinocuts->Print("plots/Pipi/deltat_pipi_nocuts.png");
   canva_efficiency->Print("efficiency_mc_cs.png");

   canva_eff_signal->Print("efficiency_mc.png");

   std::cout << "Efficiency three: " << 100*entries_sel[2][5]/(Float_t)entries[5] << "%" << std::endl; 
   std::cout << "Purity three: " << 100*entries_sel[2][5]/(Float_t)(entries_sel[2][4] + entries_sel[2][0] + entries_sel[2][1] + entries_sel[2][2] + entries_sel[2][3] + entries_sel[2][5] + entries_sel[2][6])<< "%" << std::endl << std::endl;

   std::cout << "Signal: " << entries[0] << std::endl;
   std::cout << "Regen: " << entries[1] << std::endl;
   std::cout << "Omega: " << entries[2] << std::endl;
   std::cout << "Three: " << entries[3] << std::endl;
   std::cout << "Semi: " << entries[4] << std::endl;
   std::cout << "Pipi: " << entries[5] << std::endl;
   std::cout << "Else: " << entries[6] << std::endl;
   std::cout << "DATA: " << entries[7] << std::endl;

   TLegend *legend_signal;

   legend_signal = new TLegend(0.15,0.75,0.65,0.9);

   legend_signal->AddEntry(pEff_signal, "Efficiency for time of flight method", "PE1");
   legend_signal->AddEntry(pEff_signal_tri, "Efficiency for trilateration method", "PE1");

   TCanvas *canvaefficiency = new TCanvas("canvaefficiency","", 750, 750);

   pEff_signal->SetStatisticOption(TEfficiency::kFNormal);

   pEff_signal->SetLineWidth(3);
   pEff_signal->SetLineColor(kRed);
   pEff_signal->Draw("APE1");

   legend_signal->Draw();

   pEff_signal_tri->SetStatisticOption(TEfficiency::kFNormal);

   pEff_signal_tri->SetLineWidth(3);
   pEff_signal_tri->SetLineColor(kBlue);
   pEff_signal_tri->Draw("PE1SAME");

   gPad->Update();
   pEff_signal->GetPaintedGraph()->GetXaxis()->SetLimits(-100.0, 100.0);
   pEff_signal->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 0.3);
   pEff_signal->GetPaintedGraph()->GetXaxis()->CenterTitle(1);
   pEff_signal->GetPaintedGraph()->GetYaxis()->CenterTitle(1);
   gPad->Update();

   canvaefficiency->Print("efficiency_tri.png");

   TLegend *legend_length;

   legend_length = new TLegend(0.15,0.75,0.65,0.9);

   legend_length->AddEntry(pEff_three_length, "Efficiency for K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}", "PE1");
   legend_length->AddEntry(pEff_pipi_length, "Efficiency for K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{+}#pi^{-}", "PE1");

   TCanvas *canvaefficiencylength = new TCanvas("canvaefficiencylength","", 750, 750);

   pEff_three_length->SetStatisticOption(TEfficiency::kFNormal);

   pEff_three_length->SetLineWidth(3);
   pEff_three_length->SetLineColor(kCyan);
   pEff_three_length->Draw("APE1");

   pEff_pipi_length->SetStatisticOption(TEfficiency::kFNormal);

   pEff_pipi_length->SetLineWidth(3);
   pEff_pipi_length->SetLineColor(kYellow);
   pEff_pipi_length->Draw("PE1SAME");

   legend_length->Draw();

   gPad->Update();
   pEff_three_length->GetPaintedGraph()->GetXaxis()->SetLimits(0.0, 50.0);
   pEff_three_length->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 0.3);
   pEff_three_length->GetPaintedGraph()->GetXaxis()->CenterTitle(1);
   pEff_three_length->GetPaintedGraph()->GetYaxis()->CenterTitle(1);
   gPad->Update();

   canvaefficiencylength->Print("efficiency_length.png");

   TH1 *passed_three, *passed_pipi, *passed_semi, *total_three, *total_pipi, *total_semi;

   TH1 *passed_three_mc, *passed_pipi_mc, *passed_semi_mc, *total_three_mc, *total_pipi_mc, *total_semi_mc; 
   
   passed_pipi = pEff_pipi->GetCopyPassedHisto();
   passed_three = pEff_three->GetCopyPassedHisto();
   passed_semi = pEff_semi->GetCopyPassedHisto();
   total_pipi = pEff_pipi->GetCopyTotalHisto();
   total_three = pEff_three->GetCopyTotalHisto();
   total_semi = pEff_semi->GetCopyTotalHisto();

   passed_three->Add(passed_pipi);
   passed_three->Add(passed_semi);
   total_three->Add(total_pipi);
   total_three->Add(total_semi);

   passed_pipi_mc = pEff_pipimc->GetCopyPassedHisto();
   passed_three_mc = pEff_threemc->GetCopyPassedHisto();
   passed_semi_mc = pEff_semimc->GetCopyPassedHisto();
   total_pipi_mc = pEff_pipimc->GetCopyTotalHisto();
   total_three_mc = pEff_threemc->GetCopyTotalHisto();
   total_semi_mc = pEff_semimc->GetCopyTotalHisto();

   passed_three_mc->Add(passed_pipi_mc);
   passed_three_mc->Add(passed_semi_mc);
   total_three_mc->Add(total_pipi_mc);
   total_three_mc->Add(total_semi_mc);

   TEfficiency *total = new TEfficiency(*passed_three, *total_three);
   TEfficiency *total_mc = new TEfficiency(*passed_three_mc, *total_three_mc);

   total->SetStatisticOption(TEfficiency::kFNormal);
   total_mc->SetStatisticOption(TEfficiency::kFNormal);

   TLegend *legend_eff;

   legend_eff = new TLegend(0.15,0.75,0.65,0.9);

   legend_eff->AddEntry(total, "Efficiency of control samples from DATA", "PE1");
   legend_eff->AddEntry(total_mc, "Efficiency of control samples from MC", "PE1");

   TCanvas *canvaefficiencyresult = new TCanvas("canvaefficiencyresult","", 750, 750);

   total->SetLineWidth(3);
   total->SetLineColor(kBlack);
   total->Draw("APE1");

   total_mc->SetLineWidth(3);
   total_mc->SetLineColor(kGray+1);
   total_mc->Draw("PE1SAME");

   legend_eff->Draw();

   gPad->Update();
   total->GetPaintedGraph()->GetXaxis()->SetLimits(-100.0, 100.0);
   total->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 0.3);
   total->GetPaintedGraph()->GetXaxis()->CenterTitle(1);
   total->GetPaintedGraph()->GetYaxis()->CenterTitle(1);
   gPad->Update();

   canvaefficiencyresult->Print("efficiency_result.png");

   //Division of plots of efficiency

   Double_t deltat[201] = {0}, div[201] = {0}, err[2][201] = {0};
   Double_t average = 0., denom = 0., nomin = 0., average_err = 0.;

   for(Int_t i = 1; i <= 201; i++)
   {
      deltat[i] = -100. + i*200/201;
      div[i] = total->GetEfficiency(i)/total_mc->GetEfficiency(i);
      err[0][i] = 0;
      err[1][i] = sqrt(pow(total->GetEfficiencyErrorUp(i)/total_mc->GetEfficiency(i),2) + pow(total_mc->GetEfficiencyErrorUp(i)*total->GetEfficiency(i)/pow(total_mc->GetEfficiency(i),2),2));

      denom += 1/pow((Float_t)err[1][i],2);
      nomin += div[i]/pow((Float_t)err[1][i],2);

   }

   average = nomin/denom;
   average_err = sqrt(1/denom);

   TLatex text_average;
   text_average.SetTextSize(0.035);
   TLine *line = new TLine(-100.0, average, 100.0, average);

   TGraphErrors *division = new TGraphErrors(201, deltat, div, err[0], err[1]);

   TCanvas *canvadivisionresult = new TCanvas("canvadivisionresult","", 750, 750);

   division->SetLineWidth(3);
   division->SetTitle("");
   division->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
   division->GetYaxis()->SetTitle("Correction factor");
   division->GetXaxis()->CenterTitle(1);
   division->GetYaxis()->CenterTitle(1);
   division->SetLineWidth(3);
   division->SetLineColor(kBlack);
   division->Draw("APE1");

   division->GetXaxis()->SetLimits(-100.0, 100.0);
   division->GetYaxis()->SetRangeUser(0.0, 2);

   line->SetLineWidth(3);
   line->SetLineColor(kRed);
   line->SetLineStyle(9);
   line->Draw();
   text_average.DrawLatex(-90., 1.8, Form("Average: %.3g#pm%.1g", average, average_err));

   canvadivisionresult->Print("division_result.png");


   TCanvas *canvasignal_after = new TCanvas("canvasignal_after","", 750, 750);
   TCanvas *canvasignal_before = new TCanvas("canvasignal_before","", 750, 750);

   signal_after->Scale((Float_t)signal_after->GetEntries()/signal_after->Integral(0,signal_after->GetNbinsX()+1));
   signal_before->Scale((Float_t)signal_before->GetEntries()/signal_before->Integral(0,signal_before->GetNbinsX()+1));

   signal_after->GetXaxis()->SetMaxDigits(3);
   signal_after->GetYaxis()->SetMaxDigits(3);
   signal_after->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
   signal_after->GetYaxis()->SetTitle("Counts");
   signal_after->GetXaxis()->CenterTitle(1);
   signal_after->GetYaxis()->CenterTitle(1);
   signal_after->GetYaxis()->SetRangeUser(0,1.3*signal_after->GetMaximum());
   signal_after->SetLineColor(kRed);

   signal_before->GetXaxis()->SetMaxDigits(3);
   signal_before->GetYaxis()->SetMaxDigits(3);
   signal_before->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
   signal_before->GetYaxis()->SetTitle("Counts");
   signal_before->GetXaxis()->CenterTitle(1);
   signal_before->GetYaxis()->CenterTitle(1);
   signal_before->GetYaxis()->SetRangeUser(0,1.3*signal_before->GetMaximum());
   signal_before->SetLineColor(kRed);

   canvasignal_after->cd();
   signal_after->Draw("HIST");

   canvasignal_before->cd();
   signal_before->Draw("HIST");

   canvasignal_after->Print("deltat_signal_after.png");
   canvasignal_before->Print("deltat_signal_before.png");

}
