#define correction_factor_cxx
// The class definition in correction_factor.h has been generated automatically
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
// root> T->Process("correction_factor.C")
// root> T->Process("correction_factor.C","some options")
// root> T->Process("correction_factor.C+")
//


#include "correction_factor.h"
#include <TH2.h>
#include <TStyle.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "../Include/const.h"
#include <kloe_class.h>
#include <interference.h>

auto *corr_file = new TFile("corr_file.root","RECREATE");
//auto *corr_tree = new TTree("corr","Correction factor histo");
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

const Double_t x_min = -90.0, x_max = 90.0;
const UInt_t nbins = 1 + (x_max - x_min)/2.;

using namespace KLOE;

interference event;

void correction_factor::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   pEff_signal = new TEfficiency("eff_signal",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);
   pEff_signal_tri = new TEfficiency("eff_signal_tri",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);
   pEff_semi = new TEfficiency("eff_semi",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);
   pEff_three = new TEfficiency("eff_three",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);
   pEff_pipi = new TEfficiency("eff_pipi",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);

   pEff_three_length = new TEfficiency("eff_three_length",";Length of path of Kaon [cm];Efficiency",100,0.,50.);
   pEff_pipi_length = new TEfficiency("eff_pipi_length",";Length of path of Kaon [cm];Efficiency",100,0.,50.);

   pEff_semimc = new TEfficiency("eff_semimc",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);
   pEff_threemc = new TEfficiency("eff_threemc",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);
   pEff_pipimc = new TEfficiency("eff_pipimc",";#Deltat [#tau_{S}];Efficiency", nbins, x_min, x_max);

   TString option = GetOption();
}

void correction_factor::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t correction_factor::Process(Long64_t entry)
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

   Double_t tot_br_semi = ((br_kl_piele + br_kl_pimu)/(br_kl_piele + br_kl_pimu + br_ks_piele + br_ks_pimu))*br_ks_pi0pi0 + ((br_ks_piele + br_ks_pimu)/(br_kl_piele + br_kl_pimu + br_ks_piele + br_ks_pimu))*br_kl_pi0pi0;
Double_t tot_br_pipi = ((br_kl_pippim)/(br_kl_pippim + br_ks_pippim))*br_ks_pippim + ((br_ks_pippim)/(br_kl_pippim + br_ks_pippim))*br_kl_pippim;

   fReader.SetLocalEntry(entry);

   //Momentum and position of Phi in LAB
   momPhi_mc(0) = KchmcOld[0] + KnemcOld[0];
   momPhi_mc(1) = KchmcOld[1] + KnemcOld[1];
   momPhi_mc(2) = KchmcOld[2] + KnemcOld[2];
   momPhi_mc(3) = KchmcOld[3] + KnemcOld[3];

   //Momentum and position of Phi in LAB
   momPhi(0) = *Bpx;
   momPhi(1) = *Bpy;
   momPhi(2) = *Bpz;
   momPhi(3) = *Broots;

   phi_vel = PhysicsConstants::cVel*(sqrt(pow(momPhi(0),2) + pow(momPhi(1),2) + pow(momPhi(2),2))/momPhi(3));

   pos_phi(0) = 0.;
   pos_phi(1) = 0.;
   pos_phi(2) = 0.;
   pos_phi(3) = PhysicsConstants::cVel*(sqrt(pow(pos_phi(0),2) + pow(pos_phi(1),2) + pow(pos_phi(2),2))/phi_vel);

   //Momentum and position of charged kaon MC
   mom_kch_mc(0) = KchmcOld[0];
   mom_kch_mc(1) = KchmcOld[1];
   mom_kch_mc(2) = KchmcOld[2];
   mom_kch_mc(3) = KchmcOld[3];

   kaon_pm_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kch_mc(0),2) + pow(mom_kch_mc(1),2) + pow(mom_kch_mc(2),2))/mom_kch_mc(3));

   pos_kch_mc(0) = KchmcOld[6] - ipmcOld[0];
   pos_kch_mc(1) = KchmcOld[7] - ipmcOld[1];
   pos_kch_mc(2) = KchmcOld[8] - ipmcOld[2];
   pos_kch_mc(3) = PhysicsConstants::cVel*(sqrt(pow(pos_kch_mc(0),2) + pow(pos_kch_mc(1),2) + pow(pos_kch_mc(2),2))/kaon_pm_vel);

   //Momentum and position of neutral kaon MC
   mom_kne_mc(0) = KnemcOld[0];
   mom_kne_mc(1) = KnemcOld[1];
   mom_kne_mc(2) = KnemcOld[2];
   mom_kne_mc(3) = KnemcOld[3];

   kaon_pm_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kne_mc(0),2) + pow(mom_kne_mc(1),2) + pow(mom_kne_mc(2),2))/mom_kne_mc(3));

   pos_kne_mc(0) = KnemcOld[6] - ipmcOld[0];
   pos_kne_mc(1) = KnemcOld[7] - ipmcOld[1];
   pos_kne_mc(2) = KnemcOld[8] - ipmcOld[2];
   pos_kne_mc(3) = PhysicsConstants::cVel*(sqrt(pow(pos_kne_mc(0),2) + pow(pos_kne_mc(1),2) + pow(pos_kne_mc(2),2))/kaon_pm_vel);

   //Momentum and position of charged kaon
   mom_kaon_pm(0) = Kchboost[0];
   mom_kaon_pm(1) = Kchboost[1];
   mom_kaon_pm(2) = Kchboost[2];
   mom_kaon_pm(3) = Kchboost[3];

   kaon_pm_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kaon_pm(0),2) + pow(mom_kaon_pm(1),2) + pow(mom_kaon_pm(2),2))/mom_kaon_pm(3));

   pos_kaon_pm(0) = Kchboost[6] - *Bx;
   pos_kaon_pm(1) = Kchboost[7] - *By;
   pos_kaon_pm(2) = Kchboost[8] - *Bz;
   pos_kaon_pm(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_pm(0),2) + pow(pos_kaon_pm(1),2) + pow(pos_kaon_pm(2),2))/kaon_pm_vel);

   //Momentum and position of charged kaon
   mom_kaon_pm_alt(0) = sqrt(pow(mom_kaon_pm(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_pm(0) - pos_phi(0))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(1) = sqrt(pow(mom_kaon_pm(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_pm(1) - pos_phi(1))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(2) = sqrt(pow(mom_kaon_pm(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_pm(2) - pos_phi(2))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(3) = mom_kaon_pm(3);

   if(sqrt(pow(Kchrec1[6] - *Bx,2) + pow(Kchrec1[7] - *By,2) + pow(Kchrec1[8] - *Bz,2)) < sqrt(pow(Kchrec2[6] - *Bx,2) + pow(Kchrec2[7] - *By,2) + pow(Kchrec2[8] - *Bz,2)))
   {
      //Momentum and position of charged kaon_s
      mom_kaon_s(0) = Kchrec1[0];
      mom_kaon_s(1) = Kchrec1[1];
      mom_kaon_s(2) = Kchrec1[2];
      mom_kaon_s(3) = Kchrec1[3];

      kaon_s_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kaon_s(0),2) + pow(mom_kaon_s(1),2) + pow(mom_kaon_s(2),2))/mom_kaon_s(3));

      pos_kaon_s(0) = Kchrec1[6] - *Bx;
      pos_kaon_s(1) = Kchrec1[7] - *By;
      pos_kaon_s(2) = Kchrec1[8] - *Bz;
      pos_kaon_s(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_s(0),2) + pow(pos_kaon_s(1),2) + pow(pos_kaon_s(2),2))/kaon_s_vel);

      //Momentum and position of charged kaon_s alternative
      mom_kaon_s_alt(0) = sqrt(pow(mom_kaon_s(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_s(0) - pos_phi(0))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));
      mom_kaon_s_alt(1) = sqrt(pow(mom_kaon_s(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_s(1) - pos_phi(1))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(2) = sqrt(pow(mom_kaon_s(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_s(2) - pos_phi(2))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(3) = mom_kaon_s(3);

      //Momentum and position of charged kaon_l
      mom_kaon_l(0) = Kchrec2[0];
      mom_kaon_l(1) = Kchrec2[1];
      mom_kaon_l(2) = Kchrec2[2];
      mom_kaon_l(3) = Kchrec2[3];

      kaon_l_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kaon_l(0),2) + pow(mom_kaon_l(1),2) + pow(mom_kaon_l(2),2))/mom_kaon_l(3));

      pos_kaon_l(0) = Kchrec2[6] - *Bx;
      pos_kaon_l(1) = Kchrec2[7] - *By;
      pos_kaon_l(2) = Kchrec2[8] - *Bz;
      pos_kaon_l(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_l(0),2) + pow(pos_kaon_l(1),2) + pow(pos_kaon_l(2),2))/kaon_l_vel);

      //Momentum and position of charged kaon_l alternative
      mom_kaon_l_alt(0) = sqrt(pow(mom_kaon_l(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_l(0) - pos_phi(0))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));
      mom_kaon_l_alt(1) = sqrt(pow(mom_kaon_l(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_l(1) - pos_phi(1))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(2) = sqrt(pow(mom_kaon_l(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_l(2) - pos_phi(2))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(3) = mom_kaon_l(3);
   }
   else
   {
      //Momentum and position of charged kaon_s
      mom_kaon_s(0) = Kchrec2[0];
      mom_kaon_s(1) = Kchrec2[1];
      mom_kaon_s(2) = Kchrec2[2];
      mom_kaon_s(3) = Kchrec2[3];

      kaon_s_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kaon_s(0),2) + pow(mom_kaon_s(1),2) + pow(mom_kaon_s(2),2))/mom_kaon_s(3));

      pos_kaon_s(0) = Kchrec2[6] - *Bx;
      pos_kaon_s(1) = Kchrec2[7] - *By;
      pos_kaon_s(2) = Kchrec2[8] - *Bz;
      pos_kaon_s(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_s(0),2) + pow(pos_kaon_s(1),2) + pow(pos_kaon_s(2),2))/kaon_s_vel);

      //Momentum and position of charged kaon_s alternative
      mom_kaon_s_alt(0) = sqrt(pow(mom_kaon_s(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_s(0) - pos_phi(0))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));
      mom_kaon_s_alt(1) = sqrt(pow(mom_kaon_s(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_s(1) - pos_phi(1))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(2) = sqrt(pow(mom_kaon_s(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_s(2) - pos_phi(2))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(3) = mom_kaon_s(3);

      //Momentum and position of charged kaon_l
      mom_kaon_l(0) = Kchrec1[0];
      mom_kaon_l(1) = Kchrec1[1];
      mom_kaon_l(2) = Kchrec1[2];
      mom_kaon_l(3) = Kchrec1[3];

      kaon_l_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kaon_l(0),2) + pow(mom_kaon_l(1),2) + pow(mom_kaon_l(2),2))/mom_kaon_l(3));

      pos_kaon_l(0) = Kchrec1[6] - *Bx;
      pos_kaon_l(1) = Kchrec1[7] - *By;
      pos_kaon_l(2) = Kchrec1[8] - *Bz;
      pos_kaon_l(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_l(0),2) + pow(pos_kaon_l(1),2) + pow(pos_kaon_l(2),2))/kaon_l_vel);

      //Momentum and position of charged kaon_l alternative
      mom_kaon_l_alt(0) = sqrt(pow(mom_kaon_l(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_l(0) - pos_phi(0))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));
      mom_kaon_l_alt(1) = sqrt(pow(mom_kaon_l(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_l(1) - pos_phi(1))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(2) = sqrt(pow(mom_kaon_l(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_l(2) - pos_phi(2))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(3) = mom_kaon_l(3);
   }

   //Momentum and position of neutral standard kaon
   mom_kaon_00_std(0) = Knereclor[0];
   mom_kaon_00_std(1) = Knereclor[1];
   mom_kaon_00_std(2) = Knereclor[2];
   mom_kaon_00_std(3) = Knereclor[3];

   kaon_00_std_vel = PhysicsConstants::cVel*(sqrt(pow(mom_kaon_00_std(0),2) + pow(mom_kaon_00_std(1),2) + pow(mom_kaon_00_std(2),2))/mom_kaon_00_std(3));

   pos_kaon_00_std(0) = Knereclor[6] - *Bx;
   pos_kaon_00_std(1) = Knereclor[7] - *By;
   pos_kaon_00_std(2) = Knereclor[8] - *Bz;
   pos_kaon_00_std(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_00_std(0),2) + pow(pos_kaon_00_std(1),2) + pow(pos_kaon_00_std(2),2))/kaon_00_std_vel);

   //Momentum and position of neutral standard kaon
   mom_kaon_00_std_alt(0) = sqrt(pow(mom_kaon_00_std(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_00_std(0) - pos_phi(0))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));
   mom_kaon_00_std_alt(1) = sqrt(pow(mom_kaon_00_std(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_00_std(1) - pos_phi(1))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));;
   mom_kaon_00_std_alt(2) = sqrt(pow(mom_kaon_00_std(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_00_std(2) - pos_phi(2))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));;
   mom_kaon_00_std_alt(3) = mom_kaon_00_std(3);

   //Momentum and position of neutral trilateration kaon
   mom_kaon_00_tri(0) = fourKnetri[0];
   mom_kaon_00_tri(1) = fourKnetri[1];
   mom_kaon_00_tri(2) = fourKnetri[2];
   mom_kaon_00_tri(3) = fourKnetri[3];

   kaon_00_tri_vel = PhysicsConstants::cVel*mom_kaon_00_tri.BoostVector().Mag();

   pos_kaon_00_tri(0) = fourKnetri[6] - *Bx;
   pos_kaon_00_tri(1) = fourKnetri[7] - *By;
   pos_kaon_00_tri(2) = fourKnetri[8] - *Bz;
   pos_kaon_00_tri(3) = PhysicsConstants::cVel*(pos_phi(3) + sqrt(pow(pos_kaon_00_tri(0),2) + pow(pos_kaon_00_tri(1),2) + pow(pos_kaon_00_tri(2),2))/kaon_00_tri_vel);

   //Momentum and position of neutral standard kaon
   mom_kaon_00_tri_alt(0) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_00_tri(0) - pos_phi(0))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));
   mom_kaon_00_tri_alt(1) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_00_tri(1) - pos_phi(1))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));;
   mom_kaon_00_tri_alt(2) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(PhysicsConstants::mK0,2))*(pos_kaon_00_tri(2) - pos_phi(2))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));;
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

   tpm = pos_kaon_pm(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);
   t00_std = pos_kaon_00_std(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);
   t00_tri = pos_kaon_00_tri(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);
   ts = pos_kaon_s(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);
   tl = pos_kaon_l(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);
   tpm_mc = pos_kch_mc(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);
   t00_mc = pos_kne_mc(3)/(PhysicsConstants::tau_S_nonCPT*PhysicsConstants::cVel);  

   DeltaT_signal = (tpm - t00_std);
   DeltaT_control = (tpm - t00_tri);
   DeltaT_pipi = (tl - ts);
   DeltaT_mc = (tpm_mc - t00_mc);

   for(Int_t i = 0; i < 4; i++) TRCV[i] = TclOld[fourg4taken[i]] - (sqrt(pow(Xcl[fourg4taken[i]] - fourKnetri[6],2) + pow(Ycl[fourg4taken[i]] - fourKnetri[7],2) + pow(Zcl[fourg4taken[i]] - fourKnetri[8],2))/PhysicsConstants::cVel) - (k_path00_tri/(k_beta00_tri*PhysicsConstants::cVel));
   trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

   if(*mctruth == 1)
	{
      trcv_sum_signal = trcv[g4taken[0]-1] + trcv[g4taken[1]-1] + trcv[g4taken[2]-1] + trcv[g4taken[3]-1];
      pEff_signal->FillWeighted(trcv_sum_signal > -1 && abs(*minv4gam - PhysicsConstants::mK0) < 76 && abs(Kchrec[5] - PhysicsConstants::mK0) < 1.2 && *Qmiss_inv < 3.75 && cos(M_PI*(*anglepipi_CM_kch/180.)) < -0.8, event.interf_function(DeltaT_mc,1,0), DeltaT_signal);
      pEff_signal_tri->FillWeighted(trcv_sum > -1 && abs(fourKnetri[5] - PhysicsConstants::mK0) < 76 && abs(Kchrec[5] - PhysicsConstants::mK0) < 1.2 && *Qmiss_inv < 3.75 && cos(M_PI*(*anglepipi_CM_kch/180.)) < -0.8, event.interf_function(DeltaT_mc,1,0), DeltaT_signal);

   }

   cuts_semi[0] = abs(*Qmiss_inv - 71.13) < 25;
   cuts_semi[1] = abs(*anglepipi_CM_kch - 145.8) < 10;

   cuts_signal_neutral_cs[0] = (abs(fourKnetri[5] - PhysicsConstants::mK0) < 76);
   cuts_signal_neutral_cs[1] = (trcv_sum > -1.0);

   cuts_signal_charged_cs[0] = (abs(Kchrec[5] - PhysicsConstants::mK0) < 1.2);
   cuts_signal_charged_cs[1] = (*Qmiss_inv < 3.75);
   cuts_signal_charged_cs[2] = (cos(M_PI*(*anglepipi_CM_kch)/180.) < -0.8);

   tot_cuts_semi = cuts_semi[0] && cuts_semi[1];
   tot_cuts_neutral = cuts_signal_neutral_cs[0] && cuts_signal_neutral_cs[1];
   tot_cuts_charged = cuts_signal_charged_cs[0] && cuts_signal_charged_cs[1] && cuts_signal_charged_cs[2];

   if(tot_cuts_semi && *mcflag == 1 && *mctruth != 0 && *done4 == 1) pEff_semimc->FillWeighted(tot_cuts_neutral, tot_br_semi, DeltaT_control);
   if(*done == 1 && *mcflag == 1 && *mctruth != 0 && *done4 == 1) pEff_threemc->FillWeighted(tot_cuts_charged, br_ks_pippim, DeltaT_control);
   if(donepipi[0] == 1 && donepipi[1] == 1 && abs(Kchrec1[5] - PhysicsConstants::mK0) < 2 && *mcflag == 1 && *mctruth != 0) pEff_pipimc->FillWeighted(tot_cuts_charged, tot_br_pipi, DeltaT_pipi);

   if(tot_cuts_semi && *mcflag == 0 && *done4 == 1) pEff_semi->FillWeighted(tot_cuts_neutral, tot_br_semi, DeltaT_control);
   if(*done == 1 && *mcflag == 0 && *done4 == 1) pEff_three->FillWeighted(tot_cuts_charged, br_ks_pippim, DeltaT_control);
   if(donepipi[0] == 1 && donepipi[1] == 1 && abs(Kchrec1[5] - PhysicsConstants::mK0) < 2 && *mcflag == 0) pEff_pipi->FillWeighted(tot_cuts_charged, tot_br_pipi, DeltaT_pipi);

   if(*done == 1 && *mcflag == 0 && *done4 == 1) pEff_three_length->FillWeighted(tot_cuts_charged, br_ks_pippim, sqrt(pow(Kchrec[6] - *Bx,2) + pow(Kchrec[7] - *By,2) + pow(Kchrec[8] - *Bz,2)));
   if(donepipi[0] == 1 && donepipi[1] == 1 && abs(Kchrec1[5] - PhysicsConstants::mK0) < 2 && *mcflag == 0) pEff_pipi_length->FillWeighted(tot_cuts_charged, tot_br_pipi, sqrt(pow(Kchrec2[6] - *Bx,2) + pow(Kchrec2[7] - *By,2) + pow(Kchrec2[8] - *Bz,2)));

   return kTRUE;
}

void correction_factor::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void correction_factor::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   pEff_signal->SetStatisticOption(TEfficiency::kFNormal);

   pEff_signal->SetLineWidth(3);
   pEff_signal->SetLineColor(kRed);
   pEff_signal->Draw("APE1");

   pEff_signal_tri->SetStatisticOption(TEfficiency::kFNormal);

   pEff_signal_tri->SetLineWidth(3);
   pEff_signal_tri->SetLineColor(kBlue);
   pEff_signal_tri->Draw("PE1SAME");

   gPad->Update();
   pEff_signal->GetPaintedGraph()->GetXaxis()->SetLimits(x_min, x_max);
   pEff_signal->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 0.3);
   pEff_signal->GetPaintedGraph()->GetXaxis()->CenterTitle(1);
   pEff_signal->GetPaintedGraph()->GetYaxis()->CenterTitle(1);
   gPad->Update();

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

   //corr_tree->Branch("total", "TEfficiency", &total);

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
   total->GetPaintedGraph()->GetXaxis()->SetLimits(x_min, x_max);
   total->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 0.3);
   total->GetPaintedGraph()->GetXaxis()->CenterTitle(1);
   total->GetPaintedGraph()->GetYaxis()->CenterTitle(1);
   gPad->Update();

   Double_t deltat[nbins], div[nbins], err[2][nbins];
   Double_t average = 0., denom = 0., nomin = 0., average_err = 0.;

   for(Int_t i = 0; i < nbins; i++)
   {
      deltat[i] = x_min + 2*i;
      div[i] = total->GetEfficiency(i)/total_mc->GetEfficiency(i);
      err[0][i] = 0;
      err[1][i] = sqrt(pow(total->GetEfficiencyErrorUp(i)/total_mc->GetEfficiency(i),2) + pow(total_mc->GetEfficiencyErrorUp(i)*total->GetEfficiency(i)/pow(total_mc->GetEfficiency(i),2),2));

      denom += 1/pow((Float_t)err[1][i],2);
      nomin += div[i]/pow((Float_t)err[1][i],2);

   }

   average = nomin/denom;
   average_err = sqrt(1/denom);

   TGraphErrors *division = new TGraphErrors(nbins, deltat, div, err[0], err[1]);

   division->SetLineWidth(3);
   division->SetTitle("");
   division->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
   division->GetYaxis()->SetTitle("Correction factor");
   division->GetXaxis()->CenterTitle(1);
   division->GetYaxis()->CenterTitle(1);
   division->SetLineWidth(3);
   division->SetLineColor(kBlack);
   division->Draw("APE1");

   division->GetXaxis()->SetLimits(x_min, x_max);
   division->GetYaxis()->SetRangeUser(0.0, 2);

   total->Write("cs_eff_total_data");
   total_mc->Write("cs_eff_total_mc");
   pEff_signal->Write("eff_signal");
   pEff_signal_tri->Write("eff_signal_tri");
   division->Write("correction_factor");

   corr_file->Write("");
   corr_file->Close();
}