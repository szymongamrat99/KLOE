#ifndef CONST_H
#define CONST_H

#include <json.hpp>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <boost/progress.hpp> // for loading bar display
#include <omp.h>              // for multi-threading using OpenMP

#include <TString.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include <ErrorLogs.h>
#include <MainMenu.h>
#include <RealTimeIntegration.h>
#include <SystemPaths.h> // for system paths

using json = nlohmann::json;

// Get the env variable for properties
static std::ifstream propertyFile(SystemPath::generalPropertiesPath.c_str());
static json properties = json::parse(propertyFile);

// General
const int T0 = 2.715; // ns
const unsigned int MaxNumtrkv = 200;
const unsigned int MIN_CLU_ENE = 20;

const unsigned int channNum = 6;

const TString
    channName[channNum] = {"K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                           "Regeneration",
                           "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                           "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}",
                           "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}",
                           "Other bcg"},
    channelInt[channNum] = {"1", "3", "4", "5", "6", "7"};
const TString dataName = "DATA";
const TString mcSumName = "MC sum";

const Color_t channColor[channNum] = {kRed, kGreen, kViolet, kCyan, kBlue, kGreen - 1};
const Color_t dataColor = kBlack;
const Color_t mcSumColor = kOrange;

const TString
    path_tmp = (std::string)properties["variables"]["rootFiles"]["path"],
    path_cs = (std::string)properties["variables"]["rootFiles"]["pathControlSample"],
    prod2root_path_v26 = "/data/k2/DBV-26/DK0",
    ext_root = ".root",
    ext_img = ".png",
    ext_csv = ".csv";

const TString
    gen_vars_filename = "gen_vars_",
    mctruth_filename = "mctruth_",
    neu_tri_filename = "neuvtx_tri_rec_",
    neu_triangle_filename = "neuvtx_triangle_rec_",
    neu_trilateration_kin_fit_filename = "neuvtx_tri_kin_fit_",
    omega_rec_filename = "omega_rec_",
    omega_rec_kin_fit_filename = "omega_rec_kin_fit_",
    cut_vars_filename = "cut_vars_";

const TString
    omegaRecPath = (std::string)properties["variables"]["tree"]["filename"]["omegarec"],
    mctruthPath = (std::string)properties["variables"]["tree"]["filename"]["mctruth"],
    genvarsPath = (std::string)properties["variables"]["tree"]["filename"]["generatedvars"],
    trianglePath = (std::string)properties["variables"]["tree"]["filename"]["trianglefinal"];

const TString
    gen_vars_tree = (std::string)properties["variables"]["tree"]["treename"]["mctruth"],
    neutrec_triangle_tree = (std::string)properties["variables"]["tree"]["treename"]["trianglefinal"],
    neutrec_tri_tree = "h_tri",
    neutrec_kin_fit_tree = "h_tri_kin_fit",
    omegarec_tree = (std::string)properties["variables"]["tree"]["treename"]["omegarec"],
    omegarec_kin_fit_tree = (std::string)properties["variables"]["tree"]["treename"]["omegarec"],
    cut_vars_tree = "h_cut_vars";

const int
    firstFileMax = properties["variables"]["rootFiles"]["firstFileMax"],
    lastFileMax = properties["variables"]["rootFiles"]["lastFileMax"],
    numOfThreads = properties["variables"]["parallelization"]["numOfThreads"];

struct NeuPart
{
    TLorentzVector
        vtxLAB,
        MomLAB,
        vtxPhiCM,
        MomPhiCM,
        vtxMotherCM,
        MomMotherCM;

    Double_t
        InvMass,
        TotMomLAB,
        TotMomPhiCM,
        TotMomMotherCM;

    Int_t
        CluNum[2];
};

struct ChPart
{
    TLorentzVector
        vtxLAB,
        MomLAB,
        vtxPhiCM,
        MomPhiCM,
        vtxMotherCM,
        MomMotherCM;

    Double_t
        InvMass,
        TotMomLAB,
        TotMomPhiCM,
        TotMomMotherCM;

    Int_t
        TrkNum,
        VtxNum;
};

struct Phi
{
    TLorentzVector
        vtxLAB,
        momentumLAB,
        vtxPhiCM,
        momentumPhiCM;

    Double_t
        InvMass,
        TotMomentum;

    Int_t
        TrkNum,
        VtxNum;
};

struct BaseKinematics
{
    Float_t
        Kchboost[9],
        KchboostKS[9],
        KchboostKL[9],
        Knereclor[9],
        Knerec[9],
        Kchrec[9],
        Kchmc[9],
        Knemc[9],
        ip[3],
        ipmc[3],
        phi_mom[4],
        Dtmc,
        Dtrec,
        Dtboostlor,
        Tcl[50],
        cluster[5][200],
        bhabha_vtx[3],
        T0step1,
        Chi2,
        minv4gam,
        Qmiss,
        Pgamrec[4][4],
        omega[9],
        trk[2][4],
        pi0[2][6],
        Curv[200],
        Phiv[200],
        Cotv[200],
        xv[200],
        yv[200],
        zv[200];

    UChar_t
        mctruth,
        mcisr,
        g4taken[4],
        mcflag,
        pidmc[MaxNumtrkv],
        vtxmc[MaxNumtrkv],
        mother[MaxNumtrkv],
        Vtx[MaxNumtrkv],
        ncll[200],
        VtxCh[3][MaxNumtrkv],
        iv[MaxNumtrkv];

    Int_t
        nevent,
        ntmc,
        nvtxmc,
        nclu,
        nv,
        ntv,
        mctruth_int,
        errFlag,
        errFlagKS,
        errFlagKL,
        errFlagClosest;

    std::vector<Float_t>
        Kchrecnew,
        KchrecKS,
        KchrecKL,
        KchrecClosest,
        trknew[2],
        trkKS[2],
        trkKL[2],
        trkClosest[2],
        KchrecKLTwoBody,
        trkKLTwoBody[2];

    std::vector<Int_t>
        vtaken,
        vtakenKS,
        vtakenKL,
        vtakenClosest;

    /**
     * @brief Clear all data members of BaseKinematics (scalars to zero, vectors cleared/resized).
     */
    void clear() {
        // Zero all fixed-size arrays and scalars
        memset(Kchboost, 0, sizeof(Kchboost));
        memset(KchboostKS, 0, sizeof(KchboostKS));
        memset(KchboostKL, 0, sizeof(KchboostKL));
        memset(Knereclor, 0, sizeof(Knereclor));
        memset(Knerec, 0, sizeof(Knerec));
        memset(Kchrec, 0, sizeof(Kchrec));
        memset(Kchmc, 0, sizeof(Kchmc));
        memset(Knemc, 0, sizeof(Knemc));
        memset(ip, 0, sizeof(ip));
        memset(ipmc, 0, sizeof(ipmc));
        memset(phi_mom, 0, sizeof(phi_mom));
        Dtmc = Dtrec = Dtboostlor = T0step1 = Chi2 = minv4gam = Qmiss = 0.0f;
        memset(Tcl, 0, sizeof(Tcl));
        memset(cluster, 0, sizeof(cluster));
        memset(bhabha_vtx, 0, sizeof(bhabha_vtx));
        memset(Pgamrec, 0, sizeof(Pgamrec));
        memset(omega, 0, sizeof(omega));
        memset(trk, 0, sizeof(trk));
        memset(pi0, 0, sizeof(pi0));
        memset(Curv, 0, sizeof(Curv));
        memset(Phiv, 0, sizeof(Phiv));
        memset(Cotv, 0, sizeof(Cotv));
        memset(xv, 0, sizeof(xv));
        memset(yv, 0, sizeof(yv));
        memset(zv, 0, sizeof(zv));
        mctruth = mcisr = mcflag = 0;
        memset(g4taken, 0, sizeof(g4taken));
        memset(pidmc, 0, sizeof(pidmc));
        memset(vtxmc, 0, sizeof(vtxmc));
        memset(mother, 0, sizeof(mother));
        memset(Vtx, 0, sizeof(Vtx));
        memset(ncll, 0, sizeof(ncll));
        memset(VtxCh, 0, sizeof(VtxCh));
        memset(iv, 0, sizeof(iv));
        nevent = ntmc = nvtxmc = nclu = nv = ntv = mctruth_int = 0;
        errFlag = errFlagKS = errFlagKL = errFlagClosest = 0;
        // Clear all std::vectors
        Kchrecnew.clear();
        KchrecKS.clear();
        KchrecKL.clear();
        KchrecClosest.clear();
        KchrecKLTwoBody.clear();
        trknew[0].clear(); trknew[1].clear();
        trkKS[0].clear(); trkKS[1].clear();
        trkKL[0].clear(); trkKL[1].clear();
        trkClosest[0].clear(); trkClosest[1].clear();
        trkKLTwoBody[0].clear(); trkKLTwoBody[1].clear();
        vtaken.clear();
        vtakenKS.clear();
        vtakenKL.clear();
        vtakenClosest.clear();
    }
};

struct NeutRec4
{
    Float_t
        Knerec[10],
        Photons[4][9],
        chi2min;
    Int_t
        gtaken[4],
        done;
};

inline void setGlobalStyle()
{
    // Global Style of histograms, pads, etc.

    gStyle->SetOptStat("iouMn");

    gStyle->SetFitFormat("6.2g");
    gStyle->SetStatFormat("6.2g");

    gStyle->SetCanvasDefH(750);
    gStyle->SetCanvasDefW(750);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetOptLogz(1);
    gStyle->SetPalette(1);

    gStyle->SetHistLineWidth(3);
    gStyle->SetLineWidth(2);

    gStyle->SetLabelSize(0.04, "X");
    gStyle->SetLabelSize(0.04, "Y");
    gStyle->SetLabelSize(0.04, "Z");
    gStyle->SetLabelOffset(0.02, "X");
    gStyle->SetLabelOffset(0.02, "Y");
    gStyle->SetLabelOffset(0.02, "Z");
    gStyle->SetLabelFont(62, "X");
    gStyle->SetLabelFont(62, "Y");
    gStyle->SetLabelFont(62, "Z");

    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.05, "Y");
    gStyle->SetTitleSize(0.05, "Z");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.5, "Y");
    gStyle->SetTitleOffset(1.2, "Z");
    gStyle->SetTitleFont(62, "X");
    gStyle->SetTitleFont(62, "Y");
    gStyle->SetTitleFont(62, "Z");

    gStyle->SetTitle("");

    gStyle->cd();
}

inline TString elapsedTimeHMS(double totalSeconds)
{
    int elapsedMinutes, elapsedHours;
    double elapsedSeconds;
    TString elapsedHMS;

    elapsedHours = int(totalSeconds / 3600.);
    elapsedMinutes = int(((totalSeconds / 3600.) - elapsedHours) * 60.);
    elapsedSeconds = ((((totalSeconds / 3600.) - elapsedHours) * 60.) - elapsedMinutes) * 60.;

    elapsedHMS = std::to_string(elapsedHours) + "h " + std::to_string(elapsedMinutes) + "min " + std::to_string(elapsedSeconds) + "s";

    return elapsedHMS;
};

#endif
