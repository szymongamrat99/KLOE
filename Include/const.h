#ifndef CONST_H
#define CONST_H

#include <json.hpp>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <boost/progress.hpp> // for loading bar display
// #include <omp.h>              // for multi-threading using OpenMP

#include <TString.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include <ErrorLogs.h>
#include <MainMenu.h>
#include <RealTimeIntegration.h>

using json = nlohmann::json;

// Get the env variable for properties

const std::string
    kloedataPath = getenv("KLOE_DBV26_DK0"),
    kloeMCPath = getenv("KLOE_DBV26_MK0"),
    workdirPath = getenv("WORKDIR"),
    chainDataFiles = kloedataPath + "/*.root",
    chainMCFiles = kloeMCPath + "/*.root",
    pdgConstFilePath = (std::string)getenv("PDGAPI") + "/pdg_const.json",
    propertiesPath = getenv("PROPERTIESKLOE"),
    propName = propertiesPath + "/properties.json",
    rootfilesName = propertiesPath + "/root-files.json",
    cutlimitsName = propertiesPath + "/cut-limits.json";

static std::ifstream propertyFile(propName.c_str());
static json properties = json::parse(propertyFile);

struct PDGids
{
    const TString
        Re = "/S013EPS",
        Im = "/S013EPI",
        K0mass = "/S011M",
        TauS = "/S012T",
        TauL = "/S013T",
        deltaM = "/S013D",
        modEps = "/S013EP",
        phiPM = "/S013F+-",
        phi00 = "/S013FOO";
};

static std::ifstream fconst(pdgConstFilePath);
static json constants = json::parse(fconst);

// Constants used in the analysis
// Basic quantities
const double cVel = 29.9792458;       // cm/ns
const double hBar = 6.582119569E-34;  // MeV*s
const double eleCh = 1.602176634E-19; // C

// Particles' masses
const double mPhi = 1019.461;                               // MeV/c^2
const double mK0 = (Double_t)constants["values"]["/S011M"]; // MeV/c^2
const double mPi0 = 134.9768;                               // MeV/c^2
const double mPiCh = 139.57039;                             // MeV/c^2
const double mMuon = 105.6583755;                           // MeV/c^2
const double mElec = 0.510998950;                           // MeV/c^2
const double mOmega = 782.66;                               // MeV/c^2

// Branching ratios
// Phi
const double br_phi_kskl = 0.339;
const double br_phi_omegapi0 = 4.7E-5;

// K-short
const double br_ks_pi0pi0 = 0.3069;
const double br_ks_pippim = 0.6920;
const double br_ks_pippimgamma = 1.79E-3;
const double br_ks_piele = 7.04E-4;
const double br_ks_pimu = 4.56E-4;

// K-long
const double br_kl_pi0pi0 = 8.64E-4;
const double br_kl_pippim = 1.967E-3;
const double br_kl_pippimpi0 = 0.1254;
const double br_kl_3pi0 = 0.1952;
const double br_kl_piele = 0.4055;
const double br_kl_pimu = 0.2704;

// Kaons' properties and CPV
const double tau_S_nonCPT = (Double_t)constants["values"]["/S012T"] * 1E9; // ns
const double tau_S_CPT = 0.8954E-1;                                        // ns
const double tau_L = (Double_t)constants["values"]["/S013T"] * 1E9;        // ns
const double delta_mass_nonCPT = (Double_t)constants["values"]["/S013D"];  // hbar s^-1
const double delta_mass_CPT = 0.5293E10;                                   // hbar s^-1
const double mod_epsilon = (Double_t)constants["values"]["/S013EP"];
const double Re = (Double_t)constants["values"]["/S013EPS"];
const double Im_nonCPT = (Double_t)constants["values"]["/S013EPI"] * (M_PI / 180.); // deg
const double Im_CPT = -0.002;                                                       // deg
const double phi_pm_nonCPT = (Double_t)constants["values"]["/S013F+-"];             // deg
const double phi_pm_CPT = 43.51;                                                    // deg
const double phi_00_nonCPT = (Double_t)constants["values"]["/S013FOO"];             // deg
const double phi_00_CPT = 43.52;                                                    // deg

const double TRF = 2.715; // ns - time of DAFNE bunch

// General
const int T0 = 2.715; // ns
const unsigned int MaxNumtrkv = 200;
const unsigned int MIN_CLU_ENE = 20;

const unsigned int channNum = 7;

const TString
    channName[channNum] = {"K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                           "Regeneration",
                           "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                           "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}",
                           "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}",
                           "Other bcg", "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{+}#pi^{-}"},
    channelInt[channNum] = {"1", "3", "4", "5", "6", "7", "8"};
const TString dataName = "DATA";
const TString mcSumName = "MC sum";

const Color_t channColor[channNum] = {kRed, kGreen, kViolet, kCyan, kBlue, kGreen - 1, kYellow};
const Color_t dataColor = kBlack;
const Color_t mcSumColor = kOrange;

const TString
    base_path = workdirPath + "/KLOE/",
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
    gen_vars_dir = base_path + "Subanalysis/GeneratedVars/",
    neutrec_dir = base_path + "Subanalysis/Neutrec/",
    cpfit_dir = base_path + "Subanalysis/CPFit/",
    covmatrix_dir = base_path + "Subanalysis/CovarianceMatrix/",
    initialanalysis_dir = base_path + "Subanalysis/InitialAnalysis/",
    omegarec_dir = base_path + "Subanalysis/OmegaRec/",
    efficiency_dir = base_path + "Subanalysis/EfficiencyAnalysis/",
    charged_dir = base_path + "Subanalysis/KchRec/",
    plots_dir = base_path + "Subanalysis/Plots/",
    root_files_dir = "root_files/",
    input_dir = "input/",
    logs_dir = "log/",
    result_dir = "results/",
    img_dir = "img/";

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
    BaseKinematics () : 
                       ipTriKinFit(3, 0.0),
                       KnetriKinFit(10, 0.0),
                       neuVtxTriKinFit(4, 0.0),
                       gammaMomTriKinFit1(8, 0.0),
                       gammaMomTriKinFit2(8, 0.0),
                       gammaMomTriKinFit3(8, 0.0),
                       gammaMomTriKinFit4(8, 0.0),
                       ipTriangle(3, 0.0),
                       KneTriangle(10, 0.0),
                       neuVtxTriangle(4, 0.0),
                       gammaMomTriangle1(8, 0.0),
                       gammaMomTriangle2(8, 0.0),
                       gammaMomTriangle3(8, 0.0),
                       gammaMomTriangle4(8, 0.0),
                       g4takenTriKinFit(4, 0),
                       trcfinal(4, 0.0),
                       CurvSmeared1(0.0),
                       PhivSmeared1(0.0),
                       CotvSmeared1(0.0),
                       CurvSmeared2(0.0),
                       PhivSmeared2(0.0),
                       CotvSmeared2(0.0),
                       KchrecFit(10, 0),
                       KchboostFit(10, 0),
                       ipFit(3, 0),
                       KnerecFit(10, 0),
                       KnereclorFit(10, 0)
    {};

    Float_t
        Kchboost[9],
        KchboostKSOld[9],
        KchboostKLOld[9],
        Knereclor[9],
        Knerec[9],
        Kchrec[9],
        KchmcOld[9],
        KnemcOld[9],
        ip[3],
        ipmcOld[3],
        phi_mom[4],
        Dtmc,
        Dtrec,
        Dtboostlor,
        TclOld[50],
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
        CurvOld[200],
        PhivOld[200],
        CotvOld[200],
        xvOld[200],
        yvOld[200],
        zvOld[200],
        Bx,
        By,
        Bz,
        Bsx,
        Bsy,
        Bsz,
        Bpx,
        Bpy,
        Bpz,
        Bpxerr,
        Bpyerr,
        Bpzerr,
        Broots,
        BrootsErr,
        kaonChTimeCM,
        kaonChTimeLAB,
        kaonNeTimeLAB,
        kaonNeTimeCM,
        kaonNeTimeLABMC,
        kaonNeTimeCMMC,
        kaonChTimeLABMC,
        kaonChTimeCMMC,
        Chi2TriKinFit,
        CurvSmeared1,
        PhivSmeared1,
        CotvSmeared1,
        CurvSmeared2,
        PhivSmeared2,
        CotvSmeared2,
        Chi2SignalKinFit,
        Curv1,
        Phiv1,
        Cotv1,
        Curv2,
        Phiv2,
        Cotv2;

    UChar_t
        mctruth,
        mcisr,
        g4taken[4],
        mcflag,
        pidmcOld[MaxNumtrkv],
        vtxmcOld[MaxNumtrkv],
        motherOld[MaxNumtrkv],
        Vtx[MaxNumtrkv],
        ncll[200],
        VtxCh[3][MaxNumtrkv],
        ivOld[MaxNumtrkv];

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
        errFlagClosest,
        nev,
        nrun,
        necls,
        eclfilfo,
        eclfilfoword,
        ntcl,
        bunchnum;

    std::vector<Float_t>
        Kchrecnew,
        KchrecKS,
        KchrecKL,
        KchrecClosest,
        trknew[2],
        trkKS[2],
        trkKL[2],
        trkClosest[2],
        KchrecKSTwoBody,
        KchrecKLTwoBody,
        trkKSTwoBody[2],
        trkKLTwoBody[2],
        Xcl,
        Ycl,
        Zcl,
        Tcl,
        TclCorr,
        Enecl,
        Curv,
        Phiv,
        Cotv,
        xv,
        yv,
        zv,
        xvmc,
        yvmc,
        zvmc,
        pxmc,
        pymc,
        pzmc,
        KchboostKS,
        KchboostKL,
        ipKS,
        ipKL,
        ipmc,
        Knemc,
        Kchmc,
        trkMC[2],
        trkKSmc[2],
        trkKLmc[2],
        Kchrecsmeared,
        Kchboostsmeared,
        trksmeared[2],
        Kchboostnew,
        ipnew,
        ipTriKinFit,
        KnetriKinFit,
        neuVtxTriKinFit,
        gammaMomTriKinFit1,
        gammaMomTriKinFit2,
        gammaMomTriKinFit3,
        gammaMomTriKinFit4,
        ipTriangle,
        KneTriangle,
        neuVtxTriangle,
        gammaMomTriangle1,
        gammaMomTriangle2,
        gammaMomTriangle3,
        gammaMomTriangle4,
        trcfinal,
        CurvMC,
        PhivMC,
        CotvMC,
        pullsTriKinFit,
        trkFit[2],
        KchrecFit,
        KchboostFit,
        ipFit,
        photonFit[4],
        KnerecFit,
        KnereclorFit,
        ParamSignal,
        ErrorsSignal,
        ParamSignalFit,
        ErrorsSignalFit,
        pullsSignalFit;

    TLorentzVector
        phi4Mom;

    std::vector<Int_t>
        vtaken,
        vtakenKS,
        vtakenKL,
        vtakenClosest,
        eclstream,
        Asscl,
        iv,
        pidmc,
        vtxmc,
        mother,
        goodClusIndex,
        errors,
        cuts,
        g4takenTriKinFit;

    /**
     * @brief Clear all data members of BaseKinematics (scalars to zero, vectors cleared/resized).
     */
    void clear() {
        // Zero all fixed-size arrays and scalars
        memset(Kchboost, 0, sizeof(Kchboost));
        memset(KchboostKSOld, 0, sizeof(KchboostKSOld));
        memset(KchboostKLOld, 0, sizeof(KchboostKLOld));
        memset(Knereclor, 0, sizeof(Knereclor));
        memset(Knerec, 0, sizeof(Knerec));
        memset(Kchrec, 0, sizeof(Kchrec));
        memset(KchmcOld, 0, sizeof(KchmcOld));
        memset(KnemcOld, 0, sizeof(KnemcOld));
        memset(ip, 0, sizeof(ip));
        memset(ipmcOld, 0, sizeof(ipmcOld));
        memset(phi_mom, 0, sizeof(phi_mom));
        Dtmc = Dtrec = Dtboostlor = T0step1 = Chi2 = minv4gam = Qmiss = 0.0f;
        memset(TclOld, 0, sizeof(TclOld));
        memset(cluster, 0, sizeof(cluster));
        memset(bhabha_vtx, 0, sizeof(bhabha_vtx));
        memset(Pgamrec, 0, sizeof(Pgamrec));
        memset(omega, 0, sizeof(omega));
        memset(trk, 0, sizeof(trk));
        memset(pi0, 0, sizeof(pi0));
        memset(CurvOld, 0, sizeof(CurvOld));
        memset(PhivOld, 0, sizeof(PhivOld));
        memset(CotvOld, 0, sizeof(CotvOld));
        memset(xvOld, 0, sizeof(xvOld));
        memset(yvOld, 0, sizeof(yvOld));
        memset(zvOld, 0, sizeof(zvOld));
        mctruth = mcisr = mcflag = 0;
        memset(g4taken, 0, sizeof(g4taken));
        memset(pidmcOld, 0, sizeof(pidmcOld));
        memset(vtxmcOld, 0, sizeof(vtxmcOld));
        memset(motherOld, 0, sizeof(motherOld));
        memset(Vtx, 0, sizeof(Vtx));
        memset(ncll, 0, sizeof(ncll));
        memset(VtxCh, 0, sizeof(VtxCh));
        memset(ivOld, 0, sizeof(ivOld));
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
