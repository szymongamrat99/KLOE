#pragma once

#include <json.hpp>

#include <TString.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>

#include <TMath.h>

// Get the env variable for properties

namespace Paths
{
  extern const std::string kloedataPath;
  extern const std::string kloeMCPath;
  extern const std::string workdirPath;
  extern const std::string chainDataFiles;
  extern const std::string chainMCFiles;
  extern const std::string pdgConstFilePath;
  extern const std::string propertiesPath;
  extern const std::string propName;
  extern const std::string analysisConfigPath;
  extern const std::string rootfilesName;
  extern const std::string cutlimitsName;

  extern TString base_path;
  extern TString path_tmp;
  extern TString path_cs;
  extern TString prod2root_path_v26;
  extern TString ext_root;
  extern TString ext_img;
  extern TString ext_csv;

  extern const TString gen_vars_dir;
  extern const TString neutrec_dir;
  extern const TString cpfit_dir;
  extern const TString covmatrix_dir;
  extern const TString initialanalysis_dir;
  extern const TString omegarec_dir;
  extern const TString efficiency_dir;
  extern const TString charged_dir;
  extern const TString plots_dir;
  extern const TString root_files_dir;
  extern const TString input_dir;
  extern const TString logs_dir;
  extern const TString result_dir;
  extern const TString img_dir;
}

namespace Filenames
{
  extern const TString gen_vars_filename;
  extern const TString mctruth_filename;
  extern const TString neu_tri_filename;
  extern const TString neu_triangle_filename;
  extern const TString neu_trilateration_kin_fit_filename;
  extern const TString omega_rec_filename;
  extern const TString omega_rec_kin_fit_filename;
  extern const TString cut_vars_filename;

  extern TString gen_vars_tree;
  extern TString neutrec_triangle_tree;
  extern TString neutrec_tri_tree;
  extern TString neutrec_kin_fit_tree;
  extern TString omegarec_tree;
  extern TString omegarec_kin_fit_tree;
  extern TString cut_vars_tree;

  extern TString omegaRecPath;
  extern TString mctruthPath;
  extern TString genvarsPath;
  extern TString trianglePath;
}

namespace PhysicsConstants
{
  // Constants used in the analysis
  // Fundamental quantities
  extern Double_t cVel;  // cm/ns
  extern Double_t hBar;  // MeV*s
  extern Double_t eleCh; // C

  // Particles' masses
  extern Double_t mPhi;   // MeV/c^2
  extern Double_t mK0;    // MeV/c^2
  extern Double_t mPi0;   // MeV/c^2
  extern Double_t mPiCh;  // MeV/c^2
  extern Double_t mMuon;  // MeV/c^2
  extern Double_t mElec;  // MeV/c^2
  extern Double_t mOmega; // MeV/c^2

  // Branching ratios
  // Phi
  extern Double_t br_phi_kskl;
  extern Double_t br_phi_omegapi0;

  // K-short
  extern Double_t br_ks_pi0pi0;
  extern Double_t br_ks_pippim;
  extern Double_t br_ks_pippimgamma;
  extern Double_t br_ks_piele;
  extern Double_t br_ks_pimu;

  // K-long
  extern Double_t br_kl_pi0pi0;
  extern Double_t br_kl_pippim;
  extern Double_t br_kl_pippimpi0;
  extern Double_t br_kl_3pi0;
  extern Double_t br_kl_piele;
  extern Double_t br_kl_pimu;

  // Kaons' properties and CPV
  extern Double_t tau_S_nonCPT;      // ns
  extern Double_t tau_S_CPT;         // ns
  extern Double_t tau_L;             // ns
  extern Double_t delta_mass_nonCPT; // hbar s^-1
  extern Double_t delta_mass_CPT;    // hbar s^-1
  extern Double_t mod_epsilon;
  extern Double_t Re;
  extern Double_t Im_nonCPT;     // deg
  extern Double_t Im_CPT;        // deg
  extern Double_t phi_pm_nonCPT; // deg
  extern Double_t phi_pm_CPT;    // deg
  extern Double_t phi_00_nonCPT; // deg
  extern Double_t phi_00_CPT;    // deg
}

namespace KLOE
{
  // General
  extern const Double_t T0;        // ns
  extern const UInt_t MaxNumtrkv;  //!< Maximum number of tracks in an event
  extern const UInt_t MIN_CLU_ENE; // Minimum cluster energy in MeV

  extern Int_t firstFileMax;
  extern Int_t lastFileMax;
  extern Int_t numOfThreads;

  static constexpr UInt_t channNum = 7; //!< Number of analyzed channels

  extern const std::map<Int_t, TString> channName;    //!< Map of channel names
  extern const std::map<TString, TString> channTitle; //!< Map of channel titles
  extern const std::map<TString, Color_t> channColor; //!< Map of channel colors

  extern std::map<Int_t, Int_t> channEventCount; //!< Map of event counts per channel

  Int_t TotalCountMC();

  namespace Histograms
  {
    extern const std::vector<TString> varNames; //!< Vector of variable names

    struct HistConfig1D
    {
      Int_t nBins;
      Double_t xMin;
      Double_t xMax;
      TString title;
      TString xLabel;
      TString yLabel;

      HistConfig1D(Int_t bins = 100, Double_t xMin = -10., Double_t xMax = 10., const TString &t = "", const TString &xl = "", const TString &yl = "Counts")
          : nBins(bins), xMin(xMin), xMax(xMax), title(t), xLabel(xl), yLabel(yl) {}
    }; // struct HistConfig1D

    struct HistConfig2D
    {
      Int_t nBinsX, nBinsY;
      Double_t xMin, xMax, yMin, yMax;
      TString title;
      TString xLabel, yLabel, zLabel;

      HistConfig2D(Int_t binsX = 100, Int_t binsY = 100,
                   Double_t minX = -10., Double_t maxX = 10.,
                   Double_t minY = -10., Double_t maxY = 10.,
                   const TString &t = "", const TString &xl = "",
                   const TString &yl = "", const TString &zl = "Counts")
          : nBinsX(binsX), nBinsY(binsY), xMin(minX), xMax(maxX),
            yMin(minY), yMax(maxY), title(t), xLabel(xl), yLabel(yl), zLabel(zl) {}
    };

    // Mapy konfiguracji histogramów
    extern const std::map<TString, HistConfig1D> histConfigs1D;
    extern const std::map<TString, HistConfig2D> histConfigs2D;
    extern const std::map<TString, std::pair<TString, TString>> histConfigs2D_Variables; // Para zmiennych dla 2D

    // Helper functions
    TH1F *CreateHist1D(const TString &varName, const TString &histName = "");
    TH1F *CreateHist1D(const TString &varName, const TString &channName, const TString &histName = "");

    TH2F *CreateHist2D(const TString &var1, const TString &var2, const TString &histName);
    TH2F *CreateHist2D(const TString &configName, const TString &histName);

  }

  void setGlobalStyle();

  struct BaseKinematics
  {
    BaseKinematics() : ipTriKinFit(3, 0.0),
                       KnetriKinFit(10, 0.0),
                       neuVtxTriKinFit(4, 0.0),
                       gammaMomTriKinFit1(8, 0.0),
                       gammaMomTriKinFit2(8, 0.0),
                       gammaMomTriKinFit3(8, 0.0),
                       gammaMomTriKinFit4(8, 0.0),
                       ipTriangle(3, 0.0),
                       Knerec(10, 0.0),
                       Knereclor(10, 0.0),
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
                       KnereclorFit(10, 0),
                       pi01(6, 0.0),
                       pi02(6, 0.0),
                       pi01Fit(6, 0.0),
                       pi02Fit(6, 0.0) {};

    Float_t
        Kchboost[9],
        KchboostKSOld[9],
        KchboostKLOld[9],
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
        Chi2OmegaKinFit,
        Curv1,
        Phiv1,
        Cotv1,
        Curv2,
        Phiv2,
        Cotv2;

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
        bunchnum,
        errorCode,
        doneTriKinFit,
        mcflag;

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
        pullsSignalFit,
        ParamOmega,
        ErrorsOmega,
        ParamOmegaFit,
        ErrorsOmegaFit,
        pullsOmegaFit,
        pi01,
        pi02,
        pi01Fit,
        pi02Fit,
        pi0OmegaFit[2],
        omegaFit,
        phiOmegaFit,
        trkOmegaFit[2],
        Knerec,
        Knereclor;

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
        g4takenTriKinFit,
        goodClustersTriKinFit,
        ncll;

    void resize()
    {
      clear();
      // Resize all std::vectors
      Kchrecnew.resize(10, 0.0);
      Kchboostnew.resize(10, 0.0);
      Kchboostsmeared.resize(10, 0.0);
      Kchrecsmeared.resize(10, 0.0);
      ipnew.resize(3, 0.0);
      KchrecKS.resize(10, 0.0);
      KchrecKL.resize(10, 0.0);
      KchrecClosest.resize(10, 0.0);
      KchrecKSTwoBody.resize(10, 0.0);
      KchrecKLTwoBody.resize(10, 0.0);
      Knerec.resize(10, 0.0);
      Knereclor.resize(10, 0.0);
      g4takenTriKinFit.resize(4, 0);
      goodClustersTriKinFit.resize(4, 0);
      pi01.resize(6, 0.0);
      pi02.resize(6, 0.0);
      pi01Fit.resize(6, 0.0);
      pi02Fit.resize(6, 0.0);
      ipTriKinFit.resize(3, 0.0);
      KnetriKinFit.resize(10, 0.0);
      neuVtxTriKinFit.resize(4, 0.0);
      gammaMomTriKinFit1.resize(8, 0.0);
      gammaMomTriKinFit2.resize(8, 0.0);
      gammaMomTriKinFit3.resize(8, 0.0);
      gammaMomTriKinFit4.resize(8, 0.0);
      ipTriangle.resize(3, 0.0);
      neuVtxTriangle.resize(4, 0.0);
      gammaMomTriangle1.resize(8, 0.0);
      gammaMomTriangle2.resize(8, 0.0);
      gammaMomTriangle3.resize(8, 0.0);
      gammaMomTriangle4.resize(8, 0.0);
      trcfinal.resize(4, 0.0);
      pullsTriKinFit.resize(0, 0.0);
      trkFit[0].resize(4, 0.0);
      trkFit[1].resize(4, 0.0);
      KchrecFit.resize(10, 0);
      KchboostFit.resize(10, 0);
      ipFit.resize(3, 0);
      photonFit[0].resize(4, 0.0);
      photonFit[1].resize(4, 0.0);
      photonFit[2].resize(4, 0.0);
      photonFit[3].resize(4, 0.0);
      KnerecFit.resize(10, 0.0);
      KnereclorFit.resize(10, 0.0);
      ParamSignal.resize(5, 0.0);
      ErrorsSignal.resize(5, 0.0);
      ParamSignalFit.resize(5, 0.0);
      ErrorsSignalFit.resize(5, 0.0);
      pullsSignalFit.resize(0, 0.0);
      ParamOmega.resize(5, 0.0);
      ErrorsOmega.resize(5, 0.0);
      ParamOmegaFit.resize(5, 0.0);
      ErrorsOmegaFit.resize(5, 0.0);
      pullsOmegaFit.resize(0, 0.0);
      pi0OmegaFit[0].resize(6, 0.0);
      pi0OmegaFit[1].resize(6, 0.0);
      omegaFit.resize(9, 0.0);
      phiOmegaFit.resize(4, 0.0);
      trkOmegaFit[0].resize(4, 0.0);
      trkOmegaFit[1].resize(4, 0.0);
      // Resize trk vectors
      trknew[0].resize(4, 0.0);
      trknew[1].resize(4, 0.0);
      trkKS[0].resize(4, 0.0);
      trkKS[1].resize(4, 0.0);
      trkKL[0].resize(4, 0.0);
      trkKL[1].resize(4, 0.0);
      trkClosest[0].resize(4, 0.0);
      trkClosest[1].resize(4, 0.0);
      trkKSTwoBody[0].resize(4, 0.0);
      trkKSTwoBody[1].resize(4, 0.0);
      trkKLTwoBody[0].resize(4, 0.0);
      trkKLTwoBody[1].resize(4, 0.0);
      trkMC[0].resize(4, 0.0);
      trkMC[1].resize(4, 0.0);
      trkKSmc[0].resize(4, 0.0);
      trkKSmc[1].resize(4, 0.0);
      trkKLmc[0].resize(4, 0.0);
      trkKLmc[1].resize(4, 0.0);
      trksmeared[0].resize(4, 0.0);
      trksmeared[1].resize(4, 0.0);
      CurvMC.resize(200, 0.0);
      PhivMC.resize(200, 0.0);
      CotvMC.resize(200, 0.0);
      vtaken.resize(3, 0);
      vtakenKS.resize(3, 0);
      vtakenKL.resize(3, 0);
      vtakenClosest.resize(3, 0);
    };

    /**
     * @brief Clear all data members of BaseKinematics (scalars to zero, vectors cleared/resized).
     */
    void clear()
    {
      // Clear all std::vectors
      Kchrecnew.clear();
      KchrecKS.clear();
      KchrecKL.clear();
      KchrecClosest.clear();
      KchrecKLTwoBody.clear();
      KchrecKSTwoBody.clear();
      Knerec.clear();
      Knereclor.clear();
      g4takenTriKinFit.clear();
      goodClustersTriKinFit.clear();
      pi01.clear();
      pi02.clear();
      pi01Fit.clear();
      pi02Fit.clear();
      ipTriKinFit.clear();
      KnetriKinFit.clear();
      neuVtxTriKinFit.clear();
      gammaMomTriKinFit1.clear();
      gammaMomTriKinFit2.clear();
      gammaMomTriKinFit3.clear();
      gammaMomTriKinFit4.clear();
      ipTriangle.clear();
      neuVtxTriangle.clear();
      gammaMomTriangle1.clear();
      gammaMomTriangle2.clear();
      gammaMomTriangle3.clear();
      gammaMomTriangle4.clear();
      trcfinal.clear();
      pullsTriKinFit.clear();
      trkFit[0].clear();
      trkFit[1].clear();
      KchrecFit.clear();
      KchboostFit.clear();
      ipFit.clear();
      photonFit[0].clear();
      photonFit[1].clear();
      photonFit[2].clear();
      photonFit[3].clear();
      KnerecFit.clear();
      KnereclorFit.clear();
      ParamSignal.clear();
      ErrorsSignal.clear();
      ParamSignalFit.clear();
      ErrorsSignalFit.clear();
      pullsSignalFit.clear();
      ParamOmega.clear();
      ErrorsOmega.clear();
      ParamOmegaFit.clear();
      ErrorsOmegaFit.clear();
      pullsOmegaFit.clear();
      pi0OmegaFit[0].clear();
      pi0OmegaFit[1].clear();
      omegaFit.clear();
      phiOmegaFit.clear();
      trkOmegaFit[0].clear();
      trkOmegaFit[1].clear();
      // Clear trk vectors
      trknew[0].clear();
      trknew[1].clear();
      trkKS[0].clear();
      trkKS[1].clear();
      trkKL[0].clear();
      trkKL[1].clear();
      trkClosest[0].clear();
      trkClosest[1].clear();
      trkKLTwoBody[0].clear();
      trkKLTwoBody[1].clear();
      vtaken.clear();
      vtakenKS.clear();
      vtakenKL.clear();
      vtakenClosest.clear();
    }
  };
}

namespace Utils
{
  extern nlohmann::json constants;  //!< JSON object containing analysis constants
  extern nlohmann::json properties; //!< JSON object containing analysis properties

  void InitializeVariables();
  TString elapsedTimeHMS(double totalSeconds);
}
