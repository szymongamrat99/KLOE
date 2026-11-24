#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TTree.h>
#include <string>

void cpfit_test_data(const char *directory)
{
  const std::string base(directory);
  {
    TFile file((base + "/events.root").c_str(), "RECREATE");
    TTree tree("analysis", "analysis");
    int mcflag = 0, mctruth = 0;
    double minv4gam = 0, bx = 0, by = 0, bz = 0;
    double tchmc = 0, tnemc = 0, tchfit = 0, tnefit = 0, tchtri = 0, tnetri = 0;
    double kne[9] = {}, kne6[9] = {}, kchclosest[9] = {}, knefit[9] = {}, kchfit[9] = {};
    double kch[9] = {}, trk1[3] = {0, 1, 0}, trk2[3] = {0, 1, 0};
    double ip[3] = {};
    tree.Branch("mcflag", &mcflag, "mcflag/I"); tree.Branch("mctruth", &mctruth, "mctruth/I");
    tree.Branch("minv4gam", &minv4gam, "minv4gam/D");
    tree.Branch("Bx", &bx, "Bx/D"); tree.Branch("By", &by, "By/D"); tree.Branch("Bz", &bz, "Bz/D");
    tree.Branch("KaonChTimeCMMC", &tchmc, "KaonChTimeCMMC/D");
    tree.Branch("KaonNeTimeCMMC", &tnemc, "KaonNeTimeCMMC/D");
    tree.Branch("KaonChTimeCMSignalFit", &tchfit, "KaonChTimeCMSignalFit/D");
    tree.Branch("KaonNeTimeCMSignalFit", &tnefit, "KaonNeTimeCMSignalFit/D");
    tree.Branch("KaonChTimeCMBoostTriFit", &tchtri, "KaonChTimeCMBoostTriFit/D");
    tree.Branch("KaonNeTimeCMBoostTriFit", &tnetri, "KaonNeTimeCMBoostTriFit/D");
    tree.Branch("Knerec", kne, "Knerec[9]/D"); tree.Branch("KnerecSix", kne6, "KnerecSix[9]/D");
    tree.Branch("KchrecClosest", kchclosest, "KchrecClosest[9]/D");
    tree.Branch("KnerecFit", knefit, "KnerecFit[9]/D"); tree.Branch("KchrecFit", kchfit, "KchrecFit[9]/D");
    tree.Branch("Kchrec", kch, "Kchrec[9]/D");
    tree.Branch("ip", ip, "ip[3]/D");
    tree.Branch("trk1Fit", trk1, "trk1Fit[3]/D"); tree.Branch("trk2Fit", trk2, "trk2Fit[3]/D");
    kchfit[0] = 1;
    kne6[5] = 497.611 + 7.1644;
    ip[0] = 1.0;
    for (int bin = 0; bin < 7; ++bin)
    {
      tchfit = -2.57 + bin * (5.14 / 6);
      tchtri = tchfit;
      tchmc = tchfit + 4;
      for (int channel = 0; channel <= 6; ++channel)
      {
        mcflag = channel == 0 ? 0 : 1;
        mctruth = channel;
        kch[6] = channel == 1 ? 3.0 : (channel == 2 ? 4.8 : 0.0);
        kne[6] = channel == 2 ? 4.8 : 0.0;
        const int count = channel == 0 ? 120 : (channel == 1 ? 80 : 8);
        for (int i = 0; i < count; ++i)
        {
          minv4gam = 497.611 - 92.9524 + (i % 2 == 0 ? 0 : 200);
          tree.Fill();
        }
      }
    }
    tree.Write();
  }
  {
    TFile file((base + "/signal.root").c_str(), "RECREATE");
    TH1D weights("weights", "weights", 10, -5, 5);
    for (int i = 1; i <= 10; ++i)
    {
      weights.SetBinContent(i, 0.96);
      weights.SetBinError(i, i > 5 ? 0.5 : 0.01);
    }
    weights.Write();

    TH1D ksPM("ks_pm_weights", "ks_pm_weights", 10, -5, 5);
    TH1D klPM("kl_pm_weights", "kl_pm_weights", 10, -5, 5);
    TH1D ks00("ks_00_weights", "ks_00_weights", 10, -5, 5);
    for (int i = 1; i <= 10; ++i)
    {
      ksPM.SetBinContent(i, 0.96);
      ksPM.SetBinError(i, i > 5 ? 0.5 : 0.01);
      klPM.SetBinContent(i, 0.8);
      klPM.SetBinError(i, 0.01);
      ks00.SetBinContent(i, 0.7);
      ks00.SetBinError(i, 0.01);
    }
    ksPM.Write();
    klPM.Write();
    ks00.Write();
  }
  {
    TFile file((base + "/signal_ratios.root").c_str(), "RECREATE");
    const double x[] = {0.0, 10.0, 20.0};
    const double ksValues[] = {1.0, 1.5, 2.0};
    const double klValues[] = {0.5, 1.0, 1.5};
    TGraph ksGraph(3, x, ksValues);
    TGraph klGraph(3, x, klValues);
    TSpline3 ksSpline("spline_1", &ksGraph);
    TSpline3 klSpline("spline_2", &klGraph);
    ksSpline.Write("spline_1");
    klSpline.Write("spline_2");
  }
  {
    TFile file((base + "/regen.root").c_str(), "RECREATE");
    for (const char *name : {"rhoCharged", "rhoNeutral", "RCharged", "RNeutral"})
      for (const char *prefix : {"FuncDC_", "FuncCylBP_", "FuncSphBP_"})
      {
        TF1 function((std::string(prefix) + name).c_str(), "pol1", 0, 50);
        function.SetParameters(0.8, 0);
        function.Write();
      }
  }
}
