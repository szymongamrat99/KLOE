#include "../inc/regen_analysis.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <const.h>

#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using json = nlohmann::json;

// Konstruktor: Inicjalizuje TTreeReader i mapę histogramów
RegenAnalysis::RegenAnalysis(TChain *chain) : fReader(chain)
{
  LoadConfig();
  InitHistograms();
}

RegenAnalysis::~RegenAnalysis()
{
  // ROOT automatycznie zarządza obiektami przypisanymi do plików,
  // ale warto jawnie wyczyścić mapę dla bezpieczeństwa pamięci.
  fHistos.clear();
}

void RegenAnalysis::LoadConfig()
{
  std::ifstream configFile(Paths::regen_analysis_dir + "config/config.json");
  if (!configFile.is_open())
  {
    std::cerr << "Unable to open configuration file!" << std::endl;
    return;
  }

  json config;
  configFile >> config;

  beamPipeBound = config.value("beamPipeBound", 8.0);
  driftChamberBound = config.value("driftChamberBound", 20.0);

  histRanges = config.value("histogramRanges", std::map<std::string, std::pair<double, double>>{
                                                   {"CylindricalBeamPipe", {0.0, 10.0}},
                                                   {"SphericalBeamPipe", {0.0, 30.0}}});

  fitRanges = config.value("fitRanges", std::map<std::string, std::pair<double, double>>{
                                               {"CylindricalBeamPipe", {3.0, 5.0}},
                                               {"SphericalBeamPipe", {8.0, 13.0}}});

  resolutions = config.value("resolutions", std::map<std::string, Double_t>{
                                           {"CylindricalBeamPipe", 0.5},
                                           {"SphericalBeamPipe", 1.0}});

  // --- Histogram binning and fit parameters can also be loaded here if needed ---

  std::cout << "Configuration loaded: " << std::endl;
  std::cout << "Beam Pipe Bound: " << beamPipeBound << std::endl;
  std::cout << "Drift Chamber Bound: " << driftChamberBound << std::endl;
  for (const auto &entry : histRanges)
  {
    std::cout << "Histogram: " << entry.first << " | Range: [" << entry.second.first << ", " << entry.second.second << "]" << std::endl;
  }

  for (const auto &entry : fitRanges)
  {
    std::cout << "Fit: " << entry.first << " | Range: [" << entry.second.first << ", " << entry.second.second << "]" << std::endl;
  }

  configLoaded = true;

  configFile.close();
}

void RegenAnalysis::InitHistograms()
{
  // Automatyczne generowanie histogramów dla różnych stref regeneracji
  for (const auto &type : types)
    for (const auto &reg : regions)
    {
      Int_t nBins = (reg == "CylindricalBeamPipe") ? 100 : 150; // Przykładowa liczba binów, można dostosować

      std::string histName = "Radius_" + type + "_" + reg;
      fHistos[histName] = new TH1D(Form("h_%s", histName.c_str()),
                                   Form("Radius distribution - %s %s;Radius [cm];Events", type.c_str(), reg.c_str()),
                                   nBins, histRanges[reg].first, histRanges[reg].second);

      fHistos[histName]->GetYaxis()->SetMaxDigits(3);
    }
}

void RegenAnalysis::Loop()
{
  TTreeReaderValue<Int_t> mcflag(fReader, "mcflag");
  TTreeReaderValue<Int_t> mctruth(fReader, "mctruth");

  TTreeReaderArray<Double_t> Kchboost(fReader, "Kchboost"); // Non fitted
  TTreeReaderArray<Double_t> Knerec(fReader, "Knerec");     // Non fitted
  TTreeReaderArray<Double_t> ip(fReader, "ip");

  TTreeReaderArray<Double_t> KchrecFit(fReader, "KchrecFit"); // Fitted
  TTreeReaderArray<Double_t> KnerecFit(fReader, "KnerecFit");     // Fitted
  TTreeReaderArray<Double_t> KnerecSix(fReader, "KnerecSix");     // Fitted with 6 parameters
  TTreeReaderArray<Double_t> trk1Fit(fReader, "trk1Fit");
  TTreeReaderArray<Double_t> trk2Fit(fReader, "trk2Fit");
  TTreeReaderValue<Double_t> minv4gam(fReader, "minv4gam");
  TTreeReaderValue<Double_t> KaonChTimeCMSignalFit(fReader, "KaonChTimeCMSignalFit");
  TTreeReaderValue<Double_t> KaonNeTimeCMSignalFit(fReader, "KaonNeTimeCMSignalFit");

  TTreeReaderValue<Double_t> Bx(fReader, "Bx");
  TTreeReaderValue<Double_t> By(fReader, "By");
  TTreeReaderValue<Double_t> Bz(fReader, "Bz");
  TTreeReaderArray<Double_t> KchrecClosest(fReader, "KchrecClosest");

  std::cout << "Starting event loop..." << std::endl;

  double u0 = -92.9524, v0 = 7.1644, su = 33.7633, sv = 42.0163, rho = 0.556;
  double n_sigma_cut = 2.5;

  // Definicja lambdy
  auto ellipse_cut = [=](double u, double v)
  {
    double du = u - u0;
    double dv = v - v0;
    double inv_rho2 = 1.0 / (1.0 - rho * rho);

    double dist2 = inv_rho2 * ((du * du) / (su * su) +
                               (dv * dv) / (sv * sv) -
                               2.0 * rho * (du * dv) / (su * sv));

    return std::sqrt(std::max(0.0, dist2)) > n_sigma_cut;
  };

  while (fReader.Next())
  {
    if (*mctruth != 2 || *mcflag != 1) // Interesują nas tylko zdarzenia regeneracji
      continue;

    // Cut setting
    std::array<Double_t, 3> distNeutralCharged = {KchrecFit[6] - KnerecFit[6],
                                                  KchrecFit[7] - KnerecFit[7],
                                                  KchrecFit[8] - KnerecFit[8]},
                            distNeutralIP = {Knerec[6] - *Bx,
                                             Knerec[7] - *By,
                                             Knerec[8] - KchrecClosest[8]},
                            distChargedIP = {KchrecClosest[6] - *Bx,
                                             KchrecClosest[7] - *By,
                                             KchrecClosest[8] - *Bz};

    Float_t
        rho_pm = std::sqrt(distChargedIP[0] * distChargedIP[0] + distChargedIP[1] * distChargedIP[1]),
        rho_00 = std::sqrt(distNeutralIP[0] * distNeutralIP[0] + distNeutralIP[1] * distNeutralIP[1]),
        rho = std::sqrt(std::pow(rho_pm, 2) + std::pow(rho_00, 2));

    Double_t radius00 = std::sqrt(std::pow(Knerec[6] - *Bx, 2) +
                                  std::pow(Knerec[7] - *By, 2)),
             radiuspm = std::sqrt(std::pow(KchrecClosest[6] - *Bx, 2) +
                                  std::pow(KchrecClosest[7] - *By, 2)),
             zdist00 = std::abs(Knerec[8] - KchrecClosest[8]),
             zdistpm = std::abs(KchrecClosest[8] - *Bz);

    Double_t fiducialVolume = std::sqrt(std::pow(distNeutralCharged[0], 2) + std::pow(distNeutralCharged[1], 2)) < 2.05 && std::abs(distNeutralCharged[2]) < 2.45,
             fiducialVolumeClose = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

    TVector3 KchrecVec = {KchrecFit[0], KchrecFit[1], KchrecFit[2]};
    TVector3 trk1VecFit = {trk1Fit[0], trk1Fit[1], trk1Fit[2]};
    TVector3 trk2VecFit = {trk2Fit[0], trk2Fit[1], trk2Fit[2]};

    Double_t phiTrk1Angle = cos(trk1VecFit.Angle(KchrecVec)),
             phiTrk2Angle = cos(trk2VecFit.Angle(KchrecVec));

    Bool_t global_cut = ((fiducialVolume && (abs(phiTrk1Angle) < 0.8 || abs(phiTrk2Angle) < 0.8)) || !fiducialVolume) && ((fiducialVolumeClose && rho > 1.5) || !fiducialVolumeClose) && ellipse_cut(*minv4gam - PhysicsConstants::mK0, KnerecSix[5] - PhysicsConstants::mK0);

    Double_t Dtboostlor = *KaonChTimeCMSignalFit - *KaonNeTimeCMSignalFit;

    if (!global_cut)
      continue;

    // Set up geometry vectors
    ipVector.SetXYZ(ip[0], ip[1], ip[2]);
    neutralVtxVector.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);
    chargedVtxVector.SetXYZ(Kchboost[6], Kchboost[7], Kchboost[8]);

    // Calculate radius vectors
    TVector3 ipToCharged = chargedVtxVector - ipVector;
    TVector3 ipToNeutral = neutralVtxVector - ipVector;

    // Cylindrical radius (rho) to charged and neutral vertices
    double rhoCharged = ipToCharged.Perp();
    double rhoNeutral = ipToNeutral.Perp();

    // Spherical radius to charged and neutral vertices
    double rCharged = ipToCharged.Mag();
    double rNeutral = ipToNeutral.Mag();

    // Klasyfikacja zdarzenia do odpowiedniego histogramu na podstawie geometrii detektora

    if (*KaonChTimeCMSignalFit > 7 * 1.0)
    {
      // if (rhoCharged < beamPipeBound)
      fHistos["Radius_Charged_CylindricalBeamPipe"]->Fill(rhoCharged);
      // else if (rhoCharged >= beamPipeBound && rhoCharged < driftChamberBound)
      fHistos["Radius_Charged_SphericalBeamPipe"]->Fill(rCharged);
    }

    if (*KaonNeTimeCMSignalFit > 7 * 1.0)
    {
      // if (rhoNeutral < beamPipeBound)
      fHistos["Radius_Neutral_CylindricalBeamPipe"]->Fill(rhoNeutral);
      // else if (rhoNeutral >= beamPipeBound && rhoNeutral < driftChamberBound)
      fHistos["Radius_Neutral_SphericalBeamPipe"]->Fill(rNeutral);
    }
  }
}

void RegenAnalysis::FitResults()
{
  if (!configLoaded)
  {
    std::cerr << "Configuration not loaded. Cannot perform fits." << std::endl;
    return;
  }

  TCanvas *c = new TCanvas("c_fits", "Fit Results", 800, 600);
  c->Print(Paths::regen_analysis_dir + Paths::img_dir + "Regen_Fits_Report.pdf["); // Otwarcie wielostronicowego PDF

  TFitResultPtr r;

  for (const auto &type : types)
    for (const auto &reg : regions)
    {
      std::string histName = "Radius_" + type + "_" + reg;
      if (fHistos[histName]->GetEntries() > 0)
      {
        Double_t fitMin = fitRanges[reg].first;
        Double_t fitMax = fitRanges[reg].second;

        r = fHistos[histName]->Fit("gaus", "SMQ", "", fitMin, fitMax);

        if (r.Get())
        {
          std::cout << "Fit successful for " << histName << ": Mean = " << r->Parameter(1) << ", Sigma = " << r->Parameter(2) << std::endl;

          fHistos[histName]->Draw();
          c->Print(Paths::regen_analysis_dir + Paths::img_dir + "Regen_Fits_Report.pdf"); // Dodanie strony do PDF
        }
        else
        {
          std::cerr << "Fit failed for " << histName << std::endl;
        }
      }
    }

  c->Print(Paths::regen_analysis_dir + Paths::img_dir + "Regen_Fits_Report.pdf]"); // Zamknięcie wielostronicowego PDF
  delete c;
}