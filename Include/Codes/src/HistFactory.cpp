#include "TFractionFitter.h"
#include "TObjArray.h"

#include "HistFactory.h"
#include "const.h"

HistFactory::HistFactory(const HistConfig &config)
{
  addHistogram(config);
}

HistFactory::HistFactory(const std::vector<HistConfig> &configs)
{
  for (const auto &config : configs)
    addHistogram(config);
}

void HistFactory::fill(const TString &id, Double_t x, Double_t weight)
{
  TH1 *hist = get1DHistogram(id);
  if (hist)
    hist->Fill(x, weight);
}

void HistFactory::fill(const TString &id, const Int_t mctruth, Double_t x, Double_t weight)
{
  auto it = KLOE::channName.find(mctruth);
  if (it == KLOE::channName.end())
    return;

  TH1 *hist = get1DHistogram(id, it->second);
  if (hist)
    hist->Fill(x, weight);
}

void HistFactory::fill(const TString &id, Double_t x, Double_t y, Double_t weight)
{
  TH2 *hist = get2DHistogram(id);
  if (hist)
    hist->Fill(x, y, weight);
}

void HistFactory::fill(const TString &id, const Int_t mctruth, Double_t x, Double_t y, Double_t weight)
{
  auto it = KLOE::channName.find(mctruth);
  if (it == KLOE::channName.end())
    return;

  TH2 *hist = get2DHistogram(id, it->second);
  if (hist)
    hist->Fill(x, y, weight);
}

std::map<TString, FitResult> HistFactory::fitComponentsToData(const TString &id, const std::map<TString, std::pair<Double_t, Double_t>> &customLimits)
{
  std::map<TString, FitResult> results;

  // 1. Pobieramy dane
  TH1 *histData = get1DHistogram(id, "Data");
  if (!histData || histData->Integral() <= 0)
  {
    std::cout << "[HistFactory::fitComponentsToData] Error: No data for id " << id << std::endl;
    return results;
  }

  // 2. Przygotowujemy listę szablonów MC
  TObjArray *mcTemplates = new TObjArray();
  std::vector<TString> activeChannelNames;

  TString idStr = id.Data();
  auto itId = histograms1DByChannel.find(idStr);
  if (itId == histograms1DByChannel.end())
    return results;

  for (auto &chPair : itId->second)
  {
    // Pomijamy Dane i stary MC sum
    if (chPair.first == "Data" || chPair.first == "MC sum")
      continue;

    TH1 *histMC = chPair.second;
    if (histMC && histMC->Integral() > 0)
    {
      mcTemplates->Add(histMC);
      activeChannelNames.push_back(chPair.first); // Zapamiętujemy kolejność
    }
  }

  if (mcTemplates->GetEntries() == 0)
  {
    delete mcTemplates;
    return results;
  }

  // 3. Konfiguracja i uruchomienie Fittera
  TFractionFitter *fitter = new TFractionFitter(histData, mcTemplates, "Q"); // "Q" - quiet mode

  for (Int_t i = 0; i < mcTemplates->GetEntries(); ++i)
  {
    TString chName = activeChannelNames[i];

    // Szukamy, czy użytkownik podał własne limity dla tego kanału
    auto it = customLimits.find(chName);

    if (it != customLimits.end())
    {
      Double_t minLimit = it->second.first;
      Double_t maxLimit = it->second.second;
      fitter->Constrain(i, minLimit, maxLimit);
      std::cout << "  [Fit Config] Set custom limit for " << chName
                << ": [" << minLimit << ", " << maxLimit << "]" << std::endl;
    }
    else
    {
      // Domyślny, bezpieczny zakres dla fizycznych frakcji
      fitter->Constrain(i, 0.0, 1.0);
    }
  }

  std::cout << "[HistFactory] Fiting components for ID: " << id << "..." << std::endl;
  Int_t status = fitter->Fit(); // Wykonanie dopasowania

  if (status == 0) // Status 0 oznacza sukces MINUIT-a
  {
    // Całkowita liczba zliczeń w danych
    Double_t nData = histData->Integral();

    // 4. Pobieranie wyników i skalowanie histogramów
    for (Int_t i = 0; i < mcTemplates->GetEntries(); ++i)
    {
      Double_t fraction, errorFraction;
      fitter->GetResult(i, fraction, errorFraction);

      TH1 *histMC = (TH1 *)mcTemplates->At(i);
      TString chName = activeChannelNames[i];

      // TFractionFitter zwraca FRAKCJĘ (udział w danych).
      // Czynnik skalujący dla danego histogramu to: (Dane_Total * Frakcja) / MC_Component_Total
      Double_t scaleFactor = (nData * fraction) / histMC->Integral();
      Double_t scaleError = (nData * errorFraction) / histMC->Integral();

      // Skalujemy komponent w fabryce nowym współczynnikiem z fita
      histMC->Scale(scaleFactor, "I");

      // Zapisujemy wyniki do zwrócenia użytkownikowi
      results[chName] = {scaleFactor, scaleError};

      std::cout << "  -> Channel [" << chName << "]: Scale Factor = "
                << scaleFactor << " +/- " << scaleError << " (Frac: " << fraction * 100 << "%)" << std::endl;
    }

    // 5. Na koniec przebudowujemy MC sum, żeby zawierał już dofitowane komponenty
    buildMCSum(id);
  }
  else
  {
    std::cout << "[HistFactory::fitComponentsToData] Error: Fit failed with status " << status << std::endl;
  }

  // Czyszczenie pamięci obiektów pomocniczych
  delete fitter;
  delete mcTemplates;

  return results;
}

// ---

void HistFactory::addHistogram(const HistConfig &config)
{
  if (config.is2D)
  {
    if (config.byChannel)
      for (const auto &channel : KLOE::channName)
        histograms2DByChannel[config.id][channel.second] = create2DHistogram(config, channel.second);
    else
      simple2DHistograms[config.id] = create2DHistogram(config);
  }
  else if (config.byChannel)
    for (const auto &channel : KLOE::channName)
    {
      histograms1DByChannel[config.id][channel.second] = create1DHistogram(config, channel.second);
      setChannColors(config.id, channel.second);
    }
  else
    simple1DHistograms[config.id] = create1DHistogram(config);
}

TH1 *HistFactory::create1DHistogram(const HistConfig &config, const TString &channel)
{
  TString id = channel.IsNull() ? config.id : config.id + "_" + channel;
  TH1 *hist = new TH1D(id, config.title.c_str(),
                       config.nBinsX, config.xMin, config.xMax);
  hist->GetXaxis()->SetTitle(config.xAxisTitle.c_str());
  hist->GetYaxis()->SetTitle(config.yAxisTitle.c_str());
  return hist;
}

TH2 *HistFactory::create2DHistogram(const HistConfig &config, const TString &channel)
{
  TString id = channel.IsNull() ? config.id : config.id + "_" + channel;
  TH2 *hist = new TH2D(id, config.title.c_str(),
                       config.nBinsX, config.xMin, config.xMax,
                       config.nBinsY, config.yMin, config.yMax);
  hist->GetXaxis()->SetTitle(config.xAxisTitle.c_str());
  hist->GetYaxis()->SetTitle(config.yAxisTitle.c_str());
  return hist;
}

TH1 *HistFactory::get1DHistogram(const TString &id, const TString &channel)
{
  std::string idStr = id.Data();

  if (channel.IsNull())
  {
    auto it = simple1DHistograms.find(idStr);
    if (it != simple1DHistograms.end())
      return it->second;
  }
  else
  {
    auto itId = histograms1DByChannel.find(idStr);
    if (itId != histograms1DByChannel.end())
    {
      auto itCh = itId->second.find(channel.Data());
      if (itCh != itId->second.end())
        return itCh->second;
    }
  }
  return nullptr;
}

TH2 *HistFactory::get2DHistogram(const TString &id, const TString &channel)
{
  std::string idStr = id.Data();

  if (channel.IsNull())
  {
    auto it = simple2DHistograms.find(idStr);
    if (it != simple2DHistograms.end())
      return it->second;
  }
  else
  {
    auto itId = histograms2DByChannel.find(idStr);
    if (itId != histograms2DByChannel.end())
    {
      auto itCh = itId->second.find(channel.Data());
      if (itCh != itId->second.end())
        return itCh->second;
    }
  }
  return nullptr;
}

void HistFactory::setChannColors(const TString &id, const TString &channel)
{
  auto it = KLOE::channColor.find(channel.Data());
  Color_t color = (it != KLOE::channColor.end()) ? it->second : kBlack;

  TH1 *hist = get1DHistogram(id, channel);
  if (hist)
  {
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
  }
}

void HistFactory::buildMCSum(const TString &id)
{
  TH1 *histMCSum = get1DHistogram(id);
  histMCSum->Reset(); // Clear the histogram before summing
  if (!histMCSum)
    return;

  for (const auto &channel : KLOE::channName)
  {
    if (channel.second == "Data" || channel.second == "MC sum")
      continue;

    TH1 *hist = get1DHistogram(id, channel.second);
    if (!hist)
      continue;

    histMCSum->Add(hist);
  }
}

void HistFactory::buildMCSum2D(const TString &id)
{
  TH2 *histMCSum = get2DHistogram(id);
  histMCSum->Reset(); // Clear the histogram before summing
  if (!histMCSum)
    return;

  for (const auto &channel : KLOE::channName)
  {
    if (channel.second == "Data" || channel.second == "MC sum")
      continue;

    TH2 *hist = get2DHistogram(id, channel.second);
    if (!hist)
      continue;

    histMCSum->Add(hist);
  }
}

void HistFactory::normalize(Normalization norm)
{
  for (auto &idPair : histograms1DByChannel)
  {
    TH1 *histData = get1DHistogram(idPair.first, "Data");
    TH1 *histMCSum = get1DHistogram(idPair.first, "MC sum");

    Double_t scaleFactor = histData && histData->Integral() > 0 ? histData->Integral() / histMCSum->Integral() : 1.0;

    for (auto &chPair : idPair.second)
    {
      if (chPair.first == "Data")
        continue;

      TH1 *hist = chPair.second;
      if (!hist)
        continue;

      if (norm == Normalization::MCToData)
        hist->Scale(scaleFactor, "I");
    }

    if (histMCSum)
      if (norm == Normalization::MCToData)
        histMCSum->Scale(scaleFactor, "I");
  }

  for (auto &idPair : histograms2DByChannel)
  {
    TH2 *histData = get2DHistogram(idPair.first, "Data");
    TH2 *histMCSum = get2DHistogram(idPair.first, "MC sum");

    Double_t scaleFactor = histData && histData->Integral() > 0 ? histData->Integral() / histMCSum->Integral() : 1.0;

    for (auto &chPair : idPair.second)
    {
      if (chPair.first == "Data")
        continue;

      TH2 *hist = chPair.second;
      if (!hist)
        continue;

      if (norm == Normalization::MCToData)
        hist->Scale(scaleFactor, "I");
    }

    if (histMCSum)
      if (norm == Normalization::MCToData)
        histMCSum->Scale(scaleFactor, "I");
  }
}