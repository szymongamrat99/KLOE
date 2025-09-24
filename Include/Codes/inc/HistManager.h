#pragma once
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TPad.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TFractionFitter.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TF1.h>
#include <TLatex.h>
#include "triple_gaus.h"
#include <vector>
#include <map>

/**
 * @class HistManager
 * @brief Klasa zarządzająca histogramami 1D i 2D z obsługą FractionFit
 * 
 * Klasa umożliwia tworzenie, wypełnianie, konfigurowanie i rysowanie histogramów.
 * Obsługuje zarówno dane Monte Carlo (MC) jak i dane eksperymentalne.
 * Zawiera zaawansowaną funkcjonalność FractionFit do dopasowywania składowych MC do danych.
 * 
 * @section example Przykład użycia:
 * @code
 * // Inicjalizacja z 7 kanałami MC
 * Color_t colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet};
 * std::vector<TString> names = {"Signal", "Background1", "Background2", "Background3", 
 *                               "Background4", "Background5", "Background6"};
 * HistManager* histMgr = new HistManager(7, colors, names);
 * 
 * // Konfiguracja histogramu
 * HistManager::HistConfig config;
 * config.name = "deltaT";
 * config.title = "Time difference";
 * config.bins = 100;
 * config.xmin = -10.0;
 * config.xmax = 10.0;
 * config.xtitle = "#Delta t [ns]";
 * config.ytitle = "Events";
 * 
 * histMgr->CreateHistSet1D("deltaT", config);
 * 
 * // Wypełnianie histogramów MC
 * for(int i = 0; i < nevents; ++i) {
 *     histMgr->Fill1D("deltaT", mctruth, deltaT_value, weight);
 * }
 * 
 * // Wypełnianie danych eksperymentalnych
 * for(int i = 0; i < ndata; ++i) {
 *     histMgr->FillData1D("deltaT", deltaT_data);
 * }
 * 
 * // Włączenie FractionFit
 * histMgr->SetUseFractionFitter(true);
 * // LUB proste skalowanie
 * histMgr->SetNormalizationType(HistManager::NormalizationType::SIMPLE_SCALE);
 * 
 * // Opcjonalne ustawienie ograniczeń
 * HistManager::FitConstraints constraints(7);
 * constraints.lowerBounds = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
 * constraints.upperBounds = {0.5, 0.3, 0.8, 0.6, 0.4, 0.2, 0.1};
 * constraints.fitRangeMin = 10;  // Użyj binów 10-90 do fitu
 * constraints.fitRangeMax = 90;
 * histMgr->SetFitConstraints("deltaT", constraints);
 * 
 * // Rysowanie z automatycznym fitem
 * histMgr->DrawSet1D("deltaT", "HIST", true);
 * 
 * // Sprawdzenie wyników fitu
 * auto fitResult = histMgr->GetLastFitResult();
 * if(fitResult.converged) {
 *     std::cout << "Fit converged with chi2/ndf = " << fitResult.chi2/fitResult.ndf << std::endl;
 *     for(int i = 0; i < 7; ++i) {
 *         std::cout << names[i] << ": " << fitResult.fractions[i] 
 *                   << " +/- " << fitResult.errors[i] << std::endl;
 *     }
 * }
 * 
 * // Ręczne wykonanie fitu (opcjonalne)
 * fitResult = histMgr->PerformFractionFit("deltaT", false); // bez zapisanych ograniczeń
 * @endcode
 */
class HistManager {
public:
    /**
     * @struct HistConfig
     * @brief Konfiguracja histogramu 1D
     */
    struct HistConfig {
        TString name;                ///< Nazwa histogramu
        TString title;               ///< Tytuł histogramu
        Int_t bins;                  ///< Liczba binów
        Double_t xmin;               ///< Minimalna wartość osi X
        Double_t xmax;               ///< Maksymalna wartość osi X
        TString xtitle;              ///< Tytuł osi X
        TString ytitle = "Counts";   ///< Tytuł osi Y
        Bool_t logy = false;         ///< Skala logarytmiczna na osi Y
        Bool_t showStats = true;     ///< Wyświetlanie statbox
    };

    /**
     * @struct Hist2DConfig
     * @brief Konfiguracja histogramu 2D (rozszerza HistConfig)
     */
    struct Hist2DConfig : public HistConfig {
        Int_t binsy;                 ///< Liczba binów na osi Y
        Double_t ymin;               ///< Minimalna wartość osi Y
        Double_t ymax;               ///< Maksymalna wartość osi Y
    };

    /**
     * @struct FitConstraints
     * @brief Ograniczenia dla parametrów FractionFit
     */
    struct FitConstraints {
        std::vector<Double_t> lowerBounds;  ///< Dolne granice parametrów (0.0-1.0)
        std::vector<Double_t> upperBounds;  ///< Górne granice parametrów (0.0-1.0)
        std::vector<Double_t> initialValues; ///< Wartości początkowe parametrów
        Int_t fitRangeMin = -1;             ///< Minimalny bin dla fitu (-1 = cały zakres)
        Int_t fitRangeMax = -1;             ///< Maksymalny bin dla fitu (-1 = cały zakres)
        
        /**
         * @brief Konstruktor z domyślnymi wartościami
         * @param nChannels Liczba kanałów MC
         */
        FitConstraints(Int_t nChannels = 0) {
            if(nChannels > 0) {
                lowerBounds.assign(nChannels, 0.0);
                upperBounds.assign(nChannels, 2.0);
                initialValues.assign(nChannels, 1.0/nChannels); // Równomierne wartości początkowe
            }
        }
    };

    /**
     * @struct FitResult
     * @brief Wyniki FractionFit
     */
    struct FitResult {
        std::vector<Double_t> fractions;    ///< Frakcje składowych MC
        std::vector<Double_t> errors;       ///< Błędy frakcji
        Double_t chi2;                      ///< Wartość chi-kwadrat
        Int_t ndf;                          ///< Liczba stopni swobody
        Int_t status;                       ///< Status fitu (0 = sukces)
        Bool_t converged;                   ///< Czy fit się zbiegł
        TString fitInfo;                    ///< Dodatkowe informacje o ficie
        
        FitResult() : chi2(0), ndf(0), status(-1), converged(false) {}
    };

    /**
     * @struct ArrayConfig
     * @brief Konfiguracja dla zestawów histogramów typu array (np. składowe pędu, pulls[i])
     */
    struct ArrayConfig {
        TString baseName;               ///< Bazowa nazwa zestawu (np. "momentum", "pulls")
        TString baseTitle;              ///< Bazowy tytuł (np. "Momentum components", "Pull values")
        Int_t arraySize;                ///< Rozmiar tablicy/wektora
        std::vector<TString> varNames;  ///< Nazwy zmiennych (opcjonalne, np. {"px", "py", "pz", "E"})
        std::vector<TString> varTitles; ///< Tytuły zmiennych (opcjonalne, np. {"p_{x}", "p_{y}", "p_{z}", "E"})
        HistConfig commonConfig;        ///< Wspólna konfiguracja dla wszystkich histogramów w zestawie
        
        // Indywidualne konfiguracje dla każdego sub-histogramu (opcjonalne)
        std::vector<Int_t> individualBins;      ///< Liczba binów dla każdego sub-histogramu
        std::vector<Double_t> individualXmin;   ///< Minimalne wartości X dla każdego sub-histogramu  
        std::vector<Double_t> individualXmax;   ///< Maksymalne wartości X dla każdego sub-histogramu
        std::vector<TString> individualXtitle;  ///< Tytuły osi X dla każdego sub-histogramu
        std::vector<TString> individualYtitle;  ///< Tytuły osi Y dla każdego sub-histogramu
        
        ArrayConfig() : arraySize(0) {}
        
        ArrayConfig(const TString& name, const TString& title, Int_t size, const HistConfig& config) 
            : baseName(name), baseTitle(title), arraySize(size), commonConfig(config) {}
            
        /**
         * @brief Ustawia indywidualne parametry binowania dla sub-histogramu
         * @param index Indeks sub-histogramu (0-based)
         * @param bins Liczba binów
         * @param xmin Minimalna wartość X
         * @param xmax Maksymalna wartość X
         */
        void SetIndividualBinning(Int_t index, Int_t bins, Double_t xmin, Double_t xmax) {
            if(index < 0) return;
            
            // Rozszerz wektory jeśli potrzeba
            if(index >= static_cast<Int_t>(individualBins.size())) {
                individualBins.resize(index + 1, commonConfig.bins);
                individualXmin.resize(index + 1, commonConfig.xmin);
                individualXmax.resize(index + 1, commonConfig.xmax);
            }
            
            individualBins[index] = bins;
            individualXmin[index] = xmin;
            individualXmax[index] = xmax;
        }
        
        /**
         * @brief Ustawia indywidualne tytuły osi dla sub-histogramu
         * @param index Indeks sub-histogramu (0-based)
         * @param xtitle Tytuł osi X
         * @param ytitle Tytuł osi Y (opcjonalny)
         */
        void SetIndividualAxisTitles(Int_t index, const TString& xtitle, const TString& ytitle = "") {
            if(index < 0) return;
            
            // Rozszerz wektory jeśli potrzeba
            if(index >= static_cast<Int_t>(individualXtitle.size())) {
                individualXtitle.resize(index + 1, commonConfig.xtitle);
            }
            if(index >= static_cast<Int_t>(individualYtitle.size())) {
                individualYtitle.resize(index + 1, commonConfig.ytitle);
            }
            
            individualXtitle[index] = xtitle;
            if(!ytitle.IsNull()) {
                individualYtitle[index] = ytitle;
            }
            // Jeśli ytitle jest puste, pozostaw dotychczasową wartość (domyślną lub wcześniej ustawioną)
        }
        
        /**
         * @brief Pobiera konfigurację binowania dla danego indeksu
         * @param index Indeks sub-histogramu
         * @param bins Wyjście: liczba binów
         * @param xmin Wyjście: minimalna wartość X
         * @param xmax Wyjście: maksymalna wartość X
         */
        void GetBinning(Int_t index, Int_t& bins, Double_t& xmin, Double_t& xmax) const {
            if(index >= 0 && index < static_cast<Int_t>(individualBins.size())) {
                bins = individualBins[index];
                xmin = individualXmin[index];
                xmax = individualXmax[index];
            } else {
                bins = commonConfig.bins;
                xmin = commonConfig.xmin;
                xmax = commonConfig.xmax;
            }
        }
        
        /**
         * @brief Pobiera tytuły osi dla danego indeksu
         * @param index Indeks sub-histogramu
         * @param xtitle Wyjście: tytuł osi X
         * @param ytitle Wyjście: tytuł osi Y
         */
        void GetAxisTitles(Int_t index, TString& xtitle, TString& ytitle) const {
            // Pobierz tytuł X
            if(index >= 0 && index < static_cast<Int_t>(individualXtitle.size())) {
                xtitle = individualXtitle[index];
            } else {
                xtitle = commonConfig.xtitle;
            }
            
            // Pobierz tytuł Y (może mieć inny rozmiar niż X titles)
            if(index >= 0 && index < static_cast<Int_t>(individualYtitle.size())) {
                ytitle = individualYtitle[index];
            } else {
                ytitle = commonConfig.ytitle;
            }
        }
    };

    /**
     * @struct FitParams3Gauss
     * @brief Struktura przechowująca parametry i wyniki fitowania Triple Gaussian
     */
    struct FitParams3Gauss {
        // Parametry Triple Gaussian: 9 parametrów (A1,μ1,σ1, A2,μ2,σ2, A3,μ3,σ3)
        std::vector<Double_t> initialParams;       ///< Wartości początkowe parametrów [9]
        std::vector<Double_t> lowerBounds;         ///< Dolne granice parametrów [9]
        std::vector<Double_t> upperBounds;         ///< Górne granice parametrów [9]
        std::vector<TString> paramNames;           ///< Nazwy parametrów
        
        // Wyniki fitu
        std::vector<Double_t> fittedParams;        ///< Wyniki fitu [9]
        std::vector<Double_t> paramErrors;         ///< Błędy parametrów [9]
        Double_t chi2;                             ///< Chi-kwadrat
        Int_t ndf;                                 ///< Liczba stopni swobody
        Int_t status;                              ///< Status fitu
        Bool_t converged;                          ///< Czy fit się zbiegł
        
        // Obliczone wartości kombinowane (z metod triple_gaus)
        Double_t combinedMean;                     ///< Średnia kombinowana
        Double_t combinedMeanErr;                  ///< Błąd średniej kombinowanej
        Double_t combinedStdDev;                   ///< Odchylenie standardowe kombinowane
        Double_t combinedStdDevErr;                ///< Błąd odchylenia standardowego
        
        // Zakres fitu
        Double_t fitRangeMin;                      ///< Dolna granica zakresu fitu
        Double_t fitRangeMax;                      ///< Górna granica zakresu fitu
        
        TF1* fitFunction;                          ///< Wskaźnik na funkcję fitu
        
        FitParams3Gauss() {
            InitializeDefaults();
        }
        
        ~FitParams3Gauss() {
            if(fitFunction) {
                delete fitFunction;
                fitFunction = nullptr;
            }
        }
        
        /**
         * @brief Inicjalizuje domyślne wartości parametrów
         */
        void InitializeDefaults() {
            // Nazwy parametrów
            paramNames = {"A1", "μ1", "σ1", "A2", "μ2", "σ2", "A3", "μ3", "σ3"};
            
            // Domyślne wartości początkowe (będą aktualizowane na podstawie histogramu)
            initialParams.assign(9, 0.0);
            
            // Domyślne granice (będą aktualizowane na podstawie histogramu)
            lowerBounds.assign(9, 0.0);
            upperBounds.assign(9, 1000.0);
            
            // Inicjalizuj wyniki
            fittedParams.assign(9, 0.0);
            paramErrors.assign(9, 0.0);
            chi2 = 0.0;
            ndf = 0;
            status = -1;
            converged = false;
            
            combinedMean = 0.0;
            combinedMeanErr = 0.0;
            combinedStdDev = 0.0;
            combinedStdDevErr = 0.0;
            
            fitRangeMin = -999;
            fitRangeMax = -999;
            fitFunction = nullptr;
        }
        
        /**
         * @brief Ustawia wartości początkowe parametrów
         * @param params Wektor 9 parametrów [A1,μ1,σ1, A2,μ2,σ2, A3,μ3,σ3]
         */
        void SetInitialParams(const std::vector<Double_t>& params) {
            if(params.size() == 9) {
                initialParams = params;
            }
        }
        
        /**
         * @brief Ustawia granice parametrów
         * @param lower Dolne granice [9]
         * @param upper Górne granice [9]
         */
        void SetParamBounds(const std::vector<Double_t>& lower, const std::vector<Double_t>& upper) {
            if(lower.size() == 9 && upper.size() == 9) {
                lowerBounds = lower;
                upperBounds = upper;
            }
        }
        
        /**
         * @brief Ustawia zakres fitu
         * @param xmin Dolna granica
         * @param xmax Górna granica
         */
        void SetFitRange(Double_t xmin, Double_t xmax) {
            fitRangeMin = xmin;
            fitRangeMax = xmax;
        }
    };

    enum class ImageFormat { PNG, SVG, PDF };

    /**
     * @brief Konstruktor
     * @param channNum Liczba kanałów MC (mctruth 1-7)
     * @param channColors Tablica kolorów dla kanałów MC
     * @param channelNames Nazwy kanałów (opcjonalne)
     * @param dataStyle Styl markera dla danych eksperymentalnych
     * @param dataColor Kolor dla danych eksperymentalnych
     * @param sumColor Kolor dla sumy MC
     */
    HistManager(Int_t channNum, const Color_t* channColors, 
                const std::vector<TString>& channelNames = std::vector<TString>(),
                Int_t dataStyle = kFullCircle, Int_t dataColor = kBlack, 
                Float_t dataSize = 1.0, Int_t sumColor = kOrange);

    /// Destruktor
    ~HistManager();

    // Tworzenie zestawów histogramów
    void CreateHistSet1D(const TString& setName, const HistConfig& config);
    void CreateHistSet2D(const TString& setName, const Hist2DConfig& config);

    // Wypełnianie histogramów MC
    void Fill1D(const TString& setName, Int_t mctruth, Double_t value, Double_t weight = 1.0);
    void Fill2D(const TString& setName, Int_t mctruth, Double_t x, Double_t y, Double_t weight = 1.0);
    
    // Wypełnianie histogramów z danymi eksperymentalnymi
    void FillData1D(const TString& setName, Double_t value, Double_t weight = 1.0);
    void FillData2D(const TString& setName, Double_t x, Double_t y, Double_t weight = 1.0);

    // ===== METODY ARRAY HISTOGRAMÓW =====
    
    /**
     * @brief Tworzy zestaw histogramów 1D dla zmiennych typu array (np. pulls[10], momentum[4])
     * @param config Konfiguracja zestawu array
     * @return true jeśli sukces
     */
    Bool_t CreateHistArray1D(const ArrayConfig& config);
    
    /**
     * @brief Wypełnia histogram z zestawu array dla MC
     * @param baseName Bazowa nazwa zestawu
     * @param index Indeks histogramu w zestawie (0-based)
     * @param mctruth Numer kanału MC (1-N)
     * @param value Wartość do wypełnienia
     * @param weight Waga zdarzenia (domyślnie 1.0)
     * @return true jeśli sukces
     */
    Bool_t FillArray1D(const TString& baseName, Int_t index, Int_t mctruth, Double_t value, Double_t weight = 1.0);
    
    /**
     * @brief Wypełnia histogram z zestawu array dla danych eksperymentalnych
     * @param baseName Bazowa nazwa zestawu
     * @param index Indeks histogramu w zestawie (0-based)
     * @param value Wartość do wypełnienia
     * @param weight Waga zdarzenia (domyślnie 1.0)
     * @return true jeśli sukces
     */
    Bool_t FillArrayData1D(const TString& baseName, Int_t index, Double_t value, Double_t weight = 1.0);
    
    /**
     * @brief Wypełnia wszystkie histogramy w zestawie array jednym wywołaniem (MC)
     * @param baseName Bazowa nazwa zestawu
     * @param mctruth Numer kanału MC (1-N)
     * @param values Wektor wartości do wypełnienia (rozmiar musi odpowiadać arraySize)
     * @param weight Waga zdarzenia (domyślnie 1.0)
     * @return true jeśli sukces
     */
    Bool_t FillArrayAll1D(const TString& baseName, Int_t mctruth, const std::vector<Double_t>& values, Double_t weight = 1.0);
    
    /**
     * @brief Wypełnia wszystkie histogramy w zestawie array jednym wywołaniem (dane)
     * @param baseName Bazowa nazwa zestawu
     * @param values Wektor wartości do wypełnienia (rozmiar musi odpowiadać arraySize)
     * @param weight Waga zdarzenia (domyślnie 1.0)
     * @return true jeśli sukces
     */
    Bool_t FillArrayAllData1D(const TString& baseName, const std::vector<Double_t>& values, Double_t weight = 1.0);
    
    /**
     * @brief Rysuje wszystkie histogramy z zestawu array na jednej kanwie
     * @param baseName Bazowa nazwa zestawu
     * @param drawData Czy rysować dane eksperymentalne
     * @param canvasName Nazwa kanwy (domyślnie auto-generowana)
     * @param nCols Liczba kolumn w podziale kanwy (auto = -1)
     * @param nRows Liczba wierszy w podziale kanwy (auto = -1)
     * @return Wskaźnik na utworzoną kanwę
     */
    TCanvas* DrawArray1D(const TString& baseName, Bool_t drawData = false, const TString& canvasName = "", 
                         Int_t nCols = -1, Int_t nRows = -1);
    
    /**
     * @brief Pobiera histogram z zestawu array
     * @param baseName Bazowa nazwa zestawu
     * @param index Indeks histogramu (0-based)
     * @param mctruth Numer kanału MC (1-N) lub -1 dla danych, 0 dla sumy MC
     * @return Wskaźnik na histogram lub nullptr
     */
    TH1D* GetArrayHist1D(const TString& baseName, Int_t index, Int_t mctruth = -1);
    
    /**
     * @brief Normalizuje histogram kanału w array tak, aby jego całka była równa GetEntries
     * @param baseName Bazowa nazwa zestawu array
     * @param index Indeks histogramu (0-based)
     * @param mctruth Numer kanału MC do znormalizowania (1-N)
     * @return Czynnik skalowania użyty do normalizacji, lub -1 w przypadku błędu
     */
    Double_t ScaleArrayChannelByEntries(const TString& baseName, Int_t index, Int_t mctruth);

    /**
     * @brief Skaluje wszystkie histogramy MC w array tak, aby ich suma miała tę samą całkę co dane
     * @param baseName Bazowa nazwa array histogramów  
     * @param index Indeks histogramu w array (0-based)
     * @return Czynnik skalowania użyty do normalizacji, lub -1 w przypadku błędu
     */
    Double_t PerformArraySimpleScaling(const TString& baseName, Int_t index);

    // Rysowanie histogramów
    void DrawSet1D(const TString& setName, const TString& drawOpt = "", Bool_t drawData = false);
    void DrawSet2D(const TString& setName, const TString& drawOpt = "COLZ", Bool_t drawData = false);
    
    /**
     * @brief Rysuje histogram 1D z opcjonalnym fitem Triple Gaussian
     * @param setName Nazwa zestawu histogramów
     * @param mctruth Kanał do narysowania (0=suma MC, 1-N=kanały, -1=dane)
     * @param drawOpt Opcje rysowania histogramu
     * @param performFit Czy wykonać fit Triple Gaussian
     * @param fitParams Parametry fitu (jeśli nullptr, użyje domyślnych)
     * @param drawData Czy dodatkowo narysować dane
     * @return Wskaźnik na canvas lub nullptr
     */
    TCanvas* DrawHistogram1DWithFit(const TString& setName, Int_t mctruth, 
                                   const TString& drawOpt = "HIST",
                                   Bool_t performFit = true,
                                   FitParams3Gauss* fitParams = nullptr,
                                   Bool_t drawData = false);
    
    /**
     * @brief Rysuje histogram z array z opcjonalnym fitem Triple Gaussian
     * @param baseName Bazowa nazwa zestawu array
     * @param index Indeks histogramu w zestawie
     * @param mctruth Kanał do narysowania (0=suma MC, 1-N=kanały, -1=dane)
     * @param drawOpt Opcje rysowania histogramu
     * @param performFit Czy wykonać fit Triple Gaussian
     * @param fitParams Parametry fitu (jeśli nullptr, użyje domyślnych)
     * @return Wskaźnik na canvas lub nullptr
     */
    TCanvas* DrawArrayHistogramWithFit(const TString& baseName, Int_t index, Int_t mctruth,
                                     const TString& drawOpt = "HIST", 
                                     Bool_t performFit = true,
                                     FitParams3Gauss* fitParams = nullptr);

    // Skalowanie histogramów
    /**
     * @brief Normalizuje histogram danego kanału MC tak, aby jego całka była równa GetEntries
     * @param setName Nazwa zestawu histogramów
     * @param mctruth Numer kanału MC do znormalizowania (1-N)
     * @return Czynnik skalowania użyty do normalizacji, lub -1 w przypadku błędu
     */
    Double_t ScaleChannelByEntries(const TString& setName, Int_t mctruth);
    
    /**
     * @brief Normalizuje histogram 2D danego kanału MC tak, aby jego całka była równa GetEntries
     * @param setName Nazwa zestawu histogramów 2D
     * @param mctruth Numer kanału MC do znormalizowania (1-N)
     * @return Czynnik skalowania użyty do normalizacji, lub -1 w przypadku błędu
     */
    Double_t ScaleChannel2DByEntries(const TString& setName, Int_t mctruth);

    // ==================== TRIPLE GAUSSIAN FITTING ====================
    
    /**
     * @brief Przygotowuje domyślne parametry dla fitu Triple Gaussian na podstawie histogramu
     * @param hist Histogram do analizy
     * @param params Struktura parametrów do wypełnienia
     * @return true jeśli sukces
     */
    Bool_t PrepareDefaultTripleGaussParams(TH1* hist, FitParams3Gauss& params);
    
    /**
     * @brief Wykonuje fit Triple Gaussian dla histogramu z zestawu 1D
     * @param setName Nazwa zestawu histogramów
     * @param mctruth Numer kanału MC (1-N), 0 dla sumy MC, -1 dla danych
     * @param params Parametry fitu
     * @return true jeśli fit się udał
     */
    Bool_t FitTripleGauss1D(const TString& setName, Int_t mctruth, FitParams3Gauss& params);
    
    /**
     * @brief Wykonuje fit Triple Gaussian dla histogramu z array
     * @param baseName Bazowa nazwa zestawu array
     * @param index Indeks histogramu w zestawie (0-based)
     * @param mctruth Numer kanału MC (1-N), 0 dla sumy MC, -1 dla danych
     * @param params Parametry fitu
     * @return true jeśli fit się udał
     */
    Bool_t FitTripleGaussArray1D(const TString& baseName, Int_t index, Int_t mctruth, FitParams3Gauss& params);
    
    /**
     * @brief Rysuje wyniki fitu Triple Gaussian na histogramie
     * @param hist Histogram na którym rysować
     * @param params Parametry z wynikami fitu
     * @param drawComponents Czy rysować poszczególne komponenty Gaussa
     */
    void DrawTripleGaussFit(TH1* hist, const FitParams3Gauss& params, Bool_t drawComponents = false);
    
    /**
     * @brief Rysuje poszczególne komponenty Triple Gaussian (pomocnicza dla użytkowników)
     * @param hist Histogram na którym rysować komponenty
     * @param params Parametry fitu z wynikami
     */
    void DrawTripleGaussComponents(TH1* hist, const FitParams3Gauss& params);
    
    /**
     * @brief Wyświetla wyniki fitu Triple Gaussian na terminalu i na canvie
     * @param params Parametry z wynikami fitu
     * @param hist Histogram dla kontekstu wyświetlania
     * @param addToCanvas Czy dodać tekst do aktualnego canvasu
     */
    void DisplayTripleGaussResults(const FitParams3Gauss& params, const TH1* hist, Bool_t addToCanvas = true);
    
    /**
     * @brief Pobiera histogram do fitowania Triple Gaussian
     * @param setName Nazwa zestawu lub baseName dla array
     * @param mctruth Numer kanału MC lub -1 dla danych
     * @param arrayIndex Indeks array (-1 jeśli nie array)
     * @return Wskaźnik na histogram lub nullptr
     */
    TH1* GetHistogramForTripleGaussFit(const TString& setName, Int_t mctruth, Int_t arrayIndex = -1);
    
    /**
     * @brief Pobiera wyniki fitu Triple Gaussian dla zestawu 1D
     * @param setName Nazwa zestawu
     * @param mctruth Numer kanału MC
     * @return Wskaźnik na wyniki fitu lub nullptr jeśli nie znaleziono
     */
    const FitParams3Gauss* GetTripleGaussResults1D(const TString& setName, Int_t mctruth);
    
    /**
     * @brief Pobiera wyniki fitu Triple Gaussian dla array
     * @param baseName Bazowa nazwa zestawu array
     * @param index Indeks histogramu
     * @param mctruth Numer kanału MC
     * @return Wskaźnik na wyniki fitu lub nullptr jeśli nie znaleziono
     */
    const FitParams3Gauss* GetTripleGaussResultsArray(const TString& baseName, Int_t index, Int_t mctruth);
    
    /**
     * @brief Sprawdza czy istnieją wyniki fitu Triple Gaussian
     * @param setName Nazwa zestawu
     * @param mctruth Numer kanału MC
     * @param arrayIndex Indeks array (-1 jeśli nie array)
     * @return true jeśli wyniki istnieją
     */
    Bool_t HasTripleGaussResults(const TString& setName, Int_t mctruth, Int_t arrayIndex = -1);

    // Zapisywanie
    void SaveSet(const TString& setName, const TString& filePattern);
    void SaveToRoot(const TString& filename);
    void SaveSetToRoot(const TString& setName, const TString& filename);
    void ExportSet(const TString& setName, const TString& filePattern, ImageFormat format);

    // ==================== FRACTIONFIT FUNCTIONALITY ====================
    
    /**
     * @brief Typ normalizacji MC do danych
     */
    enum class NormalizationType {
        NONE,           ///< Brak normalizacji
        SIMPLE_SCALE,   ///< Proste skalowanie na podstawie całkowitej liczby zdarzeń
        FRACTION_FIT    ///< Fraction fit z optymalizacją składowych
    };
    
    /**
     * @brief Włącza/wyłącza automatyczne wykonywanie FractionFit podczas rysowania
     * @param use Czy używać FractionFit
     */
    void SetUseFractionFitter(Bool_t use = true) { 
        fUseFractionFitter = use; 
        if(use) fNormalizationType = NormalizationType::FRACTION_FIT;
        else fNormalizationType = NormalizationType::NONE;
    }
    
    /**
     * @brief Ustawia typ normalizacji MC do danych
     * @param type Typ normalizacji
     */
    void SetNormalizationType(NormalizationType type) { 
        fNormalizationType = type; 
        fUseFractionFitter = (type == NormalizationType::FRACTION_FIT);
    }
    
    /**
     * @brief Pobiera aktualny typ normalizacji
     * @return Typ normalizacji
     */
    NormalizationType GetNormalizationType() const { return fNormalizationType; }
    
    /**
     * @brief Ustawia ograniczenia dla FractionFit
     * @param setName Nazwa zestawu histogramów
     * @param constraints Ograniczenia fitu
     */
    void SetFitConstraints(const TString& setName, const FitConstraints& constraints);
    
    /**
     * @brief Wykonuje FractionFit dla danego zestawu histogramów
     * @param setName Nazwa zestawu histogramów
     * @param useStoredConstraints Czy używać wcześniej zapisanych ograniczeń
     * @return Wyniki fitu
     */
    FitResult PerformFractionFit(const TString& setName, Bool_t useStoredConstraints = true);
    
    /**
     * @brief Wykonuje proste skalowanie MC na podstawie całkowitej liczby zdarzeń
     * @param setName Nazwa zestawu histogramów
     * @return Czynnik skalujący zastosowany do MC
     */
    Double_t PerformSimpleScaling(const TString& setName);
    
    /**
     * @brief Pobiera wyniki ostatniego fitu
     * @return Referencja do wyników ostatniego fitu
     */
    const FitResult& GetLastFitResult() const { return fLastFitResult; }
    
    /**
     * @brief Pobiera wyniki fitu dla konkretnego zestawu
     * @param setName Nazwa zestawu
     * @return Wyniki fitu (może być pusty jeśli fit nie był wykonywany)
     */
    FitResult GetFitResult(const TString& setName) const;
    
    /**
     * @brief Sprawdza czy fit został wykonany dla danego zestawu
     * @param setName Nazwa zestawu
     * @return true jeśli fit był wykonywany
     */
    Bool_t HasFitResult(const TString& setName) const;
    
    // Metody kontroli rozmiaru markerów danych
    /**
     * @brief Ustawia rozmiar markera dla danych eksperymentalnych
     * @param size Rozmiar markera (wartość Float_t)
     */
    void SetDataMarkerSize(Float_t size);
    
    /**
     * @brief Pobiera aktualny rozmiar markera dla danych
     * @return Rozmiar markera danych
     */
    Float_t GetDataMarkerSize() const { return fDataSize; }
    
    /**
     * @brief Aktualizuje rozmiar markerów w istniejących histogramach danych
     * Przydatne do natychmiastowego zastosowania nowego rozmiaru bez ponownego tworzenia histogramów
     */
    void UpdateExistingDataHistograms();

private:
    Int_t fChannNum;                                    ///< Liczba kanałów MC
    std::vector<Int_t> fChannColors;                    ///< Kolory kanałów MC
    std::vector<TString> fChannelNames;                 ///< Nazwy kanałów MC
    Int_t fDataStyle;                                   ///< Styl markera dla danych
    Int_t fDataColor;                                   ///< Kolor danych
    Float_t fDataSize;                                  ///< Rozmiar markera dla danych
    Int_t fSumColor;                                    ///< Kolor sumy MC
    Bool_t fUseFractionFitter{false};                   ///< Czy używać FractionFit
    NormalizationType fNormalizationType{NormalizationType::NONE}; ///< Typ normalizacji
    
    // Kontenery danych
    std::map<TString, std::vector<TH1*>> fHists1D;      ///< Histogramy 1D MC
    std::map<TString, std::vector<TH2*>> fHists2D;      ///< Histogramy 2D MC
    std::map<TString, TH1*> fData1D;                    ///< Histogramy 1D danych
    std::map<TString, TH2*> fData2D;                    ///< Histogramy 2D danych
    std::map<TString, std::vector<TCanvas*>> fCanvases; ///< Canvasy
    std::map<TString, HistConfig> fConfigs1D;           ///< Konfiguracje 1D
    std::map<TString, Hist2DConfig> fConfigs2D;         ///< Konfiguracje 2D
    
    // Array histogramy
    std::map<TString, std::vector<std::vector<TH1D*>>> fArrayHists1D; ///< Array histogramy 1D [baseName][index][mctruth+1]
    std::map<TString, std::vector<TH1D*>> fArrayData1D;               ///< Array histogramy danych [baseName][index]
    std::map<TString, ArrayConfig> fArrayConfigs;                     ///< Konfiguracje array histogramów
    
    // FractionFit related
    std::map<TString, FitConstraints> fFitConstraints;  ///< Ograniczenia fitu dla każdego zestawu
    std::map<TString, FitResult> fFitResults;           ///< Wyniki fitów dla każdego zestawu
    FitResult fLastFitResult;                           ///< Wyniki ostatniego fitu
    
    // Triple Gaussian fitting related
    std::map<TString, FitParams3Gauss> f3GaussFitResults;    ///< Wyniki fitów Triple Gaussian dla zestawów 1D [setName_mctruth]
    std::map<TString, FitParams3Gauss> f3GaussArrayResults;  ///< Wyniki fitów Triple Gaussian dla array [baseName_index_mctruth]

    // Metody pomocnicze
    void ConfigureHistogram(TH1* hist, Int_t color, Bool_t showStats);
    Double_t CalculateYRange(const std::vector<TH1*>& hists, const TString& setName, Bool_t logy);
    void CleanupSet(const TString& setName);
    TLegend* CreateLegend(const std::vector<TH1*>& hists, const TH1* dataHist = nullptr);
    
    /**
     * @brief Oblicza optymalne rozmiary kanwy dla N paneli
     * @param nPanels Liczba paneli
     * @param nCols Kolumny (wejście/wyjście)
     * @param nRows Wiersze (wejście/wyjście)  
     */
    void CalculateCanvasLayout(Int_t nPanels, Int_t& nCols, Int_t& nRows);
    
    /**
     * @brief Czyści zestaw array histogramów
     * @param baseName Bazowa nazwa zestawu
     */
    void CleanupArraySet(const TString& baseName);
    
    /**
     * @brief Wykonuje FractionFit z daną konfiguracją
     * @param mcHists Wektor histogramów MC (bez sumy)
     * @param dataHist Histogram z danymi
     * @param constraints Ograniczenia fitu
     * @return Wyniki fitu
     */
    FitResult DoFractionFit(const std::vector<TH1*>& mcHists, TH1* dataHist, 
                           const FitConstraints& constraints);
    
    /**
     * @brief Aktualizuje histogramy MC na podstawie wyników fitu
     * @param mcHists Wektor histogramów MC do aktualizacji
     * @param sumHist Histogram sumy MC do aktualizacji
     * @param fitResult Wyniki fitu
     */
    void UpdateMCHistograms(const std::vector<TH1*>& mcHists, TH1* sumHist, 
                           const FitResult& fitResult);
    
    /**
     * @brief Wykonuje proste skalowanie wszystkich histogramów MC
     * @param mcHists Wektor histogramów MC do przeskalowania
     * @param sumHist Histogram sumy MC do przeskalowania
     * @param scaleFactor Czynnik skalujący
     */
    void ApplySimpleScaling(const std::vector<TH1*>& mcHists, TH1* sumHist, Double_t scaleFactor);
    
    /**
     * @brief Sprawdza stabilność fitu i próbuje poprawić jeśli to konieczne
     * @param fitter Obiekt TFractionFitter
     * @param constraints Ograniczenia fitu
     * @param maxRetries Maksymalna liczba prób
     * @return Status fitu
     */
    Int_t EnsureFitStability(TFractionFitter* fitter, const FitConstraints& constraints, 
                            Int_t maxRetries = 3);
    
    // ==================== TRIPLE GAUSSIAN PRIVATE HELPERS ====================
    
    /**
     * @brief Wykonuje właściwy fit Triple Gaussian na histogramie
     * @param hist Histogram do fitowania
     * @param params Parametry fitu (wejście/wyjście)
     * @return true jeśli fit się udał
     */
    Bool_t DoTripleGaussFit(TH1* hist, FitParams3Gauss& params);
    
    /**
     * @brief Oblicza wartości kombinowane z wyników fitu
     * @param params Parametry z wynikami fitu (do aktualizacji)
     */
    void CalculateCombinedValues(FitParams3Gauss& params);
    
    /**
     * @brief Generuje klucz do przechowywania wyników fitu
     * @param setName Nazwa zestawu
     * @param mctruth Numer kanału MC
     * @param arrayIndex Indeks array (-1 jeśli nie array)
     * @return Klucz string
     */
    TString GenerateTripleGaussKey(const TString& setName, Int_t mctruth, Int_t arrayIndex = -1);
};
