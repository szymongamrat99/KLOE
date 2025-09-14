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
        
        ArrayConfig() : arraySize(0) {}
        
        ArrayConfig(const TString& name, const TString& title, Int_t size, const HistConfig& config) 
            : baseName(name), baseTitle(title), arraySize(size), commonConfig(config) {}
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
                Int_t sumColor = kOrange);

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

    // Rysowanie histogramów
    void DrawSet1D(const TString& setName, const TString& drawOpt = "", Bool_t drawData = false);
    void DrawSet2D(const TString& setName, const TString& drawOpt = "COLZ", Bool_t drawData = false);

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

private:
    Int_t fChannNum;                                    ///< Liczba kanałów MC
    std::vector<Int_t> fChannColors;                    ///< Kolory kanałów MC
    std::vector<TString> fChannelNames;                 ///< Nazwy kanałów MC
    Int_t fDataStyle;                                   ///< Styl markera dla danych
    Int_t fDataColor;                                   ///< Kolor danych
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
};
