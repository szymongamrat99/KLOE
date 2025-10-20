#pragma once

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TMath.h>
#include <TFitResult.h>
#include <iostream>
#include <vector>

// RooFit includes
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooAbsReal.h>

namespace KLOE {

/**
 * @class TripleGaussFitter
 * @brief Klasa do fitowania potrójnego Gaussa do histogramów
 * 
 * Klasa umożliwia łatwe dopasowanie funkcji triple_gaus do histogramów,
 * wyświetlanie statystyk fitu oraz wizualizację wyników.
 */
class TripleGaussFitter {
public:
    /**
     * @struct FitParams
     * @brief Parametry fitu potrójnego Gaussa
     */
    struct FitParams {
        std::vector<Double_t> initialValues;   ///< Wartości początkowe [A1,μ1,σ1,A2,μ2,σ2,A3,μ3,σ3]
        std::vector<Double_t> lowerBounds;     ///< Dolne granice parametrów
        std::vector<Double_t> upperBounds;     ///< Górne granice parametrów
        Double_t fitRangeMin = -999;           ///< Dolna granica fitu (-999 = auto)
        Double_t fitRangeMax = -999;           ///< Górna granica fitu (-999 = auto)
        TString fitOptions = "RBQ";            ///< Opcje fitu ROOT
        
        FitParams() {
            initialValues.resize(9, 0.0);
            lowerBounds.resize(9, 0.0);
            upperBounds.resize(9, 1000000.0);
        }
        
        void SetFitRange(Double_t xmin, Double_t xmax) {
            fitRangeMin = xmin;
            fitRangeMax = xmax;
        }
        
        void SetInitialValues(const std::vector<Double_t>& values) {
            if(values.size() == 9) {
                initialValues = values;
            }
        }
        
        void SetBounds(const std::vector<Double_t>& lower, const std::vector<Double_t>& upper) {
            if(lower.size() == 9 && upper.size() == 9) {
                lowerBounds = lower;
                upperBounds = upper;
            }
        }
    };
    
    /**
     * @struct FitResult
     * @brief Wyniki fitu potrójnego Gaussa
     */
    struct FitResult {
        Bool_t converged = false;              ///< Czy fit się zbiegł
        Int_t status = -1;                     ///< Status fitu ROOT
        Double_t chi2 = 0.0;                   ///< Chi-kwadrat
        Int_t ndf = 0;                         ///< Liczba stopni swobody
        Double_t probability = 0.0;            ///< Prawdopodobieństwo chi2
        
        // Parametry fitowane
        std::vector<Double_t> parameters;      ///< Parametry fitu [9]
        std::vector<Double_t> errors;          ///< Błędy parametrów [9]
        
        // Parametry kombinowane
        Double_t combinedMean = 0.0;           ///< Średnia kombinowana
        Double_t combinedMeanErr = 0.0;        ///< Błąd średniej kombinowanej
        Double_t combinedSigma = 0.0;          ///< Sigma kombinowana
        Double_t combinedSigmaErr = 0.0;       ///< Błąd sigmy kombinowanej
        
        // Funkcja fitu (dla rysowania)
        TF1* fitFunction = nullptr;            ///< Funkcja fitu (backward compatibility)
        
        // RooFit objects
        RooFitResult* rooFitResult = nullptr;  ///< RooFit result object
        RooAddPdf* rooPdf = nullptr;           ///< Combined PDF
        RooDataHist* rooDataHist = nullptr;    ///< Data histogram
        
        FitResult() {
            parameters.resize(9, 0.0);
            errors.resize(9, 0.0);
        }
        
        ~FitResult() {
            if(fitFunction) {
                delete fitFunction;
                fitFunction = nullptr;
            }
            if(rooFitResult) {
                delete rooFitResult;
                rooFitResult = nullptr;
            }
            if(rooPdf) {
                delete rooPdf;
                rooPdf = nullptr;
            }
            if(rooDataHist) {
                delete rooDataHist;
                rooDataHist = nullptr;
            }
        }
        
        // Copy constructor
        FitResult(const FitResult& other) : 
            converged(other.converged), status(other.status), 
            chi2(other.chi2), ndf(other.ndf), probability(other.probability),
            parameters(other.parameters), errors(other.errors),
            combinedMean(other.combinedMean), combinedMeanErr(other.combinedMeanErr),
            combinedSigma(other.combinedSigma), combinedSigmaErr(other.combinedSigmaErr) {
            
            fitFunction = other.fitFunction ? (TF1*)other.fitFunction->Clone() : nullptr;
        }
        
        // Assignment operator
        FitResult& operator=(const FitResult& other) {
            if(this != &other) {
                converged = other.converged;
                status = other.status;
                chi2 = other.chi2;
                ndf = other.ndf;
                probability = other.probability;
                parameters = other.parameters;
                errors = other.errors;
                combinedMean = other.combinedMean;
                combinedMeanErr = other.combinedMeanErr;
                combinedSigma = other.combinedSigma;
                combinedSigmaErr = other.combinedSigmaErr;
                
                if(fitFunction) {
                    delete fitFunction;
                }
                fitFunction = other.fitFunction ? (TF1*)other.fitFunction->Clone() : nullptr;
            }
            return *this;
        }
    };

public:
    /// Konstruktor
    TripleGaussFitter();
    
    /// Destruktor
    ~TripleGaussFitter();

    // ==================== GŁÓWNE METODY FITOWANIA ====================
    
    /**
     * @brief Fituje potrójnego Gaussa do histogramu z automatycznymi parametrami
     * @param hist Histogram do fitowania
     * @param setTitle Opcjonalny tytuł dla wyników
     * @return true jeśli fit się udał
     */
    Bool_t FitHistogram(TH1* hist, const TString& setTitle = "");
    
    /**
     * @brief Fituje potrójnego Gaussa do histogramu z niestandardowymi parametrami
     * @param hist Histogram do fitowania
     * @param params Parametry fitu
     * @param setTitle Opcjonalny tytuł dla wyników
     * @return true jeśli fit się udał
     */
    Bool_t FitHistogram(TH1* hist, const FitParams& params, const TString& setTitle = "");
    
    /**
     * @brief Przygotowuje domyślne parametry na podstawie histogramu
     * @param hist Histogram do analizy
     * @return Struktura z domyślnymi parametrami
     */
    FitParams PrepareDefaultParams(TH1* hist);
    
    /**
     * @brief Przygotowuje restrykcyjne parametry dla histogramów z długimi ogonami
     * @param hist Histogram do analizy
     * @param meanWindow Okno dla średniej w jednostkach RMS (domyślnie 3.0)
     * @param sigmaMax Maksymalna sigma w jednostkach RMS (domyślnie 3.0)
     * @return Struktura z ograniczonymi parametrami
     */
    FitParams PrepareConstrainedParams(TH1* hist, 
                                      Double_t meanWindow = 3.0,
                                      Double_t sigmaMax = 3.0);

    // ==================== DOSTĘP DO WYNIKÓW ====================
    
    /**
     * @brief Pobiera ostatnie wyniki fitu
     * @return Referencja do wyników
     */
    const FitResult& GetLastResults() const { return fLastResult; }
    
    /**
     * @brief Sprawdza czy ostatni fit się udał
     */
    Bool_t IsLastFitSuccessful() const { return fLastResult.converged; }
    
    /**
     * @brief Pobiera chi2/NDF ostatniego fitu
     */
    Double_t GetChi2NDF() const { 
        return fLastResult.ndf > 0 ? fLastResult.chi2 / fLastResult.ndf : -1.0; 
    }

    // ==================== WYŚWIETLANIE WYNIKÓW ====================
    
    /**
     * @brief Wypisuje wyniki fitu na terminal
     * @param detailed Czy wypisać szczegółowe informacje
     */
    void PrintResults(Bool_t detailed = true) const;
    
    /**
     * @brief Wypisuje skrócone podsumowanie fitu
     */
    void PrintSummary() const;
    
    /**
     * @brief Pobiera tekstowe podsumowanie wyników
     * @return String z podsumowaniem
     */
    TString GetSummaryText() const;

    // ==================== RYSOWANIE ====================
    
    /**
     * @brief Rysuje histogram z dopasowanym fitem
     * @param canvasName Nazwa canvasu (auto jeśli puste)
     * @param canvasTitle Tytuł canvasu (auto jeśli pusty)
     * @param width Szerokość canvasu
     * @param height Wysokość canvasu
     * @param drawComponents Czy rysować poszczególne komponenty Gaussa
     * @param showStats Czy pokazać statystyki na canvasie
     * @param logScale Czy używać skali logarytmicznej (auto jeśli nie podano)
     * @return Wskaźnik na canvas
     */
    TCanvas* DrawFit(const TString& canvasName = "", const TString& canvasTitle = "",
                     Int_t width = 800, Int_t height = 600,
                     Bool_t drawComponents = true, Bool_t showStats = true,
                     Bool_t logScale = false);
    
    /**
     * @brief Rysuje fit na istniejącym canvasie/padzie
     * @param drawComponents Czy rysować komponenty
     * @param showStats Czy dodać statystyki
     */
    void DrawFitOnCurrentPad(Bool_t drawComponents = true, Bool_t showStats = true);
    
    /**
     * @brief Zapisuje canvas z fitem do pliku
     * @param filename Nazwa pliku (z rozszerzeniem)
     * @return true jeśli udało się zapisać
     */
    Bool_t SaveFitToFile(const TString& filename);

    // ==================== USTAWIENIA ====================
    
    /**
     * @brief Ustawia kolory linii dla komponentów
     * @param mainColor Kolor głównej funkcji fitu
     * @param comp1Color Kolor pierwszego komponentu
     * @param comp2Color Kolor drugiego komponentu  
     * @param comp3Color Kolor trzeciego komponentu
     */
    void SetComponentColors(Color_t mainColor = kRed, Color_t comp1Color = kBlue,
                           Color_t comp2Color = kGreen, Color_t comp3Color = kMagenta);
    
    /**
     * @brief Ustawia style linii dla komponentów
     * @param mainStyle Styl głównej funkcji
     * @param componentStyle Styl komponentów
     */
    void SetLineStyles(Style_t mainStyle = 1, Style_t componentStyle = 2);
    
    /**
     * @brief Ustawia szerokość linii
     * @param mainWidth Szerokość głównej funkcji
     * @param componentWidth Szerokość komponentów
     */
    void SetLineWidths(Width_t mainWidth = 2, Width_t componentWidth = 1);
    
    /**
     * @brief Włącza/wyłącza tryb verbose
     */
    void SetVerbose(Bool_t verbose = true) { fVerbose = verbose; }
    
    /**
     * @brief Wybiera metodę fitowania: true = RooFit, false = TF1 (legacy)
     */
    void UseRooFit(Bool_t useRooFit = true) { fUseRooFit = useRooFit; }

private:
    // Wyniki ostatniego fitu
    FitResult fLastResult;
    TH1* fLastHistogram = nullptr;           ///< Ostatnio fitowany histogram
    TString fLastTitle;                      ///< Tytuł ostatniego fitu
    
    // Ustawienia graficzne
    Color_t fMainColor = kRed;
    Color_t fComp1Color = kBlue;
    Color_t fComp2Color = kGreen;
    Color_t fComp3Color = kMagenta;
    Style_t fMainLineStyle = 1;
    Style_t fComponentLineStyle = 2;
    Width_t fMainLineWidth = 2;
    Width_t fComponentLineWidth = 1;
    
    // Ustawienia
    Bool_t fVerbose = false;
    Bool_t fUseRooFit = false;  // Domyślnie używaj starego TF1
    
    // Canvas dla rysowania
    TCanvas* fCanvas = nullptr;

    // ==================== METODY POMOCNICZE ====================
    
    /**
     * @brief Wykonuje właściwy fit
     * @param hist Histogram
     * @param params Parametry fitu
     * @return true jeśli sukces
     */
    Bool_t DoFit(TH1* hist, const FitParams& params);
    
    /**
     * @brief Wykonuje fit używając TF1 (stara metoda)
     */
    Bool_t DoFitTF1(TH1* hist, const FitParams& params);
    
    /**
     * @brief Wykonuje fit używając RooFit (nowa metoda)
     */
    Bool_t DoFitRooFit(TH1* hist, const FitParams& params);
    
    /**
     * @brief Oblicza wartości kombinowane
     */
    void CalculateCombinedValues();
    
    /**
     * @brief Tworzy funkcje komponentów dla rysowania
     * @return Wektor funkcji TF1*
     */
    std::vector<TF1*> CreateComponentFunctions();
    
    /**
     * @brief Tworzy tekst ze statystykami
     * @return TPaveText z wynikami
     */
    TPaveText* CreateStatsBox();
    
    /**
     * @brief Czyści poprzednie wyniki
     */
    void CleanupPreviousResults();
    
    /**
     * @brief Automatycznie ustawia zakres osi Y na podstawie zawartości histogramu
     * @param hist Histogram do analizy
     * @param isLogScale Czy używać skali logarytmicznej
     */
    void SetOptimalYAxisRange(TH1* hist, Bool_t isLogScale = false);
    
    /**
     * @brief Tworzy TF1 dla backward compatibility z wynikami RooFit
     * @param xmin Dolna granica funkcji
     * @param xmax Górna granica funkcji
     */
    void CreateBackwardCompatibleTF1(Double_t xmin, Double_t xmax);
};

} // namespace KLOE