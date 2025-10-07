#pragma once

#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TString.h>
#include <map>
#include <memory>
#include <vector>
#include <iostream>
#include <type_traits>
#include "VariableConfig.h"

namespace KLOE {

/**
 * @class DataAccessWrapper
 * @brief Wrapper umożliwiający jednolity dostęp do plików v1.root i v2.root
 * 
 * Klasa automatycznie rozpoznaje typ pliku na podstawie nazwy i używa odpowiedniej
 * metody dostępu: TTreeReader dla v2.root, SetBranchAddress dla v1.root
 */
class DataAccessWrapper {
public:
    enum class FileVersion {
        V1,      ///< Stare pliki *v1.root - używaj SetBranchAddress
        V2,      ///< Nowe pliki *v2.root - używaj TTreeReader
        UNKNOWN  ///< Nieznany typ pliku
    };

    /**
     * @brief Konstruktor
     * @param chain Referencja do TChain
     */
    explicit DataAccessWrapper(TChain& chain);
    
    /// Destruktor
    ~DataAccessWrapper();

    /**
     * @brief Inicjalizacja wrapper'a - musi być wywołana przed użyciem
     * @return true jeśli sukces
     */
    Bool_t Initialize();
    
    /**
     * @brief Przechodzi do następnego wydarzenia
     * @return true jeśli są jeszcze wydarzenia
     */
    Bool_t Next();
    
    /**
     * @brief Pobiera całkowitą liczbę wydarzeń
     * @return Liczba wydarzeń w chain'ie
     */
    Long64_t GetEntries() const { return fChain.GetEntries(); }
    
    /**
     * @brief Pobiera aktualny numer wydarzenia
     * @return Numer aktualnego wydarzenia
     */
    Long64_t GetCurrentEntry() const { return fCurrentEntry; }

    // ==================== GETTERY DLA ZMIENNYCH SKALARNYCH ====================
    
    Int_t GetNRun() const;
    Int_t GetNEv() const;
    Int_t GetNClu() const;
    Int_t GetNTCl() const;
    Int_t GetNV() const;
    Int_t GetNTV() const;
    Int_t GetNTMC() const;
    Int_t GetNVtxMC() const;
    Int_t GetEclFilfo() const;
    Int_t GetEclFilfoWord() const;
    Int_t GetBunchNum() const;
    Int_t GetNECls() const;
    
    Float_t GetT0Step1() const;
    Float_t GetBx() const;
    Float_t GetBy() const;
    Float_t GetBz() const;
    Float_t GetBxErr() const;
    Float_t GetByErr() const;
    Float_t GetBzErr() const;
    Float_t GetBpx() const;
    Float_t GetBpy() const;
    Float_t GetBpz() const;
    Float_t GetBpxErr() const;
    Float_t GetBpyErr() const;
    Float_t GetBpzErr() const;
    Float_t GetBRoots() const;
    Float_t GetBRootsErr() const;

    // ==================== GETTERY DLA TABLIC ====================
    
    const std::vector<Int_t>& GetEclStream() const;
    const std::vector<Int_t>& GetAssCl() const;
    const std::vector<Int_t>& GetIv() const;
    const std::vector<Int_t>& GetVtxMC() const;
    const std::vector<Int_t>& GetPidMC() const;
    const std::vector<Int_t>& GetMother() const;
    
    const std::vector<Float_t>& GetXCl() const;
    const std::vector<Float_t>& GetYCl() const;
    const std::vector<Float_t>& GetZCl() const;
    const std::vector<Float_t>& GetTCl() const;
    const std::vector<Float_t>& GetEneCl() const;
    const std::vector<Float_t>& GetCurv() const;
    const std::vector<Float_t>& GetPhiv() const;
    const std::vector<Float_t>& GetCotv() const;
    const std::vector<Float_t>& GetXv() const;
    const std::vector<Float_t>& GetYv() const;
    const std::vector<Float_t>& GetZv() const;
    const std::vector<Float_t>& GetXvMC() const;
    const std::vector<Float_t>& GetYvMC() const;
    const std::vector<Float_t>& GetZvMC() const;
    const std::vector<Float_t>& GetPxMC() const;
    const std::vector<Float_t>& GetPyMC() const;
    const std::vector<Float_t>& GetPzMC() const;

    // ==================== UTILITY METHODS ====================
    
    /**
     * @brief Sprawdza aktualną wersję pliku
     * @return Wersja aktualnie przetwarzanego pliku
     */
    FileVersion GetCurrentFileVersion() const { return fCurrentFileVersion; }
    
    /**
     * @brief Wypisuje statystyki typów plików
     */
    void PrintFileTypeStats() const;
    
    /**
     * @brief Wypisuje konfigurację zmiennych
     */
    void PrintVariableConfiguration() const { fVariableConfig.PrintConfiguration(); }

    /**
     * @brief Pobiera wartość zmiennej skalarnej (metoda publiczna dla template specializations)
     * @param key Klucz zmiennej w konfiguracji
     * @return Wartość zmiennej
     */
    template<typename T>
    T GetScalarValue(const TString& key) const;
    
    /**
     * @brief Pobiera wartość tablicy (metoda publiczna dla template specializations)
     * @param key Klucz zmiennej w konfiguracji
     * @return Referencja do wektora z danymi
     */
    template<typename T>
    const std::vector<T>& GetArrayValue(const TString& key) const;
    
    /**
     * @brief Konwertuje C array do std::vector dla plików v1 (metoda publiczna dla template specializations)
     * @param key Klucz zmiennej
     * @param size Rozmiar tablicy
     * @return Referencja do wektora
     */
    template<typename T>
    const std::vector<T>& ConvertArrayToVector(const TString& key, Int_t size) const;

    // Maksymalny rozmiar tablic dla v1 (musi być publiczny dla friend classes)
    static const Int_t kMaxArraySize = 500;

    // Pola publiczne dla template specializations
    TChain& fChain;                          ///< Referencja do chain'a
    VariableConfig fVariableConfig;          ///< Konfiguracja mapowania zmiennych
    Long64_t fCurrentEntry = -1;             ///< Aktualny numer wydarzenia
    FileVersion fCurrentFileVersion = FileVersion::UNKNOWN; ///< Wersja aktualnego pliku
    
    // Mapa nazw plików do wersji
    std::map<TString, FileVersion> fFileVersionMap;
    
    // Statystyki
    Long64_t fV1Events = 0;                  ///< Liczba wydarzeń z plików v1
    Long64_t fV2Events = 0;                  ///< Liczba wydarzeń z plików v2

    // ==================== TREEREADER (dla v2) ====================
    std::unique_ptr<TTreeReader> fReader;
    
    // Mapy dla TTreeReaderValue i TTreeReaderArray
    std::map<TString, std::unique_ptr<TTreeReaderValue<Int_t>>> fIntValues_v2;
    std::map<TString, std::unique_ptr<TTreeReaderValue<Float_t>>> fFloatValues_v2;
    std::map<TString, std::unique_ptr<TTreeReaderValue<UInt_t>>> fUIntValues_v2;
    
    std::map<TString, std::unique_ptr<TTreeReaderArray<Int_t>>> fIntArrays_v2;
    std::map<TString, std::unique_ptr<TTreeReaderArray<Float_t>>> fFloatArrays_v2;
    std::map<TString, std::unique_ptr<TTreeReaderArray<UInt_t>>> fUIntArrays_v2;

    // ==================== SETBRANCHADDRESS (dla v1) ====================
    
    // Mapy dla zmiennych skalarnych v1
    std::map<TString, Int_t> fIntValues_v1;
    std::map<TString, Float_t> fFloatValues_v1;
    std::map<TString, UInt_t> fUIntValues_v1;
    
    // Mapy dla tablic v1 (fixed size arrays)
    std::map<TString, std::vector<Int_t>> fIntArrays_v1;
    std::map<TString, std::vector<Float_t>> fFloatArrays_v1;
    std::map<TString, std::vector<UInt_t>> fUIntArrays_v1;
    
    // Bufory dla konwersji z C arrays do std::vector
    mutable std::map<TString, std::vector<Int_t>> fIntVectorCache;
    mutable std::map<TString, std::vector<Float_t>> fFloatVectorCache;
    mutable std::map<TString, std::vector<UInt_t>> fUIntVectorCache;

private:    // ==================== METODY POMOCNICZE ====================
    
    /**
     * @brief Rozpoznaje wersję pliku na podstawie nazwy
     * @param filename Nazwa pliku
     * @return Wersja pliku
     */
    FileVersion DetermineFileVersion(const TString& filename);
    
    /**
     * @brief Mapuje wszystkie pliki w chain'ie do wersji
     */
    void MapFileVersions();
    
    /**
     * @brief Inicjalizuje TTreeReader dla plików v2
     * @return true jeśli sukces
     */
    Bool_t InitializeTreeReader();
    
    /**
     * @brief Inicjalizuje SetBranchAddress dla plików v1
     * @return true jeśli sukces
     */
    Bool_t InitializeBranchAddress();
    
    /**
     * @brief Pobiera aktualną nazwę pliku z chain'a
     * @return Nazwa aktualnego pliku
     */
    TString GetCurrentFileName() const;
    
    /**
     * @brief Aktualizuje wersję aktualnego pliku
     */
    void UpdateCurrentFileVersion();
};

} // namespace KLOE