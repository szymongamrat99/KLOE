#pragma once

#include <TString.h>
#include <map>
#include <vector>

namespace KLOE {

/**
 * @class VariableConfig
 * @brief Konfiguracja mapowania nazw zmiennych między plikami v1 i v2
 * 
 * Klasa pozwala na mapowanie nazw zmiennych, które różnią się między
 * starymi (v1) a nowymi (v2) plikami ROOT.
 */
class VariableConfig {
public:
    /**
     * @brief Struktura opisująca zmienną w drzewie
     */
    struct VariableInfo {
        TString nameV1;         ///< Nazwa zmiennej w plikach v1
        TString nameV2;         ///< Nazwa zmiennej w plikach v2
        TString type;           ///< Typ zmiennej (Int_t, Float_t, UInt_t)
        Bool_t isArray;         ///< Czy zmienna jest tablicą
        TString sizeVariable;   ///< Nazwa zmiennej określającej rozmiar tablicy
        Bool_t isMC;          ///< Czy zmienna dotyczy danych MC
        
        VariableInfo() : isArray(false), isMC(false) {}
        
        VariableInfo(const TString& v1, const TString& v2, const TString& t, 
                    Bool_t array = false, const TString& sizeVar = "", Bool_t mc = false) 
            : nameV1(v1), nameV2(v2), type(t), isArray(array), sizeVariable(sizeVar), isMC(mc) {}

        VariableInfo(const TString& v1, const TString& v2, const TString& t, 
                    Bool_t array = false, const TString& sizeVar = "", Bool_t mc = false) 
            : nameV1(v1), nameV2(v2), type(t), isArray(array), sizeVariable(sizeVar), isMC(mc) {}
    };

private:
    std::map<TString, VariableInfo> fVariableMap;
    
public:
    VariableConfig();
    
    /**
     * @brief Dodaj mapowanie zmiennej
     */
    void AddVariable(const TString& key, const VariableInfo& info);
    
    /**
     * @brief Pobierz informacje o zmiennej
     */
    const VariableInfo* GetVariableInfo(const TString& key) const;
    
    /**
     * @brief Pobierz nazwę zmiennej dla danej wersji pliku
     */
    TString GetVariableName(const TString& key, Bool_t isV2) const;
    
    /**
     * @brief Sprawdź czy zmienna istnieje w konfiguracji
     */
    Bool_t HasVariable(const TString& key) const;
    
    /**
     * @brief Pobierz listę wszystkich kluczy zmiennych
     */
    std::vector<TString> GetAllKeys() const;
    
    /**
     * @brief Wypisz wszystkie mapowania (debug)
     */
    void PrintConfiguration() const;

    std::vector<TString> GetMCBranches() const;
    
private:
    /**
     * @brief Inicjalizuj domyślne mapowania zmiennych
     */
    void InitializeDefaultMappings();
};

} // namespace KLOE