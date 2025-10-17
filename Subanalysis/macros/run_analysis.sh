#!/bin/bash

# Ustawienie katalogu bazowego (macros/) jako bieżącego katalogu roboczego
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR"

# --- Konfiguracja Ścieżek ---
INC_PATH="inc" 
LIB_PATH="lib" # Nowy katalog na biblioteki

DEF_FILE="$INC_PATH/HistoGeneral.cxx"
SELECTOR_FILE="plots/signal/signal_vs_bcg/signal_vs_bcg_v2.C"

# --- Tworzenie katalogu wyjściowego, jeśli nie istnieje ---
mkdir -p "$LIB_PATH"

# --- Ustawienie Zmiennej Środowiskowej ROOT-a ---
# Ta zmienna mówi ROOT-owi, gdzie ma zapisywać tymczasowe pliki kompilacji (.so)
export ROOT_ACLIC_CACHE="$LIB_PATH"

# --- Weryfikacja plików (pominięta dla zwięzłości) ---
# ...

# --- Uruchomienie ROOT-a w Trybie Batch ---

root -b << EOF
{
    // 1. Ustawienie ścieżki do nagłówków (gdzie są .h)
    gSystem->AddIncludePath("-I./$INC_PATH");

    // 2. Skompilowanie i załadowanie definicji (Mapy, Kolory, Histogramy)
    // Plik .so zostanie utworzony w katalogu $LIB_PATH dzięki ROOT_ACLIC_CACHE
    gROOT->LoadMacro("$DEF_FILE+"); 

    // 3. Skompilowanie i załadowanie TSelectora
    gROOT->LoadMacro("$SELECTOR_FILE+");
    
    // 4. Uruchomienie analizy
    TChain chain("h1");
    for (Int_t i = 1; i <= 6; i++)\
    {\
      chain.Add(Form("../../../../InitialAnalysis/root_files/2025-10-16/mk0*all_phys_SIGNAL_MIXED_%d.root",i));\
    }\
    
    if (gROOT->GetClass("signal_vs_bcg_v2")) {
        // Uruchomienie TSelector (kluczowy krok)
        chain.Process("signal_vs_bcg_v2.C"); 
        
    } else {
        std::cerr << "BŁĄD: Nie udało się załadować klasy signal_vs_bcg_v2." << std::endl;
    }
}
EOF

# Sprawdzenie kodu wyjścia
if [ $? -eq 0 ]; then
    echo "Analiza zakończona pomyślnie."
else
    echo "BŁĄD: Analiza zakończona niepowodzeniem."
fi