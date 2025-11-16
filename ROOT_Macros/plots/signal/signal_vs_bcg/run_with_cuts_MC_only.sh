#!/bin/bash

CUTS_FILE="${1:-cut-limits.json}"
HYPOTHESIS="${2:-SIMONA_ANALYSIS}"
PREDEFINED_CUTS="${3:-}"

RESULT_FILE="/tmp/selected_cuts_result.txt"

if [[ -n "$PREDEFINED_CUTS" ]]; then
    CUTS_STRING="$PREDEFINED_CUTS"
    echo "Użycie predefiniowanych cięć: $CUTS_STRING"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../" && pwd)/menu"
    
    # Uruchom interaktywnie
    bash "$SCRIPT_DIR/select_cuts.sh" "$CUTS_FILE" "$HYPOTHESIS"
    
    # Przeczytaj wynik z pliku
    if [[ -f "$RESULT_FILE" ]]; then
        CUTS_STRING=$(cat "$RESULT_FILE")
        rm -f "$RESULT_FILE"
    else
        echo "Błąd: Nie znaleziono pliku wynikowego" >&2
        exit 1
    fi
fi

echo ""
echo "Uruchamianie analizy z opcjami: $CUTS_STRING"
echo ""

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 73; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-10/mk0*all_phys_SIGNAL_MIXED_Signal_%d.root",i));\
}
chain->Process("signal_vs_bcg_v2.C", "$CUTS_STRING");
.q
EOF