#!/bin/bash

rootfilename=$1
foldername=$2

if [[ -z "$rootfilename" || -z "$foldername" ]]; then
    echo "Usage: $0 <root_file_name> <folder_name>"
    exit 1
fi

# Definiuje opcje menu (bez EXIT jako scenariusza fizycznego)
options=("NO_CUTS" "SHORTER_KAON_PATHS" "OLD_CHI2_CUT" "OLD_TRCSUM_CUT" "OLD_COMBINED_MASS_PI0_CUT" "OLD_MASS_KCH_CUT" "OLD_MASS_KNE_CUT" "OLD_QMISS_CUT" "OLD_OPENING_ANGLE_CUT" "OLD_OMEGA_GEOMETRICAL_CUT" "OLD_OMEGA_FIDUCIAL_VOLUME" "SIMONA_CHI2_CUT" "BAD_CLUS_SIMONA" "SIMONA_KIN_CUTS" "SIMONA_ALL_CUTS" "OMEGA_MASS_T0_CUT" "BLOB" "NO_BLOB")

echo "Dostępne scenariusze:"
for i in "${!options[@]}"; do
    printf "%2d) %s\n" "$((i + 1))" "${options[$i]}"
done
echo ""
echo "Podaj wiele pozycji naraz:"
echo "- numery: np. 1 4 8"
echo "- lub nazwy: np. NO_CUTS SIMONA_CHI2_CUT"
echo "- możesz mieszać i rozdzielać spacją/przecinkiem"
echo "- wpisz ALL aby wybrać wszystko, EXIT aby wyjść"

read -r -p "Twój wybór: " selection_raw

if [[ -z "$selection_raw" ]]; then
    echo "Brak wyboru. Kończę."
    exit 1
fi

selection_upper=$(echo "$selection_raw" | tr '[:lower:]' '[:upper:]')
if [[ "$selection_upper" == "EXIT" ]]; then
    echo "Exiting..."
    exit 0
fi

selected=()

append_unique() {
    local candidate="$1"
    local existing
    for existing in "${selected[@]}"; do
        if [[ "$existing" == "$candidate" ]]; then
            return
        fi
    done
    selected+=("$candidate")
}

if [[ "$selection_upper" == "ALL" ]]; then
    selected=("${options[@]}")
else
    normalized=$(echo "$selection_raw" | tr ',' ' ')
    for token in $normalized; do
        token_upper=$(echo "$token" | tr '[:lower:]' '[:upper:]')

        if [[ "$token_upper" == "EXIT" ]]; then
            echo "Exiting..."
            exit 0
        fi

        if [[ "$token" =~ ^[0-9]+$ ]]; then
            idx=$((token - 1))
            if (( idx >= 0 && idx < ${#options[@]} )); then
                append_unique "${options[$idx]}"
            else
                echo "Pomijam niepoprawny numer: $token"
            fi
            continue
        fi

        matched=0
        for opt_name in "${options[@]}"; do
            if [[ "$token_upper" == "$opt_name" ]]; then
                append_unique "$opt_name"
                matched=1
                break
            fi
        done

        if [[ $matched -eq 0 ]]; then
            echo "Pomijam nieznaną opcję: $token"
        fi
    done
fi

if (( ${#selected[@]} == 0 )); then
    echo "Nie wybrano żadnego poprawnego scenariusza."
    exit 1
fi

opt=$(IFS=';'; echo "${selected[*]}")
echo "Chosen scenarios: $opt"

opt="$opt;ADDNAME=${foldername};ROOTFILE=${rootfilename}"

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 70; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Signal_%d.root",i));\
}
for (Int_t i = 1; i <= 24; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Signal_%d.root",i));\
}
for (Int_t i = 1; i <= 29; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Signal_%d.root",i));\
}

for (Int_t i = 1; i <= 56; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Omega_%d.root",i));\
}
for (Int_t i = 1; i <= 56; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Omega_%d.root",i));\
}
for (Int_t i = 1; i <= 12; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Omega_%d.root",i));\
}

chain->Process("signal_vs_bcg_v3.C", "$opt");
.q
EOF

