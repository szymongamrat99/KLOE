#!/bin/bash


# Ustawia monit, który wyświetla się po menu
PS3='Wybierz opcję (podaj numer): '

# Definiuje opcje menu
options=("NO_CUTS" "SHORTER_KAON_PATHS" "OLD_CUTS" "SIMONA_CHI2_CUT" "BAD_CLUS_SIMONA" "SIMONA_KIN_CUTS" "SIMONA_ALL_CUTS" "OMEGA_MASS_T0_CUT" "BLOB" "NO_BLOB" "EXIT")
# select wyświetla menu i czeka na wybór
select opt in "${options[@]}"
do
    # case obsługuje wybraną opcję
    case $opt in
        "NO_CUTS")
            echo "Chosen option: NO_CUTS"
            break
            ;;
        "SHORTER_KAON_PATHS")
            echo "Chosen option: SHORTER_KAON_PATHS"
            break
            ;;
        "OLD_CUTS")
            echo "Chosen option: OLD_CUTS"
            break
            ;;
        "SIMONA_CHI2_CUT")
            echo "Chosen option: SIMONA_CHI2_CUT"
            break
            ;;
        "BAD_CLUS_SIMONA")
            echo "Chosen option: BAD_CLUS_SIMONA"
            break
            ;;
        "SIMONA_KIN_CUTS")
            echo "Chosen option: SIMONA_KIN_CUTS"
            break
            ;;
        "SIMONA_ALL_CUTS")
            echo "Chosen option: SIMONA_ALL_CUTS"
            break
            ;;
        "OMEGA_MASS_T0_CUT")
            echo "Chosen option: OMEGA_MASS_T0_CUT"
            break
            ;;
        "BLOB")
            echo "Chosen option: BLOB"
            break
            ;;
        "NO_BLOB")
            echo "Chosen option: NO_BLOB"
            break
            ;;
        "EXIT")
            echo "Exiting..."
            break
            ;;
        *)
            echo "Invalid choice: $REPLY" # $REPLY zawiera wpisany numer
            ;;
    esac
    echo # Dodaje pustą linię dla czytelności
done

if [ "$opt" == "EXIT" ]; then
    exit 0
fi

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 68; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Signal_%d.root",i));\
}
for (Int_t i = 1; i <= 24; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Signal_%d.root",i));\
}
for (Int_t i = 1; i <= 15; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Signal_%d.root",i));\
}

for (Int_t i = 1; i <= 55; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Regeneration_%d.root",i));\
}
for (Int_t i = 1; i <= 54; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Regeneration_%d.root",i));\
}
for (Int_t i = 1; i <= 12; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Regeneration_%d.root",i));\
}

for (Int_t i = 1; i <= 76; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Semileptonic_%d.root",i));\
}
for (Int_t i = 1; i <= 74; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Semileptonic_%d.root",i));\
}
for (Int_t i = 1; i <= 17; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Semileptonic_%d.root",i));\
}

for (Int_t i = 1; i <= 24; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_3pi0_%d.root",i));\
}
for (Int_t i = 1; i <= 25; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_3pi0_%d.root",i));\
}
for (Int_t i = 1; i <= 5; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_3pi0_%d.root",i));\
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

for (Int_t i = 1; i <= 30; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys_SIGNAL_MIXED_Other_%d.root",i));\
}
for (Int_t i = 1; i <= 31; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-24/mk0*all_phys2_SIGNAL_MIXED_Other_%d.root",i));\
}
for (Int_t i = 1; i <= 7; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-26/mk0*all_phys3_SIGNAL_MIXED_Other_%d.root",i));\
}

chain->Process("signal_vs_bcg_v2.C", "$opt");
.q
EOF