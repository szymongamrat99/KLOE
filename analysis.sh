#!/bin/bash

# Sprawdź, czy podano wymagane argumenty
if [ "$#" -lt 2 ]; then
    echo "Użycie: $0 <file type> <smearing option>"
    exit 1
fi

# Pobierz argumenty
arg1=$1
arg2=$2

# Pobierz aktualną datę w formacie YYYY-MM-DD (bez timestampu)
date_folder=$(date +"%Y-%m-%d")

# Pobierz aktualną datę i czas w formacie YYYY-MM-DD_HH-MM-SS
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# Utwórz folder log/<data> jeśli nie istnieje
log_dir="log/${date_folder}"
mkdir -p "${log_dir}"

# Skonstruuj nazwę pliku logu
log_file="${log_dir}/nohup_${arg1}_${arg2}_${timestamp}.log"

# Uruchom komendę z nohup
nohup ./execute_analysis.sh < input_init_analysis_${arg1}.txt > "${log_file}" &

# Wyświetl informację o uruchomieniu
echo "Skrypt został uruchomiony w tle. Log: ${log_file}"