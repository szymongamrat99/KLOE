#!/bin/bash

ANALYSIS_TYPE=$1
DATE=$2
FIRST_IDX=$3
LAST_IDX=$4
NUM_PROC=$5

# Poprawione sprawdzanie argumentów
if [ "$#" -lt 5 ]; then
    echo "Error: Usage $0 <analysis_type> <date> <first_idx> <last_idx> <num_proc>"
    exit 1
fi

RESULT_DIR=parallel_logs_${ANALYSIS_TYPE}_${FIRST_IDX}_${LAST_IDX}

# Używamy seq zamiast {..}, aby zmienne zadziałały
# Przekierowujemy wyjście do logu, aby nohup miał co zapisywać
nohup parallel -j $NUM_PROC --results $RESULT_DIR/ "cat parameters.txt | ./execute_analysis.sh 8 /data/ssd/gamrat/DBV-26/$ANALYSIS_TYPE/$DATE/job_v26_${ANALYSIS_TYPE}_4_inv_pb_{}.txt" ::: $(seq $FIRST_IDX $LAST_IDX) > "nohup_${ANALYSIS_TYPE}_${FIRST_IDX}_${LAST_IDX}.log" 2>&1 &

echo "Analiza puszczona w tle (nohup). Logi w: nohup_${ANALYSIS_TYPE}_${FIRST_IDX}_${LAST_IDX}.log"
