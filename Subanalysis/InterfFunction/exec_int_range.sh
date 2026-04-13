#!/bin/bash

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 sampling_method directory img_folder [t1Min] [t1Max] [t2Min] [t2Max] [steps_t2]"
    exit 1
fi

sampling_method="$1"
directory="$2"
img_folder="$3"

if [ "$#" -ge 4 ]; then
t1Min="$4"
fi

if [ "$#" -ge 5 ]; then
t1Max="$5"
fi

if [ "$#" -ge 6 ]; then
t2Min="$6"
fi

if [ "$#" -ge 7 ]; then
t2Max="$7"
fi

if [ "$#" -ge 8 ]; then
steps_t2="$8"
fi

nohup ./hist_fits_integral_range.exe "$sampling_method" "$directory" "$img_folder" "$t1Min" "$t1Max" "$t2Min" "$t2Max" "$steps_t2" > nohup_$img_folder.log &
