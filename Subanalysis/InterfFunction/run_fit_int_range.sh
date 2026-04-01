#!/bin/bash

options=("free" "fixed_re" "fixed_im")

select opt in "${options[@]}"; do
    case $opt in
        "free")
            echo "You chose free fit range."
            break
            ;;
        "fixed_re")
            echo "You chose fixed Re fit range."
            break
            ;;
        "fixed_im")
            echo "You chose fixed Im fit range."
            break
            ;;
        *) echo "Invalid option $REPLY";;
    esac
done

fitMax=(10 20 40 60 80 100 120 150 200 250 280 300)

for z in "${fitMax[@]}"; do
    nohup bash -c "
    ./hist_fits_integral_range << EOF
    150_points_${opt}_fit_range_${z}
    1
    0.5
    300
    150
    ${z}
    EOF
    " > "log_${opt}_fit_range_${z}.txt" 2>&1 &
done