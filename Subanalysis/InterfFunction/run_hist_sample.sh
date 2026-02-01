#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <samples>"
    exit 1
fi

echo "2 0 1 $1" | ./hist_sample_draw
