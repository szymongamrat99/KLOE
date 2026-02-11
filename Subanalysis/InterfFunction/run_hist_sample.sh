#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <samples> <method>"
    exit 1
fi

echo "2 0 1 $1 $2" | ./hist_sample_draw
