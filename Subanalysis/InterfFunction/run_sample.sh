#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <integral_limit> <samples>"
    exit 1
fi

echo "0 $1 $2" | ./bttf_sample
