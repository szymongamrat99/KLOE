#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <integral_limit> <samples>"
    exit 1
fi

echo "1 0 0 $1 $2 3" | ./bttf_sample
