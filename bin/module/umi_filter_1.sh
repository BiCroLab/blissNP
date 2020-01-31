#!/usr/bin/env bash

input=$1
output=$2

#python3 ../python/umi_filtering.py "$input" "$output"
python "$BLISS_PATH"/../python/umi_filtering.py "$input" "$output"
