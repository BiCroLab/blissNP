#!/usr/bin/env bash

input=$1
output=$2

/usr/local/share/anaconda3/bin/python3 ~/Work/pipelines/BLISS/python/umi_filtering.py "$input" "$output"
