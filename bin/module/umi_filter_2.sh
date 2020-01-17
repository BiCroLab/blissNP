#! /usr/bin/env bash

input=$1
output=$2

cat $input | tr -d "," | tr -d "'" | tr -d "[" | tr -d "]" | tr " " "\t" > $output
