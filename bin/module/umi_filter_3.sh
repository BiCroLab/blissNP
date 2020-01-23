#! /usr/bin/env bash

input=$1
output=$2

cat $input |
    awk -F'\t' '{ if ( $4 == "-" ) {$2=$2-1;$3=$3-1;print $0} else {print $0} }'| #correct for the bedtools/bamfile mismatch on negative strands end location 
    tr ' ' '\t' |
    grep -v "_" | sed -e 's/chrX/chr23/g' | sed -e 's/chrY/chr24/g' | cut -f-3 | LC_ALL=C uniq -c | awk '{OFS="\t";print $2,$3,$4,$1}' > $output
