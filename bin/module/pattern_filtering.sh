#!/usr/bin/env bash
in=$1
out=$2
patfile=$3

echo "Filterig reads based on patterns ..."
cat $in/r1.fa | parallel --tmpdir $HOME/tmp --block 100M -k --pipe -L 2 "scan_for_matches $patfile - " > $out/filtered.r1.fa
echo "Done"
