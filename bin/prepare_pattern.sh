#!/bin/usr/env bash

inputfile=$1			# samplesheet.csv

while read -r line; do
    f1=`echo $line|cut -d',' -f1`
    f2=`echo $line|cut -d',' -f2|tr -d '_/'`
    f3=`echo $line|cut -d',' -f3`
    ##############################################################################################
    # HERE CHANGE THE NUMBER OF MISMATCHES IN [1,0,0] FROM 1 TO OTHER VALUES, IE [2,0,0] OR [0,0,0]:
    #  echo ^ 8...8 "$f3"[numb_of_mismatches,0,0] 1...1000 $ > "$f1"_"$f2"_"$f3"
    ###############################################################################################
    echo ^ 8...8 "$f3"[1,0,0] 1...1000 $ > "$f1"_"$f2"_"$f3"
    ###############################
done < $inputfile
