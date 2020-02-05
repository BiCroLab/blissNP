#!/bin/usr/env bash

inputfile=$1			# samplesheet.csv
inputdir=`dirname $inputfile`

# Check fields
awk -F ',' '(NF!=5){print "Error: wrong number of columns. Expected 5 and found", NF > "/dev/stderr"; exit 1 }' $inputfile

# Check identifier uniqueness
awk -F ',' 'BEGIN{n=0;}{vec[n+1] = $2; n++}END{cnt=0;for(i=1; i<=n; i++){dup=0;for(ii=1; ii<=n; ii++){if(vec[i]==vec[ii] && i!=ii){dup=1}} if(dup==0){cnt++}} if(cnt!=n){print "Error: the sample identifier column (2) contains duplicates." > "/dev/stderr"; exit 1 }}' $inputfile

if [[ $? -ne 0 ]]; then
    exit 1
fi

while read -r line; do
    f1=`echo $line|cut -d',' -f1`
    f2=`echo $line|cut -d',' -f2|tr -d '_/'`
    f3=`echo $line|cut -d',' -f3`
    f6=`echo $line|cut -d',' -f5`
    ##############################################################################################
    # HERE CHANGE THE NUMBER OF MISMATCHES IN [1,0,0] FROM 1 TO OTHER VALUES, IE [2,0,0] OR [0,0,0]:
    #  echo ^ 8...8 "$f3"[numb_of_mismatches,0,0] 1...1000 $ > "$f1"_"$f2"_"$f3"
    ###############################################################################################
    echo ^ 8...8 "$f3"["$f6",0,0] 1...1000 $ > "$inputdir"/"$f1"_"$f2"_"$f3"
    ###############################
done < $inputfile
