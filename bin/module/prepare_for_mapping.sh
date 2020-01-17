#!/usr/bin/env bash
numb_of_files=$1
out=$2
aux=$3
in=$4

echo 'Parse the fastq files, filtering and trimming ...'
if [[ $numb_of_files == 1 ]]; then # IF SE READS
    cat $out/filtered.r1.fa | cut -d':' -f-7 | tr '>' '@' | paste - -|awk '{print $1,$NF}'|
    LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -k1,1 > $aux/ID_genomic
    LC_ALL=C join $aux/ID_genomic $in/r1oneline.fq | tr " " "\n" | grep -v "1:[YN]:0:" | paste - - - - -| 
    awk '{print $1,$2,$4,substr($5,length($5)-length($2)+1,length($5))}' | tr " " "\n" > $aux/r1.2b.aln.fq
fi
if [[ $numb_of_files == 2 ]]; then # PE READS
    cat $out/filtered.r1.fa | cut -d':' -f-7 | tr '>' '@' | paste - -|awk '{print $1,$NF}'|
    LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -k1,1 > $aux/ID_genomic
    LC_ALL=C join $aux/ID_genomic $in/r1oneline.fq | tr " " "\n" | grep -v "1:[YN]:0:" | paste - - - - -| 
    awk '{print $1,$2,$4,substr($5,length($5)-length($2)+1,length($5))}' | tr " " "\n" > $aux/r1.2b.aln.fq & pid1=$!
    cat $aux/ID_genomic | cut -d' ' -f1 | LC_ALL=C join - $in/r2oneline.fq | tr " " "\n" | grep -v "2:N:0:" > $aux/r2.2b.aln.fq & pid2=$!
    wait $pid1
    wait $pid2
fi
echo 'Done! Ready to be aligned to the reference genome!'

