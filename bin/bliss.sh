#!/usr/bin/env bash

# THIS SCRIPT CAN BE CALLED LIKE THIS: bash bliss.sh expName human ../patterns/expName_*_"$barcode" 30 /home/garner1/Work/dataset
################################################################################
# clear
# DEFINING VARIABLES
experiment=$1			# experiment ID found in fastq filename: expID_R1.fastq.gz
genome=$2			# human or mus
patfile=$3			# is the UMI+barcode pattern file used in the linker
quality=$4			# minimum mapping quality desired
fastqDir=$5			# full/path/to/directory containing the compressed fastq file
numbproc=4			# set the desired number of threads to be used while mapping to reference genome
################################################################################
# PREPARE DIRECTORY STRUCTURE
datadir=$PWD/../dataset && mkdir -p $datadir/$experiment
bin=$PWD
python=$PWD/../python #HOME/Work/pipelines/BLISS/python
in=$datadir/$experiment/indata && mkdir -p $in
out=$datadir/$experiment/outdata && mkdir -p $out
aux=$datadir/$experiment/auxdata && mkdir -p $aux
if [ $genome == human ]; then
    refgen=$HOME/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.fa #modify if necessary
fi
if [ $genome == mus ]; then
    #CHANGE THIS TO THE CORRECT LOCATION
    refgen=$HOME/Work/genomes/mm10/mm10.fa # modify if necessary
fi
################################################################################
find $fastqDir -maxdepth 1 -type f -iname "${experiment}*.fastq.gz" | sort > filelist_"$experiment"

numb_of_files=`cat filelist_"$experiment" | wc -l`
r1=`cat filelist_"$experiment" | head -n1`
echo "R1 is " $r1
if [ $numb_of_files == 2 ]; then
    r2=`cat filelist_"$experiment" | tail -n1`
    echo "R2 is " $r2
fi
if [ $numb_of_files == 0 ]; then
    echo "R1 does not exist "
    exit 1
fi
rm filelist_"$experiment"
################################################################################
if [ ! -f $in/r1oneline.fq ]; then
    "$bin"/module/prepare_files.sh  $r1 $in $numb_of_files $r2
fi
"$bin"/module/pattern_filtering.sh $in $out $patfile 
"$bin"/module/prepare_for_mapping.sh $numb_of_files $out $aux $in 
"$bin"/module/mapping.sh $numb_of_files $numbproc $refgen $aux $out $experiment $quality
"$bin"/module/umi_joining.sh $numb_of_files $out $experiment $aux $quality
cat "$datadir"/"$experiment"/outdata/_q"$quality".bed | cut -f-5 |LC_ALL=C uniq -c | awk '{print $2,$3,$4,$5,$6,$1}' | tr " " "," > "$datadir"/"$experiment"/auxdata/aux
#####UMI filtering
cp "$datadir"/"$experiment"/auxdata/aux "$datadir"/"$experiment"/outdata/pre_umi_filtering.csv

"$bin"/module/umi_filter_1.sh "$datadir"/"$experiment"/outdata/pre_umi_filtering.csv "$datadir"/"$experiment"/outdata/q"$quality"_aux
"$bin"/module/umi_filter_2.sh "$datadir"/"$experiment"/outdata/q"$quality"_aux "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-strand-umi-pcr
if [ $genome == human ]; then
    "$bin"/module/umi_filter_3.sh "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-strand-umi-pcr "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-countDifferentUMI.bed
fi
if [ $genome == mus ]; then
    input="$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-strand-umi-pcr
    output="$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-countDifferentUMI.bed
    cat $input | grep -v "_" | sed -e 's/chrX/chr21/g' | sed -e 's/chrY/chr22/g' | cut -f-3 | LC_ALL=C uniq -c | awk '{OFS="\t";print $2,$3,$4,$1}' > $output
fi

echo "Alignment statistics:" >> "$datadir"/"$experiment"/outdata/summary.txt
samtools flagstat "$datadir"/"$experiment"/outdata/*.all.bam >> "$datadir"/"$experiment"/outdata/summary.txt

rm -r "$datadir"/"$experiment"/auxdata* #clean

echo "Number of left and right cuts:" >> "$datadir"/"$experiment"/outdata/summary.txt
cat "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-strand-umi-pcr | grep -v "_" | cut -f4 | sort | uniq -c >> "$datadir"/"$experiment"/outdata/summary.txt
echo "Number of DSB locations:" >> "$datadir"/"$experiment"/outdata/summary.txt
cat "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-countDifferentUMI.bed | grep -v "_" | wc -l >> "$datadir"/"$experiment"/outdata/summary.txt
echo "Number of UMIs:" >> "$datadir"/"$experiment"/outdata/summary.txt
cat "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-strand-umi-pcr | grep -v "_" | wc -l >> "$datadir"/"$experiment"/outdata/summary.txt

name=`echo $patfile|rev|cut -d'/' -f1|rev`
cat "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-strand-umi-pcr |
    awk -F'\t' '{ if ( $4 == "-" ) {$2=$2-1;$3=$3-1;print $0} else {print $0} }'| #correct for the bedtools/bamfile mismatch on negative strands end location 
    tr ' ' '\t' > "$datadir"/"$experiment"/outdata/"$name"__q"$quality"_chr-loc-strand-umi-pcr.tsv
mv "$datadir"/"$experiment"/outdata/q"$quality"_chr-loc-countDifferentUMI.bed "$datadir"/"$experiment"/outdata/"$name"_chr-loc-countDifferentUMI.bed
mv "$datadir"/"$experiment"/outdata/summary.txt "$datadir"/"$experiment"/outdata/"$name"__summary.txt


