#!/usr/bin/env bash

################################################################################
# clear
# DEFINING VARIABLES
experiment=$1		# experiment ID found in fastq filename: expID_R1.fastq.gz
sample=$2			# sample ID
genome=$3			# human or mouse
patfile=$4			# is the UMI+barcode pattern file used in the linker
quality=$5			# minimum mapping quality desired
fastqDir=$6			# full/path/to/directory containing the compressed fastq file
numbproc=4			# set the desired number of threads to be used while mapping to reference genome
################################################################################
# PREPARE DIRECTORY STRUCTURE
datadir=./dataset && mkdir -p $datadir/$experiment/$sample
in=$datadir/$experiment/$sample/indata && mkdir -p $in
out=$datadir/$experiment/$sample/outdata && mkdir -p $out
aux=$datadir/$experiment/$sample/auxdata && mkdir -p $aux

if [ $genome == human ]; then
    #CHANGE THIS TO THE CORRECT LOCATION/INDEX_BASENAME
    refgen=$HOME/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.fa # modify if necessary
fi
if [ $genome == mouse ]; then
    #CHANGE THIS TO THE CORRECT LOCATION/INDEX_BASENAME
    refgen=$HOME/Work/genomes/mm10/mm10.fa # modify if necessary
fi
if [ ! -f "$refgen".sa ]; then # Check the last index file generated by bwa
    echo "Index not found! Please generate one or modify the 'refgen' value in the 'bliss.sh' file."
    exit 1
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
    "$BLISS_PATH"/module/prepare_files.sh  $r1 $in $numb_of_files $r2
fi
"$BLISS_PATH"/module/pattern_filtering.sh $in $out $patfile
if [[ $? -eq 0 ]]; then rm $in/r1.fa; fi 
"$BLISS_PATH"/module/prepare_for_mapping.sh $numb_of_files $out $aux $in
if [[ $? -eq 0 ]]; then rm $in/r1oneline.fq $in/r2oneline.fq >& /dev/null; fi 
"$BLISS_PATH"/module/mapping.sh $numb_of_files $numbproc $refgen $aux $out $sample $quality
if [[ $? -eq 0 ]]; then rm $aux/r1.2b.aln.fq $aux/r2.2b.aln.fq >& /dev/null; fi 
"$BLISS_PATH"/module/umi_joining.sh $numb_of_files $out $sample $aux $quality
if [[ $? -eq 0 ]]; then rm $out/filtered.r1.fa; fi 
cat $out/_q"$quality".bed | cut -f-5 |LC_ALL=C uniq -c | awk '{print $2,$3,$4,$5,$6,$1}' | tr " " "," > $aux/aux
if [[ $? -eq 0 ]]; then rm $out/_q"$quality".bed; fi 
#####UMI filtering
cp "$datadir"/"$experiment"/"$sample"/auxdata/aux $out/pre_umi_filtering.csv

"$BLISS_PATH"/module/umi_filter_1.sh $out/pre_umi_filtering.csv $out/q"$quality"_aux
"$BLISS_PATH"/module/umi_filter_2.sh $out/q"$quality"_aux $out/q"$quality"_chr-loc-strand-umi-pcr
if [[ $? -eq 0 ]]; then rm $out/q"$quality"_aux; fi 

if [ $genome == human ]; then
    "$BLISS_PATH"/module/umi_filter_3.sh $out/q"$quality"_chr-loc-strand-umi-pcr $out/q"$quality"_chr-loc-countDifferentUMI.bed
fi
if [ $genome == mouse ]; then
    input=$out/q"$quality"_chr-loc-strand-umi-pcr
    output=$out/q"$quality"_chr-loc-countDifferentUMI.bed
    cat $input | grep -v "_" | sed -e 's/chrX/chr21/g' | sed -e 's/chrY/chr22/g' | awk '{if($3=="+"){printf "%s\t%d\t%d\n", $1, $2, $3+1}else{printf "%s\t%d\t%d\n", $1, $3, $3+1}}' | sort -k1,1 -k2,2g -k3,3g | LC_ALL=C uniq -c | awk '{OFS="\t";print $2,$3,$4,$1}' > $output
fi

echo "Alignment statistics:" >> $out/summary.txt
samtools flagstat --threads $numbproc $out/*.all.bam >> $out/summary.txt

rm -r "$aux"* #clean

echo "Number of reads on plus strand and on minus strand:" >> $out/summary.txt
cat $out/q"$quality"_chr-loc-strand-umi-pcr | grep -v "_" | cut -f4 | sort | uniq -c >> $out/summary.txt
echo "Number of DSB locations:" >> $out/summary.txt
cat $out/q"$quality"_chr-loc-countDifferentUMI.bed | grep -v "_" | wc -l >> $out/summary.txt
echo "Number of UMIs:" >> $out/summary.txt
cat $out/q"$quality"_chr-loc-strand-umi-pcr | grep -v "_" | wc -l >> $out/summary.txt

name=`echo $patfile|rev|cut -d'/' -f1|rev`
cat $out/q"$quality"_chr-loc-strand-umi-pcr |
    awk -F'\t' '{ if ( $4 == "-" ) {$2=$2-1;$3=$3-1;print $0} else {print $0} }'| #correct for the bedtools/bamfile mismatch on negative strands end location 
    tr ' ' '\t' > $out/"$name"__q"$quality"_chr-loc-strand-umi-pcr.tsv
mv $out/q"$quality"_chr-loc-countDifferentUMI.bed $out/"$name"_chr-loc-countDifferentUMI.bed
mv $out/summary.txt $out/"$name"__summary.txt
