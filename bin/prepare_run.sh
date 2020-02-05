#!/bin/usr/env bash

# RUN AS: bash prepare_run.sh samplesheet.csv BICRO66 fastqdir
# head -1 samplesheet.csv as: BB78b,Caco2_Ctrl_2A,CATCAATC,Homo sapiens,5

inputfile=$1 #samplesheet.csv
inputdir=`dirname $inputfile`
run=$2
fastqdir=$3 # fullpath to dir containing FASTQ files

if [[ ! -d runs ]]; then
    mkdir runs
fi

echo "#!/bin/usr/env bash" > ./runs/run_pipeline_$run.sh
echo >> ./runs/run_pipeline_$run.sh

while read -r line; do
    code=`echo $line|cut -d',' -f1`
    sample=`echo $line|cut -d',' -f2`
    barcode=`echo $line|cut -d',' -f3`
    species=`echo $line|cut -d',' -f4 | awk '{if($0~/[Hh]omo|hs|sapiens|human/){print "human"}else if($0~/[Mm]us|mm|musculus|mouse/){print "mouse"}else{print "Error: species not recognised:", $0 > "/dev/stderr"; exit 1}}'`
    echo bash $BLISS_PATH/bliss.sh "$code" "$sample" "$species" $inputdir/"$code"_*_"$barcode" 60 "$fastqdir" >> ./runs/run_pipeline_$run.sh
    echo >> ./runs/run_pipeline_$run.sh
done < $inputfile
chmod +x ./runs/run_pipeline_$run.sh
