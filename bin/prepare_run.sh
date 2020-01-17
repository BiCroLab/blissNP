#!/bin/usr/env bash

# RUN AS: bash prepare_run.sh samplesheet.csv BICRO66 fastqdir
# head -1 samplesheet.csv as: BB78b,Caco2_Ctrl_2A,CATCAATC,Homo sapiens,5

inputfile=$1 #samplesheet.csv
run=$2
dir=$3 # fullpath to dir containing FASTQ files
rm -f ./runs/run_pipeline_$run.sh

echo "#!/bin/usr/env bash" >> ./runs/run_pipeline_$run.sh
echo >> ./runs/run_pipeline_$run.sh

while read -r line; do
    code=`echo $line|cut -d',' -f1`
    sample=`echo $line|cut -d',' -f2`
    barcode=`echo $line|cut -d',' -f3`
    echo bash $PWD/bliss.sh "$code" human $PWD/../patterns/"$code"_*_"$barcode" 60 "$dir" >> ./runs/run_pipeline_$run.sh
    echo >> ./runs/run_pipeline_$run.sh
done < $inputfile
chmod +x ./runs/run_pipeline_$run.sh
