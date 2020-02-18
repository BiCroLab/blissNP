# (s)BLISS is a versatile and quantitative method for genome-wide profiling of DNA double-strand breaks
This repository provides step-by-step instructions on how to process [BLISS](https://www.nature.com/articles/ncomms15058) or sBLISS data.

We assume that the user is using the Linux operating system (modifications to the scripts might be necessary if using Mac OS X).

Commands starting with ```$``` are executed in the command line.

The specifics of the machines necessary to run the pipeline will depend mainly on the available time and space resources. We recommend at least 15Gb of RAM and dual core processors.  

If you use this sofware please cite the original manuscript [1].

[1]: [Yan, W.X., Mirzazadeh, R., Garnerone, S., Scott, D., Schneider, M.W., Kallas, T., Custodio, J., Wernersson, E., Li, Y., Gao, L. and Federova, Y., 2017. BLISS is a versatile and quantitative method for genome-wide profiling of DNA double-strand breaks. Nature communications, 8(1), pp.1-9.](https://www.nature.com/articles/ncomms15058) 

## Setting up the pipeline

* Follow [these instructions](https://docs.conda.io/en/latest/miniconda.html) to install Miniconda (providing you with conda, Python and the basic packages they require)
* Follow [these instructions](http://blog.theseed.org/servers/2010/07/scan-for-matches.html) to install scan_for_matches

* Create a dedicated conda environment
```
$ conda create --name sbliss
```
* Activate the environment
```
$ conda activate sbliss
```
* Install bedtools (the version used in testing the pipeline is v2.29.2)
```
$ conda install -c bioconda bedtools
```
* Install samtools (the version used in testing the pipeline is v1.9)
```
$ conda install -c bioconda samtools
```
* Install bwa (the version used in testing the pipeline is v0.7.17-r1188)
```
$ conda install -c bioconda bwa
```
* Install gnu parallel (the version used in testing the pipeline is v20190722):
```
$ sudo apt install parallel
```
* Make sure that the installed software executables can be found in their conda path directory (if necessary edit your .bashrc, or the equivalent for a different shell): open the file ```~/.bashrc``` in your favorite editor and add ```export PATH="~/miniconda3/bin:$PATH"``` at the end of the file (if miniconda was installed in a different location you have to modify the path provided in .bashrc). Save and close the file. Then source your .bashrc file:
```
$ source .bashrc
```
* Make an index of the reference genome of interest using bwa:
```
$ bwa index /path/to/genome_of_interest.fa
```
* Clone or download this repository:
```
$ git clone https://github.com/BiCroLab/blissNP.git
$ cd ./blissNP
```
* Configure the BLISS\_PATH by opening the file ```~/.bashrc``` in your favorite editor and add ```export BLISS_PATH="<path/to/blissNP/bin>"``` at the end of the file. Save and close the file. Then source your .bashrc file:
```
$ source .bashrc
```
__IMPORTANT:__ To configure the pipeline for general usage you should:
* In ```blissNP/bin/bliss.sh``` set the number of threads (default to 4) used during alignment on [this line](https://github.com/BiCroLab/blissNP/blob/f1aec60e1c4d2631fb4add82505deb06598c0017/bin/bliss.sh#L12) 
* In ```blissNP/bin/bliss.sh``` set the location of the human reference genome on [this line](https://github.com/BiCroLab/blissNP/blob/f1aec60e1c4d2631fb4add82505deb06598c0017/bin/bliss.sh#L22) or of the mouse reference genome [on this line](https://github.com/BiCroLab/blissNP/blob/f1aec60e1c4d2631fb4add82505deb06598c0017/bin/bliss.sh#L26)
* Prepare a sample sheet configuration file (CSV format) with five fields:
  1. FASTQ file base name
  2. sample ID
  3. sample barcode
  4. organism of interest (must be one of: homo sapiens, hs or human, mus musculus, mm, or mouse)
  5. number of mismatches allowed in the sample barcode

To run the pipeline on your dataset, use the following command:
```
$ bash "$BLISS_PATH"/prepare_pattern.sh <sample sheet>
$ bash "$BLISS_PATH"/prepare_run.sh <sample sheet> <run name> <full/path/to/directory_with_fastq_files>
$ bash ./runs/run_pipeline_<run name>
```

## Test demonstration
For demonstration and testing purposes, we prepared a small dataset contained in the ```blissNP/test``` directory. The directory contains a fastq file in the *blissNP/test/fastq* directory and a configuration csv file in the *bliss_NP/test/samplesheet* directory. The configuration file has 5 fields: *experiment ID, sample ID, sample barcode,  genome of interest, number of mismatches allowed in the sample barcode*. The *genome of interest* has to be one of: [Hh]omo|hs|sapiens|human|[Mm]us|mm|musculus|mouse.

* Move into the  *blissNP/test* directory and run the *test_pipeline* script
```
$ cd ./test
$ bash test_pipeline
```
* In the *blissNP/test/runs* directory a file named ```run_pipeline_TEST.sh``` is generated which contains the command line that runs the pipeline:
```
bash /home/garner1/Work/pipelines/blissNP/bin/bliss.sh TEST DMSO human samplesheet/TEST_DMSO_AGCCATCA 60 fastq
```
The arguments passed to the *bliss.sh* script are: ```experiment ID```, ```sample ID```, ```genome ID```, ```path/to/UMI-barcode/pattern/file```, ```threshold on the quality of alignment``` and ```/path/to/dir/containing/fastq```.

* The screen output should look like this
```
R1 is  ../test/fastq/test.fastq.gz
Filtering reads based on patterns ...
Done
Parse the fastq files, filtering and trimming ...
Done! Ready to be aligned to the reference genome!
Aligning reads to the reference genome ...
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -v 1 -t 4 /path/to/reference/genome.fa /path/to/downloaded/repo/blissNP/bin/../dataset/TEST/auxdata/r1.2b.aln.fq
[main] Real time: 3.023 sec; CPU: 3.093 sec
Done
Selecting unique UMIs
Done
Done with filtering UMIs!
```

## Output from the test
The output from the test dataset is located in ```blissNP/dataset/TEST/outdata```.
Relevant files are (for the test dataset: experiment ID = TEST, sample ID = DMSO, sample barcode = AGCCATCA):
* ```sampleID.all.bam```: bam file before UMI deduplication
* ```sampleID.q60.bam```: bam file after UMI deduplication
* ```experimentID_sampleID_samplebarcode__summary.txt```: summary of the analysis
* ```experimentID_sampleID_samplebarcode_chr-loc-countDifferentUMI.bed```: list of (chr,start,end,number of unique DSB at this location) for each DSB location
* ```experimentID_sampleID_samplebarcode__q60_chr-loc-strand-umi-pcr.tsv```: list of (chr,start,end,strand,UMI,number of PCR duplicates) for each DSB


