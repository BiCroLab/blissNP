# BLISS is a versatile and quantitative method for genome-wide profiling of DNA double-strand breaks
This repository provides step-by-step instructions on how to process [BLISS](https://www.nature.com/articles/ncomms15058) data.
We assume that the user is using the Linux operating system (modifications to the scripts might be necessary if using Mac OS X). Commands starting with ```$``` are executed in the command line.

## Setting up the pipeline

* Follow [these instructions](https://docs.conda.io/en/latest/miniconda.html) to install Miniconda (providing you with conda, Python and the basic packages they require)
* Follow [these instructions](http://blog.theseed.org/servers/2010/07/scan-for-matches.html) to install scan_for_matches
* Install bedtools (the version used in testing the pipeline is v2.29.2), samtools (the version used in testing the pipeline is v1.9) and gnu parallel (the version used in testing the pipeline is v20190722):
```
$ conda install -c bioconda bedtools
$ conda install -c bioconda samtools
$ sudo apt install parallel
```
* Clone or download this repository:
```
$ git clone https://github.com/garner1/bliss_NP.git
$ cd ./bliss_NP
```
* Make an index of the reference genome of interest using bwa:
```
$ bwa index /path/to/genome_of_interest.fa
```
* For demonstration and testing purposes, we prepared a small dataset contained in the *bliss_NP/test* directory. The directory contains a fastq file in the *bliss_NP/test/fastq* directory and a configuration file in the *bliss_NP/test/samplesheet* directory
* Move into the  *bliss_NP/test* directory and run the *test_pipeline* script
```
$ cd ./test
$ bash test_pipeline
```
* The screen output should look like this
```
R1 is  ../test/fastq/test.fastq.gz
Filterig reads based on patterns ...
Done
Parse the fastq files, filtering and trimming ...
Done! Ready to be aligned to the reference genome!
Aligning reads to the reference genome ...
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -v 1 -t 4 /home/garner1/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.fa /home/garner1/Work/pipelines/bliss_NP/bin/../dataset/TEST/auxdata/r1.2b.aln.fq
[main] Real time: 3.023 sec; CPU: 3.093 sec
Done
Selecting unique UMIs
Done
Done with filtering UMIs!
```
* In ./bin/bliss.sh set the number of threads used during alignment in line 13 (default 24)
* In ./bin/bliss.sh set the location of the human or mus reference genome (fasta file) in lines 24 and 28
* In ./bin/prepare_pattern.sh set the number of mismatches allowed in the barcode in line 10
