# BLISS is a versatile and quantitative method for genome-wide profiling of DNA double-strand breaks
This repository provides step-by-step instructions on how to process [BLISS](https://www.nature.com/articles/ncomms15058) data.
We assume that the user is using the Linux operating system (modifications to the scripts might be necessary if using Mac OS X). Commands starting with ```$``` are executed in the command line.

## Setting up the pipeline

* Follow [these instruction](https://docs.conda.io/en/latest/miniconda.html) to install Miniconda (providing you with conda, Python and the basic packages they require)
* Clone or download this repository
```
$ git clone https://github.com/garner1/bliss_NP.git
```

* In ./bin/bliss.sh set the number of threads used during alignment in line 13 (default 24)
* In ./bin/bliss.sh set the location of the human or mus reference genome (fasta file) in lines 24 and 28
* In ./bin/prepare_pattern.sh set the number of mismatches allowed in the barcode in line 10
