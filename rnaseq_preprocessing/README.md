# Snakemake workflow: rna-seq-star
## This workflow aligns RNA-seq data using STAR and performs QC steps  

## Adapted from the following sources for use in the Philpott lab

https://github.com/snakemake-workflows/rna-seq-star-deseq2  

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2)
[![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)](https://cdn.rawgit.com/snakemake-workflows/rna-seq-star-deseq2/master/.test/report.html)

## Authors

* Johannes Köster (@johanneskoester), https://koesterlab.github.io
* Sebastian Schmeier (@sschmeier), https://sschmeier.com
* Jose Maturana (@matrs)
* Laura Woods (@lmw48)

## Usage

### Simple

#### Step 1: Install workflow

Download the rnaseq_preprocessing directory

#### Step 2: Configure workflow

##### sample and units tables  
Example samples and units tables can be found in the config directory  
  
Fastq file names should be in the format basename.unit.R{1,2}.fq.gz where basename is the sample name that will be appended to all output files, and unit is the sequencing run the sample belongs to. In the units sheet the filename should be appended with /fastq/ to point to the appropriate directory.  
   
In the example below the sample with basename BE2C_ASCL1_WTA1_R1_uninduced_S11 was sequenced paired-end in 2 sequencing runs (s_1 and s_2)  
  
e.g.  
- /fastq/BE2C_ASCL1_WTA1_R1_uninduced_S11.s_1.R1.fq.gz  
- /fastq/BE2C_ASCL1_WTA1_R1_uninduced_S11.s_1.R2.fq.gz  
- /fastq/BE2C_ASCL1_WTA1_R1_uninduced_S11.s_2.R1.fq.gz  
- /fastq/BE2C_ASCL1_WTA1_R1_uninduced_S11.s_2.R2.fq.gz  

#### Step 2a: Build and run docker image
A docker file is provided with the workflow to set up a working environment with all required software. To build the image run the following command.  

```
docker build -t rnaseq /path/to/dockerfile
```

Set up a docker container using the following command, which runs the docker container interactively, mounts: the workflow directory, location of fastq files, output directory, and location of genome indexes and annotations:  

```
docker run -it --rm -v <pipeline-directory>:/work -v <ref-directory>:/ref -v /path/to/fastq:/fastq -v /path/to/output/:/output -v /path/to/config/:/work/config rnaseq:latest
```

The ref directory should contain:  
- prebuilt star index e.g. ~/ref/starIndex/hg38  
- chrom sizes file which can be downloaded from UCSC  
- annotation gtf file which can be downloaded from UCSC  

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via  

```
snakemake -np --cores 4  
```

Run the workflow with:  

```
snakemake -p --cores 20  
```

You can change the number of cores depending on your compute power and the size of you dataset  
