# Run alignment and processing pipeline for transcription factor ChIP-seq
===========================================================
## This pipeline runs QC, alignment, and peak calling of transcription factor ChIP-seq data  

Adapted from https://github.com/snakemake-workflows/chipseq for use in the Philpott lab  

## Authors  
- Antonie Vietor (@AntonieV)
- Laura Woods (@lmw48)

## Usage instructions  
1. Check docker is installed, if not then install docker

2. Build docker image in current directory: docker build --tag chipseq:2.0 .
- NB: if running on ap05 this step is unnecessary as the image already exists  

3. Clone this repository with
```
git clone https://gitlab.developers.cam.ac.uk/jcbc/csci/philpottlab/chip-seq_2021
```

4. Modify samples.tsv and units.tsv in config folder to reflect your samples. Examples are included in config. Fastq file names should be in the format basename_unit_R{1,2}.fq.gz. In the units sheet fq file names should be appended with /fastq/.  
- Basename is the sample name that will be appended to output files  
- Unit is the sequencing run the file comes from. This is to allow merging of technical replicates
- e.g. paired-end fq files for sample BE2C-Rep-1-Reverted_S4 sequenced in runs 1 and 2:  
	- /fastq/BE2C-Rep-1-Reverted_S4_1_R1.fq.gz  
	- /fastq/BE2C-Rep-1-Reverted_S4_1_R2.fq.gz  
	- /fastq/BE2C-Rep-1-Reverted_S4_2_R1.fq.gz  
	- /fastq/BE2C-Rep-1-Reverted_S4_2_R2.fq.gz  

5. Get bowtie2 genome index 
- either build with bowtie2-build, or download from http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer  
- Genome index may already be built for you - if working on Philpott lab Linux ap05 check in /mnt/tank/genomeIndex  

6. Modify config.yaml to specify
- genome build  
- path to bowtie2 index  
- peak calling parameters  

7. Set up docker container, use -v to mount current directory containing Dockerfile, Snakefile, and data/fastq/*.fq.gz
Mount pipeline directory at /work
Mount fastq files at /fastq
Mount output directory at /output
Mount config files at /work/config
Mount bowtie2 genome index at /genomeIndex

```
docker run --rm -it -v <pipeline-directory>:/work -v <fastq-directory>:/fastq -v <output-directory>:/output -v <config-files>/:/work/config -v <genome-index-files>:/genomeIndex chipseq:2.0
```

8. Run snakemake

Dry run to show steps that will be run:
	snakemake -np --cores 4

Run pipeline, print shell commands:  

```
snakemake -p --cores 20
```
20 cores recommended, can use "all" if available  

If output directory is on FUSE mounted google drive, FastQC and trimgalore will throw errors.  
