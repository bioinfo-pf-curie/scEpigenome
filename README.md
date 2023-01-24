# scEpigenome

**Institut Curie - single-cell Epigenomics analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.11-blue.svg)](https://multiqc.info/)
[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7443721.svg)](https://doi.org/10.5281/zenodo.7443721)-->

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.

The goal of this pipeline is to process multiple type of single-cell epigenomics profiles, including scCut&Tag (10X, CellenOne, inDrop), and scChIP-seq.

### Pipline summary

TODO

### Quick help

TODO

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda
```

#### Run the pipeline from a sample plan

```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --aligner 'star' --counts 'star' --genome 'hg38' --outDir MY_OUTPUT_DIR -profile conda
```

#### Run the pipeline on a computational cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --aligner 'star' --counts 'star' --genome 'hg19' --outDir MY_OUTPUT_DIR -profile singularity,cluster" | qsub -N rnaseq
```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option. See the [full documentation](docs/profiles.md) for details.

```
## Run the pipeline locally, using the paths defined in the configuration for each tool (see conf/path.config)
-profile path --globalPath INSTALLATION_PATH 

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityImagePath SINGULARITY_IMAGE_PATH 

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda --condaCacheDir CONDA_CACHE 
```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs, **with no header**.


SAMPLE_ID,SAMPLE_NAME,PATH_TO_R1_FASTQ,[PATH_TO_R2_FASTQ]


### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/referenceGenomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the single-cell custom and bioinformatics facilities of the Institut Curie (L. Hadj-Abed, P. Prompsy, C. Vallot, N. Servant)

#### Citation

If you use this pipeline for your project, please cite it using the following doi: **TODO**
Do not hesitate to use the Zenodo doi corresponding to the version you used !

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.
