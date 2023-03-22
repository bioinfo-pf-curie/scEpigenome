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

This pipeline process 3 types of epigenomics data : i) scChIPseq, ii) scCUT&Tag done by an indrop fashion using a microfluidics device & iii) scCUT&Tag done in a 10X like fashion using 10XGenomics device. 

The pipeline goes from raw reads (fastq, paired end) to exploitable count matrices. The multiple steps involved in the pipeline are :

1. Align barcode read parts on barcode index libraries
3. Align genomic read parts on the genome
4. Assignation of cell barcodes to aligned read
5. Removal of duplicates (PCR & RT)
6. Removal of reads based on window screening (if Read2 was unmapped)
7. Removal of black regions (repeated regions, low mappability regions)
8. Counting (Generation of count matrix) in bins or by TSS (transcription start sites) as an approximation of genes 
9. Generation of coverage file (bigwig) (CPM normalization)
10. Reporting

### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `main.nf` [stupefied_darwin] - revision: aa905ab621
=======================================================

Usage:

Mandatory arguments:
--reads [file]                   Path to input data (must be surrounded with quotes)
--samplePlan [file]              Path to sample plan file if '--reads' is not specified
--genome [str]                   Name of the reference genome. See the `--genomeAnnotationPath` to defined the annotation path
-profile [str]                   Configuration profile to use (multiple profiles can be specified with comma separated values)
--protocol [str]                 Chose between: 'scchip_indrop', 'sccuttag_indrop', 'sccuttag_10X'

Skip options: All are false by default
--skipSoftVersion [bool]         Do not report software version
--skipMultiQC [bool]             Skip MultiQC
--skipBigWig                     Skip bigwig file generation (decrease the running time)

Other options:
--outDir [dir]                  The output directory where the results will be saved
-w/--work-dir [dir]             The temporary directory where intermediate data will be saved

--removeBlackRegion [bool]       Remove black region. Default is true
--binSize [str]                  Bin size to use (in base pairs). Default is '50000,250'
--tssWindow [int]                Number of base pairs arround TSS. Default is 5000
--minReadsPerCellmatrix [int]    Minimum number of reads per cell for the matrices. Default is 100
--minReadsPerCellmqc [int]       Minimum number of reads to account for a cell in the multiqc report. Default is 1000

=======================================================
Available profiles
-profile test                    Run the test dataset
-profile conda                   Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile multiconda              Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile path                    Use the installation path defined for all tools. Use `--globalPath` to define the installation path
-profile multipath               Use the installation paths defined for each tool. Use `--globalPath` to define the installation path
-profile docker                  Use the Docker images for each process
-profile singularity             Use the Singularity images for each process. Use `--singularityPath` to define the path of the singularity containers
-profile cluster                 Run the workflow on the cluster, instead of locally

```


### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow:

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda --protocol 'scchip_indrop' 
```

#### Run the pipeline from a sample plan

```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --genome 'hg38' --outDir MY_OUTPUT_DIR -profile conda
```

#### Run the pipeline on a computational cluster

```
echo "nextflow run main.nf --reads '*.R{1,2}.fastq.gz' --genome 'hg19' --outDir MY_OUTPUT_DIR -profile singularity,cluster" | qsub -N scchip
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
