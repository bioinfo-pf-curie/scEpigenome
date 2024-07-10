# scEpigenome

**Institut Curie - single-cell Epigenomics analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.22-blue.svg)](https://multiqc.info/)
[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7443721.svg)](https://doi.org/10.5281/zenodo.7443721)-->

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.

The goal of this pipeline is to process multiple type of single-cell epigenomics profiles, including scCut&Tag (10X, inDrop), and scChIP-seq.

### Pipline summary

This pipeline process 2 types of epigenomics data : scChIPseq and scCUT&Tag, generated with various protocoles including 10X barcoding, indrop microfluidics protocols or plate systems.

The pipeline goes from raw reads (fastq, paired end) to exploitable count matrices as follow:

1. Align barcode read parts on barcode index libraries
3. Align genomic read parts on the genome
4. Assignation of cell barcodes to aligned read
5. Removal of duplicates (PCR & extra duplicates)
7. Removal of black regions (repeated regions, low mappability regions)
8. Counting (Generation of count matrix) in bins or by TSS (transcription start sites) as an approximation of genes 
9. Generation of coverage file (bigwig) (CPM normalization)
10. Peak Calling (pseudo-bulk)
11. Reporting

### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 22.10.6
Launching `main.nf` [sad_magritte] DSL2 - revision: 59d670a8d1
=======================================================


Usage:
	
The typical command for running the pipeline is as follows:
	
nextflow run main.nf --reads PATH --samplePlan PATH --genome STRING --protocol STRING
			
MANDATORY ARGUMENTS:
--genome     STRING                                                                Name of the reference genome.
--protocol   STRING [scchip_indrop, sccuttag_indrop, sccuttag_10X, sccutag_plate]  Specify which protocol to run
--reads      PATH                                                                  Path to input data (must be surrounded with quotes)
--samplePlan PATH                                                                  Path to sample plan (csv format) with raw reads (if `--reads` is not specified)

REFERENCES:
--genomeAnnotationPath   PATH      Path to genome annotations folder
--effGenomeSize          INTEGER   Effective genome size
--fasta                  PATH      Path to genome fasta file
--geneBed                PATH      Path to gene file (BED)
--genomeAnnotationPath   PATH      Path to genome annotations folder
--gtf                    PATH      Path to GTF annotation file. Used in HOMER peak annotation
--starIndex              PATH      Indexes for STAR aligner
--bwaIndex               PATH      Indexes for Bwa-mem aligner

INPUTS:
--batchSize              INTEGER   Number of cells to merge together to work in batch (only for plate protocols)

ALIGNMENT:
--aligner                STRING    Aligner to use ('star' or 'bwa-mem2')

BARCODES:
--mapqBarcode            INTEGER   Mapping quality for the barcode alignment (40)
--barcodeTag             STRING    Barcode tag ('XB')

FILTERING:
--blackList              PATH      Path to black list regions (.bed). See the genome.config for details
--mapq                   INTEGER   Minimum mapping quality after reads alignment (20)
--rmSingleton                      Remove singleton
--extraDup                         Remove extra duplicates (RT and window)
--keepRTdup                        Keep RT duplicates (if --extraDup is specified)
--keepDups                         Keep all duplicated reads
--keepBlackList                    Keep reads in blacklist regions
--distDup                INTEGER   Genomic distance to consider a read as a window duplicate (for scChIP only)

MATRICES
--minReadsPerCellmatrix  INTEGER   Cells having less than this number are removed from final matrices
--binSize                INTEGER [50000]  Size of bins to create matrices

PEAK CALLING
--peakCalling                                  Run bulk peak calling analysis
--macs2Opts              STRING                MACS2 parameters
--peakDist               INTEGER               Maximum distance between peaks to be merged
--tssWindow              INTEGER               Distance (upstream/downstream) to transcription start point to consider

SKIP OPTIONS:
--skipBigWig                 Disable BigWig
--skipMultiQC                Disable MultiQC
--skipSoftVersions           Disable Soft Versions

OTHER OPTIONS:
--metadata      PATH     Specify a custom metadata file for MultiQC
--multiqcConfig PATH     Specify a custom config file for MultiQC
--name          STRING   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

OUTPUTs:
--outDir             PATH     The output directory where the results will be saved
--saveIntermediates           Save intermediates files
--cleanup            STRING [none, auto, success]  Cleaning strategy of the work/ directory

======================================================
Available Profiles
-profile test                        Run the test dataset
-profile conda                       Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile multiconda                  Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile path                        Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
-profile multipath                   Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
-profile docker                      Use the Docker images for each process
-profile singularity                 Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
-profile cluster                     Run the workflow on the cluster, instead of locally
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
