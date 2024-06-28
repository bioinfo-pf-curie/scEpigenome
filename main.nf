#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2022
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
========================================================================================
                         DSL2 Template
========================================================================================
Analysis Pipeline DSL2 template.
https://patorjk.com/software/taag/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Initialize lintedParams and paramsWithUsage
NFTools.welcome(workflow, params)

// Use lintedParams as default params object
paramsWithUsage = NFTools.readParamsFromJsonSettings("${projectDir}/parameters.settings.json")
params.putAll(NFTools.lint(params, paramsWithUsage))

// Run name
customRunName = NFTools.checkRunName(workflow.runName, params.protocol)

// Custom functions/variables
mqcReport = []
include {checkAlignmentPercent} from './lib/functions'

/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/

// Genome-based variables
if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory. genome params : '${params.genome}'"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$projectDir/docs/output.md")
chOutputDocsImages = file("$projectDir/docs/images/", checkIfExists: true)

// Initialize variable from the genome.conf file
params.starIndex = NFTools.getGenomeAttribute(params, 'starIndex')
params.blackList = NFTools.getGenomeAttribute(params, 'blackList')
params.gtf = NFTools.getGenomeAttribute(params, 'gtf')
params.effGenomeSize = NFTools.getGenomeAttribute(params, 'effGenomeSize')
params.fasta = NFTools.getGenomeAttribute(params, 'fasta')
params.geneBed = NFTools.getGenomeAttribute(params, 'geneBed')

/*
==========================
 VALIDATE INPUTS
==========================
*/

if ((params.reads && params.samplePlan) || (params.readPaths && params.samplePlan)){
  exit 1, "Input reads must be defined using either '--reads' or '--samplePlan' parameter. Please choose one way"
}

/*
==========================
 BUILD CHANNELS
==========================
*/

chStarIndex     = params.starIndex     ? Channel.fromPath(params.starIndex, checkIfExists: true).collect()  : Channel.empty()
chBlackList     = params.blackList     ? Channel.fromPath(params.blackList, checkIfExists: true).collect()  : Channel.empty()
chGtf           = params.gtf           ? Channel.fromPath(params.gtf, checkIfExists: true).collect()        : Channel.empty()
chGeneBed       = params.geneBed       ? Channel.fromPath(params.geneBed, checkIfExists: true).collect()    : channel.empty()
chFasta         = params.fasta         ? Channel.fromPath(params.fasta, checkIfExists: true).collect()      : Channel.empty()
chEffGenomeSize = params.effGenomeSize ? Channel.of(params.effGenomeSize)                                   : Channel.value([])
chGeneBed       = params.geneBed       ? Channel.fromPath(params.geneBed, checkIfExists: true).collect()    : channel.empty()
chMetadata      = params.metadata      ? Channel.fromPath(params.metadata, checkIfExists: true).collect()   : channel.empty()
chBinSize       = Channel.from(params.binSize).splitCsv().flatten().toInteger()

chBinSize.view()

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline' : workflow.manifest.name ?: null,
  'Protocol': params.protocol,
  'Version': workflow.manifest.version ?: null,
  'DOI': workflow.manifest.doi ?: null,
  'Run Name': customRunName,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'Genome' : params.genome,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Profile' : workflow.profile,
  'OutDir' : params.outDir,
  'WorkDir': workflow.workDir,
  'CommandLine': workflow.commandLine
].findAll{ it.value != null }

workflowSummaryCh = NFTools.summarize(summary, workflow, params)

/*
==============================
  LOAD INPUT DATA
==============================
*/

// Load raw reads
chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.protocol, params)

// Make samplePlan if not available
sPlanCh = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.protocol) 

/*
==================================
           INCLUDE
==================================
*/ 

include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { intersectBed as rmBlackList } from './nf-modules/common/process/bedtools/intersectBed'
include { samtoolsIndex } from './nf-modules/common/process/samtools/samtoolsIndex'
include { samtoolsSort as samtoolsSortByName } from './nf-modules/common/process/samtools/samtoolsSort'
include { bamToFrag } from './nf-modules/local/process/bamToFrag'
include { multiqc } from './nf-modules/local/process/multiqc'

include { scchipFlow } from './nf-modules/local/subworkflow/scchipFlow'
include { sccuttagIndropFlow } from './nf-modules/local/subworkflow/sccuttagIndropFlow' 
include { sccuttag10XFlow } from './nf-modules/local/subworkflow/sccuttag10XFlow' 
include { bigwigFlow } from './nf-modules/local/subworkflow/bigwigFlow'
include { countMatricesFlow } from './nf-modules/local/subworkflow/countMatricesFlow'
include { peakCallingFlow } from './nf-modules/local/subworkflow/peakCallingFlow'

/*
=====================================
            WORKFLOW 
=====================================
*/

workflow {

  main:
    chVersions = Channel.empty()

    // subroutines
    outputDocumentation(
      chOutputDocs,
      chOutputDocsImages
    )

    // subworkflow: From fastq to BAM
    if (params.protocol=='sccuttag_10X'){
      sccuttag10XFlow(
        chRawReads,
        chStarIndex
      )
      chBam = sccuttag10XFlow.out.bam
      chBcLogs = sccuttag10XFlow.out.bcLogs
      chStarLogs = sccuttag10XFlow.out.starLogs
      chMdLogs = sccuttag10XFlow.out.mdLogs
      chStats = sccuttag10XFlow.out.stats
      chBarcodes = sccuttag10XFlow.out.barcodes
      chBarcodesCounts = sccuttag10XFlow.out.counts
      chWhist = sccuttag10XFlow.out.whist
      chVersions = sccuttag10XFlow.out.versions
    }

    if (params.protocol=='sccuttag_indrop'){
      sccuttagIndropFlow(
        chRawReads,
        chStarIndex
      )
      chBam = sccuttagIndropFlow.out.bam
      chBcLogs = sccuttagIndropFlow.out.bcLogs
      chStarLogs = sccuttagIndropFlow.out.starLogs
      chMdLogs = sccuttagIndropFlow.out.mdLogs
      chStats = sccuttagIndropFlow.out.stats
      chBarcodes = sccuttagIndropFlow.out.barcodes
      chBarcodesCounts = sccuttagIndropFlow.out.counts
      chWhist =	 sccuttagIndropFlow.out.whist
      chVersions = sccuttagIndropFlow.out.versions
    }

    if (params.protocol=='scchip_indrop'){
      scchipFlow(
        chRawReads,
        chStarIndex,
      )
      chBam = scchipFlow.out.bam
      chBcLogs = scchipFlow.out.bcLogs
      chStarLogs = scchipFlow.out.starLogs
      chMdLogs = scchipFlow.out.mdLogs
      chStats = scchipFlow.out.stats
      chBarcodes = scchipFlow.out.barcodes
      chBarcodesCounts = scchipFlow.out.counts
      chWhist =	 scchipFlow.out.whist
      chVersions = scchipFlow.out.versions
    }

    // process: remove blacklist regions
    rmBlackList(
      chBam,
      chBlackList
    )
    chVersions = chVersions.mix(rmBlackList.out.versions)

    samtoolsIndex(
      rmBlackList.out.bam
    )
    chBamBai = rmBlackList.out.bam.join(samtoolsIndex.out.bai)

    // subworkflow: generate bigwig files
    chBigWig = Channel.empty()
    if (!params.skipBigwig){
      bigwigFlow(
        chBamBai,
        chEffGenomeSize,
        chGeneBed
      )
      chVersions = chVersions.mix(bigwigFlow.out.versions)
      chBigWig = bigwigFlow.out.mqc.collect()
    }
 
    //subworflow: get count matrices
    countMatricesFlow(
      chBamBai,
      chBarcodes,
      chBinSize,
      chGeneBed
    )
    chVersions = chVersions.mix(countMatricesFlow.out.versions)

    //process: generate fragment data
    samtoolsSortByName(
      chBamBai.map{meta, bam, bai -> [meta, bam]}
    )
    chVersions = chVersions.mix(samtoolsSortByName.out.versions)

    bamToFrag(
      samtoolsSortByName.out.bam 
    )
    chVersions = chVersions.mix(bamToFrag.out.versions)

    // subWorkflow: Pseudo Bulk peak calling
    chPeaksCountsMqc = Channel.empty()
    chPeaksFrip = Channel.empty()
    chPeaksQC = Channel.empty()
    if (params.peakCalling){
      peakCallingFlow(
        chBamBai,
        chEffGenomeSize,
        chGtf,
        chFasta
      )
      chPeaksCountsMqc = peakCallingFlow.out.peaksCountsMqc.collect().ifEmpty([]),
      chPeaksFrip = peakCallingFlow.out.frip.collect().ifEmpty([]),
      chPeaksQC = peakCallingFlow.out.qc.collect().ifEmpty([]),
      chVersions = chVersions.mix(peakCallingFlow.out.versions)
    }

    //*******************************************
    // MULTIQC
    if (!params.skipMultiQC){
       getSoftwareVersions(
          chVersions.unique().collectFile()
       )

       chWarn = Channel.empty()
       
       multiqc(
         customRunName,
         sPlanCh.collect(),
         chMetadata.ifEmpty([]),
         chMultiqcConfig.ifEmpty([]),
         chBcLogs.collect().ifEmpty([]),
    //     joinBcIndexesLogsCollected.collect().ifEmpty([])
         chBarcodesCounts.map{meta,counts->counts}.collect().ifEmpty([]),
         chWhist.collect().ifEmpty([]),
         chStarLogs.collect().ifEmpty([]),
	 chStats.map{it->[it[1]]}.collect().ifEmpty([]),
	 chMdLogs.collect().ifEmpty([]),
         chPeaksCountsMqc.collect().ifEmpty([]),
         chPeaksFrip.collect().ifEmpty([]),
         chPeaksQC.collect().ifEmpty([]),
         //chPeaksSizesMqc.collect().ifEmpty([]),
         chBigWig.collect().ifEmpty([]),
         getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
         workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
         chWarn.collect().ifEmpty([])
       )
    }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
