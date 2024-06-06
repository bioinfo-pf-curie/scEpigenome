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
multiqcConfigCh = Channel.fromPath(params.multiqcConfig)
outputDocsCh = Channel.fromPath("$projectDir/docs/output.md")
outputDocsImagesCh = file("$projectDir/docs/images/", checkIfExists: true)

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
chFasta         = params.fasta         ? Channel.fromPath(params.fasta, checkIfExists: true).collect()      : Channel.empty()
chEffGenomeSize = params.effGenomeSize ? Channel.of(params.effGenomeSize)                                   : Channel.value([])
chGeneBed       = params.geneBed       ? Channel.fromPath(params.geneBed, checkIfExists: true).collect()    : channel.empty()
chMetadata      = params.metadata      ? Channel.fromPath(params.metadata, checkIfExists: true).collect()   : channel.empty()
chBinsize       = Channel.from(params.binSize).splitCsv().flatten().toInteger()

//------- Custom barcode indexes--------
//--------------------------------------
//for ( idx in params.barcodesIndrop.keySet() ){
//  if ( params.barcodesIndrop[ idx ].bwt2 ){
//    lastPath = params.barcodesIndrop[ idx ].bwt2.lastIndexOf(File.separator)
//    bt2Dir = params.barcodesIndrop[ idx ].bwt2.substring(0,lastPath+1)
//    bt2Base = params.barcodesIndrop[ idx ].bwt2.substring(lastPath+1)
//    params.barcodesIndrop[ idx ].base = bt2Base
//    params.barcodesIndrop[ idx ].dir = bt2Dir
//  }
//}

//Channel
//   .from(params.barcodesIndrop)
//   .flatMap()
//   .map { it -> [ it.key, file(it.value['dir']) ] } 
//    // [indexB, /data/annotations/pipelines/tools/scRNA_LBC_bowtie2_indexes]
//    // [indexC, /data/annotations/pipelines/tools/scRNA_LBC_bowtie2_indexes]
//    // [indexD, /data/annotations/pipelines/tools/scRNA_LBC_bowtie2_indexes]
//   .ifEmpty { exit 1, "Bowtie2 index not found" }
//   .set { chIndexBwt2 } 

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline' : workflow.manifest.name ?: null,
  'Version': workflow.manifest.version ?: null,
  'DOI': workflow.manifest.doi ?: null,
  'Run Name': customRunName,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'Genome' : params.genome,
  'Protocol': params.protocol,
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

// Workflows and processes

include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { intersectBed as rmBlackList } from './nf-modules/common/process/bedtools/intersectBed'
include { samtoolsIndex } from './nf-modules/common/process/samtools/samtoolsIndex'
include { bamToFrag } from './nf-modules/local/process/bamToFrag'

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
      outputDocsCh,
      outputDocsImagesCh
    )

    // subworkflow: From fastq to BAM
    if (params.protocol=='sccuttag_10X'){
      sccuttag10XFlow(
        chRawReads,
        chStarIndex
      )
      chBam = sccuttag10XFlow.out.bam
      chBarcodes = sccuttag10XFlow.out.barcodes
      chVersions = sccuttag10XFlow.out.versions
    }

    if (params.protocol=='sccuttag_indrop'){
      sccuttagIndropFlow(
        chRawReads,
        chStarIndex
      )
      chBam = sccuttagIndropFlow.out.bam
      chBarcodes = sccuttagIndropFlow.out.barcodes
      chVersions = sccuttagIndropFlow.out.versions
    }

    if (params.protocol=='scchip_indrop'){
      scchipFlow(
        chRawReads,
        chStarIndex,
      )
      chBam = scchipFlow.out.bam
      chBarcodes = scchipFlow.out.barcodes
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
    bigwigFlow(
      chBamBai,
      chEffGenomeSize,
      chGeneBed
    )
    chVersions = chVersions.mix(bigwigFlow.out.versions)

    //subworflow: get count matrices
    // To modify - based on XB tag
    //countMatricesFlow(
    //  chBamBai,
    //  chBarcodes,
    //  chBinSize,
    //  chGtf
    //)
    //chVersions = chVersions.mix(countMatricesFlow.out.versions)

    //process: generate fragment data
    // To modify - based on XB tag
    //bamToFrag(
    //  chBamBai
    //)
    //chVersions = chVersions.mix(bamToFrag.out.versions)

    // subWorkflow: Pseudo Bulk peak calling
    peakCallingFlow(
       chBamBai,
       chEffGenomeSize,
       chGtf,
       chFasta
    )
    chVersions = chVersions.mix(peakCallingFlow.out.versions)
    //  peaksPseudoBulkBed = peaksPseudoBulk.out.mergedPeaks
    //  chPeaksCountsMqc = peaksPseudoBulk.out.peaksCountsMqc
    //  chPeaksSizesMqc = peaksPseudoBulk.out.peaksSizesMqc
    //  chFripResults = peaksPseudoBulk.out.fripResults
    //  chPeaksQCMqc = peaksPseudoBulk.out.peaksQCMqc

    //*******************************************
    // MULTIQC
    if (!params.skipMultiQC){
       getSoftwareVersions(
          chVersions.unique().collectFile()
       )

    //   multiqc(
    //     customRunName,
    //     sPlanCh.collect(),
    //     metadataCh.ifEmpty([]),
    //     multiqcConfigCh.ifEmpty([]),
    //     getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
    //     workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
    //     //warnCh.collect().ifEmpty([]),
    //     chAlignedLogs.collect().ifEmpty([]), //star = *Log.final.out
    //     // bcAlign:
    //     chIndexBowtie2Logs.collect().ifEmpty([]),//index/${sample}_Bowtie2.log
    //     // bcSubset:
    //     joinBcIndexesLogsCollected.collect().ifEmpty([]),//bowtie2/${sample}_bowtie2.log
    //     // countSummary:
    //     chDedupCountSummary.collect().ifEmpty([]),//allDup/${sample}_allDup.log
    //     // countSummary:
    //     chfinalBClistCollected.collect().ifEmpty([]),//cellThresholds/${sample}_rmDup.txt
    //     // removeWindowDup:
    //     chRemoveDupLog.collect().ifEmpty([]),//removeWindowDup/${sample}_removeWindowDup.log (#Number of duplicates: nnnn)
    //     //distribUMIs
    //     chMqcDistribUMI.collect().ifEmpty([]), //pour config graph
    //     // macs2
    //     chPeaksCountsMqc.collect().ifEmpty([]),
    //     //peakQC
    //     chFripResults.collect().ifEmpty([]),
    //     chPeaksQCMqc.collect().ifEmpty([]),
    //     // from deeptools matrix
    //     chDeeptoolsProfileMqc.collect().ifEmpty([]),
    //     // macs2
    //     chPeaksSizesMqc.collect().ifEmpty([])
    //   )
    //   chMqcReport = multiqc.out.report.toList()
    }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
