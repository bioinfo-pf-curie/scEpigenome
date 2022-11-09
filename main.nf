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
customRunName = NFTools.checkRunName(workflow.runName, params.name)

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
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Initialize variable from the genome.conf file
//params.bowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2')

// Stage config files
multiqcConfigCh = Channel.fromPath(params.multiqcConfig)
outputDocsCh = Channel.fromPath("$projectDir/docs/output.md")
outputDocsImagesCh = file("$projectDir/docs/images/", checkIfExists: true)

// Initialize variable from the genome.conf file
params.starIndex = NFTools.getGenomeAttribute(params, 'starIndex')
params.blackList = NFTools.getGenomeAttribute(params, 'blackList')
params.gtf = NFTools.getGenomeAttribute(params, 'gtf')

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

chStarIndex    = params.starIndex   ? Channel.fromPath(params.starIndex, checkIfExists: true).collect()         : Channel.empty()
chBlackList    = params.blackList   ? Channel.fromPath(params.blackList, checkIfExists: true).collect()         : Channel.empty()
chGtf          = params.gtf         ? Channel.fromPath(params.gtf, checkIfExists: true).collect()               : Channel.empty()
chBinSize      = Channel.from(params.binSize).splitCsv().flatten().toInteger()

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { metadataCh }
}

//------- Custom barcode indexes--------
//--------------------------------------
for ( idx in params.barcodes.keySet() ){
  if ( params.barcodes[ idx ].bwt2 ){
    lastPath = params.barcodes[ idx ].bwt2.lastIndexOf(File.separator)
    bt2Dir = params.barcodes[ idx ].bwt2.substring(0,lastPath+1)
    bt2Base = params.barcodes[ idx ].bwt2.substring(lastPath+1)
    params.barcodes[ idx ].base = bt2Base
    params.barcodes[ idx ].dir = bt2Dir
  }
}

Channel
   .from(params.barcodes)
   .flatMap()
   .map { it -> [ it.key, file(it.value['dir']) ] }
   .ifEmpty { exit 1, "Bowtie2 index not found" }
   .set { chIndexBwt2 } 

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline Release': workflow.revision ?: null,
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
// R3 added :
chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)

// Make samplePlan if not available
// R3 added :
sPlanCh = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd)

/*
==================================
           INCLUDE
==================================
*/ 

// Workflows

// Processes
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { starAlign } from './nf-modules/common/process/star/starAlign'
//local
include { multiqc } from './nf-modules/local/process/multiqc'
include { bcAlign } from './nf-modules/local/process/bcAlign'
include { bcSubset } from './nf-modules/local/process/bcSubset'
include { bcTrim } from './nf-modules/local/process/bcTrim'
include { addBarcodeTag } from './nf-modules/local/process/addBarcodeTag'
    // remove duplicates
include { removePCRdup } from './nf-modules/local/process/removePCRdup'
include { removeRTdup } from './nf-modules/local/process/removeRTdup'
include { removeWindoWdup } from './nf-modules/local/process/removeWindoWdup'
    //------
  // blackRegions
include { gtfToTSSBed } from './nf-modules/local/process/gtfToTSSBed'
include { removeBlackRegions } from './nf-modules/local/process/removeBlackRegions'
  //--------
include { countSummary } from './nf-modules/local/process/countSummary'
include { distribUMIs } from './nf-modules/local/process/distribUMIs'
include { bigwig } from './nf-modules/local/process/bigwig'
include { bamToFrag } from './nf-modules/local/process/bamToFrag'
//subworkflow
include { countMatricesPerBin } from './nf-modules/local/subworkflow/countMatricesPerBin' 

/*
=====================================
            WORKFLOW 
=====================================
*/

workflow {
  chVersions = Channel.empty()

  main:
    // Init Channels
    chAlignedLogs = Channel.empty()
    outputDocsImagesCh = Channel.empty()

    // subroutines
    outputDocumentation(
      outputDocsCh,
      outputDocsImagesCh
    )

    // 1) Barcode alignement and extrcation part
    bcAlign(
      chRawReads.combine(chIndexBwt2)
    )
    chReadsMatchingIndex = bcAlign.out.results
    chIndexCount = bcAlign.out.counts
    chIndexBowtie2Logs = bcAlign.out.logs
    chVersions = chVersions.mix(bcAlign.out.versions)

    bcSubset(
      chReadsMatchingIndex.groupTuple(),
      chIndexCount.groupTuple()
    )
    chReadBcNames = bcSubset.out.results
    chBowtie2Logs = bcSubset.out.logs

    // 2) DNA alignment part
    bcTrim(
      chRawReads
    )
    chTrimmedReads = bcTrim.out.reads
    chTrimmedReadsLogs = bcTrim.out.logs
    chVersions = chVersions.mix(bcTrim.out.versions)

    chRawReads
      .join(chTrimmedReads)
      .map{ it -> [it[0], it[1][0], it[2]]}
      .view()
      .set{chReads}
    
    starAlign(
      //inputs
      chReads,
      chStarIndex
      //parameters to add in conf/modules
    )
    //outputs
    chAlignedBam = starAlign.out.bam
    chAlignedLogs = starAlign.out.logs
    chVersions = chVersions.mix(starAlign.out.versions)

    // Add barcode info into dna info
    addBarcodeTag(
      chAlignedBam.join(chReadBcNames)
    )
    chTaggedBam=addBarcodeTag.out.bam

    removePCRdup(
      //inputs
      chTaggedBam
    )
    //outputs
    chRemovePCRdupBam = removePCRdup.out.bam
    chRemovePCRdupSam = removePCRdup.out.sam
    chPCRdupCount = removePCRdup.out.count
    chR1unmappedR2Count = removePCRdup.out.countR1unmapped
    chRemovePcrRtDup_Log = removePCRdup.out.logs

    removeRTdup(
      //inputs
      chTaggedBam,
      chRemovePCRdupBam,
      chRemovePCRdupSam
    )
    //outputs
    chRemovePcrRtDup = removeRTdup.out.bam
    chRTdupCount = removeRTdup.out.logs
    
    removeWindoWdup(
      //inputs
      chRemovePcrRtDup
    )
    //outputs
    chRemoveBlackReg = removeWindoWdup.out.bam
    chRemoveDupLog = removeWindoWdup.out.logs

    removeBlackRegions(
      //inputs
      chRemoveBlackReg,
      chBlackList.collect()
    )
    chVersions = chVersions.mix(removeBlackRegions.out.versions)
    chNoDupBam = removeBlackRegions.out.bam
    chNoDupBai = removeBlackRegions.out.bai
    chBlackRegionsCount = removeBlackRegions.out.sam

    countSummary(
      //inputs
      chRemovePcrRtDup_Log,
      chPCRdupCount,
      chRTdupCount,
      chR1unmappedR2Count,
      chBlackRegionsCount
    )
    chDedupCountSummary = countSummary.out.logs
    chfinalBClist = countSummary.out.result
    chfinalBCcounts = countSummary.out.count

    countMatricesPerBin(
      chfinalBClist,
      chNoDupBam,
      chNoDupBai,
      chBinSize
    )
    chMatrices=countMatricesPerBin.out.matrix

    distribUMIs(
      //inputs
      chfinalBClist
    )
    chMqcDistribUMI = distribUMIs.out.mqc
    chPdfDist = distribUMIs.out.pdf
    chVersions = chVersions.mix(removeBlackRegions.out.versions)

    bigwig(
      //inputs
      chNoDupBam.join(chNoDupBai),
      chBlackList.collect()
    )
    //outputs
    chBigWig = bigwig.out.bigwig
    chBigWigLogs = bigwig.out.logs
    chVersions = chVersions.mix(bigwig.out.versions)

    bamToFrag(
      //inputs
      chNoDupBam.join(chNoDupBai)
    )
    //outputs
    chFragmentFiles = bamToFrag.out.gz

    gtfToTSSBed(
      //inputs
      chGtf
    )
    //outputs
    chGtfToTSSBed= gtfToTSSBed.out.bed

    //*******************************************
    // MULTIQC
  
    // Warnings that will be printed in the mqc report
    warnCh = Channel.empty()

    if (!params.skipMultiQC){

      getSoftwareVersions(
        chVersions.unique().collectFile()
      )

      multiqc(
        customRunName,
        sPlanCh.collect(),
        metadataCh.ifEmpty([]),
        multiqcConfigCh.ifEmpty([]),
        getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
        workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
        warnCh.collect().ifEmpty([])
      )
      mqcReport = multiqc.out.report.toList()
    }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
