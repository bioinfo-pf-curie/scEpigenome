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

chStarIndex          = params.starIndex                ? Channel.fromPath(params.starIndex, checkIfExists: true).collect()         : Channel.empty()
chBlackList          = params.blackList                ? Channel.fromPath(params.blackList, checkIfExists: true).collect()         : Channel.empty()
chGtf                = params.gtf                      ? Channel.fromPath(params.gtf, checkIfExists: true).collect()               : Channel.empty()
chBowtie2_10Xbc      = params.barcodes10X_bwt2         ? Channel.fromPath(params.barcodes10X_bwt2, checkIfExists: true).collect()  : Channel.empty()
chBinSize            = Channel.from(params.binSize).splitCsv().flatten().toInteger()

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { metadataCh }
}

//------- Custom barcode indexes--------
//--------------------------------------
for ( idx in params.barcodesIndrop.keySet() ){
  if ( params.barcodesIndrop[ idx ].bwt2 ){
    lastPath = params.barcodesIndrop[ idx ].bwt2.lastIndexOf(File.separator)
    bt2Dir = params.barcodesIndrop[ idx ].bwt2.substring(0,lastPath+1)
    bt2Base = params.barcodesIndrop[ idx ].bwt2.substring(lastPath+1)
    params.barcodesIndrop[ idx ].base = bt2Base
    params.barcodesIndrop[ idx ].dir = bt2Dir
  }
}

Channel
   .from(params.barcodesIndrop)
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
chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.libType, params)

// Make samplePlan if not available
// R3 added :
sPlanCh = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd) // FUTURE : modif params.singleEnd by libType !!!!!!!!!!!!!!!!!!!!!!!!!!!

/*
==================================
           INCLUDE
==================================
*/ 

// Workflows


// countMatricesPerPeak to be create

include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { scchip } from './nf-modules/local/subworkflow/scchip'
include { sccuttag_indrop } from './nf-modules/local/subworkflow/sccuttag_indrop' 
include { sccuttag_10X } from './nf-modules/local/subworkflow/sccuttag_10X' 
include { starAlign } from './nf-modules/common/process/star/starAlign'


/*
=====================================
            WORKFLOW 
=====================================
*/

workflow {

  main:

    // subroutines
    outputDocumentation(
      outputDocsCh,
      outputDocsImagesCh
    )

    chRawReads
        .collect()
        .filter('/data/users/lhadjabe/scChIP/testData_scEpigenome/testdata_scCut/scCUT10X/H01_1/T_AE4734_hu_nuclei_m06y22_H3K4me1_S1_R1_001.fastq.gz')
        .set{test}
      
      test.view() 

    if (params.protocol=='sccuttag_10X'){

      sccuttag_10X(
        chRawReads,
        workflowSummaryCh,
        multiqcConfigCh,
        metadataCh,
        sPlanCh,
        customRunName,
        chBowtie2_10Xbc,
        chStarIndex,
        chBlackList,
        chGtf,
        chBinSize
      )
      chBam = sccuttag_10X.out.bam
      chBai = sccuttag_10X.out.bai
      chBw = sccuttag_10X.out.bigwig
      chTSSmat  = sccuttag_10X.out.matrixTSS
      chBinmat = sccuttag_10X.out.matrixBin 
      chMQChtml = sccuttag_10X.out.mqcreport 
    }

    if (params.protocol=='sccuttag_indrop'){
      // want to select only id, R2 == BC
      chRawReads
        .collect() {item -> [item[0], item[1][1]] }
        .set{chBarcodeRead}
      // want to select only id, R1 and R3 == DNA
      chRawReads
        .collect() {item -> [item[0], [item[1][0], item[1][2]]] }
        .set{chDNAreads}
      sccuttag_indrop(
        chBarcodeRead,
        chDNAreads,
        workflowSummaryCh,
        multiqcConfigCh,
        metadataCh,
        sPlanCh,
        customRunName,
        chIndexBwt2,
        chStarIndex,
        chBlackList,
        chGtf,
        chBinSize
      )
      chBam = sccuttag_indrop.out.bam
      chBai = sccuttag_indrop.out.bai
      chBw = sccuttag_indrop.out.bigwig
      chTSSmat  = sccuttag_indrop.out.matrixTSS
      chBinmat = sccuttag_indrop.out.matrixBin 
      chMQChtml = sccuttag_indrop.out.mqcreport 
    }

    if (params.protocol=='scchip_indrop'){
      chRawReads
        .collect() {item -> [item[0], item[1][1]] }
        .set{chBarcodeRead}
      scchip(
        chRawReads,
        chBarcodeRead,
        workflowSummaryCh,
        multiqcConfigCh,
        metadataCh,
        sPlanCh,
        customRunName,
        chIndexBwt2,
        chStarIndex,
        chBlackList,
        chGtf,
        chBinSize
      )
      chBam = scchip.out.bam
      chBai = scchip.out.bai
      chBw = scchip.out.bigwig
      chTSSmat  = scchip.out.matrixTSS
      chBinmat = scchip.out.matrixBin 
      chMQChtml = scchip.out.mqcreport 
    }
    
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
