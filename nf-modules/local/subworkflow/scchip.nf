//common
include { getSoftwareVersions } from '../../common/process/utils/getSoftwareVersions'
include { starAlign } from '../../common/process/star/starAlign'
include { deeptoolsBamCoverage } from '../../common/process/deeptools/deeptoolsBamCoverage'
include { deeptoolsComputeMatrix } from '../../common/process/deeptools/deeptoolsComputeMatrix'

//include { bigwig } from '../../local/process/bigwig' // move to common one condition a mettre dans modules pour les args
// add preseq
//local
include { multiqc } from '../../local/process/multiqc'
include { bcAlign } from '../../local/process/bcAlign'
include { joinBcIndexes } from '../../local/process/joinBcIndexes'
include { bcTrim } from '../../local/process/bcTrim'
include { addFlags } from '../../local/process/addFlags'
  // remove duplicates
include { removePCRdup } from '../../local/process/removePCRdup' 
include { removeRTdup } from '../../local/process/removeRTdup'
include { removeWindowDup } from '../../local/process/removeWindowDup'
  // blackRegions
include { removeBlackRegions } from '../../local/process/removeBlackRegions'
  //--------
include { countSummary } from '../../local/process/countSummary' // empty channels pour éviter bug car pas de RT ni Window?
include { distribUMIs } from '../../local/process/distribUMIs'
include { bamToFrag } from '../../local/process/bamToFrag'
include { countMatricesPerBin } from '../../local/process/countMatricesPerBin'
//subworkflow
include { countMatricesPerTSSFlow } from '../../local/subworkflow/countMatricesPerTSSFlow' 
include { peaksPseudoBulk } from '../../local/subworkflow/peaksPseudoBulk' 


workflow scchip {

  take:
  reads
  barcodeRead
  workflowSummaryCh
  multiqcConfigCh
  metadataCh
  sPlanCh
  customRunName
  bowtie2Index
  starIndex
  blackList
  gtf
  fasta
  binsize
  effGenomeSize ////////////////////////////////////////////////////////////////
  geneBed ////////////////////////////////////////////////////////////////

  main:
    // Init Channels
    // channels never filled
    chStarGtf  = Channel.value([])
    // channels filled
    chRemoveDupLog = Channel.empty()
    chBigWig= Channel.empty()
    chAlignedLogs = Channel.empty()
    chRemoveRtSummary = Channel.empty()
    warnCh = Channel.empty()
    chVersions = Channel.empty()

    // if BigWig
    chDeeptoolsProfileMqc = Channel.empty()

    // 1) Barcode alignement and extrcation part
    bcAlign(
      barcodeRead.combine(bowtie2Index)
    )
    chReadsMatchingIndex = bcAlign.out.results
    chIndexCount = bcAlign.out.counts
    chIndexBowtie2Logs = bcAlign.out.logs
    chVersions = chVersions.mix(bcAlign.out.versions)

    joinBcIndexes(
      chReadsMatchingIndex.groupTuple(),
      chIndexCount.groupTuple()
    )
    chReadBcNames = joinBcIndexes.out.results
    joinBcIndexesLogs = joinBcIndexes.out.logs

    // 2) DNA alignment part
    bcTrim(
      reads
    )
    chTrimmedReads = bcTrim.out.reads
    chTrimmedReadsLogs = bcTrim.out.logs
    chVersions = chVersions.mix(bcTrim.out.versions)

    reads
      .join(chTrimmedReads)
      .map{ it -> [it[0], [it[1][0], it[2]]]}
      .set{chReads}

    starAlign(
      //inputs
      chReads,
      starIndex,
      chStarGtf
      //parameters to add in conf/modules
    )
    //outputs
    chAlignedBam = starAlign.out.bam
    chAlignedLogs = starAlign.out.logs
    chVersions = chVersions.mix(starAlign.out.versions)

    // Add barcode info into dna info
    addFlags(
      chAlignedBam.join(chReadBcNames)
    )
    chTaggedBam=addFlags.out.bam

    removePCRdup(
      //inputs
      chTaggedBam
    )
    //outputs
    chRemovePCRdupBam = removePCRdup.out.bam
    chRemovePCRdupSam = removePCRdup.out.sam
    chRemovePCRdupSummary = removePCRdup.out.count
    chR1unmappedR2Summary = removePCRdup.out.countR1unmapped
    chRemovePcrBamSummary = removePCRdup.out.bamLogs

    removeRTdup(
      //inputs
      chTaggedBam,
      chRemovePCRdupBam,
      chRemovePCRdupSam
    )
    //outputs
    chRemovePcrRtBam = removeRTdup.out.bam
    chRemoveRtSummary = removeRTdup.out.logs
    
    removeWindowDup(
      //inputs
      chRemovePcrRtBam
    )
    //outputs
    chRemoveBlackReg = removeWindowDup.out.bam
    chRemoveDupLog = removeWindowDup.out.logs

    //Chout.map{meta, table -> [meta, table, []]}

    removeBlackRegions(
      //inputs
      chRemoveBlackReg,
      blackList.collect()
    )
    chVersions = chVersions.mix(removeBlackRegions.out.versions)
    chNoDupBam = removeBlackRegions.out.bam
    chNoDupBai = removeBlackRegions.out.bai
    chfinalBClist = removeBlackRegions.out.list

    peaksPseudoBulk( ////////////////////////////////////////////////////////////////
      chNoDupBam,
      chNoDupBai,
      effGenomeSize,
      gtf,
      fasta
    )
    peaksPseudoBulkBed = peaksPseudoBulk.out.mergedPeaks
    chPeaksCountsMqc = peaksPseudoBulk.out.peaksCountsMqc
    chFripResults = peaksPseudoBulk.out.fripResults
    chPeaksQCMqc = peaksPseudoBulk.out.peaksQCMqc
    chVersions = chVersions.mix(peaksPseudoBulk.out.versions)

    countSummary(
      //inputs
      chRemovePCRdupSummary, // pcr
      chRemovePcrBamSummary, // pcr
      chR1unmappedR2Summary, // pcr
      chRemoveRtSummary // faire des empty channels 
    )
    chDedupCountSummary = countSummary.out.logs

    // Subworkflow
    countMatricesPerBin( ////////////////////////////////////////////////////////////////
      binsize,
      chNoDupBam.join(chNoDupBai),
      chfinalBClist
    )
    chBinMatrices=countMatricesPerBin.out.matrix
    chVersions = chVersions.mix(countMatricesPerBin.out.versions)

    // Subworkflow
    countMatricesPerTSSFlow( ////////////////////////////////////////////////////////////////
      chNoDupBam.join(chNoDupBai),
      chfinalBClist,
      gtf
    )
    chTssMatrices=countMatricesPerTSSFlow.out.matrix
    chVersions = chVersions.mix(countMatricesPerTSSFlow.out.versions)

    distribUMIs(
      //inputs
      chfinalBClist
    )
    chMqcDistribUMI = distribUMIs.out.mqc
    chPdfDist = distribUMIs.out.pdf
    chVersions = chVersions.mix(removeBlackRegions.out.versions)

    if (!params.skipBigWig){
      deeptoolsBamCoverage(
        //inputs
        chNoDupBam.join(chNoDupBai),
        blackList.collect(),
        effGenomeSize
      )
      //outputs
      chBigWig = deeptoolsBamCoverage.out.bigwig
      chVersions = chVersions.mix(deeptoolsBamCoverage.out.versions)
      
      deeptoolsComputeMatrix( ////////////////////////////////////////////////////////////////
        chBigWig,
        geneBed.collect()
      )
      chDeeptoolsProfileMqc = deeptoolsComputeMatrix.out.mqc
      chVersions = chVersions.mix(deeptoolsComputeMatrix.out.versions) 
    }

    bamToFrag(
      //inputs
      chNoDupBam.join(chNoDupBai)
    )
    //outputs
    chFragmentFiles = bamToFrag.out.gz

    // delete $meta for mqc input
    chfinalBClist
    .map{it -> it[1]}
    .set{chfinalBClistCollected}

    //*******************************************
    // MULTIQC
  
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
        warnCh.collect().ifEmpty([]),
        chAlignedLogs.collect().ifEmpty([]), //star
        // bcAlign:
        chIndexBowtie2Logs.collect().ifEmpty([]),//index/${sample}_indexBBowtie2.log
        // joinBcIndexes:
        joinBcIndexesLogs.collect().ifEmpty([]),//bowtie2/${sample}_bowtie2.log
        // countSummary:
        chDedupCountSummary.collect().ifEmpty([]),//removeRtPcr/${sample}_removePcrRtDup.log
        // countSummary:
        chfinalBClistCollected.collect().ifEmpty([]),//cellThresholds/${sample}_rmDup.txt
        // removeWindowDup:
        chRemoveDupLog.collect().ifEmpty([]),//removeWindowDup/${sample}_removeWindowDup.log (#Number of duplicates: nnnn)
        //distribUMIs
        chMqcDistribUMI.collect().ifEmpty([]), //pour config graph
        chPeaksCountsMqc.collect().ifEmpty([]),
        chFripResults.collect().ifEmpty([]),
        chPeaksQCMqc.collect().ifEmpty([]),
        chDeeptoolsProfileMqc.collect().ifEmpty([])
      )
      chMqcReport = multiqc.out.report.toList()
    }

  emit:
  bam = chNoDupBam
  bai = chNoDupBai
  bigwig = chBigWig
  matrixTSS = chTssMatrices
  matrixBin = chBinMatrices
  mqcreport = chMqcReport
}