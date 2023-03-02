//common
include { getSoftwareVersions } from '../../common/process/utils/getSoftwareVersions'
include { starAlign } from '../../common/process/star/starAlign'
include { deeptoolsBamCoverage } from '../../common/process/deeptools/deeptoolsBamCoverage'
// add preseq
//local
include { concatenate_fastqs_from_10X } from '../../local/process/concatenate_fastqs_from_10X'
include { multiqc } from '../../local/process/multiqc'
include { bcAlign10X } from '../../local/process/bcAlign10X'
include { addFlags } from '../../local/process/addFlags'
  // remove duplicates
include { removePCRdup } from '../../local/process/removePCRdup' // je les passe dans common ?? Non
  // blackRegions
include { removeBlackRegions } from '../../local/process/removeBlackRegions'
  //--------
include { countSummary } from '../../local/process/countSummary' // empty channels pour é ter bug car pas de RT ni Window?
include { distribUMIs } from '../../local/process/distribUMIs'
include { bamToFrag } from '../../local/process/bamToFrag'
include { reverseComplement } from '../../local/process/reverseComplement'
include { countMatricesPerBin } from '../../local/process/countMatricesPerBin'
//subworkflow
include { countMatricesPerTSSFlow } from '../../local/subworkflow/countMatricesPerTSSFlow' 
include { peaksPseudoBulk } from '../../local/subworkflow/peaksPseudoBulk' 
include { deeptoolsComputeMatrix } from '../../common/process/deeptools/deeptoolsComputeMatrix'


workflow sccuttag_10X {

  take:
  reads
  workflowSummaryCh
  multiqcConfigCh
  metadataCh
  sPlanCh
  customRunName
  starIndex
  blackList
  gtf
  fasta
  binsize
  effGenomeSize 
  geneBed 

  main:
    // Init Channels
    // channels never filled
    chStarGtf  = Channel.value([])
    chEffGenomeSize = Channel.value([])
    // channels filled
    // filled avec ifEmpty([]) in mqc
    joinBcIndexesLogs = Channel.empty()
    chRemoveDupLog = Channel.empty()
    // filled by process
    chBigWig= Channel.empty()
    chAlignedLogs = Channel.empty()
    warnCh = Channel.empty()
    chVersions = Channel.empty()
    chNoDupBam  = Channel.empty()
    chNoDupBai  = Channel.empty()
    chTssMatrices  = Channel.empty()
    chBinMatrices  = Channel.empty()
    chMqcReport  = Channel.empty()

    // if BigWig
    chDeeptoolsProfileMqc = Channel.empty()

    /*reads
      .groupTuple()
      .flatten()
      .toList()
      .set{allSamples}*/ //ne marche pas 
  
    reads
      .groupTuple()
      .map{ it -> [it[0], [it[1][0][0], it[1][0][1], it[1][0][2], it[1][1][0], it[1][1][1], it[1][1][2], it[1][2][0], it[1][2][1], it[1][2][2], it[1][3][0], it[1][3][1], it[1][3][2]]]}
      .set{allSamples}

    concatenate_fastqs_from_10X(
      allSamples
    )
    barcodeRead=concatenate_fastqs_from_10X.out.barcodeRead
    dnaRead=concatenate_fastqs_from_10X.out.dnaRead

    // 1) Barcode alignement and extrcation part
    bcAlign10X(
      barcodeRead
    )
    chReadsMatchingIndex = bcAlign10X.out.results
    chReadBcNames = bcAlign10X.out.bcNames
    chIndexBowtie2Logs = bcAlign10X.out.logs
    joinBcIndexesLogs = bcAlign10X.out.counts
    //chRemoveDupLog = bcAlign10X.out.mqc
    chVersions = chVersions.mix(bcAlign10X.out.versions)

    starAlign(
      dnaRead,
      starIndex,
      chStarGtf
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

    chRemovePCRdupSummary
    .map { meta, val -> [ meta, []] }
    .set { chRemoveRtSummary }

    countSummary(
      //inputs
      chRemovePCRdupSummary.join(chRemovePcrBamSummary).join(chR1unmappedR2Summary),
      chRemoveRtSummary
    )
    chDedupCountSummary = countSummary.out.logs

    removeBlackRegions(
      //inputs
      chRemovePCRdupBam,
      blackList.collect()
    )
    chVersions = chVersions.mix(removeBlackRegions.out.versions)
    chNoDupBam = removeBlackRegions.out.bam
    chNoDupBai = removeBlackRegions.out.bai
    chfinalBClist = removeBlackRegions.out.list

    peaksPseudoBulk( 
      chNoDupBam,
      chNoDupBai,
      effGenomeSize,
      gtf,
      fasta
    )
    peaksPseudoBulkBed = peaksPseudoBulk.out.mergedPeaks
    chPeaksCountsMqc = peaksPseudoBulk.out.peaksCountsMqc
    chPeaksSizesMqc = peaksPseudoBulk.out.peaksSizesMqc
    chFripResults = peaksPseudoBulk.out.fripResults
    chPeaksQCMqc = peaksPseudoBulk.out.peaksQCMqc
    chVersions = chVersions.mix(peaksPseudoBulk.out.versions)

    // Subworkflow
    countMatricesPerBin(
      chNoDupBam.join(chNoDupBai).combine(binsize),
      chfinalBClist
    )
    chBinMatrices=countMatricesPerBin.out.matrix
    chVersions = chVersions.mix(countMatricesPerBin.out.versions)

    // Subworkflow
    countMatricesPerTSSFlow(
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
        chNoDupBam.join(chNoDupBai).map{it->[it[0],it[1],it[2],[]]},
        blackList.collect(),
        chEffGenomeSize
      )
      //outputs
      chBigWig = deeptoolsBamCoverage.out.bigwig
      chVersions = chVersions.mix(deeptoolsBamCoverage.out.versions)

      deeptoolsComputeMatrix( 
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

    // delete $meta to input to mqc a path
    chfinalBClist
      .map{it -> it[1]}
      .set{chfinalBClistCollected}
    joinBcIndexesLogs
      .map{it -> it[1]}
      .set{joinBcIndexesLogsCollected}

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
        chAlignedLogs.collect().ifEmpty([]), //star = *Log.final.out
        // bcAlign:
        chIndexBowtie2Logs.collect().ifEmpty([]),//index/${sample}_Bowtie2.log
        // bcSubset:
        joinBcIndexesLogsCollected.collect().ifEmpty([]),//bowtie2/${sample}_bowtie2.log
        // countSummary:
        chDedupCountSummary.collect().ifEmpty([]),//allDup/${sample}_allDup.log
        // countSummary:
        chfinalBClistCollected.collect().ifEmpty([]),//cellThresholds/${sample}_rmDup.txt
        // removeWindowDup:
        chRemoveDupLog.collect().ifEmpty([]),//removeWindowDup/${sample}_removeWindowDup.log (#Number of duplicates: nnnn)
        //distribUMIs
        chMqcDistribUMI.collect().ifEmpty([]), //pour config graph
        // macs2
        chPeaksCountsMqc.collect().ifEmpty([]),
        //peakQC
        chFripResults.collect().ifEmpty([]),
        chPeaksQCMqc.collect().ifEmpty([]),
        // from deeptools matrix
        chDeeptoolsProfileMqc.collect().ifEmpty([]),
        // macs2
        chPeaksSizesMqc.collect().ifEmpty([])
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