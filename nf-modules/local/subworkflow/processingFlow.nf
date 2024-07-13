include { starAlign } from '../../common/process/star/starAlign'
include { bwaMem2 } from '../../common/process/bwaMem2/bwaMem2'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'
include { samtoolsFlagstat as markdupStat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFixmate } from '../../common/process/samtools/samtoolsFixmate'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsMarkdup } from '../../common/process/samtools/samtoolsMarkdup'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { barcode2tag } from '../../local/process/barcode2tag'
include { removeExtraDup } from '../../local/process/removeExtraDup'
include { intersectBed as rmBlackList } from '../../common/process/bedtools/intersectBed'
include { getTagValues } from '../../local/process/getTagValues'
include { weightedDistrib } from '../../local/process/weightedDistrib'

workflow processingFlow {

  take:
  reads
  index
  blackList

  main:

  chVersions = Channel.empty()

  // Alignment on reference genome
  if (params.aligner == "star"){
    starAlign(
      reads,
      index,
      Channel.value([])
    )
    chVersions = chVersions.mix(starAlign.out.versions)
    chBams = starAlign.out.bam
  }else if (params.aligner = "bwa-mem2"){
    bwaMem2(
      chReads,
      index,
      Channel.of([])
    )
    chVersions = chVersions.mix(bwaMem2.out.versions)
    chBams = bwaMem2.out.bam
  }

  chBams.view()

  // Add barcodes as read tag
  barcode2tag(
    chBams.map{meta, bam -> [meta, bam, []]}
  )
  chVersions = chVersions.mix(barcode2tag.out.versions)

  barcode2tag.out.bam.view()

  // Merge multiple BAM files from the same sample
  chTaggedBams = barcode2tag.out.bam
    .map{meta, bam ->
       def newMeta = [ id: meta.id, name: meta.name, protocol: meta.protocol, part:meta.part ]
       [ groupKey(newMeta, meta.part), bam ]
     }.groupTuple()
     .branch {
       single: it[0].part <= 1
       multiple: it[0].part > 1
     }

  chTaggedBams.multiple.view()

  samtoolsMerge(
    chTaggedBams.multiple
  )

  samtoolsStats(
    samtoolsMerge.out.bam.mix(chTaggedBams.single),
    Channel.value([])
  )
  chVersions = chVersions.mix(samtoolsStats.out.versions)

  //********************************************************
  // Mark PCR reads duplicates

  samtoolsFixmate(
    samtoolsMerge.out.bam.mix(chTaggedBams.single)
  )
  chVersions = chVersions.mix(samtoolsFixmate.out.versions)

  samtoolsSort(
    samtoolsFixmate.out.bam
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsMarkdup(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsMarkdup.out.versions)

  // Addition duplicates (RT / window)
  removeExtraDup(
    samtoolsMarkdup.out.bam
  )
  chMdBam = params.extraDup ? removeExtraDup.out.bam : samtoolsMarkdup.out.bam

  // Stats on mapped reads including duplicates
  markdupStat(
    chMdBam
  )
  chVersions = chVersions.mix(markdupStat.out.versions)

  //********************************************************
  // Filter out aligned reads
  
  samtoolsFilter(
    chMdBam
  )
  chVersions = chVersions.mix(samtoolsFilter.out.versions)

  rmBlackList(
    samtoolsFilter.out.bam,
    blackList
  )
  chVersions = chVersions.mix(rmBlackList.out.versions)
  chFinalBam = params.keepBlackList ? samtoolsFilter.out.bam : rmBlackList.out.bam
                                                                                                                                                                                                           
  samtoolsIndex(
    chFinalBam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    chFinalBam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  //*********************************************************
  // Get barcodes information from the final BAM file

  getTagValues(
    chFinalBam.join(samtoolsIndex.out.bai)
  )
  chVersions = chVersions.mix(getTagValues.out.versions)

  weightedDistrib(
    getTagValues.out.counts
  )
  chVersions = chVersions.mix(weightedDistrib.out.versions)

  emit:
  bam = chFinalBam.join(samtoolsIndex.out.bai)
  mdLogs = samtoolsMarkdup.out.logs.mix(removeExtraDup.out.logs)
  stats = samtoolsFlagstat.out.stats.mix(markdupStat.out.stats).mix(samtoolsStats.out.stats)
  barcodes = getTagValues.out.barcodes 
  counts = getTagValues.out.counts
  whist = weightedDistrib.out.mqc
  versions = chVersions
}