include { starAlign } from '../../common/process/star/starAlign'
include { bwaMem2 } from '../../common/process/bwaMem2/bwaMem2'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex as samtoolsIndexMarkdup } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexFilter } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat as filterAlignedStat} from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'
include { samtoolsFlagstat as markdupStat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFixmate } from '../../common/process/samtools/samtoolsFixmate'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsMarkdup } from '../../common/process/samtools/samtoolsMarkdup'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { barcode2tag } from '../../local/process/barcode2tag'
include { removeExtraDup } from '../../local/process/removeExtraDup'
include { pairToBed as rmBlackList } from '../../common/process/bedtools/pairToBed'
include { getTagfragmentCounts } from '../../local/process/getTagfragmentCounts'
include { weightedDistrib } from '../../local/process/weightedDistrib'

include {checkAlignmentPercent} from '../../../lib/functions'

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
    chStar = starAlign.out.bam
    chStarLogs = starAlign.out.logs

  // Filter removes all 'aligned' channels that fail the check
  chStarLogs.join(chStar)
    .filter { meta, logs, bam -> checkAlignmentPercent(meta, logs) }
    .map { meta, logs, bam -> [ meta, bam ] }
    .set { chBams }

  }else if (params.aligner = "bwa-mem2"){
    bwaMem2(
      reads,
      index,
      Channel.value([])
    )
    chVersions = chVersions.mix(bwaMem2.out.versions)
    chBams = bwaMem2.out.bam
  }


  // Add barcodes as read tag
  barcode2tag(
    chBams.map{meta, bam -> [meta, bam, []]}
  )
  chVersions = chVersions.mix(barcode2tag.out.versions)

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

  // merge batches
  samtoolsMerge(
    chTaggedBams.multiple
  )

  samtoolsStats(
    samtoolsMerge.out.bam.mix(chTaggedBams.single),
    Channel.value([])
  )
  chVersions = chVersions.mix(samtoolsStats.out.versions)

  //*******************************************************
  // Remove blacklist using name sorted BAM

  rmBlackList(
    samtoolsMerge.out.bam.mix(chTaggedBams.single),
    blackList
  )
  chVersions = chVersions.mix(rmBlackList.out.versions)
  chNameSortedBam = params.keepBlackList ? samtoolsMerge.out.bam.mix(chTaggedBams.single) : rmBlackList.out.bam
 
  //********************************************************
  // Mark PCR reads duplicates

  samtoolsFixmate(
    chNameSortedBam
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

  // Stats on mapped+unmapped reads, including duplicates
  markdupStat(
    chMdBam
  )
  chVersions = chVersions.mix(markdupStat.out.versions)

  // Index markdup file
  samtoolsIndexMarkdup(
    chMdBam
  )
  chVersions = chVersions.mix(samtoolsIndexMarkdup.out.versions)

  //********************************************************
  // Filter out aligned reads
  
  samtoolsFilter(
    chMdBam
  )
  chVersions = chVersions.mix(samtoolsFilter.out.versions)
                                                                                                                                                                                                       
  samtoolsIndexFilter(
    samtoolsFilter.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndexFilter.out.versions)

  filterAlignedStat(
    samtoolsFilter.out.bam
  )
  chVersions = chVersions.mix(filterAlignedStat.out.versions)

  //*********************************************************
  // Get barcodes information from the final BAM file

  getTagfragmentCounts(
    samtoolsFilter.out.bam
  )
  chVersions = chVersions.mix(getTagfragmentCounts.out.versions)

  weightedDistrib(
    getTagfragmentCounts.out.counts
  )
  chVersions = chVersions.mix(weightedDistrib.out.versions)

  emit:
  bam = samtoolsFilter.out.bam.join(samtoolsIndexFilter.out.bai)
  mdLogs = samtoolsMarkdup.out.logs.mix(removeExtraDup.out.logs)
  stats = filterAlignedStat.out.stats.mix(markdupStat.out.stats).mix(samtoolsStats.out.stats)
  barcodes = getTagfragmentCounts.out.barcodes 
  counts = getTagfragmentCounts.out.counts
  whist = weightedDistrib.out.mqc
  versions = chVersions
}