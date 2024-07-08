include { starAlign } from '../../common/process/star/starAlign'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFlagstat as mappingStat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsMarkdup } from '../../common/process/samtools/samtoolsMarkdup'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { samtoolsFixmate } from '../../common/process/samtools/samtoolsFixmate'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { mergeBarcodes } from '../../local/process/mergeBarcodes'
include { barcode2tag } from '../../local/process/barcode2tag'
include { reverseComplement } from '../../local/process/reverseComplement'
include { extractBarcodeFlow } from '../../local/subworkflow/extractBarcodeFlow'
include { intersectBed as rmBlackList } from './nf-modules/common/process/bedtools/intersectBed'

// Set the meta.chunk value in case of multiple sequencing lanes
def setMetaChunk(row){
  def map = []
  row[1].eachWithIndex() { file, i ->
    meta = row[0].clone()
    meta.chunk = i+1
    meta.part = row[1].size()
    map += [meta, file]
  }
  return map
}

workflow sccuttagIndropFlow {

  take:
  reads
  starIndex

  main:
  chVersions = Channel.empty()
   
  // Manage multiple fastq files for the same sample
  chReads = reads.groupTuple()
    .flatMap { it -> setMetaChunk(it) }
    .collate(2)

  // DNA reads
  chDNAReads = chReads.map{it -> [it[0], [it[1][0],it[1][2]]]}

  // The cell barcode is always on R2
  chBarcodeReads = chReads.map{ it -> [it[0], it[1][1]] }

  // Get barcodes information
  // /!\ +1 for the cut&tag indrop protocol
  if (params.darkCycleDesign){
    chBcInfo = Channel.from(params.barcodesIndrop)
      .flatMap()
      .map { it -> [ [id:it.key], it.value['start_darkcycle'] +1, it.value['size'], file(it.value['bwt2']) ] }
  }else{
    chBcInfo = Channel.from(params.barcodesIndrop)
      .flatMap()
      .map { it -> [ [id:it.key], it.value['start_nodarkcycles'] +1, it.value['size'], file(it.value['bwt2']) ] }
  }
  chBcMapQ = Channel.of(params.mapqBarcode).collect()

  // Reverse Complement the barcode reads
  reverseComplement(
    chBarcodeReads
  )
  chVersions = chVersions.mix(reverseComplement.out.versions)

  // Extrat barcode sequence
  extractBarcodeFlow(
    chDNAReads,
    reverseComplement.out.reads,
    chBcInfo,
    chBcMapQ
  )
  chVersions = chVersions.mix(extractBarcodeFlow.out.versions)

  // Alignment on reference genome
  starAlign(
    extractBarcodeFlow.out.fastq,
    starIndex,
    Channel.value([])
  )
  chVersions = chVersions.mix(starAlign.out.versions)

  // Add barcodes as read groups
  barcode2tag(
    starAlign.out.bam.join(extractBarcodeFlow.out.barcodes)
  )
  chVersions = chVersions.mix(barcode2tag.out.versions)                                                                                                                                                 

  // Merge multiple BAM files from the same sample
  chBams = barcode2tag.out.bam
    .map{meta, bam ->
      def newMeta = [ id:meta.id, name:meta.name, protocol:meta.protocol, part:meta.part ]
      [ groupKey(newMeta, meta.part), bam ]
    }.groupTuple()
     .branch {
       single: it[0].part <= 1
       multiple: it[0].part > 1
     }

  samtoolsMerge(
    chBams.multiple
  )

  mappingStat(
    samtoolsMerge.out.bam.mix(chBams.single)
  )
  chVersions = chVersions.mix(mappingStat.out.versions)

  //********************************************************
  // Mark reads duplicates

  samtoolsFixmate(
    samtoolsMerge.out.bam.mix(chBams.single)
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

  //********************************************************
  // Filter out aligned reads

  samtoolsFilter(
    samtoolsMarkdup.out.bam
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

  emit:
  bam = samtoolsFilter.out.bam.join(samtoolsIndex.out.bai)
  bcLogs = extractBarcodeFlow.out.mappingLogs.mix(extractBarcodeFlow.out.stats)
  starLogs = starAlign.out.logs
  mdLogs = samtoolsMarkdup.out.logs
  stats = mappingStat.out.stats.mix(samtoolsFlagstat.out.stats)
  barcodes = getTagValues.out.barcodes
  counts = getTagValues.out.counts
  versions = chVersions
}