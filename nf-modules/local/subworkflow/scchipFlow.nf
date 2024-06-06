include { starAlign } from '../../common/process/star/starAlign'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as indexFilter } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { markDuplicates } from '../../common/process/picard/markDuplicates'
include { mergeSamFiles } from '../../common/process/picard/mergeSamFiles'
include { catTxt as mergeBarcodes } from '../../common/process/cat/catTxt'
include { cutadapt } from '../../common/process/cutadapt/cutadapt'
include { barcode2rg } from '../../local/process/barcode2rg'
include { extractBarcodeFlow } from '../../local/subworkflow/extractBarcodeFlow'


// ?? remove duplicates ??
include { removePCRdup } from '../../local/process/removePCRdup' 
include { removeRTdup } from '../../local/process/removeRTdup'
include { removeWindowDup } from '../../local/process/removeWindowDup'

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

workflow scchipFlow {

  take:
  reads
  starIndex

  main:

  chVersions = Channel.empty()

  // Manage multiple fastq files for the same sample
  chReads = reads.groupTuple()
    .flatMap { it -> setMetaChunk(it) }
    .collate(2)

  // The cell barcode is always on R2
  chBarcodeReads = chReads.map{ it -> [it[0], it[1][1]] }

  // Get barcodes information
  if (params.darkCycleDesign){
    chBcInfo = Channel.from(params.barcodesIndrop)
      .flatMap()
      .map { it -> [ [id:it.key], it.value['start_darkcycle'], it.value['size'], file(it.value['bwt2']) ] }
  }else{
    chBcInfo = Channel.from(params.barcodesIndrop)
      .flatMap()
      .map { it -> [ [id:it.key], it.value['start_nodarkcycles'], it.value['size'], file(it.value['bwt2']) ] }
  }

  // Extrat barcode sequence
  extractBarcodeFlow(
    chReads,
    chBarcodeReads,
    chBcInfo
  )
  chVersions = chVersions.mix(extractBarcodeFlow.out.versions)                                                                                                                                              
  // Trim R2 reads
  chR2reads = chReads.map{ meta, reads ->
    newMeta = meta.clone()
    newMeta.singleEnd = true
    [newMeta, reads[1]]
  }

  cutadapt(
    chR2reads
  )

  chTrimReads = cutadapt.out.fastq
    .map{ meta, r2 ->
      newMeta = [id: meta.id, name: meta.name, protocol: meta.protocol, chunk: meta.chunk, part:meta.part]
      [newMeta, r2]
    }.join(chReads)
    .map{meta, r2trim, reads ->
      [meta, [reads[0], r2trim] ]
    }

  // Alignment on reference genome
  starAlign(
    chTrimReads,
    starIndex,
    Channel.value([])
  )
  chVersions = chVersions.mix(starAlign.out.versions)

  samtoolsIndex(
    starAlign.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  // Add barcodes as read groups
  barcode2rg(
    starAlign.out.bam.join(samtoolsIndex.out.bai).join(extractBarcodeFlow.out.barcodes)
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  // Merge multiple BAM files from the same sample
  chBams = barcode2rg.out.bam
    .map{meta, bam, bai ->
       def newMeta = [ id: meta.id, name: meta.name, protocol: meta.protocol, part:meta.part ]
       [ groupKey(newMeta, meta.part), bam, bai ]
     }.groupTuple()
     .branch {
       single: it[0].part <= 1
       multiple: it[0].part > 1
     }

  mergeSamFiles(
    chBams.multiple
  )

  // Merge multiple Barcode files from the same sample
  chAllBarcodes = extractBarcodeFlow.out.barcodes
    .map{meta, bc ->
       def newMeta = [ id: meta.id, name: meta.name, protocol: meta.protocol, part:meta.part ]
       [ groupKey(newMeta, meta.part), bc ]
     }.groupTuple()
     .branch {
       single: it[0].part <= 1
       multiple: it[0].part > 1
     }

  mergeBarcodes(
    chAllBarcodes.multiple,
    Channel.of(true).collect()
  )

  // Mark reads duplicates
  markDuplicates(
    mergeSamFiles.out.bam.mix(chBams.single)
  )
  chVersions = chVersions.mix(markDuplicates.out.versions)

  // Filter out aligned reads
  samtoolsFilter(
    markDuplicates.out.bam
  )
  chVersions = chVersions.mix(samtoolsFilter.out.versions)

  indexFilter(
    samtoolsFilter.out.bam
  )

  samtoolsFlagstat(
    samtoolsFilter.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  //  removePCRdup(
  //    //inputs
  //    chTaggedBam
  //  )
  //  chRemovePCRdupBam = removePCRdup.out.bam
  //  chRemovePCRdupSam = removePCRdup.out.sam
  //  chRemovePCRdupSummary = removePCRdup.out.count
  //  chR1unmappedR2Summary = removePCRdup.out.countR1unmapped

  //  removeRTdup(
  //    chTaggedBam.join(chRemovePCRdupBam).join(chRemovePCRdupSam)
  //  )
  //  chRemovePcrRtBam = removeRTdup.out.bam
  //  chRemoveRtSummary = removeRTdup.out.logs
    
  //  removeWindowDup(
  //    chRemovePcrRtBam
  //  )
  //  chRemoveBlackReg = removeWindowDup.out.bam
  //  removeWindowDup = removeWindowDup.out.logs

  //  countSummary(
  //    chRemovePCRdupSummary.join(chTaggedBam).join(chR1unmappedR2Summary).join(chRemoveRtSummary).join(removeWindowDup), 
  //  )
  //  chDedupCountSummary = countSummary.out.logs

  emit:
  bam = samtoolsFilter.out.bam.join(indexFilter.out.bai)
  barcodes = mergeBarcodes.out.txt.mix(chAllBarcodes.single)
  versions = chVersions
}