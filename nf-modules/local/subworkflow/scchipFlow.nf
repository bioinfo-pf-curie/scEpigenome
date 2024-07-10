include { cutadapt } from '../../common/process/cutadapt/cutadapt'
include { barcode2tag } from '../../local/process/barcode2tag'
include { extractBarcodeFlow } from '../../local/subworkflow/extractBarcodeFlow'

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
  chBcMapQ = Channel.of(params.mapqBarcode).collect()

  // Extrat barcode sequence
  extractBarcodeFlow(
    chReads,
    chBarcodeReads,
    chBcInfo,
    chBcMapQ
  )
  chVersions = chVersions.mix(extractBarcodeFlow.out.versions)

  // Trim R2 reads
  chDNAreads = extractBarcodeFlow.out.fastq
  chR2reads = chDNAreads.map{ meta, reads ->
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
    }.join(chDNAreads)
    .map{meta, r2trim, reads ->
      [meta, [reads[0], r2trim] ]
    }

  emit:
  reads = chTrimReads
  barcodes = extractBarcodeFlow.out.barcodes
  logs = extractBarcodeFlow.out.mappingLogs
  versions = chVersions
}