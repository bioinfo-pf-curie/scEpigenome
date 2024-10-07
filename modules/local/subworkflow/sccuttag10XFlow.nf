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

workflow sccuttag10XFlow {

  take:
  reads

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

  // Barcode information
  chBcInfo = Channel.of([[id:'barcode10X'], params.barcodes10X['start'], 
    params.barcodes10X['len'], file(params.barcodes10X['bwt2'])])
  chBcMapQ = Channel.of(params.mapqBarcode).collect()
  
  // Extrat barcode sequence
  extractBarcodeFlow(
    chDNAReads,
    chBarcodeReads,
    chBcInfo,
    chBcMapQ
  )
  chVersions = chVersions.mix(extractBarcodeFlow.out.versions)

  emit:
  reads = extractBarcodeFlow.out.fastq
  barcodes = extractBarcodeFlow.out.barcodes
  logs = extractBarcodeFlow.out.mappingLogs
  stats = extractBarcodeFlow.out.stats
  versions = chVersions 
}