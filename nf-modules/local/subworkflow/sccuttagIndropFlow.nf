include { reverseComplement } from '../../local/process/reverseComplement'
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

workflow sccuttagIndropFlow {

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

  emit:
  reads = extractBarcodeFlow.out.reads
  barcodes = extractBarcodeFlow.out.barcodes
  versions = chVersions
}