include { getBarcodes } from '../../local/process/getBarcodes'
include { alignBarcodes } from '../../local/process/alignBarcodes'
include { addBarcodes } from '../../local/process/addBarcodes'
include { joinBcIndexes } from '../../local/process/joinBcIndexes'
include { getBarcodesCounts } from '../../local/process/getBarcodesCounts'
//include { weightedDistrib } from '../../local/process/weightedDistrib'

/*
 * Extract Barcodes from reads and align them using a barcode reference
 */

workflow extractBarcodeFlow {

  take:
  dnaReads // [meta,reads1, reads2]
  barcodeReads // [meta, reads]
  barcodeInfo // [meta, start, len, index]
  barcodeMapQ // val

  main:
  chVersions = Channel.empty()

  // Extract barcode sequences
  chBarcodeReads = barcodeReads
    .combine(barcodeInfo)
    .map{ meta, reads, meta_bc, start, end, index -> 
      def newMeta = meta.clone()
      newMeta.index_id = meta_bc.id
      [ newMeta, reads, start, end ]
    }

  getBarcodes(
    chBarcodeReads
  )
  chVersions = chVersions.mix(getBarcodes.out.versions)

  // Align them to the bowtie2 indexes
  chBarcode2align = getBarcodes.out.fasta
    .combine(barcodeInfo)
    .filter { it[0].index_id == it[2].id }
    .map{ meta, reads, meta_bc, start, end, index ->
      [ meta, reads, meta_bc, index]
    }

  alignBarcodes(
    chBarcode2align,
    barcodeMapQ
  )
  chVersions = chVersions.mix(alignBarcodes.out.versions)

  // Remove index ids and join barcodes
  chBarcodes2merge = alignBarcodes.out.readBarcodes
    .map{ meta, fasta ->
      def newMeta = [id: meta.id, name: meta.name, protocol: meta.protocol, chunk: meta.chunk, part:meta.part]
      [newMeta, fasta]
    }.groupTuple()
    .branch {
       single: it[0].protocol == 'sccuttag_10X'
       multiple: it[0].protocol != 'sccuttag_10X'
    }

  joinBcIndexes(
    chBarcodes2merge.multiple
  )

  // add barcode info in reads' name
  chFinalBarcodes = chBarcodes2merge.single.mix(joinBcIndexes.out.results)

  dnaReadsWithName=dnaReads.map{ meta, reads ->
      def newMeta = [id: meta.id, name: meta.name, protocol: meta.protocol, chunk: meta.chunk, part:meta.part]
      [newMeta, reads]
    }

  addBarcodes(
    dnaReadsWithName.join(chFinalBarcodes)
  )
  chVersions = chVersions.mix(addBarcodes.out.versions)

  getBarcodesCounts(
    chFinalBarcodes
  )

  emit:
  versions = chVersions
  barcodes = getBarcodesCounts.out.barcodes
  mappingLogs = alignBarcodes.out.logs
  stats = addBarcodes.out.log
  counts = getBarcodesCounts.out.counts
  fastq = addBarcodes.out.fastq
}