include { getBarcodes } from '../../local/process/getBarcodes'
include { alignBarcodes } from '../../local/process/alignBarcodes'
include { addBarcodes } from '../../local/process/addBarcodes'
include { joinBcIndexes } from '../../local/process/joinBcIndexes'
include { getUniqueBarcodes } from '../../local/process/getUniqueBarcodes'

/*
 * Extract Barcodes from reads and align them using a barcode reference
 */

workflow extractBarcodeFlow {

  take:
  dnaReads // [meta,reads1, reads2]
  barcodeReads // [meta, reads]
  barcodeInfo // [meta, start, len, index]

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
    chBarcode2align
  )
  chVersions = chVersions.mix(alignBarcodes.out.versions)

  // Remove index ids and join barcodes
  chBarcodes2merge = alignBarcodes.out.readBarcodes
    .map { meta, fasta ->
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
  chFinalBarcodes = chBarcodes2merge.single
    .mix(joinBcIndexes.out.results)

  addBarcodes(
    dnaReads.join(chFinalBarcodes)
  )
  chVersions = chVersions.mix(addBarcodes.out.versions)

  getUniqueBarcodes(
    chFinalBarcodes
  )

  // What is the INPUT ? is it after all filtering ?
  //    distribUMIs(
  //      chfinalBClist
  //    )
  //    chMqcDistribUMI = distribUMIs.out.mqc
  //    chPdfDist = distribUMIs.out.pdf
  //    chVersions = chVersions.mix(removeBlackRegions.out.versions)

  emit:
  versions = chVersions
  barcodes = getUniqueBarcodes.out.barcodes
  fastq = addBarcodes.out.fastq
}