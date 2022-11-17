/* 
 * Create sparse binned matrix
 */

include { nbBarcodes } from '../../local/process/nbBarcodes'
include { createBinMatrices } from '../../local/process/createBinMatrices'

workflow countMatricesPerBin {

  take:
  bins
  bam 
  bai
  bcList

  main:

  nbBarcodes(
    bcList
  )

  chVersions = Channel.empty()

  createBinMatrices(
    nbBarcodes.out.count,
    bam.join(bai).combine(bins)
  )
  chVersions = chVersions.mix(createBinMatrices.out.versions)

  emit:
  matrix = createBinMatrices.out.matrix
  versions = chVersions
}