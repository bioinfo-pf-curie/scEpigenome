/* 
 * Create sparse binned matrix
 */

include { nbBarcodes } from '../../local/process/nbBarcodes'
include { createMatrices } from '../../local/process/createMatrices'

workflow countMatricesPerBin {

  take:
  bam.join(bai).combine(bins)
  bcList

  main:
  
  chVersions = Channel.empty()

  nbBarcodes(
    bcList
  )

  createMatrices(
    nbBarcodes.out.count
    bam.join(bai).combine(bins)
  )
  chVersions = chVersions.mix(createMatrices.out.versions)

  emit:
  matrix = createMatrices.out.matrix
  versions = chVersions
}