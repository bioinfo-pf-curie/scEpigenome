/* 
 * Create sparse binned matrix
 */

include { nbBarcodes } from '../../local/process/nbBarcodes'
include { createMatrices } from '../../local/process/createMatrices'

workflow countMatricesPerBin {

  take:
  bins
  bam 
  bai
  bcList

  main:
  
  chVersions = Channel.empty()
  chInputs = bam.join(bai).combine(bins)

  nbBarcodes(
    bcList
  )

  createMatrices(
    nbBarcodes.out.count
    chInputs
  )
  chVersions = chVersions.mix(createMatrices.out.versions)

  emit:
  matrix = createMatrices.out.matrix
  versions = chVersions
}