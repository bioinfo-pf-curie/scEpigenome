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

  nbBarcodes(
    bcList
  )
  chNbBc=nbBarcodes.out.count

  chNbBc.view()

  createMatrices(
    chNbBc
    bam.join(bai).combine(bins)
  )
  chVersions = chVersions.mix(createMatrices.out.versions)

  emit:
  matrix = createMatrices.out.matrix
  versions = chVersions
}