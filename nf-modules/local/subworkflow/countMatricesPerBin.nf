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

  bam.join(bai).combine(bins).set{chBC}

  chBC.view()

  createMatrices(
    nbBarcodes.out.count
    //chBC
  )
  chVersions = chVersions.mix(createMatrices.out.versions)

  emit:
  matrix = createMatrices.out.matrix
  versions = chVersions
}