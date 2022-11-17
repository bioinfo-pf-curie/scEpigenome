/* 
 * Create sparse binned matrix
 */

include { gtfToTSSBed } from '../../local/process/gtfToTSSBed'
include { nbBarcodes } from '../../local/process/nbBarcodes'
include { createTssMatrices } from '../../local/process/createTssMatrices'

workflow countMatricesPerTSS {

  take:
  bam 
  bai
  bcList
  gtf

  main:

  gtfToTSSBed(
    gtf
  )

  nbBarcodes(
    bcList
  )
  chNbBarcodes=nbBarcodes.out.count

  chVersions = Channel.empty()

  createTssMatrices(
    gtfToTSSBed.out.bed,
    bam.join(bai).join(chNbBarcodes)
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  emit:
  matrix = createTssMatrices.out.matrix
  versions = chVersions
}