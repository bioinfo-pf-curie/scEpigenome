/* 
 * Create sparse binned matrix
 */

include { gtfToTSSBed } from '../../local/process/gtfToTSSBed'
include { countMatricesPerTSS } from '../../local/process/countMatricesPerTSS'

workflow countMatricesPerTSSFlow {

  take:
  bamBai
  bcList
  gtf

  main:

  gtfToTSSBed(
    gtf
  )

  chVersions = Channel.empty()

  countMatricesPerTSS(
    gtfToTSSBed.out.bed,
    bamBai
  )
  chVersions = chVersions.mix(countMatricesPerTSS.out.versions)

  emit:
  matrix = countMatricesPerTSS.out.matrix
  versions = chVersions
}