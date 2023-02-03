/* 
 * Create sparse binned matrix
 */

include { gtfToTSSBed } from '../../local/process/gtfToTSSBed'
include { countMatricesPerTSS } from '../../local/process/countMatricesPerTSS'

workflow countMatricesPerTSSFlow {

  take:
  bam 
  bai
  bcList
  gtf

  main:

  gtfToTSSBed(
    gtf
  )

  chVersions = Channel.empty()

  countMatricesPerTSS(
    gtfToTSSBed.out.bed,
    bam.join(bai)
  )
  chVersions = chVersions.mix(countMatricesPerTSS.out.versions)

  emit:
  matrix = countMatricesPerTSSFlow.out.matrix
  versions = chVersions
}