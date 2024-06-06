/* 
 * Create sparse binned matrix
 */

include { gtfToTSSBed } from '../../local/process/gtfToTSSBed'
include { countMatricesPerTSS } from '../../local/process/countMatricesPerTSS'
include { countMatricesPerBin } from '../../local/process/countMatricesPerBin'

workflow countMatricesFlow {

  take:
  bam //[meta, bam, bai]
  bcList
  bins
  gtf

  main:
  chVersions = Channel.empty()

  // Counts per genomic bins
  countMatricesPerBin(
    bam.join(bcList).combine(bins)
  )
  chVersions = chVersions.mix(countMatricesPerBin.out.versions)

  // Counts on TSS
  gtfToTSSBed(
    gtf
  )
  //chVersions = chVersions.mix(gtfToTSSBed.out.versions)

  countMatricesPerTSS(
    gtfToTSSBed.out.bed,
    bam,
    bcList
  )
  chVersions = chVersions.mix(countMatricesPerTSS.out.versions)

  emit:
  bins = countMatricesPerBin.out.matrix
  tss = countMatricesPerTSS.out.matrix
  versions = chVersions
}