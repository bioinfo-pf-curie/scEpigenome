/* 
 * Create sparse binned matrix
 */

include { gtfToTSSBed } from '../../local/process/gtfToTSSBed'
include { countMatricesPerFeature } from '../../local/process/countMatricesPerFeature'
include { countMatricesPerBin } from '../../local/process/countMatricesPerBin'

workflow countMatricesFlow {

  take:
  bam //[meta, bam, bai]
  //bcList
  bins
  gtf

  main:
  chVersions = Channel.empty()

  // Counts per genomic bins
  countMatricesPerBin(
    bam.map{meta, bam, bai -> [ meta, bam, bai, [] ]}.combine(bins)
  )
  chVersions = chVersions.mix(countMatricesPerBin.out.versions)

  // Counts on TSS
  gtfToTSSBed(
    gtf
  )
  //chVersions = chVersions.mix(gtfToTSSBed.out.versions)

  countMatricesPerFeature(
    bam.map{meta, bam, bai -> [ meta, bam, bai, [] ]},
    gtfToTSSBed.out.bed,
  )
  chVersions = chVersions.mix(countMatricesPerFeature.out.versions)

  emit:
  bins = countMatricesPerBin.out.matrix
  tss = countMatricesPerFeature.out.matrix
  versions = chVersions
}