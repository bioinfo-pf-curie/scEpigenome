/* 
 * Create sparse binned matrix
 */

include { extractTSS } from '../../local/process/extractTSS'
include { countMatricesPerFeature } from '../../local/process/countMatricesPerFeature'
include { countMatricesPerBin } from '../../local/process/countMatricesPerBin'
include { weightedDistrib } from '../../local/process/weightedDistrib'

workflow countMatricesFlow {

  take:
  bam //[meta, bam, bai]
  bcList
  bins
  geneBed

  main:
  chVersions = Channel.empty()

  chCountMatricesPerBin = Channel.empty()
  if (params.createBinMatrix){
  // Counts per genomic bins
  countMatricesPerBin(
    bam.join(bcList).combine(bins)
  )
  chCountMatricesPerBin= countMatricesPerBin.out.matrix
  chCountMatricesPerBinLogs = countMatricesPerBin.out.logs
  chVersions = chVersions.mix(countMatricesPerBin.out.versions)
  }

  // Counts on TSS
  extractTSS(
    geneBed
  )

  countMatricesPerFeature(
    bam.join(bcList),
    extractTSS.out.tss,
  )
  chVersions = chVersions.mix(countMatricesPerFeature.out.versions)

  emit:
  bins = chCountMatricesPerBin.ifEmpty([])
  binLogs= chCountMatricesPerBinLogs.ifEmpty([])
  tss = countMatricesPerFeature.out.matrix
  tssLogs= countMatricesPerFeature.out.logs
  versions = chVersions
}