/* 
 * Manage bigwig files
 */

include { deeptoolsBamCoverage } from '../../common/process/deeptools/deeptoolsBamCoverage'
include { deeptoolsComputeMatrix } from '../../common/process/deeptools/deeptoolsComputeMatrix'

workflow bigwigFlow {

  take:
  bam // [meta, bam, bai]
  effGenomeSize
  geneBed

  main:
  chVersions = Channel.empty()
  bam.view()
  deeptoolsBamCoverage(
    bam.map(it -> [it[0], it[1], it[2], []]),
    Channel.value([]),
    effGenomeSize
  )
  chVersions = chVersions.mix(deeptoolsBamCoverage.out.versions)

  deeptoolsComputeMatrix(
    deeptoolsBamCoverage.out.bigwig,
    geneBed.collect()
  )
  chVersions = chVersions.mix(deeptoolsComputeMatrix.out.versions)

  emit:
  bigwig = deeptoolsBamCoverage.out.bigwig
  matrix = deeptoolsComputeMatrix.out.mqc
  versions = chVersions
}