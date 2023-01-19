/* 
 * merge all peak beds into one == pseudo bulk
 */

include { bedtoolsMergePeaks } from '../../local/process/bedtoolsMergePeaks'
include { bedtoolsSort } from '../../common/process/bedtools/bedtoolsSort'

workflow peaksPseudoBulk {

  take:
  peaks

  main:

  chVersions = Channel.empty()

  bedtoolsMergePeaks(
    peaks
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  bedtoolsSort(
    bedtoolsMerge.out.bed,
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  emit:
  peaksBulkBed = bedtoolsSort.out.bed
  versions = chVersions
}