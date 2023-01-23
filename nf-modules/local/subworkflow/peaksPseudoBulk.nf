/* 
 * merge all peak beds into one == pseudo bulk
 */

include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

include { macs2 as macs2Sharp} from '../../common/process/macs2/macs2'
include { bedtoolsMergePeaks as mergePeaksSharp} from '../../local/process/bedtoolsMergePeaks'

include { macs2 as macs2Broad} from '../../local/process/macs2'
include { bedtoolsMergePeaks as mergePeaksBroad} from '../../local/process/bedtoolsMergePeaks'


Channel
  .fromPath("$projectDir/assets/peak_count_header.txt")
  .set { chPeakCountHeader }
Channel
  .fromPath("$projectDir/assets/frip_score_header.txt")
  .set { chFripScoreHeader }

workflow peaksPseudoBulk {

  take:
  bam
  effgsize

  main:

  chVersions = Channel.empty()

  samtoolsFlagstat(
    bam.map{it -> [it[0], it[1]]}
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)


  /*******************************
   * Macs2  - Sharp mode
   */

  macs2Sharp(
    bam
    effgsize.first(),
    chPeakCountHeader.collect()
    //ext.args = "--nomodel --extsize 200 --keep-dup all -f BAM ${params.macs_sharp}"
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  mergePeaksSharp(
    macs2Sharp.out.peaks
    //ext.args = " ${params.max_feature_dist_sharp} "
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  /*********************************
   * Macs2 - Broad mode
   */ 

  macs2Broad(
    bam
    effgsize.first(),
    chPeakCountHeader.collect()
    //ext.args = "--nomodel --extsize 200 --keep-dup all -f BAM ${params.macs_broad}
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  mergePeaksBroad(
    macs2Broad.out.peaks
    //ext.args = " ${params.max_feature_dist_board} "
  )
  chVersions = chVersions.mix(createTssMatrices.out.versions)

  emit:
  peaksPseudoBulkBed = bedtoolsMergePeaks.out.bed
  versions = chVersions
}