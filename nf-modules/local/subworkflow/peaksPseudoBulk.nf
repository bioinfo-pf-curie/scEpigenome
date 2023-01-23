/* 
 * merge all peak beds into one == pseudo bulk
 */

include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

include { macs2 as macs2Sharp} from '../../common/process/macs2/macs2'
include { bedtoolsMergePeaks as mergePeaksSharp} from '../../local/process/bedtoolsMergePeaks'

include { macs2 as macs2Broad} from '../../common/process/macs2/macs2'
include { bedtoolsMergePeaks as mergePeaksBroad} from '../../local/process/bedtoolsMergePeaks'

Channel
  .fromPath("$projectDir/assets/peak_count_header.txt")
  .set { chPeakCountHeader }

workflow peaksPseudoBulk {

  take:
  bam
  bai
  effgsize

  main:

  chVersions = Channel.empty()

  /*******************************
   * Macs2  - Sharp mode
   */

  macs2Sharp(
    bam.join(bai),
    effgsize.first(),
    chPeakCountHeader.collect()
  )
  chXls = macs2Sharp.out.outputXls
  chSharpPeaks = macs2Sharp.out.peaks
  chPeaksMqc = macs2Sharp.out.mqc
  chVersions = chVersions.mix(macs2Sharp.out.versions)

  mergePeaksSharp(
    chSharpPeaks
  )
  chVersions = chVersions.mix(mergePeaksSharp.out.versions)

  /*********************************
   * Macs2 - Broad mode
   */ 

  macs2Broad(
    bam.join(bai),
    effgsize.first(),
    chPeakCountHeader.collect()
  )
  chVersions = chVersions.mix(macs2Broad.out.versions)

  mergePeaksBroad(
    macs2Broad.out.peaks
  )
  chVersions = chVersions.mix(mergePeaksBroad.out.versions)

  /********************************
   * FRIP
   */

  samtoolsFlagstat(
    bam.map{it -> [it[0], it[1]]}
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  frip(
    bam.join(samtoolsFlagstat.out.stats).join(),
    chFripScoreHeader.collect()
  )
  chVersions = chVersions.mix(frip.out.versions)


  emit:
  peaksPseudoBulkBed = mergePeaksSharp.out.bed
  peaksPseudoBulkBed = mergePeaksSharp.out.bed
  versions = chVersions
}