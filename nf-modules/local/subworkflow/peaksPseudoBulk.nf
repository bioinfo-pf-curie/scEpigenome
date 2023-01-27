/* 
 * merge all peak beds into one == pseudo bulk
 */

include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

include { macs2 } from '../../common/process/macs2/macs2'
include { bedtoolsMergePeaks } from '../../local/process/bedtoolsMergePeaks'

include { frip} from '../../local/process/frip'
include { annotatePeaks } from '../../common/process/homer/annotatePeaks'


Channel
  .fromPath("$projectDir/assets/peak_count_header.txt")
  .set { chPeakCountHeader }

Channel
  .fromPath("$projectDir/assets/frip_score_header.txt")
  .set { chFripScoreHeader }

Channel
  .fromPath("$projectDir/assets/peak_annotation_header.txt")
  .set{ chPeakAnnotationHeader }


workflow peaksPseudoBulk {

  take:
  bam
  bai
  effgsize
  gtf
  fasta

  main:

  chVersions = Channel.empty()

  /*******************************
   * Macs2  - sharp mode
   */

  macs2(
    bam.join(bai).join(Channel.value([])).join(Channel.value([])),
    effgsize.first(),
    chPeakCountHeader.collect()
  )
  chXls = macs2.out.outputXls
  chPeaks = macs2.out.peaks
  chPeaksMqc = macs2.out.mqc
  chVersions = chVersions.mix(macs2.out.versions)

  bedtoolsMergePeaks(
    chPeaks
  )
  chMergePeaksBed = bedtoolsMergePeaks.out.bed  
  chMergePeaksLogs = bedtoolsMergePeaks.out.logs
  chVersions = chVersions.mix(bedtoolsMergePeaks.out.versions)

  /********************************
   * FRIP
   */

  samtoolsFlagstat(
    bam.map{it -> [it[0], it[1]]}
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  frip(
    bam.join(samtoolsFlagstat.out.stats).join(chMergePeaksBed),
    chFripScoreHeader.collect()
  )
  chFripTsv = frip.out.fripTsv
  chVersions = chVersions.mix(frip.out.versions)

  /********************************
   * Peaks Annotation
   */       
  
  annotatePeaks(
    chMergePeaksBed,
    gtf.collect(),
    fasta.collect()
  )

  peakQC(
    chMergePeaksBed.map{it->it[1]}.collect(),
    annotatePeaks.out.output.map{it->it[1]}.collect(),
    chPeakAnnotationHeader
  )
  chVersions = chVersions.mix(peakQC.out.versions)
  chPeakQC = peakQC.out.mqc

  emit:
  //macs2
  peaksOutput = chXls 
  peaksCountsMqc = chPeaksMqc
  // merged peaks
  mergedPeaks = chMergePeaksBed
  // QC
  peaksQCMqc = chPeakQC
  fripResults = chFripTsv
  versions = chVersions
}