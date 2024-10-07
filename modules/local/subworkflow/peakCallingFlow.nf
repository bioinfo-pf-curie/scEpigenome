/* 
 * Run peak calling analysis using pseudo-bulk profiles
 */

include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { macs2 } from '../../common/process/macs2/macs2'
include { annotatePeaks } from '../../common/process/homer/annotatePeaks'
include { bedtoolsMerge } from '../../common/process/bedtools/bedtoolsMerge'
include { frip} from '../../local/process/frip'
include { peakQC } from '../process/peakQC'

chPeakCountHeader = Channel.fromPath("$projectDir/assets/peak_count_header.txt")
chFripScoreHeader = Channel.fromPath("$projectDir/assets/frip_score_header.txt")
chPeakAnnotationHeader = Channel.fromPath("$projectDir/assets/peak_annotation_header.txt")

// Create special channel to deal with no input cases
chNoInput = Channel.from( [ [], [] ] ).toList()

workflow peakCallingFlow {

  take:
  bam //[meta, bam, bai]
  effgsize
  gtf
  fasta

  main:

  chVersions = Channel.empty()

  // Macs2 peak calling
  macs2(
    bam.combine(chNoInput),
    effgsize.collect(),
    chPeakCountHeader.collect()
  )
  chVersions = chVersions.mix(macs2.out.versions)

  bedtoolsMerge(
    macs2.out.peaks
  )
  chMergePeaksBed = bedtoolsMerge.out.bed  
  chVersions = chVersions.mix(bedtoolsMerge.out.versions)

  // FRIP
  samtoolsFlagstat(
    bam.map{it -> [it[0], it[1]]}
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  chFripInput = bam.map{it -> [it[0], it[1]]}
    .join(samtoolsFlagstat.out.stats).join(chMergePeaksBed)

  frip(
    chFripInput,
    chFripScoreHeader.collect()
  )
  chVersions = chVersions.mix(frip.out.versions)

  // Annotate peaks
  annotatePeaks(
    chMergePeaksBed,
    gtf.collect(),
    fasta.collect()
  )

  // Peaks quality controls
  peakQC(
    chMergePeaksBed.map{it->it[1]}.collect(),
    annotatePeaks.out.output.map{it->it[1]}.collect(),
    chPeakAnnotationHeader
  )
  chVersions = chVersions.mix(peakQC.out.versions)

  emit:
  peaksOutput = macs2.out.outputXls 
  peaksCountsMqc = macs2.out.mqc
  peaksMerged = chMergePeaksBed
  frip = frip.out.frip
  qc = peakQC.out.mqc
  versions = chVersions
}
