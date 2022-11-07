/* 
 * samtools PCR deduplaication Workflow
 */

include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsSort2 } from '../../common/process/samtools/samtoolsSort2'
include { samtoolsFixmate } from '../../common/process/samtools/samtoolsFixmate'
include { samtoolsMarkdup } from '../../common/process/samtools/samtoolsMarkdup'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'

workflow removePCRdup {

  take:
  bam 

  main:
  
  chVersions = Channel.empty()

  samtoolsSort(
    bam
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsFixmate(
    samtoolsSort.out.bam
  )
  
  samtoolsSort2(
    samtoolsFixmate.out.bam
  )

  samtoolsMarkdup(
    samtoolsSort2.out.bam
  )

  samtoolsIndex(
    samtoolsMarkdup.out.bam
  )

  emit:
  bam = samtoolsMarkdup.out.bam
  logs = samtoolsMarkdup.out.logs
  bai = samtoolsIndex.out.bai
  versions = chVersions
}