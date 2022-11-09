/*
 * Bam to BigWig files
 */

process bigwig {
  tag "${meta.id}"
  label 'deeptools'
  label 'extraCpu'
  label 'extraMem'

  when:
  !params.skipBigWig

  input:
  tuple val(meta), path(bam), path(bai)
  path (blackListBed)

  output:
  tuple val(meta), path("*.bw"), emit: bigwig
  tuple val(meta), path("*_bamToBigWig.log"), emit: logs
  path ("versions.txt"), emit: versions
  
  script:
  def prefix = task.ext.prefix ?: "${bam.baseName}"
  """
  if [[ "${params.removeBlackRegions}" == "true" ]]
  then
      bamCoverage --bam ${bam} --outFileName ${prefix}.bw --numberOfProcessors ${task.cpus} --normalizeUsing CPM --ignoreForNormalization chrX --binSize 50 --smoothLength 500 --extendReads 150 --blackListFileName ${blackListBed} &> ${prefix}_bamToBigWig.log
  else
      bamCoverage --bam ${bam} --outFileName ${prefix}.bw --numberOfProcessors ${task.cpus} --normalizeUsing CPM --ignoreForNormalization chrX --binSize 50 --smoothLength 500 --extendReads 150 &> ${prefix}_bamToBigWig.log
  fi

  deeptools --version &> versions.txt
  """
}
