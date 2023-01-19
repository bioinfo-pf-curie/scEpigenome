/*
 * Bowtie2 index barcode alignment
 */

process bcAlign {
  tag "$meta.id - ${index}"
  label 'bedtools'
  label 'highCpu'
  label 'highMem'

  input:
  peaks

  output:
  tuple val(meta), path ("*_merged_peaks.bed"), emit: bed  
  path ("*_macs2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  cut -f1-3 ${peaks} | bedtools merge -d ${params.max_features_dist} -i /dev/stdin > ${prefix}_merged_peaks.bed 2>> ${prefix}_macs2.log
  """
}
