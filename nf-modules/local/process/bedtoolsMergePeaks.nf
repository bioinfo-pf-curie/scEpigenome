/*
 * Bowtie2 index barcode alignment
 */

process bedtoolsMergePeaks {
  tag "$meta.id"
  label 'macs2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(peaks)

  output:
  tuple val(meta), path ("*_merged_peaks_sorted.bed"), emit: bed  
  path ("*_macs2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  bedtools merge ${args} -i ${peaks} > ${prefix}_merged_peaks.bed 2>> ${prefix}_macs2.log
  bedtools sort -i ${prefix}_merged_peaks.bed > ${prefix}_merged_peaks_sorted.bed
  echo \$(bedtools --version echo) &> versions.txt
  """
}
