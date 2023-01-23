/*
 * Bowtie2 index barcode alignment
 */

process bedtoolsMergePeaks {
  tag "$meta.id"
  label 'macs2'
  label 'highCpu'
  label 'highMem'

  input:
  bed

  output:
  tuple val(meta), path ("*_merged_peaks.bed"), emit: bed  
  path ("*_macs2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  cut -f1-3 ${bed} | bedtools merge ${args} -i /dev/stdin > ${prefix}_merged_peaks.bed 2>> ${prefix}_macs2.log
  bedtools sort ${prefix}_merged_peaks.bed > ${prefix}_merged_peaks_sorted.bed
  """
}
