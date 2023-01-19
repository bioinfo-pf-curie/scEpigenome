/*
 * Add barcode tag into bam files
 */

process Macs2peakCalling {
  tag "$meta.id"
  label 'macs2'
  label 'highCpu'
  label 'highMem'

  input:
  bam

  output:
  broad_peaks_bed
  sharp_peaks_bed
  tuple val(meta), path (".*Peak")
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: '' 
  """
  macs2 callpeak -t ${bam} \
  ${args} \
  -n ${prefix} 2> ${prefix}_macs2.log
  """
}


# pseudo bulk

cut -f1-3 ${prefix}.*Peak > ${prefix}_coulmns.*Peaks

bedtools merge -d ${params.max_features_dist} -i /dev/stdin > ${prefix}_merged_peaks.bed 2>> ${prefix}_macs2.log

bedtools sort ${prefix}_merged_peaks.bed > ${prefix}_merged_peaks_sorted.bed