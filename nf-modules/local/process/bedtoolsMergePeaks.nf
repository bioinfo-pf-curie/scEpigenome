/*
 * Bowtie2 index barcode alignment
 */

process bedtoolsMergePeaks {
  tag "$meta.id"
  label 'bedtools'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(peaks)

  output:
  tuple val(meta), path ("*.bed"), emit: bed  
  path ("*_macs2.log"), emit: logs
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  peak_type=\$(echo ${peaks} | cut -f2 -d.)
  cut -f1-3 ${peaks} | bedtools merge ${args} | bedtools sort > ${prefix}_merged_peaks_sorted."\$peak_type".bed 2>> ${prefix}_"\$peak_type"_macs2.log
  echo \$(bedtools --version echo) &> versions.txt
  """
}
