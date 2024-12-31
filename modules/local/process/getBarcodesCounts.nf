/*
 * Get counts per barcode
 */

process getTagfragmentCounts {
  tag "$meta.id"
  label 'onlyLinux'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(readsBarcodes)

  output:
  tuple val(meta), path ("*_counts.txt"), emit: counts
  tuple val(meta), path ("*_list.txt"), emit: barcodes
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  awk -F"\t" '\$2!="None"{print \$2}' ${readsBarcodes} | sort | uniq -c > ${prefix}_barcodes_counts.txt
  awk '{print \$2}' ${prefix}_barcodes_counts.txt > ${prefix}_barcodes_list.txt
  """
}
