/*
 * Get unique barcode sequences
 */

process getUniqueBarcodes {
  tag "$meta.id"
  label 'unix'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(readsBarcodes)

  output:
  tuple val(meta), path ("*_unique_barcodes.txt"), emit: barcodes
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  awk -F"\t" '\$2!=None{print \$2}' ${readsBarcodes} | sort -u > ${prefix}_unique_barcodes.txt
  """
}
