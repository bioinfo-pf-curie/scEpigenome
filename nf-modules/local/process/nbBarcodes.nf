/*
 * Get number of barcodes
 */

process nbBarcodes {
  tag "${meta.id}"
  label 'python'
  label 'medCpu'
  label 'medMem'

  input:
  path(bcList)

  output:
  tuple val(meta), path("*_nbBarcodes.txt"), emit: count
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  barcodes=\$(wc -l ${bcList} | awk '{print \$1}')
  echo \$barcodes > ${prefix}_nbBarcodes.txt
  """
}
