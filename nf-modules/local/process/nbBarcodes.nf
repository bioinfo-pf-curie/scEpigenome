/*
 * Get number of barcodes
 */

process nbBarcodes {
  tag "$meta.id"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  path (bcList)

  output:
  tuple val(meta), path("*_nbBarcodes.txt"), emit: count
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  barcodes=\$(wc -l ${bcList} | awk '{print \$1}')
  echo \$barcodes > ${prefix}_nbBarcodes.txt
  """
}
