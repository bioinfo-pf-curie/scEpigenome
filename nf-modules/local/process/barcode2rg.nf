/*
 * Add barcodes info to read group information
 */

process barcode2rg {
  tag "$meta.id"
  label 'python'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai), path(barcodes)

  output:
  tuple val(meta), path("*wRG.bam"), path("*wRG.bam.bai"), emit: bam
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  barcode2rg.py -i ${bam} -o ${prefix}_wRG.bam --barcode ${barcodes} -SM ${meta.id} 2> ${prefix}_barcode2rg.log
  echo \$(python --version 2>&1) > versions.txt 
  """
}
