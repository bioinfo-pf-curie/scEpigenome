/*
 * Add barcodes info in reads names of the DNA part of the fastqs
 */

process addBarcodes {
  tag "$meta.id"
  label 'python'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(reads), path(barcodes)

  output:
  tuple val(meta), path ("*_barcoded.fastq.gz"), emit: fastq
  path ("*.log"), emit: log
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  addBarcode.py --fastq1 ${reads[0]} --fastq2 ${reads[1]} --barcode ${barcodes} 2> ${prefix}_addbarcodes.log
  echo \$(python --version 2>&1) > versions.txt 
  """
}
