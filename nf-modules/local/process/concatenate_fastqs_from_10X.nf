/*
 * Concatenate 4 samples per reads 
 */

process concatenate_fastqs_from_10X{ 
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads)

  output:
  //tuple val(meta), path ("*_R2.fastq.gz"), emit: barcodeRead
  //tuple val(meta), path ("*_R1.fastq.gz"), path ("*_R3.fastq.gz"), emit: dnaRead
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo $prefix
  """
}