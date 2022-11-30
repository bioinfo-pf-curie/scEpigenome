/*
 * Reverse Complement R2 
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process reverseComplement {
  tag "$meta.id"
  label 'fastx'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(r2barcode) // take only R2

  output:
  tuple val(meta), path('*_reverseComp.R2.fastq.gz'), emit: reads
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  fastx_reverse_complement -Q33 -i <(gzip -cd ${r2barcode}) -z -o ${prefix}_reverseComp.R2.fastq.gz
  fastx_toolkit --version &> versions.txt
  """
}