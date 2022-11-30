/*
 * Reverse Complement R2 
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process reverseComplement {
  label 'fastx'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(r2barcode) // take only R2

  output:
  tuple val(meta), path('*_reverseComp.R2.fastq.gz'), emit: reads
  path ("versions.txt"), emit: versions

  script:
  """
  fastx_reverse_complement -Q33 -i <(gzip -cd ${r2barcode}) -z -o ${meta}_reverseComp.R2.fastq.gz
  fastx_toolkit --version &> versions.txt
  """
}