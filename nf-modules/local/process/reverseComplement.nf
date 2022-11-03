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
  tuple val(meta), path(reads) // take only R2

  output:
  tuple val(meta), path('*_reverseComp.R2.fastq.gz')

  script:
  """
  fastx_reverse_complement -Q33 -i <(gzip -cd ${reads[1]}) -z -o ${meta}_reverseComp.R2.fastq.gz
  """
}