/*
 * Reverse Complement R2 
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process reverseComplement {
  tag "$meta.id"
  label 'seqkit'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(r2barcode) // take only R2

  output:
  tuple val(meta), path('*_reverseComp.R2.fastq.gz'), emit: reads
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  seqkit seq ${r2barcode} --threads ${task.cpus} --reverse --complement --seq-type dna -o ${prefix}_reverseComp.R2.fastq.gz
  seqkit version  &> versions.txt
  """
}

