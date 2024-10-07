/*
 * Reverse complement a sequence
 */

process reverseComplement {
  tag "$meta.id"
  label 'seqkit'
  label 'highCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*_reverseComp.fastq.gz'), emit: reads
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  seqkit seq ${reads} --threads ${task.cpus} --reverse --complement --seq-type dna -o ${prefix}_reverseComp.fastq.gz
  seqkit version  &> versions.txt
  """
}

