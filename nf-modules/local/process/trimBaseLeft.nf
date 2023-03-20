/*
 * Seqtk trim R2 first base
 */

process trimBaseLeft {
  tag "$meta.id"
  label 'seqtk'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*_firstBaseTrim.R2.fastq.gz'), emit: reads
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  # Trim first base after reverseComp which is not part of the barcode
  seqtk trimfq -b 1 ${reads} > ${prefix}_firstBaseTrim.R2.fastq
  """
}
