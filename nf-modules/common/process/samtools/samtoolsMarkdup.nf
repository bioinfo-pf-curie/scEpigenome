/*
 * Samtools - Markdup
 */

process samtoolsMarkdup {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path ("*_markdup.bam"), emit: bam
  path ("*_markdup.log"), emit: logs

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${bam.baseName}"
  """
  # remove duplicates 
  samtools markdup -s -r \\
    ${bam} \\
    ${prefix}_markdup.bam &> ${prefix}_markdup.log
  """
}

