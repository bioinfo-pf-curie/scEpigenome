/*
 * Macs2 - peak calling
 */

process test{
  tag "$meta.id"
  label 'macs2'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(fq), path(fq2)
  
  output:
  tuple val(meta), path('*.txt'), emit: txt
  
  errorStrategy 'ignore'

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo ${fq} > ${prefix}".txt"
  echo ${fq2} >> ${prefix}".txt"
  """
}


