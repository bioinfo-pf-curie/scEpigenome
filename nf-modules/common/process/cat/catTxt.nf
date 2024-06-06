/*
 * Cat text files
 */

process catTxt {
  tag "$meta.id"
  label 'unix'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(inputs)
  val(unique)

  output:
  tuple val(meta), path ("*txt"), emit: txt
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def opts = unique ? " | sort -u" : ''
  """
  cat ${inputs} ${opts} > ${prefix}.txt 
  """
}
