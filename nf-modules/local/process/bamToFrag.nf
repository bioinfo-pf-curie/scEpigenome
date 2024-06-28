/*
 * Bam to fragment files files
 */

process bamToFrag {
  tag "${meta.id}"
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam)
  
  output:
  tuple val(meta), path("*.fragments.tsv*"), emit: tsv
  path("*.log"), emit: log
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  bamToFrag.py ${args} --input ${bam} --output ${prefix}.fragments.tsv > ${prefix}_bam2frag.log 2>&1
  echo \$(python --version 2>&1) > versions.txt
  """
}
