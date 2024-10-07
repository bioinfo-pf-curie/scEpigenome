/* 
 *  Weighted histogram to highlight empty droplets 
 */

process weightedDistrib {
  tag "$meta.id"
  label 'R'
  label 'medCpu'
  label 'lowMem'
  
  input:
  tuple val(meta), path(counts)

  output:
  path("*distDF.mqc"), emit: mqc
  tuple val(meta), path("*distribution.pdf"), emit: pdf
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  weightedDistribution.r ${counts} ${prefix}
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  """
}