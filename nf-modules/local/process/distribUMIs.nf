/* 
 *  Weighted histogram to highlight empty droplets 
 */

process distribUMIs {
  tag "$meta.id"
  label 'R'
  label 'medCpu'
  label 'lowMem'
  
  input:
  tuple val(meta), path(countedReadsPerCell_matrix) 

  output:
  path("*distDF.mqc"), emit: mqc
  tuple val(meta), path("*distribution.pdf"), emit: pdf
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  umisDistribution.r ${countedReadsPerCell_matrix} ${prefix}
  R --version &> versions.txt
  """
}