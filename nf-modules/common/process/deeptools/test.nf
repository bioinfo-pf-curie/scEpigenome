/*
 * Macs2 - peak calling
 */

process test{
  tag "$meta.id"
  label 'macs2'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(fq), val(effGenomeSize)
  val(sf)
  path(blacklistBed)

  output:
  tuple val(meta), path('*.txt'), emit: txt
  
  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  effGsizeOpts = effGenomeSize.size() ? "--effectiveGenomeSize ${effGenomeSize[0]}" : ""
  sfOpts = sf.size() ? "--scaleFactor $sf" : ""
  """
  echo ${fq} > ${prefix}".txt"
  """
}


