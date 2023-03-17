process test{
  tag "$meta.id"
  label 'deeptools'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(fq), path(fq2)
  val(sf)
  val(effGenomeSize)
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

  echo ${fq2} >> ${prefix}".txt"

  """
}