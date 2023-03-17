process test {
  tag "${meta.id}"
  label 'deeptools'
  label 'highCpu'
  label 'medMem'

  input:
  tuple val(meta), path(fq), path(fq2)
  val(sf)
  val(effGenomeSize)
  path(blacklistBed)
  
  output:
  tuple val(meta), path('*.txt'), emit: txt

  when:
  task.ext.when == null || task.ext.when

  script:
  blacklistOpts = blacklistBed.size() ? "--blackListFileName ${blacklistBed}" : ""
  effGsizeOpts = effGenomeSize.size() ? "--effectiveGenomeSize ${effGenomeSize[0]}" : ""
  sfOpts = sf.size() ? "--scaleFactor $sf" : ""
  strandOpts = meta.strandness == 'forward' ? '--filterRNAstrand forward' : meta.strandness == 'reverse' ? '--filterRNAstrand reverse' : ''
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo $prefix > $prefix".txt"

  """
}