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
  echo $fq > $prefix".txt"

  echo $fq2 >> $prefix".txt"

  """
}