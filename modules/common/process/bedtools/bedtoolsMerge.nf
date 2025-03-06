/*
 * Fusionner des intervalles BED chevauchants
 */

process bedtoolsMerge {
  tag "$meta.id"
  label 'bedtools'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bed)

  output:
  tuple val(meta), path ("*_merged.bed"), emit: bed  
  path ("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  bedtools merge ${args} -i ${bed} > ${prefix}_merged.bed
  echo \$(bedtools --version echo) &> versions.txt
  """
}
