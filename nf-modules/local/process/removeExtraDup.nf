/* 
 *  Remove RT and window duplicates
 */

process removeExtraDup {
  tag "$meta.id"
  label 'python'
  label 'medCpu'
  label 'lowMem'
  
  input:
  tuple val(meta), path(bam)
  
  output:
  tuple val(meta), path("*_extraDup.bam"), emit: bam
  path("*_removeExtraDup.log"), emit: logs
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  rmDup.py -v ${args} -i ${bam} -o ${prefix}_extraDup.bam > ${prefix}_removeExtraDup.log
  echo \$(python --version 2>&1) > versions.txt
  """
}