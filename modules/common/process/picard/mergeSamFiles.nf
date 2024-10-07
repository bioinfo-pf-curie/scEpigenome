/*
 * Picard Merge SAM/BAM files
 * Compared to samtools merge, it alsom merge the RG tags
 */

process mergeSamFiles {
  tag "${meta.id}"
  label 'picard'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bams), path(bai)

  output:
  tuple val(meta), path('*.bam'), path('*.bai'), emit: bam
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def javaArgs = task.ext.args ?: ''
  def args = task.ext.args2 ?: ''
  memOption = "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
  """
  picard ${memOption} ${javaArgs} MergeSamFiles \\
    --CREATE_INDEX true \\
    ${args} \\
    ${'--INPUT '+bams.join(' --INPUT ')} \\
    --OUTPUT ${prefix}.bam

  echo \$(picard CollectInsertSizeMetrics --version 2>&1 | grep Version | sed -e 's/.*Version:/picard /') > versions.txt
  """
}
