/*
 * Peaks Annotation with HOMER
 */

process annotatePeaks {
  tag "${meta.id}"
  label 'homer'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(peaks)
  path gtf
  path fasta

  output:
  tuple val(meta), path("*.txt"), emit: output

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  annotatePeaks.pl ${peaks} \\
        ${fasta} \\
        -gtf ${gtf} \\
        -cpu ${task.cpus} \\
        > ${prefix}_annotHomer.txt
  """
}


