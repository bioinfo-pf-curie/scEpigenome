/*
 * Create sparse matrices on features
 */

process countMatricesPerFeature {
  tag "${meta.id}"
  label 'python'
  label 'medCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai), path(bcList)
  path(bed)

  output:
  tuple val(meta), path ("*.tar.gz"), emit: matrix
  tuple val(meta), path ("*_counts.logs"), emit:logs
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}_TSS"
  def args = task.ext.args ?: ''
  def bcCmd = bcList ? "nbbarcodes=\$(wc -l < ${bcList})" : ''
  def bcOpts = bcList ? "-s \$nbbarcodes" : ""
  """
  ${bcCmd}
  sc2sparsecounts.py -i ${bam} -o ${prefix}_counts -B ${bed} ${bcOpts} ${args} 2> ${prefix}_counts.logs
  tar -zcvf ${prefix}_counts.tar.gz ${prefix}_counts
  rm -rf ${prefix}_counts

  python --version &> versions.txt
  """
}
