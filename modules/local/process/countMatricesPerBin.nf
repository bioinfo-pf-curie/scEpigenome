/*
 * Create binned sparse matrices
 */

process countMatricesPerBin {
  tag "${meta.id} [${bins}]"
  label 'python'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(bam), path(bai), path (bcList), val(bins)

  output:
  tuple val(meta), path ("*.tar.gz"), emit: matrix
  path ("versions.txt"), emit: versions
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}_bin_${bins}"
  def args = task.ext.args ?: ''
  def bcCmd = bcList ? "nbbarcodes=\$(wc -l < ${bcList})" : ''
  def bcOpts = bcList ? "-s \$nbbarcodes" : ""
  """
  ${bcCmd}
  sc2sparsecounts.py -i ${bam} -o ${prefix}_counts -b ${bins} ${bcOpts} ${args}
  tar -zcvf ${prefix}_counts.tar.gz ${prefix}_counts
  rm -rf ${prefix}_counts

  echo \$(python --version 2>&1) > versions.txt
  """
}