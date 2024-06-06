/*
 * Create binned sparse matrices
 */

process countMatricesPerBin {
  tag "${meta.id}"
  label 'python'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(bam), path(bai), path (bcList), val(bins)

  output:
  tuple val(meta), path ("*.tar.gz"), emit: matrix
  path ("versions.txt"), emit: versions
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  nbbarcodes=\$(wc -l < ${bcList})

  # Counts per bin (--bin)
  sc2sparsecounts.py -i ${bam} -o ${prefix}_counts_bin_${bins} -b ${bins} -s \$nbbarcodes -v ${args}

  tar -zcvf ${prefix}_counts_bin_${bins}.tar.gz ${prefix}_counts_bin_${bins}
  rm -rf ${prefix}_counts_bin_${bins}

  python --version &> versions.txt
  """
}