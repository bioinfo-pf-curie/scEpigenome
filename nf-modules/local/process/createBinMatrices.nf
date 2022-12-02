/*
 * Create binned sparse matrices
 */

process createBinMatrices {
  tag "$meta.id"
  label 'python'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai), path(nbBc), val(bins)

  output:
  tuple val(meta), path ("*.zip"), emit: matrix
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  # Counts per bin (--bin)
  nbbarcodes=\$(awk '{print \$1}' ${nbBc})
  sc2sparsecounts.py -i ${bam} -o ${prefix}_counts_bin_${bins} -b ${bins} -s \$nbbarcodes -v ${args}

  zip -r ${prefix}_counts_bin_${bins}.zip ${prefix}_counts_bin_${bins}
  rm -rf ${prefix}_counts_bin_${bins}

  python --version &> versions.txt
  """
}
