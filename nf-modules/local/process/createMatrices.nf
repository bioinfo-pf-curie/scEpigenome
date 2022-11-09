/*
 * Create binned sparse matrices
 */

process createMatrices {
  tag "$meta.id"
  label 'python'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(nbBc)
  tuple val(meta), path(bam), path(bai), val(bins)

  output:
  tuple val(meta), path ("*.zip"), emit: matrix
  path ("versions.txt"), emit: versions

 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  # Counts per bin (--bin)
  sc2sparsecounts.py -i ${rmDupBam} -o ${prefix}_counts_bin_${bins} -b ${bins} -s \$barcodes -v ${args}

  zip -r ${prefix}_counts_bin_${bins}.zip ${prefix}_counts_bin_${bins}
  rm -rf ${prefix}_counts_bin_${bins}

  python --version &> versions.txt
  """
}
