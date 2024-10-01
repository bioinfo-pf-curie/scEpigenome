/*
 * Extract barcode sequences from reads
 */

process getBarcodes {
  tag "$meta.id"
  label 'onlyLinux'
  label 'minCpu'
  label 'minMem'

  input:
  tuple val(meta), path(reads), val(start), val(len)

  output:
  tuple val(meta), path ("*barcodeInReads.fasta"), emit: fasta
  path ("versions.txt"), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  gzip -cd  ${reads} | \
    awk -v start=${start} -v size=${len} 'NR%4==1{print ">"\$0}; NR%4==2{print substr(\$0,start,size)}' > ${prefix}_barcodeInReads.fasta

  ## version
  echo \$(awk --version | awk 'NR==1{print "awk "\$3}') > versions.txt
  """
}
