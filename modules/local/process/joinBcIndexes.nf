/*
 * Bowtie2 index barcode alignment
 */

process joinBcIndexes {
  tag "$meta.id"
  label 'onlyLinux'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(barcodeReads) 

  output:
  tuple val(meta), path("*_read_barcodes.txt"), emit: results

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputs = barcodeReads.sort()
  """
  paste ${inputs} |\
     awk '\$1!=\$3 || \$1!=\$5{exit -1}\
     { if(\$2=="None"||\$4=="None"||\$6=="None"){print \$1"\tNone"}else{print \$1"\t"\$2"-"\$4"-"\$6} }' > ${prefix}_read_barcodes.txt
  """
}