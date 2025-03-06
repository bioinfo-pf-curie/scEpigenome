/*
 * Merge all barcodes info
 */

process mergeBarcodes {
  tag "$meta.id"
  label 'onlyLinux'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bclist), path(bccounts)

  output:
  tuple val(meta), path ("*barcodes.txt"), emit: barcodes
  tuple val(meta), path ("*counts.txt"), emit: counts
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  cat ${bclist} | sort -u -T ${params.tmpDir} --parallel=${task.cpus} > ${prefix}_barcodes.txt 
  cat ${bccounts} | sort -T ${params.tmpDir} --parallel=${task.cpus} -k2,2 | \
    awk 'NR==1{bc=\$2;c=\$1} NR>1{if(\$2==bc){c+=\$1}else{print c,bc; bc=\$2; c=\$1}}END{print c,bc}' > ${prefix}_counts.txt
  """
}
