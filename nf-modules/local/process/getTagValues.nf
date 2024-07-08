/*
 * Extract a tag from a BAM file
 */

process getTagValues {
  tag "$meta.id"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path ("*barcodes.txt"), emit: barcodes
  tuple val(meta), path ("*counts.txt"), emit: counts
  path('versions.txt'), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  samtools view ${bam} |\
    awk -v tag=XB '{for(i=1;i<=NF;i++){if (\$i~tag){split(\$i,bc,":");print bc[length(bc)]}}}' |\
     sort | uniq -c > ${prefix}_final_barcodes_counts.txt
  awk '{print \$2}' ${prefix}_final_barcodes_counts.txt > ${prefix}_final_barcodes.txt
  echo \$(samtools --version | head -1) > versions.txt
  """
}
