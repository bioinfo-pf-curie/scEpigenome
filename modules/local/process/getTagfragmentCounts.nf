/*
 * Extract a tag from a BAM file and count the number of fragment
 */

process getTagfragmentCounts {
  tag "$meta.id"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path ("*barcodes.txt"), emit: barcodes
  tuple val(meta), path ("*counts.txt"), emit: counts
  path('versions.txt'), emit: versions
 
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  samtools view ${bam} |\
    awk -v tag=XB 'NR==1{for(i=1;i<=NF;i++){if (\$i~tag) print \$i; read=\$1} \$1!=read;print \$i;read=\$1}' |\
    cut -f3 -d: | sort | uniq -c > ${prefix}_final_barcodes_counts.txt
  awk '{print \$2}' ${prefix}_final_barcodes_counts.txt > ${prefix}_final_barcodes.txt
  echo \$(samtools --version | head -1) > versions.txt
  """
}
