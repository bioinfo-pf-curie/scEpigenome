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
    awk -v tag=${params.barcodeTag} '
  {
      for(i=1; i<=NF; i++){
          if (\$i ~ tag":"){  # Vérifie que le champ contient bien "XB:"
              split(\$i, bc, ":"); 
              print bc[length(bc)];  # Récupère la dernière partie (valeur du tag XB)
          }
      }
  }' | \
  sort | uniq -c > ${prefix}_final_barcodes_counts.txt
  awk '{print \$2}' ${prefix}_final_barcodes_counts.txt > ${prefix}_final_barcodes.txt
  echo \$(samtools --version | head -1) > versions.txt
  """
}
