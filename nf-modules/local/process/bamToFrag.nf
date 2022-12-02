/*
 * Bam to fragment files (~scbed) files
 */

process bamToFrag {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)
  
  output:
  tuple val(meta), path("*.fragments.tsv.gz"), emit: gz

  script:
  def prefix = task.ext.prefix ?: "${bam.baseName}"
  """
  ##Sort by barcode then chromosome then position R2
  #Find the column containing the barcode tag XB
  barcode_field=\$(samtools view ${bam} |sed -n "1 s/XB.*//p" |sed 's/[^\t]//g' | wc -c)

  #Sort by barcode then chromosome then read position
  samtools view ${bam} | grep -E "XB:Z" | awk -v bc_field=\$barcode_field -v OFS="\t" '{gsub("XB:Z:","",\$bc_field); print \$3,\$4,\$4+100,\$bc_field,1}' > ${prefix}.fragments.tsv

  #Compress
  bgzip -@ 8 -f -l 9 ${prefix}.fragments.tsv

  ## Index flagged_rmPCR_RT file
  tabix -p bed ${prefix}.fragments.tsv.gz
  """
}
