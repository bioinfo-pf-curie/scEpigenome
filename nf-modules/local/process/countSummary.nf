/* 
 *  Summary of all logs and counts
 */

process countSummary {
  tag "$meta.id"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(flaggedBam) from chRemovePcrRtDup_Log
  tuple val(meta), path(pcrDup) from chPCRdupCount
  tuple val(meta), path(rtDup) from chRTdupCount
  tuple val(meta), path(r1UnmappedR2) from chR1unmappedR2Count
  tuple val(meta), path(rmDupSam) from chCountSummary

  output:
  tuple val(meta), path("*_removePcrRtDup.log"), emit: logs //into chPcrRtCountsLog
  tuple val(meta), path("*_rmDup.txt"), emit: result //into chDistribUMIs, chRemoveDupBarcodeLog, chPerBin, chPerTSS
  tuple val(meta), path("*.count"), emit: count

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  ### 1. Count duplicate proportions
  n_mapped_barcoded=\$(samtools view -c ${flaggedBam})
  n_pcr_duplicates=\$(cat ${pcrDup})
  n_rt_duplicates=\$(cat ${rtDup})
  n_R1_mapped_R2_unmapped=\$(cat ${r1UnmappedR2})
  n_unique_except_R1_unmapped_R2=\$((\$n_mapped_barcoded - \$n_pcr_duplicates - \$n_rt_duplicates))

  echo "## Number of reads mapped and barcoded: \$n_mapped_barcoded" > ${prefix}_removePcrRtDup.log
  echo "## Number of pcr duplicates: \$n_pcr_duplicates" >> ${prefix}_removePcrRtDup.log
  echo "## Number of rt duplicates: \$n_rt_duplicates" >> ${prefix}_removePcrRtDup.log
  echo "## Number of R1 mapped but R2 unmapped: \$n_R1_mapped_R2_unmapped" >> ${prefix}_removePcrRtDup.log
  echo "## Number of reads after PCR and RT removal (not R1 unmapped R2): \$n_unique_except_R1_unmapped_R2" >> ${prefix}_removePcrRtDup.log

  ### 2. Count final number of barcoded cells 
  # Count nb barcodes from flagged - PCR, RT & window dups  (need to sort by barcode)
  barcode_field=\$( cat ${rmDupSam} | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)
  cat ${rmDupSam} | awk -v bc_field=\$barcode_field '{print substr(\$bc_field,6)}' | sort | uniq -c > ${prefix}_rmDup.txt
  barcodes=\$(wc -l ${prefix}_rmDup.txt | awk '{print \$1}')
  echo "Barcodes found = \$barcodes" > ${prefix}.count  
  """
}