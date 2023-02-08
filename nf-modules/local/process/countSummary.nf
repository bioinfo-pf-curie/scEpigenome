/* 
 *  Summary of all logs and counts
 */

process countSummary {
  tag "$meta.id"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  
  input:
  tuple val(meta), path(pcrDup)
  tuple val(meta), path(flaggedBam)
  tuple val(meta), path(r1UnmappedR2)
  tuple val(meta), path(rtDup)

  output:
  tuple val(meta), path ("*_allDup.log"), emit: logs 

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  ### 1. Count duplicate proportions
  n_mapped_barcoded=\$(samtools view -c ${flaggedBam})
  n_pcr_duplicates=\$(cat ${pcrDup})
  n_R1_mapped_R2_unmapped=\$(cat ${r1UnmappedR2})

  echo "## Number of reads mapped and barcoded: \$n_mapped_barcoded" > ${prefix}_allDup.log
  echo "## Number of pcr duplicates: \$n_pcr_duplicates" >> ${prefix}_allDup.log
  echo "## Number of R1 mapped but R2 unmapped: \$n_R1_mapped_R2_unmapped" >> ${prefix}_allDup.log

  if [[ "${params.protocol}" == "scchip_indrop" ]] 
  then
    n_rt_duplicates=\$(cat ${rtDup})
    n_unique_except_R1_unmapped_R2=\$((\$n_mapped_barcoded - \$n_pcr_duplicates - \$n_rt_duplicates))
    echo "## Number of rt duplicates: \$n_rt_duplicates" >> ${prefix}_allDup.log
    echo "## Number of reads after PCR and RT removal (not R1 unmapped R2): \$n_unique_except_R1_unmapped_R2" >> ${prefix}_allDup.log
  fi
  """
}