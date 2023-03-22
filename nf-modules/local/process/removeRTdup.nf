/* 
 *  RT deduplaication Workflow
 */

process removeRTdup {
  tag "$meta.id"
  label 'samtools'
  label 'medCpu'
  label 'lowMem'
  
  input:
  tuple val(meta), path(flaggedBam)
  tuple val(meta), path(flaggedRmPCRbam)
  tuple val(meta), path(flaggedRmPCRsam)

  output:
  tuple val(meta), path("*_flagged_rmPCR_RT.bam"), emit: bam
  // For countSummary
  tuple val(meta), path("*_count_RT_duplicates.txt"), emit: logs

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  #Find the column containing the barcode tag XB
  barcode_field=\$(samtools view ${flaggedBam} | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)
  #Find the column containing the position R2 tag XS
  posR2_field=\$(samtools view ${flaggedBam} | sed -n \"1 s/XS.*//p\" | sed 's/[^\t]//g' | wc -c)

  ## Index flagged_rmPCR file
  samtools index ${flaggedRmPCRbam}

  if [ ${params.keepRTdup} == 'false' ] 
  then
    #Remove RT duplicates (if two consecutive reads have the same barcode and same R2 chr&start) but not same R1 
    cat ${flaggedRmPCRsam} | awk -v bc_field=\$barcode_field -v R2_field=\$posR2_field 'BEGIN{count=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\")} ; NR>=2{split( \$R2_field,R2Pos,\":\");if((R2Pos[3]==lastR2Pos[3]) && (R2Pos[3]!=2147483647) && (lastR2Pos[3]!=2147483647)  && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){count++;next} {print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print count > \"count_RT_duplicates\"}' > ${prefix}_flagged_rmPCR_RT.sam
    # sam to bam
    samtools view -H ${flaggedBam}  | sed '/^@CO/ d' > ${prefix}_header.sam
    cat ${prefix}_flagged_rmPCR_RT.sam >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged_rmPCR_RT.bam 
    ## Sort flagged_rmPCR_RT file
    samtools sort -@ ${task.cpus} ${prefix}_flagged_rmPCR_RT.bam > ${prefix}_flagged_rmPCR_RT_sorted.bam
    ## Rename flagged_rmPCR_RT file
    mv ${prefix}_flagged_rmPCR_RT_sorted.bam ${prefix}_flagged_rmPCR_RT.bam

    # If no RT duplicates removing:
    else
      ## Copy flagged_rmPCR to flagged_rmPCR_RT
      cp ${flaggedRmPCRbam} ${prefix}_flagged_rmPCR_RT.bam
      ## Set RT duplicate count to 0
      echo 0 > count_RT_duplicates
    fi

    # Save counts
    mv count_RT_duplicates ${prefix}_count_RT_duplicates.txt
  """
}