/* 
 *  PCR deduplaication Workflow
 */

process removePCRdup_cut {
  tag "$meta.id"
  label 'samtools'
  label 'medCpu'
  label 'lowMem'
  
  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*_flagged_rmPCR.bam"), emit: bam
  tuple val(meta), path("*_flagged_rmPCR.sam"), emit: sam
  // For countSummary
  tuple val(meta), path("*_count_PCR_duplicates.txt"), emit: count
  tuple val(meta), path("*_countR1unmappedR2.txt"), emit: countR1unmapped

  errorStrategy 'ignore'

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  #Find the column containing the barcode tag XB
  barcode_field=\$(samtools view ${bam} | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)
  #Find the column containing the position R2 tag XS
  posR2_field=\$(samtools view ${bam} | sed -n \"1 s/XS.*//p\" | sed 's/[^\t]//g' | wc -c)

  #Sort by barcode then chromosome then position R2 then Position R1 (for the PCR removal) 
  #It is important to sort by R1 pos also	because the removal is done by comparing consecutive lines
  printf '@HD\tVN:1.4\tSO:unsorted\n' > ${prefix}_header.sam
  samtools view -H ${bam} | sed '/^@HD/ d' >> ${prefix}_header.sam
  samtools view ${bam} | LC_ALL=C sort -T ${params.tmpDir} --parallel=${task.cpus} -t \$'\t' -k \"\$barcode_field.6" -k 3.4,3g -k 4,4n >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged.sorted.bam
  
  # counts and remove PCR duplicates
  samtools view ${prefix}_flagged.sorted.bam | awk -v bc_field=\$barcode_field 'BEGIN {countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; lastR1Pos=\$4} ; NR>=2{R1Pos=\$4; if( (R1Pos==lastR1Pos) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field}} END {print countPCR > \"count_PCR_duplicates\"}' > ${prefix}_flagged_rmPCR.sam

  # Save counts
  mv count_PCR_duplicates ${prefix}_count_PCR_duplicates.txt
  # no RT duplicates for cut&tag
  echo 0 > ${prefix}_countR1unmappedR2.txt

  # sam to bam
  samtools view -H ${prefix}_flagged.sorted.bam  | sed '/^@CO/ d' > ${prefix}_header.sam
  cat ${prefix}_flagged_rmPCR.sam >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged_rmPCR.bam

  ## Sort flagged_rmPCR file
  samtools sort -@ ${task.cpus} ${prefix}_flagged_rmPCR.bam > ${prefix}_flagged_rmPCR_sorted.bam

  ## Rename flagged_rmPCR file
  mv ${prefix}_flagged_rmPCR_sorted.bam ${prefix}_flagged_rmPCR.bam

  rm ${prefix}_header.sam 
  """
}