/* 
 *  PCR deduplaication Workflow
 */

process removePCRdup {
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
  tuple val(meta), path("*_flagged.sorted.bam"), emit: bamLogs

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
  samtools view ${bam} | LC_ALL=C sort -T ${params.tmpDir} --parallel=${task.cpus} -t \$'\t' -k \"\$barcode_field.8,\$barcode_field\"n -k 3.4,3g -k \"\$posR2_field.6,\$posR2_field\"n -k 4,4n >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged.sorted.bam
  
  # counts
  samtools view ${prefix}_flagged.sorted.bam | awk -v bc_field=\$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${prefix}_flagged.count
  samtools view ${prefix}_flagged.sorted.bam | awk -v bc_field=\$barcode_field -v R2_field=\$posR2_field 'BEGIN {countR1unmappedR2=0;countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\"); lastR1Pos=\$4} ; NR>=2{split( \$R2_field,R2Pos,\":\");R1Pos=\$4; if(R2Pos[3]==2147483647){print \$0;countR1unmappedR2++; next}; if( (R1Pos==lastR1Pos) && (R2Pos[3]==lastR2Pos[3]) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print countPCR > \"count_PCR_duplicates\";print countR1unmappedR2 > \"countR1unmappedR2\"}' > ${prefix}_flagged_rmPCR.sam

  # Save counts
  mv count_PCR_duplicates ${prefix}_count_PCR_duplicates.txt
  mv countR1unmappedR2 ${prefix}_countR1unmappedR2.txt

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